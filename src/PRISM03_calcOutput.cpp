// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// Prismatic is distributed under the GNU General Public License (GPL)
// If you use Prismatic, we kindly ask that you cite the following papers:

// 1. Ophus, C.: A fast image simulation algorithm for scanning
//    transmission electron microscopy. Advanced Structural and
//    Chemical Imaging 3(1), 13 (2017)

// 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
//    Implementation of Image Simulation Algorithms for Scanning
//	  Transmission Electron Microscopy. arXiv:1706.08563 (2017)

#include "PRISM03_calcOutput.h"
#include "params.h"
#include <iostream>
#include <algorithm>
#include <string>
#include <cmath>
#include <thread>
#include <mutex>
#include <numeric>
#include <vector>
#include "fftw3.h"
#include "utility.h"
#include "WorkDispatcher.h"
#include "ArrayND.h"
#include "fileIO.h"

#ifdef PRISMATIC_BUILDING_GUI
#include "prism_progressbar.h"
#endif

namespace Prismatic
{
extern std::mutex fftw_plan_lock; // for synchronizing access to shared FFTW resources
// extern std::mutex HDF5_lock;
using namespace std;
const static std::complex<PRISMATIC_FLOAT_PRECISION> i(0, 1);
// this might seem a strange way to get pi, but it's slightly more future proof
const static PRISMATIC_FLOAT_PRECISION pi = std::acos(-1);
Array2D<PRISMATIC_FLOAT_PRECISION> array2D_subset(const Array2D<PRISMATIC_FLOAT_PRECISION> &arr,
												  const size_t &starty, const size_t &stepy, const size_t &stopy,
												  const size_t &startx, const size_t &stepx, const size_t &stopx)
{
	vector<PRISMATIC_FLOAT_PRECISION> _d;
	size_t dimx = 0;
	size_t dimy = 0;
	for (auto y = starty; y < stopy; y += stepy)
	{
		for (auto x = startx; x < stopx; x += stepx)
		{
			_d.push_back(arr.at(y, x));
		}
		++dimy;
	}
	for (auto x = startx; x < stopx; x += stepx)
		++dimx;
	Array2D<PRISMATIC_FLOAT_PRECISION> result(_d, {{dimy, dimx}});
	return result;
}

void setupCoordinates_2(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	std::vector<PRISMATIC_FLOAT_PRECISION> xp_d;
	std::vector<PRISMATIC_FLOAT_PRECISION> yp_d;

	if(pars.meta.arbitraryProbes)
	{
		xp_d = pars.meta.probes_x;
		yp_d = pars.meta.probes_y;
		pars.numXprobes = xp_d.size();
		pars.numYprobes = 1;
		pars.numProbes = xp_d.size();
	}
	else
	{
		Array1D<PRISMATIC_FLOAT_PRECISION> xR = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{2}});
		xR[0] = pars.scanWindowXMin * pars.tiledCellDim[2];
		xR[1] = pars.scanWindowXMax * pars.tiledCellDim[2];
		Array1D<PRISMATIC_FLOAT_PRECISION> yR = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{2}});
		yR[0] = pars.scanWindowYMin * pars.tiledCellDim[1];
		yR[1] = pars.scanWindowYMax * pars.tiledCellDim[1];

		PRISMATIC_FLOAT_PRECISION probeStepX;
		PRISMATIC_FLOAT_PRECISION probeStepY;
		if (pars.meta.nyquistSampling)
		{
			int numX = nyquistProbes(pars, 2); //x is dim 2
			int numY = nyquistProbes(pars, 1); //y is dim 1
			probeStepX = pars.tiledCellDim[2] / numX;
			probeStepY = pars.tiledCellDim[1] / numY;
		}
		else
		{
			probeStepX = pars.meta.probeStepX;
			probeStepY = pars.meta.probeStepY;
		}

		xp_d = vecFromRange(xR[0], probeStepX, xR[1]);
		yp_d = vecFromRange(yR[0], probeStepY, yR[1]);
		pars.numXprobes = xp_d.size();
		pars.numYprobes = yp_d.size();
		pars.numProbes = xp_d.size()*yp_d.size();
	}

	Array1D<PRISMATIC_FLOAT_PRECISION> xp(xp_d, {{xp_d.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> yp(yp_d, {{yp_d.size()}});

	pars.xp = xp;
	pars.yp = yp;
}

void setupDetector(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	// setup coordinates related to detector size, angles, and image output

	pars.alphaMax = pars.qMax * pars.lambda;

	vector<PRISMATIC_FLOAT_PRECISION> detectorAngles_d = vecFromRange(pars.meta.detectorAngleStep * 1000 / 2,
																	  pars.meta.detectorAngleStep * 1000,
																	  (pars.alphaMax - pars.meta.detectorAngleStep / 2) * 1000);
	Array1D<PRISMATIC_FLOAT_PRECISION> detectorAngles(detectorAngles_d, {{detectorAngles_d.size()}});
	pars.detectorAngles = detectorAngles;
	PRISMATIC_FLOAT_PRECISION r_0 = pars.imageSizeOutput[0] / pars.meta.interpolationFactorY / 2;
	PRISMATIC_FLOAT_PRECISION r_1 = pars.imageSizeOutput[1] / pars.meta.interpolationFactorX / 2;
	vector<PRISMATIC_FLOAT_PRECISION> yVec_d = vecFromRange(-r_0, (PRISMATIC_FLOAT_PRECISION)1.0, r_0 - 1);
	vector<PRISMATIC_FLOAT_PRECISION> xVec_d = vecFromRange(-r_1, (PRISMATIC_FLOAT_PRECISION)1.0, r_1 - 1);
	Array1D<PRISMATIC_FLOAT_PRECISION> yVec(yVec_d, {{yVec_d.size()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> xVec(xVec_d, {{xVec_d.size()}});
	pars.yVec = yVec;
	pars.xVec = xVec;
	pars.Ndet = pars.detectorAngles.size();
}

void setupBeams_2(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	// setup some coordinates for the beams

	Array2D<PRISMATIC_FLOAT_PRECISION> beamsReduce = array2D_subset(pars.beamsOutput,
																	0, pars.meta.interpolationFactorY,
																	pars.beamsOutput.get_dimj(),
																	0, pars.meta.interpolationFactorX,
																	pars.beamsOutput.get_dimi());

	vector<size_t> imageSizeReduce{beamsReduce.get_dimj(), beamsReduce.get_dimi()};
	pars.imageSizeReduce = imageSizeReduce;
	pars.xyBeams = zeros_ND<2, long>({{pars.beamsIndex.size(), 2}});

	for (auto a0 = 1; a0 <= pars.beamsIndex.size(); ++a0)
	{
		for (auto y = 0; y < beamsReduce.get_dimj(); ++y)
		{
			for (auto x = 0; x < beamsReduce.get_dimi(); ++x)
			{
				if (beamsReduce.at(y, x) == a0)
				{
					pars.xyBeams.at(a0 - 1, 0) = y;
					pars.xyBeams.at(a0 - 1, 1) = x;
				}
			}
		}
	}
}

void createStack_integrate(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	// create output of a size corresponding to 3D mode (integration)
	pars.numLayers = 1;
	pars.output = zeros_ND<4, PRISMATIC_FLOAT_PRECISION>({{1, pars.numYprobes, pars.numXprobes, pars.Ndet}});
	if (pars.meta.saveDPC_CoM)
		pars.DPC_CoM = zeros_ND<4, PRISMATIC_FLOAT_PRECISION>({{1, pars.numYprobes, pars.numXprobes, 2}});
	
	std::vector<PRISMATIC_FLOAT_PRECISION> depths(1);
	depths[0] = pars.numPlanes*pars.meta.sliceThickness;
	pars.depths = depths;
	if (pars.meta.save4DOutput)
	{
		if(pars.fpFlag == 0 or pars.meta.saveComplexOutputWave) setup4DOutput(pars);
		
		if(pars.meta.saveComplexOutputWave)
		{
			pars.cbed_buffer_c = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>({{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
			if(pars.meta.crop4DOutput) pars.cbed_buffer_c = cropOutput(pars.cbed_buffer_c, pars);
		}
		else
		{
			pars.cbed_buffer = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
			if(pars.meta.crop4DOutput) pars.cbed_buffer = cropOutput(pars.cbed_buffer, pars);
		}
	}
}

void setupFourierCoordinates(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	// create Fourier space coordinates
	pars.qxaReduce = array2D_subset(pars.qxaOutput,
									0, pars.meta.interpolationFactorY, pars.qxaOutput.get_dimj(),
									0, pars.meta.interpolationFactorX, pars.qxaOutput.get_dimi());
	pars.qyaReduce = array2D_subset(pars.qyaOutput,
									0, pars.meta.interpolationFactorY, pars.qyaOutput.get_dimj(),
									0, pars.meta.interpolationFactorX, pars.qyaOutput.get_dimi());

	//grabbing 1D vector of fourier coordinates for output storage
	Array1D<PRISMATIC_FLOAT_PRECISION> qx = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.qxaReduce.get_dimi()}});
	Array1D<PRISMATIC_FLOAT_PRECISION> qy = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.qyaReduce.get_dimj()}});
	for (auto i = 0; i < pars.qxaReduce.get_dimi(); i++)
		qx[i] = pars.qxaReduce.at(0, i);
	for (auto j = 0; j < pars.qyaReduce.get_dimj(); j++)
		qy[j] = pars.qyaReduce.at(j, 0);

	pars.qx = qx;
	pars.qy = qy;

	pars.q1 = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
	pars.q2 = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>({{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
}

std::pair<Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>>, Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>>>
getSinglePRISMProbe_CPU(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const PRISMATIC_FLOAT_PRECISION xp, const PRISMATIC_FLOAT_PRECISION yp)
{

	// compute a single probe (for comparison with multislice in the GUI)

	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> realspace_probe;
	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> kspace_probe;

	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> psi = Prismatic::zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>(
		{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
	unique_lock<mutex> gatekeeper(fftw_plan_lock);
	PRISMATIC_FFTW_PLAN plan = PRISMATIC_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
														  reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi[0]),
														  reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi[0]),
														  FFTW_FORWARD, FFTW_ESTIMATE);
	gatekeeper.unlock();
	const static std::complex<PRISMATIC_FLOAT_PRECISION> i(0, 1);
	const static PRISMATIC_FLOAT_PRECISION pi = std::acos(-1);

	// setup some coordinates
	PRISMATIC_FLOAT_PRECISION x0 = xp / pars.pixelSizeOutput[1];
	PRISMATIC_FLOAT_PRECISION y0 = yp / pars.pixelSizeOutput[0];

	Array1D<PRISMATIC_FLOAT_PRECISION> x = pars.xVec + round(x0);
	// the second call to fmod here is to make sure the result is positive
	transform(x.begin(), x.end(), x.begin(), [&pars](PRISMATIC_FLOAT_PRECISION &a) {
		return fmod((PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[1] +
						fmod(a, (PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[1]),
					(PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[1]);
	});
	Array1D<PRISMATIC_FLOAT_PRECISION> y = pars.yVec + round(y0);
	transform(y.begin(), y.end(), y.begin(), [&pars](PRISMATIC_FLOAT_PRECISION &a) {
		return fmod((PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[0] +
						fmod(a, (PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[0]),
					(PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[0]);
	});
	Array2D<PRISMATIC_FLOAT_PRECISION> intOutput = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>(
		{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});

	memset(&psi[0], 0, sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>) * psi.size());
	for (auto a4 = 0; a4 < pars.beamsIndex.size(); ++a4)
	{
		PRISMATIC_FLOAT_PRECISION yB = pars.xyBeams.at(a4, 0);
		PRISMATIC_FLOAT_PRECISION xB = pars.xyBeams.at(a4, 1);

		if (abs(pars.psiProbeInit.at(yB, xB)) > 0)
		{
			PRISMATIC_FLOAT_PRECISION q0_0 = pars.qxaReduce.at(yB, xB);
			PRISMATIC_FLOAT_PRECISION q0_1 = pars.qyaReduce.at(yB, xB);
			std::complex<PRISMATIC_FLOAT_PRECISION> phaseShift = exp(
				-2 * pi * i * (q0_0 * (xp + pars.xTiltShift) + q0_1 * (yp + pars.yTiltShift)));
			const std::complex<PRISMATIC_FLOAT_PRECISION> tmp_const = pars.psiProbeInit.at(yB, xB) * phaseShift;
			auto psi_ptr = psi.begin();
			for (auto j = 0; j < y.size(); ++j)
			{
				for (auto i = 0; i < x.size(); ++i)
				{
					*psi_ptr++ += (tmp_const * pars.Scompact.at(a4, y[j], x[i]));
				}
			}
		}
	}
	realspace_probe = psi;
	PRISMATIC_FFTW_EXECUTE(plan);
	kspace_probe = psi;
	gatekeeper.lock();
	PRISMATIC_FFTW_DESTROY_PLAN(plan);
	gatekeeper.unlock();
	return std::make_pair(realspace_probe, kspace_probe);
}

void buildPRISMOutput_CPUOnly(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{

	// launch threads to compute results for batches of xp, yp
	// I do this by dividing the xp points among threads, and each computes
	// all of the relevant yp for each of its xp. This seems an okay strategy
	// as long as the number of xp and yp are similar.
	// If that is not the case
	// this may need to be adapted

	// initialize FFTW threads
	PRISMATIC_FFTW_INIT_THREADS();
	PRISMATIC_FFTW_PLAN_WITH_NTHREADS(pars.meta.numThreads);
	vector<thread> workers;
	workers.reserve(pars.meta.numThreads);																  // prevents multiple reallocations
	const size_t PRISMATIC_PRINT_FREQUENCY_PROBES = max((size_t)1, pars.numProbes / 10); // for printing status
	WorkDispatcher dispatcher(0, pars.numProbes);
	for (auto t = 0; t < pars.meta.numThreads; ++t)
	{
		cout << "Launching CPU worker thread #" << t << " to compute partial PRISM result\n";
		workers.push_back(thread([&pars, &dispatcher, &PRISMATIC_PRINT_FREQUENCY_PROBES]() {
			size_t Nstart, Nstop, ay, ax;
			Nstart = Nstop = 0;
			if (dispatcher.getWork(Nstart, Nstop))
			{ // synchronously get work assignment
				Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> psi = Prismatic::zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>(
					{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});

				unique_lock<mutex> gatekeeper(fftw_plan_lock);
				PRISMATIC_FFTW_PLAN plan = PRISMATIC_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
																	  reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi[0]),
																	  reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&psi[0]),
																	  FFTW_FORWARD, FFTW_MEASURE);
				gatekeeper.unlock();

				// main work loop
				do
				{
					while (Nstart < Nstop)
					{
						if (Nstart % PRISMATIC_PRINT_FREQUENCY_PROBES == 0 | Nstart == 100)
						{
							cout << "Computing Probe Position #" << Nstart << "/" << pars.numProbes << endl;
						}
						ay = (pars.meta.arbitraryProbes) ? Nstart : Nstart / pars.numXprobes;
						ax = (pars.meta.arbitraryProbes) ? Nstart : Nstart % pars.numXprobes;
						buildSignal_CPU(pars, ay, ax, plan, psi);
#ifdef PRISMATIC_BUILDING_GUI
						pars.progressbar->signalOutputUpdate(Nstart, pars.numProbes);
#endif
						++Nstart;
					}
				} while (dispatcher.getWork(Nstart, Nstop));
				gatekeeper.lock();
				PRISMATIC_FFTW_DESTROY_PLAN(plan);
				gatekeeper.unlock();
			}
		}));
	}
	// synchronize
	cout << "Waiting for threads...\n";
	for (auto &t : workers)
		t.join();
	PRISMATIC_FFTW_CLEANUP_THREADS();
}

void buildSignal_CPU(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
					 const size_t &ay,
					 const size_t &ax,
					 PRISMATIC_FFTW_PLAN &plan,
					 Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> &psi)
{
	// build the output for a single probe position using CPU resources

	const static std::complex<PRISMATIC_FLOAT_PRECISION> i(0, 1);
	const static PRISMATIC_FLOAT_PRECISION pi = std::acos(-1);

	// setup some coordinates
	PRISMATIC_FLOAT_PRECISION x0 = pars.xp[ax] / pars.pixelSizeOutput[1];
	PRISMATIC_FLOAT_PRECISION y0 = pars.yp[ay] / pars.pixelSizeOutput[0];
	Array1D<PRISMATIC_FLOAT_PRECISION> x = pars.xVec + round(x0);

	// the second call to fmod here is to make sure the result is positive
	transform(x.begin(), x.end(), x.begin(), [&pars](PRISMATIC_FLOAT_PRECISION &a) {
		return fmod((PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[1] +
						fmod(a, (PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[1]),
					(PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[1]);
	});

	Array1D<PRISMATIC_FLOAT_PRECISION> y = pars.yVec + round(y0);
	transform(y.begin(), y.end(), y.begin(), [&pars](PRISMATIC_FLOAT_PRECISION &a) {
		return fmod((PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[0] +
						fmod(a, (PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[0]),
					(PRISMATIC_FLOAT_PRECISION)pars.imageSizeOutput[0]);
	});

	Array2D<PRISMATIC_FLOAT_PRECISION> intOutput = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>(
		{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});

	memset(&psi[0], 0, sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>) * psi.size());

	for (auto a4 = 0; a4 < pars.beamsIndex.size(); ++a4)
	{
		PRISMATIC_FLOAT_PRECISION yB = pars.xyBeams.at(a4, 0);
		PRISMATIC_FLOAT_PRECISION xB = pars.xyBeams.at(a4, 1);

		if (abs(pars.psiProbeInit.at(yB, xB)) > 0)
		{
			PRISMATIC_FLOAT_PRECISION q0_0 = pars.qxaReduce.at(yB, xB);
			PRISMATIC_FLOAT_PRECISION q0_1 = pars.qyaReduce.at(yB, xB);
			std::complex<PRISMATIC_FLOAT_PRECISION> phaseShift = exp(
				-2 * pi * i * (q0_0 * (pars.xp[ax] + pars.xTiltShift) + q0_1 * (pars.yp[ay] + pars.yTiltShift)));

			const std::complex<PRISMATIC_FLOAT_PRECISION> tmp_const = pars.psiProbeInit.at(yB, xB) * phaseShift;
			for (auto j = 0; j < y.size(); ++j)
			{
				for (auto i = 0; i < x.size(); ++i)
				{
					psi.at(j, i) += (tmp_const * pars.Scompact.at(a4, y[j], x[i]));
				}
			}
		}
	}

	PRISMATIC_FFTW_EXECUTE(plan);

	for (auto jj = 0; jj < intOutput.get_dimj(); ++jj)
	{
		for (auto ii = 0; ii < intOutput.get_dimi(); ++ii)
		{
			intOutput.at(jj, ii) += pow(abs(psi.at(jj, ii)), 2) * pars.scale;
		}
	}

	size_t write_ay = (pars.meta.arbitraryProbes) ? 0 : ay;
	if (pars.meta.saveDPC_CoM)
	{
		//calculate center of mass; qxa, qya are the fourier coordinates, should have 0 components at boundaries
		for (long y = 0; y < intOutput.get_dimj(); ++y)
		{
			for (long x = 0; x < intOutput.get_dimi(); ++x)
			{
				pars.DPC_CoM.at(0, write_ay, ax, 0) += pars.qxaReduce.at(y, x) * intOutput.at(y, x);
				pars.DPC_CoM.at(0, write_ay, ax, 1) += pars.qyaReduce.at(y, x) * intOutput.at(y, x);
			}
		}
		//divide by sum of intensity
		PRISMATIC_FLOAT_PRECISION intensitySum = 0;
		for (auto iter = intOutput.begin(); iter != intOutput.end(); ++iter)
		{
			intensitySum += *iter;
		}
		pars.DPC_CoM.at(0, write_ay, ax, 0) /= intensitySum;
		pars.DPC_CoM.at(0, write_ay, ax, 1) /= intensitySum;
	}

	//         update output -- ax,ay are unique per thread so this write is thread-safe without a lock
	auto idx = pars.alphaInd.begin();
	for (auto counts = intOutput.begin(); counts != intOutput.end(); ++counts)
	{
		if (*idx <= pars.Ndet)
		{
			pars.output.at(0, write_ay, ax, (*idx) - 1) += *counts;
		}
		++idx;
	};

	//save 4D output if applicable
	if (pars.meta.save4DOutput)
	{
		std::string nameString ="4DSTEM_simulation/data/datacubes/CBED_array_depth" + getDigitString(0)+pars.currentTag;


		PRISMATIC_FLOAT_PRECISION numFP = pars.meta.numFP;
		hsize_t offset[4] = {ax, write_ay, 0, 0}; //order by ax, ay so that aligns with py4DSTEM
		if(pars.meta.saveComplexOutputWave)
		{
			nameString += "_fp" + getDigitString(pars.meta.fpNum);
			Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> finalOutput;
			if(pars.meta.crop4DOutput)
			{
				finalOutput = cropOutput(psi, pars);
				finalOutput *= sqrt(pars.scale);
				hsize_t mdims[4] = {1, 1, finalOutput.get_dimj(), finalOutput.get_dimi()};
				writeDatacube4D(pars, &finalOutput[0], &pars.cbed_buffer_c[0], mdims, offset, numFP, nameString.c_str());
			}
			else
			{
				finalOutput = fftshift2_flip(psi);
				finalOutput *= sqrt(pars.scale);
				hsize_t mdims[4] = {1, 1, finalOutput.get_dimj(), finalOutput.get_dimi()};
				writeDatacube4D(pars, &finalOutput[0], &pars.cbed_buffer_c[0], mdims, offset, numFP, nameString.c_str());
			}
		}
		else
		{
			if(pars.meta.crop4DOutput)
			{
				Array2D<PRISMATIC_FLOAT_PRECISION> croppedOutput = cropOutput(intOutput, pars);
				hsize_t mdims[4] = {1, 1, croppedOutput.get_dimj(), croppedOutput.get_dimi()};
				writeDatacube4D(pars, &croppedOutput[0],  &pars.cbed_buffer[0], mdims, offset, numFP, nameString.c_str());
			}
			else
			{
				intOutput = fftshift2_flip(intOutput);
				hsize_t mdims[4] = {1, 1, intOutput.get_dimj(), intOutput.get_dimi()};
				writeDatacube4D(pars, &intOutput[0],  &pars.cbed_buffer[0], mdims, offset, numFP, nameString.c_str());
			}
		}

	}
}

void transformIndices(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	// setup some relevant coordinates

	pars.dq = (pars.qxaReduce.at(0, 1) + pars.qyaReduce.at(1, 0)) / 2;
	PRISMATIC_FLOAT_PRECISION scale = pow(pars.meta.interpolationFactorX, 2) * pow(pars.meta.interpolationFactorY, 2);

	pars.scale = scale;

	//		 The operators +, -, /, * return PRISM arrays by value, so to avoid unnecessary memory
	//		 allocations/copies for chained operations I try to do things like create variables
	//		 initially with at most one operation, and then perform in-place transforms if more is needed
	Array2D<PRISMATIC_FLOAT_PRECISION> qxaShift = pars.qxaReduce - (pars.meta.probeXtilt / pars.lambda);
	Array2D<PRISMATIC_FLOAT_PRECISION> qyaShift = pars.qyaReduce - (pars.meta.probeYtilt / pars.lambda);
	transform(qxaShift.begin(), qxaShift.end(),
			  qyaShift.begin(), pars.q2.begin(),
			  [](const PRISMATIC_FLOAT_PRECISION &a, const PRISMATIC_FLOAT_PRECISION &b) { return a * a + b * b; });
	transform(pars.q2.begin(), pars.q2.end(),
			  pars.q1.begin(),
			  [](const PRISMATIC_FLOAT_PRECISION &a) { return sqrt(a); });
	//						Array2D<PRISMATIC_FLOAT_PRECISION> alphaInd(pars.q1); // copy constructor more efficient than assignment

	
	pars.qTheta = pars.q1;
	std::transform(qxaShift.begin(), qxaShift.end(),
				   qyaShift.begin(), pars.qTheta.begin(), [](const PRISMATIC_FLOAT_PRECISION&a, const PRISMATIC_FLOAT_PRECISION& b){
						return atan2(b,a);
					});

	pars.alphaInd = pars.q1;
	transform(pars.alphaInd.begin(), pars.alphaInd.end(),
			  pars.alphaInd.begin(),
			  // pars.detectorAngles is in mrad
			  [&pars](const PRISMATIC_FLOAT_PRECISION &a) {
				  return 1 + round((a * pars.lambda - pars.detectorAngles[0] / 1000) / pars.meta.detectorAngleStep);
			  });
	transform(pars.alphaInd.begin(), pars.alphaInd.end(),
			  pars.alphaInd.begin(),
			  [](const PRISMATIC_FLOAT_PRECISION &a) { return a < 1 ? 1 : a; });
	Prismatic::ArrayND<2, std::vector<unsigned short>> alphaMask(
		std::vector<unsigned short>(pars.alphaInd.size(), 0),
		{{pars.alphaInd.get_dimj(), pars.alphaInd.get_dimi()}});
	transform(pars.alphaInd.begin(), pars.alphaInd.end(),
			  alphaMask.begin(),
			  [&pars](const PRISMATIC_FLOAT_PRECISION &a) { return (a < pars.Ndet) ? 1 : 0; });
}

void initializeProbes(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	// initialize the probe

	pars.psiProbeInit = zeros_ND<2, complex<PRISMATIC_FLOAT_PRECISION>>(
		{{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
	PRISMATIC_FLOAT_PRECISION qProbeMax = pars.meta.probeSemiangle / pars.lambda;

	// erf probe is deprecated, but keeping the source here in case we ever want to flexibly switch
	// transform(pars.psiProbeInit.begin(), pars.psiProbeInit.end(),
	// 		  pars.q1.begin(), pars.psiProbeInit.begin(),
	// 		  [&pars, &qProbeMax](std::complex<PRISMATIC_FLOAT_PRECISION> &a, PRISMATIC_FLOAT_PRECISION &q1_t) {
	// 			  a.real(erf((qProbeMax - q1_t) / (0.5 * pars.dq)) * 0.5 + 0.5);
	// 			  a.imag(0);
	// 			  return a;
	// 		  });

	PRISMATIC_FLOAT_PRECISION dqx = pars.qxaReduce.at(0,1);
	PRISMATIC_FLOAT_PRECISION dqy = pars.qyaReduce.at(1,0);
	for(auto j = 0; j < pars.q1.get_dimj(); j++)
	{
		for(auto i = 0; i < pars.q1.get_dimi(); i++)
		{
			PRISMATIC_FLOAT_PRECISION tmp_val = (qProbeMax*pars.q1.at(j,i) - pars.q2.at(j,i));
			tmp_val /= sqrt(dqx*dqx*pow(pars.qxaReduce.at(j,i),2.0)+dqy*dqy*pow(pars.qyaReduce.at(j,i),2.0));					
			tmp_val += 0.5; 
			tmp_val = std::max(tmp_val, (PRISMATIC_FLOAT_PRECISION) 0.0);
			tmp_val = std::min(tmp_val, (PRISMATIC_FLOAT_PRECISION) 1.0);
			pars.psiProbeInit.at(j,i).real(tmp_val);
		}
	}

	pars.psiProbeInit.at(0,0).real(1.0);


	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> chi = getChi(pars.q1, pars.qTheta, pars.lambda, pars.meta.aberrations);
	transform(pars.psiProbeInit.begin(), pars.psiProbeInit.end(),
				chi.begin(), pars.psiProbeInit.begin(),
				[](std::complex<PRISMATIC_FLOAT_PRECISION> &a, std::complex<PRISMATIC_FLOAT_PRECISION> &b) {
					a = a * exp(-i * b);
					return a;
				});

	PRISMATIC_FLOAT_PRECISION norm_constant = sqrt(accumulate(pars.psiProbeInit.begin(), pars.psiProbeInit.end(),
															  (PRISMATIC_FLOAT_PRECISION)0.0,
															  [](PRISMATIC_FLOAT_PRECISION accum,
																 std::complex<PRISMATIC_FLOAT_PRECISION> &a) {
																  return accum + abs(a) * abs(a);
															  })); // make sure to initialize with 0.0 and NOT 0 or it won't be a float and answer will be wrong
	PRISMATIC_FLOAT_PRECISION a = 0;
	for (auto &i : pars.psiProbeInit)
	{
		a += i.real();
	};
	transform(pars.psiProbeInit.begin(), pars.psiProbeInit.end(),
			  pars.psiProbeInit.begin(), [&norm_constant](std::complex<PRISMATIC_FLOAT_PRECISION> &a) {
				  return a / norm_constant;
			  });

	if(pars.meta.saveProbe && pars.fpFlag == 0)
	{
		setupProbeOutput(pars);
        saveProbe(pars);
	}
}

void PRISM03_calcOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	// compute final image

	cout << "Entering PRISM03_calcOutput" << endl;
	// setup necessary coordinates
	setupCoordinates_2(pars);
	
	// setup angles of detector and image sizes
	setupDetector(pars);

	// setup coordinates and indices for the beams
	setupBeams_2(pars);
	
	// setup Fourier coordinates for the S-matrix
	setupFourierCoordinates(pars);

	// initialize the output to the correct size for the output mode
	createStack_integrate(pars);

	// perform some necessary setup transformations of the data
	transformIndices(pars);

	// initialize/compute the probes
	initializeProbes(pars);

#ifdef PRISMATIC_BUILDING_GUI
	pars.progressbar->signalDescriptionMessage("Computing final output (PRISM)");
	pars.progressbar->signalOutputUpdate(0, pars.numProbes);
#endif

	// compute the final PRISM output
	buildPRISMOutput(pars);
}
} // namespace Prismatic
