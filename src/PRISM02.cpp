// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include "PRISM02.h"
#include "params.h"
#include <iostream>
#include <vector>
#include <thread>
#include "fftw3.h"
#include <mutex>
#include "ArrayND.h"
#include <complex>
#include "utility.h"
#include "configure.h"
#include "WorkDispatcher.h"
#ifdef PRISM_BUILDING_GUI
#include "prism_progressbar.h"
#endif

namespace PRISM {

	using namespace std;
	const PRISM_FLOAT_PRECISION pi = acos(-1);
	const std::complex<PRISM_FLOAT_PRECISION> i(0, 1);

	void setupCoordinates(Parameters<PRISM_FLOAT_PRECISION> &pars) {

		// setup some Fourier coordinates and propagators
		pars.imageSize[0] = pars.pot.get_dimj();
		pars.imageSize[1] = pars.pot.get_dimi();
		Array1D<PRISM_FLOAT_PRECISION> qx = makeFourierCoords(pars.imageSize[1], pars.pixelSize[1]);
		Array1D<PRISM_FLOAT_PRECISION> qy = makeFourierCoords(pars.imageSize[0], pars.pixelSize[0]);

		pair<Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > mesh = meshgrid(qy, qx);
		pars.qya = mesh.first;
		pars.qxa = mesh.second;
		Array2D<PRISM_FLOAT_PRECISION> q2(pars.qya);
		transform(pars.qxa.begin(), pars.qxa.end(),
		          pars.qya.begin(), q2.begin(), [](const PRISM_FLOAT_PRECISION &a, const PRISM_FLOAT_PRECISION &b) {
					return a * a + b * b;
				});
		pars.q2 = q2;

		// get qMax
		pars.qMax = 0;
		{
			PRISM_FLOAT_PRECISION qx_max = 0;
			PRISM_FLOAT_PRECISION qy_max = 0;
			for (auto i = 0; i < qx.size(); ++i) {
				qx_max = (abs(qx[i]) > qx_max) ? abs(qx[i]) : qx_max;
				qy_max = (abs(qy[i]) > qy_max) ? abs(qy[i]) : qy_max;
			}
			pars.qMax = min(qx_max, qy_max) / 2;
		}

		pars.qMask = zeros_ND<2, unsigned int>({{pars.imageSize[0], pars.imageSize[1]}});
		{
			long offset_x = pars.qMask.get_dimi() / 4;
			long offset_y = pars.qMask.get_dimj() / 4;
			long ndimy = (long) pars.qMask.get_dimj();
			long ndimx = (long) pars.qMask.get_dimi();
			for (long y = 0; y < pars.qMask.get_dimj() / 2; ++y) {
				for (long x = 0; x < pars.qMask.get_dimi() / 2; ++x) {
					pars.qMask.at(((y - offset_y) % ndimy + ndimy) % ndimy,
					              ((x - offset_x) % ndimx + ndimx) % ndimx) = 1;
				}
			}
		}

		// build propagators
		pars.prop = zeros_ND<2, std::complex<PRISM_FLOAT_PRECISION> >({{pars.imageSize[0], pars.imageSize[1]}});
		pars.propBack = zeros_ND<2, std::complex<PRISM_FLOAT_PRECISION> >({{pars.imageSize[0], pars.imageSize[1]}});
		for (auto y = 0; y < pars.qMask.get_dimj(); ++y) {
			for (auto x = 0; x < pars.qMask.get_dimi(); ++x) {
				if (pars.qMask.at(y, x) == 1) {
					pars.prop.at(y, x) = exp(-i * pi * complex<PRISM_FLOAT_PRECISION>(pars.lambda, 0) *
					                         complex<PRISM_FLOAT_PRECISION>(pars.meta.sliceThickness, 0) *
					                         complex<PRISM_FLOAT_PRECISION>(pars.q2.at(y, x), 0));
					pars.propBack.at(y, x) = exp(i * pi * complex<PRISM_FLOAT_PRECISION>(pars.lambda, 0) *
					                             complex<PRISM_FLOAT_PRECISION>(pars.meta.cellDim[0], 0) *
					                             complex<PRISM_FLOAT_PRECISION>(pars.q2.at(y, x), 0));
				}
			}
		}
	}

	inline void setupBeams(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		Array1D<PRISM_FLOAT_PRECISION> xv = makeFourierCoords(pars.imageSize[1],
		                                                      (PRISM_FLOAT_PRECISION) 1 / pars.imageSize[1]);
		Array1D<PRISM_FLOAT_PRECISION> yv = makeFourierCoords(pars.imageSize[0],
		                                                      (PRISM_FLOAT_PRECISION) 1 / pars.imageSize[0]);
		pair<Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > mesh_a = meshgrid(yv, xv);

		// create beam mask and count beams
		PRISM::Array2D<unsigned int> mask;
		mask = zeros_ND<2, unsigned int>({{pars.imageSize[0], pars.imageSize[1]}});
		pars.numberBeams = 0;
		long interp_fx = (long) pars.meta.interpolationFactorX;
		long interp_fy = (long) pars.meta.interpolationFactorY;
		for (auto y = 0; y < pars.qMask.get_dimj(); ++y) {
			for (auto x = 0; x < pars.qMask.get_dimi(); ++x) {
				if (pars.q2.at(y, x) < pow(pars.meta.alphaBeamMax / pars.lambda, 2) &&
				    pars.qMask.at(y, x) == 1 &&
				    (long) round(mesh_a.first.at(y, x)) % interp_fy == 0 &&
				    (long) round(mesh_a.second.at(y, x)) % interp_fx == 0) {
					mask.at(y, x) = 1;
					++pars.numberBeams;
				}
			}
		}
		pars.beams = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.imageSize[0], pars.imageSize[1]}});
		{
			int beam_count = 1;
			for (auto y = 0; y < pars.qMask.get_dimj(); ++y) {
				for (auto x = 0; x < pars.qMask.get_dimi(); ++x) {
					if (mask.at(y, x) == 1) {
						pars.beamsIndex.push_back((size_t) y * pars.qMask.get_dimi() + (size_t) x);
						pars.beams.at(y, x) = beam_count++;
					}
				}
			}
		}
	}

	inline void setupSMatrixCoordinates(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		// TODO: ensure this block is correct for arbitrary dimension
		// get the indices for the compact S-matrix
		pars.qxInd = zeros_ND<1, size_t>({{pars.imageSize[1] / 2}});
		pars.qyInd = zeros_ND<1, size_t>({{pars.imageSize[0] / 2}});
		{
			long n_0 = pars.imageSize[0];
			long n_1 = pars.imageSize[1];
			long n_half0 = pars.imageSize[0] / 2;
			long n_half1 = pars.imageSize[1] / 2;
			long n_quarter0 = pars.imageSize[0] / 4;
			long n_quarter1 = pars.imageSize[1] / 4;
			for (auto i = 0; i < n_quarter0; ++i) {
				pars.qyInd[i] = i;
				pars.qyInd[i + n_quarter0] = (i - n_quarter0) + n_0;
			}
			for (auto i = 0; i < n_quarter1; ++i) {
				pars.qxInd[i] = i;
				pars.qxInd[i + n_quarter1] = (i - n_quarter1) + n_1;
			}
		}
	}

	inline void downsampleFourierComponents(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		// downsample Fourier components by x2 to match output
		pars.imageSizeOutput = pars.imageSize;
		pars.imageSizeOutput[0] /= 2;
		pars.imageSizeOutput[1] /= 2;
		pars.pixelSizeOutput = pars.pixelSize;
		pars.pixelSizeOutput[0] *= 2;
		pars.pixelSizeOutput[1] *= 2;

		pars.qxaOutput = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.qyInd.size(), pars.qxInd.size()}});
		pars.qyaOutput = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.qyInd.size(), pars.qxInd.size()}});
		pars.beamsOutput = zeros_ND<2, PRISM_FLOAT_PRECISION>({{pars.qyInd.size(), pars.qxInd.size()}});
		for (auto y = 0; y < pars.qyInd.size(); ++y) {
			for (auto x = 0; x < pars.qxInd.size(); ++x) {
				pars.qxaOutput.at(y, x) = pars.qxa.at(pars.qyInd[y], pars.qxInd[x]);
				pars.qyaOutput.at(y, x) = pars.qya.at(pars.qyInd[y], pars.qxInd[x]);
				pars.beamsOutput.at(y, x) = pars.beams.at(pars.qyInd[y], pars.qxInd[x]);
			}
		}
	}

	void propagatePlaneWave_CPU(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                            size_t currentBeam,
	                            Array2D<complex<PRISM_FLOAT_PRECISION> > &psi,
	                            const PRISM_FFTW_PLAN &plan_forward,
	                            const PRISM_FFTW_PLAN &plan_inverse,
	                            mutex &fftw_plan_lock) {
		// propagates a single plan wave and fills in the corresponding section of compact S-matrix

		psi[pars.beamsIndex[currentBeam]] = 1;
		const PRISM_FLOAT_PRECISION slice_size= (PRISM_FLOAT_PRECISION) psi.size();


		PRISM_FFTW_EXECUTE(plan_inverse);
		for (auto &i : psi)i /= slice_size; // fftw scales by N, need to correct
		const complex<PRISM_FLOAT_PRECISION> *trans_t = &pars.transmission[0];
		for (auto a2 = 0; a2 < pars.numPlanes; ++a2) {

			for (auto &p:psi)p *= (*trans_t++); // transmit
			PRISM_FFTW_EXECUTE(plan_forward); // FFT
			for (auto i = psi.begin(), j = pars.prop.begin(); i != psi.end(); ++i, ++j)*i *= (*j); // propagate
			PRISM_FFTW_EXECUTE(plan_inverse); // IFFT
			for (auto &i : psi)i /= slice_size; // fftw scales by N, need to correct
		}
		PRISM_FFTW_EXECUTE(plan_forward);

		Array2D<complex<PRISM_FLOAT_PRECISION> > psi_small = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.qyInd.size(), pars.qxInd.size()}});


		unique_lock<mutex> gatekeeper(fftw_plan_lock);
		PRISM_FFTW_PLAN plan_final = PRISM_FFTW_PLAN_DFT_2D(psi_small.get_dimj(), psi_small.get_dimi(),
		                                                    reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_small[0]),
		                                                    reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_small[0]),
		                                                    FFTW_BACKWARD, FFTW_ESTIMATE);
		gatekeeper.unlock();
		for (auto y = 0; y < pars.qyInd.size(); ++y) {
			for (auto x = 0; x < pars.qxInd.size(); ++x) {
				psi_small.at(y, x) = psi.at(pars.qyInd[y], pars.qxInd[x]);
			}
		}
		PRISM_FFTW_EXECUTE(plan_final);
		gatekeeper.lock();
		PRISM_FFTW_DESTROY_PLAN(plan_final);
		gatekeeper.unlock();


		complex<PRISM_FLOAT_PRECISION> *S_t = &pars.Scompact[currentBeam * pars.Scompact.get_dimj() * pars.Scompact.get_dimi()];
		const PRISM_FLOAT_PRECISION N_small = (PRISM_FLOAT_PRECISION) psi_small.size();
		for (auto &jj:psi_small) {
			*S_t++ = jj / N_small;
		}
	}



	void propagatePlaneWave_CPU_batch(Parameters<PRISM_FLOAT_PRECISION> &pars,
	                                  size_t currentBeam,
	                                  size_t stopBeam,
	                                  Array1D<complex<PRISM_FLOAT_PRECISION> > &psi_stack,
	                                  const PRISM_FFTW_PLAN &plan_forward,
	                                  const PRISM_FFTW_PLAN &plan_inverse,
	                                  mutex &fftw_plan_lock) {
		// propagates a batch of plane waves and fills in the corresponding sections of compact S-matrix
		const size_t slice_size                  =  pars.imageSize[0]*pars.imageSize[1];
		const PRISM_FLOAT_PRECISION slice_size_f = (PRISM_FLOAT_PRECISION) slice_size;
		{
			int beam_count = 0;
			for (auto jj = currentBeam; jj < stopBeam; ++jj) {
				psi_stack[beam_count*slice_size + pars.beamsIndex[jj]] = 1;
				++beam_count;
			}
		}


		PRISM_FFTW_EXECUTE(plan_inverse);
		for (auto &i : psi_stack)i /= slice_size_f; // fftw scales by N, need to correct
		complex<PRISM_FLOAT_PRECISION>* slice_ptr = &pars.transmission[0];
		for (auto a2 = 0; a2 < pars.numPlanes; ++a2) {

			// transmit each of the probes in the batch
			for (auto batch_idx = 0; batch_idx < min(pars.meta.batch_size_CPU, stopBeam - currentBeam); ++batch_idx){
				auto t_ptr   = slice_ptr; // start at the beginning of the current slice
				auto psi_ptr = &psi_stack[batch_idx * slice_size];
				for (auto jj = 0; jj < slice_size; ++jj){
					*psi_ptr++ *= (*t_ptr++);// transmit
				}
			}
			slice_ptr += slice_size; // advance to point to the beginning of the next potential slice
//			for (auto &p:psi)p *= (*trans_t++); // transmit
			PRISM_FFTW_EXECUTE(plan_forward); // FFT

			// propagate each of the probes in the batch
			for (auto batch_idx = 0; batch_idx < min(pars.meta.batch_size_CPU, stopBeam - currentBeam); ++batch_idx) {
				auto p_ptr = pars.prop.begin();
				auto psi_ptr = &psi_stack[batch_idx * slice_size];
				for (auto jj = 0; jj < slice_size; ++jj){
					*psi_ptr++ *= (*p_ptr++);// propagate
				}
			}
//			for (auto i = psi.begin(), j = pars.prop.begin(); i != psi.end(); ++i, ++j)*i *= (*j); // propagate
			PRISM_FFTW_EXECUTE(plan_inverse); // IFFT
			for (auto &i : psi_stack)i /= slice_size_f; // fftw scales by N, need to correct
		}
		PRISM_FFTW_EXECUTE(plan_forward);

		Array2D<complex<PRISM_FLOAT_PRECISION> > psi_small = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.qyInd.size(), pars.qxInd.size()}});
		const PRISM_FLOAT_PRECISION N_small = (PRISM_FLOAT_PRECISION) psi_small.size();
		unique_lock<mutex> gatekeeper(fftw_plan_lock);
		PRISM_FFTW_PLAN plan_final = PRISM_FFTW_PLAN_DFT_2D(psi_small.get_dimj(), psi_small.get_dimi(),
		                                                    reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_small[0]),
		                                                    reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_small[0]),
		                                                    FFTW_BACKWARD, FFTW_ESTIMATE);
		gatekeeper.unlock();
		int batch_idx = 0;
		while (currentBeam < stopBeam) {

			// can later remove this and batch this section as well
//			Array2D<complex<PRISM_FLOAT_PRECISION> > psi = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >(
//					{{pars.imageSize[0], pars.imageSize[1]}});
//			auto psi_ptr = &psi_stack[batch_idx * slice_size];
//			for (auto&i:psi)i=*psi_ptr++;
			for (auto y = 0; y < pars.qyInd.size(); ++y) {
				for (auto x = 0; x < pars.qxInd.size(); ++x) {
//					psi_small.at(y, x) = psi.at(pars.qyInd[y], pars.qxInd[x]);
					psi_small.at(y, x) = psi_stack[batch_idx*slice_size + pars.qyInd[y]*pars.imageSize[1] +  pars.qxInd[x]];
				}
			}
			PRISM_FFTW_EXECUTE(plan_final);
			complex<PRISM_FLOAT_PRECISION> *S_t = &pars.Scompact[currentBeam * pars.Scompact.get_dimj() * pars.Scompact.get_dimi()];
			for (auto &jj:psi_small) {
				*S_t++ = jj / N_small;
			}
			++currentBeam;
			++batch_idx;
		}
		gatekeeper.lock();
		PRISM_FFTW_DESTROY_PLAN(plan_final);
		gatekeeper.unlock();
	}


	void fill_Scompact_CPUOnly(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		// populates the compact S-matrix using CPU resources
#ifdef PRISM_BUILDING_GUI
        pars.progressbar->signalDescriptionMessage("Computing compact S-matrix");
		pars.progressbar->signalScompactUpdate(-1, pars.numberBeams);
#endif
        extern mutex fftw_plan_lock;
		pars.Scompact = zeros_ND<3, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.numberBeams, pars.imageSize[0] / 2, pars.imageSize[1] / 2}});
		pars.transmission = zeros_ND<3, complex<PRISM_FLOAT_PRECISION> >(
				{{pars.pot.get_dimk(), pars.pot.get_dimj(), pars.pot.get_dimi()}});
		{
			auto p = pars.pot.begin();
			for (auto &j:pars.transmission)j = exp(i * pars.sigma * (*p++));
		}

		vector<thread> workers;
		workers.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
//		setWorkStartStop(0, pars.numberBeams, 1);
		const size_t PRISM_PRINT_FREQUENCY_BEAMS = pars.numberBeams / 10; // for printing status
		WorkDispatcher dispatcher(0, pars.numberBeams);
		pars.meta.batch_size_CPU = min(pars.meta.batch_size_target_CPU, pars.numberBeams / pars.meta.NUM_THREADS);
		cout << "PRISM02 pars.meta.batch_size_CPU = " << pars.meta.batch_size_CPU << endl;

		PRISM_FFTW_INIT_THREADS();
		PRISM_FFTW_PLAN_WITH_NTHREADS(pars.meta.NUM_THREADS);
		for (auto t = 0; t < pars.meta.NUM_THREADS; ++t) {

			cout << "Launching thread #" << t << " to compute beams\n";
            workers.push_back(thread([&pars, &dispatcher, &PRISM_PRINT_FREQUENCY_BEAMS]() {


	            // allocate array for psi just once per thread
//				Array2D<complex<PRISM_FLOAT_PRECISION> > psi = zeros_ND<2, complex<PRISM_FLOAT_PRECISION> >(
//						{{pars.imageSize[0], pars.imageSize[1]}});
	            Array1D<complex<PRISM_FLOAT_PRECISION> > psi_stack = zeros_ND<1, complex<PRISM_FLOAT_PRECISION> >({{pars.imageSize[0]*pars.imageSize[1] * pars.meta.batch_size_CPU}});
//				PRISM_FFTW_PLAN plan_forward = PRISM_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
//				                                                      reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
//				                                                      reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
//				                                                      FFTW_FORWARD, FFTW_MEASURE);
//				PRISM_FFTW_PLAN plan_inverse = PRISM_FFTW_PLAN_DFT_2D(psi.get_dimj(), psi.get_dimi(),
//				                                                      reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
//				                                                      reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi[0]),
//				                                                      FFTW_BACKWARD, FFTW_MEASURE);


	            // setup batch FFTW parameters
	            constexpr int rank = 2;
	            int n[] = {pars.imageSize[0], pars.imageSize[1]};
	            const int howmany = pars.meta.batch_size_CPU;
	            int idist = n[0]*n[1];
	            int odist = n[0]*n[1];
	            int istride = 1;
	            int ostride = 1;
	            int *inembed = n;
	            int *onembed = n;

	            unique_lock<mutex> gatekeeper(fftw_plan_lock);
	            PRISM_FFTW_PLAN plan_forward = PRISM_FFTW_PLAN_DFT_BATCH(rank, n, howmany,
	                                                                     reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), inembed,
	                                                                     istride, idist,
	                                                                     reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), onembed,
	                                                                     ostride, odist,
	                                                                     FFTW_FORWARD, FFTW_MEASURE);
	            PRISM_FFTW_PLAN plan_inverse = PRISM_FFTW_PLAN_DFT_BATCH(rank, n, howmany,
	                                                                     reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), inembed,
	                                                                     istride, idist,
	                                                                     reinterpret_cast<PRISM_FFTW_COMPLEX *>(&psi_stack[0]), onembed,
	                                                                     ostride, odist,
	                                                                     FFTW_BACKWARD, FFTW_MEASURE);
				gatekeeper.unlock(); // unlock it so we only block as long as necessary to deal with plans
				size_t currentBeam, stopBeam;
                currentBeam=stopBeam=0;
//				while (getWorkID(pars, currentBeam, stopBeam)) { // synchronously get work assignment
                while (dispatcher.getWork(currentBeam, stopBeam, pars.meta.batch_size_CPU)) { // synchronously get work assignment
                    while (currentBeam < stopBeam) {
	                    if (currentBeam % PRISM_PRINT_FREQUENCY_BEAMS < pars.meta.batch_size_CPU| currentBeam == 100){
		                    cout << "Computing Plane Wave #" << currentBeam << "/" << pars.numberBeams << endl;
	                    }

						// re-zero psi each iteration
						memset((void *) &psi_stack[0], 0, psi_stack.size() * sizeof(complex<PRISM_FLOAT_PRECISION>));
//						propagatePlaneWave_CPU(pars, currentBeam, psi, plan_forward, plan_inverse, fftw_plan_lock);
	                    propagatePlaneWave_CPU_batch(pars, currentBeam, stopBeam, psi_stack, plan_forward, plan_inverse, fftw_plan_lock);
#ifdef PRISM_BUILDING_GUI
                        pars.progressbar->signalScompactUpdate(currentBeam, pars.numberBeams);
#endif
                        currentBeam=stopBeam;
					}
				}
				// clean up
				gatekeeper.lock();
				PRISM_FFTW_DESTROY_PLAN(plan_forward);
				PRISM_FFTW_DESTROY_PLAN(plan_inverse);
				gatekeeper.unlock();
			}));
		}
		cout << "Waiting for threads...\n";
		for (auto &t:workers)t.join();
		PRISM_FFTW_CLEANUP_THREADS();
#ifdef PRISM_BUILDING_GUI
        pars.progressbar->setProgress(100);
        pars.progressbar->signalCalcStatusMessage(QString("Plane Wave ") +
                                                  QString::number(pars.numberBeams) +
                                                  QString("/") +
                                                  QString::number(pars.numberBeams));
#endif //PRISM_BUILDING_GUI
	}

	void PRISM02(Parameters<PRISM_FLOAT_PRECISION> &pars) {
		// propagate plane waves to construct compact S-matrix

		cout << "Entering PRISM02" << endl;

		// setup some coordinates
		setupCoordinates(pars);

		// setup the beams and their indices
		setupBeams(pars);

		// setup coordinates for nonzero values of compact S-matrix
		setupSMatrixCoordinates(pars);

		cout << "Computing compact S matrix" << endl;

		// populate compact-S matrix
		fill_Scompact(pars);

		downsampleFourierComponents(pars);
	}
}
