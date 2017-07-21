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

#include "configure.h"
#include "PRISM_entry.h"
#include "Multislice_entry.h"
#include "Multislice_calcOutput.h"
#include <iostream>
#include "PRISM01_calcPotential.h"
#include "PRISM02_calcSMatrix.h"
#include "PRISM03_calcOutput.h"
//#define PRISMATIC_ENABLE_GPU



#ifdef PRISMATIC_ENABLE_GPU
#include "Multislice_calcOutput.cuh"
#include "PRISM02_calcSMatrix.cuh"
#include "PRISM03_calcOutput.cuh"
#include "utility.cuh"
#include "Multislice_entry.h"
#endif //PRISMATIC_ENABLE_GPU
namespace Prismatic {
	entry_func execute_plan;
	ms_output_func buildMultisliceOutput;
	prism_output_func buildPRISMOutput;
	format_output_func formatOutput_CPU;
	fill_Scompact_func fill_Scompact;

#ifdef PRISMATIC_ENABLE_GPU
    template <class T>
    StreamingMode transferMethodAutoChooser(Prismatic::Metadata<T>& meta) {
		// query all devices and choose based on the minimum compute capability
		// get the total memory on device
		// based on whether the algo is prism or multislice, estimate the largest array and trigger streaming
		// safely far from this limit
		// For PRISM: limit should be the maximum size of the potential or Scompact
		// which is imageSize/2 x imageSize/2 * numberBeams of type complex float/double
		// numberBeams can be calculated again here... not a big deal
		// For Multislice: limit is the potential array which is imageSize x imageSize x numPlanes of type float/double
		constexpr double memoryThreshholdFraction = 0.5; // if estimated size of largest array is greater than this fraction
		// times the amount of available memory on a GPU, streaming mode will be triggered
		size_t estimatedMaxMemoryUsage;

	    T f_x = 4 * meta.interpolationFactorX;
	    T f_y = 4 * meta.interpolationFactorY;
	    Array1D<size_t> imageSize({{(size_t)(meta.cellDim[1] * meta.tileY), (size_t)(meta.cellDim[2] * meta.tileX)}}, {{2}});
	    imageSize[0] = (size_t)std::max((PRISMATIC_FLOAT_PRECISION)4.0,  (f_y * round(((T)imageSize[0]) / meta.realspacePixelSize[0] / f_y)));
	    imageSize[1] = (size_t)std::max((PRISMATIC_FLOAT_PRECISION)4.0,  (f_x * round(((T)imageSize[1]) / meta.realspacePixelSize[1] / f_x)));



		size_t estimatedPotentialSize = (meta.cellDim[0] * meta.tileZ / meta.sliceThickness) * imageSize[0] * imageSize[1] *
		                                sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>);


		// Estimate the amount of memory needed for the various buffers. This is affected by the batch size, which is inputted
		// by the user but will be adjusted if it is inappropriate (i.e. not enough work per thread).
		// Figure out the scan configuration to determine how many probes there are to compute
		Array1D<PRISMATIC_FLOAT_PRECISION> xR = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{2}});
		xR[0] = meta.scanWindowXMin * meta.cellDim[2] * meta.tileX;
		xR[1] = meta.scanWindowXMax * meta.cellDim[2] * meta.tileX;
		Array1D<PRISMATIC_FLOAT_PRECISION> yR = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{2}});
		yR[0] = meta.scanWindowYMin * meta.cellDim[1] * meta.tileY;
		yR[1] = meta.scanWindowYMax * meta.cellDim[1] * meta.tileY;
		vector<PRISMATIC_FLOAT_PRECISION> xp_d = vecFromRange(xR[0], meta.probeStepX, xR[1]);
		vector<PRISMATIC_FLOAT_PRECISION> yp_d = vecFromRange(yR[0], meta.probeStepY, yR[1]);

		// determine the batch size
		size_t batch_size = std::min(meta.batchSizeTargetGPU, max((size_t)1, xp_d.size()*yp_d.size()/ max((size_t)1,(meta.numStreamsPerGPU*meta.numGPUs)))); // make sure the batch is small enough to spread work to all threads

		// estimate the amount of buffer memory needed. The factor of 3 is because there are two arrays that must be allocated space that scales
		// with the batch size, and the cuFFT plans also allocate internal buffers.
		size_t estimatedBatchBufferSize =  meta.numStreamsPerGPU*3*batch_size*imageSize[0]*imageSize[1]*sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>);
		estimatedMaxMemoryUsage = 3*estimatedPotentialSize + estimatedBatchBufferSize; // factor of 3 is because there is a complex array of the same size created

		cout << "Estimated potential array size = " << estimatedPotentialSize << '\n';
		cout << "Estimated buffer memory needed = " << estimatedBatchBufferSize << '\n';
		cout << "meta.numStreamsPerGPU*2*batch_size*imageSize[0]*imageSize[1]= " << meta.numStreamsPerGPU*2*batch_size*imageSize[0]*imageSize[1] << '\n';

		if (meta.algorithm == Prismatic::Algorithm::PRISM) {
			Array1D<PRISMATIC_FLOAT_PRECISION> xv = makeFourierCoords(imageSize[1],
			                                                      (PRISMATIC_FLOAT_PRECISION) 1 / imageSize[1]);
			Array1D<PRISMATIC_FLOAT_PRECISION> yv = makeFourierCoords(imageSize[0],
			                                                      (PRISMATIC_FLOAT_PRECISION) 1 / imageSize[0]);
			pair<Array2D<PRISMATIC_FLOAT_PRECISION>, Array2D<PRISMATIC_FLOAT_PRECISION> > mesh_a = meshgrid(yv, xv);

			Array1D<PRISMATIC_FLOAT_PRECISION> qx = makeFourierCoords(imageSize[1], meta.realspacePixelSize[1]);
			Array1D<PRISMATIC_FLOAT_PRECISION> qy = makeFourierCoords(imageSize[0], meta.realspacePixelSize[0]);

			pair<Array2D<PRISMATIC_FLOAT_PRECISION>, Array2D<PRISMATIC_FLOAT_PRECISION> > mesh = meshgrid(qy, qx);
			Array2D<PRISMATIC_FLOAT_PRECISION> q2(mesh.first);
			transform(mesh.second.begin(), mesh.second.end(),
			          mesh.first.begin(), q2.begin(),
			          [](const PRISMATIC_FLOAT_PRECISION &a, const PRISMATIC_FLOAT_PRECISION &b) {
				          return a * a + b * b;
			          });


			Array2D<unsigned int> qMask = zeros_ND<2, unsigned int>({{imageSize[0], imageSize[1]}});
			{
				long offset_x = qMask.get_dimi() / 4;
				long offset_y = qMask.get_dimj() / 4;
				long ndimy = (long) qMask.get_dimj();
				long ndimx = (long) qMask.get_dimi();
				for (long y = 0; y < qMask.get_dimj() / 2; ++y) {
					for (long x = 0; x < qMask.get_dimi() / 2; ++x) {
						qMask.at(((y - offset_y) % ndimy + ndimy) % ndimy,
						         ((x - offset_x) % ndimx + ndimx) % ndimx) = 1;
					}
				}
			}

			constexpr double m = 9.109383e-31;
			constexpr double e = 1.602177e-19;
			constexpr double c = 299792458;
			constexpr double h = 6.62607e-34;
			T lambda = (T) (h / sqrt(2 * m * e * meta.E0) / sqrt(1 + e * meta.E0 / 2 / m / c / c) * 1e10);

			// create beam mask and count beams
			Prismatic::Array2D<unsigned int> mask;
			mask = zeros_ND<2, unsigned int>({{imageSize[0], imageSize[1]}});
			size_t numberBeams = 0;
			long interp_fx = (long) meta.interpolationFactorX;
			long interp_fy = (long) meta.interpolationFactorY;
			for (auto y = 0; y < qMask.get_dimj(); ++y) {
				for (auto x = 0; x < qMask.get_dimi(); ++x) {
					if (q2.at(y, x) < pow(meta.alphaBeamMax / lambda, 2) &&
					    qMask.at(y, x) == 1 &&
					    (long) round(mesh_a.first.at(y, x)) % interp_fy == 0 &&
					    (long) round(mesh_a.second.at(y, x)) % interp_fx == 0) {
						mask.at(y, x) = 1;
						++numberBeams;
					}
				}
			}

			size_t estimatedSMatrixSize =
			numberBeams * imageSize[0] * imageSize[1] / 4 * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>);

			estimatedMaxMemoryUsage = std::max(estimatedSMatrixSize, estimatedMaxMemoryUsage);
		}

#ifdef PRISMATIC_ENABLE_GPU

		size_t available_memory;
		cudaErrchk(cudaMemGetInfo(&available_memory, NULL));
		cout << "Available GPU memory = " << available_memory << '\n';
		cout << "Estimated GPU memory usage for single transfer method = " << estimatedMaxMemoryUsage  << '\n';
		return (estimatedMaxMemoryUsage > memoryThreshholdFraction * available_memory) ?
		       Prismatic::StreamingMode::Stream : Prismatic::StreamingMode::SingleXfer;
#else
        return Prismatic::StreamingMode::SingleXfer;
#endif //PRISMATIC_ENABLE_GPU
    }
	format_output_func_GPU formatOutput_GPU;

#endif
	void configure(Metadata<PRISMATIC_FLOAT_PRECISION>& meta) {
		std::cout << "Formatting" << std::endl;
		formatOutput_CPU = formatOutput_CPU_integrate;
#ifdef PRISMATIC_ENABLE_GPU
		formatOutput_GPU = formatOutput_GPU_integrate;
#endif
		if (meta.algorithm == Algorithm::PRISM) {
			std::cout << "Execution plan: PRISM\n";
			execute_plan = PRISM_entry;
#ifdef PRISMATIC_ENABLE_GPU
            if (meta.transferMode == Prismatic::StreamingMode::Auto){
            	meta.transferMode = transferMethodAutoChooser(meta);
            }
			std::cout << "Using GPU codes" << '\n';
			if (meta.transferMode == Prismatic::StreamingMode::Stream) {
				cout << "Using streaming method\n";
				fill_Scompact = fill_Scompact_GPU_streaming;
				buildPRISMOutput = buildPRISMOutput_GPU_streaming;
			} else {
				cout << "Using single transfer method\n";
				fill_Scompact = fill_Scompact_GPU_singlexfer;
				buildPRISMOutput = buildPRISMOutput_GPU_singlexfer;
			}

#else
			fill_Scompact = fill_Scompact_CPUOnly;
			buildPRISMOutput = buildPRISMOutput_CPUOnly;
#endif //PRISMATIC_ENABLE_GPU
		} else if (meta.algorithm == Algorithm::Multislice) {
			std::cout << "Execution plan: Multislice\n";
			execute_plan = Multislice_entry;
#ifdef PRISMATIC_ENABLE_GPU
			std::cout << "Using GPU codes" << '\n';
			if (meta.transferMode == Prismatic::StreamingMode::Auto){
				meta.transferMode = transferMethodAutoChooser(meta);
			}
			if (meta.transferMode == Prismatic::StreamingMode::Stream) {
				cout << "Using streaming method\n";
				buildMultisliceOutput = buildMultisliceOutput_GPU_streaming;
			} else {
				cout << "Using single transfer method\n";
				buildMultisliceOutput = buildMultisliceOutput_GPU_singlexfer;
			}

#else
			buildMultisliceOutput = buildMultisliceOutput_CPUOnly;
#endif //PRISMATIC_ENABLE_GPU
		}
	}
}
