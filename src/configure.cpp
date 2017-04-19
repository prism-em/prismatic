// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:


#include "configure.h"
#include "PRISM_entry.h"
#include "Multislice_entry.h"
#include "Multislice.h"
#include <iostream>
#include "PRISM01.h"
#include "PRISM02.h"
#include "PRISM03.h"
//#define PRISM_ENABLE_GPU



#ifdef PRISM_ENABLE_GPU
#include "Multislice.cuh"
#include "PRISM02.cuh"
#include "PRISM03.cuh"
#include "utility.cuh"
#include "Multislice_entry.h"
#endif //PRISM_ENABLE_GPU
namespace PRISM {
	entry_func execute_plan;
	ms_output_func buildMultisliceOutput;
	prism_output_func buildPRISMOutput;
	format_output_func formatOutput_CPU;
	fill_Scompact_func fill_Scompact;

//#ifdef PRISM_ENABLE_GPU
    template <class T>
    StreamingMode transferMethodAutoChooser(PRISM::Metadata<T>& meta) {
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
		T f = 4 * meta.interpolationFactor;
		Array1D<size_t> imageSize({{meta.cellDim[1], meta.cellDim[2]}}, {{2}});
		std::transform(imageSize.begin(), imageSize.end(), imageSize.begin(),
		               [&f, &meta](size_t &a) {
			               return (size_t) (f * round(((T) a) / meta.realspace_pixelSize / f));
		               });
		size_t estimatedPotentialSize = (meta.cellDim[0] / meta.sliceThickness) * imageSize[0] * imageSize[1] *
		                                sizeof(std::complex<PRISM_FLOAT_PRECISION>);
		estimatedMaxMemoryUsage = estimatedPotentialSize;
		if (meta.algorithm == PRISM::Algorithm::PRISM) {
			Array1D<PRISM_FLOAT_PRECISION> xv = makeFourierCoords(imageSize[1],
			                                                      (PRISM_FLOAT_PRECISION) 1 / imageSize[1]);
			Array1D<PRISM_FLOAT_PRECISION> yv = makeFourierCoords(imageSize[0],
			                                                      (PRISM_FLOAT_PRECISION) 1 / imageSize[0]);
			pair<Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > mesh_a = meshgrid(yv, xv);

			Array1D<PRISM_FLOAT_PRECISION> qx = makeFourierCoords(imageSize[1], meta.realspace_pixelSize);
			Array1D<PRISM_FLOAT_PRECISION> qy = makeFourierCoords(imageSize[0], meta.realspace_pixelSize);

			pair<Array2D<PRISM_FLOAT_PRECISION>, Array2D<PRISM_FLOAT_PRECISION> > mesh = meshgrid(qy, qx);
			Array2D<PRISM_FLOAT_PRECISION> q2(mesh.first);
			transform(mesh.second.begin(), mesh.second.end(),
			          mesh.first.begin(), q2.begin(),
			          [](const PRISM_FLOAT_PRECISION &a, const PRISM_FLOAT_PRECISION &b) {
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
			PRISM::Array2D<unsigned int> mask;
			mask = zeros_ND<2, unsigned int>({{imageSize[0], imageSize[1]}});
			size_t numberBeams = 0;
			long interp_f = (long) meta.interpolationFactor;
			for (auto y = 0; y < qMask.get_dimj(); ++y) {
				for (auto x = 0; x < qMask.get_dimi(); ++x) {
					if (q2.at(y, x) < pow(meta.alphaBeamMax / lambda, 2) &&
					    qMask.at(y, x) == 1 &&
					    (long) round(mesh_a.first.at(y, x)) % interp_f == 0 &&
					    (long) round(mesh_a.second.at(y, x)) % interp_f == 0) {
						mask.at(y, x) = 1;
						++numberBeams;
					}
				}
			}

			size_t estimatedSMatrixSize =
					numberBeams * imageSize[0] * imageSize[1] / 4 * sizeof(std::complex<PRISM_FLOAT_PRECISION>);
			estimatedMaxMemoryUsage = std::max(estimatedSMatrixSize, estimatedPotentialSize);
		}

		size_t available_memory;
		cudaErrchk(cudaMemGetInfo(NULL, &available_memory));
		cout << "available_memory = " << available_memory << endl;
		return (estimatedMaxMemoryUsage > memoryThreshholdFraction * available_memory) ?
		       PRISM::StreamingMode::Stream : PRISM::StreamingMode::SingleXfer;
//        return PRISM::StreamingMode::SingleXfer;
    }
	format_output_func_GPU formatOutput_GPU;

//#endif
	void configure(Metadata<PRISM_FLOAT_PRECISION>& meta) {
		formatOutput_CPU = formatOutput_CPU_integrate;
#ifdef PRISM_ENABLE_GPU
		formatOutput_GPU = formatOutput_GPU_integrate;
#endif
		if (meta.algorithm == Algorithm::PRISM) {
			std::cout << "Execution plan: PRISM w/ single FP configuration" << std::endl;
			execute_plan = PRISM_entry;
#ifdef PRISM_ENABLE_GPU
            if (meta.transfer_mode == PRISM::StreamingMode::Auto){
            	meta.transfer_mode = transferMethodAutoChooser(meta);
            }
			if (meta.transfer_mode == PRISM::StreamingMode::Stream) {
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
#endif //PRISM_ENABLE_GPU
		} else if (meta.algorithm == Algorithm::Multislice) {
			std::cout << "Execution plan: Multislice w/ single FP configuration" << std::endl;
			execute_plan = Multislice_entry;
#ifdef PRISM_ENABLE_GPU
			std::cout << "Using GPU codes" << std::endl;
			if (meta.transfer_mode == PRISM::StreamingMode::Auto){
				meta.transfer_mode = transferMethodAutoChooser(meta);
			}
			if (meta.transfer_mode == PRISM::StreamingMode::Stream) {
				cout << "Using streaming method\n";
				buildMultisliceOutput = buildMultisliceOutput_GPU_streaming;
			} else {
				cout << "Using single transfer method\n";
				buildMultisliceOutput = buildMultisliceOutput_GPU_singlexfer;
			}

#else
			buildMultisliceOutput = buildMultisliceOutput_CPUOnly;
#endif //PRISM_ENABLE_GPU
		}
	}
}
