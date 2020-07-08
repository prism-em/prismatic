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

#ifndef PRISM_PARAMS_H
#define PRISM_PARAMS_H
#include <vector>
#include <string>
#include <algorithm>
#include <mutex>
#include <complex>
#include "ArrayND.h"
#include "atom.h"
#include "meta.h"
#include "H5Cpp.h"
// #include "utility.h"

#ifdef PRISMATIC_BUILDING_GUI
class prism_progressbar;
#endif

namespace Prismatic{

	template <class T>
	using Array1D = Prismatic::ArrayND<1, std::vector<T> >;
	template <class T>
	using Array2D = Prismatic::ArrayND<2, std::vector<T> >;
	template <class T>
	using Array3D = Prismatic::ArrayND<3, std::vector<T> >;
	template <class T>
	using Array4D = Prismatic::ArrayND<4, std::vector<T> >;

	// for monitoring memory consumption on GPU
	static std::mutex memLock;

	// for threadsafe HDF5 writing to 4D output
	// static std::mutex HDF5_lock;
	
    template <class T>
    class Parameters {

    public:

	    void calculateLambda();
		void calculateFileSize();
	    Metadata<T> meta;
	    Array3D< std::complex<T>  > Scompact;
	    Array4D<T> output;
		Array4D<std::complex<T>> output_c;
	    Array4D<T> net_output;
		Array4D<std::complex<T>> net_output_c;
		Array4D<T> DPC_CoM;
		Array4D<T> net_DPC_CoM;
		Array3D<T> pot;
	    Array3D<std::complex<T> > transmission;
	    Array2D< std::complex<T>  > prop;
	    Array2D< std::complex<T> > propBack;
	    Array2D< std::complex<T> > psiProbeInit;
		Array2D<T> cbed_buffer;
		Array2D<std::complex<T>> cbed_buffer_c;
	    Array2D<unsigned int> qMask;
	    T zTotal;
	    T xTiltShift;
	    T yTiltShift;
		T minXtilt_tem;
		T minYtilt_tem;
		T maxXtilt_tem;
		T maxYtilt_tem;
		T xTiltOffset_tem;
		T yTiltOffset_tem;
		T xTiltStep_tem;
		T yTiltStep_tem;
		std::vector<T> xTilts_tem;
		std::vector<T> yTilts_tem;
		std::vector<int> xTiltsInd_tem;
		std::vector<int> yTiltsInd_tem;
	    Array2D<T> qxa;
	    Array2D<T> qya;
	    Array2D<T> qxaOutput;
	    Array2D<T> qyaOutput;
        Array2D<T> qxaReduce;
        Array2D<T> qyaReduce;
	    Array2D<T> alphaInd;
	    Array2D<T> q2;
	    Array2D<T> q1;
        Array1D<T> xp;
        Array1D<T> yp;
		Array1D<T> qx;
		Array1D<T> qy;
        std::vector<size_t> beamsIndex;
		std::vector<size_t> HRTEMbeamOrder;
	    Prismatic::ArrayND<2, std::vector<long> > xyBeams;
		Array2D<T> beams;
	    Array2D<T> beamsOutput;
        Array1D<T> xVec;
        Array1D<T> yVec;
        Array1D<T> detectorAngles;
	    std::vector<atom> atoms;
	    std::vector<T> pixelSize;
	    std::vector<T> pixelSizeOutput;
		T dzPot;
	    Array1D<size_t> imageSize;
	    std::vector<size_t> imageSizeReduce;
	    Array1D<size_t> imageSizeOutput;
	    Array1D<size_t> qxInd;
	    Array1D<size_t> qyInd;
	    std::vector<T> tiledCellDim;
	    T scale;
        T lambda;
        T dr;
        T dq;
	    T sigma;
	    T qMax;
	    T alphaMax;
		T scanWindowXMin;
		T scanWindowXMax;
		T scanWindowYMin;
		T scanWindowYMax;
        size_t Ndet;
        size_t numPlanes;
		size_t numSlices;
		size_t zStartPlane;
		size_t numLayers;
		size_t numXprobes;
		size_t numYprobes;
		size_t numProbes;
		std::vector<T> depths;
	    size_t numberBeams;
		H5::H5File outputFile;
		size_t fpFlag; //flag to prevent creation of new HDF5 files

		#ifdef PRISMATIC_ENABLE_GPU
				cudaDeviceProp deviceProperties;
		//#ifndef NDEBUG
				// for monitoring memory consumption on GPU
				size_t maxGPUMem;
				size_t targetNumBlocks; // estimate for a good number of blocks to launch on GPU so that enough are made to fill the device without incurring too much overhead unnecessarily
		//#endif //NDEBUG
		#endif // PRISMATIC_ENABLE_GPU
		#ifdef PRISMATIC_BUILDING_GUI
				prism_progressbar *progressbar;
		#endif
		Parameters(){};
		#ifdef PRISMATIC_BUILDING_GUI
		Parameters(Metadata<T> _meta, prism_progressbar* _progressbar = NULL) : meta(_meta), progressbar(_progressbar){
		#else
	    Parameters(Metadata<T> _meta) : meta(_meta){
		#endif

		    constexpr double m = 9.109383e-31;
		    constexpr double e = 1.602177e-19;
		    constexpr double c = 299792458;
		    //constexpr double h = 6.62607e-34;
		    const double pi = std::acos(-1);

			try {
				atoms = tileAtoms(meta.tileX, meta.tileY, meta.tileZ, readAtoms_xyz(meta.filenameAtoms));
				if (!meta.userSpecifiedCelldims){
					std::array<double, 3> dims = peekDims_xyz(meta.filenameAtoms);
					meta.cellDim[0] = dims[0];
					meta.cellDim[1] = dims[1];
					meta.cellDim[2] = dims[2];
				}
			}
			catch (const std::runtime_error &e) {
				std::cout << "Prismatic: Error opening " << meta.filenameAtoms << std::endl;
				std::cout << e.what();
				throw;
			}
			catch (const std::domain_error &e) {
				std::cout << "Prismatic: Error extracting atomic data from " << meta.filenameAtoms << "!" << std::endl;
				std::cout << e.what();
				throw;
			}

		    // tile the cell dimension.
		    tiledCellDim     = meta.cellDim;
		    tiledCellDim[2] *= meta.tileX;
		    tiledCellDim[1] *= meta.tileY;
		    tiledCellDim[0] *= meta.tileZ;
			//	std::cout << "tiledCellDim[0]= " << tiledCellDim[0]<< std::endl;
		    zTotal = tiledCellDim[0];
		    xTiltShift = -zTotal * tan(meta.probeXtilt);
		    yTiltShift = -zTotal * tan(meta.probeYtilt);

			if(meta.realSpaceWindow_x){
				scanWindowXMin = std::min(meta.scanWindowXMin_r, tiledCellDim[2]) / tiledCellDim[2]; //default to max size if dimension exceeds tiled cell
				scanWindowXMax = std::min(meta.scanWindowXMax_r, tiledCellDim[2]) / tiledCellDim[2];
			}else{
				scanWindowXMin = meta.scanWindowXMin;
				scanWindowXMax = meta.scanWindowXMax;
			}

			if(meta.realSpaceWindow_y){
				scanWindowYMin = std::min(meta.scanWindowYMin_r, tiledCellDim[1]) / tiledCellDim[1]; //default to max size if dimension exceeds tiled cell
				scanWindowYMax = std::min(meta.scanWindowYMax_r, tiledCellDim[1]) / tiledCellDim[1];
			}else{
				scanWindowYMin = meta.scanWindowYMin;
				scanWindowYMax = meta.scanWindowYMax;
			}

			calculateLambda();
		    sigma = (T)((2 * pi / lambda / meta.E0) * (m * c * c + e * meta.E0) /
		                       (2 * m * c * c + e * meta.E0));

		    T f_x = 4 * meta.interpolationFactorX;
		    T f_y = 4 * meta.interpolationFactorY;
		    //std::cout << "f_x = " << f_x << std::endl;
		    //std::cout << "f_y = " << f_y << std::endl;
		    //std::cout << "tiledCellDim[1] = " << tiledCellDim[1] << std::endl;
		    //std::cout << "tiledCellDim[2] = " << tiledCellDim[2] << std::endl;
		    Array1D<size_t> _imageSize({{(size_t)tiledCellDim[1], (size_t)tiledCellDim[2]}}, {{2}});
		    _imageSize[0] = (size_t)std::max((PRISMATIC_FLOAT_PRECISION)4.0,  (PRISMATIC_FLOAT_PRECISION)(f_y * round((tiledCellDim[1]) / meta.realspacePixelSize[0] / f_y)));
		    _imageSize[1] = (size_t)std::max((PRISMATIC_FLOAT_PRECISION)4.0,  (PRISMATIC_FLOAT_PRECISION)(f_x * round((tiledCellDim[2]) / meta.realspacePixelSize[1] / f_x)));

		    // std::cout << "(f_y * round((tiledCellDim[1]) / meta.realspacePixelSize[0] / f_y) = " << (f_y * round((tiledCellDim[1]) / meta.realspacePixelSize[0] / f_y)) << std::endl;
		    // std::cout << "_imageSize[0] = " << _imageSize[0] << std::endl;
		    // std::cout << "_imageSize[1] = " << _imageSize[1] << std::endl;
			//		    std::transform(_imageSize.begin(), _imageSize.end(), _imageSize.begin(),
			//		                   [&f, this](size_t &a) {
			//			                   return (size_t)std::max(4.0,  (f * round(((T)a) / meta.realspacePixelSize / f)));
			//		                   });
		    this->imageSize = _imageSize;

		    std::vector<T> _pixelSize{(T) tiledCellDim[1], (T) tiledCellDim[2]};
		    pixelSize = _pixelSize;
		    pixelSize[0] /= (T)imageSize[0];
		    pixelSize[1] /= (T)imageSize[1];

			numSlices = meta.numSlices;
			zStartPlane = (size_t) std::ceil(meta.zStart / meta.sliceThickness);

			//set tilt properties to prevent out of bound access
			if(meta.algorithm == Algorithm::HRTEM)
			{
				T maxXtilt = lambda / (4.0 * pixelSize[1]);
				T maxYtilt = lambda / (4.0 * pixelSize[0]);
				minXtilt_tem = std::min(meta.minXtilt, maxXtilt);
				maxXtilt_tem = std::min(meta.maxXtilt, maxXtilt);
				minYtilt_tem = std::min(meta.minYtilt, maxYtilt);
				maxYtilt_tem = std::min(meta.maxYtilt, maxYtilt);
				xTiltOffset_tem = meta.xTiltOffset;
				yTiltOffset_tem = meta.yTiltOffset;
				xTiltStep_tem = std::min(meta.xTiltStep, maxXtilt_tem);
				yTiltStep_tem = std::min(meta.yTiltStep, maxYtilt_tem);
				

				if(maxXtilt_tem + meta.xTiltOffset > maxXtilt)
				{
					std::cout << "Requested tilt series lies outside of antialiasing aperture in X" << std::endl;
					std::cout << "Resetting X offset to 0.0 mrad" << std::endl;
					xTiltOffset_tem = 0.0;
				}

				if(maxYtilt_tem + meta.yTiltOffset > maxYtilt)
				{
					std::cout << "Requested tilt series lies outside of antialiasing aperture in Y" << std::endl;
					std::cout << "Resetting Y offset to 0.0 mrad" << std::endl;
					yTiltOffset_tem = 0.0;
				}

				std::cout << "tilt settings" << std::endl;
				std::cout << "min tilt step x: " << lambda / (imageSize[1] * pixelSize[1]) << std::endl;
				std::cout << "min tilt step x: " << lambda / (tiledCellDim[2]) << std::endl;
				std::cout << "min tilt step y: " << lambda / (imageSize[0] * pixelSize[0]) << std::endl;
				std::cout << "min tilt step y: " << lambda / (tiledCellDim[1]) << std::endl;
				std::cout << "minXtilt: " << minXtilt_tem << std::endl;
				std::cout << "maxXtilt: " << maxXtilt_tem << std::endl;
				std::cout << "xTiltOffset: " << xTiltOffset_tem << std::endl;
				std::cout << "minYtilt: " << minYtilt_tem << std::endl;
				std::cout << "maxYtilt: " << maxYtilt_tem << std::endl;
				std::cout << "yTiltOffset: " << yTiltOffset_tem << std::endl;
			}

			#ifdef PRISMATIC_ENABLE_GPU
			#ifndef NDEBUG
					// for monitoring memory consumption on GPU
					maxGPUMem = 0;
			#endif //NDEBUG
			// query GPU properties
		    int nDevices;
		    cudaGetDeviceCount(&nDevices);
		    if (nDevices < meta.numGPUs){
			    std::cout << "Warning: User requested " << meta.numGPUs << " GPUs but only " << nDevices << " were found. Proceeding with " << nDevices << ".\n";
			    meta.numGPUs = nDevices;
		    }

		    // Check the properties of each GPU, which is used to choose kernel launch configurations. .
		    // Most machines with multiple GPUs will have identical or similar models, so we keep just one copy of the metadata from the device
		    // with the smallest compute capability
		    cudaErrchk(cudaGetDeviceProperties(&deviceProperties, 0));
		    if (nDevices > 1) {
			    for (auto g = 1; g < nDevices; ++g) {
				    cudaDeviceProp tmp_prop;
				    cudaErrchk(cudaGetDeviceProperties(&tmp_prop, g));
				    if (tmp_prop.major < deviceProperties.major){
					    deviceProperties = tmp_prop;
				    }
			    }
		    }
		    targetNumBlocks = deviceProperties.multiProcessorCount * deviceProperties.maxThreadsPerBlock / BLOCK_SIZE1D *4; // the 4 is a fudge factor
		    std::cout << "deviceProperties.major = " << deviceProperties.major << std::endl;
		    std::cout << "deviceProperties.maxThreadsPerBlock = " << deviceProperties.maxThreadsPerBlock << std::endl;
		    std::cout << "targetNumBlocks = " << targetNumBlocks << std::endl;

			#endif //PRISMATIC_ENABLE_GPU
	    };

    };

	template <class T>
	void Parameters<T>::calculateLambda(){
		constexpr double m = 9.109383e-31;
		constexpr double e = 1.602177e-19;
		constexpr double c = 299792458;
		constexpr double h = 6.62607e-34;
		lambda = (T)(h / sqrt(2 * m * e * meta.E0) / sqrt(1 + e * meta.E0 / 2 / m / c / c) * 1e10);
	}

	// template <class T>
	// void Parameters<T>::calculateFileSize(){

	// 	//calc detector size
	// 	long long ncx = (long long) floor((PRISMATIC_FLOAT_PRECISION) imageSize[1] / 2);
	// 	PRISMATIC_FLOAT_PRECISION dpx = 1.0 / ((PRISMATIC_FLOAT_PRECISION)imageSize[1] * meta.realspacePixelSize[1]);
	// 	long long ncy = (long long) floor((PRISMATIC_FLOAT_PRECISION) imageSize[0] / 2);
	// 	PRISMATIC_FLOAT_PRECISION dpy = 1.0 / ((PRISMATIC_FLOAT_PRECISION)imageSize[0] * meta.realspacePixelSize[0]);
	// 	PRISMATIC_FLOAT_PRECISION qmax_tmp = std::min(dpx*(ncx), dpy*(ncy)) / 2;

	// 	PRISMATIC_FLOAT_PRECISION alphaMax_tmp = qmax_tmp*lambda;
	// 	std::vector<PRISMATIC_FLOAT_PRECISION> detectorAngles_d = vecFromRange(meta.detectorAngleStep / 2, meta.detectorAngleStep, alphaMax_tmp - meta.detectorAngleStep / 2);
	// 	size_t Ndet_tmp = detectorAngles_d.size();

	// 	//calc num probes
	// 	Array1D<PRISMATIC_FLOAT_PRECISION> xR = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{2}});
	// 	xR[0] = scanWindowXMin * tiledCellDim[2];
	// 	xR[1] = scanWindowXMax * tiledCellDim[2];
	// 	Array1D<PRISMATIC_FLOAT_PRECISION> yR = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{2}});
	// 	yR[0] = scanWindowYMin * tiledCellDim[1];
	// 	yR[1] = scanWindowYMax * tiledCellDim[1];

	// 	std::vector<PRISMATIC_FLOAT_PRECISION> xp_d = vecFromRange(xR[0], meta.probeStepX, xR[1]);
	// 	std::vector<PRISMATIC_FLOAT_PRECISION> yp_d = vecFromRange(yR[0], meta.probeStepY, yR[1]);
	// 	size_t numProbes = xp_d.size() + yp_d.size();

	// 	long long uint numElems = 0;
	// 	long long uint mult = (meta.saveComplexOutputWave) ? 2 : 1;
	// 	if(meta.save2DOutput) numElems += numProbes*mult;
	// 	if(meta.save3DOutput) numElems += numProbes*Ndet_tmp*mult;
	// 	if(meta.saveDPC_CoM) numElems += 2*numProbes*(mult % 2); 
	// 	if(meta.savePotentialSlices) numElems += imageSize[0]*imageSize[1]*std::ceil(tiledCellDim[0]/meta.sliceThickness);

	// 	if(meta.algorithm == Algorithm::Multislice)
	// 	{
	// 		//TODO: num depth outputs
	// 		if(meta.save4DOutput) numElems += numProbes*imageSize[0]*imageSize[1]*mult/4;
	// 	}
	// 	else if(meta.algorithm == Algorithm::PRISM)
	// 	{
	// 		size_t numBeams = imageSize[0]*imageSize[1]/(4*meta.interpolationFactorX*meta.interpolationFactorY);
	// 		if(meta.save4DOutput) numElems += numProbes*imageSize[0]*imageSize[1]*mult/(4*meta.interpolationFactorX*meta.interpolationFactorY);
	// 		if(meta.saveSMatrix) numElems += numBeams*imageSize[0]*imageSize[1]/4; 
	// 	}
	// 	else if(meta.algorithm == Algorithm::HRTEM)
	// 	{
	// 		//TODO: calculate number of tilts here
	// 		numElems = 0;
	// 		size_t numBeams = 1;
	// 		if(meta.savePotentialSlices) numElems += imageSize[0]*imageSize[1]*std::ceil(tiledCellDim[0]/meta.sliceThickness);
	// 		numElems += imageSize[0]*imageSize[1]*numBeams*mult/2; 
	// 	}

	// 	std::cout << "app. output file size is: " << numElems*sizeof(PRISMATIC_FLOAT_PRECISION) << std::endl;
	// }

	template <class T>
	class Parameters_C : Parameters<T> {
	
	public:
		Array4D<std::complex<T>> output;
	};

}
#endif //PRISM_PARAMS_H
