//
// Created by AJ Pryor on 2/13/17.
//
#ifndef PRISM_PARAMS_H
#define PRISM_PARAMS_H
#include <vector>
#include <string>
#include <algorithm>
#include "ArrayND.h"
#include <complex>
#include "atom.h"
#include "meta.h"
namespace PRISM{
    template <class T>
    class Parameters {
    public:
	    using Array1D        = PRISM::ArrayND<1, std::vector<T> >;
	    using Array1D_dims   = PRISM::ArrayND<1, std::vector<size_t> >;
	    using Array2D        = PRISM::ArrayND<2, std::vector<T> >;
	    using Array2D_cx     = PRISM::ArrayND<2, std::vector< std::complex<T> > >;
	    using Array2D_mask   = PRISM::ArrayND<2, std::vector<unsigned int> >;
	    using Array2D_dims   = PRISM::ArrayND<2, std::vector<size_t> >;
	    using Array3D        = PRISM::ArrayND<3, std::vector<T> >;
	    using Array3D_cx     = PRISM::ArrayND<3, std::vector< std::complex<T> > >;
		using Array4D        = PRISM::ArrayND<4, std::vector<T> >;
	    Metadata<T> meta;
	    Array3D_cx Scompact;
	    Array4D stack;
		Array3D pot;


	    Array2D_cx prop;
	    Array2D_cx propBack;
//	    size_t interpolationFactor;
	    Array2D_mask qMask;
        Array1D probeDefocusArray;
        Array1D probeSemiangleArray;
        Array1D probeXtiltArray;
        Array1D probeYtiltArray;
	    Array2D qxa;
	    Array2D qya;
	    Array2D qxaOutput;
	    Array2D qyaOutput;
        Array2D qxaReduce;
        Array2D qyaReduce;
        Array1D xp;
        Array1D yp;
        std::vector<size_t> beamsIndex;
	    PRISM::ArrayND<2, std::vector<long> > xyBeams;
		Array2D beams;
	    Array2D beamsOutput;
        Array1D xVec;
        Array1D yVec;
        Array1D detectorAngles;
	    Array1D u;
	    std::vector<atom> atoms;
	    std::vector<T> pixelSize;
	    std::vector<T> pixelSizeOutput;
	    Array1D_dims imageSize;
	    std::vector<size_t> imageSizeReduce;
	    Array1D_dims imageSizeOutput;
	    Array1D_dims qxInd;
	    Array1D_dims qyInd;

	    T scale;
        T lambda;
        T dr;
        T dq;
	    T sigma;
	    T qMax;
	    T alphaMax;
        size_t Ndet;
        size_t numPlanes;
	    size_t numberBeams;
		Parameters(){};
	    Parameters(Metadata<T> _meta) : meta(_meta){

		    constexpr double m = 9.109383e-31;
		    constexpr double e = 1.602177e-19;
		    constexpr double c = 299792458;
		    constexpr double h = 6.62607e-34;
		    const double pi = std::acos(-1);
		    lambda = h / sqrt(2 * m * e * meta.E0) / sqrt(1 + e * meta.E0 / 2 / m / c / c) * 1e10;
		    sigma = (2 * pi / lambda / meta.E0) * (m * c * c + e * meta.E0) /
		                       (2 * m * c * c + e * meta.E0);

		    T f = 4 * meta.interpolationFactor;
		    Array1D_dims _imageSize({{meta.cellDim[1], meta.cellDim[2]}}, {{2}});
		    std::transform(_imageSize.begin(), _imageSize.end(), _imageSize.begin(),
		                   [&f, this](size_t &a) {
			                   return (size_t) (f * round(((T)a) / meta.realspace_pixelSize / f));
		                   });
		    this->imageSize = _imageSize;

		    std::vector<T> _pixelSize{(T) meta.cellDim[1], (T) meta.cellDim[2]};
		    pixelSize = _pixelSize;
		    pixelSize[0] /= (T)imageSize[0];
		    pixelSize[1] /= (T)imageSize[1];
		    try {
			    atoms = readAtoms(meta.filename_atoms);
		    }
		    catch (const std::runtime_error &e) {
			    std::cout << "PRISM: Error opening " << meta.filename_atoms << std::endl;
			    std::cout << e.what();
			    std::cout << "Terminating" << std::endl;
//			    return -1;
		    }
		    catch (const std::domain_error &e) {
			    std::cout << "PRISM: Error extracting atomic data from " << meta.filename_atoms << "!" << std::endl;
			    std::cout << e.what();
			    std::cout << "Terminating" << std::endl;
//			    return -2;
		    }

		    this->u = ones_ND<1, double>({{118}}) * 0.08;
//		    u = u;
		    std::cout << " prism_pars.pixelSize[1] = " << pixelSize[1] << std::endl;
		    std::cout << " prism_pars.pixelSize[0] = " << pixelSize[0] << std::endl;
	    };
    };


}
#endif //PRISM_PARAMS_H
