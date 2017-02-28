//
// Created by AJ Pryor on 2/13/17.
//
#ifndef PRISM_EMDSTEM_H
#define PRISM_EMDSTEM_H

#include "ArrayND.h"
#include <vector>
#include <complex>
#include "atom.h"

namespace PRISM{
    template <class T>
    struct emdSTEM {
	    using Array1D        = PRISM::ArrayND<1, std::vector<T> >;
	    using Array1D_dims   = PRISM::ArrayND<1, std::vector<size_t> >;
	    using Array2D        = PRISM::ArrayND<2, std::vector<T> >;
	    using Array2D_cx     = PRISM::ArrayND<2, std::vector< std::complex<T> > >;
	    using Array2D_mask   = PRISM::ArrayND<2, std::vector<unsigned int> >;
	    using Array3D        = PRISM::ArrayND<3, std::vector<T> >;
	    using Array3D_cx     = PRISM::ArrayND<3, std::vector< std::complex<T> > >;

	    Array3D_cx Scompact;
	    Array3D stack;
		Array3D pot;

	    Array2D_cx prop;
	    Array2D_cx propBack;

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
        Array2D xyBeams;
		Array2D beams;
	    Array2D beamsOutput;
        Array2D xVec;
        Array2D yVec;
        Array1D detectorAngles;
	    Array1D u;
	    std::vector<atom> atoms;
	    Array1D pixelSize;
        Array1D pixelSizeOutput;
	    Array1D_dims cellDim;
	    Array1D_dims imageSize;
	    Array1D_dims imageSizeReduce;
	    Array1D_dims imageSizeOutput;
	    Array1D_dims qxInd;
	    Array1D_dims qyInd;

	    Array2D_cx PsiProbeInit;
        Array2D q1;
        Array2D q2;

	    T scale;
        T lambda;
        T dr;
        T dq;
        T potBound;
	    T E0;
	    T sigma;
	    T alphaBeamMax;
	    T qMax;
        size_t Ndet;
        size_t numFP;
        size_t sliceThickness;
        size_t interpolationFactor;
        size_t numPlanes;
	    size_t numberBeams;
		size_t NUM_THREADS;
	    size_t NUM_GPUS;
		emdSTEM(){};
    };


}
#endif //PRISM_EMDSTEM_H
