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
	    using Array1D      = PRISM::ArrayND<1, std::vector<T> >;
	    using Array1D_dims = PRISM::ArrayND<1, std::vector<size_t> >;
	    using Array2D      = PRISM::ArrayND<2, std::vector<T> >;
	    using Array2D_cx   = PRISM::ArrayND<2, std::vector< std::complex<T> > >;
	    using Array3D      = PRISM::ArrayND<3, std::vector<T> >;
	    using Array3D_cx   = PRISM::ArrayND<3, std::vector< std::complex<T> > >;

	    Array3D_cx Scompact;
	    Array3D stack;
		Array3D pot;

        Array2D probeDefocusArray;
        Array2D probeSemiangleArray;
        Array2D probeXtiltArray;
        Array2D probeYtiltArray;
        Array2D qxaReduce;
        Array2D qyaReduce;
        Array2D xp;
        Array2D yp;
        Array2D beamsIndex;
        Array2D xyBeams;
        Array2D xVec;
        Array2D yVec;
        Array2D detectorAngles;
	    Array2D u;
	    std::vector<atom> atoms;
	    Array1D pixelSize;
        Array1D pixelSizeOutput;
	    Array1D_dims cellDim;
	    Array1D_dims imageSize;
	    Array1D_dims imageSizeReduce;
	    Array1D_dims imageSizeOutput;

	    Array2D_cx PsiProbeInit;
        Array2D q1;
        Array2D q2;

	    T scale;
        T lambda;
        T dr;
        T dq;
        T potBound;
        size_t Ndet;
        size_t numFP;
        size_t sliceThickness;
        size_t interpolationFactor;
        size_t numPlanes;
	    emdSTEM(){};
    };


}
#endif //PRISM_EMDSTEM_H
