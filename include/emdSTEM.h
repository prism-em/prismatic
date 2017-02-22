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
        PRISM::ArrayND<3, std::vector< std::complex<T> > > Scompact;
        PRISM::ArrayND<3, std::vector<T> > stack;

        PRISM::ArrayND<2, std::vector<T> > probeDefocusArray;
        PRISM::ArrayND<2, std::vector<T> > probeSemiangleArray;
        PRISM::ArrayND<2, std::vector<T> > probeXtiltArray;
        PRISM::ArrayND<2, std::vector<T> > probeYtiltArray;
        PRISM::ArrayND<2, std::vector<T> > qxaReduce;
        PRISM::ArrayND<2, std::vector<T> > qyaReduce;
        PRISM::ArrayND<2, std::vector<T> > xp;
        PRISM::ArrayND<2, std::vector<T> > yp;
        PRISM::ArrayND<2, std::vector<T> > beamsIndex;
        PRISM::ArrayND<2, std::vector<T> > xyBeams;
        PRISM::ArrayND<2, std::vector<T> > xVec;
        PRISM::ArrayND<2, std::vector<T> > yVec;
        PRISM::ArrayND<2, std::vector<T> > detectorAngles;

	    std::vector<atom> atoms;
	    PRISM::ArrayND<1, std::vector<T> > pixelSize;
        PRISM::ArrayND<1, std::vector<T> > pixelSizeOutput;
	    PRISM::ArrayND<1, std::vector<size_t> > cellDim;
	    PRISM::ArrayND<1, std::vector<size_t> > imageSize;
	    PRISM::ArrayND<1, std::vector<size_t> > imageSizeReduce;
	    PRISM::ArrayND<1, std::vector<size_t> > imageSizeOutput;

        PRISM::ArrayND<2, std::vector< std::complex<T> > > PsiProbeInit;
        PRISM::ArrayND<2, std::vector<T> > q1;
        PRISM::ArrayND<2, std::vector<T> > q2;

	    T scale;
        T lambda;
        T dr;
        T dq;
        T potBound;
        size_t Ndet;
        size_t numFP;
        size_t sliceThickness;
        size_t interpolationFactor;
        emdSTEM(){};
    };


}
#endif //PRISM_EMDSTEM_H
