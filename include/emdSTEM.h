//
// Created by AJ Pryor on 2/13/17.
//
#ifndef PRISM_EMDSTEM_H
#define PRISM_EMDSTEM_H

#include "Array2D.h"
#include "Array3D.h"
#include <vector>
#include <complex>

namespace PRISM{
    template <class T>
    struct emdSTEM {
        PRISM::Array3D< std::vector< std::complex<T> > > Scompact;
        PRISM::Array3D< std::vector<T> > stack;
        PRISM::Array2D< std::vector<T> > probeDefocusArray;
        PRISM::Array2D< std::vector<T> > probeSemiangleArray;
        PRISM::Array2D< std::vector<T> > probeXtiltArray;
        PRISM::Array2D< std::vector<T> > probeYtiltArray;
        PRISM::Array2D< std::vector<T> > qxaReduce;
        PRISM::Array2D< std::vector<T> > qyaReduce;
        PRISM::Array2D< std::vector<T> > xp;
        PRISM::Array2D< std::vector<T> > yp;
        PRISM::Array2D< std::vector<T> > beamsIndex;
        PRISM::Array2D< std::vector<T> > xyBeams;
        PRISM::Array2D< std::vector<T> > xVec;
        PRISM::Array2D< std::vector<T> > yVec;
        PRISM::Array2D< std::vector<T> > imageSizeReduce;
        PRISM::Array2D< std::vector<T> > imageSizeOutput;
        PRISM::Array2D< std::vector<T> > detectorAngles;
        PRISM::Array2D< std::vector<T> > cellDim;
        PRISM::Array2D< std::vector<T> > pixelSizeOutput;

        PRISM::Array2D< std::vector< std::complex<T> > > PsiProbeInit;
        PRISM::Array2D< std::vector<T> > q1;
        PRISM::Array2D< std::vector<T> > q2;

        T scale;
        T lambda;
        T dr;
        T dq;
        T Ndet;
        T numFP;
        emdSTEM(){};
    };


}
#endif //PRISM_EMDSTEM_H
