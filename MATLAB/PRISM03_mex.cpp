#include <vector>
#include <stdexcept>
#include <string>
#include "fftw3.h"
#include "mex.h"
#include "Array2D.h"
#include "Array3D.h"

// define the position of each input here for clarity, since there are so many
#define POS_Scompact            0
#define POS_stack               1
#define POS_probeDefocusArray   2
#define POS_probeSemiangleArray 3
#define POS_probeXtiltArray     4
#define POS_probeYtiltArray     5
#define POS_qxaReduce           6
#define POS_qyaReduce           7
#define POS_xp                  8
#define POS_yp                  9
#define POS_beamsIndex          10
#define POS_xyBeams             11
#define POS_xVec                12
#define POS_yVec                13
#define POS_imageSizeOutput     14
#define POS_detectorAngles      15
#define POS_cellDim             16
#define POS_scale               17
#define POS_lambda              18
#define POS_dr                  19
#define POS_dq                  20
#define POS_Ndet                21
#define POS_numFP               22
#define NUM_INPUTS              22


template <typename T>
PRISM::Array3D< std::vector<T> > mat3DtoPRISM3D(const mxArray *array){
    const mwSize *dims   = mxGetDimensions(array);
    const size_t nrows   = (size_t)dims[0];
    const size_t ncols   = (size_t)dims[1];
    const size_t nlayers = (size_t)dims[2];
    mexPrintf("nrows = %i\n",nrows);
    mexPrintf("ncols = %i\n",ncols);
    mexPrintf("nlayers = %i\n",nlayers);
    const size_t N       = nrows*ncols*nlayers;
    double *ptr  = mxGetPr(array);
    PRISM::Array3D< std::vector<T> > arr( std::vector<T>(N, 0), nrows, ncols, nlayers);
    //for (auto& i:arr)i=(T)*ptr++;
    for (auto& i:arr){i=(T)*ptr++;mexPrintf("%f\n",*ptr);}
    return arr;
}


template <typename T>
PRISM::Array2D< std::vector<T> > mat2DtoPRISM2D(const mxArray *array){
    const mwSize *dims = mxGetDimensions(array);
    const size_t nrows = (size_t)dims[0];
    const size_t ncols = (size_t)dims[1];
    const size_t N     = nrows*ncols;

    double *ptr  = mxGetPr(array);
    PRISM::Array2D< std::vector<T> > arr( std::vector<T>(N, 0), nrows, ncols);
    for (auto& i:arr)i=(T)*ptr++;
    return arr;
}

template <typename T>
T matGetScalar(const mxArray *array){
    double *ptr  = mxGetPr(array);
    return (T)*ptr;
}

template <typename T>
PRISM::Array2D< std::vector<T> > ones_2D(const size_t& nrows, const size_t& ncols){
    return PRISM::Array2D< std::vector<T> >(std::vector<T>(nrows*ncols,1),nrows, ncols);
}

template <typename T>
PRISM::Array2D< std::vector<T> > zeros_2D(const size_t& nrows, const size_t& ncols){
    return PRISM::Array2D< std::vector<T> >(std::vector<T>(nrows*ncols,0),nrows, ncols);
}

template <typename T>
PRISM::Array3D< std::vector<T> > ones_3D(const size_t& nrows, const size_t& ncols, const size_t& nlayers){
    return PRISM::Array3D< std::vector<T> >(std::vector<T>(nrows*ncols,1),nrows, ncols, nlayers);
}

template <typename T>
PRISM::Array3D< std::vector<T> > zeros_3D(const size_t& nrows, const size_t& ncols, const size_t& nlayers){
    return PRISM::Array3D< std::vector<T> >(std::vector<T>(nrows*ncols,0),nrows, ncols, nlayers);
}




void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    
    // Check we for correct number of arguments
//     if (nrhs != NUM_INPUTS)mexErrMsgTxt("Incorrect number of inputs to mex function, should be 22");
    using PRISM_FLOAT_TYPE = double;

    PRISM::Array3D< std::vector<PRISM_FLOAT_TYPE> > Scompact = mat3DtoPRISM3D<PRISM_FLOAT_TYPE>(prhs[0]); // change from 0
   // PRISM::Array2D< std::vector<PRISM_FLOAT_TYPE> > qxaReduce = mat2DtoPRISM2D<double>(prhs[1]); // change from 0


    for (auto &i : Scompact)mexPrintf("%f\n",i);
    //plhs[0] = mxCreateDoubleMatrix(qxaReduce.get_nrows(), qxaReduce.get_ncols(), mxREAL);
    mexPrintf("size = %i\n",Scompact.size());
    plhs[0] = mxCreateDoubleMatrix(Scompact.size(), 1, mxREAL);

    double * ptr = mxGetPr(plhs[0]);
    for (auto &i:Scompact)*ptr++=i;
}