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
#define POS_imageSizeReduce     14
#define POS_imageSizeOutput     15
#define POS_detectorAngles      16
#define POS_cellDim             17
#define POS_pixelSizeOutput     18
#define POS_scale               19
#define POS_lambda              20
#define POS_dr                  21
#define POS_dq                  22
#define POS_Ndet                23
#define POS_numFP               24
#define NUM_INPUTS              25


template <typename T>
PRISM::Array3D< std::vector<T> > mat3DtoPRISM3D(const mxArray *array){
    // convert 3D MATLAB array to 3D PRISM array and convert from F to C order
    const mwSize ndims   = mxGetNumberOfDimensions(array);
    if (ndims!=3)mexErrMsgTxt("Array is not 3D.\n");
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
    //for (auto& i:arr){i=(T)*ptr++;}
    for (auto k = 0; k < nlayers; ++k){
        for (auto j = 0; j < ncols; ++j){
            for (auto i = 0; i < nrows; ++i){
                arr.at(i,j,k) = (T)*ptr++;
            }
        }
    }
    return arr;
}

template <typename T>
PRISM::Array2D< std::vector<T> > mat2DtoPRISM2D(const mxArray *array){
    const mwSize ndims   = mxGetNumberOfDimensions(array);
    if (ndims!=2)mexErrMsgTxt("Array is not 2D.\n");
    // convert 3D MATLAB array to 3D PRISM array and convert from F to C order
    const mwSize *dims   = mxGetDimensions(array);
    const size_t nrows   = (size_t)dims[0];
    const size_t ncols   = (size_t)dims[1];
    mexPrintf("nrows = %i\n",nrows);
    mexPrintf("ncols = %i\n",ncols);
    const size_t N       = nrows*ncols;
    double *ptr  = mxGetPr(array);
    PRISM::Array2D< std::vector<T> > arr( std::vector<T>(N, 0), nrows, ncols);
    //for (auto& i:arr)i=(T)*ptr++;
    //for (auto& i:arr){i=(T)*ptr++;}
        for (auto j = 0; j < ncols; ++j){
            for (auto i = 0; i < nrows; ++i){
                arr.at(i,j) = (T)*ptr++;
            }
        }
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
//    if (nrhs != NUM_INPUTS){
//        std::string what = "Incorrect number of inputs to mex function (expected " + std::to_string(NUM_INPUTS) + ')';
//        mexErrMsgTxt(what.c_str());
//    }
    using PRISM_FLOAT_TYPE = double;
    using Array3D = PRISM::Array3D< std::vector<PRISM_FLOAT_TYPE> >;
    using Array2D = PRISM::Array2D< std::vector<PRISM_FLOAT_TYPE> >;

    Array3D Scompact = mat3DtoPRISM3D<PRISM_FLOAT_TYPE>(prhs[POS_Scompact]);
    Array3D stack    = mat3DtoPRISM3D<PRISM_FLOAT_TYPE>(prhs[POS_stack]);

        Array2D probeDefocusArray   = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_probeDefocusArray]);
    Array2D probeSemiangleArray = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_probeSemiangleArray]);
    Array2D probeXtiltArray     = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_probeXtiltArray]);
    Array2D probeYtiltArray     = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_probeYtiltArray]);
    Array2D qxaReduce           = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_qxaReduce]);
    Array2D qyaReduce           = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_qyaReduce]);
    Array2D xp                  = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_xp]);
    Array2D yp                  = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_yp]);
    Array2D beamsIndex          = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_beamsIndex]);
    Array2D xyBeams             = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_xyBeams]);
    Array2D xVec                = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_xVec]);
    Array2D yVec                = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_yVec]);
    Array2D imageSizeReduce     = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_imageSizeReduce]);
    Array2D imageSizeOutput     = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_imageSizeOutput]);
    Array2D detectorAngles      = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_detectorAngles]);
    Array2D cellDim             = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_cellDim]);
    Array2D pixelSizeOutput     = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_pixelSizeOutput]);

    PRISM_FLOAT_TYPE scale           = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_scale]);
    PRISM_FLOAT_TYPE lambda          = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_lambda]);
    PRISM_FLOAT_TYPE dr              = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_dr]);
    PRISM_FLOAT_TYPE dq              = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_dq]);
    PRISM_FLOAT_TYPE Ndet            = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_Ndet]);
    PRISM_FLOAT_TYPE numFP           = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_numFP]);


   // PRISM::Array2D< std::vector<PRISM_FLOAT_TYPE> > qxaReduce = mat2DtoPRISM2D<double>(prhs[1]); // change from 0


//    for (auto &i : Scompact)mexPrintf("%f\n",i);
//    for (auto &i : probeDefocusArray)mexPrintf("%f\n",i);
    //plhs[0] = mxCreateDoubleMatrix(qxaReduce.get_nrows(), qxaReduce.get_ncols(), mxREAL);
    mexPrintf("Scompact size = %i\n",Scompact.size());
    mexPrintf("stack size = %i\n",stack.size());


    plhs[0] = mxCreateDoubleMatrix(Scompact.size(), 1, mxREAL);

    double * ptr = mxGetPr(plhs[0]);
//    for (auto &i:Scompact)*ptr++=i;
}