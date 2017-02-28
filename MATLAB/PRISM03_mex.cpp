#include <vector>
#include <stdexcept>
#include <string>
#include <complex>
#include "fftw3.h"
#include "mex.h"
#include "ArrayND.h"
#include "emdSTEM.h"
#include "PRISM03.h"

// define the position of each input here for clarity, since there are many
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
PRISM::ArrayND<3,  std::vector<T> > mat3DtoPRISM3D(const mxArray *array){
    // convert 3D MATLAB array to 3D PRISM array and convert from F to C order
    // TODO: implement reorder parameter like in mat3DtoPRISM3D_cx
    
    const mwSize ndims   = mxGetNumberOfDimensions(array);
    if (ndims!=3)mexErrMsgTxt("Array is not 3D.\n");
    const mwSize *dims   = mxGetDimensions(array);
    const size_t nrows   = (size_t)dims[0];
    const size_t ncols   = (size_t)dims[1];
    const size_t nlayers = (size_t)dims[2];
    const size_t N       = nrows*ncols*nlayers;
    double *ptr  = mxGetPr(array);
    PRISM::ArrayND<3, std::vector<T> > arr( std::vector<T>(N, 0), {nrows, ncols, nlayers});
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
PRISM::ArrayND<2, std::vector<T> > mat2DtoPRISM2D(const mxArray *array){
    // convert 2D MATLAB array to 2D PRISM array and convert from F to C order
    // TODO: implement reorder parameter like in mat3DtoPRISM3D_cx

    const mwSize ndims   = mxGetNumberOfDimensions(array);
    if (ndims!=2)mexErrMsgTxt("Array is not 2D.\n");
    const mwSize *dims   = mxGetDimensions(array);
    const size_t nrows   = (size_t)dims[0];
    const size_t ncols   = (size_t)dims[1];
    const size_t N       = nrows*ncols;
    double *ptr  = mxGetPr(array);
    PRISM::ArrayND<2, std::vector<T> > arr( std::vector<T>(N, 0), {nrows, ncols});
        for (auto j = 0; j < ncols; ++j){
            for (auto i = 0; i < nrows; ++i){
                arr.at(i,j) = (T)*ptr++;
            }
        }
    return arr;
}


template <typename T>
PRISM::ArrayND<1, std::vector<T> > mat1DtoPRISM1D(const mxArray *array){
    // convert 2D MATLAB array to 2D PRISM array and convert from F to C order
    // TODO: implement reorder parameter like in mat3DtoPRISM3D_cx

    const mwSize ndims   = mxGetNumberOfDimensions(array);
    if (ndims!=2)mexErrMsgTxt("Array is not 2D.\n");
    const mwSize *dims   = mxGetDimensions(array);
    const size_t nrows   = (size_t)dims[0];
    const size_t ncols   = (size_t)dims[1];
    const size_t N       = nrows*ncols;
    double *ptr  = mxGetPr(array);
    PRISM::ArrayND<1, std::vector<T> > arr( std::vector<T>(N, 0), {nrows*ncols});
    for (auto i = 0; i < nrows*ncols; ++i){
            arr.at(i) = (T)*ptr++;
    }
    return arr;
}


template <typename T>
T matGetScalar(const mxArray *array){
    // Get a scalar value
    
    double *ptr  = mxGetPr(array);
    return (T)*ptr;
}


template <typename T>
PRISM::ArrayND<3, std::vector< std::complex<T> > > mat3DtoPRISM3D_cx(const mxArray *array, bool reorder=1){
    // convert 3D MATLAB complex array to 3D PRISM complex array and convert from F to C order if reorder == 1
    
    const mwSize ndims   = mxGetNumberOfDimensions(array);
    if (ndims!=3)mexErrMsgTxt("Array is not 3D.\n");
    const mwSize *dims   = mxGetDimensions(array);
    const size_t nrows   = (size_t)dims[0];
    const size_t ncols   = (size_t)dims[1];
    const size_t nlayers = (size_t)dims[2];
    const size_t N       = nrows*ncols*nlayers;
    double *ptr_r  = mxGetPr(array);
    double *ptr_i  = mxGetPi(array);

    if (reorder) {
        PRISM::ArrayND<3, std::vector< std::complex<T> > > arr( std::vector< std::complex<T> >(N, 0), {nrows, ncols, nlayers});
        for (auto k = 0; k < nlayers; ++k) {
            for (auto j = 0; j < ncols; ++j) {
                for (auto i = 0; i < nrows; ++i) {
                    arr.at(i, j, k).real((T) *ptr_r++);
                    arr.at(i, j, k).imag((T) *ptr_i++);
                }
            }
        }
        return arr;
    } else
    {
        PRISM::ArrayND<3, std::vector< std::complex<T> > > arr( std::vector< std::complex<T> >(N, 0), {nlayers, ncols, nrows});
        for (auto &i:arr){
            i.real((T)*ptr_r++);
            i.imag((T)*ptr_i++);
        }
        return arr;
    }

}

template <typename T>
PRISM::ArrayND<2, std::vector< std::complex<T> > > mat2DtoPRISM2D_cx(const mxArray *array){
    // convert 2D MATLAB complex array to 2D PRISM complex array and convert from F to C order
    // TODO: implement reorder parameter like in mat3DtoPRISM3D_cx
    
    const mwSize ndims   = mxGetNumberOfDimensions(array);
    if (ndims!=2)mexErrMsgTxt("Array is not 2D.\n");
    const mwSize *dims   = mxGetDimensions(array);
    const size_t nrows   = (size_t)dims[0];
    const size_t ncols   = (size_t)dims[1];
    const size_t N       = nrows*ncols;
    double *ptr_r  = mxGetPr(array);
    double *ptr_i  = mxGetPi(array);
    PRISM::ArrayND<2, std::vector< std::complex<T> > > arr( std::vector< std::complex<T> >(N, 0), {nrows, ncols});
    for (auto j = 0; j < ncols; ++j){
        for (auto i = 0; i < nrows; ++i){
            arr.at(i,j).real((T)*ptr_r++);
            arr.at(i,j).imag((T)*ptr_i++);
        }
    }
    return arr;
}



void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]){
    
    // Check we for correct number of arguments
    if (nrhs != NUM_INPUTS){
        std::string what = "Incorrect number of inputs to mex function (expected " + std::to_string(NUM_INPUTS) + ')';
        mexErrMsgTxt(what.c_str());
    }
    
    // Define convenience aliases
    using PRISM_FLOAT_TYPE = double;
    using PRISM_COMPLEX_TYPE = std::complex<PRISM_FLOAT_TYPE>;
    using Array3D_r = PRISM::ArrayND<3, std::vector<PRISM_FLOAT_TYPE> >;
    using Array2D_r = PRISM::ArrayND<2, std::vector<PRISM_FLOAT_TYPE> >;
    using Array3D_cx = PRISM::ArrayND<3, std::vector<PRISM_COMPLEX_TYPE> >;
    using Array2D_cx = PRISM::ArrayND<2, std::vector<PRISM_COMPLEX_TYPE> >;
    using Array2D_dims = PRISM::ArrayND<2, std::vector<size_t> >;
    using Array1D      = PRISM::ArrayND<1, std::vector<PRISM_FLOAT_TYPE> >;
    using Array1D_dims = PRISM::ArrayND<1, std::vector<size_t> >;

    // extract all the data from MATLAB into PRISM data types
    Array3D_cx Scompact = mat3DtoPRISM3D_cx<PRISM_FLOAT_TYPE>(prhs[POS_Scompact], 0);
    Array3D_r stack     = mat3DtoPRISM3D<PRISM_FLOAT_TYPE>(prhs[POS_stack]);

    Array2D_r probeDefocusArray   = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_probeDefocusArray]);
    Array2D_r probeSemiangleArray = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_probeSemiangleArray]);
    Array2D_r probeXtiltArray     = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_probeXtiltArray]);
    Array2D_r probeYtiltArray     = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_probeYtiltArray]);
    Array2D_r qxaReduce           = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_qxaReduce]);
    Array2D_r qyaReduce           = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_qyaReduce]);
    Array2D_r xp                  = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_xp]);
    Array2D_r yp                  = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_yp]);
    Array2D_r beamsIndex          = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_beamsIndex]);
    Array2D_r xyBeams             = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_xyBeams]);
    Array2D_r xVec                = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_xVec]);
    Array2D_r yVec                = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_yVec]);
    Array1D_dims imageSizeReduce  = mat1DtoPRISM1D<size_t>(prhs[POS_imageSizeReduce]);
    Array1D_dims imageSizeOutput  = mat1DtoPRISM1D<size_t>(prhs[POS_imageSizeOutput]);
    Array2D_r detectorAngles      = mat2DtoPRISM2D<PRISM_FLOAT_TYPE>(prhs[POS_detectorAngles]);
    Array1D_dims cellDim          = mat1DtoPRISM1D<size_t>(prhs[POS_cellDim]);
    Array1D pixelSizeOutput       = mat1DtoPRISM1D<PRISM_FLOAT_TYPE>(prhs[POS_pixelSizeOutput]);

    PRISM_FLOAT_TYPE scale        = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_scale]);
    PRISM_FLOAT_TYPE lambda       = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_lambda]);
    PRISM_FLOAT_TYPE dr           = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_dr]);
    PRISM_FLOAT_TYPE dq           = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_dq]);
    PRISM_FLOAT_TYPE Ndet         = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_Ndet]);
    PRISM_FLOAT_TYPE numFP        = matGetScalar<PRISM_FLOAT_TYPE>(prhs[POS_numFP]);

    // an emdSTEM object holds all of this data in one place for convenience
    PRISM::emdSTEM<PRISM_FLOAT_TYPE> PRISM_pars;

    PRISM_pars.Scompact = Scompact;
    PRISM_pars.stack    = stack;

    PRISM_pars.probeDefocusArray   = probeDefocusArray;
    PRISM_pars.probeSemiangleArray = probeSemiangleArray;
    PRISM_pars.probeXtiltArray     = probeXtiltArray;
    PRISM_pars.probeYtiltArray     = probeYtiltArray;
    PRISM_pars.qxaReduce           = qxaReduce;
    PRISM_pars.qyaReduce           = qyaReduce;
    PRISM_pars.xp                  = xp;
    PRISM_pars.yp                  = yp;
    PRISM_pars.beamsIndex          = beamsIndex;
    PRISM_pars.xyBeams             = xyBeams;
    PRISM_pars.xVec                = xVec;
    PRISM_pars.yVec                = yVec;
    PRISM_pars.imageSizeReduce     = imageSizeReduce;
    PRISM_pars.imageSizeOutput     = imageSizeOutput;
    PRISM_pars.detectorAngles      = detectorAngles;
    PRISM_pars.cellDim             = cellDim;
    PRISM_pars.pixelSizeOutput     = pixelSizeOutput;

    PRISM_pars.PsiProbeInit = PRISM::zeros_ND<2, PRISM_COMPLEX_TYPE>({imageSizeReduce[0], imageSizeReduce[1]});
    PRISM_pars.q1           = PRISM::zeros_ND<2, PRISM_FLOAT_TYPE>({imageSizeReduce[0], imageSizeReduce[1]});
    PRISM_pars.q2           = PRISM::zeros_ND<2, PRISM_FLOAT_TYPE>({imageSizeReduce[0], imageSizeReduce[1]});

    PRISM_pars.scale  = scale;
    PRISM_pars.lambda = lambda;
    PRISM_pars.dr     = dr;
    PRISM_pars.dq     = dq;
    PRISM_pars.Ndet   = Ndet;
    PRISM_pars.numFP  = numFP;
    
    // run PRISM03
    PRISM::PRISM03<PRISM_FLOAT_TYPE>(PRISM_pars);

    // create output. TODO: make it an nD output instead of reshaping a 1D
    plhs[0] = mxCreateDoubleMatrix(stack.size(), 1, mxREAL);

    double * ptr_r = mxGetPr(plhs[0]);

    for (auto k = 0; k < PRISM_pars.stack.get_nlayers(); ++k){
        for (auto j = 0; j < PRISM_pars.stack.get_ncols(); ++j){
            for (auto i = 0; i < PRISM_pars.stack.get_nrows(); ++i){
                *ptr_r++ = PRISM_pars.stack.at(i,j,k);
            }
        }
    }
}


using namespace std;
namespace PRISM {

    // forward declare the helper function
    template<class T>
    void buildSignal(const emdSTEM<T> &pars, const size_t &ax, const size_t &ay);

    template<class T>
    void PRISM03(emdSTEM<T> &pars) {
////	    emdSTEM.probeDefocusArray = 0;%[-40 0 40];%[-80 -60 -40 -20 0 20];%[-40 -30 -20 -10 0];%[-40 -20 0];%[-100 -50 0];  % in Angstroms      % dim 4
////	    emdSTEM.probeSemiangleArray = 20/1000;%[30 20 10]/1000;%fliplr([25 50 75])/1000;%25/1000;  % rads      % dim 5
////	    emdSTEM.probeXtiltArray = 0/1000;  % rads           % dim 6
////	    emdSTEM.probeYtiltArray = 0/1000;  % rads           % dim 7
////	    dxy = 0.25 * 2;
////	    xR = [0.1 0.9]*emdSTEM.cellDim(1);
////	    yR = [0.1 0.9]*emdSTEM.cellDim(2);
////	    emdSTEM.xp = (xR(1)+dxy/2):dxy:(xR(2)-dxy/2);
////	    emdSTEM.yp = (yR(1)+dxy/2):dxy:(yR(2)-dxy/2);
////	    dr = 2.5 / 1000;
////	    alphaMax = emdSTEM.qMax * emdSTEM.lambda;
////	    emdSTEM.detectorAngles = (dr/2):dr:(alphaMax-dr/2);
////	    flag_plot = 0;
////	    flag_keep_beams = 0;
////	    r = emdSTEM.imageSizeOutput / emdSTEM.interpolationFactor / 2;
////	    xVec = ((-r(1)):(r(1)-1));
////	    yVec = ((-r(2)):(r(2)-1));
////	    % Downsampled beams
////	    emdSTEM.beamsReduce = emdSTEM.beamsOutput( ...
////	    1:(emdSTEM.interpolationFactor):end,...
////	    1:(emdSTEM.interpolationFactor):end);
////	    imageSizeReduce = size(emdSTEM.beamsReduce);
////	    xyBeams = zeros(length(emdSTEM.beamsIndex),2);
////	    for a0 = 1:emdSTEM.numberBeams;
////	    [~,ind] = min(abs(emdSTEM.beamsReduce(:) - a0));
////	    [xx,yy] = ind2sub(imageSizeReduce,ind);
////	    xyBeams(a0,:) = [xx yy];
////	    end
////	    // alias some types to avoid so much text
////        using Array3D = PRISM::ArrayND<3, std::vector<T> >;
////        using Array2D = PRISM::ArrayND<2, std::vector<T> >;
////	    qxaReduce = emdSTEM.qxaOutput( ...
////	    1:emdSTEM.interpolationFactor:end,...
////	    1:emdSTEM.interpolationFactor:end);
////	    qyaReduce = emdSTEM.qyaOutput( ...
////	    1:emdSTEM.interpolationFactor:end,...
////	    1:emdSTEM.interpolationFactor:end);
////	    Ndet = length(emdSTEM.detectorAngles);
////	    % Initialize pieces
////	    emdSTEM.stackSize = [ ...
////	    length(emdSTEM.xp) ...
////	    length(emdSTEM.yp) ...
////	    length(emdSTEM.detectorAngles) ...
////	    length(emdSTEM.probeDefocusArray) ...
////	    length(emdSTEM.probeSemiangleArray) ...
////	    length(emdSTEM.probeXtiltArray) ...
////	    length(emdSTEM.probeYtiltArray)];
////	    emdSTEM.stack = zeros(emdSTEM.stackSize,'single');
////	    q1 = zeros(imageSizeReduce);
////	    q2 = zeros(imageSizeReduce);
////	    dq = mean([qxaReduce(2,1) qyaReduce(1,2)]);
////	    PsiProbeInit = zeros(imageSizeReduce);
////	    psi = zeros(imageSizeReduce);
//	    intOutput = zeros(imageSizeReduce);
//	    % Main loops
//	    scale = emdSTEM.interpolationFactor^4;

        // Most of this is transcribed directly from the original MATLAB version.
        // The operators +, -, /, * return PRISM arrays by value, so to avoid unnecessary memory
        // allocations/copies for chained operations I try to do things like create variables
        // initially with at most one operation, and then perform in-place transforms if more is needed
        for (auto a0 = 0; a0 < pars.probeDefocusArray.size(); ++a0) {
            for (auto a1 = 0; a1 < pars.probeSemiangleArray.size(); ++a1) {
                T qProbeMax = pars.probeSemiangleArray[a0] / pars.lambda;
                for (auto a2 = 0; a2 < pars.probeXtiltArray.size(); ++a2) {
                    for (auto a3 = 0; a3 < pars.probeYtiltArray.size(); ++a3) {
                        Array2D qxaShift = pars.qxaReduce - (pars.probeXtiltArray[a2] / pars.lambda);
                        Array2D qyaShift = pars.qyaReduce - (pars.probeYtiltArray[a2] / pars.lambda);
                        transform(qxaShift.begin(), qxaShift.end(),
                                  qyaShift.begin(), pars.q2.begin(),
                                  [](const T &a, const T &b) { return a * a + b * b; });
                        transform(pars.q2.begin(), pars.q2.end(),
                                  pars.q1.begin(),
                                  [](const T &a) { return sqrt(a); });
                        Array2D alphaInd(pars.q1); // copy constructor more efficient than assignment
                        transform(alphaInd.begin(), alphaInd.end(),
                                  alphaInd.begin(),
                                  [&pars](const T &a) {
                                      return 1 + round((a * pars.lambda - pars.detectorAngles[0]) / pars.dr);
                                  });
                        transform(alphaInd.begin(), alphaInd.end(),
                                  alphaInd.begin(),
                                  [](const T &a) { return a < 1 ? 1 : a; });
                        PRISM::ArrayND<2, std::vector<unsigned short> > alphaMask(
                                std::vector<unsigned short>(alphaInd.size(), 0),
                                {alphaInd.get_nrows(), alphaInd.get_ncols()});
                        transform(alphaInd.begin(), alphaInd.end(),
                                  alphaMask.begin(),
                                  [&pars](const T &a) { return (a < pars.Ndet) ? 1 : 0; });

                        transform(pars.PsiProbeInit.begin(), pars.PsiProbeInit.end(),
                                  pars.q1.begin(), pars.PsiProbeInit.begin(),
                                  [&pars, &qProbeMax](std::complex<T> &a, T &q1_t) {
                                      a.real(erf((qProbeMax - q1_t) / (0.5 * pars.dq)) * 0.5 + 0.5);
                                      a.imag(0);
                                      return a;
                                  });

                        const static std::complex<T> i(0, 1);
                        // this might seem a strange way to get pi, but it's slightly more future proof
                        const static double pi = std::acos(-1);
                        transform(pars.PsiProbeInit.begin(), pars.PsiProbeInit.end(),
                                  pars.q2.begin(), pars.PsiProbeInit.begin(),
                                  [&pars, &a0](std::complex<T> &a, T &q2_t) {
                                      a = a * exp(-i * pi * pars.lambda * pars.probeDefocusArray[a0] * q2_t);
                                      return a;
                                  });
                        T norm_constant = sqrt(accumulate(pars.PsiProbeInit.begin(), pars.PsiProbeInit.end(),
                                                          0.0, [](T accum, std::complex<T> &a) {
                                    return accum + abs(a) * abs(a);
                                })); // make sure to initialize with 0.0 and NOT 0 or it won't be a float and answer will be wrong

                        double a = 0;
                        for (auto &i : pars.PsiProbeInit) { a += i.real(); };
                        transform(pars.PsiProbeInit.begin(), pars.PsiProbeInit.end(),
                                  pars.PsiProbeInit.begin(), [&norm_constant](std::complex<T> &a) {
                                    return a / norm_constant;
                                });

                        T zTotal = pars.cellDim[2];
                        T xTiltShift = -zTotal * tan(pars.probeXtiltArray[a3]);
                        T yTiltShift = -zTotal * tan(pars.probeYtiltArray[a3]);

                        // launch threads to compute results for batches of xp, yp
                        // I do this by dividing the xp points among threads, and each computes
                        // all of the relevant yp for each of its xp. This seems an okay strategy
                        // as long as the number of xp and yp are similar. If that is not the case
                        // this may need to be adapted
                        vector<thread> workers;
                        workers.reserve(pars.NUM_THREADS); // prevents multiple reallocations
                        auto WORK_CHUNK_SIZE = ( (pars.xp.size()-1) / pars.NUM_THREADS) + 1;
                        auto start = 0;
                        auto stop = start + WORK_CHUNK_SIZE;
                        while (start < pars.xp.size()) {
                            // emplace_back is better whenever constructing a new object
                            workers.emplace_back(thread([&pars, &xTiltShift, &yTiltShift, &alphaInd, start,stop]() {
                                for (auto ax = start; ax < min((size_t)stop, pars.xp.size()); ++ax) {
                                    for (auto ay = 0; ay < pars.yp.size(); ++ay) {
                                        buildSignal(pars, ax, ay, xTiltShift, yTiltShift, alphaInd);
                                    }
                                }
                            }));
                            start += WORK_CHUNK_SIZE;
                            if (start >= pars.xp.size())break;
                            stop  += WORK_CHUNK_SIZE;
                        }

                        // synchronize
                        for (auto &t:workers)t.join();
                    }
                }
            }
        }
    }



    template<class T>
    void buildSignal(emdSTEM<T> &pars, const size_t &ax, const size_t &ay, const T &xTiltShift, const T &yTiltShift, PRISM::ArrayND<2, std::vector<T> > &alphaInd) {
        static mutex fftw_plan_lock; // for synchronizing access to shared FFTW resources
        const static std::complex<T> i(0, 1);
        const static double pi = std::acos(-1);

        // convenience aliases
        using Array3D = PRISM::ArrayND<3, std::vector<T> >;
        using Array2D = PRISM::ArrayND<2, std::vector<T> >;
        using Array2D_cx = PRISM::ArrayND<2, std::vector<std::complex<T> > >;

        T x0 = pars.xp[ax] / pars.pixelSizeOutput[0];
        T y0 = pars.yp[ay] / pars.pixelSizeOutput[1];
        Array2D x = pars.xVec + round(x0);
        transform(x.begin(), x.end(), x.begin(), [&pars](T &a) { return fmod(a, pars.imageSizeOutput[0]); });
        Array2D y = pars.yVec + round(y0);
        transform(y.begin(), y.end(), y.begin(), [&pars](T &a) { return fmod(a, pars.imageSizeOutput[1]); });
        Array2D intOutput = PRISM::zeros_ND<2, T>({pars.imageSizeReduce[0], pars.imageSizeReduce[1]});
        for (auto a5 = 0; a5 < pars.numFP; ++a5) {
            Array2D_cx psi = PRISM::zeros_ND<2, std::complex<T> >({pars.imageSizeReduce[0], pars.imageSizeReduce[1]});
            for (auto a4 = 0; a4 < pars.beamsIndex.size(); ++a4) {
                T xB = pars.xyBeams.at(a4, 0) - 1;
                T yB = pars.xyBeams.at(a4, 1) - 1;
                if (abs(pars.PsiProbeInit.at(xB, yB)) > 0) {
                    T q0_0 = pars.qxaReduce.at(xB, yB);
                    T q0_1 = pars.qyaReduce.at(xB, yB);
                    std::complex<T> phaseShift = exp(-2 * pi * i * (q0_0 * (pars.xp[ax] + xTiltShift) +
                                                                    q0_1 * (pars.yp[ay] + yTiltShift)));

                    // caching this constant made a 5x performance improvement even with
                    // full compiler optimization turned on. Optimizing compilers aren't perfect...
                    const std::complex<T> tmp_const = pars.PsiProbeInit.at(xB, yB) * phaseShift;
                    auto psi_ptr = psi.begin();
                    for (auto j = 0; j < y.size(); ++j) {
                        for (auto i = 0; i < x.size(); ++i) {
                            // access contiguously for performance
                            *psi_ptr++ +=  (tmp_const * pars.Scompact.at(a4, y[j], x[i]));
                        }
                    }
                }
            }

            // fftw_execute is the only thread-safe function in the library, so we need to synchronize access
            // to the plan creation methods
            unique_lock<mutex> gatekeeper(fftw_plan_lock);
            fftw_plan plan = fftw_plan_dft_2d(psi.get_nrows(), psi.get_ncols(),
                                              reinterpret_cast<fftw_complex *>(&psi[0]),
                                              reinterpret_cast<fftw_complex *>(&psi[0]),
                                              FFTW_FORWARD, FFTW_ESTIMATE);

            gatekeeper.unlock(); // unlock it so we only block as long as necessary to deal with plans
            fftw_execute(plan);
            gatekeeper.lock();
            fftw_destroy_plan(plan);
            gatekeeper.unlock();

            // since I'm transcribing from MATLAB, here there is an F to C ordering change. Tiny cost for convenience
            for (auto ii = 0; ii < intOutput.get_nrows(); ++ii){
                for (auto jj = 0; jj < intOutput.get_ncols(); ++jj){
                    intOutput.at(ii,jj) += pow(abs(psi.at(jj,ii)),2);
                }
            }
        }

        // update emdSTEM.stack -- ax,ay are unique per thread so this write is thread-safe without a lock
        auto idx = alphaInd.begin();
        for (auto counts = intOutput.begin(); counts != intOutput.end(); ++counts){
            if (*idx <= pars.Ndet){
                pars.stack.at(ax,ay,(*idx)-1) += *counts * pars.scale;
            }
            ++idx;
        };
    }
}


