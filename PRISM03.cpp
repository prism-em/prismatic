//
// Created by AJ Pryor on 2/13/17.
//

#include "include/PRISM03.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>
#include "fftw3.h"

#ifndef NUM_THREADS
#define NUM_THREADS 12
#endif //NUM_THREADS

using namespace std;
namespace PRISM {

    // forward declare the helper function
    template<class T>
    void buildSignal(const emdSTEM<T> &pars, const size_t &ax, const size_t &ay);

    template<class T>
    void PRISM03(emdSTEM<T> &pars) {

        // alias some types to avoid so much text
        using Array3D = PRISM::Array3D<std::vector<T> >;
        using Array2D = PRISM::Array2D<std::vector<T> >;

        // most of this is transcribed directly from the original MATLAB version. The operators +, -, /, * return
        // arrays by value, so to avoid unnecessary memory allocations/copies for chained operations I try to do things
        // like create variables initially with at most one operation, and then perform in-place transforms if more is needed
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
                        PRISM::Array2D<std::vector<unsigned short> > alphaMask(
                                std::vector<unsigned short>(alphaInd.size(), 0),
                                alphaInd.get_nrows(), alphaInd.get_ncols());
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

                        // constexpr means these are evaluated at compile time and no memory/overhead occurs at runtime
                        constexpr static std::complex<T> i(0, 1);
                        // this might seem a strange way to get pi, but it's slightly more future proof
                        constexpr static double pi = std::acos(-1);
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
                        workers.reserve(NUM_THREADS); // prevents multiple reallocations
                        auto WORK_CHUNK_SIZE = ( (pars.xp.size()-1) / NUM_THREADS) + 1;
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
                            stop  += WORK_CHUNK_SIZE;
                            if (start >= pars.xp.size())break;
                        }

                        // synchronize
                        for (auto &t:workers)t.join();
                    }
                }
            }
        }
    }

    static mutex fftw_plan_lock; // for synchronizing access to shared FFTW resources

    template<class T>
    void buildSignal(emdSTEM<T> &pars, const size_t &ax, const size_t &ay, const T &xTiltShift, const T &yTiltShift, PRISM::Array2D<std::vector<T> > &alphaInd) {
        // constexpr explained above
        constexpr static std::complex<T> i(0, 1);
        constexpr static double pi = std::acos(-1);

        // convenience aliases
        using Array3D = PRISM::Array3D<std::vector<T> >;
        using Array2D = PRISM::Array2D<std::vector<T> >;
        using Array2D_cx = PRISM::Array2D<std::vector<std::complex<T> > >;

        T x0 = pars.xp[ax] / pars.pixelSizeOutput[0];
        T y0 = pars.yp[ay] / pars.pixelSizeOutput[1];
        Array2D x = pars.xVec + round(x0);
        transform(x.begin(), x.end(), x.begin(), [&pars](T &a) { return fmod(a, pars.imageSizeOutput[0]); });
        Array2D y = pars.yVec + round(y0);
        transform(y.begin(), y.end(), y.begin(), [&pars](T &a) { return fmod(a, pars.imageSizeOutput[1]); });
        Array2D intOutput = PRISM::zeros_2D<T>(pars.imageSizeReduce[0], pars.imageSizeReduce[1]);
        for (auto a5 = 0; a5 < pars.numFP; ++a5) {
            Array2D_cx psi = PRISM::zeros_2D<std::complex<T> >(pars.imageSizeReduce[0], pars.imageSizeReduce[1]);
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
            fftw_plan plan = fftw_plan_dft_2d(psi.get_nrows(), psi.get_nrows(),
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