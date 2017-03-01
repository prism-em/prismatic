//
// Created by AJ Pryor on 2/13/17.
//

#include "include/PRISM03.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <mutex>
#include <numeric>
#include "fftw3.h"


using namespace std;
namespace PRISM {
    template<class T>
    using Array3D = PRISM::ArrayND<3, std::vector<T> >;
    template<class T>
    using Array2D = PRISM::ArrayND<2, std::vector<T> >;
    template<class T>
    using Array1D = PRISM::ArrayND<1, std::vector<T> >;

    template<class T>
    vector<T> vecFromRange(const T &start, const T &step, const T &stop) {
        vector<T> result;
        for (auto i = start; i <= stop; i += step) {
            result.push_back(i);
        }
        return result;
    };

    template<class T>
    Array2D<T> array2D_subset(const Array2D<T> &arr,
                              const size_t &starty, const size_t &stepy, const size_t &stopy,
                              const size_t &startx, const size_t &stepx, const size_t &stopx) {
        vector<T> _d;
        size_t dimx = 0;
        size_t dimy = 0;
        for (auto y = starty; y < stopy; y += stepy) {
            for (auto x = startx; x < stopx; x += stepx) {
                _d.push_back(arr.at(y, x));
            }
            ++dimy;
        }
        for (auto x = startx; x < stopx; x += stepx)++dimx;
        Array2D<T> result(_d, {{dimy, dimx}});
        return result;
    }

    template<class T>
    void PRISM03(emdSTEM<T> &pars) {
		// compute final image

	    cout << "Entering PRISM02" << endl;
        pars.probeDefocusArray = zeros_ND<1, T>({{1}});
        pars.probeSemiangleArray = zeros_ND<1, T>({{1}});
        pars.probeXtiltArray = zeros_ND<1, T>({{1}});
        pars.probeYtiltArray = zeros_ND<1, T>({{1}});
        pars.probeDefocusArray[0] = 0.0;
        pars.probeSemiangleArray[0] = 20.0 / 1000;
        pars.probeXtiltArray[0] = 0.0 / 1000;
        pars.probeYtiltArray[0] = 0.0 / 1000;

        T dxy = 0.25 * 2;

        Array1D<T> xR = zeros_ND<1, T>({{2}});
        xR[0] = 0.1 * pars.cellDim[2];
        xR[1] = 0.9 * pars.cellDim[2];
        Array1D<T> yR = zeros_ND<1, T>({{2}});
        yR[0] = 0.1 * pars.cellDim[1];
        yR[1] = 0.9 * pars.cellDim[1];

        vector<T> xp_d = vecFromRange(xR[0] + dxy / 2, dxy, xR[1] - dxy / 2);
        vector<T> yp_d = vecFromRange(yR[0] + dxy / 2, dxy, yR[1] - dxy / 2);

        Array1D<T> xp(xp_d, {{xp_d.size()}});
        Array1D<T> yp(yp_d, {{yp_d.size()}});
        pars.xp = xp;
        pars.yp = yp;
        pars.dr = 2.5 / 1000;
        pars.alphaMax = pars.qMax * pars.lambda;

        vector<T> detectorAngles_d = vecFromRange(pars.dr / 2, pars.dr, pars.alphaMax - pars.dr / 2);
        Array1D<T> detectorAngles(detectorAngles_d, {{detectorAngles_d.size()}});
        pars.detectorAngles = detectorAngles;
//        bool flag_plot = 0;
//        bool flag_keep_beams = 0;
        T r_0 = pars.imageSizeOutput[0] / pars.interpolationFactor / 2;
        T r_1 = pars.imageSizeOutput[1] / pars.interpolationFactor / 2;
        vector<T> yVec_d = vecFromRange(-r_0, 1.0, r_0 - 1);
        vector<T> xVec_d = vecFromRange(-r_1, 1.0, r_1 - 1);
        Array1D<T> yVec(yVec_d,{{yVec_d.size()}});
        Array1D<T> xVec(xVec_d,{{xVec_d.size()}});
        pars.yVec = yVec;
        pars.xVec = xVec;

        Array2D<T> beamsReduce = array2D_subset(pars.beamsOutput,
                                                0, pars.interpolationFactor, pars.beamsOutput.get_dimj(),
                                                0, pars.interpolationFactor, pars.beamsOutput.get_dimi());

        vector<size_t> imageSizeReduce{beamsReduce.get_dimj(), beamsReduce.get_dimi()};
        pars.xyBeams = zeros_ND<2, long>({pars.beamsIndex.size(), 2});

        for (auto a0 = 1; a0 <= pars.beamsIndex.size(); ++a0) {
            for (auto y = 0; y < beamsReduce.get_dimj(); ++y) {
                for (auto x = 0; x < beamsReduce.get_dimi(); ++x) {
                    if (beamsReduce.at(y, x) == a0) {
                        pars.xyBeams.at(a0 - 1, 0) = y;
                        pars.xyBeams.at(a0 - 1, 1) = x;
                    }
                }
            }
        }

        pars.qxaReduce = array2D_subset(pars.qxaOutput,
                                        0, pars.interpolationFactor, pars.qxaOutput.get_dimj(),
                                        0, pars.interpolationFactor, pars.qxaOutput.get_dimi());
        pars.qyaReduce = array2D_subset(pars.qyaOutput,
                                        0, pars.interpolationFactor, pars.qyaOutput.get_dimj(),
                                        0, pars.interpolationFactor, pars.qyaOutput.get_dimi());
        pars.Ndet = pars.detectorAngles.size();
        pars.stack = zeros_ND<4, T>({{pars.yp.size(), pars.xp.size(), pars.Ndet, 1}});

        Array2D<T> q1 = zeros_ND<2, T>({{imageSizeReduce[0], imageSizeReduce[1]}});
        Array2D<T> q2 = zeros_ND<2, T>({{imageSizeReduce[0], imageSizeReduce[1]}});
        Array2D<T> intOutput = zeros_ND<2, T>({{imageSizeReduce[0], imageSizeReduce[1]}});
        Array2D<complex<T> > PsiProbeInit, psi;
        PsiProbeInit = zeros_ND<2, complex<T> >({{imageSizeReduce[0], imageSizeReduce[1]}});
        psi = zeros_ND<2, complex<T> >({{imageSizeReduce[0], imageSizeReduce[1]}});
        pars.imageSizeReduce = imageSizeReduce;
        pars.dq = (pars.qxaReduce.at(0, 1) + pars.qyaReduce.at(1, 0)) / 2;
        T scale = pow(pars.interpolationFactor, 4);
        pars.scale = scale;

        // Most of this is transcribed directly from the original MATLAB version.
        // The operators +, -, /, * return PRISM arrays by value, so to avoid unnecessary memory
        // allocations/copies for chained operations I try to do things like create variables
        // initially with at most one operation, and then perform in-place transforms if more is needed
        for (auto a0 = 0; a0 < pars.probeDefocusArray.size(); ++a0) {
            for (auto a1 = 0; a1 < pars.probeSemiangleArray.size(); ++a1) {
                T qProbeMax = pars.probeSemiangleArray[a1] / pars.lambda;
                for (auto a2 = 0; a2 < pars.probeXtiltArray.size(); ++a2) {
                    for (auto a3 = 0; a3 < pars.probeYtiltArray.size(); ++a3) {
                        Array2D<T> qxaShift = pars.qxaReduce - (pars.probeXtiltArray[a2] / pars.lambda);
                        Array2D<T> qyaShift = pars.qyaReduce - (pars.probeYtiltArray[a2] / pars.lambda);
                        transform(qxaShift.begin(), qxaShift.end(),
                                  qyaShift.begin(), q2.begin(),
                                  [](const T &a, const T &b) { return a * a + b * b; });
                        transform(q2.begin(), q2.end(),
                                  q1.begin(),
                                  [](const T &a) { return sqrt(a); });
                        Array2D<T> alphaInd(q1); // copy constructor more efficient than assignment
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
                                {alphaInd.get_dimj(), alphaInd.get_dimi()});
                        transform(alphaInd.begin(), alphaInd.end(),
                                  alphaMask.begin(),
                                  [&pars](const T &a) { return (a < pars.Ndet) ? 1 : 0; });

                        transform(PsiProbeInit.begin(), PsiProbeInit.end(),
                                  q1.begin(), PsiProbeInit.begin(),
                                  [&pars, &qProbeMax](std::complex<T> &a, T &q1_t) {
                                      a.real(erf((qProbeMax - q1_t) / (0.5 * pars.dq)) * 0.5 + 0.5);
                                      a.imag(0);
                                      return a;
                                  });


                        const static std::complex<T> i(0, 1);
                        // this might seem a strange way to get pi, but it's slightly more future proof
                        const static double pi = std::acos(-1);
                        transform(PsiProbeInit.begin(), PsiProbeInit.end(),
                                  q2.begin(), PsiProbeInit.begin(),
                                  [&pars, &a0](std::complex<T> &a, T &q2_t) {
                                      a = a * exp(-i * pi * pars.lambda * pars.probeDefocusArray[a0] * q2_t);
                                      return a;
                                  });
                        T norm_constant = sqrt(accumulate(PsiProbeInit.begin(), PsiProbeInit.end(),
                                                          0.0, [](T accum, std::complex<T> &a) {
                                    return accum + abs(a) * abs(a);
                                })); // make sure to initialize with 0.0 and NOT 0 or it won't be a float and answer will be wrong
                        double a = 0;
                        for (auto &i : PsiProbeInit) { a += i.real(); };
                        transform(PsiProbeInit.begin(), PsiProbeInit.end(),
                                  PsiProbeInit.begin(), [&norm_constant](std::complex<T> &a) {
                                    return a / norm_constant;
                                });


                        T zTotal = pars.cellDim[0];
                        T xTiltShift = -zTotal * tan(pars.probeXtiltArray[a3]);
                        T yTiltShift = -zTotal * tan(pars.probeYtiltArray[a3]);

                        // launch threads to compute results for batches of xp, yp
                        // I do this by dividing the xp points among threads, and each computes
                        // all of the relevant yp for each of its xp. This seems an okay strategy
                        // as long as the number of xp and yp are similar. If that is not the case
                        // this may need to be adapted
                        vector<thread> workers;
                        workers.reserve(pars.NUM_THREADS); // prevents multiple reallocations
                        auto WORK_CHUNK_SIZE = ((pars.yp.size() - 1) / pars.NUM_THREADS) + 1;
                        auto start = 0;
                        auto stop = start + WORK_CHUNK_SIZE;
                        while (start < pars.yp.size()) {
	                        cout << "Launching thread to compute all x-probe positions for y-probes "
	                             << start << "/" << min(stop,pars.yp.size()) << '\n';
                            // emplace_back is better whenever constructing a new object
                            workers.emplace_back(thread([&pars, &xTiltShift, &yTiltShift,
                                                                &alphaInd, &PsiProbeInit,
                                                                start, stop]() {
                                for (auto ay = start; ay < min((size_t) stop, pars.yp.size()); ++ay) {
                                    for (auto ax = 0; ax < pars.xp.size(); ++ax) {
                                        buildSignal(pars, ay, ax, yTiltShift, xTiltShift, alphaInd, PsiProbeInit);
                                    }
                                }
                            }));
                            start += WORK_CHUNK_SIZE;
                            if (start >= pars.yp.size())break;
                            stop += WORK_CHUNK_SIZE;
                        }

                        // synchronize
	                    cout << "Waiting for threads...\n";
                        for (auto &t:workers)t.join();
                    }
                }
            }
        }
    }


    template<class T>
    void buildSignal(emdSTEM<T> &pars,
                     const size_t &ay,
                     const size_t &ax,
                     const T &yTiltShift,
                     const T &xTiltShift,
                     PRISM::ArrayND<2, std::vector<T> > &alphaInd,
                     PRISM::ArrayND<2, std::vector<std::complex<T> > > &PsiProbeInit) {
        static mutex fftw_plan_lock; // for synchronizing access to shared FFTW resources
        const static std::complex<T> i(0, 1);
        const static double pi = std::acos(-1);

        // convenience aliases
        using Array3D    = PRISM::ArrayND<3, std::vector<T> >;
        using Array2D    = PRISM::ArrayND<2, std::vector<T> >;
        using Array2D_cx = PRISM::ArrayND<2, std::vector<std::complex<T> > >;
        using Array1D    = PRISM::ArrayND<1, std::vector<T> >;

        T x0 = pars.xp[ax] / pars.pixelSizeOutput[1];
        T y0 = pars.yp[ay] / pars.pixelSizeOutput[0];
        Array1D x = pars.xVec + round(x0);
        transform(x.begin(), x.end(), x.begin(), [&pars](T &a) { return fmod(a, (T)pars.imageSizeOutput[1]); });
        Array1D y = pars.yVec + round(y0);
        transform(y.begin(), y.end(), y.begin(), [&pars](T &a) { return fmod(a, (T)pars.imageSizeOutput[0]); });
        Array2D intOutput = PRISM::zeros_ND<2, T>({{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
        for (auto a5 = 0; a5 < pars.numFP; ++a5) {
            Array2D_cx psi = PRISM::zeros_ND<2, std::complex<T> >({{pars.imageSizeReduce[0], pars.imageSizeReduce[1]}});
            for (auto a4 = 0; a4 < pars.beamsIndex.size(); ++a4) {
                T yB = pars.xyBeams.at(a4, 0);
                T xB = pars.xyBeams.at(a4, 1);
                if (abs(PsiProbeInit.at(yB, xB)) > 0) {
                    T q0_0 = pars.qxaReduce.at(yB, xB);
                    T q0_1 = pars.qyaReduce.at(yB, xB);
                    std::complex<T> phaseShift = exp(-2 * pi * i * (q0_0 * (pars.xp[ax] + xTiltShift) +
                                                                    q0_1 * (pars.yp[ay] + yTiltShift)));
//                    std::complex<T> phaseShift = exp(-2 * pi * i);
                    // caching this constant made a 5x performance improvement even with
                    // full compiler optimization turned on. Optimizing compilers aren't perfect...
                    const std::complex<T> tmp_const = PsiProbeInit.at(yB, xB) * phaseShift;
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
            fftw_plan plan = fftw_plan_dft_2d(psi.get_dimj(), psi.get_dimi(),
                                                  reinterpret_cast<fftw_complex *>(&psi[0]),
                                                  reinterpret_cast<fftw_complex *>(&psi[0]),
                                                  FFTW_FORWARD, FFTW_ESTIMATE);

            gatekeeper.unlock(); // unlock it so we only block as long as necessary to deal with plans
            fftw_execute(plan);
            gatekeeper.lock();
            fftw_destroy_plan(plan);
            gatekeeper.unlock();

            for (auto jj = 0; jj < intOutput.get_dimj(); ++jj){
                for (auto ii = 0; ii < intOutput.get_dimi(); ++ii){
                    intOutput.at(jj,ii) += pow(abs(psi.at(jj,ii)),2);
//	                intOutput.at(ii,jj) += pow(abs(psi.at(jj,ii)),2);
                }
            }
        }

//         update emdSTEM.stack -- ax,ay are unique per thread so this write is thread-safe without a lock
        auto idx = alphaInd.begin();
        for (auto counts = intOutput.begin(); counts != intOutput.end(); ++counts){
            if (*idx <= pars.Ndet){
                pars.stack.at(ay,ax,(*idx)-1, 1) += *counts * pars.scale;
            }
            ++idx;
        };
    }
}

