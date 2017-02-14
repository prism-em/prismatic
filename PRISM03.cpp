//
// Created by AJ Pryor on 2/13/17.
//

#include "include/PRISM03.h"
#include <iostream>
#include <algorithm>
#include <cmath>

using namespace std;
namespace PRISM {


    template <class T>
    void buildSignal(const emdSTEM<T>& pars, const size_t& ax, const size_t& ay);

    template <class T>
    void PRISM03(emdSTEM<T>& pars) {

        using Array3D = PRISM::Array3D< std::vector<T> >;
        using Array2D = PRISM::Array2D< std::vector<T> >;
        cout << "pars.scale = " << pars.scale << endl;
        cout << "pars.q1.nrows = " << pars.q1.get_nrows() << endl;
        cout << "pars.q1.nrows = " << pars.q1.get_ncols() << endl;
        cout << "pars.q1.size = " << pars.q1.size() << endl;

        for (auto a0 = 0; a0 < pars.probeDefocusArray.size(); ++a0){
            for (auto a1 = 0; a1 < pars.probeSemiangleArray.size(); ++a1){
                T qProbeMax = pars.probeSemiangleArray[a0] / pars.lambda;
                cout << "qProbeMax = " << qProbeMax << endl;
                for (auto a2 = 0; a2 < pars.probeXtiltArray.size(); ++a2){
                    for (auto a3 = 0; a3 < pars.probeYtiltArray.size(); ++a3){
                        Array2D qxaShift = pars.qxaReduce - (pars.probeXtiltArray[a2] / pars.lambda);
                        Array2D qyaShift = pars.qyaReduce - (pars.probeYtiltArray[a2] / pars.lambda);
                        transform(qxaShift.begin(),qxaShift.end(),
                                       qyaShift.begin(), pars.q2.begin(),
                                       [](const T& a, const T& b){return a*a + b*b;});
                        transform(pars.q2.begin(),pars.q2.end(),
                                       pars.q1.begin(),
                                       [](const T& a){return sqrt(a);});
                        Array2D alphaInd(pars.q1);
                        transform(alphaInd.begin(), alphaInd.end(),
                                  alphaInd.begin(),
                                  [&pars](const T& a){return 1 + round( (a*pars.lambda - pars.detectorAngles[0]) / pars.dr);});
                        transform(alphaInd.begin(), alphaInd.end(),
                                  alphaInd.begin(),
                                  [](const T& a){return a<1?1:a;});
                        PRISM::Array2D< std::vector<unsigned short> > alphaMask(std::vector<unsigned short>(alphaInd.size(),0),
                                                                      alphaInd.get_nrows(), alphaInd.get_ncols());
                        transform(alphaInd.begin(), alphaInd.end(),
                                  alphaMask.begin(),
                                  [&pars](const T& a){return (a<pars.Ndet)?1:0;});

                        cout << "PsiProbeInit.at(5,5) = " << pars.PsiProbeInit.at(5,5) << endl;
                        transform(pars.PsiProbeInit.begin(), pars.PsiProbeInit.end(),
                                  pars.q1.begin(),pars.PsiProbeInit.begin(),
                                  [&pars, &qProbeMax](std::complex<T>& a, T& q1_t){
                                      a.real(erf( (qProbeMax - q1_t) / (0.5 * pars.dq) ) *0.5 + 0.5);
                                      a.imag(0);
//                                      a.imag(erf( (qProbeMax - q1_t) / (0.5 * pars.dq) ) *0.5 + 0.5);
                                      return a;
                                  });
                        constexpr static std::complex<T> i(0,1);
                        constexpr double pi = std::acos(-1);
                        transform(pars.PsiProbeInit.begin(), pars.PsiProbeInit.end(),
                                  pars.q2.begin(),pars.PsiProbeInit.begin(),
                                  [&pars, &a0](std::complex<T>& a, T& q2_t){
                                      a = a * exp(-i * pi * pars.lambda * pars.probeDefocusArray[a0]*q2_t);
//                                      a.real(1);
//                                      a.imag(1);
                                      return a;
                                  });
                        cout << "PsiProbeInit.at(5,5) = " << pars.PsiProbeInit.at(5,5) << endl;
//                        T norm_constant = sqrt(accumulate(pars.PsiProbeInit.begin(),pars.PsiProbeInit.end(),
//                                                     0, [](T accum, std::complex<T>& a){return accum + abs(a) * abs(a);}));
                        T norm_constant = sqrt(accumulate(pars.PsiProbeInit.begin(),pars.PsiProbeInit.end(),
                                                          0.0, [](T accum, std::complex<T>& a){return accum + abs(a)*abs(a);}));

                        cout << " norm_constant = " << norm_constant << endl;
                        cout << "  PsiProbeInit.size()= " << pars.PsiProbeInit.size() << endl;
                        double a = 0;
                        for (auto&i : pars.PsiProbeInit){a+=i.real();};cout<<"a= " << a << endl;
                        transform(pars.PsiProbeInit.begin(), pars.PsiProbeInit.end(),
                                  pars.PsiProbeInit.begin(), [&norm_constant](std::complex<T>& a){
                                    return a / norm_constant;
                                });

                        T zTotal     = pars.cellDim[2];
                        T xTiltShift = -zTotal*tan(pars.probeXtiltArray[a3]);
                        T yTiltShift = -zTotal*tan(pars.probeYtiltArray[a3]);
                        for (auto ax = 0; ax < pars.xp.size(); ++ax){
                            cout << " ax  = " << ax << endl;
                            for (auto ay = 0; ay < pars.yp.size(); ++ay){
//                        for (auto ax = 0; ax < 1; ++ax){
//                            for (auto ay = 0; ay < 1; ++ay){
                                buildSignal(pars, ax, ay, xTiltShift, yTiltShift);
                            }
                        }
                        cout << "pars.Scompact.at(0,0,0) = " << pars.Scompact.at(0,0,0) << endl;
                        cout << "pars.Scompact.at(0,0,1) = " << pars.Scompact.at(0,0,1) << endl;

                        cout << "PsiProbeInit.at(5,5) = " << pars.PsiProbeInit.at(5,5) << endl;
                        cout << "alphaInd.at(5,5) = " << alphaInd.at(5,5) << endl;
                        cout << "alphaInd.at(0,0) = " << alphaInd.at(0,0) << endl;
                        cout << "qxaShift[101] = " << qxaShift[101] << endl;
                        cout << "pars.q2[1] = " << pars.q2.at(0,0) << endl;
                        cout << "pars.q2[101] = " << pars.q2.at(1,2) << endl;
                        cout << "pars.probeXtiltArray[a2] = " << pars.probeXtiltArray[a2] << endl;
                        cout << "pars.lambda = " << pars.lambda << endl;
                        cout << "pars.qxaReduce = " << pars.qxaReduce[1] << endl;
                        cout << "alphaMask.at(44,44) = " << alphaMask.at(44,44) << endl;
                        cout << "alphaMask.at(0,0) = " << alphaMask.at(0,0) << endl;
                        cout << "zTotal= " << zTotal << endl;
                        cout << "xTiltShift = " << xTiltShift << endl;
                        cout << " yTiltShift= " << yTiltShift << endl;
                    }
                }
            }
        }
    }

    template <class T>
    void buildSignal(emdSTEM<T>& pars, const size_t& ax, const size_t& ay, const T& xTiltShift, const T& yTiltShift){
//        cout << " ax  = " << ax << endl;
//        cout << " ay  = " << ay << endl;
        constexpr static std::complex<T> i(0,1);
        constexpr double pi = std::acos(-1);
        using Array3D = PRISM::Array3D< std::vector<T> >;
        using Array2D = PRISM::Array2D< std::vector<T> >;
        using Array2D_cx = PRISM::Array2D< std::vector< std::complex<T> > >;
        T x0 = pars.xp[ax] / pars.pixelSizeOutput[0]; // could not create these
        T y0 = pars.yp[ay] / pars.pixelSizeOutput[1];
        Array2D x = pars.xVec + round(x0) - 1;
        transform(x.begin(),x.end(),x.begin(),[&pars](T& a){return fmod(a, pars.imageSizeOutput[0]) + 1;});
        Array2D y = pars.yVec + round(y0) - 1;
        transform(y.begin(),y.end(),y.begin(),[&pars](T& a){return fmod(a, pars.imageSizeOutput[1]) + 1;});
        Array2D intOutput = PRISM::zeros_2D<T>(pars.imageSizeReduce[0], pars.imageSizeReduce[1]);

        for (auto a5 = 0; a5 < pars.numFP; ++a5){
            Array2D_cx psi = PRISM::zeros_2D< std::complex<T> >(pars.imageSizeReduce[0], pars.imageSizeReduce[1]);
            for (auto a4 = 0; a4 < pars.beamsIndex.size(); ++a4){
                T xB = pars.xyBeams.at(a4,0) - 1;
                T yB = pars.xyBeams.at(a4,1) - 1;
//                if (a5 ==0 & a4 ==0 & ax==0 & ay==0) {
////                    cout << "xB = " << xB << endl;
//                    cout << "yB = " << yB << endl;
//                }

                if (abs(pars.PsiProbeInit.at(xB, yB)) > 0){
                    T q0_0 = pars.qxaReduce.at(xB, yB);
                    T q0_1 = pars.qyaReduce.at(xB, yB);
                    std::complex<T> phaseShift = exp(-2*pi*i*(q0_0 * (pars.xp[ax] + xTiltShift) +
                                                q0_1 * (pars.yp[ax] + yTiltShift)));
//                    auto x_ptr = x.begin();
//                    auto y_ptr = y.begin();
                    std::complex<T> t = pars.PsiProbeInit.at(xB,yB) * phaseShift;
                    auto psi_ptr = psi.begin();
//                    if (a4==0 && a5==0)cout << "phase shift = " << phaseShift << endl;


                    for (auto j = 0; j < y.size(); ++j){
                        for (auto i = 0; i < x.size(); ++i) {
//                            if (i==0 && a4==0 && a5==0){
//                                cout << "val = " << t * pars.Scompact.at(a4, y[j], x[i]) << endl;
//                                cout << " t = " << t << endl;
//                                cout << "pars.Scompact.at(a4, y[j], x[i])  = " << pars.Scompact.at(a4, y[j], x[i])<< endl;
//                            }
//                            *psi_ptr += t;
//                            *psi_ptr += (pars.PsiProbeInit.at(xB,yB) * phaseShift);
                            *psi_ptr = *psi_ptr +  (t * pars.Scompact.at(a4, y[j], x[i]));
//                            *psi_ptr = *psi_ptr +  (t * pars.Scompact.at(x[i], y[j], a4));
                            ++psi_ptr;
// i+=std::complex<T>(1,1);
//                            i += (pars.PsiProbeInit.at(xB, yB) * phaseShift);
//                        i += (pars.PsiProbeInit.at(xB,yB) * phaseShift * pars.Scompact.at(*x_ptr++, *y_ptr++, a4));

                        }
                    }





//                    for (auto i = 0; i < x.size(); ++i){
//                        for (auto j = 0; j < y.size(); ++j) {
////                            if (i==0 && a4==0 && a5==0){
////                                cout << "val = " << t * pars.Scompact.at(a4, y[j], x[i]) << endl;
////                                cout << " t = " << t << endl;
////                                cout << "pars.Scompact.at(a4, y[j], x[i])  = " << pars.Scompact.at(a4, y[j], x[i])<< endl;
////                            }
////                            *psi_ptr += t;
////                            *psi_ptr += (pars.PsiProbeInit.at(xB,yB) * phaseShift);
//                            *psi_ptr = *psi_ptr +  (t * pars.Scompact.at(a4, y[j], x[i]));
////                            *psi_ptr = *psi_ptr +  (t * pars.Scompact.at(x[i], y[j], a4));
//                            ++psi_ptr;
//// i+=std::complex<T>(1,1);
////                            i += (pars.PsiProbeInit.at(xB, yB) * phaseShift);
////                        i += (pars.PsiProbeInit.at(xB,yB) * phaseShift * pars.Scompact.at(*x_ptr++, *y_ptr++, a4));
//
//                        }
//                    }
//                    if (a5 ==0 & a4 ==0 & ax==0 & ay==0) {
//                        cout << "q0_0 = " << q0_0 << endl;
//                        cout << "q0_1 = " << q0_1 << endl;
////                        cout << "phaseShift = " << phaseShift << endl;
//                    }
                }


            }
//            if (a5 ==0 && ax==0 && ay==0) {
//                cout << "DEBUG psi[0] = " << psi[0] << endl;
//
//                cout << "psi.at(0,0) = " << psi.at(0,0) << endl;
//                cout << "psi.at(0,1) = " << psi.at(0,1) << endl;
//                cout << "psi.at(1,0) = " << psi.at(1,0) << endl;
//                cout << "psi.at(1,1) = " << psi.at(1,1) << endl;
//                cout << "psi.at(6,5) = " << psi.at(6,5) << endl;
//            }
        }
//        if (ax==0 && ay==0) {
//            cout << "x[0] = " << x[0] << endl;
//
//        }
//        cout << "pars.Scompact.at(0,0,0) = " << pars.Scompact.at(0,0,0) << endl;
//        cout << "pars.Scompact.at(0,0,1) = " << pars.Scompact.at(0,0,1) << endl;
//        cout << "pars.Scompact.at(0,0,2) = " << pars.Scompact.at(0,0,2) << endl;
//        cout << "pars.Scompact.at(1,1,0) = " << pars.Scompact.at(1,1,0) << endl;
//        cout << "pars.Scompact.at(4,3,0) = " << pars.Scompact.at(4,3,0) << endl;
//        cout << "pars.Scompact.get_nrows() = " << pars.Scompact.get_nrows() << endl;
//        cout << "pars.Scompact.get_ncols() = " << pars.Scompact.get_ncols() << endl;
//        cout << "pars.Scompact.get_nlayers() = " << pars.Scompact.get_nlayers() << endl;


//        cout << "ax = " << ax << endl = ;

    }
}