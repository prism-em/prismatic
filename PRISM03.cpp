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
                    }
                }
            }
        }
    }
}