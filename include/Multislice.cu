#include "Multislice.cuh"
#include "Multislice.h"
#include "cuComplex.h"
#include "cufft.h"

#include "fftw3.h"
#include "utility.h"
#define NX 64
#define NY 64
#define NZ 128

namespace PRISM{
    __host__ void getMultisliceProbe_gpu(){
cufftHandle plan;
cufftComplex *data1, *data2;
cudaMalloc((void**)&data1, sizeof(cufftComplex)*NX*NY*NZ);
cudaMalloc((void**)&data2, sizeof(cufftComplex)*NX*NY*NZ);
/* Create a 3D FFT plan. */
cufftPlan3d(&plan, NX, NY, NZ, CUFFT_C2C);
cufftDestroy(plan);
}
    __host__ void buildMultisliceOutput_gpu(Parameters <PRISM_FLOAT_PRECISION> &pars,
                                   Array3D <std::complex<PRISM_FLOAT_PRECISION>> &trans,
                                   Array2D <std::complex<PRISM_FLOAT_PRECISION>> &PsiProbeInit,
                                   Array2D <PRISM_FLOAT_PRECISION> &alphaInd) {
        using namespace std;
        cout << "Test GPU function from CUDA host" << endl;
getMultisliceProbe_gpu();
    /*
static mutex fftw_plan_lock; // for synchronizing access to shared FFTW resources
		static const PRISM_FLOAT_PRECISION pi = acos(-1);
		static const std::complex<PRISM_FLOAT_PRECISION> i(0, 1);
		// populates the output stack for Multislice simulation using the CPU. The number of
		// threads used is determined by pars.meta.NUM_THREADS

		Array2D<complex<PRISM_FLOAT_PRECISION> > psi(PsiProbeInit);

		{
			auto qxa_ptr = pars.qxa.begin();
			auto qya_ptr = pars.qya.begin();
            for (auto& p:psi)p*=exp(-2 * pi * i * ( (*qxa_ptr++)*pars.xp[ax] +
                                                    (*qya_ptr++)*pars.yp[ay]));
		}

        // fftw_execute is the only thread-safe function in the library, so we need to synchronize access
        // to the plan creation methods
        unique_lock<mutex> gatekeeper(fftw_plan_lock);
        fftwf_plan plan_forward = fftwf_plan_dft_2d(psi.get_dimj(), psi.get_dimi(),
                                                    reinterpret_cast<fftwf_complex *>(&psi[0]),
                                                    reinterpret_cast<fftwf_complex *>(&psi[0]),
                                                    FFTW_FORWARD, FFTW_ESTIMATE);
        fftwf_plan plan_inverse = fftwf_plan_dft_2d(psi.get_dimj(), psi.get_dimi(),
                                                    reinterpret_cast<fftwf_complex *>(&psi[0]),
                                                    reinterpret_cast<fftwf_complex *>(&psi[0]),
                                                    FFTW_BACKWARD, FFTW_ESTIMATE);
        gatekeeper.unlock(); // unlock it so we only block as long as necessary to deal with plans

        for (auto a2 = 0; a2 < pars.numPlanes; ++a2){
            fftwf_execute(plan_inverse);
            complex<PRISM_FLOAT_PRECISION>* t_ptr = &trans[a2 * trans.get_dimj() * trans.get_dimi()];
            for (auto& p:psi)p *= (*t_ptr++); // transmit
            fftwf_execute(plan_forward);
            auto p_ptr = pars.prop.begin();
            for (auto& p:psi)p *= (*p_ptr++); // propagate
            for (auto& p:psi)p /= psi.size(); // scale FFT

        }

        Array2D<PRISM_FLOAT_PRECISION> intOutput = zeros_ND<2, PRISM_FLOAT_PRECISION>({{psi.get_dimj(), psi.get_dimi()}});
        auto psi_ptr = psi.begin();
        for (auto& j:intOutput) j = pow(abs(*psi_ptr++),2);
        //update stack -- ax,ay are unique per thread so this write is thread-safe without a lock
        auto idx = alphaInd.begin();
        for (auto counts = intOutput.begin(); counts != intOutput.end(); ++counts){
            if (*idx <= pars.Ndet){
                pars.stack.at(ay,ax,(*idx)-1, 1) += *counts;
            }
            ++idx;
        };
*/



};
}
