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
    __host__ void getMultisliceProbe_gpu(Parameters<PRISM_FLOAT_PRECISION>& pars,
                                Array3D<complex<PRISM_FLOAT_PRECISION> >& trans,
                                const Array2D<complex<PRISM_FLOAT_PRECISION> >& PsiProbeInit,
                                const size_t& ay,
                                const size_t& ax,
                                Array2D<PRISM_FLOAT_PRECISION> &alphaInd){
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
	const PRISM_FLOAT_PRECISION cpu_stop = std::floor(pars.meta.cpu_gpu_ratio*pars.yp.size());
		vector<thread> workers_gpu;
		vector<thread> workers_cpu;
		workers_gpu.reserve(pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU); // prevents multiple reallocations
		workers_cpu.reserve(pars.meta.NUM_THREADS); // prevents multiple reallocations
		auto WORK_CHUNK_SIZE_GPU = std::floor( (((pars.yp.size() - cpu_stop - 1)) / (pars.meta.NUM_GPUS * pars.meta.NUM_STREAMS_PER_GPU)) + 1); //TODO: divide work more generally than just splitting up by yp. If input isn't square this might not do a good job
		cout << "WORK_CHUNK_SIZE_GPU = " << WORK_CHUNK_SIZE_GPU << endl;
		auto start = cpu_stop;// gpu work starts where the cpu work will stop
		auto stop = start + WORK_CHUNK_SIZE_GPU;
		while (start < pars.yp.size()) {
			cout << "Launching thread to compute all x-probe positions for y-probes "
				 << start << "/" << std::min((size_t)stop,pars.yp.size()) << " on GPU\n";
			// emplace_back is better whenever constructing a new object
			workers_gpu.emplace_back(thread([&pars, &trans,
												&alphaInd, &PsiProbeInit,
												start, stop]() {
				for (auto ay = start; ay < std::min((size_t) stop, pars.yp.size()); ++ay) {
					for (auto ax = 0; ax < pars.xp.size(); ++ax) {
						getMultisliceProbe_gpu(pars, trans, PsiProbeInit, ay, ax, alphaInd);
					}
				}
			}));
			start += WORK_CHUNK_SIZE_GPU;
			if (start >= pars.yp.size())break;
			stop += WORK_CHUNK_SIZE_GPU;
		}


		// now launch CPU work
		auto WORK_CHUNK_SIZE_CPU = std::floor(((cpu_stop - 1) / pars.meta.NUM_THREADS) + 1); //TODO: divide work more generally than just splitting up by yp. If input isn't square this might not do a good job
		cout << "WORK_CHUNK_SIZE_CPU = " << WORK_CHUNK_SIZE_CPU << endl;
                start = 0;// cpu work starts at beginning
                stop = start + WORK_CHUNK_SIZE_CPU;
                while (start < cpu_stop) {
                        cout << "Launching thread to compute all x-probe positions for y-probes "
                                 << start << "/" << std::min(stop,cpu_stop) << " on CPU\n";
                        // emplace_back is better whenever constructing a new object
                        workers_gpu.emplace_back(thread([&pars, &trans,
                                                                                                &alphaInd, &PsiProbeInit,
                                                                                                start, cpu_stop,stop]() {
                                for (auto ay = start; ay < std::min(stop, cpu_stop); ++ay) {
                                        for (auto ax = 0; ax < pars.xp.size(); ++ax) {
                                                getMultisliceProbe_cpu(pars, trans, PsiProbeInit, ay, ax, alphaInd);
                                        }
                                }
                        }));
                        start += WORK_CHUNK_SIZE_CPU;
                        if (start >= cpu_stop)break;
                        stop += WORK_CHUNK_SIZE_CPU;
                }

		// synchronize threads
		for (auto& t:workers_gpu)t.join();
		for (auto& t:workers_cpu)t.join();

	}
}
