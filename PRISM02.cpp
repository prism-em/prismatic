//
// Created by AJ Pryor on 2/24/17.
//

#include "PRISM02.h"
#include <iostream>
#include <vector>
#include <thread>
#include "fftw3.h"
#include <mutex>
using namespace std;
namespace PRISM {

	template<class T>
	using Array3D = PRISM::ArrayND<3, std::vector<T> >;
	template<class T>
	using Array2D = PRISM::ArrayND<2, std::vector<T> >;
	template<class T>
	using Array1D = PRISM::ArrayND<1, std::vector<T> >;

	template<class T>
	Array1D<T> makeFourierCoords(const size_t &N, const T &pixel_size) {
		Array1D<T> result = zeros_ND<1, T>({N});
		long long nc = (size_t) floor((double) N / 2);

		T dp = 1 / (N * pixel_size);
		for (auto i = 0; i < N; ++i) {
			result[(nc + (size_t) i) % N] = (i - nc) * dp;
		}
		return result;
	};

	template<class T>
	void propagatePlaneWave(emdSTEM<T> &pars,
	                        Array3D<complex<T> >& trans,
	                        size_t a0,
	                        Array2D<complex<T> > &psi,
	                        const fftw_plan &plan_forward,
	                        const fftw_plan &plan_inverse,
							mutex& fftw_plan_lock){
		psi[pars.beamsIndex[a0]] = 1;
		const T N = (T)psi.size();
		fftw_execute(plan_inverse);
		for (auto &i : psi)i /= N; // fftw scales by N, need to correct
		const complex<T>* trans_t = &trans[0];
		for (auto a2 = 0; a2 < pars.numPlanes; ++a2){

			for (auto& p:psi)p*=(*trans_t++); // transmit
			fftw_execute(plan_forward); // FFT
			for (auto i = psi.begin(), j = pars.prop.begin(); i != psi.end();++i,++j)*i*=(*j); // propagate
			fftw_execute(plan_inverse); // IFFT
			for (auto &i : psi)i /= N; // fftw scales by N, need to correct
		}
		fftw_execute(plan_forward);

		Array2D< complex<T> > psi_small = zeros_ND<2, complex<T> >({{pars.qyInd.size(), pars.qxInd.size()}});
		for (auto y = 0; y < pars.qyInd.size(); ++y){
			for (auto x = 0; x < pars.qxInd.size(); ++x){
				psi_small.at(y,x) = psi.at(pars.qyInd[y], pars.qxInd[x]);
			}
		}

		unique_lock<mutex> gatekeeper(fftw_plan_lock);
		fftw_plan plan_final = fftw_plan_dft_2d(psi_small.get_nrows(), psi_small.get_ncols(),
		                                        reinterpret_cast<fftw_complex *>(&psi_small[0]),
		                                        reinterpret_cast<fftw_complex *>(&psi_small[0]),
		                                        FFTW_BACKWARD, FFTW_ESTIMATE);

		gatekeeper.unlock();
		fftw_execute(plan_final);
		gatekeeper.lock();
		fftw_destroy_plan(plan_final);
		gatekeeper.unlock();
	complex<T>* S_t = &pars.Scompact[a0 * pars.Scompact.get_nlayers() * pars.Scompact.get_ncols()];
	//cout << "pars.Scompact.get_nlayers() = " << pars.Scompact.get_nlayers() << endl;
	//	cout << "pars.Scompact.get_nrows() = " << pars.Scompact.get_nrows() << endl;
	//	cout << "pars.Scompact.get_ncols() = " << pars.Scompact.get_ncols() << endl;
		const T N_small = (T)psi_small.size();
		for (auto& i:psi_small) {
			*S_t++ = i/N_small;
		}
	};

	template<class T>
	void fill_Scompact(emdSTEM<T> &pars) {
		mutex fftw_plan_lock;

		const double pi = acos(-1);
		const std::complex<double> i(0, 1);
		pars.Scompact = zeros_ND<3, complex<T> > ({{pars.numberBeams,pars.imageSize[1], pars.imageSize[0]}});
		Array3D<complex<T> > trans = zeros_ND<3, complex<T> >(
				{{pars.pot.get_nrows(), pars.pot.get_ncols(), pars.pot.get_nlayers()}});
		{
			auto p = pars.pot.begin();
			for (auto &j:trans)j = exp(i * pars.sigma * (*p++));
		}

		for (auto& i:trans){i.real(1);i.imag(2);};
//
//
//		Array2D<complex<T> > psi = zeros_ND<2, complex<T> >({{pars.imageSize[1], pars.imageSize[0]}});
//		unique_lock<mutex> gatekeeper(fftw_plan_lock);
//		fftw_plan plan_forward = fftw_plan_dft_2d(psi.get_nrows(), psi.get_ncols(),
//		                                          reinterpret_cast<fftw_complex *>(&psi[0]),
//		                                          reinterpret_cast<fftw_complex *>(&psi[0]),
//		                                          FFTW_FORWARD, FFTW_ESTIMATE);
//		fftw_plan plan_inverse = fftw_plan_dft_2d(psi.get_nrows(), psi.get_ncols(),
//		                                          reinterpret_cast<fftw_complex *>(&psi[0]),
//		                                          reinterpret_cast<fftw_complex *>(&psi[0]),
//		                                          FFTW_BACKWARD, FFTW_ESTIMATE);
//		gatekeeper.unlock();
//		//for (auto a0 = 0; a0 < pars.numberBeams; ++a0){
//		for (auto a0 = 0; a0 < 1; ++a0) {
//			cout << "a0 = " << a0 << endl;
//			// re-zero psi each iteration
//			memset((void *) &psi[0], 0, psi.size() * sizeof(complex<T>));
//			propagatePlaneWave(pars, trans, a0, psi, plan_forward, plan_inverse, fftw_plan_lock);
//		}
//		gatekeeper.lock();
//		fftw_destroy_plan(plan_forward);
//		fftw_destroy_plan(plan_inverse);
//		gatekeeper.unlock();










		vector<thread> workers;
		workers.reserve(pars.NUM_THREADS); // prevents multiple reallocations
		auto WORK_CHUNK_SIZE = ( (pars.numberBeams-1) / pars.NUM_THREADS) + 1;
		auto start = 0;
		auto stop = start + WORK_CHUNK_SIZE;
		while (start < pars.numberBeams){

		workers.emplace_back([&pars, start, stop, &fftw_plan_lock, &trans](){
				// allocate array for psi just once per thread
				Array2D< complex<T> > psi = zeros_ND<2, complex<T> >({{pars.imageSize[1], pars.imageSize[0]}});

				unique_lock<mutex> gatekeeper(fftw_plan_lock);
				fftw_plan plan_forward = fftw_plan_dft_2d(psi.get_nrows(), psi.get_ncols(),
				                                          reinterpret_cast<fftw_complex *>(&psi[0]),
				                                          reinterpret_cast<fftw_complex *>(&psi[0]),
				                                          FFTW_FORWARD, FFTW_ESTIMATE);
				fftw_plan plan_inverse = fftw_plan_dft_2d(psi.get_nrows(), psi.get_ncols(),
				                                          reinterpret_cast<fftw_complex *>(&psi[0]),
				                                          reinterpret_cast<fftw_complex *>(&psi[0]),
				                                          FFTW_BACKWARD, FFTW_ESTIMATE);
				gatekeeper.unlock(); // unlock it so we only block as long as necessary to deal with plans
				for (auto a0 = start; a0 < min(stop, pars.numberBeams); ++a0){
					// re-zero psi each iteration
					memset((void*)&psi[0], 0, psi.size()*sizeof(complex<T>));
					propagatePlaneWave(pars, trans, a0, psi, plan_forward, plan_inverse, fftw_plan_lock);
				}
				// clean up
				gatekeeper.lock();
				fftw_destroy_plan(plan_forward);
				fftw_destroy_plan(plan_inverse);
				gatekeeper.unlock();
			});

			start += WORK_CHUNK_SIZE;
			if (start >= pars.numberBeams)break;
			stop  += WORK_CHUNK_SIZE;
		}
		for (auto& t:workers)t.join();


		cout << "pars.Scompact.at(0,0,0)= " << pars.Scompact.at(0,0,0) << endl;
		cout << "pars.Scompact.at(0,0,1)= " << pars.Scompact.at(0,0,1) << endl;
		cout << "pars.Scompact.at(0,0,2)= " << pars.Scompact.at(0,0,2) << endl;
		cout << "pars.Scompact.at(0,0,3)= " << pars.Scompact.at(0,0,3) << endl;
		cout << "pars.Scompact.at(1,0,0)= " << pars.Scompact.at(1,0,0) << endl;
		cout << "pars.Scompact.at(1,0,1)= " << pars.Scompact.at(1,0,1) << endl;
		cout << "pars.Scompact.at(1,0,2)= " << pars.Scompact.at(1,0,2) << endl;
		cout << "pars.Scompact.at(1,0,3)= " << pars.Scompact.at(1,0,3) << endl;
//		cout << "pars.Scompact.at(2,0,0)= " << pars.Scompact.at(2,0,0) << endl;
//		cout << "pars.Scompact.at(2,0,1)= " << pars.Scompact.at(2,0,1) << endl;
//		cout << "pars.Scompact.at(2,0,2)= " << pars.Scompact.at(2,0,2) << endl;
//		cout << "pars.Scompact.at(2,0,3)= " << pars.Scompact.at(2,0,3) << endl;
//		cout << "pars.Scompact.at(10,0,0)= " << pars.Scompact.at(10,0,0) << endl;
//		cout << "pars.Scompact.at(10,0,1)= " << pars.Scompact.at(10,0,1) << endl;
//		cout << "pars.Scompact.at(10,0,2)= " << pars.Scompact.at(10,0,2) << endl;
//		cout << "pars.Scompact.at(10,0,3)= " << pars.Scompact.at(10,0,3) << endl;
//		cout << "pars.Scompact.at(100,0,0)= " << pars.Scompact.at(100,0,0) << endl;
//		cout << "pars.Scompact.at(100,0,1)= " << pars.Scompact.at(100,0,1) << endl;
//		cout << "pars.Scompact.at(100,0,2)= " << pars.Scompact.at(100,0,2) << endl;
//		cout << "pars.Scompact.at(100,0,3)= " << pars.Scompact.at(100,0,3) << endl;

	}

	template <class T>
	void PRISM02(emdSTEM<T>& pars){
		// propogate plane waves to construct compact S-matrix

		// constants
		constexpr double m = 9.109383e-31;
		constexpr double e = 1.602177e-19;
		constexpr double c = 299792458;
		constexpr double h = 6.62607e-34;
        const double pi = acos(-1);
        const std::complex<double> i(0, 1);

		// setup some coordinates
		pars.imageSize[0] = pars.pot.get_nlayers();
		pars.imageSize[1] = pars.pot.get_ncols();
		Array1D<T> qx = makeFourierCoords(pars.imageSize[0], pars.pixelSize[0]);
		Array1D<T> qy = makeFourierCoords(pars.imageSize[1], pars.pixelSize[1]);

		pair< Array2D<T>, Array2D<T> > mesh = meshgrid(qx,qy);
		pars.qxa = mesh.first;
		pars.qya = mesh.second;
		Array2D<T> q2(pars.qya);
		transform(pars.qxa.begin(), pars.qxa.end(),
				  pars.qya.begin(), q2.begin(), [](const T& a, const T& b){
				return a*a + b*b;
				});

		// get qMax
		pars.qMax = 0;
		{
			T qx_max;
			T qy_max;
			for (auto i = 0; i < qx.size(); ++i) {
				qx_max = ( abs(qx[i]) > qx_max) ? abs(qx[i]) : qx_max;
				qy_max = ( abs(qy[i]) > qy_max) ? abs(qy[i]) : qy_max;
			}
			pars.qMax = min(qx_max, qy_max) / 2;
		}

		pars.qMask = zeros_ND<2, unsigned int>({pars.imageSize[1], pars.imageSize[0]});
		{
			long offset_x = pars.qMask.get_ncols()/4;		cout << "offset_x = " << offset_x << endl;
			long offset_y = pars.qMask.get_nrows()/4;
			long ndimy = (long)pars.qMask.get_nrows();
			long ndimx = (long)pars.qMask.get_nrows();
			for (long y = 0; y < pars.qMask.get_nrows() / 2; ++y) {
				for (long x = 0; x < pars.qMask.get_ncols() / 2; ++x) {
					pars.qMask.at( ((y-offset_y) % ndimy + ndimy) % ndimy,
								   ((x-offset_x) % ndimx + ndimx) % ndimx) = 1;
				}
			}
		}

        // build propagators
        pars.prop     = zeros_ND<2, std::complex<T> >({pars.imageSize[1], pars.imageSize[0]});
        pars.propBack = zeros_ND<2, std::complex<T> >({pars.imageSize[1], pars.imageSize[0]});
        cout << "exp(-i * pi) << endl =" << exp(-i * pi) << endl;
        for (auto y = 0; y < pars.qMask.get_nrows(); ++y) {
            for (auto x = 0; x < pars.qMask.get_ncols(); ++x) {
                if (pars.qMask.at(y,x)==1)
                {
                    pars.prop.at(y,x)     = exp(-i * pi * complex<T>(pars.lambda, 0) *
                                                          complex<T>(pars.sliceThickness, 0) *
                                                          complex<T>(q2.at(y, x), 0));
                    pars.propBack.at(y,x) = exp(i * pi * complex<T>(pars.lambda, 0) *
                                                         complex<T>(pars.cellDim[2], 0) *
                                                         complex<T>(q2.at(y, x), 0));
                }
            }
        }
        cout << "pars.propBack.at(3,4) = " << pars.propBack.at(3,4) << endl;
        cout << "pars.prop.at(33,44) = " << pars.prop.at(33,44) << endl;


        Array1D<T> xv = makeFourierCoords(pars.imageSize[1], (double)1/pars.imageSize[1]);
        Array1D<T> yv = makeFourierCoords(pars.imageSize[0], (double)1/pars.imageSize[0]);
        pair< Array2D<T>, Array2D<T> > mesh_a = meshgrid(xv, yv);

        // create beam mask and count beams
        PRISM::ArrayND<2, std::vector<unsigned int> > mask;
        mask = zeros_ND<2, unsigned int>({pars.imageSize[1], pars.imageSize[0]});
        pars.numberBeams = 0;
        long interp_f = (long)pars.interpolationFactor;
        for (auto y = 0; y < pars.qMask.get_ncols(); ++y) {
            for (auto x = 0; x < pars.qMask.get_nrows(); ++x) {
                if (q2.at(y,x) < pow(pars.alphaBeamMax / pars.lambda,2) &&
                    pars.qMask.at(y,x)==1 &&
                    (long)round(mesh_a.first.at(y,x))  % interp_f == 0 &&
                    (long)round(mesh_a.second.at(y,x)) % interp_f == 0){
                    mask.at(y,x)=1;
                    ++pars.numberBeams;
                }
            }
        }
		cout << "pars.alphaBeamMax = " << pars.alphaBeamMax<< endl;
		cout << "pars.lambda = " << pars.lambda<< endl;

        pars.beams     = zeros_ND<2, T>({{pars.imageSize[1], pars.imageSize[0]}});
		{
			int beam_count = 1;
			for (auto y = 0; y < pars.qMask.get_ncols(); ++y) {
				for (auto x = 0; x < pars.qMask.get_nrows(); ++x) {
					if (mask.at(y,x)==1){
						pars.beamsIndex.push_back((size_t)y*pars.qMask.get_ncols() + (size_t)x);
						pars.beams.at(y,x) = beam_count++;
					}
				}
			}
		}


		// TODO: ensure this block is correct for arbitrary dimension
		// get the indices for the compact S-matrix
		pars.qxInd = zeros_ND<1, size_t>({{pars.imageSize[1]/2}});
		pars.qyInd = zeros_ND<1, size_t>({{pars.imageSize[0]/2}});
		{
			long n_0        = pars.imageSize[0];
			long n_1        = pars.imageSize[1];
			long n_half0    = pars.imageSize[0]/2;
			long n_half1    = pars.imageSize[1]/2;
			long n_quarter0 = pars.imageSize[0]/4;
			long n_quarter1 = pars.imageSize[1]/4;
			for (auto i = 0; i < n_quarter0; ++i) {
				pars.qyInd[i] = i;
				pars.qyInd[i + n_quarter0] = (i-n_quarter0) + n_0;
			}
			for (auto i = 0; i < n_quarter1; ++i) {
				pars.qxInd[i] = i;
				pars.qxInd[i + n_quarter1] = (i-n_quarter1) + n_1;
			}
		}

		cout << "Fetching compact S" << endl;

		// populate compact-S matrix
		fill_Scompact(pars);

		// downsample Fourier components by x2 to match output
		pars.imageSizeOutput = pars.imageSize;
		pars.imageSizeOutput[0]/=2;
		pars.imageSizeOutput[1]/=2;
		pars.pixelSizeOutput = pars.pixelSize;
		pars.pixelSizeOutput[0]/=2;
		pars.pixelSizeOutput[1]/=2;

		pars.qxaOutput   = zeros_ND<2, T>({{pars.qyInd.size(), pars.qxInd.size()}});
		pars.qyaOutput   = zeros_ND<2, T>({{pars.qyInd.size(), pars.qxInd.size()}});
		pars.beamsOutput = zeros_ND<2, T>({{pars.qyInd.size(), pars.qxInd.size()}});
		for (auto y = 0; y < pars.qyInd.size(); ++y){
			for (auto x = 0; x < pars.qxInd.size(); ++x){
				pars.qxaOutput.at(y,x)   = pars.qxa.at(pars.qyInd[y],pars.qxInd[x]);
				pars.qyaOutput.at(y,x)   = pars.qya.at(pars.qyInd[y],pars.qxInd[x]);
				pars.beamsOutput.at(y,x) = pars.beams.at(pars.qyInd[y],pars.qxInd[x]);
			}
		}

	/*
		cout << "pars.qxaOutput.at(0,1) = " << pars.qxaOutput.at(0,1) << endl;
		cout << "pars.qyaOutput.at(4,3) = " << pars.qyaOutput.at(4,3) << endl;
		cout << "pars.beamsOutput.at(4,3) = " << pars.beamsOutput.at(4,3) << endl;
		pars.qxaOutput.toMRC_f("/mnt/spareA/clion/PRISM/MATLAB/qxa.mrc");
		pars.qyaOutput.toMRC_f("/mnt/spareA/clion/PRISM/MATLAB/qya.mrc");
		pars.beamsOutput.toMRC_f("/mnt/spareA/clion/PRISM/MATLAB/beamsOutput.mrc");
		cout << "pars.qyInd[2] = " << pars.qyInd[2]<< endl;
		cout << "pars.qxInd[2] = " << pars.qxInd[2]<< endl;
		cout << "pars.qyInd[350] = " << pars.qyInd[350]<< endl;
		cout << "pars.qxInd[350] = " << pars.qxInd[350]<< endl;
        q2.toMRC_f("/Users/ajpryor/Documents/MATLAB/multislice/PRISM/q2.mrc");
        mesh_a.first.toMRC_f("/Users/ajpryor/Documents/MATLAB/multislice/PRISM/xa.mrc");
        mesh_a.second.toMRC_f("/Users/ajpryor/Documents/MATLAB/multislice/PRISM/ya.mrc");

        mask.toMRC_f("/Users/ajpryor/Documents/MATLAB/multislice/PRISM/mask.mrc");
        pars.beams.toMRC_f("/Users/ajpryor/Documents/MATLAB/multislice/PRISM/beams.mrc");

*/
		/*

        //pars.qMask.toMRC_f("/mnt/spareA/clion/PRISM/MATLABdebug.mrc");
		pars.qMask.toMRC_f("/Users/ajpryor/Documents/MATLAB/multislice/PRISM/debug.mrc");
		cout << "pars.qMax = " << pars.qMax << endl;
		cout << "pars.qya.at(1,1) = " << pars.qya.at(1,1) << endl;
		cout << "pars.qya.at(0,1) = " << pars.qya.at(0,1) << endl;
		cout << "q2.at(3,4) = " << q2.at(3,4) << endl;
		cout << "q2.at(5,5) = " << q2.at(5,5) << endl;

		cout << "pars.pixelSize[0] = " << pars.pixelSize[0]<< endl;
		cout << "pars.pixelSize[1] = " << pars.pixelSize[1]<< endl;
		for (auto i = 0; i < 10; ++i){
			cout << "qx[" << i << "] = " << qx[i] << endl;
			cout << "qy[" << i << "] = " << qy[i] << endl;
		}
		cout << "qx[499] = " << qx[499] << endl;
		cout << "qy[499] = " << qy[499] << endl;
		cout << "qx[500] = " << qx[500] << endl;
		cout << "qy[500] = " << qy[500] << endl;
*/
	}
}