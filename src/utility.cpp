// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// Prismatic is distributed under the GNU General Public License (GPL)
// If you use Prismatic, we kindly ask that you cite the following papers:

// 1. Ophus, C.: A fast image simulation algorithm for scanning
//    transmission electron microscopy. Advanced Structural and
//    Chemical Imaging 3(1), 13 (2017)

// 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
//    Implementation of Image Simulation Algorithms for Scanning
//	  Transmission Electron Microscopy. arXiv:1706.08563 (2017)

#include "utility.h"
#include <complex>
#include "defines.h"
#include "configure.h"
namespace Prismatic {



	std::pair<Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> >, Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > >
	upsamplePRISMProbe(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > probe,
	                   const long dimj, const long dimi, long ys, long xs) {
		Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > realspace_probe;
		Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > buffer_probe;
		Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > kspace_probe;

		buffer_probe = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION> >({{(size_t)dimj, (size_t)dimi}});
//		std::cout << "dimj = " << dimj << std::endl;
		long ncy = probe.get_dimj() / 2;
		long ncx = probe.get_dimi() / 2;
		for (auto j = 0; j < probe.get_dimj(); ++j) {
			for (auto i = 0; i < probe.get_dimi(); ++i) {
                buffer_probe.at( (dimj + ((j - ncy + ys) % dimj)) % dimj,
                                 (dimi + ((i - ncx + xs) % dimi)) % dimi) = probe.at(j, i);
//				std::cout << "(dimj + ((j - ncy) % dimj)) % dimj= " << (dimj + ((j - ncy) % dimj)) % dimj<< std::endl;
//				std::cout << "(j - ncy)= " << (j - ncy) << std::endl;
//				std::cout << "(j - ncy) % dimj)= " << (j - ncy) % dimj<< std::endl;

//				buffer_probe.at( (dimj + ((j - ncy) % dimj)) % dimj,
//				                 (dimi + ((i - ncx) % dimi)) % dimi) = probe.at(j, i);
			}
		}
		std::unique_lock<std::mutex> gatekeeper(fftw_plan_lock);
		PRISMATIC_FFTW_PLAN plan = PRISMATIC_FFTW_PLAN_DFT_2D(buffer_probe.get_dimj(), buffer_probe.get_dimi(),
		                                              reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&buffer_probe[0]),
		                                              reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&buffer_probe[0]),
		                                              FFTW_FORWARD, FFTW_ESTIMATE);
		gatekeeper.unlock();
		realspace_probe = buffer_probe;
		PRISMATIC_FFTW_EXECUTE(plan);
		kspace_probe = buffer_probe;
		gatekeeper.lock();
		PRISMATIC_FFTW_DESTROY_PLAN(plan);
		gatekeeper.unlock();
		return std::make_pair(realspace_probe, kspace_probe);
	}

	PRISMATIC_FLOAT_PRECISION computePearsonCorrelation(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > left,
	                                                Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > right){
		PRISMATIC_FLOAT_PRECISION m1, m2, sigma1, sigma2, R;
		m1=m2=sigma1=sigma2=R=0;

		for (auto &i:left) m1 += std::abs(i);
		for (auto &i:right)m2 += std::abs(i);

		m1 /= (left.size());
		m2 /= (right.size());

		for (auto &i:left)sigma1  += pow(std::abs(i)-m1, 2);
		for (auto &i:right)sigma2 += pow(std::abs(i)-m2, 2);

		sigma1 /= (left.size());
		sigma2 /= (right.size());

		sigma1 = sqrt(sigma1);
		sigma2 = sqrt(sigma2);
		for (auto i = 0; i < std::min(left.size(), right.size()); ++i){
			R = R + (std::abs(left[i]) - m1) * (std::abs(right[i]) - m2);
		}
		R/=sqrt(left.size()*right.size());
		return R / (sigma1 * sigma2);
	}
	PRISMATIC_FLOAT_PRECISION computeRfactor(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > left,
	                                     Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > right){
		PRISMATIC_FLOAT_PRECISION accum, diffs;
		accum = diffs = 0;
		for (auto i = 0; i < std::min(left.size(), right.size()); ++i){
			diffs += std::abs(left[i] - right[i]);
			accum += std::abs(left[i]);
		}
		return diffs / accum;
	}



}
