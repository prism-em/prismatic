// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include "utility.h"
#include <complex>
#include "defines.h"
#include "configure.h"
namespace PRISM {
	std::pair<PRISM::Array2D<std::complex<PRISM_FLOAT_PRECISION> >, PRISM::Array2D<std::complex<PRISM_FLOAT_PRECISION> > >
	upsamplePRISMProbe(PRISM::Array2D<std::complex<PRISM_FLOAT_PRECISION> > probe, const size_t dimj, const size_t dimi) {
		Array2D<std::complex<PRISM_FLOAT_PRECISION> > realspace_probe;
		Array2D<std::complex<PRISM_FLOAT_PRECISION> > buffer_probe;
		Array2D<std::complex<PRISM_FLOAT_PRECISION> > kspace_probe;

		buffer_probe = zeros_ND<2, std::complex<PRISM_FLOAT_PRECISION> >({{dimj, dimi}});

		for (auto j = 0; j < buffer_probe.get_dimj(); ++j) {
			for (auto i = 0; i < buffer_probe.get_dimi(); ++i) {
				buffer_probe.at(j, i) = probe.at(j, i);
			}
		}
		std::unique_lock<std::mutex> gatekeeper(fftw_plan_lock);
		PRISM_FFTW_PLAN plan = PRISM_FFTW_PLAN_DFT_2D(buffer_probe.get_dimj(), buffer_probe.get_dimi(),
		                                              reinterpret_cast<PRISM_FFTW_COMPLEX *>(&buffer_probe[0]),
		                                              reinterpret_cast<PRISM_FFTW_COMPLEX *>(&buffer_probe[0]),
		                                              FFTW_FORWARD, FFTW_ESTIMATE);
		gatekeeper.unlock();
		realspace_probe = buffer_probe;
		PRISM_FFTW_EXECUTE(plan);
		kspace_probe = buffer_probe;
		gatekeeper.lock();
		PRISM_FFTW_DESTROY_PLAN(plan);
		gatekeeper.unlock();
		return std::make_pair(realspace_probe, kspace_probe);
	}
}