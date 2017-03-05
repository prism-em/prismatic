//
// Created by AJ Pryor on 3/5/17.
//

#ifndef PRISM_MULTISLICE_ENTRY_H
#define PRISM_MULTISLICE_ENTRY_H
#include "meta.h"
#include "params.h"
#include "ArrayND.h"
#include <algorithm>

namespace PRISM{
	template <class T>
	int Multislice_entry(Metadata<T>& meta){

		using PRISM_FLOAT_TYPE = double;
		using vec_d = std::vector<PRISM_FLOAT_TYPE>;
		using Array3D = ArrayND<3, vec_d>;
		using Array2D = ArrayND<2, vec_d>;
		using Array1D = ArrayND<1, vec_d>;
		using Array1D_dims = ArrayND<1, std::vector<size_t> >;


		Parameters<PRISM_FLOAT_TYPE> prism_pars;
		prism_pars.meta = meta;
		std::cout<<"Dummy code for Multislice entrypoint" << std::endl;
//		constexpr double m = 9.109383e-31;
//		constexpr double e = 1.602177e-19;
//		constexpr double c = 299792458;
//		constexpr double h = 6.62607e-34;
//		const double pi = std::acos(-1);
//		prism_pars.lambda = h / sqrt(2 * m * e * prism_pars.meta.E0) / sqrt(1 + e * prism_pars.meta.E0 / 2 / m / c / c) * 1e10;
//		prism_pars.sigma = (2 * pi / prism_pars.lambda / prism_pars.meta.E0) * (m * c * c + e * prism_pars.meta.E0) /
//		                   (2 * m * c * c + e * prism_pars.meta.E0);
//
//		PRISM_FLOAT_TYPE f = 4 * prism_pars.meta.interpolationFactor;
//		Array1D_dims imageSize({{meta.cellDim[1], meta.cellDim[2]}}, {{2}});
//		std::transform(imageSize.begin(), imageSize.end(), imageSize.begin(),
//		               [&f, &prism_pars](size_t &a) {
//			               return (size_t) (f * round((PRISM_FLOAT_TYPE) a / prism_pars.meta.realspace_pixelSize / f));
//		               });
//		prism_pars.imageSize = imageSize;
//
//		Array1D pixelSize({{(PRISM_FLOAT_TYPE) meta.cellDim[1], (PRISM_FLOAT_TYPE) meta.cellDim[2]}}, {{2}});
//		prism_pars.pixelSize = pixelSize;
//		prism_pars.pixelSize[0] /= (PRISM_FLOAT_TYPE)prism_pars.imageSize[0];
//		prism_pars.pixelSize[1] /= (PRISM_FLOAT_TYPE)prism_pars.imageSize[1];
//		try {
//			prism_pars.atoms = readAtoms(prism_pars.meta.filename_atoms);
//		}
//		catch (const std::runtime_error &e) {
//			std::cout << "PRISM: Error opening " << prism_pars.meta.filename_atoms << std::endl;
//			std::cout << e.what();
//			std::cout << "Terminating" << std::endl;
//			return -1;
//		}
//		catch (const std::domain_error &e) {
//			std::cout << "PRISM: Error extracting atomic data from " << prism_pars.meta.filename_atoms << "!" << std::endl;
//			std::cout << e.what();
//			std::cout << "Terminating" << std::endl;
//			return -2;
//		}
//
//		Array1D u = ones_ND<1, double>({{118}}) * 0.08;
//		prism_pars.u = u;
//		PRISM01(prism_pars);
//
//
//		size_t lower = 13;
//		size_t upper = 18;
//		Array2D prism_image;
//		prism_image = zeros_ND<2, PRISM_FLOAT_TYPE>({{prism_pars.stack.get_diml(), prism_pars.stack.get_dimk()}});
//		for (auto y = 0; y < prism_pars.stack.get_diml(); ++y){
//			for (auto x = 0; x < prism_pars.stack.get_dimk(); ++x){
//				for (auto b = lower; b < upper; ++b){
//					prism_image.at(y,x) += prism_pars.stack.at(y,x,b,1);
//				}
//			}
//		}
//
//		prism_image.toMRC_f(prism_pars.meta.filename_output.c_str());
//		std::cout << "Calculation complete.\n" << std::endl;
		return 0;
	}

}
#endif //PRISM_MULTISLICE_ENTRY_H
