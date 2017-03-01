#include "include/fparams.h"
#include "include/ArrayND.h"
#include "include/PRISM01.h"
#include "include/PRISM02.h"
#include "include/PRISM03.h"
#include "include/atom.h"
#include "include/emdSTEM.h"
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include <string>

using namespace std;
int main(int argc, const char** argv) {
	using PRISM_FLOAT_TYPE = double;
	using vec_d = std::vector<PRISM_FLOAT_TYPE>;
	using Array3D = PRISM::ArrayND<3, vec_d>;
	using Array2D = PRISM::ArrayND<2, vec_d>;
	using Array1D = PRISM::ArrayND<1, vec_d>;
	using Array1D_dims = PRISM::ArrayND<1, std::vector<size_t> >;
	PRISM::emdSTEM<PRISM_FLOAT_TYPE> prism_pars;

	std::string filename;
	if (argc>1) {
		filename = std::string(argv[1]);
	} else{
		cout << "PRISM: Correct syntax is prism filename [interpolation factor]" << endl;
		return 0;
	}
	prism_pars.interpolationFactor = (argc>2) ? atoi(argv[2]) : 50;

	PRISM_FLOAT_TYPE one_pixel_size = 100.0 / 1000.0;
	prism_pars.potBound = 1.0;
	prism_pars.numFP = 8.0 / 8.0;
	prism_pars.sliceThickness = 2;
	Array1D_dims cellDim({80, 100, 100}, {3}); // this is z,y,x format
	prism_pars.cellDim = cellDim;
	prism_pars.E0 = 80e3;
	prism_pars.alphaBeamMax = 24 / 1000.0;
	prism_pars.NUM_GPUS = 1;
	prism_pars.NUM_THREADS = 12;

	constexpr double m = 9.109383e-31;
	constexpr double e = 1.602177e-19;
	constexpr double c = 299792458;
	constexpr double h = 6.62607e-34;
	const double pi = std::acos(-1);
	prism_pars.lambda = h / sqrt(2 * m * e * prism_pars.E0) / sqrt(1 + e * prism_pars.E0 / 2 / m / c / c) * 1e10;
	prism_pars.sigma = (2 * pi / prism_pars.lambda / prism_pars.E0) * (m * c * c + e * prism_pars.E0) /
	                   (2 * m * c * c + e * prism_pars.E0);

	PRISM_FLOAT_TYPE f = 4 * prism_pars.interpolationFactor;
	Array1D_dims imageSize({{cellDim[1], cellDim[2]},
	                        {2}});
	std::transform(imageSize.begin(), imageSize.end(), imageSize.begin(),
	               [&f, &prism_pars, &one_pixel_size](size_t &a) {
		               return (size_t) (f * round((PRISM_FLOAT_TYPE) a / one_pixel_size / f));
	               });
	prism_pars.imageSize = imageSize;

	Array1D pixelSize({{(PRISM_FLOAT_TYPE) cellDim[1], (PRISM_FLOAT_TYPE) cellDim[2]},
	                   {2}});
	prism_pars.pixelSize = pixelSize;
	prism_pars.pixelSize[0] /= (PRISM_FLOAT_TYPE)prism_pars.imageSize[0];
	prism_pars.pixelSize[1] /= (PRISM_FLOAT_TYPE)prism_pars.imageSize[1];
	try {
		prism_pars.atoms = PRISM::readAtoms(filename);
	}
	catch (const std::runtime_error &e) {
		cout << "PRISM: Error opening " << filename << endl;
		cout << e.what();
		cout << "Terminating" << endl;
		return 1;
	}
	catch (const std::domain_error &e) {
		cout << "PRISM: Error extracting atomic data from " << filename << "!" << endl;
		cout << e.what();
		cout << "Terminating" << endl;
		return 1;
	}

	Array1D u = PRISM::ones_ND<1, double>({{118}}) * 0.08;
	prism_pars.u = u;
	PRISM::PRISM01(prism_pars);
	PRISM::PRISM02(prism_pars);
	PRISM::PRISM03(prism_pars);

	size_t lower = 13;
	size_t upper = 18;
	Array2D multislice_image;
	multislice_image = PRISM::zeros_ND<2, PRISM_FLOAT_TYPE>({{prism_pars.stack.get_diml(), prism_pars.stack.get_dimk()}});
	for (auto y = 0; y < prism_pars.stack.get_diml(); ++y){
		for (auto x = 0; x < prism_pars.stack.get_dimk(); ++x){
			for (auto b = lower; b < upper; ++b){
				multislice_image.at(y,x) += prism_pars.stack.at(y,x,b,1);
			}
		}
	}
	multislice_image.toMRC_f("multislice_image.mrc");
	cout << "PRISM: Calculation complete. Exiting." << endl;
	return 0;
}
