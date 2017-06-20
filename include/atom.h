// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_ATOM_H
#define PRISM_ATOM_H
#include <string>
#include <vector>
#include <array>
#include <iostream>
struct atom{
	double x;
	double y;
	double z;
	size_t species;
	double sigma;
	void to_string(){
		std::cout << "x = " << x << std::endl;
		std::cout << "y = " << y << std::endl;
		std::cout << "z = " << z << std::endl;
		std::cout << "Debye-Waller thermal displacement standard deviation = " << sigma << std::endl;
		std::cout << "Z = " << species << std::endl;
	}
};

namespace Prismatic {
	std::vector<atom> tileAtoms(const size_t tileX, const size_t tileY, const size_t tileZ, std::vector<atom> atoms);

	std::vector<atom> readAtoms(const std::string& filename);

//	std::array<double, 3> peekDims(const std::string& filename);

	std::array<double, 3> peekDims_xyz(const std::string& filename);

//	std::vector<atom> readAtoms_csv(const std::string& filename);

	std::vector<atom> readAtoms_xyz(const std::string& filename);

    std::string getLowercaseExtension(const std::string filename);

	std::vector<atom> defaultAtoms();

}
#endif //PRISM_ATOM_H