// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_ATOM_H
#define PRISM_ATOM_H
#include <string>
#include <vector>
//#include <stdlib.h>
//#include <sstream>
//#include <fstream>
//#include <stdexcept>
#include <iostream>
struct atom{
	double x,y,z;
	size_t species;
	void to_string(){
		std::cout << "x = " << x << std::endl;
		std::cout << "y = " << y << std::endl;
		std::cout << "z = " << z << std::endl;
		std::cout << "Z = " << species << std::endl;
	}
};

namespace PRISM {
	std::vector<atom> tileAtoms(const size_t tileX, const size_t tileY, const size_t tileZ, std::vector<atom> atoms);

	std::vector<atom> readAtoms(const std::string& filename);

	std::vector<atom> readAtoms_csv(const std::string& filename);

	std::vector<atom> readAtoms_XYZ(const std::string& filename);

}
#endif //PRISM_ATOM_H
