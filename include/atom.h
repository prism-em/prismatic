//
// Created by AJ Pryor on 2/21/17.
//

#ifndef PRISM_ATOM_H
#define PRISM_ATOM_H
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>

struct atom{
	double x,y,z;
	size_t species;
};

namespace PRISM{
	std::vector<atom> readAtoms(const std::string& filename){
		std::vector<atom> atoms;
		std::ifstream f(filename);
		std::string line;
		std::cout << "opening file " << filename << std::endl;
		while (std::getline(f,line)){
			double tx, ty, tz;
			size_t tspecies;
			std::stringstream ss;
			ss << line;
			std::cout <<  ss.str() << std::endl;
		}
		atom a{1,2,3,4};
		atoms.push_back(a);
		return atoms;
	};

}
#endif //PRISM_ATOM_H
