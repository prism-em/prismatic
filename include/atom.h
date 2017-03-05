//
// Created by AJ Pryor on 2/21/17.
//

#ifndef PRISM_ATOM_H
#define PRISM_ATOM_H
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <stdexcept>

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

namespace PRISM{
	inline std::vector<atom> readAtoms(const std::string& filename){
		std::vector<atom> atoms;
		std::ifstream f(filename);
		if (!f)throw std::runtime_error("Unable to open file.\n");
		std::string line;
		std::string token;
		size_t line_num = 0;
		size_t atom_count = 0;
		while (std::getline(f,line)){
			++atom_count;
			++line_num;
			double tx, ty, tz;
			size_t tspecies;
			std::stringstream ss(line);
			if(!(ss >> tx))throw std::domain_error("Bad input data for X. The txt file should continue 4 comma separated values per line (x,y,z,species).\n");
			if(ss.peek()==',')ss.ignore();
			if(!(ss >> ty))throw std::domain_error("Bad input data for Y. The txt file should continue 4 comma separated values per line (x,y,z,species).\n");
			if(ss.peek()==',')ss.ignore();
			if(!(ss >> tz))throw std::domain_error("Bad input data for Z. The txt file should continue 4 comma separated values per line (x,y,z,species).\n");
			if(ss.peek()==',')ss.ignore();
			if(!(ss >> tspecies))throw std::domain_error("Bad input data for atomic species. The txt file should continue 4 comma separated values per line (x,y,z,species).\n");
			if(ss.peek()==',')ss.ignore();
			atoms.emplace_back(atom{tx,ty,tz,tspecies});
		}
		if (atom_count == 0){
			std::domain_error("Bad input data. No atoms were found in this file.\n");
		} else {
			std::cout << "extracted " << atom_count << " atoms from " << line_num << " lines in " << filename
			          << std::endl;
		}
			return atoms;
	};

}
#endif //PRISM_ATOM_H
