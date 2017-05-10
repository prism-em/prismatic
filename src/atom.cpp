// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include "atom.h"
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iostream>
//
//struct atom{
//	double x,y,z;
//	size_t species;
//	void to_string(){
//		std::cout << "x = " << x << std::endl;
//		std::cout << "y = " << y << std::endl;
//		std::cout << "z = " << z << std::endl;
//		std::cout << "Z = " << species << std::endl;
//	}
//};

namespace PRISM {
	std::vector<atom> tileAtoms(const size_t tileX, const size_t tileY, const size_t tileZ, std::vector<atom> atoms) {
		if (tileX == 1 & tileY == 1 & tileZ == 1)return atoms; // case where no tiling is necessary
		std::vector<atom> tiled_atoms;
		tiled_atoms.reserve(atoms.size() * tileX * tileY * tileZ);
		for (auto tz = 0; tz < tileZ; ++tz) {
			for (auto ty = 0; ty < tileY; ++ty) {
				for (auto tx = 0; tx < tileX; ++tx) {
					for (auto i = 0; i < atoms.size(); ++i){
						tiled_atoms.emplace_back(atom{(atoms[i].x + tx) / tileX, (atoms[i].y + ty) / tileY, (atoms[i].z + tz) / tileZ, atoms[i].species});
					}
				}
			}
		}
		return tiled_atoms;
	}

	std::vector<atom> readAtoms_csv(const std::string& filename){
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

	std::vector<atom> readAtoms_xyz(const std::string& filename){
		std::vector<atom> atoms;
		std::ifstream f(filename);
		if (!f)throw std::runtime_error("Unable to open file.\n");
		std::string line;
		std::string token;
		size_t line_num = 0;
		size_t atom_count = 0;
		if (!std::getline(f,line)) throw std::runtime_error("Error reading comment line.\n");
		if (!std::getline(f,line)) throw std::runtime_error("Error reading unit cell params.\n");
        double a,b,c; // unit cell params
		{
			std::stringstream ss(line);
			if (!(ss >> a))
				throw std::domain_error(
						"Bad input data for unit cell dimension a.\n");
			if (!(ss >> b))
				throw std::domain_error(
						"Bad input data for unit cell dimension b.\n");
			if (!(ss >> c))
				throw std::domain_error(
						"Bad input data for unit cell dimension c.\n");
		}
        std::cout << "Unit cell a, b, c = " << a << ", " << b << ", " << c << std::endl;
			while (std::getline(f, line)) {
                line = line.substr(line.find_first_not_of(" \n\t"), line.find_last_not_of(" \n\t"));
				if (line.size() <=3){
					break;
				}
                ++atom_count;
                ++line_num;
//                std::cout << "atom_count = " << atom_count << " line = " << line <<  std::endl;
				double tx, ty, tz;
				size_t tspecies;
				std::stringstream ss(line);
				if (!(ss >> tspecies))
					throw std::domain_error(
							"Bad input data for atomic species. The txt file should continue 4 comma separated values per line (x,y,z,species).\n");
				if (ss.peek() == ',')ss.ignore();
				std::cout << ss.str() << std::endl;
				if (!(ss >> tx))
					throw std::domain_error(
							"Bad input data for X. The txt file should continue 4 comma separated values per line (x,y,z,species).\n");
				if (ss.peek() == ',')ss.ignore();
				if (!(ss >> ty))
					throw std::domain_error(
							"Bad input data for Y. The txt file should continue 4 comma separated values per line (x,y,z,species).\n");
				if (ss.peek() == ',')ss.ignore();
				if (!(ss >> tz))
					throw std::domain_error(
							"Bad input data for Z. The txt file should continue 4 comma separated values per line (x,y,z,species).\n");
				if (ss.peek() == ',')ss.ignore();

				atoms.emplace_back(atom{tx / a, ty / b , tz / c, tspecies});
			}
			if (atom_count == 0) {
				std::domain_error("Bad input data. No atoms were found in this file.\n");
			} else {
				std::cout << "extracted " << atom_count << " atoms from " << line_num << " lines in " << filename
						  << std::endl;
			}

		return atoms;
	};

	std::vector<atom> readAtoms(const std::string& filename){
//		constexpr std::string supported_filetypes[] = {".csv",".xyz"};
		std::string::size_type idx;
		idx = filename.rfind('.');
		if(idx != std::string::npos) {
			std::string ext = filename.substr(idx+1);
			std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
			std::cout << "extension = " << ext << std::endl;
			if (ext == "csv"){
				return readAtoms_csv(filename);
			}
			if (ext == "xyz"){
				return readAtoms_xyz(filename);
			}
			throw std::domain_error("Unsupported file extension. Make sure the input data is of a supported filetype and named accordingly.\n");
		} else {
			throw std::domain_error("Unable to determine file extension. Make sure an extension exists.\n");
		}

	}

}
