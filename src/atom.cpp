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

#include "atom.h"
#include <array>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include "kirkland_params.h"

namespace Prismatic
{
std::string atomReadError(size_t line_num, const std::string str)
{
	std::string msg(" \n\nPrismatic: Error getting atomic species from");
	std::stringstream ssError;
	msg += " line ";
	ssError << line_num;
	msg += ssError.str();
	msg += ":\n ";
	msg += str;
	msg += " \n";
	return msg;
}
std::vector<atom> tileAtoms(const size_t tileX, const size_t tileY, const size_t tileZ, std::vector<atom> atoms)
{
	if (tileX == 1 & tileY == 1 & tileZ == 1)
		return atoms; // case where no tiling is necessary
	std::vector<atom> tiled_atoms;
	tiled_atoms.reserve(atoms.size() * tileX * tileY * tileZ);
	for (auto tz = 0; tz < tileZ; ++tz)
	{
		for (auto ty = 0; ty < tileY; ++ty)
		{
			for (auto tx = 0; tx < tileX; ++tx)
			{
				for (auto i = 0; i < atoms.size(); ++i)
				{
					tiled_atoms.emplace_back(atom{(atoms[i].x + tx) / tileX, (atoms[i].y + ty) / tileY, (atoms[i].z + tz) / tileZ, atoms[i].species, atoms[i].sigma, atoms[i].occ});
				}
			}
		}
	}
	return tiled_atoms;
}

void to_xyz(const std::vector<atom> atoms, const std::string filename, const std::string comment, double a, double b, double c)
{
	std::ofstream f(filename, std::ios::out);
	if (f)
	{
		std::stringstream ss;
		ss.precision(8);
		f << comment << '\n';
		ss << '\t';
		ss << a << '\t';
		ss << b << '\t';
		ss << c << '\n';
		f << ss.str();
		for (auto &atom : atoms)
		{
			ss.str("\t");
			ss << atom.species << '\t';
			ss << atom.x * a << '\t';
			ss << atom.y * b << '\t';
			ss << atom.z * c << '\t';
			ss << atom.occ << '\t';
			ss << atom.sigma << '\n';
			f << ss.str();
		}
		f << "-1\n";
	}
}

std::vector<atom> readAtoms_xyz(const std::string &filename)
{
	std::vector<atom> atoms;
	std::ifstream f(filename);
	if (!f)
		throw std::runtime_error("Unable to open file.\n");
	std::string line;
	std::string token;
	size_t line_num = 2;
	size_t atom_count = 0;
	if (!std::getline(f, line))
		throw std::runtime_error("Error reading comment line.\n");
	if (!std::getline(f, line))
		throw std::runtime_error("Error reading unit cell params.\n");
	double a, b, c; // unit cell params
	{
		std::stringstream ss;
		ss.precision(8);
		ss << line;
		if (!(ss >> a) || (a <= 0))
			throw std::domain_error(
				"Bad input data for unit cell dimension a.\n");
		if (!(ss >> b) || (b <= 0))
			throw std::domain_error(
				"Bad input data for unit cell dimension b.\n");
		if (!(ss >> c) || (c <= 0))
			throw std::domain_error(
				"Bad input data for unit cell dimension c.\n");
	}
	while (std::getline(f, line))
	{
		line = line.substr(line.find_first_not_of(" \n\t"), line.find_last_not_of(" \n\t") + 1);
		if (line.size() <= 3)
		{
			break;
		}
		++atom_count;
		++line_num;
		double tx, ty, tz, occ, sigma;
		size_t tspecies;
		std::stringstream ss;
		ss.precision(8);
		ss << line;
		if (!(ss >> tspecies) || (tspecies > NUM_SPECIES_KIRKLAND))
		{
			throw std::domain_error(atomReadError(line_num, line));
		}
		if (ss.peek() == ',')
			ss.ignore();
		if (!(ss >> tx))
		{
			throw std::domain_error(atomReadError(line_num, line));
		}
		if (ss.peek() == ',')
			ss.ignore();
		if (!(ss >> ty))
		{
			throw std::domain_error(atomReadError(line_num, line));
		}
		if (ss.peek() == ',')
			ss.ignore();
		if (!(ss >> tz))
		{
			throw std::domain_error(atomReadError(line_num, line));
		}
		if (ss.peek() == ',')
			ss.ignore();
		if (!(ss >> occ))
		{
			throw std::domain_error(atomReadError(line_num, line));
		}
		if (ss.peek() == ',')
			ss.ignore();
		if (!(ss >> sigma))
		{
			throw std::domain_error(atomReadError(line_num, line));
		}
		if (ss.peek() == ',')
			ss.ignore();
		atoms.emplace_back(atom{tx / a, ty / b, tz / c, tspecies, sigma, occ});
		//				atoms.emplace_back(atom{tx / a, ty / b , tz / c, tspecies});
	}
	if (atom_count == 0)
	{
		std::domain_error("Bad input data. No atoms were found in this file.\n");
	}
	else
	{
		std::cout << "extracted " << atom_count << " atoms from " << line_num << " lines in " << filename
				  << std::endl;
	}
	return atoms;
};

std::array<double, 3> peekDims_xyz(const std::string &filename)
{
	std::ifstream f(filename);
	if (!f)
		throw std::runtime_error("Unable to open file.\n");
	std::string line;
	std::string token;
	if (!std::getline(f, line))
		throw std::runtime_error("Error reading comment line.\n");
	if (!std::getline(f, line))
		throw std::runtime_error("Error reading unit cell params.\n");
	double a, b, c; // unit cell params
	{
		std::stringstream ss;
		ss.precision(8);
		ss << line;
		if (!(ss >> a))
			throw std::domain_error("Bad input data for unit cell dimension a.\n");
		if (!(ss >> b))
			throw std::domain_error("Bad input data for unit cell dimension b.\n");
		if (!(ss >> c))
			throw std::domain_error("Bad input data for unit cell dimension c.\n");
	}
	//		return {a,b,c};
	return {c, b, a};
}

std::string getLowercaseExtension(const std::string filename)
{
	std::string::size_type idx;
	idx = filename.rfind('.');
	if (idx != std::string::npos)
	{
		std::string ext = filename.substr(idx + 1);
		std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
		return ext;
	}
	else
	{
		return "";
	}
}

std::vector<atom> defaultAtoms()
{
	// returns the unit cell of 100 Silicon from the file SI100.XYZ. This is sometimes used as a default input in case
	// the user hasn't provided one, for example in the GUI

	//		one unit cell of 100 silicon
	//		5.43    5.43    5.43
	//		14  0.0000  0.0000  0.0000  1.0  0.076
	//		14  2.7150  2.7150  0.0000  1.0  0.076
	//		14  1.3575  4.0725  1.3575  1.0  0.076
	//		14  4.0725  1.3575  1.3575  1.0  0.076
	//		14  2.7150  0.0000  2.7150  1.0  0.076
	//		14  0.0000  2.7150  2.7150  1.0  0.076
	//		14  1.3575  1.3575  4.0725  1.0  0.076
	//		14  4.0725  4.0725  4.0725  1.0  0.076
	//		-1

	std::vector<atom> result;
	result.resize(8);
	result.emplace_back(atom{0.0000 / 5.43, 0.0000 / 5.43, 0.0000 / 5.43, 14, 0.076});
	result.emplace_back(atom{2.7150 / 5.43, 2.7150 / 5.43, 0.0000 / 5.43, 14, 0.076});
	result.emplace_back(atom{1.3575 / 5.43, 4.0725 / 5.43, 1.3575 / 5.43, 14, 0.076});
	result.emplace_back(atom{4.0725 / 5.43, 1.3575 / 5.43, 1.3575 / 5.43, 14, 0.076});
	result.emplace_back(atom{2.7150 / 5.43, 0.0000 / 5.43, 2.7150 / 5.43, 14, 0.076});
	result.emplace_back(atom{0.0000 / 5.43, 2.7150 / 5.43, 2.7150 / 5.43, 14, 0.076});
	result.emplace_back(atom{1.3575 / 5.43, 1.3575 / 5.43, 4.0725 / 5.43, 14, 0.076});
	result.emplace_back(atom{4.0725 / 5.43, 4.0725 / 5.43, 4.0725 / 5.43, 14, 0.076});
	return result;
}
} // namespace Prismatic
