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

#include "aberration.h"
#include <array>
#include <algorithm>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iostream>

namespace Prismatic
{

std::string aberrationReadError(size_t line_num, const std::string str)
{
	std::string msg(" \n\nPrismatic: Error getting aberration from");
	std::stringstream ssError;
	msg += " line ";
	ssError << line_num;
	msg += ssError.str();
	msg += ":\n ";
	msg += str;
	msg += " \n";
	return msg;
}

std::vector<aberration> readAberrations(const std::string &filename)
{
	std::vector<aberration> aberrations;
	std::ifstream f(filename);
	if (!f)
		throw std::runtime_error("Unable to open file.\n");
	std::string line;
	std::string token;
	size_t line_num = 1;
	size_t aberration_count = 0;
	if (!std::getline(f, line))
		throw std::runtime_error("Error reading comment line.\n");
	while (std::getline(f, line))
	{
		line = line.substr(line.find_first_not_of(" \n\t"), line.find_last_not_of(" \n\t") + 1);
		if (line.size() <= 3)
		{
			break;
		}
		++aberration_count;
		++line_num;
		int m, n;
        double mag, angle;
		std::stringstream ss;
		ss.precision(8);
		ss << line;
		if (!(ss >> m))
		{
			throw std::domain_error(aberrationReadError(line_num, line));
		}
		if (ss.peek() == ',')
			ss.ignore();
		if (!(ss >> n))
		{
			throw std::domain_error(aberrationReadError(line_num, line));
		}
		if (ss.peek() == ',')
			ss.ignore();
		if (!(ss >> mag))
		{
			throw std::domain_error(aberrationReadError(line_num, line));
		}
		if (ss.peek() == ',')
			ss.ignore();
		if (!(ss >> angle))
		{
			throw std::domain_error(aberrationReadError(line_num, line));
		}
		if (ss.peek() == ',')
			ss.ignore();
		aberrations.emplace_back(aberration{m,n,mag,angle});
	}
	if (aberration_count == 0)
	{
		std::domain_error("Bad input data. No aberrations were found in this file.\n");
	}
	else
	{
		std::cout << "extracted " << aberration_count << " aberrations from " << line_num << " lines in " << filename
				  << std::endl;
	}
	return aberrations;
};

} //namespace Prismatic