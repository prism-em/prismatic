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
#include "utility.h"
#include <vector>

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
        PRISMATIC_FLOAT_PRECISION mag, angle;
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

Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> getChi(Array2D<PRISMATIC_FLOAT_PRECISION> &q,
                                                        Array2D<PRISMATIC_FLOAT_PRECISION> &qTheta,
                                                        PRISMATIC_FLOAT_PRECISION &lambda, 
                                                        std::vector<aberration> &ab)
{

    Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> chi = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>({{q.get_dimj(), q.get_dimi()}});
    const PRISMATIC_FLOAT_PRECISION pi = acos(-1);
    for(auto n = 0; n < ab.size(); n++)
    {
        for(auto j = 0; j < chi.get_dimj(); j++)
        {
            for(auto i = 0; i < chi.get_dimi(); i++)
            {
				PRISMATIC_FLOAT_PRECISION rad = ab[n].angle * pi / 180.0;
                PRISMATIC_FLOAT_PRECISION cx = ab[n].mag * cos(ab[n].n*rad);
                PRISMATIC_FLOAT_PRECISION cy = ab[n].mag * sin(ab[n].n*rad);
                PRISMATIC_FLOAT_PRECISION tmp = chi.at(j,i).real();
				tmp += cx*pow(lambda*q.at(j,i), ab[n].m)*cos(ab[n].n * qTheta.at(j,i));
                tmp += cy*pow(lambda*q.at(j,i), ab[n].m)*sin(ab[n].n * qTheta.at(j,i));
                chi.at(j,i).real(tmp);
            }
        }
    }

    return chi;

};

std::vector<aberration> updateAberrations(std::vector<aberration> ab, 
                                        PRISMATIC_FLOAT_PRECISION C1, 
                                        PRISMATIC_FLOAT_PRECISION C3, 
                                        PRISMATIC_FLOAT_PRECISION C5,
                                        const PRISMATIC_FLOAT_PRECISION &lambda)
{
    //prune for unique aberrations
    if(ab.size() > 0)
    {
        std::sort(ab.begin(), ab.end(), 
                    [](aberration &a, aberration &b){
                        return std::make_pair(a.m, a.n) < std::make_pair(b.m, b.n);
                    });
        ab = getUnique(ab);

        //m < n and m+n % 2 == 1 aren't valid components of basis set
        std::vector<aberration> tmp;
        for(auto i = 0; i < ab.size(); i++)
        {
            bool check = (ab[i].m  >= ab[i].n) && !((ab[i].m + ab[i].n) % 2);
            if(check)
            {
                tmp.push_back(ab[i]);
            }
        }

        ab = tmp;
    }

    //input C1, C3, C5 aassumed to all be in angstrom
    //first, check if aberration exists
    //if exists and CX is not nan, update; else, don't update
    //if does not exist and CX is not nan, set to CX; else if CX is nan, do nothing
    PRISMATIC_FLOAT_PRECISION pi = acos(-1);

    bool C1_exists = false;
    for(auto i = 0; i < ab.size(); i++)
    {
        if(ab[i].m ==0 and ab[i].n == 2)
        {
            C1_exists = true;
            if(not std::isnan(C1)) ab[i].mag = C1 * pi / lambda;
        }
    }

    if((not C1_exists) and (not std::isnan(C1)))
    {
        PRISMATIC_FLOAT_PRECISION mag = C1 * pi / lambda;
        aberration new_C1 = aberration{2, 0, mag, 0.0};
        ab.push_back(new_C1);
    }

    bool C3_exists = false;
    for(auto i = 0; i < ab.size(); i++)
    {
        if(ab[i].m ==0 and ab[i].n == 4)
        {
            C3_exists = true;
            if(not std::isnan(C3)) ab[i].mag = C3 * pi / (2.0*lambda);
        }
    }

    if((not C3_exists) and (not std::isnan(C3)))
    {
        PRISMATIC_FLOAT_PRECISION mag = C3 * pi / (2.0*lambda);
        aberration new_C3 = aberration{4, 0, mag, 0.0};
        ab.push_back(new_C3);
    }

    bool C5_exists = false;
    for(auto i = 0; i < ab.size(); i++)
    {
        if(ab[i].m ==0 and ab[i].n == 6)
        {
            C5_exists = true;
            if(not std::isnan(C5)) ab[i].mag = C5 * pi / (3.0*lambda);
        }
    }

    if((not C5_exists) and (not std::isnan(C5)))
    {
        PRISMATIC_FLOAT_PRECISION mag = C5 * pi / (3.0*lambda);
        aberration new_C5 = aberration{6, 0, mag, 0.0};
        ab.push_back(new_C5);
    }

    return ab;
};

} //namespace Prismatic