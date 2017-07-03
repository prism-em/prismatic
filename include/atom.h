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
    double occ;
    void to_string() const{
        std::cout << "x = " << x << std::endl;
        std::cout << "y = " << y << std::endl;
        std::cout << "z = " << z << std::endl;
        std::cout << "Debye-Waller thermal displacement standard deviation = " << sigma << std::endl;
        std::cout << "Z = " << species << std::endl;
        std::cout << "Occupancy = " << occ << std::endl;
    }
};

namespace Prismatic {

    void to_xyz(const std::vector<atom> atoms, const std::string filename, const std::string comment, double a, double b, double c);

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