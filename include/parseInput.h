// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_PARSEINPUT_H
#define PRISM_PARSEINPUT_H
#include "defines.h"
#include "meta.h"

namespace Prismatic{
    bool parseInputs(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                            int& argc, const char*** argv);
    bool parseInput(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                            int& argc, const char*** argv);
    void printHelp();
}
#endif //PRISM_PARSEINPUT_H
