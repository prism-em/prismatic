//
// Created by Alan Pryor on 4/12/17.
//

#ifndef PRISM_PARSEINPUT_H
#define PRISM_PARSEINPUT_H
#include "defines.h"
#include "meta.h"

namespace PRISM{
    bool parseInputs(Metadata<PRISM_FLOAT_PRECISION>& meta,
                            int& argc, const char*** argv);
    bool parseInput(Metadata<PRISM_FLOAT_PRECISION>& meta,
                            int& argc, const char*** argv);
    void printHelp();
}
#endif //PRISM_PARSEINPUT_H
