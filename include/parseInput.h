//
// Created by Alan Pryor on 4/12/17.
//

#ifndef PRISM_PARSEINPUT_H
#define PRISM_PARSEINPUT_H
#include "defines.h"
#include "meta.h"

namespace PRISM{
    enum class ParseResult{Success, Failure, Help};
    bool parseInputs(Metadata<PRISM_FLOAT_PRECISION>& meta,
                            int& argc, const char*** argv);
    ParseResult parseInput(Metadata<PRISM_FLOAT_PRECISION>& meta,
                            int& argc, const char*** argv);
    void printHelp();
}
#endif //PRISM_PARSEINPUT_H
