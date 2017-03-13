//
// Created by AJ Pryor on 3/13/17.
//

#ifndef PRISM_GETWORKID_H
#define PRISM_GETWORKID_H
#include "params.h"
#include "configure.h"
//#include <mutex>
//template <class T>
//class PRISM::Parameters;
// helper function for dispatching work
bool getWorkID(const PRISM::Parameters<PRISM_FLOAT_PRECISION>& pars, long long& Nstart, long long& Nstop);
#endif //PRISM_GETWORKID_H
