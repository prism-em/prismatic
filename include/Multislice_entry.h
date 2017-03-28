//
// Created by AJ Pryor on 3/5/17.
//

#ifndef PRISM_MULTISLICE_ENTRY_H
#define PRISM_MULTISLICE_ENTRY_H
#include "meta.h"
#include "params.h"
#include "ArrayND.h"
#include "configure.h"
#include "Multislice.h"
#include "PRISM01.h"
#include "PRISM02.h"
#include <algorithm>


namespace PRISM{
	Parameters<PRISM_FLOAT_PRECISION> Multislice_entry(Metadata<PRISM_FLOAT_PRECISION>& meta);
}
#endif //PRISM_MULTISLICE_ENTRY_H
