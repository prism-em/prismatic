#ifndef PRISMATIC_FILEIO_H
#define PRISMATIC_FILEIO_H
#include "H5Cpp.h"
#include "params.h"

namespace Prismatic{

void setupOutputFile(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> prismatic_pars);

void setup4DOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const float dummy);

void setup4DOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const double dummy);

void setupVDOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const float dummy);

void setupVDOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const double dummy);

void setup2DOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const float dummy);

void setup2DOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const double dummy);

void setupDPCOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const float dummy);

void setupDPCOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const double dummy);

void writeRealSlice(H5::DataSet dataset, const float *buffer, const hsize_t *mdims);

void writeRealSlice(H5::DataSet dataset, const double *buffer, const hsize_t *mdims);

void writeDatacube3D(H5::DataSet dataset, const float *buffer, const hsize_t *mdims);

void writeDatacube3D(H5::DataSet dataset, const double *buffer, const hsize_t *mdims);

void writeDatacube4D(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, float *buffer, const hsize_t *mdims, const hsize_t *offset, const float numFP, const std::string nameString);

void writeDatacube4D(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, double *buffer, const hsize_t *mdims, const hsize_t *offset, const double numFP, const std::string nameString);

void writeStringArray(H5::DataSet dataset,H5std_string * string_array, hsize_t elements);

std::string getDigitString(int digit);

void writeMetadata(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, float dummy);

void writeMetadata(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, double dummy);

} //namespace Prismatic

#endif //PRISMATIC_FILEIO_H