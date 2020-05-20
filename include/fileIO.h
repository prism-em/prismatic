#ifndef PRISMATIC_FILEIO_H
#define PRISMATIC_FILEIO_H
#include "H5Cpp.h"
#include "params.h"

namespace Prismatic{

void setupOutputFile(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);

void setup4DOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers, const float dummy);

void setup4DOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers, const double dummy);

void setupVDOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers, const float dummy);

void setupVDOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers, const double dummy);

void setup2DOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers, const float dummy);

void setup2DOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers, const double dummy);

void setupDPCOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers, const float dummy);

void setupDPCOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers, const double dummy);

void writeRealSlice(H5::DataSet dataset, const float *buffer, const hsize_t *mdims);

void writeRealSlice(H5::DataSet dataset, const double *buffer, const hsize_t *mdims);

void writeDatacube3D(H5::DataSet dataset, const float *buffer, const hsize_t *mdims);

void writeDatacube3D(H5::DataSet dataset, const double *buffer, const hsize_t *mdims);

void writeDatacube4D(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, float *buffer, const hsize_t *mdims, const hsize_t *offset, const float numFP, const std::string nameString);

void writeDatacube4D(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, double *buffer, const hsize_t *mdims, const hsize_t *offset, const double numFP, const std::string nameString);

void writeStringArray(H5::DataSet dataset,H5std_string * string_array, hsize_t elements);

void savePotentialSlices(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);

std::string getDigitString(int digit);

void writeMetadata(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, float dummy);

void writeMetadata(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, double dummy);

Array2D<PRISMATIC_FLOAT_PRECISION> readDataset2D(const std::string &filename, const std::string &dataPath);

Array3D<PRISMATIC_FLOAT_PRECISION> readDataset3D(const std::string &filename, const std::string &dataPath);

Array4D<PRISMATIC_FLOAT_PRECISION> readDataset4D(const std::string &filename, const std::string &dataPath);

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, PRISMATIC_FLOAT_PRECISION &val);

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, int &val);

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, std::string &val);

} //namespace Prismatic

#endif //PRISMATIC_FILEIO_H