#ifndef PRISMATIC_FILEIO_H
#define PRISMATIC_FILEIO_H
#include "H5Cpp.h"
#include "params.h"

struct complex_float_t
{
	PRISMATIC_FLOAT_PRECISION re;
	PRISMATIC_FLOAT_PRECISION im;
};


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

void setupSMatrixOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const int FP, const float dummy);

void setupSMatrixOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const int FP, const double dummy);

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

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, PRISMATIC_FLOAT_PRECISION *val);

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, int &val);

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, std::string &val);

void writeComplexDataset(H5::Group group, const std::string &dsetname, const std::complex<float> *buffer, const hsize_t *mdims, const size_t &rank);

void writeComplexDataset(H5::Group group, const std::string &dsetname, const std::complex<double> *buffer, const hsize_t *mdims, const size_t &rank);

template <size_t N>
void readComplexDataset(ArrayND<N, std::vector<std::complex<PRISMATIC_FLOAT_PRECISION>>> &output, const std::string &filename, const std::string &dataPath)
{
	H5::H5File input = H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	H5::DataSet dataset = input.openDataSet(dataPath.c_str());
	H5::DataSpace dataspace = dataset.getSpace();
    H5::CompType complex_type = H5::CompType(sizeof(complex_float_t));
	const H5std_string re_str("r"); //using h5py default configuration
	const H5std_string im_str("i");

	if(sizeof(PRISMATIC_FLOAT_PRECISION) == 4)
	{
		complex_type.insertMember(re_str, 0, H5::PredType::NATIVE_FLOAT);
		complex_type.insertMember(im_str, 4, H5::PredType::NATIVE_FLOAT);
	}
	else{
		complex_type.insertMember(re_str, 0, H5::PredType::NATIVE_DOUBLE);
		complex_type.insertMember(im_str, 4, H5::PredType::NATIVE_DOUBLE);
	}
	

	hsize_t dims_out[N];
	int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
	H5::DataSpace mspace(N,dims_out);

    std::array<size_t, N> data_dims;
    size_t totalSize = 1;
	for(auto i = 0; i < N; i++) totalSize *= dims_out[i];
    for(auto i = 0; i < N; i++) data_dims[N-1-i] = dims_out[i];

    std::complex<PRISMATIC_FLOAT_PRECISION> data_in[totalSize];
    dataset.read(data_in, complex_type, mspace, dataspace);

    output = zeros_ND<N, std::complex<PRISMATIC_FLOAT_PRECISION>>(data_dims);
    for(auto i = 0; i < output.size(); i++) output[i] = data_in[i];

    mspace.close();
    dataspace.close();
    dataset.close();
    input.close();
};

int countDataGroups(H5::Group group, const std::string &basename);

void configureSupergroup(Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
						const std::string &sgName,
						const std::vector<std::vector<PRISMATIC_FLOAT_PRECISION>> &dims,
						const std::vector<std::string> &dims_name,
						const std::vector<std::string> &dims_units,
						const std::vector<std::vector<PRISMATIC_FLOAT_PRECISION>> &sgdims,
						const std::vector<std::string> &sgdims_name,
						const std::vector<std::string> &sgdims_units);

void writeVirtualDataSet(H5::Group group,
						const std::string &dsetName,
						std::vector<H5::DataSet> &datasets,
						std::vector<std::vector<size_t>> indices);

void depthSeriesSG(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);

std::string getDatasetName(H5::DataSet &dataset);
						
} //namespace Prismatic

#endif //PRISMATIC_FILEIO_H