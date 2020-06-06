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

static std::mutex write4D_lock;

void setupOutputFile(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);

void setup4DOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers);

void setupVDOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers);

void setup2DOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers);

void setupDPCOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const size_t numLayers);

void setupSMatrixOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const int FP);

void writeRealSlice(H5::DataSet dataset, const PRISMATIC_FLOAT_PRECISION *buffer, const hsize_t *mdims);

void writeDatacube3D(H5::DataSet dataset, const PRISMATIC_FLOAT_PRECISION *buffer, const hsize_t *mdims);

void writeDatacube3D(H5::DataSet dataset, const std::complex<PRISMATIC_FLOAT_PRECISION> *buffer, const hsize_t *mdims);

// void writeDatacube4D(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, PRISMATIC_FLOAT_PRECISION *buffer, const hsize_t *mdims, const hsize_t *offset, const PRISMATIC_FLOAT_PRECISION numFP, const std::string nameString);
template<class T>
void writeDatacube4D(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, T *buffer, const hsize_t *mdims, const hsize_t *offset, const PRISMATIC_FLOAT_PRECISION numFP, const std::string nameString)
{
	//for 4D writes, need to first read the data set and then add; this way, FP are accounted for
	//lock the whole file access/writing procedure in only one location
	std::unique_lock<std::mutex> writeGatekeeper(write4D_lock);

    H5::Group dataGroup = pars.outputFile.openGroup(nameString);
    H5::DataSet dataset = dataGroup.openDataSet("datacube");

    //set up file and memory spaces
    H5::DataSpace fspace = dataset.getSpace();

    H5::DataSpace mspace(4, mdims); //rank = 4

    fspace.selectHyperslab(H5S_SELECT_SET, mdims, offset);

    //divide by num FP
    for (auto i = 0; i < mdims[0] * mdims[1] * mdims[2] * mdims[3]; i++)
        buffer[i] /= numFP;

    //restride the dataset so that qx and qy are flipped
    T *finalBuffer = (T *)malloc(mdims[0] * mdims[1] * mdims[2] * mdims[3] * sizeof(T));
    for (auto i = 0; i < mdims[2]; i++)
    {
        for (auto j = 0; j < mdims[3]; j++)
        {
            finalBuffer[i * mdims[3] + j] = buffer[j * mdims[2] + i];
        }
    }

    //add frozen phonon set
    T *readBuffer = (T *)malloc(mdims[0] * mdims[1] * mdims[2] * mdims[3] * sizeof(T));
    dataset.read(&readBuffer[0], dataset.getDataType(), mspace, fspace);
    for (auto i = 0; i < mdims[0] * mdims[1] * mdims[2] * mdims[3]; i++)
        finalBuffer[i] += readBuffer[i];
    free(readBuffer);
		
    dataset.write(finalBuffer, dataset.getDataType(), mspace, fspace);
    free(finalBuffer);
    fspace.close();
    mspace.close();
    dataset.flush(H5F_SCOPE_LOCAL);
    dataset.close();
    dataGroup.flush(H5F_SCOPE_LOCAL);
    dataGroup.close();
    pars.outputFile.flush(H5F_SCOPE_LOCAL);

	writeGatekeeper.unlock();
};

void writeStringArray(H5::DataSet dataset,H5std_string * string_array, hsize_t elements);

void savePotentialSlices(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);

std::string getDigitString(int digit);

void writeMetadata(Parameters<PRISMATIC_FLOAT_PRECISION> &pars);

Array2D<PRISMATIC_FLOAT_PRECISION> readDataSet2D(const std::string &filename, const std::string &dataPath);

Array3D<PRISMATIC_FLOAT_PRECISION> readDataSet3D(const std::string &filename, const std::string &dataPath);

Array4D<PRISMATIC_FLOAT_PRECISION> readDataSet4D(const std::string &filename, const std::string &dataPath);

Array4D<PRISMATIC_FLOAT_PRECISION> readDataSet4D_keepOrder(const std::string &filename, const std::string &dataPath);

// template <size_t N, class T>
// ArrayND<N, T> readDataSet(const std::string &filename, const std::string &dataPath, size_t );

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, PRISMATIC_FLOAT_PRECISION &val);

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, PRISMATIC_FLOAT_PRECISION *val);

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, int &val);

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, std::string &val);

void writeComplexDataSet(H5::Group group, const std::string &dsetname, const std::complex<PRISMATIC_FLOAT_PRECISION> *buffer, const hsize_t *mdims, const size_t &rank);

void writeRealDataSet(H5::Group group, const std::string &dsetname, const PRISMATIC_FLOAT_PRECISION *buffer, const hsize_t *mdims, const size_t &rank);

void writeScalarAttribute(H5::H5Object &object, const std::string &name, const int &data);

void writeScalarAttribute(H5::H5Object &object, const std::string &name, const PRISMATIC_FLOAT_PRECISION &data);

void writeScalarAttribute(H5::H5Object &object, const std::string &name, const std::string &data);

template <size_t N>
void readComplexDataSet(ArrayND<N, std::vector<std::complex<PRISMATIC_FLOAT_PRECISION>>> &output, const std::string &filename, const std::string &dataPath)
{
	H5::H5File input = H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	H5::DataSet dataset = input.openDataSet(dataPath.c_str());
	H5::DataSpace dataspace = dataset.getSpace();

	hsize_t dims_out[N];
	int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
	H5::DataSpace mspace(N,dims_out);

    std::array<size_t, N> data_dims;
    size_t totalSize = 1;
	for(auto i = 0; i < N; i++) totalSize *= dims_out[i];
    for(auto i = 0; i < N; i++) data_dims[N-1-i] = dims_out[i];

    std::vector<std::complex<PRISMATIC_FLOAT_PRECISION>> data_in(totalSize);
    dataset.read(&data_in[0], dataset.getDataType(), mspace, dataspace);

	output = ArrayND<N, std::vector<std::complex<PRISMATIC_FLOAT_PRECISION>>>(data_in, data_dims);

    mspace.close();
    dataspace.close();
    dataset.close();
    input.close();
};

int countDataGroups(H5::Group group, const std::string &basename);

int countDimensions(H5::Group group, const std::string &basename);

void configureSupergroup(H5::Group &new_sg,
						H5::Group &sourceExample,
						const std::vector<std::vector<PRISMATIC_FLOAT_PRECISION>> &sgdims,
						const std::vector<std::string> &sgdims_name,
						const std::vector<std::string> &sgdims_units);

void writeVirtualDataSet(H5::Group group,
						const std::string &dsetName,
						std::vector<H5::DataSet> &datasets,
						std::vector<std::vector<size_t>> indices);

void depthSeriesSG(H5::H5File &file);

std::string reducedDataSetName(std::string &fullPath);

void copyDataSet(H5::Group &targetGroup, H5::DataSet &source);
						
} //namespace Prismatic

#endif //PRISMATIC_FILEIO_H