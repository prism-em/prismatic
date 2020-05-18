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

#include "utility.h"
#include <complex>
#include "defines.h"
#include "configure.h"
#include "H5Cpp.h"
#include <string>
#include <stdio.h>
#ifdef _WIN32
#include <io.h>
#define access _access_s
#else
#include <unistd.h>
#endif
#include <mutex>
#include <thread>

namespace Prismatic
{

std::mutex write4D_lock;

std::pair<Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>>, Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>>>
upsamplePRISMProbe(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> probe,
				   const long dimj, const long dimi, long ys, long xs)
{
	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> realspace_probe;
	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> buffer_probe;
	Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> kspace_probe;

	buffer_probe = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>({{(size_t)dimj, (size_t)dimi}});
	//		std::cout << "dimj = " << dimj << std::endl;
	long ncy = probe.get_dimj() / 2;
	long ncx = probe.get_dimi() / 2;
	for (auto j = 0; j < probe.get_dimj(); ++j)
	{
		for (auto i = 0; i < probe.get_dimi(); ++i)
		{
			buffer_probe.at((dimj + ((j - ncy + ys) % dimj)) % dimj,
							(dimi + ((i - ncx + xs) % dimi)) % dimi) = probe.at(j, i);
			//				std::cout << "(dimj + ((j - ncy) % dimj)) % dimj= " << (dimj + ((j - ncy) % dimj)) % dimj<< std::endl;
			//				std::cout << "(j - ncy)= " << (j - ncy) << std::endl;
			//				std::cout << "(j - ncy) % dimj)= " << (j - ncy) % dimj<< std::endl;

			//				buffer_probe.at( (dimj + ((j - ncy) % dimj)) % dimj,
			//				                 (dimi + ((i - ncx) % dimi)) % dimi) = probe.at(j, i);
		}
	}
	std::unique_lock<std::mutex> gatekeeper(fftw_plan_lock);
	PRISMATIC_FFTW_PLAN plan = PRISMATIC_FFTW_PLAN_DFT_2D(buffer_probe.get_dimj(), buffer_probe.get_dimi(),
														  reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&buffer_probe[0]),
														  reinterpret_cast<PRISMATIC_FFTW_COMPLEX *>(&buffer_probe[0]),
														  FFTW_FORWARD, FFTW_ESTIMATE);
	gatekeeper.unlock();
	realspace_probe = buffer_probe;
	PRISMATIC_FFTW_EXECUTE(plan);
	kspace_probe = buffer_probe;
	gatekeeper.lock();
	PRISMATIC_FFTW_DESTROY_PLAN(plan);
	gatekeeper.unlock();
	return std::make_pair(realspace_probe, kspace_probe);
}

PRISMATIC_FLOAT_PRECISION computePearsonCorrelation(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> left,
													Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> right)
{
	PRISMATIC_FLOAT_PRECISION m1, m2, sigma1, sigma2, R;
	m1 = m2 = sigma1 = sigma2 = R = 0;

	for (auto &i : left)
		m1 += std::abs(i);
	for (auto &i : right)
		m2 += std::abs(i);

	m1 /= (left.size());
	m2 /= (right.size());

	for (auto &i : left)
		sigma1 += pow(std::abs(i) - m1, 2);
	for (auto &i : right)
		sigma2 += pow(std::abs(i) - m2, 2);

	sigma1 /= (left.size());
	sigma2 /= (right.size());

	sigma1 = sqrt(sigma1);
	sigma2 = sqrt(sigma2);
	for (auto i = 0; i < std::min(left.size(), right.size()); ++i)
	{
		R = R + (std::abs(left[i]) - m1) * (std::abs(right[i]) - m2);
	}
	R /= sqrt(left.size() * right.size());
	return R / (sigma1 * sigma2);
}
PRISMATIC_FLOAT_PRECISION computeRfactor(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> left,
										 Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> right)
{
	PRISMATIC_FLOAT_PRECISION accum, diffs;
	accum = diffs = 0;
	for (auto i = 0; i < std::min(left.size(), right.size()); ++i)
	{
		diffs += std::abs(left[i] - right[i]);
		accum += std::abs(left[i]);
	}
	return diffs / accum;
}

int nyquistProbes(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, size_t dim)
{
	int nProbes = ceil(4 * (pars.meta.probeSemiangle / pars.lambda) * pars.tiledCellDim[dim]);
	return nProbes;
}

std::string remove_extension(const std::string &filename)
{
	size_t lastdot = filename.find_last_of(".");
	if (lastdot == std::string::npos)
		return filename;
	return filename.substr(0, lastdot);
}

int testFilenameOutput(const std::string &filename)
{
	bool exists = !testExist(filename);
	bool write_ok = !testWrite(filename);
	//Check if file already exists and if we can write to it
	if (exists && write_ok)
	{
		std::cout << "Warning " << filename << " already exists and will be overwritten" << std::endl;
		return 2;
	}
	else if (exists && !write_ok)
	{
		std::cout << filename << " isn't an accessible write destination" << std::endl;
		return 0;
	}
	else
	{
		//If the file does not exist, check to see if we can open a file of that name
		std::ofstream f(filename, std::ios::binary | std::ios::out);
		if (f)
		{
			//If we can open such a file, close the file and delete it.
			f.close();
			std::remove(filename.c_str());
			return 1;
		}
		else
		{
			std::cout << filename << " isn't an accessible write destination" << std::endl;
			return 0;
		}
	}
}

int testWrite(const std::string &filename)
{
	int answer = access(filename.c_str(), 02); //W_OK = 02
	return answer;
}

int testExist(const std::string &filename)
{
	int answer = access(filename.c_str(), 00); //F_OK == 00
	return answer;
}

void setupOutputFile(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars)
{
	//create main groups
	H5::Group simulation(pars.outputFile.createGroup("/4DSTEM_simulation"));

	//set version attributes
	int maj_data = 0;
	int min_data = 5;
	int group_type = 2;

	H5::DataSpace attr_dataspace(H5S_SCALAR);

	H5::Attribute maj_attr = simulation.createAttribute("version_major", H5::PredType::NATIVE_INT, attr_dataspace);
	H5::Attribute min_attr = simulation.createAttribute("version_minor", H5::PredType::NATIVE_INT, attr_dataspace);
	H5::Attribute emd_group_type_attr = simulation.createAttribute("emd_group_type", H5::PredType::NATIVE_INT, attr_dataspace);

	maj_attr.write(H5::PredType::NATIVE_INT, &maj_data);
	min_attr.write(H5::PredType::NATIVE_INT, &min_data);
	emd_group_type_attr.write(H5::PredType::NATIVE_INT, &group_type);

	//data groups
	H5::Group data(simulation.createGroup("data"));
	H5::Group datacubes(data.createGroup("datacubes"));
	H5::Group dslices(data.createGroup("diffractionslices"));
	H5::Group rslices(data.createGroup("realslices"));
	H5::Group pointlists(data.createGroup("pointlists"));	//point lists and point list arrays are not used in prismatic
	H5::Group plarrays(data.createGroup("pointlistarrays")); //included here to maintain consistency with format

	//log group
	H5::Group log(simulation.createGroup("log"));

	//metadata groups
	H5::Group metadata(simulation.createGroup("metadata"));
	H5::Group metadata_0(metadata.createGroup("metadata_0")); //for consistency with py4DSTEM v0.4

	H5::Group original(metadata_0.createGroup("original"));
	H5::Group shortlist(original.createGroup("shortlist"));
	H5::Group all(original.createGroup("all"));
	H5::Group microscope(metadata_0.createGroup("microscope"));
	H5::Group sample(metadata_0.createGroup("sample"));
	H5::Group user(metadata_0.createGroup("user"));
	H5::Group calibration(metadata_0.createGroup("calibration"));
	H5::Group comments(metadata_0.createGroup("comments"));
}

//use dummy variable to overload float/double dependence
void setup4DOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const float dummy)
{
	H5::Group datacubes = pars.outputFile.openGroup("4DSTEM_simulation/data/datacubes");

	//shared properties
	std::string base_name = "CBED_array_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[4];
	data_dims[0] = {pars.xp.size()};
	data_dims[1] = {pars.yp.size()};
	hsize_t chunkDims[4];
	chunkDims[0] = chunkDims[1] = {1};
	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};
	hsize_t qx_dim[1];
	hsize_t qy_dim[1];

	Prismatic::Array1D<PRISMATIC_FLOAT_PRECISION> qx;
	Prismatic::Array1D<PRISMATIC_FLOAT_PRECISION> qy;
	long offset_qx;
	long offset_qy;

    size_t qxInd_max = 0;
    size_t qyInd_max = 0;

    if(pars.meta.crop4DOutput)
    {

        PRISMATIC_FLOAT_PRECISION qMax = pars.meta.crop4Damax / pars.lambda;

        for(auto i = 0; i < pars.qx.get_dimi(); i++)
        {
            if(pars.qx.at(i) < qMax)
            {
                qxInd_max++;
            }
            else
            {
                break;
            }
        }

        for(auto j = 0; j < pars.qy.get_dimi(); j++)
        {
            if(pars.qy.at(j) < qMax)
            {
                qyInd_max++;
            }
            else
            {
                break;
            }
        }

        qxInd_max *= 2;
        qyInd_max *= 2;
        //TODO: Figure out how to correctly sit the dimension offset for a cropped image
        offset_qx = 0;
        offset_qy = 0;
    }
    else
    {
        if (pars.meta.algorithm == Prismatic::Algorithm::Multislice)
        {
            qxInd_max = pars.psiProbeInit.get_dimi() / 2;
            qyInd_max = pars.psiProbeInit.get_dimj() / 2;
            offset_qx = pars.psiProbeInit.get_dimi() / 4;
            offset_qy = pars.psiProbeInit.get_dimj() / 4;
        }
        else
        {
            qxInd_max = pars.qx.get_dimi();
            qyInd_max = pars.qy.get_dimi();
            offset_qx = 0;
            offset_qy = 0;
        }
    }

	if (pars.meta.algorithm == Prismatic::Algorithm::Multislice)
	{
		data_dims[2] = {qxInd_max};
		data_dims[3] = {qyInd_max};
		qx_dim[0] = {qxInd_max};
		qy_dim[0] = {qyInd_max};
		qx = fftshift(pars.qx);
		qy = fftshift(pars.qy);
		chunkDims[2] = {qxInd_max};
		chunkDims[3] = {qyInd_max};
	}
	else
	{
		data_dims[2] = {qxInd_max};
		data_dims[3] = {qyInd_max};
		qx_dim[0] = {qxInd_max};
		qy_dim[0] = {qyInd_max};
		qx = pars.qx;
		qy = pars.qy;
		chunkDims[2] = {qxInd_max};
		chunkDims[3] = {qyInd_max};
	}

	for (auto n = 0; n < numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n);
		H5::Group CBED_slice_n(datacubes.createGroup(nth_name.c_str()));

		//write group type attribute
		H5::DataSpace attr1_dataspace(H5S_SCALAR);
		H5::Attribute emd_group_type = CBED_slice_n.createAttribute("emd_group_type", H5::PredType::NATIVE_INT, attr1_dataspace);
		int group_type = 1;
		emd_group_type.write(H5::PredType::NATIVE_INT, &group_type);

		//write metadata attribute
		H5::DataSpace attr2_dataspace(H5S_SCALAR);
		H5::Attribute metadata_group = CBED_slice_n.createAttribute("metadata", H5::PredType::NATIVE_INT, attr2_dataspace);
		int mgroup = 0;
		metadata_group.write(H5::PredType::NATIVE_INT, &mgroup);

		//setup data set chunking properties
		H5::DSetCreatPropList plist;
		plist.setChunk(4, chunkDims);

		//create dataset
		H5::DataSpace mspace(4, data_dims); //rank is 4
		H5::DataSet CBED_data = CBED_slice_n.createDataSet("datacube", H5::PredType::NATIVE_FLOAT, mspace, plist);
		mspace.close();

		//write dimensions
		H5::DataSpace str_name_ds(H5S_SCALAR);
		H5::StrType strdatatype(H5::PredType::C_S1, 256);

		H5::DataSpace dim1_mspace(1, rx_dim);
		H5::DataSpace dim2_mspace(1, ry_dim);
		H5::DataSpace dim3_mspace(1, qx_dim);
		H5::DataSpace dim4_mspace(1, qy_dim);

		H5::DataSet dim1 = CBED_slice_n.createDataSet("dim1", H5::PredType::NATIVE_FLOAT, dim1_mspace);
		H5::DataSet dim2 = CBED_slice_n.createDataSet("dim2", H5::PredType::NATIVE_FLOAT, dim2_mspace);
		H5::DataSet dim3 = CBED_slice_n.createDataSet("dim3", H5::PredType::NATIVE_FLOAT, dim3_mspace);
		H5::DataSet dim4 = CBED_slice_n.createDataSet("dim4", H5::PredType::NATIVE_FLOAT, dim4_mspace);

		H5::DataSpace dim1_fspace = dim1.getSpace();
		H5::DataSpace dim2_fspace = dim2.getSpace();
		H5::DataSpace dim3_fspace = dim3.getSpace();
		H5::DataSpace dim4_fspace = dim4.getSpace();

		dim1.write(&pars.xp[0], H5::PredType::NATIVE_FLOAT, dim1_mspace, dim1_fspace);
		dim2.write(&pars.yp[0], H5::PredType::NATIVE_FLOAT, dim2_mspace, dim2_fspace);
		dim3.write(&qx[offset_qx], H5::PredType::NATIVE_FLOAT, dim3_mspace, dim3_fspace);
		dim4.write(&qy[offset_qy], H5::PredType::NATIVE_FLOAT, dim4_mspace, dim4_fspace);

		//dimension attributes
		const H5std_string dim1_name_str("R_x");
		const H5std_string dim2_name_str("R_y");
		const H5std_string dim3_name_str("Q_x");
		const H5std_string dim4_name_str("Q_y");

		H5::Attribute dim1_name = dim1.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim2_name = dim2.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim3_name = dim3.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim4_name = dim4.createAttribute("name", strdatatype, str_name_ds);

		dim1_name.write(strdatatype, dim1_name_str);
		dim2_name.write(strdatatype, dim2_name_str);
		dim3_name.write(strdatatype, dim3_name_str);
		dim4_name.write(strdatatype, dim4_name_str);

		const H5std_string dim1_unit_str("[Å]");
		const H5std_string dim2_unit_str("[Å]");
		const H5std_string dim3_unit_str("[Å^-1]");
		const H5std_string dim4_unit_str("[Å^-1]");

		H5::Attribute dim1_unit = dim1.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim2_unit = dim2.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim3_unit = dim3.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim4_unit = dim4.createAttribute("units", strdatatype, str_name_ds);

		dim1_unit.write(strdatatype, dim1_unit_str);
		dim2_unit.write(strdatatype, dim2_unit_str);
		dim3_unit.write(strdatatype, dim3_unit_str);
		dim4_unit.write(strdatatype, dim4_unit_str);

		CBED_slice_n.close();
	}

	datacubes.close();
};

//use dummy variable to overload float/double dependence
void setup4DOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const double dummy)
{
	H5::Group datacubes = pars.outputFile.openGroup("4DSTEM_simulation/data/datacubes");

	//shared properties
	std::string base_name = "CBED_array_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[4];
	data_dims[0] = {pars.xp.size()};
	data_dims[1] = {pars.yp.size()};
	hsize_t chunkDims[4];
	chunkDims[0] = chunkDims[1] = {1};
	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};
	hsize_t qx_dim[1];
	hsize_t qy_dim[1];
	Prismatic::Array1D<PRISMATIC_FLOAT_PRECISION> qx;
	Prismatic::Array1D<PRISMATIC_FLOAT_PRECISION> qy;
	long offset_qx;
	long offset_qy;

	if (pars.meta.algorithm == Prismatic::Algorithm::Multislice)
	{
		data_dims[2] = {pars.psiProbeInit.get_dimi() / 2};
		data_dims[3] = {pars.psiProbeInit.get_dimj() / 2};
		qx_dim[0] = {pars.psiProbeInit.get_dimi() / 2};
		qy_dim[0] = {pars.psiProbeInit.get_dimj() / 2};
		qx = fftshift(pars.qx);
		qy = fftshift(pars.qy);
		offset_qx = pars.psiProbeInit.get_dimi() / 4;
		offset_qy = pars.psiProbeInit.get_dimj() / 4;
		chunkDims[2] = {pars.psiProbeInit.get_dimi() / 2};
		chunkDims[3] = {pars.psiProbeInit.get_dimj() / 2};
	}
	else
	{
		data_dims[2] = {pars.qx.get_dimi()};
		data_dims[3] = {pars.qy.get_dimi()};
		qx_dim[0] = {pars.qx.get_dimi()};
		qy_dim[0] = {pars.qy.get_dimi()};
		qx = pars.qx;
		qy = pars.qy;
		offset_qx = 0;
		offset_qy = 0;
		chunkDims[2] = {pars.qx.get_dimi()};
		chunkDims[3] = {pars.qy.get_dimi()};
		//std::cout << "Probe size: " << pars.psiProbeInit.get_dimi() << std::endl;
	}

	for (auto n = 0; n < numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n);
		H5::Group CBED_slice_n(datacubes.createGroup(nth_name.c_str()));

		//write group type attribute
		H5::DataSpace attr1_dataspace(H5S_SCALAR);
		H5::Attribute emd_group_type = CBED_slice_n.createAttribute("emd_group_type", H5::PredType::NATIVE_INT, attr1_dataspace);
		int group_type = 1;
		emd_group_type.write(H5::PredType::NATIVE_INT, &group_type);

		//write metadata attribute
		H5::DataSpace attr2_dataspace(H5S_SCALAR);
		H5::Attribute metadata_group = CBED_slice_n.createAttribute("metadata", H5::PredType::NATIVE_INT, attr2_dataspace);
		int mgroup = 0;
		metadata_group.write(H5::PredType::NATIVE_INT, &mgroup);

		//set chunk properties
		H5::DSetCreatPropList plist;
		plist.setChunk(4, chunkDims);

		//create dataset
		H5::DataSpace mspace(4, data_dims); //rank is 4
		H5::DataSet CBED_data = CBED_slice_n.createDataSet("datacube", H5::PredType::NATIVE_DOUBLE, mspace, plist);
		mspace.close();

		//write dimensions
		H5::DataSpace str_name_ds(H5S_SCALAR);
		H5::StrType strdatatype(H5::PredType::C_S1, 256);

		H5::DataSpace dim1_mspace(1, rx_dim);
		H5::DataSpace dim2_mspace(1, ry_dim);
		H5::DataSpace dim3_mspace(1, qx_dim);
		H5::DataSpace dim4_mspace(1, qy_dim);

		H5::DataSet dim1 = CBED_slice_n.createDataSet("dim1", H5::PredType::NATIVE_DOUBLE, dim1_mspace);
		H5::DataSet dim2 = CBED_slice_n.createDataSet("dim2", H5::PredType::NATIVE_DOUBLE, dim2_mspace);
		H5::DataSet dim3 = CBED_slice_n.createDataSet("dim3", H5::PredType::NATIVE_DOUBLE, dim3_mspace);
		H5::DataSet dim4 = CBED_slice_n.createDataSet("dim4", H5::PredType::NATIVE_DOUBLE, dim4_mspace);

		H5::DataSpace dim1_fspace = dim1.getSpace();
		H5::DataSpace dim2_fspace = dim2.getSpace();
		H5::DataSpace dim3_fspace = dim3.getSpace();
		H5::DataSpace dim4_fspace = dim4.getSpace();

		dim1.write(&pars.xp[0], H5::PredType::NATIVE_DOUBLE, dim1_mspace, dim1_fspace);
		dim2.write(&pars.yp[0], H5::PredType::NATIVE_DOUBLE, dim2_mspace, dim2_fspace);
		dim3.write(&qx[offset_qx], H5::PredType::NATIVE_DOUBLE, dim3_mspace, dim3_fspace);
		dim4.write(&qy[offset_qy], H5::PredType::NATIVE_DOUBLE, dim4_mspace, dim4_fspace);

		//dimension attributes
		const H5std_string dim1_name_str("R_x");
		const H5std_string dim2_name_str("R_y");
		const H5std_string dim3_name_str("Q_x");
		const H5std_string dim4_name_str("Q_y");

		H5::Attribute dim1_name = dim1.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim2_name = dim2.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim3_name = dim3.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim4_name = dim4.createAttribute("name", strdatatype, str_name_ds);

		dim1_name.write(strdatatype, dim1_name_str);
		dim2_name.write(strdatatype, dim2_name_str);
		dim3_name.write(strdatatype, dim3_name_str);
		dim4_name.write(strdatatype, dim4_name_str);

		const H5std_string dim1_unit_str("[Å]");
		const H5std_string dim2_unit_str("[Å]");
		const H5std_string dim3_unit_str("[Å^-1]");
		const H5std_string dim4_unit_str("[Å^-1]");

		H5::Attribute dim1_unit = dim1.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim2_unit = dim2.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim3_unit = dim3.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim4_unit = dim4.createAttribute("units", strdatatype, str_name_ds);

		dim1_unit.write(strdatatype, dim1_unit_str);
		dim2_unit.write(strdatatype, dim2_unit_str);
		dim3_unit.write(strdatatype, dim3_unit_str);
		dim4_unit.write(strdatatype, dim4_unit_str);

		CBED_slice_n.close();
	}

	datacubes.close();
};

void setupVDOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const float dummy)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	//shared properties
	std::string base_name = "virtual_detector_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[3];
	data_dims[0] = {pars.xp.size()};
	data_dims[1] = {pars.yp.size()};
	data_dims[2] = {pars.Ndet};

	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};
	hsize_t bin_dim[1] = {pars.Ndet};

	for (auto n = 0; n < numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n);
		H5::Group VD_slice_n(realslices.createGroup(nth_name.c_str()));

		//write group type attribute
		H5::DataSpace attr1_dataspace(H5S_SCALAR);
		H5::Attribute emd_group_type = VD_slice_n.createAttribute("emd_group_type", H5::PredType::NATIVE_INT, attr1_dataspace);
		int group_type = 1;
		emd_group_type.write(H5::PredType::NATIVE_INT, &group_type);

		//write metadata attribute
		H5::DataSpace attr2_dataspace(H5S_SCALAR);
		H5::Attribute metadata_group = VD_slice_n.createAttribute("metadata", H5::PredType::NATIVE_INT, attr2_dataspace);
		int mgroup = 0;
		metadata_group.write(H5::PredType::NATIVE_INT, &mgroup);

		//create datasets
		H5::DataSpace mspace(3, data_dims); //rank is 2 for each realslice
		H5::DataSet VD_data = VD_slice_n.createDataSet("realslice", H5::PredType::NATIVE_FLOAT, mspace);
		VD_data.close();
		mspace.close();

		//write dimensions
		H5::DataSpace str_name_ds(H5S_SCALAR);
		H5::StrType strdatatype(H5::PredType::C_S1, 256);

		H5::DataSpace dim1_mspace(1, rx_dim);
		H5::DataSpace dim2_mspace(1, ry_dim);
		H5::DataSpace dim3_mspace(1, bin_dim);

		H5::DataSet dim1 = VD_slice_n.createDataSet("dim1", H5::PredType::NATIVE_FLOAT, dim1_mspace);
		H5::DataSet dim2 = VD_slice_n.createDataSet("dim2", H5::PredType::NATIVE_FLOAT, dim2_mspace);
		H5::DataSet dim3 = VD_slice_n.createDataSet("dim3", H5::PredType::NATIVE_FLOAT, dim3_mspace);

		H5::DataSpace dim1_fspace = dim1.getSpace();
		H5::DataSpace dim2_fspace = dim2.getSpace();
		H5::DataSpace dim3_fspace = dim3.getSpace();

		dim1.write(&pars.xp[0], H5::PredType::NATIVE_FLOAT, dim1_mspace, dim1_fspace);
		dim2.write(&pars.yp[0], H5::PredType::NATIVE_FLOAT, dim2_mspace, dim2_fspace);
		dim3.write(&pars.detectorAngles[0], H5::PredType::NATIVE_FLOAT, dim3_mspace, dim3_fspace);
		//dimension attributes
		const H5std_string dim1_name_str("R_x");
		const H5std_string dim2_name_str("R_y");
		const H5std_string dim3_name_str("bin_outer_angle");

		H5::Attribute dim1_name = dim1.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim2_name = dim2.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim3_name = dim3.createAttribute("name", strdatatype, str_name_ds);

		dim1_name.write(strdatatype, dim1_name_str);
		dim2_name.write(strdatatype, dim2_name_str);
		dim3_name.write(strdatatype, dim3_name_str);

		const H5std_string dim1_unit_str("[Å]");
		const H5std_string dim2_unit_str("[Å]");
		const H5std_string dim3_unit_str("[rad]");

		H5::Attribute dim1_unit = dim1.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim2_unit = dim2.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim3_unit = dim3.createAttribute("units", strdatatype, str_name_ds);

		dim1_unit.write(strdatatype, dim1_unit_str);
		dim2_unit.write(strdatatype, dim2_unit_str);
		dim3_unit.write(strdatatype, dim3_unit_str);

		VD_slice_n.close();
	}

	realslices.close();
};

void setupVDOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const double dummy)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	//shared properties
	std::string base_name = "virtual_detector_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[3];
	data_dims[0] = {pars.xp.size()};
	data_dims[1] = {pars.yp.size()};
	data_dims[2] = {pars.Ndet};

	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};
	hsize_t bin_dim[1] = {pars.Ndet};

	for (auto n = 0; n < numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n);
		H5::Group VD_slice_n(realslices.createGroup(nth_name.c_str()));

		//write group type attribute
		H5::DataSpace attr1_dataspace(H5S_SCALAR);
		H5::Attribute emd_group_type = VD_slice_n.createAttribute("emd_group_type", H5::PredType::NATIVE_INT, attr1_dataspace);
		int group_type = 1;
		emd_group_type.write(H5::PredType::NATIVE_INT, &group_type);

		//write metadata attribute
		H5::DataSpace attr2_dataspace(H5S_SCALAR);
		H5::Attribute metadata_group = VD_slice_n.createAttribute("metadata", H5::PredType::NATIVE_INT, attr2_dataspace);
		int mgroup = 0;
		metadata_group.write(H5::PredType::NATIVE_INT, &mgroup);

		//create datasets
		H5::DataSpace mspace(3, data_dims); //rank is 2 for each realslice
		H5::DataSet VD_data = VD_slice_n.createDataSet("realslice", H5::PredType::NATIVE_DOUBLE, mspace);
		VD_data.close();
		mspace.close();

		//write dimensions
		H5::DataSpace str_name_ds(H5S_SCALAR);
		H5::StrType strdatatype(H5::PredType::C_S1, 256);

		H5::DataSpace dim1_mspace(1, rx_dim);
		H5::DataSpace dim2_mspace(1, ry_dim);
		H5::DataSpace dim3_mspace(1, bin_dim);

		H5::DataSet dim1 = VD_slice_n.createDataSet("dim1", H5::PredType::NATIVE_DOUBLE, dim1_mspace);
		H5::DataSet dim2 = VD_slice_n.createDataSet("dim2", H5::PredType::NATIVE_DOUBLE, dim2_mspace);
		H5::DataSet dim3 = VD_slice_n.createDataSet("dim3", H5::PredType::NATIVE_DOUBLE, dim3_mspace);

		H5::DataSpace dim1_fspace = dim1.getSpace();
		H5::DataSpace dim2_fspace = dim2.getSpace();
		H5::DataSpace dim3_fspace = dim3.getSpace();

		dim1.write(&pars.xp[0], H5::PredType::NATIVE_DOUBLE, dim1_mspace, dim1_fspace);
		dim2.write(&pars.yp[0], H5::PredType::NATIVE_DOUBLE, dim2_mspace, dim2_fspace);
		dim3.write(&pars.detectorAngles[0], H5::PredType::NATIVE_DOUBLE, dim3_mspace, dim3_fspace);

		//dimension attributes
		const H5std_string dim1_name_str("R_x");
		const H5std_string dim2_name_str("R_y");
		const H5std_string dim3_name_str("bin_outer_angle");

		H5::Attribute dim1_name = dim1.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim2_name = dim2.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim3_name = dim3.createAttribute("name", strdatatype, str_name_ds);

		dim1_name.write(strdatatype, dim1_name_str);
		dim2_name.write(strdatatype, dim2_name_str);
		dim3_name.write(strdatatype, dim3_name_str);

		const H5std_string dim1_unit_str("[Å]");
		const H5std_string dim2_unit_str("[Å]");
		const H5std_string dim3_unit_str("[rad]");

		H5::Attribute dim1_unit = dim1.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim2_unit = dim2.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim3_unit = dim3.createAttribute("units", strdatatype, str_name_ds);

		dim1_unit.write(strdatatype, dim1_unit_str);
		dim2_unit.write(strdatatype, dim2_unit_str);
		dim3_unit.write(strdatatype, dim3_unit_str);

		VD_slice_n.close();
	}

	realslices.close();
};

void setup2DOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const float dummy)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	//shared properties
	std::string base_name = "annular_detector_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[2];
	data_dims[0] = {pars.xp.size()};
	data_dims[1] = {pars.yp.size()};

	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};

	for (auto n = 0; n < numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n);
		H5::Group annular_slice_n(realslices.createGroup(nth_name.c_str()));

		//write group type attribute
		H5::DataSpace attr1_dataspace(H5S_SCALAR);
		H5::Attribute emd_group_type = annular_slice_n.createAttribute("emd_group_type", H5::PredType::NATIVE_INT, attr1_dataspace);
		int group_type = 1;
		emd_group_type.write(H5::PredType::NATIVE_INT, &group_type);

		//write metadata attribute
		H5::DataSpace attr2_dataspace(H5S_SCALAR);
		H5::Attribute metadata_group = annular_slice_n.createAttribute("metadata", H5::PredType::NATIVE_INT, attr2_dataspace);
		int mgroup = 0;
		metadata_group.write(H5::PredType::NATIVE_INT, &mgroup);

		//write depth attribute
		int depth = 1;
		H5::DataSpace attr3_dataspace(H5S_SCALAR);
		H5::Attribute depth_attr = annular_slice_n.createAttribute("depth", H5::PredType::NATIVE_INT, attr3_dataspace);
		depth_attr.write(H5::PredType::NATIVE_INT, &depth);

		//create dataset
		H5::DataSpace mspace(2, data_dims); //rank is 3
		H5::DataSet CBED_data = annular_slice_n.createDataSet("realslice", H5::PredType::NATIVE_FLOAT, mspace);
		mspace.close();

		//write dimensions
		H5::DataSpace str_name_ds(H5S_SCALAR);
		H5::StrType strdatatype(H5::PredType::C_S1, 256);

		H5::DataSpace dim1_mspace(1, rx_dim);
		H5::DataSpace dim2_mspace(1, ry_dim);

		H5::DataSet dim1 = annular_slice_n.createDataSet("dim1", H5::PredType::NATIVE_FLOAT, dim1_mspace);
		H5::DataSet dim2 = annular_slice_n.createDataSet("dim2", H5::PredType::NATIVE_FLOAT, dim2_mspace);

		H5::DataSpace dim1_fspace = dim1.getSpace();
		H5::DataSpace dim2_fspace = dim2.getSpace();

		dim1.write(&pars.xp[0], H5::PredType::NATIVE_FLOAT, dim1_mspace, dim1_fspace);
		dim2.write(&pars.yp[0], H5::PredType::NATIVE_FLOAT, dim2_mspace, dim2_fspace);

		//dimension attributes
		const H5std_string dim1_name_str("R_x");
		const H5std_string dim2_name_str("R_y");

		H5::Attribute dim1_name = dim1.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim2_name = dim2.createAttribute("name", strdatatype, str_name_ds);

		dim1_name.write(strdatatype, dim1_name_str);
		dim2_name.write(strdatatype, dim2_name_str);

		const H5std_string dim1_unit_str("[Å]");
		const H5std_string dim2_unit_str("[Å]");

		H5::Attribute dim1_unit = dim1.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim2_unit = dim2.createAttribute("units", strdatatype, str_name_ds);

		dim1_unit.write(strdatatype, dim1_unit_str);
		dim2_unit.write(strdatatype, dim2_unit_str);

		annular_slice_n.close();
	}

	realslices.close();
};

void setup2DOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const double dummy)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	//shared properties
	std::string base_name = "annular_detector_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[2];
	data_dims[0] = {pars.xp.size()};
	data_dims[1] = {pars.yp.size()};

	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};

	for (auto n = 0; n < numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n);
		H5::Group annular_slice_n(realslices.createGroup(nth_name.c_str()));

		//write group type attribute
		H5::DataSpace attr1_dataspace(H5S_SCALAR);
		H5::Attribute emd_group_type = annular_slice_n.createAttribute("emd_group_type", H5::PredType::NATIVE_INT, attr1_dataspace);
		int group_type = 1;
		emd_group_type.write(H5::PredType::NATIVE_INT, &group_type);

		//write metadata attribute
		H5::DataSpace attr2_dataspace(H5S_SCALAR);
		H5::Attribute metadata_group = annular_slice_n.createAttribute("metadata", H5::PredType::NATIVE_INT, attr2_dataspace);
		int mgroup = 0;
		metadata_group.write(H5::PredType::NATIVE_INT, &mgroup);

		//write depth attribute
		int depth = 1;
		H5::DataSpace attr3_dataspace(H5S_SCALAR);
		H5::Attribute depth_attr = annular_slice_n.createAttribute("depth", H5::PredType::NATIVE_INT, attr3_dataspace);
		depth_attr.write(H5::PredType::NATIVE_INT, &depth);

		//create dataset
		H5::DataSpace mspace(2, data_dims); //rank is 2
		H5::DataSet CBED_data = annular_slice_n.createDataSet("realslice", H5::PredType::NATIVE_DOUBLE, mspace);
		mspace.close();

		//write dimensions
		H5::DataSpace str_name_ds(H5S_SCALAR);
		H5::StrType strdatatype(H5::PredType::C_S1, 256);

		H5::DataSpace dim1_mspace(1, rx_dim);
		H5::DataSpace dim2_mspace(1, ry_dim);

		H5::DataSet dim1 = annular_slice_n.createDataSet("dim1", H5::PredType::NATIVE_DOUBLE, dim1_mspace);
		H5::DataSet dim2 = annular_slice_n.createDataSet("dim2", H5::PredType::NATIVE_DOUBLE, dim2_mspace);

		H5::DataSpace dim1_fspace = dim1.getSpace();
		H5::DataSpace dim2_fspace = dim2.getSpace();

		dim1.write(&pars.xp[0], H5::PredType::NATIVE_DOUBLE, dim1_mspace, dim1_fspace);
		dim2.write(&pars.yp[0], H5::PredType::NATIVE_DOUBLE, dim2_mspace, dim2_fspace);

		//dimension attributes
		const H5std_string dim1_name_str("R_x");
		const H5std_string dim2_name_str("R_y");

		H5::Attribute dim1_name = dim1.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim2_name = dim2.createAttribute("name", strdatatype, str_name_ds);

		dim1_name.write(strdatatype, dim1_name_str);
		dim2_name.write(strdatatype, dim2_name_str);

		const H5std_string dim1_unit_str("[Å]");
		const H5std_string dim2_unit_str("[Å]");

		H5::Attribute dim1_unit = dim1.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim2_unit = dim2.createAttribute("units", strdatatype, str_name_ds);

		dim1_unit.write(strdatatype, dim1_unit_str);
		dim2_unit.write(strdatatype, dim2_unit_str);

		annular_slice_n.close();
	}

	realslices.close();
};

void setupDPCOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const float dummy)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	//shared properties
	std::string base_name = "DPC_CoM_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[3];
	data_dims[0] = {pars.xp.size()};
	data_dims[1] = {pars.yp.size()};
	data_dims[2] = {2};

	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};
	hsize_t str_dim[1] = {2};

	for (auto n = 0; n < numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n);
		H5::Group DPC_CoM_slice_n(realslices.createGroup(nth_name.c_str()));

		//write group type attribute
		H5::DataSpace attr1_dataspace(H5S_SCALAR);
		H5::Attribute emd_group_type = DPC_CoM_slice_n.createAttribute("emd_group_type", H5::PredType::NATIVE_INT, attr1_dataspace);
		int group_type = 1;
		emd_group_type.write(H5::PredType::NATIVE_INT, &group_type);

		//write metadata attribute
		H5::DataSpace attr2_dataspace(H5S_SCALAR);
		H5::Attribute metadata_group = DPC_CoM_slice_n.createAttribute("metadata", H5::PredType::NATIVE_INT, attr2_dataspace);
		int mgroup = 0;
		metadata_group.write(H5::PredType::NATIVE_INT, &mgroup);

		//create dataset
		//create dataset
		H5::DataSpace mspace(3, data_dims); //rank is 3
		H5::DataSet DPC_data = DPC_CoM_slice_n.createDataSet("realslice", H5::PredType::NATIVE_FLOAT, mspace);
		mspace.close();

		//write dimensions
		H5::DataSpace str_name_ds(H5S_SCALAR);
		H5::StrType strdatatype(H5::PredType::C_S1, 256);

		H5::DataSpace dim1_mspace(1, rx_dim);
		H5::DataSpace dim2_mspace(1, ry_dim);
		H5::DataSpace dim3_mspace(1, str_dim);

		H5::DataSet dim1 = DPC_CoM_slice_n.createDataSet("dim1", H5::PredType::NATIVE_FLOAT, dim1_mspace);
		H5::DataSet dim2 = DPC_CoM_slice_n.createDataSet("dim2", H5::PredType::NATIVE_FLOAT, dim2_mspace);

		H5::DataSpace dim1_fspace = dim1.getSpace();
		H5::DataSpace dim2_fspace = dim2.getSpace();

		dim1.write(&pars.xp[0], H5::PredType::NATIVE_FLOAT, dim1_mspace, dim1_fspace);
		dim2.write(&pars.yp[0], H5::PredType::NATIVE_FLOAT, dim2_mspace, dim2_fspace);

		H5::DataSet dim3 = DPC_CoM_slice_n.createDataSet("dim3", strdatatype, dim3_mspace);
		H5std_string dpc_x("DPC_CoM_x");
		H5std_string dpc_y("DPC_CoM_y");
		H5std_string str_buffer_array[2];
		str_buffer_array[0] = dpc_x;
		str_buffer_array[1] = dpc_y;

		writeStringArray(dim3, str_buffer_array, 2);

		//dimension attributes
		const H5std_string dim1_name_str("R_x");
		const H5std_string dim2_name_str("R_y");

		H5::Attribute dim1_name = dim1.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim2_name = dim2.createAttribute("name", strdatatype, str_name_ds);

		dim1_name.write(strdatatype, dim1_name_str);
		dim2_name.write(strdatatype, dim2_name_str);

		const H5std_string dim1_unit_str("[Å]");
		const H5std_string dim2_unit_str("[Å]");

		H5::Attribute dim1_unit = dim1.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim2_unit = dim2.createAttribute("units", strdatatype, str_name_ds);

		dim1_unit.write(strdatatype, dim1_unit_str);
		dim2_unit.write(strdatatype, dim2_unit_str);

		DPC_CoM_slice_n.close();
	}

	realslices.close();
};

void setupDPCOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers, const double dummy)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	//shared properties
	std::string base_name = "DPC_CoM_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[3];
	data_dims[0] = {pars.xp.size()};
	data_dims[1] = {pars.yp.size()};
	data_dims[2] = {2};

	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};
	hsize_t str_dim[1] = {2};

	for (auto n = 0; n < numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n);
		H5::Group DPC_CoM_slice_n(realslices.createGroup(nth_name.c_str()));

		//write group type attribute
		H5::DataSpace attr1_dataspace(H5S_SCALAR);
		H5::Attribute emd_group_type = DPC_CoM_slice_n.createAttribute("emd_group_type", H5::PredType::NATIVE_INT, attr1_dataspace);
		int group_type = 1;
		emd_group_type.write(H5::PredType::NATIVE_INT, &group_type);

		//write metadata attribute
		H5::DataSpace attr2_dataspace(H5S_SCALAR);
		H5::Attribute metadata_group = DPC_CoM_slice_n.createAttribute("metadata", H5::PredType::NATIVE_INT, attr2_dataspace);
		int mgroup = 0;
		metadata_group.write(H5::PredType::NATIVE_INT, &mgroup);

		//create dataset

		H5::DataSpace mspace(3, data_dims); //rank is 3
		H5::DataSet DPC_data = DPC_CoM_slice_n.createDataSet("realslice", H5::PredType::NATIVE_FLOAT, mspace);
		mspace.close();

		//write dimensions
		H5::DataSpace str_name_ds(H5S_SCALAR);
		H5::StrType strdatatype(H5::PredType::C_S1, 256);

		H5::DataSpace dim1_mspace(1, rx_dim);
		H5::DataSpace dim2_mspace(1, ry_dim);
		H5::DataSpace dim3_mspace(1, str_dim);

		H5::DataSet dim1 = DPC_CoM_slice_n.createDataSet("dim1", H5::PredType::NATIVE_DOUBLE, dim1_mspace);
		H5::DataSet dim2 = DPC_CoM_slice_n.createDataSet("dim2", H5::PredType::NATIVE_DOUBLE, dim2_mspace);

		H5::DataSpace dim1_fspace = dim1.getSpace();
		H5::DataSpace dim2_fspace = dim2.getSpace();

		dim1.write(&pars.xp[0], H5::PredType::NATIVE_DOUBLE, dim1_mspace, dim1_fspace);
		dim2.write(&pars.yp[0], H5::PredType::NATIVE_DOUBLE, dim2_mspace, dim2_fspace);

		H5::DataSet dim3 = DPC_CoM_slice_n.createDataSet("dim3", strdatatype, dim3_mspace);
		H5std_string dpc_x("DPC_CoM_x");
		H5std_string dpc_y("DPC_CoM_y");
		H5std_string str_buffer_array[2];
		str_buffer_array[0] = dpc_x;
		str_buffer_array[1] = dpc_y;

		writeStringArray(dim3, str_buffer_array, 2);

		//dimension attributes
		const H5std_string dim1_name_str("R_x");
		const H5std_string dim2_name_str("R_y");

		H5::Attribute dim1_name = dim1.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim2_name = dim2.createAttribute("name", strdatatype, str_name_ds);

		dim1_name.write(strdatatype, dim1_name_str);
		dim2_name.write(strdatatype, dim2_name_str);

		const H5std_string dim1_unit_str("[Å]");
		const H5std_string dim2_unit_str("[Å]");

		H5::Attribute dim1_unit = dim1.createAttribute("units", strdatatype, str_name_ds);
		H5::Attribute dim2_unit = dim2.createAttribute("units", strdatatype, str_name_ds);

		dim1_unit.write(strdatatype, dim1_unit_str);
		dim2_unit.write(strdatatype, dim2_unit_str);

		DPC_CoM_slice_n.close();
	}

	realslices.close();
};

void writeRealSlice(H5::DataSet dataset, const float *buffer, const hsize_t *mdims)
{
	H5::DataSpace fspace = dataset.getSpace(); //all realslices have data written all at once
	H5::DataSpace mspace(2, mdims);			   //rank = 2

	dataset.write(buffer, H5::PredType::NATIVE_FLOAT, mspace, fspace);

	fspace.close();
	mspace.close();
}

void writeRealSlice(H5::DataSet dataset, const double *buffer, const hsize_t *mdims)
{
	H5::DataSpace fspace = dataset.getSpace(); //all realslices have data written all at once
	H5::DataSpace mspace(2, mdims);			   //rank = 2

	dataset.write(buffer, H5::PredType::NATIVE_DOUBLE, mspace, fspace);

	fspace.close();
	mspace.close();
}

void writeDatacube3D(H5::DataSet dataset, const float *buffer, const hsize_t *mdims)
{
	//set up file and memory spaces
	H5::DataSpace fspace = dataset.getSpace(); //all 3D cubes will write full buffer at once
	H5::DataSpace mspace(3, mdims);			   //rank = 3

	dataset.write(buffer, H5::PredType::NATIVE_FLOAT, mspace, fspace);

	fspace.close();
	mspace.close();
};

void writeDatacube3D(H5::DataSet dataset, const double *buffer, const hsize_t *mdims)
{
	//set up file and memory spaces
	H5::DataSpace fspace = dataset.getSpace(); //all 3D cubes will write full buffer at once
	H5::DataSpace mspace(3, mdims);			   //rank = 3

	dataset.write(buffer, H5::PredType::NATIVE_DOUBLE, mspace, fspace);

	fspace.close();
	mspace.close();
};

//for 4D writes, need to first read the data set and then add; this way, FP are accounted for
void writeDatacube4D(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, float *buffer, const hsize_t *mdims, const hsize_t *offset, const float numFP, const std::string nameString)
{
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
    float *finalBuffer = (float *)malloc(mdims[0] * mdims[1] * mdims[2] * mdims[3] * sizeof(float));
    for (auto i = 0; i < mdims[2]; i++)
    {
        for (auto j = 0; j < mdims[3]; j++)
        {
            finalBuffer[i * mdims[3] + j] = buffer[j * mdims[2] + i];
        }
    }

    //add frozen phonon set
    float *readBuffer = (float *)malloc(mdims[0] * mdims[1] * mdims[2] * mdims[3] * sizeof(float));
    dataset.read(&readBuffer[0], H5::PredType::NATIVE_FLOAT, mspace, fspace);
    for (auto i = 0; i < mdims[0] * mdims[1] * mdims[2] * mdims[3]; i++)
        finalBuffer[i] += readBuffer[i];
    free(readBuffer);

    dataset.write(finalBuffer, H5::PredType::NATIVE_FLOAT, mspace, fspace);
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

void writeDatacube4D(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, double *buffer, const hsize_t *mdims, const hsize_t *offset, const double numFP, const std::string nameString)
{
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
	double *finalBuffer = (double *)malloc(mdims[0] * mdims[1] * mdims[2] * mdims[3] * sizeof(double));
	for (auto i = 0; i < mdims[2]; i++)
	{
		for (auto j = 0; j < mdims[3]; j++)
		{
			finalBuffer[i * mdims[3] + j] = buffer[j * mdims[2] + i];
		}
	}

	//add frozen phonon set
	double *readBuffer = (double *)malloc(mdims[0] * mdims[1] * mdims[2] * mdims[3] * sizeof(double));
	dataset.read(&readBuffer[0], H5::PredType::NATIVE_DOUBLE, mspace, fspace);
	for (auto i = 0; i < mdims[0] * mdims[1] * mdims[2] * mdims[3]; i++)
		finalBuffer[i] += readBuffer[i];
	free(readBuffer);

	dataset.write(finalBuffer, H5::PredType::NATIVE_DOUBLE, mspace, fspace);
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

void writeStringArray(H5::DataSet dataset, H5std_string *string_array, const hsize_t elements)
{
	//assumes that we are writing a 1 dimensional array of strings- used pretty much only for DPC
	H5::StrType strdatatype(H5::PredType::C_S1, 256);
	hsize_t offset[1] = {0};
	H5::DataSpace fspace = dataset.getSpace();
	hsize_t str_write_dim[1] = {1};
	H5::DataSpace mspace(1, str_write_dim);

	for (hsize_t i = 0; i < elements; i++)
	{
		offset[0] = {i};
		fspace.selectHyperslab(H5S_SELECT_SET, str_write_dim, offset);
		dataset.write(string_array[i], strdatatype, mspace, fspace);
	}
	fspace.close();
	mspace.close();
}

std::string getDigitString(int digit)
{
	char buffer[20];
	sprintf(buffer, "%04d", digit);
	std::string output = buffer;
	return output;
};

void writeMetadata(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, float dummy)
{
	//set up group
	H5::Group metadata = pars.outputFile.openGroup("4DSTEM_simulation/metadata/metadata_0/original");
	H5::Group sim_params = metadata.createGroup("simulation_parameters");

	//write all parameters as attributes

	//create common dataspaces
	H5::DataSpace str_name_ds(H5S_SCALAR); //string dataspaces and types
	H5::StrType strdatatype(H5::PredType::C_S1, 256);
	H5::DataSpace scalar_attr(H5S_SCALAR);

	//initialize string parameter data
	H5std_string algorithm;
	if (pars.meta.algorithm == Prismatic::Algorithm::Multislice)
	{
		algorithm = "m";
	}
	else
	{
		algorithm = "p";
	}

	const H5std_string filenameAtoms(pars.meta.filenameAtoms);

	//create string attributes
	H5::Attribute atoms_attr = sim_params.createAttribute("i", strdatatype, str_name_ds);
	H5::Attribute alg_attr = sim_params.createAttribute("a", strdatatype, str_name_ds);

	//create scalar logical/integer attributes
	H5::Attribute fx_attr = sim_params.createAttribute("fx", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute fy_attr = sim_params.createAttribute("fy", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute numFP_attr = sim_params.createAttribute("F", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute numSlices_attr = sim_params.createAttribute("ns", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute te_attr = sim_params.createAttribute("te", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute oc_attr = sim_params.createAttribute("oc", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute save3D_attr = sim_params.createAttribute("3D", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute save4D_attr = sim_params.createAttribute("4D", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute saveDPC_attr = sim_params.createAttribute("DPC", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute savePS_attr = sim_params.createAttribute("ps", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute nyquist_attr = sim_params.createAttribute("nqs", H5::PredType::NATIVE_INT, scalar_attr);

	//create scalar float/double attributes (changes based on prismatic float precision)
	H5::Attribute px_attr = sim_params.createAttribute("px", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute py_attr = sim_params.createAttribute("py", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute potBound_attr = sim_params.createAttribute("P", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute sliceThickness_attr = sim_params.createAttribute("s", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute zStart_attr = sim_params.createAttribute("zs", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute E0_attr = sim_params.createAttribute("E", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute alphaMax_attr = sim_params.createAttribute("A", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute rx_attr = sim_params.createAttribute("rx", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute ry_attr = sim_params.createAttribute("ry", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute df_attr = sim_params.createAttribute("df", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute C3_attr = sim_params.createAttribute("C3", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute C5_attr = sim_params.createAttribute("C5", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute semiangle_attr = sim_params.createAttribute("sa", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute detector_attr = sim_params.createAttribute("d", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute tx_attr = sim_params.createAttribute("tx", H5::PredType::NATIVE_FLOAT, scalar_attr);
	H5::Attribute ty_attr = sim_params.createAttribute("ty", H5::PredType::NATIVE_FLOAT, scalar_attr);

	//create vector spaces
	hsize_t two[1] = {2};
	hsize_t three[1] = {3};
	H5::DataSpace v_two_dataspace(1, two);
	H5::DataSpace v_three_dataspace(1, three);

	H5::Attribute cell_dim_attr = sim_params.createAttribute("c", H5::PredType::NATIVE_FLOAT, v_three_dataspace);
	H5::Attribute tile_attr = sim_params.createAttribute("t", H5::PredType::NATIVE_FLOAT, v_three_dataspace);
	H5::Attribute scanWindow_x_attr = sim_params.createAttribute("wx", H5::PredType::NATIVE_FLOAT, v_two_dataspace);
	H5::Attribute scanWindow_y_attr = sim_params.createAttribute("wy", H5::PredType::NATIVE_FLOAT, v_two_dataspace);

	H5::Attribute scanWindow_x_r_attr;
	H5::Attribute scanWindow_y_r_attr;
	if (pars.meta.realSpaceWindow_x)
		scanWindow_x_r_attr = sim_params.createAttribute("wxr", H5::PredType::NATIVE_FLOAT, v_two_dataspace);
	if (pars.meta.realSpaceWindow_y)
		scanWindow_y_r_attr = sim_params.createAttribute("wyr", H5::PredType::NATIVE_FLOAT, v_two_dataspace);

	H5::Attribute save2D_attr;
	if (pars.meta.save2DOutput)
	{
		save2D_attr = sim_params.createAttribute("2D", H5::PredType::NATIVE_FLOAT, v_two_dataspace);
	}

	//write data
	//strings
	atoms_attr.write(strdatatype, filenameAtoms);
	alg_attr.write(strdatatype, algorithm);

	//scalar logical/integers
	fx_attr.write(H5::PredType::NATIVE_INT, &pars.meta.interpolationFactorX);
	fy_attr.write(H5::PredType::NATIVE_INT, &pars.meta.interpolationFactorY);
	numFP_attr.write(H5::PredType::NATIVE_INT, &pars.meta.numFP);
	numSlices_attr.write(H5::PredType::NATIVE_INT, &pars.meta.numSlices);

	//logicals first need to be cast to ints
	int tmp_te = {pars.meta.includeThermalEffects};
	int tmp_oc = {pars.meta.includeOccupancy};
	int tmp_3D = {pars.meta.save3DOutput};
	int tmp_4D = {pars.meta.save4DOutput};
	int tmp_DPC = {pars.meta.saveDPC_CoM};
	int tmp_PS = {pars.meta.savePotentialSlices};
	int tmp_nqs = {pars.meta.nyquistSampling};

	te_attr.write(H5::PredType::NATIVE_INT, &tmp_te);
	oc_attr.write(H5::PredType::NATIVE_INT, &tmp_oc);
	save3D_attr.write(H5::PredType::NATIVE_INT, &tmp_3D);
	save4D_attr.write(H5::PredType::NATIVE_INT, &tmp_4D);
	saveDPC_attr.write(H5::PredType::NATIVE_INT, &tmp_DPC);
	savePS_attr.write(H5::PredType::NATIVE_INT, &tmp_PS);
	nyquist_attr.write(H5::PredType::NATIVE_INT, &tmp_nqs);

	//scalar floats/doubles
	px_attr.write(H5::PredType::NATIVE_FLOAT, &pars.meta.realspacePixelSize[1]);
	py_attr.write(H5::PredType::NATIVE_FLOAT, &pars.meta.realspacePixelSize[0]);
	potBound_attr.write(H5::PredType::NATIVE_FLOAT, &pars.meta.potBound);
	sliceThickness_attr.write(H5::PredType::NATIVE_FLOAT, &pars.meta.sliceThickness);
	zStart_attr.write(H5::PredType::NATIVE_FLOAT, &pars.meta.zStart);
	rx_attr.write(H5::PredType::NATIVE_FLOAT, &pars.meta.probeStepX);
	ry_attr.write(H5::PredType::NATIVE_FLOAT, &pars.meta.probeStepY);
	df_attr.write(H5::PredType::NATIVE_FLOAT, &pars.meta.probeDefocus);
	C3_attr.write(H5::PredType::NATIVE_FLOAT, &pars.meta.C3);
	C5_attr.write(H5::PredType::NATIVE_FLOAT, &pars.meta.C5);

	//scalars with unit adjustments
	PRISMATIC_FLOAT_PRECISION tmp_tx[1] = {pars.meta.probeXtilt * 1000};
	PRISMATIC_FLOAT_PRECISION tmp_ty[1] = {pars.meta.probeYtilt * 1000};
	PRISMATIC_FLOAT_PRECISION tmp_E0[1] = {pars.meta.E0 / 1000};
	PRISMATIC_FLOAT_PRECISION tmp_alphaMax[1] = {pars.meta.alphaBeamMax * 1000};
	PRISMATIC_FLOAT_PRECISION tmp_sa[1] = {pars.meta.probeSemiangle * 1000};
	PRISMATIC_FLOAT_PRECISION tmp_d[1] = {pars.meta.detectorAngleStep * 1000};

	tx_attr.write(H5::PredType::NATIVE_FLOAT, &tmp_tx);
	ty_attr.write(H5::PredType::NATIVE_FLOAT, &tmp_ty);
	E0_attr.write(H5::PredType::NATIVE_FLOAT, &tmp_E0);
	alphaMax_attr.write(H5::PredType::NATIVE_FLOAT, &tmp_alphaMax);
	semiangle_attr.write(H5::PredType::NATIVE_FLOAT, &tmp_sa);
	detector_attr.write(H5::PredType::NATIVE_FLOAT, &tmp_d);

	//vector spaces
	PRISMATIC_FLOAT_PRECISION tmp_buffer[2];

	if (pars.meta.save2DOutput)
	{
		tmp_buffer[0] = pars.meta.integrationAngleMin * 1000;
		tmp_buffer[1] = pars.meta.integrationAngleMax * 1000;
		save2D_attr.write(H5::PredType::NATIVE_FLOAT, tmp_buffer);
	}

	if (pars.meta.realSpaceWindow_x)
	{
		tmp_buffer[0] = pars.meta.scanWindowXMin_r;
		tmp_buffer[1] = pars.meta.scanWindowXMax_r;
		scanWindow_x_r_attr.write(H5::PredType::NATIVE_FLOAT, tmp_buffer);
	}

	if (pars.meta.realSpaceWindow_y)
	{
		tmp_buffer[0] = pars.meta.scanWindowYMin_r;
		tmp_buffer[1] = pars.meta.scanWindowYMax_r;
		scanWindow_y_r_attr.write(H5::PredType::NATIVE_FLOAT, tmp_buffer);
	}

	tmp_buffer[0] = pars.meta.scanWindowXMin;
	tmp_buffer[1] = pars.meta.scanWindowXMax;
	scanWindow_x_attr.write(H5::PredType::NATIVE_FLOAT, tmp_buffer);

	tmp_buffer[0] = pars.meta.scanWindowYMin;
	tmp_buffer[1] = pars.meta.scanWindowYMax;
	scanWindow_y_attr.write(H5::PredType::NATIVE_FLOAT, tmp_buffer);

	int tile_buffer[3];
	tile_buffer[0] = pars.meta.tileX;
	tile_buffer[1] = pars.meta.tileY;
	tile_buffer[2] = pars.meta.tileZ;
	tile_attr.write(H5::PredType::NATIVE_INT, tile_buffer);

	PRISMATIC_FLOAT_PRECISION cellBuffer[3];
	cellBuffer[0] = pars.meta.cellDim[0];
	cellBuffer[1] = pars.meta.cellDim[1];
	cellBuffer[2] = pars.meta.cellDim[2];
	cell_dim_attr.write(H5::PredType::NATIVE_FLOAT, cellBuffer);

	metadata.close();
};

void writeMetadata(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, double dummy)
{
	//set up group
	H5::Group metadata = pars.outputFile.openGroup("4DSTEM_simulation/metadata/metadata_0/original");
	H5::Group sim_params = metadata.createGroup("simulation_parameters");

	//write all parameters as attributes

	//create common dataspaces
	H5::DataSpace str_name_ds(H5S_SCALAR); //string dataspaces and types
	H5::StrType strdatatype(H5::PredType::C_S1, 256);
	H5::DataSpace scalar_attr(H5S_SCALAR);

	//initialize string parameter data
	H5std_string algorithm;
	if (pars.meta.algorithm == Prismatic::Algorithm::Multislice)
	{
		algorithm = "m";
	}
	else
	{
		algorithm = "p";
	}

	const H5std_string filenameAtoms(pars.meta.filenameAtoms);

	//create string attributes
	H5::Attribute atoms_attr = sim_params.createAttribute("i", strdatatype, str_name_ds);
	H5::Attribute alg_attr = sim_params.createAttribute("a", strdatatype, str_name_ds);

	//create scalar logical/integer attributes
	H5::Attribute fx_attr = sim_params.createAttribute("fx", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute fy_attr = sim_params.createAttribute("fy", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute numFP_attr = sim_params.createAttribute("F", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute numSlices_attr = sim_params.createAttribute("ns", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute te_attr = sim_params.createAttribute("te", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute oc_attr = sim_params.createAttribute("oc", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute save3D_attr = sim_params.createAttribute("3D", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute save4D_attr = sim_params.createAttribute("4D", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute saveDPC_attr = sim_params.createAttribute("DPC", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute savePS_attr = sim_params.createAttribute("ps", H5::PredType::NATIVE_INT, scalar_attr);
	H5::Attribute nyquist_attr = sim_params.createAttribute("nqs", H5::PredType::NATIVE_INT, scalar_attr);

	//create scalar float/double attributes (changes based on prismatic float precision)
	H5::Attribute px_attr = sim_params.createAttribute("px", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute py_attr = sim_params.createAttribute("py", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute potBound_attr = sim_params.createAttribute("P", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute sliceThickness_attr = sim_params.createAttribute("s", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute zStart_attr = sim_params.createAttribute("zs", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute E0_attr = sim_params.createAttribute("E", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute alphaMax_attr = sim_params.createAttribute("A", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute rx_attr = sim_params.createAttribute("rx", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute ry_attr = sim_params.createAttribute("ry", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute df_attr = sim_params.createAttribute("df", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute C3_attr = sim_params.createAttribute("C3", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute C5_attr = sim_params.createAttribute("C5", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute semiangle_attr = sim_params.createAttribute("sa", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute detector_attr = sim_params.createAttribute("d", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute tx_attr = sim_params.createAttribute("tx", H5::PredType::NATIVE_DOUBLE, scalar_attr);
	H5::Attribute ty_attr = sim_params.createAttribute("ty", H5::PredType::NATIVE_DOUBLE, scalar_attr);

	//create vector spaces
	hsize_t two[1] = {2};
	hsize_t three[1] = {3};
	H5::DataSpace v_two_dataspace(1, two);
	H5::DataSpace v_three_dataspace(1, three);

	H5::Attribute cell_dim_attr = sim_params.createAttribute("c", H5::PredType::NATIVE_DOUBLE, v_three_dataspace);
	H5::Attribute tile_attr = sim_params.createAttribute("t", H5::PredType::NATIVE_DOUBLE, v_three_dataspace);
	H5::Attribute scanWindow_x_attr = sim_params.createAttribute("wx", H5::PredType::NATIVE_DOUBLE, v_two_dataspace);
	H5::Attribute scanWindow_y_attr = sim_params.createAttribute("wy", H5::PredType::NATIVE_DOUBLE, v_two_dataspace);

	H5::Attribute scanWindow_x_r_attr;
	H5::Attribute scanWindow_y_r_attr;
	if (pars.meta.realSpaceWindow_x)
		scanWindow_x_r_attr = sim_params.createAttribute("wxr", H5::PredType::NATIVE_DOUBLE, v_two_dataspace);
	if (pars.meta.realSpaceWindow_y)
		scanWindow_y_r_attr = sim_params.createAttribute("wyr", H5::PredType::NATIVE_DOUBLE, v_two_dataspace);

	H5::Attribute save2D_attr;
	if (pars.meta.save2DOutput)
	{
		save2D_attr = sim_params.createAttribute("2D", H5::PredType::NATIVE_DOUBLE, v_two_dataspace);
	}

	//write data
	//strings
	atoms_attr.write(strdatatype, filenameAtoms);
	alg_attr.write(strdatatype, algorithm);

	//scalar integers
	fx_attr.write(H5::PredType::NATIVE_INT, &pars.meta.interpolationFactorX);
	fy_attr.write(H5::PredType::NATIVE_INT, &pars.meta.interpolationFactorY);
	numFP_attr.write(H5::PredType::NATIVE_INT, &pars.meta.numFP);
	numSlices_attr.write(H5::PredType::NATIVE_INT, &pars.meta.numSlices);

	//logicals first need to be cast to ints
	int tmp_te = {pars.meta.includeThermalEffects};
	int tmp_oc = {pars.meta.includeOccupancy};
	int tmp_3D = {pars.meta.save3DOutput};
	int tmp_4D = {pars.meta.save4DOutput};
	int tmp_DPC = {pars.meta.saveDPC_CoM};
	int tmp_PS = {pars.meta.savePotentialSlices};
	int tmp_nqs = {pars.meta.nyquistSampling};

	te_attr.write(H5::PredType::NATIVE_INT, &tmp_te);
	oc_attr.write(H5::PredType::NATIVE_INT, &tmp_oc);
	save3D_attr.write(H5::PredType::NATIVE_INT, &tmp_3D);
	save4D_attr.write(H5::PredType::NATIVE_INT, &tmp_4D);
	saveDPC_attr.write(H5::PredType::NATIVE_INT, &tmp_DPC);
	savePS_attr.write(H5::PredType::NATIVE_INT, &tmp_PS);
	nyquist_attr.write(H5::PredType::NATIVE_INT, &tmp_nqs);

	//scalar floats/doubles
	px_attr.write(H5::PredType::NATIVE_DOUBLE, &pars.meta.realspacePixelSize[1]);
	py_attr.write(H5::PredType::NATIVE_DOUBLE, &pars.meta.realspacePixelSize[0]);
	potBound_attr.write(H5::PredType::NATIVE_DOUBLE, &pars.meta.potBound);
	sliceThickness_attr.write(H5::PredType::NATIVE_DOUBLE, &pars.meta.sliceThickness);
	zStart_attr.write(H5::PredType::NATIVE_DOUBLE, &pars.meta.zStart);
	rx_attr.write(H5::PredType::NATIVE_DOUBLE, &pars.meta.probeStepX);
	ry_attr.write(H5::PredType::NATIVE_DOUBLE, &pars.meta.probeStepY);
	df_attr.write(H5::PredType::NATIVE_DOUBLE, &pars.meta.probeDefocus);
	C3_attr.write(H5::PredType::NATIVE_DOUBLE, &pars.meta.C3);
	C5_attr.write(H5::PredType::NATIVE_DOUBLE, &pars.meta.C5);

	//scalars with unit adjustments
	PRISMATIC_FLOAT_PRECISION tmp_tx[1] = {pars.meta.probeXtilt * 1000};
	PRISMATIC_FLOAT_PRECISION tmp_ty[1] = {pars.meta.probeYtilt * 1000};
	PRISMATIC_FLOAT_PRECISION tmp_E0[1] = {pars.meta.E0 / 1000};
	PRISMATIC_FLOAT_PRECISION tmp_alphaMax[1] = {pars.meta.alphaBeamMax * 1000};
	PRISMATIC_FLOAT_PRECISION tmp_sa[1] = {pars.meta.probeSemiangle * 1000};
	PRISMATIC_FLOAT_PRECISION tmp_d[1] = {pars.meta.detectorAngleStep * 1000};

	tx_attr.write(H5::PredType::NATIVE_DOUBLE, &tmp_tx);
	ty_attr.write(H5::PredType::NATIVE_DOUBLE, &tmp_ty);
	E0_attr.write(H5::PredType::NATIVE_DOUBLE, &tmp_E0);
	alphaMax_attr.write(H5::PredType::NATIVE_DOUBLE, &tmp_alphaMax);
	semiangle_attr.write(H5::PredType::NATIVE_DOUBLE, &tmp_sa);
	detector_attr.write(H5::PredType::NATIVE_DOUBLE, &tmp_d);

	//vector spaces
	PRISMATIC_FLOAT_PRECISION tmp_buffer[2];

	if (pars.meta.save2DOutput)
	{
		tmp_buffer[0] = pars.meta.integrationAngleMin * 1000;
		tmp_buffer[1] = pars.meta.integrationAngleMax * 1000;
		save2D_attr.write(H5::PredType::NATIVE_DOUBLE, tmp_buffer);
	}

	if (pars.meta.realSpaceWindow_x)
	{
		tmp_buffer[0] = pars.meta.scanWindowXMin_r;
		tmp_buffer[1] = pars.meta.scanWindowXMax_r;
		scanWindow_x_r_attr.write(H5::PredType::NATIVE_DOUBLE, tmp_buffer);
	}

	if (pars.meta.realSpaceWindow_y)
	{
		tmp_buffer[0] = pars.meta.scanWindowYMin_r;
		tmp_buffer[1] = pars.meta.scanWindowYMax_r;
		scanWindow_y_r_attr.write(H5::PredType::NATIVE_DOUBLE, tmp_buffer);
	}

	tmp_buffer[0] = pars.meta.scanWindowXMin;
	tmp_buffer[1] = pars.meta.scanWindowXMax;
	scanWindow_x_attr.write(H5::PredType::NATIVE_DOUBLE, tmp_buffer);

	tmp_buffer[0] = pars.meta.scanWindowYMin;
	tmp_buffer[1] = pars.meta.scanWindowYMax;
	scanWindow_y_attr.write(H5::PredType::NATIVE_DOUBLE, tmp_buffer);

	int tile_buffer[3];
	tile_buffer[0] = pars.meta.tileX;
	tile_buffer[1] = pars.meta.tileY;
	tile_buffer[2] = pars.meta.tileZ;
	tile_attr.write(H5::PredType::NATIVE_INT, tile_buffer);

	PRISMATIC_FLOAT_PRECISION cellBuffer[3];
	cellBuffer[0] = pars.meta.cellDim[0];
	cellBuffer[1] = pars.meta.cellDim[1];
	cellBuffer[2] = pars.meta.cellDim[2];
	cell_dim_attr.write(H5::PredType::NATIVE_DOUBLE, cellBuffer);

	metadata.close();
};
} // namespace Prismatic
