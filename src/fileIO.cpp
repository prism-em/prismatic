#include "H5Cpp.h"
#include "params.h"
#include "fileIO.h"
#include "utility.h"
#include <mutex>

namespace Prismatic{

void setupOutputFile(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
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
	H5::Group supergroups(data.createGroup("supergroups"));

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

void setup4DOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	std::cout << "here?0" << std::endl;
	H5::Group datacubes = pars.outputFile.openGroup("4DSTEM_simulation/data/datacubes");
	std::cout << "here?0.5" << std::endl;

	//shared properties
	std::string base_name = "CBED_array_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[4];
	data_dims[0] = {pars.numXprobes};
	data_dims[1] = {pars.numYprobes};
	hsize_t chunkDims[4];
	chunkDims[0] = chunkDims[1] = {1};
	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};
	hsize_t qx_dim[1];
	hsize_t qy_dim[1];

	Array1D<PRISMATIC_FLOAT_PRECISION> qx;
	Array1D<PRISMATIC_FLOAT_PRECISION> qy;
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
        offset_qx = 0;
        offset_qy = 0;
    }
    else
    {
        if (pars.meta.algorithm == Algorithm::Multislice)
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

	if (pars.meta.algorithm == Algorithm::Multislice)
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

	H5::CompType complex_type = H5::CompType(sizeof(complex_float_t));
	const H5std_string re_str("r"); //using h5py default configuration
	const H5std_string im_str("i");
	complex_type.insertMember(re_str, 0, PFP_TYPE);
	complex_type.insertMember(im_str, 4, PFP_TYPE);

	std::cout << "here?1" << std::endl;
	for (auto n = 0; n < pars.numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n) + pars.currentTag;
		H5::Group CBED_slice_n(datacubes.createGroup(nth_name.c_str()));
		std::cout << "here?2" << std::endl;

		//write attributes
		writeScalarAttribute(CBED_slice_n, "emd_group_type", 1);
		writeScalarAttribute(CBED_slice_n, "metadata", 0);
		writeScalarAttribute(CBED_slice_n, "output_depth", pars.depths[n]);
		
		//setup data set chunking properties
		H5::DSetCreatPropList plist;
		plist.setChunk(4, chunkDims);

		//create dataset
		H5::DataSpace mspace(4, data_dims); //rank is 4
		H5::DataSet CBED_data;
		if(pars.meta.saveComplexOutputWave)
		{
			CBED_data = CBED_slice_n.createDataSet("datacube", complex_type, mspace, plist);
		}
		else
		{
			CBED_data = CBED_slice_n.createDataSet("datacube", PFP_TYPE, mspace, plist);
		}
		mspace.close();
		std::cout << "here?3" << std::endl;

		//write dimensions
		H5::DataSpace str_name_ds(H5S_SCALAR);
		H5::StrType strdatatype(H5::PredType::C_S1, 256);

		std::vector<size_t> order = {0};
		writeRealDataSet(CBED_slice_n, "dim1", &pars.xp[0], rx_dim, 1, order);
		writeRealDataSet(CBED_slice_n, "dim2", &pars.yp[0], ry_dim, 1, order);
		writeRealDataSet(CBED_slice_n, "dim3", &qx[offset_qx], qx_dim, 1, order);
		writeRealDataSet(CBED_slice_n, "dim4", &qy[offset_qy], qy_dim, 1, order);

		//dimension attributes
		H5::DataSet dim1 = CBED_slice_n.openDataSet("dim1");
		H5::DataSet dim2 = CBED_slice_n.openDataSet("dim2");
		H5::DataSet dim3 = CBED_slice_n.openDataSet("dim3");
		H5::DataSet dim4 = CBED_slice_n.openDataSet("dim4");

		writeScalarAttribute(dim1, "name", "R_x");
		writeScalarAttribute(dim2, "name", "R_y");
		writeScalarAttribute(dim3, "name", "Q_x");
		writeScalarAttribute(dim4, "name", "Q_y");

		writeScalarAttribute(dim1, "units", "[Å]");
		writeScalarAttribute(dim2, "units", "[Å]");
		writeScalarAttribute(dim3, "units", "[Å^-1]");
		writeScalarAttribute(dim4, "units", "[Å^-1]");

		dim1.close();
		dim2.close();
		dim3.close();
		dim4.close();
		CBED_slice_n.close();
	}

	datacubes.close();
};

void setupVDOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	//shared properties
	std::string base_name = "virtual_detector_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[3];
	data_dims[0] = {pars.numXprobes};
	data_dims[1] = {pars.numYprobes};
	data_dims[2] = {pars.Ndet};

	hsize_t chunkDims[3] = {1, pars.numYprobes, pars.Ndet};

	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};
	hsize_t bin_dim[1] = {pars.Ndet};

	
	H5::CompType complex_type = H5::CompType(sizeof(complex_float_t));
	const H5std_string re_str("r"); //using h5py default configuration
	const H5std_string im_str("i");
	complex_type.insertMember(re_str, 0, PFP_TYPE);
	complex_type.insertMember(im_str, 4, PFP_TYPE);

	for (auto n = 0; n < pars.numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n) + pars.currentTag;
		std::cout << "here" << std::endl;
		H5::Group VD_slice_n(realslices.createGroup(nth_name.c_str()));
		std::cout << "here 2" << std::endl;

		//write attributes
		writeScalarAttribute(VD_slice_n, "emd_group_type", 1);
		writeScalarAttribute(VD_slice_n, "metadata", 0);
		writeScalarAttribute(VD_slice_n, "output_depth", pars.depths[n]);
		if(pars.meta.simSeries) writeScalarAttribute(VD_slice_n, "output_defocus", pars.meta.probeDefocus);

		//create datasets
		H5::DataSpace mspace(3, data_dims); //rank is 2 for each realslice

		H5::DataSet VD_data;
		H5::DSetCreatPropList plist;
		plist.setChunk(3, chunkDims);
		if(pars.meta.saveComplexOutputWave)
		{
			VD_slice_n.createDataSet("realslice", complex_type, mspace, plist);
		}
		else
		{
			VD_slice_n.createDataSet("realslice", PFP_TYPE, mspace, plist);
		}

		VD_data.close();
		mspace.close();

		//write dimensions
		std::vector<size_t> order = {0};
		writeRealDataSet(VD_slice_n, "dim1", &pars.xp[0], rx_dim, 1, order);
		writeRealDataSet(VD_slice_n, "dim2", &pars.yp[0], ry_dim, 1, order);
		writeRealDataSet(VD_slice_n, "dim3", &pars.detectorAngles[0], bin_dim, 1, order);

		//dimension attribute
		H5::DataSet dim1 = VD_slice_n.openDataSet("dim1");
		H5::DataSet dim2 = VD_slice_n.openDataSet("dim2");
		H5::DataSet dim3 = VD_slice_n.openDataSet("dim3");

		writeScalarAttribute(dim1, "name", "R_x");
		writeScalarAttribute(dim2, "name", "R_y");
		writeScalarAttribute(dim3, "name", "bin_outer_angle");

		writeScalarAttribute(dim1, "units", "[Å]");
		writeScalarAttribute(dim2, "units", "[Å]");
		writeScalarAttribute(dim3, "units", "[mrad]");

		dim1.close();
		dim2.close();
		dim3.close();
		VD_slice_n.close();
	}

	realslices.close();
};

void setup2DOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	//shared properties
	std::string base_name = "annular_detector_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[2];
	data_dims[0] = {pars.numXprobes};
	data_dims[1] = {pars.numYprobes};

	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};

	H5::CompType complex_type = H5::CompType(sizeof(complex_float_t));
	const H5std_string re_str("r"); //using h5py default configuration
	const H5std_string im_str("i");
	complex_type.insertMember(re_str, 0, PFP_TYPE);
	complex_type.insertMember(im_str, 4, PFP_TYPE);

	for (auto n = 0; n < pars.numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n) + pars.currentTag;
		H5::Group annular_slice_n(realslices.createGroup(nth_name.c_str()));

		//write attributes
		writeScalarAttribute(annular_slice_n, "emd_group_type", 1);
		writeScalarAttribute(annular_slice_n, "metadata", 0);
		writeScalarAttribute(annular_slice_n, "output_depth", pars.depths[n]);
		writeScalarAttribute(annular_slice_n, "depth", 1);

		//create dataset
		H5::DataSpace mspace(2, data_dims); //rank is 2
		H5::DataSet annular_data;
		if(pars.meta.saveComplexOutputWave)
		{
			annular_data = annular_slice_n.createDataSet("realslice", complex_type, mspace);
		}
		else
		{
			annular_data = annular_slice_n.createDataSet("realslice", PFP_TYPE, mspace);
		}
		mspace.close();

		//write dimensions
		std::vector<size_t> order = {0};
		writeRealDataSet(annular_slice_n, "dim1", &pars.xp[0], rx_dim, 1, order);
		writeRealDataSet(annular_slice_n, "dim2", &pars.yp[0], ry_dim, 1, order);

		//dimension attribute
		H5::DataSet dim1 = annular_slice_n.openDataSet("dim1");
		H5::DataSet dim2 = annular_slice_n.openDataSet("dim2");

		writeScalarAttribute(dim1, "name", "R_x");
		writeScalarAttribute(dim2, "name", "R_y");

		writeScalarAttribute(dim1, "units", "[Å]");
		writeScalarAttribute(dim2, "units", "[Å]");

		dim1.close();
		dim2.close();

		annular_slice_n.close();
	}

	realslices.close();
};

void setupDPCOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	//shared properties
	std::string base_name = "DPC_CoM_depth";
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[3];
	data_dims[0] = {pars.numXprobes};
	data_dims[1] = {pars.numYprobes};
	data_dims[2] = {2};

	hsize_t rx_dim[1] = {pars.xp.size()};
	hsize_t ry_dim[1] = {pars.yp.size()};
	hsize_t str_dim[1] = {2};

	for (auto n = 0; n < pars.numLayers; n++)
	{
		//create slice group
		std::string nth_name = base_name + getDigitString(n) + pars.currentTag;
		H5::Group DPC_CoM_slice_n(realslices.createGroup(nth_name.c_str()));

		//write attributes
		writeScalarAttribute(DPC_CoM_slice_n, "emd_group_type", 1);
		writeScalarAttribute(DPC_CoM_slice_n, "metadata", 0);
		writeScalarAttribute(DPC_CoM_slice_n, "output_depth", pars.depths[n]);

		//create dataset
		H5::DataSpace mspace(3, data_dims); //rank is 3
		H5::DataSet DPC_data = DPC_CoM_slice_n.createDataSet("realslice", PFP_TYPE, mspace);
		mspace.close();

		//write dimensions
		std::vector<size_t> order = {0};
		writeRealDataSet(DPC_CoM_slice_n, "dim1", &pars.xp[0], rx_dim, 1, order);
		writeRealDataSet(DPC_CoM_slice_n, "dim2", &pars.yp[0], ry_dim, 1, order);

		H5::StrType strdatatype(H5::PredType::C_S1, 256);
		H5::DataSpace dim3_mspace(1, str_dim);
		H5::DataSet dim3 = DPC_CoM_slice_n.createDataSet("dim3", strdatatype, dim3_mspace);
		H5std_string dpc_x("DPC_CoM_x");
		H5std_string dpc_y("DPC_CoM_y");
		H5std_string str_buffer_array[2] = {dpc_x, dpc_y};
		writeStringArray(dim3, str_buffer_array, 2);

		//dimension attribute
		H5::DataSet dim1 = DPC_CoM_slice_n.openDataSet("dim1");
		H5::DataSet dim2 = DPC_CoM_slice_n.openDataSet("dim2");

		writeScalarAttribute(dim1, "name", "R_x");
		writeScalarAttribute(dim2, "name", "R_y");
		writeScalarAttribute(dim3, "name", "X/Y");

		writeScalarAttribute(dim1, "units", "[Å]");
		writeScalarAttribute(dim2, "units", "[Å]");
		writeScalarAttribute(dim3, "units", "[none]");

		dim1.close();
		dim2.close();
		dim3.close();
		DPC_CoM_slice_n.close();
	}

	realslices.close();
};

void setupSMatrixOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, const int FP)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	std::string base_name = "smatrix_fp" + getDigitString(FP);
	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[3] = {pars.Scompact.get_dimi(), pars.Scompact.get_dimj(), pars.Scompact.get_dimk()};

	hsize_t x_size[1] = {pars.imageSize[1] / 2};
	hsize_t y_size[1] = {pars.imageSize[0] / 2};
	hsize_t beams[1] = {pars.numberBeams};

	H5::CompType complex_type = H5::CompType(sizeof(complex_float_t));
	const H5std_string re_str("r"); //using h5py default configuration
	const H5std_string im_str("i");
	complex_type.insertMember(re_str, 0, PFP_TYPE);
	complex_type.insertMember(im_str, 4, PFP_TYPE);

	H5::Group smatrix_group(realslices.createGroup(base_name.c_str()));

	//write attributes
	writeScalarAttribute(smatrix_group, "emd_group_type", 1);
	writeScalarAttribute(smatrix_group, "metadata", 0);

	//create datasets
	H5::DataSpace mspace(3, data_dims); //rank is 2 for each realslice
	H5::DataSet smatrix_data = smatrix_group.createDataSet("realslice", complex_type, mspace);
	smatrix_data.close();
	mspace.close();

	//write dimensions
	Array1D<PRISMATIC_FLOAT_PRECISION> x_dim_data = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.imageSize[1] / 2}});
	Array1D<PRISMATIC_FLOAT_PRECISION> y_dim_data = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.imageSize[0] / 2}});
	for (auto i = 0; i < pars.imageSize[1] / 2; i++) x_dim_data[i] = i * pars.pixelSize[1]*2;
	for (auto i = 0; i < pars.imageSize[0] / 2; i++) y_dim_data[i] = i * pars.pixelSize[0]*2;

	std::vector<PRISMATIC_FLOAT_PRECISION> beamsIndex(pars.numberBeams); //convert to float
	for(auto i = 0; i < pars.numberBeams; i++) beamsIndex.push_back(pars.beamsIndex[i]);

	std::vector<size_t> order = {0};
	writeRealDataSet(smatrix_group, "dim1", &x_dim_data[0], x_size, 1, order);
	writeRealDataSet(smatrix_group, "dim2", &y_dim_data[0], y_size, 1, order);
	writeRealDataSet(smatrix_group, "dim3", &beamsIndex[0], beams, 1, order);

	//dimension attributes
	H5::DataSet dim1 = smatrix_group.openDataSet("dim1");
	H5::DataSet dim2 = smatrix_group.openDataSet("dim2");
	H5::DataSet dim3 = smatrix_group.openDataSet("dim3");

	writeScalarAttribute(dim1, "name", "R_x");
	writeScalarAttribute(dim2, "name", "R_y");
	writeScalarAttribute(dim3, "name", "beam_number");

	writeScalarAttribute(dim1, "units", "[Å]");
	writeScalarAttribute(dim2, "units", "[Å]");
	writeScalarAttribute(dim3, "units", "[none]");

	smatrix_group.close();
	realslices.close();

};

void setupHRTEMOutput(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");

	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[3] = {pars.Scompact.get_dimi(), pars.Scompact.get_dimj(),
						    pars.Scompact.get_dimk()};

	hsize_t x_size[1] = {pars.imageSize[1] / 2};
	hsize_t y_size[1] = {pars.imageSize[0] / 2};
	hsize_t tilt_size[2] = {pars.Scompact.get_dimk(), 2};

	H5::CompType complex_type = H5::CompType(sizeof(complex_float_t));
	const H5std_string re_str("r"); //using h5py default configuration
	const H5std_string im_str("i");
	complex_type.insertMember(re_str, 0, PFP_TYPE);
	complex_type.insertMember(im_str, 4, PFP_TYPE);

	H5::Group hrtem_group(realslices.createGroup("HRTEM"));

	//write attributes
	writeScalarAttribute(hrtem_group, "emd_group_type", 1);
	writeScalarAttribute(hrtem_group, "metadata", 0);

	//create datasets
	H5::DataSpace mspace(3, data_dims); //rank is 2 for each realslice
	H5::DataSet hrtem_data;
	if(pars.meta.saveComplexOutputWave)
	{
		hrtem_data = hrtem_group.createDataSet("realslice", complex_type, mspace);
	}
	else
	{
		hrtem_data = hrtem_group.createDataSet("realslice", PFP_TYPE, mspace);
	}
	
	hrtem_data.close();
	mspace.close();

	//write dimensions
	Array1D<PRISMATIC_FLOAT_PRECISION> x_dim_data = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.imageSize[1] / 2}});
	Array1D<PRISMATIC_FLOAT_PRECISION> y_dim_data = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.imageSize[0] / 2}});
	for (auto i = 0; i < pars.imageSize[1] / 2; i++) x_dim_data[i] = i * pars.pixelSize[1]*2;
	for (auto i = 0; i < pars.imageSize[0] / 2; i++) y_dim_data[i] = i * pars.pixelSize[0]*2;

	std::vector<PRISMATIC_FLOAT_PRECISION> xTilts_write(pars.xTilts_tem);
	std::vector<PRISMATIC_FLOAT_PRECISION> yTilts_write(pars.yTilts_tem);

	for(auto i = 0; i < xTilts_write.size(); i++) xTilts_write[i] *= 1000; //convert to mrad for writing
	for(auto i = 0; i < yTilts_write.size(); i++) yTilts_write[i] *= 1000; //convert to mrad for writing
	Array2D<PRISMATIC_FLOAT_PRECISION> dim_tilts = zeros_ND<2,PRISMATIC_FLOAT_PRECISION>({{tilt_size[0], tilt_size[1]}});
	for(auto i = 0; i < xTilts_write.size(); i++)
	{
		dim_tilts.at(i,0) = xTilts_write[i];
		dim_tilts.at(i,1) = yTilts_write[i];
	}
	
	std::vector<size_t> order_1D = {0};
	std::vector<size_t> order_2D = {0, 1};

	writeRealDataSet(hrtem_group, "dim1", &x_dim_data[0], x_size, 1, order_1D);
	writeRealDataSet(hrtem_group, "dim2", &y_dim_data[0], y_size, 1, order_1D);
	writeRealDataSet(hrtem_group, "dim3", &dim_tilts[0], tilt_size, 2, order_2D);

	//dimension attributes
	H5::DataSet dim1 = hrtem_group.openDataSet("dim1");
	H5::DataSet dim2 = hrtem_group.openDataSet("dim2");
	H5::DataSet dim3 = hrtem_group.openDataSet("dim3");

	writeScalarAttribute(dim1, "name", "R_x");
	writeScalarAttribute(dim2, "name", "R_y");
	writeScalarAttribute(dim3, "name", "Tilts");

	writeScalarAttribute(dim1, "units", "[Å]");
	writeScalarAttribute(dim2, "units", "[Å]");
	writeScalarAttribute(dim3, "units", "[mrad]");

	hrtem_group.close();
	realslices.close();

};

void setupHRTEMOutput_virtual(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	H5::Group datacubes = pars.outputFile.openGroup("4DSTEM_simulation/data/datacubes");

	//get unique vectors for qx and qy first; assumes sorting before set up
	std::vector<PRISMATIC_FLOAT_PRECISION> xTilts_unique = getUnique(pars.xTilts_tem);

	std::vector<PRISMATIC_FLOAT_PRECISION> yTilts_unique(pars.yTilts_tem);
	std::sort(yTilts_unique.begin(), yTilts_unique.end());
	yTilts_unique = getUnique(yTilts_unique);
	std::cout << "yTilts_unique size: " << yTilts_unique.size() << std::endl;

	hsize_t attr_dims[1] = {1};
	hsize_t data_dims[4] = {pars.Scompact.get_dimi(), pars.Scompact.get_dimj(),
						    xTilts_unique.size(), yTilts_unique.size()};

	hsize_t x_size[1] = {pars.imageSize[1] / 2};
	hsize_t y_size[1] = {pars.imageSize[0] / 2};
	hsize_t tiltX_size[1] = {xTilts_unique.size()};
	hsize_t tiltY_size[1] = {yTilts_unique.size()};

	H5::Group hrtem_group(datacubes.createGroup("HRTEM_virtual"));

	//write attributes
	writeScalarAttribute(hrtem_group, "emd_group_type", 1);
	writeScalarAttribute(hrtem_group, "metadata", 0);

	//create virtual dataset
	H5::DataSet src_data = pars.outputFile.openDataSet("4DSTEM_simulation/data/realslices/HRTEM/realslice");
	std::string src_path = src_data.getObjName();
	H5::DataSpace src_mspace = src_data.getSpace();

	H5::DataSpace vds_mspace(4, data_dims);
	H5::DSetCreatPropList plist;
	hsize_t dest_offset[4] = {0,0,0,0};
	hsize_t dest_mdims[4] = {data_dims[0], data_dims[1], 1, 1};
	hsize_t src_offset[3] = {0,0,0};
	hsize_t src_mdims[3] = {data_dims[0], data_dims[1], 1};

	for(auto i = 0; i < pars.xTiltsInd_tem.size(); i++)
	{
		src_offset[2] = i; //since the dimensions are alerady sorted, this does not need to be calculated
		src_mspace.selectHyperslab(H5S_SELECT_SET, src_mdims, src_offset);

		dest_offset[2] = pars.xTiltsInd_tem[i];
		dest_offset[3] = pars.yTiltsInd_tem[i];
		vds_mspace.selectHyperslab(H5S_SELECT_SET, dest_mdims, dest_offset);

		plist.setVirtual(vds_mspace, src_data.getFileName(), src_path, src_mspace);
	}

	dest_offset[2] = 0;
	dest_offset[3] = 0;
	vds_mspace.selectHyperslab(H5S_SELECT_SET, data_dims, dest_offset);
	if(pars.meta.saveComplexOutputWave)
	{
		std::complex<PRISMATIC_FLOAT_PRECISION> fillVal = {nan(""), nan("")};
		plist.setFillValue(src_data.getDataType(), &fillVal);
	}
	else
	{
		PRISMATIC_FLOAT_PRECISION fillVal = nan("");
		plist.setFillValue(src_data.getDataType(), &fillVal);
	}
	H5::DataSet hrtem_data = hrtem_group.createDataSet("datacube", src_data.getDataType(), vds_mspace, plist);
	
	src_data.close();
	src_mspace.close();
	hrtem_data.close();
	vds_mspace.close();

	//write dimensions
	Array1D<PRISMATIC_FLOAT_PRECISION> x_dim_data = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.imageSize[1] / 2}});
	Array1D<PRISMATIC_FLOAT_PRECISION> y_dim_data = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.imageSize[0] / 2}});
	for (auto i = 0; i < pars.imageSize[1] / 2; i++) x_dim_data[i] = i * pars.pixelSize[1]*2;
	for (auto i = 0; i < pars.imageSize[0] / 2; i++) y_dim_data[i] = i * pars.pixelSize[0]*2;

	std::vector<PRISMATIC_FLOAT_PRECISION> xTilts_write(xTilts_unique);
	std::vector<PRISMATIC_FLOAT_PRECISION> yTilts_write(yTilts_unique);

	for(auto i = 0; i < xTilts_write.size(); i++) xTilts_write[i] *= 1000; //convert to mrad for writing
	for(auto i = 0; i < yTilts_write.size(); i++) yTilts_write[i] *= 1000; //convert to mrad for writing
	
	std::vector<size_t> order = {0};
	writeRealDataSet(hrtem_group, "dim1", &x_dim_data[0], x_size, 1, order);
	writeRealDataSet(hrtem_group, "dim2", &y_dim_data[0], y_size, 1, order);
	writeRealDataSet(hrtem_group, "dim3", &xTilts_write[0], tiltX_size, 1, order);
	writeRealDataSet(hrtem_group, "dim4", &yTilts_write[0], tiltY_size, 1, order);

	//dimension attributes
	H5::DataSet dim1 = hrtem_group.openDataSet("dim1");
	H5::DataSet dim2 = hrtem_group.openDataSet("dim2");
	H5::DataSet dim3 = hrtem_group.openDataSet("dim3");
	H5::DataSet dim4 = hrtem_group.openDataSet("dim4");

	writeScalarAttribute(dim1, "name", "R_x");
	writeScalarAttribute(dim2, "name", "R_y");
	writeScalarAttribute(dim3, "name", "Tilt_x");
	writeScalarAttribute(dim4, "name", "Tilt_y");

	writeScalarAttribute(dim1, "units", "[Å]");
	writeScalarAttribute(dim2, "units", "[Å]");
	writeScalarAttribute(dim3, "units", "[mrad]");
	writeScalarAttribute(dim4, "units", "[mrad]");

	hrtem_group.close();
	datacubes.close();
};

void sortHRTEMbeams(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	//convert beam indices into qx, qy indices by sorted tilts
	size_t N = pars.numberBeams;
	std::vector<std::pair<PRISMATIC_FLOAT_PRECISION, PRISMATIC_FLOAT_PRECISION>> tilts(N);
	for(auto i = 0; i < N; i++) tilts[i] = std::make_pair(pars.xTilts_tem[i], pars.yTilts_tem[i]);
	std::vector<size_t> indices(N);
	for(auto i = 0; i < N; i++) indices[i] = i;
    std::sort(indices.begin(), indices.end(), [&](int i, int j){return tilts[i]<tilts[j];} );

	pars.HRTEMbeamOrder = indices;
	
	int minXtiltInd = *std::min_element(pars.xTiltsInd_tem.begin(), pars.xTiltsInd_tem.end());
	int minYtiltInd = *std::min_element(pars.yTiltsInd_tem.begin(), pars.yTiltsInd_tem.end());

	//sort tilt arrays and generate 2D array coords for each beam
	std::vector<PRISMATIC_FLOAT_PRECISION> xTilts_tmp(pars.xTilts_tem);
	std::vector<PRISMATIC_FLOAT_PRECISION> yTilts_tmp(pars.yTilts_tem);
	std::vector<int> xTiltsInd_tmp(pars.xTiltsInd_tem);
	std::vector<int> yTiltsInd_tmp(pars.yTiltsInd_tem);
	for(auto i = 0; i < N; i++)
	{
		pars.xTilts_tem[i] = xTilts_tmp[indices[i]];
		pars.yTilts_tem[i] = yTilts_tmp[indices[i]];
		pars.xTiltsInd_tem[i] = xTiltsInd_tmp[indices[i]]-minXtiltInd;
		pars.yTiltsInd_tem[i] = yTiltsInd_tmp[indices[i]]-minYtiltInd;
	}
}

//these write functions will soon be deprecated
void writeRealSlice(H5::DataSet dataset, const PRISMATIC_FLOAT_PRECISION *buffer, const hsize_t *mdims)
{
	H5::DataSpace fspace = dataset.getSpace(); //all realslices have data written all at once
	H5::DataSpace mspace(2, mdims);			   //rank = 2

	dataset.write(buffer, PFP_TYPE, mspace, fspace);

	fspace.close();
	mspace.close();
}

void writeDatacube3D(H5::DataSet dataset, const PRISMATIC_FLOAT_PRECISION *buffer, const hsize_t *mdims)
{
	//set up file and memory spaces
	H5::DataSpace fspace = dataset.getSpace(); //all 3D cubes will write full buffer at once
	H5::DataSpace mspace(3, mdims);			   //rank = 3

	dataset.write(buffer, PFP_TYPE, mspace, fspace);

	fspace.close();
	mspace.close();
};

void writeDatacube3D(H5::DataSet dataset, const std::complex<PRISMATIC_FLOAT_PRECISION> *buffer, const hsize_t *mdims)
{
	//set up file and memory spaces
	H5::DataSpace fspace = dataset.getSpace(); //all 3D cubes will write full buffer at once
	H5::DataSpace mspace(3, mdims);			   //rank = 3

	dataset.write(buffer, dataset.getDataType(), mspace, fspace);

	fspace.close();
	mspace.close();
};

void writeStringArray(H5::DataSet dataset, H5std_string *string_array, const hsize_t elements)
{
	//assumes that we are writing a 1 dimensional array of strings- used only for DPC
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

void savePotentialSlices(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	// create new datacube group
	H5::Group realslices = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices");
	std::stringstream nameString;
	nameString << "ppotential_fp" << getDigitString(pars.fpFlag);
	std::string groupName = nameString.str();
	H5::Group ppotential;

	ppotential = realslices.createGroup(groupName);

	//write attributes
	writeScalarAttribute(ppotential, "emd_group_type", 1);
	writeScalarAttribute(ppotential, "metadata", 0);

	//write dimensions
	hsize_t x_size[1] = {pars.imageSize[1]};
	hsize_t y_size[1] = {pars.imageSize[0]};
	hsize_t z_size[1] = {pars.numPlanes};

	Array1D<PRISMATIC_FLOAT_PRECISION> x_dim_data = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.imageSize[1]}});
	Array1D<PRISMATIC_FLOAT_PRECISION> y_dim_data = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.imageSize[0]}});
	Array1D<PRISMATIC_FLOAT_PRECISION> z_dim_data = zeros_ND<1, PRISMATIC_FLOAT_PRECISION>({{pars.numPlanes}});

	for (auto i = 0; i < pars.imageSize[1]; i++) x_dim_data[i] = i * pars.pixelSize[1];
	for (auto i = 0; i < pars.imageSize[0]; i++) y_dim_data[i] = i * pars.pixelSize[0];
	for (auto i = 0; i < pars.numPlanes; i++) z_dim_data[i] = i * pars.meta.sliceThickness;

	std::vector<size_t> order_1D = {0};
	writeRealDataSet(ppotential, "dim1", &x_dim_data[0], x_size, 1, order_1D);
	writeRealDataSet(ppotential, "dim2", &y_dim_data[0], y_size, 1, order_1D);
	writeRealDataSet(ppotential, "dim3", &z_dim_data[0], z_size, 1, order_1D);

	//dimension attributes
	H5::DataSet dim1 = ppotential.openDataSet("dim1");
	H5::DataSet dim2 = ppotential.openDataSet("dim2");
	H5::DataSet dim3 = ppotential.openDataSet("dim3");

	writeScalarAttribute(dim1, "name", "R_x");
	writeScalarAttribute(dim2, "name", "R_y");
	writeScalarAttribute(dim3, "name", "R_z");

	writeScalarAttribute(dim1, "units", "[Å]");
	writeScalarAttribute(dim2, "units", "[Å]");
	writeScalarAttribute(dim3, "units", "[Å]");

	//create dataset
	//first, in potential array and re-stride
	std::vector<size_t> order_3D = {0, 1, 2};
	hsize_t dataDims[3] = {pars.imageSize[1], pars.imageSize[0], pars.numPlanes};
	writeRealDataSet(ppotential, "realslice", &pars.pot[0], dataDims, 3, order_3D);

	dim1.close();
	dim2.close();
	dim3.close();
	ppotential.close();
}

void saveHRTEM(Parameters<PRISMATIC_FLOAT_PRECISION> &pars, 
				Array3D<std::complex<PRISMATIC_FLOAT_PRECISION>> &net_output_c, 
				Array3D<PRISMATIC_FLOAT_PRECISION> &net_output)
{
	//open group and write tilt images incrementally in loop
	//this avoids restriding S-matrix array all at once, which is memory intensive
	H5::Group hrtem_group = pars.outputFile.openGroup("4DSTEM_simulation/data/realslices/HRTEM");
	hsize_t mdims[3] = {pars.Scompact.get_dimi(), pars.Scompact.get_dimj(), pars.numberBeams};
	std::vector<size_t> order = {0,1,2};
	
	if(pars.meta.saveComplexOutputWave)
	{
		writeComplexDataSet(hrtem_group, "realslice", &net_output_c[0], mdims, 3, order);
	}
	else
	{
		writeRealDataSet(hrtem_group, "realslice", &net_output[0], mdims, 3, order);
	}

	hrtem_group.close();
}

void saveSTEM(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	pars.outputFile = H5::H5File(pars.meta.filenameOutput.c_str(), H5F_ACC_RDWR);
	if (pars.meta.save3DOutput)
	{
		setupVDOutput(pars);
		hsize_t mdims[3] = {pars.numXprobes, pars.numYprobes, pars.Ndet};
		size_t strides = mdims[0]*mdims[1]*mdims[2];
		std::vector<size_t> order = {0,1,2};
		for (auto j = 0; j < pars.numLayers; j++)
		{
			std::string nameString;
			nameString = "4DSTEM_simulation/data/realslices/virtual_detector_depth" + getDigitString(j);
			nameString = nameString + pars.currentTag;
			std::cout << nameString << std::endl;
			H5::Group dataGroup = pars.outputFile.openGroup(nameString);
			if(pars.meta.saveComplexOutputWave)
			{
				writeComplexDataSet(dataGroup, "realslice", &pars.net_output_c[j*strides], mdims, 3, order);
			}
			else
			{
				Array3D<PRISMATIC_FLOAT_PRECISION> tmp_array = zeros_ND<3, PRISMATIC_FLOAT_PRECISION>({{pars.net_output.get_dimi(), pars.net_output.get_dimj(), pars.net_output.get_dimk()}});
				for(auto ii = 0; ii < pars.net_output.get_dimi(); ii++)
				{
					for(auto jj = 0; jj < pars.net_output.get_dimj(); jj++)
					{
						for(auto kk = 0; kk < pars.net_output.get_dimk(); kk++)
						{
							tmp_array.at(ii,jj,kk) = pars.net_output.at(j,kk,jj,ii);
						}
					}
				}
				writeRealDataSet_inOrder(dataGroup, "realslice", &tmp_array[0], mdims, 3);
			}
			dataGroup.close();
		}
	}

	if (pars.meta.save2DOutput)
	{
		size_t lower = std::max((size_t)0, (size_t)(pars.meta.integrationAngleMin / pars.meta.detectorAngleStep));
		size_t upper = std::min(pars.detectorAngles.size(), (size_t)(pars.meta.integrationAngleMax / pars.meta.detectorAngleStep));
		setup2DOutput(pars);
		std::vector<size_t> order = {0,1};
		hsize_t mdims[2] = {pars.numXprobes, pars.numYprobes};
		for (auto j = 0; j < pars.numLayers; j++)
		{
			std::string nameString = "4DSTEM_simulation/data/realslices/annular_detector_depth" + getDigitString(j);
			nameString += pars.currentTag;
			H5::Group dataGroup = pars.outputFile.openGroup(nameString.c_str());

			if(pars.meta.saveComplexOutputWave)
			{
				Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> prism_image;
				//need to initiliaze output image at each slice to prevent overflow of value
				prism_image = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>(
					{{pars.net_output_c.get_dimj(), pars.net_output_c.get_dimk()}});

				for (auto y = 0; y < pars.net_output_c.get_dimk(); ++y)
				{
					for (auto x = 0; x < pars.net_output_c.get_dimj(); ++x)
					{
						for (auto b = lower; b < upper; ++b)
						{
							prism_image.at(x, y) += pars.net_output_c.at(j, y, x, b);
						}
					}
				}
				writeComplexDataSet(dataGroup, "realslice", &prism_image[0], mdims, 2, order);
			}
			else
			{
				Array2D<PRISMATIC_FLOAT_PRECISION> prism_image;
				//need to initiliaze output image at each slice to prevent overflow of value
				prism_image = zeros_ND<2, PRISMATIC_FLOAT_PRECISION>(
					{{pars.net_output.get_dimj(), pars.net_output.get_dimk()}});

				for (auto y = 0; y < pars.net_output.get_dimk(); ++y)
				{
					for (auto x = 0; x < pars.net_output.get_dimj(); ++x)
					{
						for (auto b = lower; b < upper; ++b)
						{
							prism_image.at(x, y) += pars.net_output.at(j, y, x, b);
						}
					}
				}
				writeRealDataSet(dataGroup, "realslice", &prism_image[0], mdims, 2, order);
			}
			dataGroup.close();
		}
	}

	if (pars.meta.saveDPC_CoM and not pars.meta.saveComplexOutputWave)
	{
		setupDPCOutput(pars);
		hsize_t mdims[3] = {pars.numXprobes, pars.numYprobes, 2};
		size_t strides = mdims[0]*mdims[1]*mdims[2];
		std::vector<size_t> order = {2,0,1};
		for (auto j = 0; j < pars.numLayers; j++)
		{
			std::string nameString = "4DSTEM_simulation/data/realslices/DPC_CoM_depth" + getDigitString(j);
			nameString += pars.currentTag;
			H5::Group dataGroup = pars.outputFile.openGroup(nameString.c_str());
			writeRealDataSet(dataGroup, "realslice", &pars.net_DPC_CoM[j*strides], mdims, 3, order);
			dataGroup.close();
		}
	}

	pars.outputFile.close();
};

void configureImportFP(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{

	if(pars.meta.importSMatrix)
	{
		std::cout << "Skipping PRISM01. Using precalculated scattering matrix from: "  << pars.meta.importFile << std::endl;
	}
	else if(pars.meta.importPotential)
	{
		std::cout << "Using precalculated potential from " << pars.meta.importFile << std::endl;
		std::string groupName = "4DSTEM_simulation/data/realslices/";
		std::string baseName = "ppotential_fp";
		H5::H5File inFile = H5::H5File(pars.meta.importFile.c_str(), H5F_ACC_RDONLY);
		H5::Group realslices = inFile.openGroup(groupName.c_str());
		int configurations = countDataGroups(realslices, baseName);

		std::cout << configurations << " frozen phonon configurations available in " << pars.meta.importFile << std::endl;
		int tmp_fp = pars.meta.numFP;

		//if user requests more than the available configurations, only run number available
		if(pars.meta.numFP > 1)
		{
			pars.meta.numFP = (tmp_fp > configurations) ? configurations : tmp_fp;
			std::cout << "User requested " << tmp_fp  << " frozen phonons." << std::endl;
			std::cout << "Running " << pars.meta.numFP << " frozen phonons out of " << configurations << " available configurations." << std::endl;
		}
		else
		{
			if( not pars.meta.userSpecifiedNumFP)
			{
				//if user specifically specifies to run a single frozen phonon, this is skipped and only the first configuration will run
				pars.meta.numFP = configurations;
			}
		}
	}

	if(pars.meta.importSMatrix)
	{
		std::string groupName = "4DSTEM_simulation/data/realslices/";
		std::string baseName = "smatrix_fp";
		H5::H5File inFile = H5::H5File(pars.meta.importFile.c_str(), H5F_ACC_RDONLY);
		H5::Group realslices = inFile.openGroup(groupName.c_str());
		int configurations = countDataGroups(realslices, baseName);

		std::cout << configurations << " frozen phonon configurations available in " << pars.meta.importFile << std::endl;
		int tmp_fp = pars.meta.numFP;

		//if user requests more than the available configurations, only run number available
		if(pars.meta.numFP > 1)
		{
			pars.meta.numFP = (tmp_fp > configurations) ? configurations : tmp_fp;
			std::cout << "User requested " << tmp_fp  << " frozen phonons." << std::endl;
			std::cout << "Running " << pars.meta.numFP << " frozen phonons out of " << configurations << " available configurations." << std::endl;
		}
		else
		{
			if( not pars.meta.userSpecifiedNumFP)
			{
				//if user specifically specifies to run a single frozen phonon, this is skipped and only the first configuration will run
				pars.meta.numFP = configurations;
			}
		}
	}

};

std::string getDigitString(int digit)
{
	char buffer[20];
	sprintf(buffer, "%04d", digit);
	std::string output = buffer;
	return output;
};

void writeMetadata(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	//set up group
	H5::Group metadata = pars.outputFile.openGroup("4DSTEM_simulation/metadata/metadata_0/original");
	H5::Group sim_params = metadata.createGroup("simulation_parameters");

	//write all parameters as attributes

	//create common dataspaces
	H5::DataSpace str_name_ds(H5S_SCALAR); //string dataspaces and types
	H5::StrType strdatatype(H5::PredType::C_S1, 256);
	H5::DataSpace scalar_attr(H5S_SCALAR);

	//create string attributes
	H5std_string algorithm = (pars.meta.algorithm == Algorithm::Multislice) ? "m" : "p";

	writeScalarAttribute(sim_params, "i", pars.meta.filenameAtoms);
	writeScalarAttribute(sim_params, "a", algorithm);

	//create scalar logical/integer attributes
	writeScalarAttribute(sim_params, "fx", (int) pars.meta.interpolationFactorX);
	writeScalarAttribute(sim_params, "fy", (int) pars.meta.interpolationFactorY);
	writeScalarAttribute(sim_params, "F", (int) pars.meta.numFP);
	writeScalarAttribute(sim_params, "ns", (int) pars.meta.numSlices);
	writeScalarAttribute(sim_params, "te", (int) pars.meta.includeThermalEffects);
	writeScalarAttribute(sim_params, "oc", (int) pars.meta.includeOccupancy);
	writeScalarAttribute(sim_params, "3D", (int) pars.meta.save3DOutput);
	writeScalarAttribute(sim_params, "4D", (int) pars.meta.save4DOutput);
	writeScalarAttribute(sim_params, "DPC", (int) pars.meta.saveDPC_CoM);
	writeScalarAttribute(sim_params, "ps", (int) pars.meta.savePotentialSlices);
	writeScalarAttribute(sim_params, "nqs", (int) pars.meta.nyquistSampling);

	//create scalar float attributes
	writeScalarAttribute(sim_params, "px", pars.meta.realspacePixelSize[1]);
	writeScalarAttribute(sim_params, "py", pars.meta.realspacePixelSize[0]);
	writeScalarAttribute(sim_params, "P", pars.meta.potBound);
	writeScalarAttribute(sim_params, "s", pars.meta.sliceThickness);
	writeScalarAttribute(sim_params, "zs", pars.meta.zStart);
	writeScalarAttribute(sim_params, "E", pars.meta.E0 / 1000);
	writeScalarAttribute(sim_params, "A", pars.meta.alphaBeamMax * 1000);
	writeScalarAttribute(sim_params, "rx", pars.meta.probeStepX);
	writeScalarAttribute(sim_params, "ry", pars.meta.probeStepY);
	writeScalarAttribute(sim_params, "df", pars.meta.probeDefocus);
	writeScalarAttribute(sim_params, "C3", pars.meta.C3);
	writeScalarAttribute(sim_params, "C5", pars.meta.C5);
	writeScalarAttribute(sim_params, "sa", pars.meta.probeSemiangle * 1000);
	writeScalarAttribute(sim_params, "d", pars.meta.detectorAngleStep * 1000);
	writeScalarAttribute(sim_params, "tx", pars.meta.probeXtilt * 1000);
	writeScalarAttribute(sim_params, "ty", pars.meta.probeYtilt * 1000);

	//create vector spaces
	hsize_t two[1] = {2};
	hsize_t three[1] = {3};
	H5::DataSpace v_two_dataspace(1, two);
	H5::DataSpace v_three_dataspace(1, three);

	H5::Attribute cell_dim_attr = sim_params.createAttribute("c", PFP_TYPE, v_three_dataspace);
	H5::Attribute tile_attr = sim_params.createAttribute("t", PFP_TYPE, v_three_dataspace);
	H5::Attribute scanWindow_x_attr = sim_params.createAttribute("wx", PFP_TYPE, v_two_dataspace);
	H5::Attribute scanWindow_y_attr = sim_params.createAttribute("wy", PFP_TYPE, v_two_dataspace);

	H5::Attribute scanWindow_x_r_attr;
	H5::Attribute scanWindow_y_r_attr;
	if (pars.meta.realSpaceWindow_x) scanWindow_x_r_attr = sim_params.createAttribute("wxr", PFP_TYPE, v_two_dataspace);
	if (pars.meta.realSpaceWindow_y) scanWindow_y_r_attr = sim_params.createAttribute("wyr", PFP_TYPE, v_two_dataspace);

	H5::Attribute save2D_attr;
	if (pars.meta.save2DOutput) save2D_attr = sim_params.createAttribute("2D", PFP_TYPE, v_two_dataspace);
	
	PRISMATIC_FLOAT_PRECISION tmp_buffer[2];

	if (pars.meta.realSpaceWindow_x)
	{
		tmp_buffer[0] = pars.meta.scanWindowXMin_r;
		tmp_buffer[1] = pars.meta.scanWindowXMax_r;
		scanWindow_x_r_attr.write(PFP_TYPE, tmp_buffer);
	}

	if (pars.meta.realSpaceWindow_y)
	{
		tmp_buffer[0] = pars.meta.scanWindowYMin_r;
		tmp_buffer[1] = pars.meta.scanWindowYMax_r;
		scanWindow_y_r_attr.write(PFP_TYPE, tmp_buffer);
	}

	if (pars.meta.save2DOutput)
	{
		tmp_buffer[0] = pars.meta.integrationAngleMin * 1000;
		tmp_buffer[1] = pars.meta.integrationAngleMax * 1000;
		save2D_attr.write(PFP_TYPE, tmp_buffer);
	}

	tmp_buffer[0] = pars.meta.scanWindowXMin;
	tmp_buffer[1] = pars.meta.scanWindowXMax;
	scanWindow_x_attr.write(PFP_TYPE, tmp_buffer);

	tmp_buffer[0] = pars.meta.scanWindowYMin;
	tmp_buffer[1] = pars.meta.scanWindowYMax;
	scanWindow_y_attr.write(PFP_TYPE, tmp_buffer);

	int tile_buffer[3] = {pars.meta.tileX, pars.meta.tileY, pars.meta.tileZ};
	tile_attr.write(H5::PredType::NATIVE_INT, tile_buffer);

	PRISMATIC_FLOAT_PRECISION cellBuffer[3] = {pars.meta.cellDim[0], pars.meta.cellDim[1], pars.meta.cellDim[2]};
	cell_dim_attr.write(PFP_TYPE, cellBuffer);

	metadata.close();
};

Array2D<PRISMATIC_FLOAT_PRECISION> readDataSet2D(const std::string &filename, const std::string &dataPath)
{
	H5::H5File input = H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	H5::DataSet dataset = input.openDataSet(dataPath.c_str());
	H5::DataSpace dataspace = dataset.getSpace();

	hsize_t dims_out[2];
	int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
	H5::DataSpace mspace(2,dims_out);

	std::vector<PRISMATIC_FLOAT_PRECISION> data_in(dims_out[0]*dims_out[1]);
	dataset.read(&data_in[0], PFP_TYPE, mspace, dataspace);

	mspace.close();
	dataspace.close();
	dataset.close();
	input.close();

	Array2D<PRISMATIC_FLOAT_PRECISION> data = ArrayND<2, std::vector<PRISMATIC_FLOAT_PRECISION>>(data_in, {{dims_out[1], dims_out[0]}});

	return data;
};

Array3D<PRISMATIC_FLOAT_PRECISION> readDataSet3D(const std::string &filename, const std::string &dataPath)
{
	H5::H5File input = H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	H5::DataSet dataset = input.openDataSet(dataPath.c_str());
	H5::DataSpace dataspace = dataset.getSpace();

	hsize_t dims_out[3];
	int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
	H5::DataSpace mspace(3,dims_out);

	std::vector<PRISMATIC_FLOAT_PRECISION> data_in(dims_out[0]*dims_out[1]*dims_out[2]);
	dataset.read(&data_in[0], PFP_TYPE, mspace, dataspace);

	mspace.close();
	dataspace.close();
	dataset.close();
	input.close();

	Array3D<PRISMATIC_FLOAT_PRECISION> data = ArrayND<3, std::vector<PRISMATIC_FLOAT_PRECISION>>(data_in, {{dims_out[2], dims_out[1], dims_out[0]}});

	return data;
};

Array4D<PRISMATIC_FLOAT_PRECISION> readDataSet4D(const std::string &filename, const std::string &dataPath)
{
	H5::H5File input = H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	H5::DataSet dataset = input.openDataSet(dataPath.c_str());
	H5::DataSpace dataspace = dataset.getSpace();

	hsize_t dims_out[4];
	int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
	H5::DataSpace mspace(4,dims_out);

	std::vector<PRISMATIC_FLOAT_PRECISION> data_in(dims_out[0]*dims_out[1]*dims_out[2]*dims_out[3]);
	dataset.read(&data_in[0], PFP_TYPE, mspace, dataspace);

	mspace.close();
	dataspace.close();
	dataset.close();
	input.close();


	//mem dims are stored 0->3 kx, ky, qx, qy
	//flipping x, y back for storage into 4D array
	Array4D<PRISMATIC_FLOAT_PRECISION> data = ArrayND<4, std::vector<PRISMATIC_FLOAT_PRECISION>>(data_in, {{dims_out[1], dims_out[0], dims_out[3], dims_out[2]}});
	std::array<size_t, 4> dims_in = {dims_out[1], dims_out[0], dims_out[3], dims_out[2]};
	std::array<size_t, 4> order = {1, 0, 3, 2};
	data = restride(data, dims_in, order);

	return data;
};

Array4D<PRISMATIC_FLOAT_PRECISION> readDataSet4D_keepOrder(const std::string &filename, const std::string &dataPath)
{
	//copy of above, but does not move any dimensions around
	H5::H5File input = H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	H5::DataSet dataset = input.openDataSet(dataPath.c_str());
	H5::DataSpace dataspace = dataset.getSpace();

	hsize_t dims_out[4];
	hsize_t dims_out_switch[4];
	int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
	dims_out_switch[0] = dims_out[1];
	dims_out_switch[1] = dims_out[2];
	dims_out_switch[2] = dims_out[3];
	dims_out_switch[3] = dims_out[0];
	H5::DataSpace mspace(4,dims_out_switch);

	std::vector<PRISMATIC_FLOAT_PRECISION> data_in(dims_out[0]*dims_out[1]*dims_out[2]*dims_out[3]);
	dataset.read(&data_in[0], PFP_TYPE, mspace, dataspace);

	mspace.close();
	dataspace.close();
	dataset.close();
	input.close();

	Array4D<PRISMATIC_FLOAT_PRECISION> data = ArrayND<4, std::vector<PRISMATIC_FLOAT_PRECISION>>(data_in, {{dims_out[3], dims_out[2], dims_out[1], dims_out[0]}});
	return data;
};

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, PRISMATIC_FLOAT_PRECISION &val)
{
	//read an attribute from a group into val
	//overloaded by expected value type

	H5::H5File input = H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	H5::Group group = input.openGroup(groupPath);
	H5::Attribute attribute = group.openAttribute(attr);
	H5::DataType type =  attribute.getDataType();
	attribute.read(type,&val);

	attribute.close();
	group.close();
	input.close();
	
}

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, PRISMATIC_FLOAT_PRECISION *val)
{
	//read an attribute from a group into val
	//overloaded by expected value type

	H5::H5File input = H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	H5::Group group = input.openGroup(groupPath);
	H5::Attribute attribute = group.openAttribute(attr);
	H5::DataType type =  attribute.getDataType();
	attribute.read(type,&val[0]);

	attribute.close();
	group.close();
	input.close();
	
}

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, int &val)
{
	//read an attribute from a group into val
	//overloaded by expected value type

	H5::H5File input = H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	H5::Group group = input.openGroup(groupPath);
	H5::Attribute attribute = group.openAttribute(attr);
	H5::DataType type =  attribute.getDataType();
	attribute.read(type,&val);

	attribute.close();
	group.close();
	input.close();

}

void readAttribute(const std::string &filename, const std::string &groupPath, const std::string &attr, std::string &val)
{
	//read an attribute from a group into val
	//overloaded by expected value type

	H5::H5File input = H5::H5File(filename.c_str(), H5F_ACC_RDONLY);
	H5::Group group = input.openGroup(groupPath);
	H5::Attribute attribute = group.openAttribute(attr);
	H5::DataType type =  attribute.getDataType();

	H5std_string tmpVal;
	attribute.read(type,tmpVal);
	val = tmpVal;

	attribute.close();
	group.close();
	input.close();

}

void writeComplexDataSet(H5::Group group, const std::string &dsetname,
							const std::complex<PRISMATIC_FLOAT_PRECISION> *buffer,
							const hsize_t *mdims,
							const size_t &rank,
							std::vector<size_t> &order)
{
	H5::CompType complex_type = H5::CompType(sizeof(complex_float_t));
	const H5std_string re_str("r"); //using h5py default configuration
	const H5std_string im_str("i");
	complex_type.insertMember(re_str, 0, PFP_TYPE);
	complex_type.insertMember(im_str, 4, PFP_TYPE);

	//create dataset and write
	H5::DataSpace mspace(rank, mdims);

	H5::DataSet complex_dset;
	if(group.nameExists(dsetname.c_str()))
	{
		complex_dset = group.openDataSet(dsetname.c_str());
	}
	else
	{
		complex_dset = group.createDataSet(dsetname.c_str(), complex_type, mspace);
	}		
	
	H5::DataSpace fspace = complex_dset.getSpace();
	if(rank > 1)
	{
		std::vector<size_t> rdims;
		for(auto i =0; i < rank; i++) rdims.push_back(mdims[i]);
		restrideElements(fspace, rdims, order);
	}
	complex_dset.write(buffer, complex_type, mspace, fspace);

	//close spaces
	fspace.close();
	mspace.close();
	complex_dset.close();

}

void writeRealDataSet(H5::Group group,
						const std::string &dsetname,
						const PRISMATIC_FLOAT_PRECISION *buffer,
						const hsize_t *mdims,
						const size_t &rank,
						std::vector<size_t> &order)
{
	//create dataset and write
	H5::DataSpace mspace(rank, mdims);

	H5::DataSet real_dset;
	if(group.nameExists(dsetname.c_str()))
	{
		real_dset = group.openDataSet(dsetname.c_str());
	}
	else
	{
		real_dset = group.createDataSet(dsetname.c_str(), PFP_TYPE, mspace);
	}		
	
	H5::DataSpace fspace = real_dset.getSpace();
	std::cout << "Organizing dataspace." << std::endl;
	if(rank > 1)
	{
		std::vector<size_t> rdims;
		for(auto i =0; i < rank; i++) rdims.push_back(mdims[i]);
		restrideElements(fspace, rdims, order);
	}
	std::cout << "Writing dataset." << std::endl;
	real_dset.write(buffer, PFP_TYPE, mspace, fspace);
	std::cout << "Write finished." << std::endl;

	//close spaces
	fspace.close();
	mspace.close();
	real_dset.close();

};

void writeRealDataSet_inOrder(H5::Group group,
						const std::string &dsetname,
						const PRISMATIC_FLOAT_PRECISION *buffer,
						const hsize_t *mdims,
						const size_t &rank)
{
	//create dataset and write
	H5::DataSpace mspace(rank, mdims);

	H5::DataSet real_dset;
	if(group.nameExists(dsetname.c_str()))
	{
		real_dset = group.openDataSet(dsetname.c_str());
	}
	else
	{
		real_dset = group.createDataSet(dsetname.c_str(), PFP_TYPE, mspace);
	}		
	
	H5::DataSpace fspace = real_dset.getSpace();
	real_dset.write(buffer, PFP_TYPE, mspace, fspace);

	//close spaces
	fspace.close();
	mspace.close();
	real_dset.close();
};

void writeScalarAttribute(H5::H5Object &object, const std::string &name, const int &data)
{
	H5::DataSpace attr_dataspace(H5S_SCALAR);
	H5::Attribute attr = object.createAttribute(name.c_str(), H5::PredType::NATIVE_INT, attr_dataspace);
	attr.write(H5::PredType::NATIVE_INT, &data);
};

void writeScalarAttribute(H5::H5Object &object, const std::string &name, const PRISMATIC_FLOAT_PRECISION &data)
{
	H5::DataSpace attr_dataspace(H5S_SCALAR);
	H5::Attribute attr = object.createAttribute(name.c_str(), PFP_TYPE, attr_dataspace);
	attr.write(PFP_TYPE, &data);
};

void writeScalarAttribute(H5::H5Object &object, const std::string &name, const std::string &data)
{
	H5::StrType strdatatype(H5::PredType::C_S1, 256);
	const H5std_string h5_str(data);
	H5::DataSpace attr_dataspace(H5S_SCALAR);
	H5::Attribute attr = object.createAttribute(name.c_str(), strdatatype, attr_dataspace);
	attr.write(strdatatype, h5_str);
};

int countDataGroups(H5::Group group, const std::string &basename)
{
	//count number of subgroups in specified group with desired basename + "####"
	int count = 0;
	std::string currentName = basename + getDigitString(count);
	while(group.nameExists(currentName.c_str()))
	{
		count += 1;
		currentName = basename + getDigitString(count);
	}

	return count;

};

int countDimensions(H5::Group group, const std::string &basename)
{
	//count number of dimensions in specified group with desired basename + "#"
	//basename used to specified dim or sgdim
	int count = 1;
	std::string currentName = basename + std::to_string(count);
	while(group.nameExists(currentName.c_str()))
	{
		count += 1;
		currentName = basename + std::to_string(count);
	}

	return count-1;

};

void configureSupergroup(H5::Group &new_sg,
						H5::Group &sourceExample,
						const std::vector<std::vector<PRISMATIC_FLOAT_PRECISION>> &sgdims,
						const std::vector<std::string> &sgdims_name,
						const std::vector<std::string> &sgdims_units)
{
	std::cout << "writing supergroup attributes" << std::endl;
	//write group type attribute
	H5::DataSpace attr1_dataspace(H5S_SCALAR);
	H5::Attribute emd_group_type = new_sg.createAttribute("emd_group_type", H5::PredType::NATIVE_INT, attr1_dataspace);
	int group_type = 3;
	emd_group_type.write(H5::PredType::NATIVE_INT, &group_type);

	//write metadata attribute
	H5::DataSpace attr2_dataspace(H5S_SCALAR);
	H5::Attribute metadata_group = new_sg.createAttribute("metadata", H5::PredType::NATIVE_INT, attr2_dataspace);
	int mgroup = 0;
	metadata_group.write(H5::PredType::NATIVE_INT, &mgroup);

	H5::DataSpace str_name_ds(H5S_SCALAR);
	H5::StrType strdatatype(H5::PredType::C_S1, 256);

	//write common dimensions
	int numDims = countDimensions(sourceExample, "dim");
	std::cout << "writing " << numDims << " common dimensions" << std::endl;
	for(auto i = 0; i < numDims; i++)
	{
		H5::DataSet dim_data = sourceExample.openDataSet("dim"+std::to_string(i+1));
		copyDataSet(new_sg, dim_data);
	}

	//write supergroup dimensions
	std::cout << "writing " << sgdims.size() << " supergroup dimensions" << std::endl;
	for(auto i = 0; i < sgdims.size(); i++)
	{
		hsize_t dim_size[1] = {sgdims[i].size()};
		H5::DataSpace dim_mspace(1, dim_size);
		H5::DataSet dim = new_sg.createDataSet("sgdim"+std::to_string(i+1), PFP_TYPE, dim_mspace);
		H5::DataSpace dim_fspace = dim.getSpace();
		dim.write(&sgdims[i][0], PFP_TYPE, dim_mspace, dim_fspace);

		//dimension attributes
		const H5std_string dim_name_str(sgdims_name[i]);
		const H5std_string dim_unit_str(sgdims_units[i]);
		H5::Attribute dim_name = dim.createAttribute("name", strdatatype, str_name_ds);
		H5::Attribute dim_unit = dim.createAttribute("units", strdatatype, str_name_ds);
		dim_name.write(strdatatype, dim_name_str);
		dim_unit.write(strdatatype, dim_unit_str);
	}

};

void writeVirtualDataSet(H5::Group group,
						const std::string &dsetName,
						std::vector<H5::DataSet> &datasets,
						std::vector<std::vector<size_t>> indices)
{
	//determine new dimensionality of VDS and get extent along each dim
	size_t new_rank = indices[0].size();
	size_t max_dims[new_rank] = {0};
	for(auto i = 0; i < indices.size(); i++)
	{
		for(auto j = 0; j < new_rank; j++)
		{
			max_dims[j] = (max_dims[j] < indices[i][j]) ? indices[i][j] : max_dims[j];
		}
	}

	//organize dimensions and create mem space for virtual dataset
	std::cout << "setting up virtual dataset memspace" << std::endl;
	H5::DataSpace sampleSpace = datasets[0].getSpace();
	int rank = sampleSpace.getSimpleExtentNdims();
	hsize_t dims_out[rank];
	int ndims = sampleSpace.getSimpleExtentDims(dims_out, NULL); //nidms and rank are redundant, but rank is not known a priori
	
	hsize_t mdims[rank+new_rank] = {0};
	hsize_t mdims_ind[rank+new_rank] = {0};

	for(auto i = 0; i < rank; i++)
	{
		mdims[i] = dims_out[i];
		mdims_ind[i] = dims_out[i];
	}

	for(auto i = rank; i < rank+new_rank; i++)
	{
		mdims[i] = max_dims[i-rank]+1;
		mdims_ind[i] = 1;
	}

	sampleSpace.close();

	//create virtual dataspace and plist mapping
	H5::DataSpace vds_mspace(rank+new_rank, mdims);
	H5::DataSpace src_mspace;
	H5::DSetCreatPropList plist;
    std::string path;
	hsize_t offset[rank+new_rank] = {0};
	
	std::cout << "mapping virtual dataset dataspace to source dataset spaces" << std::endl;
	for(auto i = 0; i < datasets.size(); i++)
	{
		path = datasets[i].getObjName();
		src_mspace = datasets[i].getSpace();

		for(auto j = rank; j < rank+new_rank; j++) offset[j] = indices[i][j-rank];
		
		vds_mspace.selectHyperslab(H5S_SELECT_SET, mdims_ind, offset);
		plist.setVirtual(vds_mspace, datasets[i].getFileName(), path, src_mspace);
	}

	for(auto i = rank; i < rank+new_rank; i++) offset[i] = 0;
	vds_mspace.selectHyperslab(H5S_SELECT_SET, mdims, offset);
	std::cout << "creating virtual dataset" << std::endl;
	H5::DataSet vds = group.createDataSet(dsetName, datasets[0].getDataType(), vds_mspace, plist);

	src_mspace.close();
	vds.close();

};

void depthSeriesSG(H5::H5File &file)
{
	H5::Group supergroups = file.openGroup("4DSTEM_simulation/data/supergroups");
	H5::Group depthSeries = supergroups.createGroup("vd_depth_series");

	//gather datasets and create mapping
	std::string basename = "virtual_detector_depth";
	H5::Group realslices = file.openGroup("4DSTEM_simulation/data/realslices");

	int numDataSets = countDataGroups(realslices, basename);
	std::vector<H5::DataSet> datasets;
	std::vector<std::vector<size_t>> indices;

	std::vector<PRISMATIC_FLOAT_PRECISION> depths;
	PRISMATIC_FLOAT_PRECISION tmp_depth;
	for(auto i = 0; i < numDataSets; i++)
	{
		std::string tmp_name = basename+getDigitString(i);
		H5::Group tmp_group = realslices.openGroup(tmp_name.c_str());
		H5::DataSet tmp_dataset = tmp_group.openDataSet("realslice");
		datasets.push_back(tmp_dataset);
		indices.push_back(std::vector<size_t>{i});
		readAttribute(tmp_group.getFileName(), tmp_group.getObjName(), "output_depth", tmp_depth);
		depths.push_back(tmp_depth);
	}

	//write dataset
	writeVirtualDataSet(depthSeries, "supergroup", datasets, indices);

	//configure supergroup
	//gather dim properties from first datagroup and write sgdims
	H5::Group firstGroup = realslices.openGroup(basename+getDigitString(0));
	std::vector<std::vector<PRISMATIC_FLOAT_PRECISION>> sgdims;
	sgdims.push_back(depths);
	
	std::vector<std::string> sgdims_name;
	sgdims_name.push_back("Depth");

	std::vector<std::string> sgdims_units;
	sgdims_units.push_back("[Å]");

	configureSupergroup(depthSeries, firstGroup, sgdims, sgdims_name, sgdims_units);

}

void CCseriesSG(H5::H5File &file)
{
	H5::Group supergroups = file.openGroup("4DSTEM_simulation/data/supergroups");
	H5::Group CC_series = supergroups.createGroup("vd_CC_series");

	//gather datasets and create mapping
	std::string basename = "virtual_detector_depth0000_df";
	H5::Group realslices = file.openGroup("4DSTEM_simulation/data/realslices");

	int numDataSets = countDataGroups(realslices, basename);
	std::vector<H5::DataSet> datasets;
	std::vector<std::vector<size_t>> indices;

	std::vector<PRISMATIC_FLOAT_PRECISION> defocii;
	PRISMATIC_FLOAT_PRECISION tmp_defocus;
	std::cout << "reading attributes from datasets" << std::endl;
	for(auto i = 0; i < numDataSets; i++)
	{
		std::string tmp_name = basename+getDigitString(i);
		H5::Group tmp_group = realslices.openGroup(tmp_name.c_str());
		H5::DataSet tmp_dataset = tmp_group.openDataSet("realslice");
		datasets.push_back(tmp_dataset);
		indices.push_back(std::vector<size_t>{i});
		readAttribute(tmp_group.getFileName(), tmp_group.getObjName(), "output_defocus", tmp_defocus);
		defocii.push_back(tmp_defocus);
	}

	//write dataset
	std::cout << "configuring virtual dataset" << std::endl;
	writeVirtualDataSet(CC_series, "supergroup", datasets, indices);

	//configure supergroup
	//gather dim properties from first datagroup and write sgdims
	H5::Group firstGroup = realslices.openGroup(basename+getDigitString(0));
	std::vector<std::vector<PRISMATIC_FLOAT_PRECISION>> sgdims;
	sgdims.push_back(defocii);
	
	std::vector<std::string> sgdims_name;
	sgdims_name.push_back("Defocus");

	std::vector<std::string> sgdims_units;
	sgdims_units.push_back("[Å]");

	std::cout << "configuring supergroup" << std::endl;
	configureSupergroup(CC_series, firstGroup, sgdims, sgdims_name, sgdims_units);
}

std::string reducedDataSetName(std::string &fullPath)
{
	size_t index = fullPath.find_last_of("/");
	return fullPath.substr(index+1);
}

void copyDataSet(H5::Group &targetGroup, H5::DataSet &source)
{
	//grab properties from source dataset
	std::string dsName = source.getObjName();
	dsName = reducedDataSetName(dsName);

	H5::DataSpace sourceSpace = source.getSpace();
	int rank = sourceSpace.getSimpleExtentNdims();
	hsize_t dims_out[rank];
	int ndims = sourceSpace.getSimpleExtentDims(dims_out, NULL); //nidms and rank are redundant, but rank is not known a priori

	//create buffer array and read data from source
	H5::DataSpace mspace(rank, dims_out);
	size_t bufferSize = source.getInMemDataSize(); //size in bytes
	unsigned char buffer[bufferSize]; //unsigned char is always byte sized, let's us be agnostic to storage type of array
	source.read(buffer, source.getDataType(), mspace, sourceSpace);	

	//create new 
	H5::DataSet target = targetGroup.createDataSet(dsName.c_str(), source.getDataType(), mspace);
	H5::DataSpace fspace = target.getSpace();
	target.write(&buffer[0], source.getDataType(), mspace, fspace);

	//look for attributes in source and copy
	for(auto i = 0; i < source.getNumAttrs(); i++)
	{
		H5::Attribute tmp_attr = source.openAttribute(i);
		H5::DataSpace attr_space = tmp_attr.getSpace();
		int attr_rank = attr_space.getSimpleExtentNdims();
		hsize_t attr_dims_out[attr_rank];
		int ndims = attr_space.getSimpleExtentDims(attr_dims_out, NULL); //nidms and rank are redundant, but rank is not known a priori
		
		H5::DataSpace t_attr_space(attr_rank, attr_dims_out);
		H5::Attribute t_attr = target.createAttribute(tmp_attr.getName(), tmp_attr.getDataType(), t_attr_space);
		
		unsigned char attr_buffer[tmp_attr.getInMemDataSize()];
		tmp_attr.read(tmp_attr.getDataType(), attr_buffer);
		t_attr.write(tmp_attr.getDataType(), &attr_buffer[0]);
	}

	target.close();


};

void restrideElements(H5::DataSpace &fspace, std::vector<size_t> &dims, std::vector<size_t> &order)
{
	//clear selection
	fspace.selectNone();

	//get rank
	size_t rank = dims.size();
	
	//block indicates which dimensions are kept contiguous in selects
	size_t block0 = dims[order[0]];
	size_t block1 = dims[order[1]];
	size_t Nblocks = 1;
	std::vector<size_t> Nblock_order;
	for(auto i = 2; i < rank; i++)	Nblocks *= dims[order[i]];

	//setup divisors and modulos
	//not generalized since will only deal with at most 4 dimension restrides
	std::vector<size_t> divs(rank);
	std::vector<size_t> mods(rank);
	divs[0] = 1; 
	mods[0] = block0;
	
	divs[1] = block0;
	mods[1] = block0*block1;

	if(rank == 3)
	{
		divs[2] = 1;
		mods[2] = Nblocks;
	}
	else if(rank == 4)
	{
		divs[2] = 1;
		mods[2] = dims[order[2]];
		
		divs[3] = dims[order[2]];
		mods[3] = Nblocks;
	}

	//fill coordinates and make element selections
	hsize_t * coords = (hsize_t*) malloc(block0*block1*rank*sizeof(hsize_t));
	for(auto n = 0; n < Nblocks; n++)
	{
		for(auto m = 0; m < block0*block1; m++)
		{
			for(auto j = 0; j < 2; j++)	coords[m*rank+order[j]] = (m / divs[j]) % mods[j];
			for(auto j = 2; j < rank; j++)	coords[m*rank+order[j]] = (n / divs[j]) % mods[j];
		}
		fspace.selectElements(H5S_SELECT_APPEND, block0*block1, coords);
	}
	
};


void restrideElements_subset(H5::DataSpace &fspace, std::vector<size_t> &dims, std::vector<size_t> &order, std::vector<size_t> &offset)
{
	//for large arrays; restrides a subset and adds offset to non-continguous dims
	//clear selection
	fspace.selectNone();

	//get rank
	size_t rank = dims.size();
	
	//block indicates which dimensions are kept contiguous in selects
	size_t block0 = dims[order[0]];
	size_t block1 = dims[order[1]];
	size_t Nblocks = 1;
	std::vector<size_t> Nblock_order;
	for(auto i = 2; i < rank; i++)	Nblocks *= dims[order[i]];

	//setup divisors and modulos
	//not generalized since will only deal with at most 4 dimension restrides
	std::vector<size_t> divs(rank);
	std::vector<size_t> mods(rank);
	divs[0] = 1; 
	mods[0] = block0;
	
	divs[1] = block0;
	mods[1] = block0*block1;

	//fill coordinates and make element selections
	hsize_t coords[block0*block1][rank];
	for(auto m = 0; m < block0*block1; m++)
	{
		for(auto j = 0; j < 2; j++)	coords[m][order[j]] = (m / divs[j]) % mods[j];
		for(auto j = 2; j < rank; j++)	coords[m][order[j]] =  offset[order[j]];
	}
	fspace.selectElements(H5S_SELECT_APPEND, block0*block1, &coords[0][0]);
	
};

hsize_t* restrideElements_subset(std::vector<size_t> &dims, std::vector<size_t> &order, std::vector<size_t> &offset)
{
	//get rank
	size_t rank = dims.size();
	
	//block indicates which dimensions are kept contiguous in selects
	size_t block0 = dims[order[0]];
	size_t block1 = dims[order[1]];

	//setup divisors and modulos
	//not generalized since will only deal with at most 4 dimension restrides
	std::vector<size_t> divs(rank);
	std::vector<size_t> mods(rank);
	divs[0] = 1; 
	mods[0] = block0;
	
	divs[1] = block0;
	mods[1] = block0*block1;

	//fill coordinates and make element selections
	hsize_t coords[block0*block1][rank];
	for(auto m = 0; m < block0*block1; m++)
	{
		for(auto j = 0; j < 2; j++)	coords[m][order[j]] = (m / divs[j]) % mods[j];
		for(auto j = 2; j < rank; j++)	coords[m][order[j]] =  offset[order[j]];
	}

	return (hsize_t *) coords;	
};

void createScratchFile(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	//TODO: decide home directory for windows
	pars.scratchFile = H5::H5File("prismatic_scratch.h5", H5F_ACC_TRUNC);
	H5::Group scratchGroup = pars.scratchFile.createGroup("scratch");

	//initialize datasets
	hsize_t mdims[4] = {pars.output.get_diml(), pars.output.get_dimk(), pars.output.get_dimj(), pars.output.get_dimi()};
	H5::DataSpace mspace(4,mdims);
	for(auto i = 0; i < pars.meta.seriesTags.size(); i++)
	{
		scratchGroup.createDataSet(pars.meta.seriesTags[i].c_str(), PFP_TYPE, mspace);
			
	}

	if(pars.meta.saveDPC_CoM)
	{
		hsize_t mdims_dpc[4] = {pars.DPC_CoM.get_diml(), pars.DPC_CoM.get_dimk(), pars.DPC_CoM.get_dimj(), pars.DPC_CoM.get_dimi()};
		std::cout << "MDIMS: " << mdims_dpc[0] << " " << mdims_dpc[1] << " " << mdims_dpc[2] << " " << mdims_dpc[3] << std::endl;
		H5::DataSpace mspace_dpc(4,mdims_dpc);
		for(auto i = 0; i < pars.meta.seriesTags.size(); i++)
		{
			std::string current_name = pars.meta.seriesTags[i]+"_DPC";
			scratchGroup.createDataSet(current_name.c_str(), PFP_TYPE, mspace_dpc);
				
		}
	}
	scratchGroup.close();
};

void removeScratchFile(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	//TODO: make better decision about where scratch file is written
	// and save name as a member of pars
	pars.scratchFile.close();
	if( remove( "prismatic_scratch.h5" ) != 0 )
        perror( "Error deleting scratch file" );
    else
        puts( "Scratch file successfully deleted" );
};

void updateScratchData(Parameters<PRISMATIC_FLOAT_PRECISION> &pars)
{
	//for series simulations, need to update the 3D output array independently with a read-add-write process
	Array4D<PRISMATIC_FLOAT_PRECISION> tmp_buffer(pars.output);
	std::string currentName = pars.currentTag;
	std::cout << "Reading scratch dataset " << currentName << " from " <<  "prismatic_scratch.h5" << std::endl;
	readRealDataSet_inOrder(tmp_buffer, "prismatic_scratch.h5", "scratch/"+currentName);
	for(auto i = 0; i < pars.output.size(); i++) tmp_buffer[i] += pars.output[i];

	H5::Group scratch = pars.scratchFile.openGroup("scratch");
	hsize_t mdims[4] = {pars.output.get_diml(), pars.output.get_dimk(), pars.output.get_dimj(), pars.output.get_dimi()};
	std::cout << "Writing scratch dataset " << currentName << " to " <<  "prismatic_scratch.h5" << std::endl;
	writeRealDataSet_inOrder(scratch, currentName, &tmp_buffer[0], mdims, 4);

	if(pars.meta.saveDPC_CoM)
	{
		Array4D<PRISMATIC_FLOAT_PRECISION> dpc_buffer(pars.DPC_CoM);
		currentName += "_DPC";
		std::cout << "Reading scratch dataset " << currentName << " from " <<  "prismatic_scratch.h5" << std::endl;
		readRealDataSet_inOrder(dpc_buffer, "prismatic_scratch.h5", "scratch/"+currentName);
		for(auto i = 0; i < pars.DPC_CoM.size(); i++) dpc_buffer[i] += pars.DPC_CoM[i];

		H5::Group scratch = pars.scratchFile.openGroup("scratch");
		hsize_t mdims_dpc[4] = {pars.DPC_CoM.get_diml(), pars.DPC_CoM.get_dimk(), pars.DPC_CoM.get_dimj(), pars.DPC_CoM.get_dimi()};
		std::cout << "MDIMS: " << mdims_dpc[0] << " " << mdims_dpc[1] << " " << mdims_dpc[2] << " " << mdims_dpc[3] << std::endl;

		std::cout << "Writing scratch dataset " << currentName << " to " <<  "prismatic_scratch.h5" << std::endl;

		writeRealDataSet_inOrder(scratch, currentName, &dpc_buffer[0], mdims_dpc, 4);

	}
	scratch.close();

};


} //namespace Prismatic