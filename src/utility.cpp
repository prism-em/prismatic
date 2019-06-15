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
#ifdef _WIN32
   #include <io.h> 
   #define access    _access_s
#else
   #include <unistd.h>
#endif


namespace Prismatic {



	std::pair<Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> >, Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > >
	upsamplePRISMProbe(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > probe,
	                   const long dimj, const long dimi, long ys, long xs) {
		Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > realspace_probe;
		Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > buffer_probe;
		Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > kspace_probe;

		buffer_probe = zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION> >({{(size_t)dimj, (size_t)dimi}});
//		std::cout << "dimj = " << dimj << std::endl;
		long ncy = probe.get_dimj() / 2;
		long ncx = probe.get_dimi() / 2;
		for (auto j = 0; j < probe.get_dimj(); ++j) {
			for (auto i = 0; i < probe.get_dimi(); ++i) {
                buffer_probe.at( (dimj + ((j - ncy + ys) % dimj)) % dimj,
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

	PRISMATIC_FLOAT_PRECISION computePearsonCorrelation(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > left,
	                                                Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > right){
		PRISMATIC_FLOAT_PRECISION m1, m2, sigma1, sigma2, R;
		m1=m2=sigma1=sigma2=R=0;

		for (auto &i:left) m1 += std::abs(i);
		for (auto &i:right)m2 += std::abs(i);

		m1 /= (left.size());
		m2 /= (right.size());

		for (auto &i:left)sigma1  += pow(std::abs(i)-m1, 2);
		for (auto &i:right)sigma2 += pow(std::abs(i)-m2, 2);

		sigma1 /= (left.size());
		sigma2 /= (right.size());

		sigma1 = sqrt(sigma1);
		sigma2 = sqrt(sigma2);
		for (auto i = 0; i < std::min(left.size(), right.size()); ++i){
			R = R + (std::abs(left[i]) - m1) * (std::abs(right[i]) - m2);
		}
		R/=sqrt(left.size()*right.size());
		return R / (sigma1 * sigma2);
	}
	PRISMATIC_FLOAT_PRECISION computeRfactor(Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > left,
	                                     Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION> > right){
		PRISMATIC_FLOAT_PRECISION accum, diffs;
		accum = diffs = 0;
		for (auto i = 0; i < std::min(left.size(), right.size()); ++i){
			diffs += std::abs(left[i] - right[i]);
			accum += std::abs(left[i]);
		}
		return diffs / accum;
	}

	std::string remove_extension(const std::string& filename) {
      size_t lastdot = filename.find_last_of(".");
      if (lastdot == std::string::npos) return filename;
      return filename.substr(0, lastdot);
	}

	bool testFilenameOutput(const std::string& filename){
		bool exists = testExist(filename);
		bool write_ok = testWrite(filename);
		//Check if file already exists and if we can write to it
		if(exists && write_ok){
		std::cout<<"Warning "<<filename<<" already exists and will be overwritten"<<std::endl;
			return true;
		}
		else if(exists && !write_ok){
			std::cout<<filename<<" isn't an accessible write destination"<<std::endl;
			return false;
		}
		else{
			//If the file does not exist, check to see if we can open a file of that name
			std::ofstream f(filename, std::ios::binary |std::ios::out);
			if(f){
				//If we can open such a file, close the file and delete it.
				f.close();
				std::remove(filename.c_str());
				return true;
			}
			else{
				std::cout<<filename<<" isn't an accessible write destination"<<std::endl;
				return false;
			}
		}


	}

	int testWrite(const std::string& filename){
        int answer = access(filename.c_str(),02); //W_OK = 02
        return answer;
    }

    int testExist(const std::string& filename){
        int answer = access(filename.c_str(),00); //F_OK == 00
        return answer;
    }

	void setupOutputFile(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars){
		//set version attributes
		int maj_data[1] = {0};
		int min_data[1] = {4};
		hsize_t attr_dims[1] = {1};

		H5::DataSpace attr_dataspace(1,attr_dims);

		H5::Attribute maj_attr = pars.outputFile.createAttribute("version_major=0",H5::PredType::NATIVE_INT,attr_dataspace);
		H5::Attribute min_attr = pars.outputFile.createAttribute("version_minor=4",H5::PredType::NATIVE_INT,attr_dataspace);

		maj_attr.write(H5::PredType::NATIVE_INT, maj_data);
		min_attr.write(H5::PredType::NATIVE_INT, min_data);

		//create main groups
		H5::Group simulation(pars.outputFile.createGroup("/4DSTEM_experiment"));

		//data groups
		H5::Group data(simulation.createGroup("data"));
		H5::Group datacubes(data.createGroup("datacubes"));
		H5::Group dslices(data.createGroup("diffractionslices"));
		H5::Group rslices(data.createGroup("realslices"));
		H5::Group pointlists(data.createGroup("pointlists")); //point lists and point list arrays are not used in prismatic
		H5::Group plarrays(data.creaeGroup("pointlistarrays")); //included here to maintain consistency with format

		//log group
		H5::Group log(simulation.createGroup("log"));

		//metadata groups
		H5::Group metadata(simulation.createGroup("metadata"));
		H5::Group original(metadata.createGroup("original"));
		H5::Group microscope(metadata.createGroup("microscope"));
		H5::Group sample(metadata.createGroup("sample"));
		H5::Group user(metadata.createGroup("user"));
		H5::Group calibration(metadata.createGroup("calibration"));
		H5::Group comments(metadata.createGroup("comments"));

	}

	void setup4DOutput(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> pars, const size_t numLayers){
		H5::Group datacubes = pars.outputFile.openGroup("data/datacubes");

		//shared properties 
		std::string base_name = "CBED_array_slice";
		hsize_t attr_dims[1] = {1};
		hsize_t data_dims[4] = {pars.xp.size(),pars.yp.size(),pars.psiProbeInit.get_dimi()/2,pars.psiProbeInit.get_dimj()/2};
		hsize_t rx_dim[1] = {pars.xp.size()};
		hsize_t ry_dim[1] = {pars.yp.size()};

		if(pars.meta.algorithm == Prismatic::Algorithm::Multislice){
			hsize_t qx_dim[1] = {pars.psiProbeInit.get_dimi()/2};
			hsize_t qy_dim[1] = {pars.psiProbeInit.get_dimj()/2};
		}else{
			hsize_t qx_dim[1] = {pars.psiProbeInit.get_dimi()};
			hsize_t qy_dim[1] = {pars.psiProbeInit.get_dimj()};			
		}

		for(auto n = 0; n < numLayers; n++){
			//create slice group
			std::string nth_name = base_name + std::to_string(n);
			H5::Group CBED_slice_n(datacubes.createGroup(nth_name.c_str()));
			
			//write group type attribute
			H5::DataSpace attr_dataspace(1,attr_dims);
			H5::Attribute emd_group_type = CBED_slice_n.createAttribute("emd_group_type",H5::PredType::NATIVE_INT,attr_dataspace);
			int group_type[1] = {1};
			emd_group_type.write(H5::PredType::NATIVE_INT, group_type);	

			//create dataset
			H5::DataSpace mspace(4,dataDims); //rank is 4
			H5::DataSet CBED_data = CBED_slice_n.createDataSet("datacube",H5::PredType::NATIVE_FLOAT,mspace);			
			mspace.close();

			//write dimensions
			//TODO: fftshift the dimensions so they are consistent with image; write the fftshift1 function to do so
			//TODO: qx, qy truncate so that match antialiasing filter if using multislice
			H5::DataSpace str_name_ds(H5S_SCALAR);
			H5::StrType strdatatype(H5::PredType::C_S1,256);

			H5::DataSpace dim1_mspace(1,rx_dim);
			H5::DataSpace dim2_mspace(1,ry_dim);
			H5::DataSpace dim3_mspace(1,qx_dim);
			H5::DataSpace dim4_mspace(1,qy_dim);

			H5::DataSet dim1 = CBED_slice_n.createDataSet("dim1",H5::PredType::NATIVE_FLOAT,dim1_mspace);
			H5::DataSet dim2 = CBED_slice_n.createDataSet("dim2",H5::PredType::NATIVE_FLOAT,dim2_mspace);
			H5::DataSet dim3 = CBED_slice_n.createDataSet("dim3",H5::PredType::NATIVE_FLOAT,dim3_mspace);
			H5::DataSet dim4 = CBED_slice_n.createDataSet("dim4",H5::PredType::NATIVE_FLOAT,dim4_mspace);

			H5::DataSpace dim1_fspace = dim1.getspace();
			H5::DataSpace dim2_fspace = dim2.getspace();
			H5::DataSpace dim3_fspace = dim3.getspace();
			H5::DataSpace dim4_fspace = dim4.getspace();

			dim1.write(&pars.xp[0],H5::PredType::NATIVE_FLOAT,dim1_mspace,dim1_fspace);
			dim2.write(&pars.yp[0],H5::PredType::NATIVE_FLOAT,dim2_mspace,dim2_fspace);
			dim3.write(&pars.qx[0],H5::PredType::NATIVE_FLOAT,dim3_mspace,dim3_fspace);
			dim4.write(&pars.qy[0],H5::PredType::NATIVE_FLOAT,dim4_mspace,dim4_fspace);
			
			//dimension attributes
			const H5std_string dim1_name_str("R_x");
			const H5std_string dim2_name_str("R_y");
			const H5std_string dim3_name_str("Q_x");
			const H5std_string dim4_name_str("Q_y");

			H5::Attribute dim1_name = dim1.createAttribute("name",strdatatype,str_name_ds);
			H5::Attribute dim2_name = dim2.createAttribute("name",strdatatype,str_name_ds);
			H5::Attribute dim3_name = dim3.createAttribute("name",strdatatype,str_name_ds);
			H5::Attribute dim4_name = dim4.createAttribute("name",strdatatype,str_name_ds);

			dim1_name.write(strdatatype,dim1_name_str);
			dim2_name.write(strdatatype,dim2_name_str);
			dim3_name.write(strdatatype,dim3_name_str);
			dim4_name.write(strdatatype,dim4_name_str);

			const H5std_string dim1_unit_str("[n_m]");
			const H5std_string dim2_unit_str("[n_m]");
			const H5std_string dim3_unit_str("[n_m^-1]");
			const H5std_string dim3_unit_str("[n_m^-1]");

			H5::Attribute dim1_unit = dim1.createAttribute("units",strdatatype,str_name_ds);
			H5::Attribute dim2_unit = dim2.createAttribute("units",strdatatype,str_name_ds);
			H5::Attribute dim3_unit = dim3.createAttribute("units",strdatatype,str_name_ds);
			H5::Attribute dim4_unit = dim4.createAttribute("units",strdatatype,str_name_ds);

			dim1_unit.write(strdatatype,dim1_unit_str);
			dim2_unit.write(strdatatype,dim2_unit_str);
			dim3_unit.write(strdatatype,dim3_unit_str);
			dim4_unit.write(strdatatype,dim4_unit_str);

			CBED_slice_n.close();
		}

		datacubes.close();
	};

	void writeDatacube3D(H5::DataSet dataset, const PRISMATIC_FLOAT_PRECISION* buffer, const hsize_t* mdims){
		//set up file and memory spaces
		H5::DataSpace fspace = dataset.getSpace(); //all 3D cubes will write full buffer at once
		H5::DataSpace mspace(3,mdims); //rank = 3

		dataset.write(buffer,H5::PredType::NATIVE_FLOAT,mspace,fspace);

		fspace.close();
		mspace.close();
	};

	void writeDatacube4D(H5::DataSet dataset, const PRISMATIC_FLOAT_PRECISION* buffer, const hsize_t* mdims, const hsize_t* offset){
		//set up file and memory spaces
		H5::DataSpace fspace = dataset.getSpace();
		H5::DataSpace mspace(4,mdims); //rank = 4
		
		fspace.selectHyperslab(H5S_SELECT_SET,mdims,offset);
		dataset.write(buffer,H5::PredType::NATIVE_FLOAT,mspace,fspace);

		fspace.close();
		mspace.close();
	};

}
