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

#include "utility.cuh"
#include "utility.h"
#include "cuComplex.h"
#include <iostream>
#include <sstream>
#include <mutex>
#include "fileIO.h"
#include "fileIO.cuh"

void formatOutput_GPU_integrate(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
                                PRISMATIC_FLOAT_PRECISION *psiIntensity_ds,
                                const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
                                PRISMATIC_FLOAT_PRECISION *output_ph,
								PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
								const PRISMATIC_FLOAT_PRECISION* qya_d,
								const PRISMATIC_FLOAT_PRECISION* qxa_d,
								const size_t currentSlice,
                                const size_t ay,
                                const size_t ax,
                                const size_t& dimj,
                                const size_t& dimi,
                                const cudaStream_t& stream,
                                const long& scale) {

	//save 4D output if applicable
    if (pars.meta.save4DOutput)
    {
		// This section could be improved. It currently makes a new 2D array, copies to it, and
		// then saves the image. This allocates arrays multiple times unneccessarily, and the allocated
		// memory isn't pinned, so the memcpy is not asynchronous.
		//std::string section4DFilename = generateFilename(pars, currentSlice, ay, ax);
		
		Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> currentImage = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>(
				{{pars.psiProbeInit.get_dimj(), pars.psiProbeInit.get_dimi()}});
		cudaErrchk(cudaMemcpyAsync(&currentImage[0],
		                           psiIntensity_ds,
		                           pars.psiProbeInit.size() * sizeof(PRISMATIC_FLOAT_PRECISION),
		                           cudaMemcpyDeviceToHost,
								   stream));
								   
		//Need to scale the output by the square of the PRISM interpolation factor 
		// std::unique_lock<std::mutex> HDF5_gatekeeper(Prismatic::HDF5_lock);

		currentImage *= pars.scale;
		std::string nameString = "/4DSTEM_simulation/data/datacubes/CBED_array_depth" + Prismatic::getDigitString(currentSlice);
		nameString += pars.currentTag;

		hsize_t offset[4] = {ax,ay,0,0}; //order by ax, ay so that aligns with py4DSTEM
        PRISMATIC_FLOAT_PRECISION numFP = pars.meta.numFP;
        
        if(pars.meta.crop4DOutput)
        {
            Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> finalImage = cropOutput(currentImage, pars);
            hsize_t mdims[4] = {1,1,finalImage.get_dimj(),finalImage.get_dimi()};
            Prismatic::writeDatacube4D(pars, &finalImage[0],&pars.cbed_buffer[0],mdims,offset,numFP,nameString.c_str());
        }
        else
        {

            if (pars.meta.algorithm == Prismatic::Algorithm::Multislice){
                Prismatic::Array2D<PRISMATIC_FLOAT_PRECISION> finalImage = Prismatic::zeros_ND<2, PRISMATIC_FLOAT_PRECISION>(
                    {{pars.psiProbeInit.get_dimi()/2,pars.psiProbeInit.get_dimj()/2}});
                    {
                        long offset_x = pars.psiProbeInit.get_dimi() / 4;
                        long offset_y = pars.psiProbeInit.get_dimj() / 4;
                        long ndimy = (long) pars.psiProbeInit.get_dimj();
                        long ndimx = (long) pars.psiProbeInit.get_dimi();
                        for (long y = 0; y < pars.psiProbeInit.get_dimj() / 2; ++y) {
                            for (long x = 0; x < pars.psiProbeInit.get_dimi() / 2; ++x) {
                                finalImage.at(x, y) = currentImage.at(((y - offset_y) % ndimy + ndimy) % ndimy,
                                ((x - offset_x) % ndimx + ndimx) % ndimx);
                            }
                        }
                    }
                    
                    hsize_t mdims[4] = {1,1, finalImage.get_dimj(), finalImage.get_dimi()};
                    Prismatic::writeDatacube4D(pars, &finalImage[0],&pars.cbed_buffer[0],mdims,offset,numFP,nameString.c_str());
                }else{                     
                    currentImage = fftshift2_flip(currentImage);
                    hsize_t mdims[4] = {1,1,currentImage.get_dimj(), currentImage.get_dimi()};
                    Prismatic::writeDatacube4D(pars, &currentImage[0],&pars.cbed_buffer[0],mdims,offset,numFP,nameString.c_str());
                }
        }
    }

	size_t num_integration_bins = pars.detectorAngles.size();
	setAll <<< (num_integration_bins - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>>
	                                                                            (integratedOutput_ds, 0, num_integration_bins);

	integrateDetector <<< (dimj * dimi - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>>
	                                                                              (psiIntensity_ds, alphaInd_d, integratedOutput_ds,
			                                                                              dimj *
			                                                                              dimi, num_integration_bins);

	multiply_arr_scalar <<< (dimj * dimi - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>>
	                                                                                (integratedOutput_ds, scale, num_integration_bins);

	cudaErrchk(cudaMemcpyAsync(output_ph, integratedOutput_ds,
	                           num_integration_bins * sizeof(PRISMATIC_FLOAT_PRECISION),
	                           cudaMemcpyDeviceToHost, stream));

	//	 wait for the copy to complete and then copy on the host. Other host threads exist doing work so this wait isn't costing anything
	cudaErrchk(cudaStreamSynchronize(stream));
	const size_t stack_start_offset =
			currentSlice * pars.output.get_dimk() * pars.output.get_dimj() * pars.output.get_dimi() + ay * pars.output.get_dimj() * pars.output.get_dimi() + ax * pars.output.get_dimi();
	memcpy(&pars.output[stack_start_offset], output_ph, num_integration_bins * sizeof(PRISMATIC_FLOAT_PRECISION));
	
    if(pars.meta.saveDPC_CoM)
    {
		//device variables
		PRISMATIC_FLOAT_PRECISION *num_qx_d;
		PRISMATIC_FLOAT_PRECISION *num_qy_d;
		PRISMATIC_FLOAT_PRECISION *denominator_d;
		cudaErrchk(cudaMallocManaged(&num_qx_d, 1*sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocManaged(&num_qy_d, 1*sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocManaged(&denominator_d, 1*sizeof(PRISMATIC_FLOAT_PRECISION)));

		//host variables
		PRISMATIC_FLOAT_PRECISION *num_qx_h = new PRISMATIC_FLOAT_PRECISION[1];
		PRISMATIC_FLOAT_PRECISION *num_qy_h = new PRISMATIC_FLOAT_PRECISION[1];
		PRISMATIC_FLOAT_PRECISION *denominator_h = new PRISMATIC_FLOAT_PRECISION[1];
		num_qx_h[0] = 0.0;
		num_qy_h[0] = 0.0;
		denominator_h[0] = 0.0;

		//initialize device variables
		cudaErrchk(cudaMemcpyAsync(num_qx_d,&num_qx_h[0],1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyHostToDevice));
		cudaErrchk(cudaMemcpyAsync(num_qy_d,&num_qy_h[0],1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyHostToDevice));
		cudaErrchk(cudaMemcpyAsync(denominator_d,&denominator_h[0],1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyHostToDevice));
		
		//reduce in X
		DPC_numerator_reduce <<< (dimj * dimi - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>>
		(psiIntensity_ds,qxa_d, num_qx_d, dimj * dimi);
		
		//reduce in Y
		DPC_numerator_reduce <<< (dimj * dimi - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>>
		(psiIntensity_ds,qya_d, num_qy_d, dimj * dimi);
		
		DPC_denominator_reduce <<< (dimj * dimi - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>> (psiIntensity_ds, denominator_d, dimj*dimi);
		
		//copy back to host
		cudaErrchk(cudaMemcpyAsync(&num_qx_h[0],num_qx_d,1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyDeviceToHost));
		cudaErrchk(cudaMemcpyAsync(&num_qy_h[0],num_qy_d,1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyDeviceToHost));
		cudaErrchk(cudaMemcpyAsync(&denominator_h[0],denominator_d,1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyDeviceToHost));

		PRISMATIC_FLOAT_PRECISION DPC_CoM[2];
		DPC_CoM[0] = num_qx_h[0]/denominator_h[0]; //measurement at ax,ay of CoM w.r.t. qx
		DPC_CoM[1] = num_qy_h[0]/denominator_h[0]; //measurement at ax,ay of CoM w.r.t. qy

		//copy to memory and free variables
		const size_t dpc_stack_offset = 
				currentSlice*pars.DPC_CoM.get_dimk() * pars.DPC_CoM.get_dimj() * pars.DPC_CoM.get_dimi() + ay * pars.DPC_CoM.get_dimj() * pars.DPC_CoM.get_dimi() + ax * pars.DPC_CoM.get_dimi();
		memcpy(&pars.DPC_CoM[dpc_stack_offset],&DPC_CoM[0],2*sizeof(PRISMATIC_FLOAT_PRECISION));
		cudaErrchk(cudaFree(num_qx_d));
		cudaErrchk(cudaFree(num_qy_d));
		cudaErrchk(cudaFree(denominator_d));
		free(num_qx_h);
		free(num_qy_h);
		free(denominator_h);
	}
}

void formatOutput_GPU_c_integrate(Prismatic::Parameters<PRISMATIC_FLOAT_PRECISION> &pars,
									PRISMATIC_CUDA_COMPLEX_FLOAT *psi,
									PRISMATIC_FLOAT_PRECISION *psiIntensity_ds,
									const PRISMATIC_FLOAT_PRECISION *alphaInd_d,
									PRISMATIC_FLOAT_PRECISION *output_ph,
									PRISMATIC_FLOAT_PRECISION *integratedOutput_ds,
									const PRISMATIC_FLOAT_PRECISION* qya_d,
									const PRISMATIC_FLOAT_PRECISION* qxa_d,
									const size_t currentSlice,
									const size_t ay,
									const size_t ax,
									const size_t& dimj,
									const size_t& dimi,
									const cudaStream_t& stream,
									const long& scale){

	//save 4D output if applicable
    if (pars.meta.save4DOutput)
    {
		// This section could be improved. It currently makes a new 2D array, copies to it, and
		// then saves the image. This allocates arrays multiple times unneccessarily, and the allocated
		// memory isn't pinned, so the memcpy is not asynchronous.
		Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> currentImage = 
						Prismatic::zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>({{pars.psiProbeInit.get_dimj(), pars.psiProbeInit.get_dimi()}});
		cudaErrchk(cudaMemcpyAsync(&currentImage[0],
									psi,
									pars.psiProbeInit.size() * sizeof(std::complex<PRISMATIC_FLOAT_PRECISION>),
									cudaMemcpyDeviceToHost,
									stream));
									
		// Need to scale the output by the square of the PRISM interpolation factor 
		currentImage *= sqrt(pars.scale);
		std::string nameString = "/4DSTEM_simulation/data/datacubes/CBED_array_depth" + Prismatic::getDigitString(currentSlice);
		nameString += pars.currentTag + "_fp" + Prismatic::getDigitString(pars.meta.fpNum);
		

		Prismatic::Array2D<std::complex<PRISMATIC_FLOAT_PRECISION>> finalImage;
		
		hsize_t offset[4] = {ax,ay,0,0}; //order by ax, ay so that aligns with py4DSTEM
		PRISMATIC_FLOAT_PRECISION numFP = pars.meta.numFP;
		
        if(pars.meta.crop4DOutput)
        {
            finalImage = cropOutput(currentImage, pars);
            hsize_t mdims[4] = {1,1,finalImage.get_dimj(),finalImage.get_dimi()};
            Prismatic::writeDatacube4D(pars, &finalImage[0],&pars.cbed_buffer_c[0],mdims,offset,numFP,nameString.c_str());
        }
        else
        {
            if (pars.meta.algorithm == Prismatic::Algorithm::Multislice){
                finalImage = Prismatic::zeros_ND<2, std::complex<PRISMATIC_FLOAT_PRECISION>>(
                    {{pars.psiProbeInit.get_dimi()/2,pars.psiProbeInit.get_dimj()/2}});
                    {
                        long offset_x = pars.psiProbeInit.get_dimi() / 4;
                        long offset_y = pars.psiProbeInit.get_dimj() / 4;
                        long ndimy = (long) pars.psiProbeInit.get_dimj();
                        long ndimx = (long) pars.psiProbeInit.get_dimi();
                        for (long y = 0; y < pars.psiProbeInit.get_dimj() / 2; ++y) {
                            for (long x = 0; x < pars.psiProbeInit.get_dimi() / 2; ++x) {
                                finalImage.at(x, y) = currentImage.at(((y - offset_y) % ndimy + ndimy) % ndimy,
                                ((x - offset_x) % ndimx + ndimx) % ndimx);
                            }
                        }
                    }
                    
                    hsize_t mdims[4] = {1,1, finalImage.get_dimj(), finalImage.get_dimi()};
                    Prismatic::writeDatacube4D(pars, &finalImage[0],&pars.cbed_buffer_c[0],mdims,offset,numFP,nameString.c_str());
                }else{                     
                    currentImage = fftshift2_flip(currentImage);
                    hsize_t mdims[4] = {1, 1, currentImage.get_dimj(), currentImage.get_dimi()};
                    Prismatic::writeDatacube4D(pars, &currentImage[0],&pars.cbed_buffer_c[0],mdims,offset,numFP,nameString.c_str());
                }
        }
	}


	size_t num_integration_bins = pars.detectorAngles.size();
	setAll <<< (num_integration_bins - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>>
	                                                                            (integratedOutput_ds, 0, num_integration_bins);

	integrateDetector <<< (dimj * dimi - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>>
	                                                                              (psiIntensity_ds, alphaInd_d, integratedOutput_ds,
			                                                                              dimj *
			                                                                              dimi, num_integration_bins);

	multiply_arr_scalar <<< (dimj * dimi - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>>
	                                                                                (integratedOutput_ds, scale, num_integration_bins);

	cudaErrchk(cudaMemcpyAsync(output_ph, integratedOutput_ds,
	                           num_integration_bins * sizeof(PRISMATIC_FLOAT_PRECISION),
	                           cudaMemcpyDeviceToHost, stream));

	//	 wait for the copy to complete and then copy on the host. Other host threads exist doing work so this wait isn't costing anything
	cudaErrchk(cudaStreamSynchronize(stream));
	const size_t stack_start_offset =
			currentSlice * pars.output.get_dimk() * pars.output.get_dimj() * pars.output.get_dimi() + ay * pars.output.get_dimj() * pars.output.get_dimi() + ax * pars.output.get_dimi();
	memcpy(&pars.output[stack_start_offset], output_ph, num_integration_bins * sizeof(PRISMATIC_FLOAT_PRECISION));
	
    if(pars.meta.saveDPC_CoM)
    {
		//device variables
		PRISMATIC_FLOAT_PRECISION *num_qx_d;
		PRISMATIC_FLOAT_PRECISION *num_qy_d;
		PRISMATIC_FLOAT_PRECISION *denominator_d;
		cudaErrchk(cudaMallocManaged(&num_qx_d, 1*sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocManaged(&num_qy_d, 1*sizeof(PRISMATIC_FLOAT_PRECISION)));
		cudaErrchk(cudaMallocManaged(&denominator_d, 1*sizeof(PRISMATIC_FLOAT_PRECISION)));

		//host variables
		PRISMATIC_FLOAT_PRECISION *num_qx_h = new PRISMATIC_FLOAT_PRECISION[1];
		PRISMATIC_FLOAT_PRECISION *num_qy_h = new PRISMATIC_FLOAT_PRECISION[1];
		PRISMATIC_FLOAT_PRECISION *denominator_h = new PRISMATIC_FLOAT_PRECISION[1];
		num_qx_h[0] = 0.0;
		num_qy_h[0] = 0.0;
		denominator_h[0] = 0.0;

		//initialize device variables
		cudaErrchk(cudaMemcpyAsync(num_qx_d,&num_qx_h[0],1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyHostToDevice));
		cudaErrchk(cudaMemcpyAsync(num_qy_d,&num_qy_h[0],1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyHostToDevice));
		cudaErrchk(cudaMemcpyAsync(denominator_d,&denominator_h[0],1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyHostToDevice));
		
		//reduce in X
		DPC_numerator_reduce <<< (dimj * dimi - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>>
		(psiIntensity_ds,qxa_d, num_qx_d, dimj * dimi);
		
		//reduce in Y
		DPC_numerator_reduce <<< (dimj * dimi - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>>
		(psiIntensity_ds,qya_d, num_qy_d, dimj * dimi);
		
		DPC_denominator_reduce <<< (dimj * dimi - 1) / BLOCK_SIZE1D + 1, BLOCK_SIZE1D, 0, stream >>> (psiIntensity_ds, denominator_d, dimj*dimi);
		
		//copy back to host
		cudaErrchk(cudaMemcpyAsync(&num_qx_h[0],num_qx_d,1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyDeviceToHost));
		cudaErrchk(cudaMemcpyAsync(&num_qy_h[0],num_qy_d,1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyDeviceToHost));
		cudaErrchk(cudaMemcpyAsync(&denominator_h[0],denominator_d,1*sizeof(PRISMATIC_FLOAT_PRECISION),cudaMemcpyDeviceToHost));

		PRISMATIC_FLOAT_PRECISION DPC_CoM[2];
		DPC_CoM[0] = num_qx_h[0]/denominator_h[0]; //measurement at ax,ay of CoM w.r.t. qx
		DPC_CoM[1] = num_qy_h[0]/denominator_h[0]; //measurement at ax,ay of CoM w.r.t. qy

		//copy to memory and free variables
		const size_t dpc_stack_offset = 
				currentSlice*pars.DPC_CoM.get_dimk() * pars.DPC_CoM.get_dimj() * pars.DPC_CoM.get_dimi() + ay * pars.DPC_CoM.get_dimj() * pars.DPC_CoM.get_dimi() + ax * pars.DPC_CoM.get_dimi();
		memcpy(&pars.DPC_CoM[dpc_stack_offset],&DPC_CoM[0],2*sizeof(PRISMATIC_FLOAT_PRECISION));
		cudaErrchk(cudaFree(num_qx_d));
		cudaErrchk(cudaFree(num_qy_d));
		cudaErrchk(cudaFree(denominator_d));
		free(num_qx_h);
		free(num_qy_h);
		free(denominator_h);
	}
}
