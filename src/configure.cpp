//
// Created by AJ Pryor on 3/5/17.
//
#include "configure.h"
#include "PRISM_entry.h"
#include "Multislice_entry.h"
#include "Multislice.h"
#include <iostream>
#include "PRISM01.h"
#include "PRISM02.h"
#include "PRISM03.h"
//#define PRISM_ENABLE_GPU



#ifdef PRISM_ENABLE_GPU
#include "Multislice.cuh"
#include "PRISM02.cuh"
#include "PRISM03.cuh"
#include "utility.cuh"
#include "Multislice_entry.h"
#endif //PRISM_ENABLE_GPU
namespace PRISM {
	entry_func execute_plan;
	ms_output_func buildMultisliceOutput;
	prism_output_func buildPRISMOutput;
	format_output_func formatOutput_CPU;
	fill_Scompact_func fill_Scompact;

#ifdef PRISM_ENABLE_GPU
    template <class T>
    StreamingMode transferMethodAutoChooser(PRISM::Metadata<T>& meta){
        // query all devices and choose based on the minimum compute capability
        // get the total memory on device
        // based on whether the algo is prism or multislice, estimate the largest array and trigger streaming
        // safely far from this limit
        // For PRISM: limit should be the maximum size of the potential or Scompact
        // which is imageSize/2 x imageSize/2 * numberBeams of type complex float/double
        // numberBeams can be calculated again here... not a big deal
        // For Multislice: limit is the potential array which is imageSize x imageSize x numPlanes of type float/double

        //placeholder for now
        return PRISM::StreamingMode::Stream;
    }
	format_output_func_GPU formatOutput_GPU;

#endif
	void configure(Metadata<PRISM_FLOAT_PRECISION>& meta) {
		formatOutput_CPU = formatOutput_CPU_integrate;
#ifdef PRISM_ENABLE_GPU
		formatOutput_GPU = formatOutput_GPU_integrate;
#endif
		if (meta.algorithm == Algorithm::PRISM) {
			std::cout << "Execution plan: PRISM w/ single FP configuration" << std::endl;
			execute_plan = PRISM_entry;
#ifdef PRISM_ENABLE_GPU
            if (meta.transfer_mode == PRISM::StreamingMode::Auto){
            	meta.transfer_mode = transferMethodAutoChooser(meta);
            }
			if (meta.transfer_mode == PRISM::StreamingMode::Stream) {
				cout << "Using streaming method\n";
				fill_Scompact = fill_Scompact_GPU_streaming;
				buildPRISMOutput = buildPRISMOutput_GPU_streaming;
			} else {
				cout << "Using single transfer method\n";
				fill_Scompact = fill_Scompact_GPU_singlexfer;
				buildPRISMOutput = buildPRISMOutput_GPU_singlexfer;
			}
#else
			fill_Scompact = fill_Scompact_CPUOnly;
			buildPRISMOutput = buildPRISMOutput_CPUOnly;
#endif //PRISM_ENABLE_GPU
		} else if (meta.algorithm == Algorithm::Multislice) {
			std::cout << "Execution plan: Multislice w/ single FP configuration" << std::endl;
			execute_plan = Multislice_entry;
#ifdef PRISM_ENABLE_GPU
			std::cout << "Using GPU codes" << std::endl;
			if (meta.transfer_mode == PRISM::StreamingMode::Auto){
				meta.transfer_mode = transferMethodAutoChooser(meta);
			}
			if (meta.transfer_mode == PRISM::StreamingMode::Stream) {
				cout << "Using streaming method\n";
				buildMultisliceOutput = buildMultisliceOutput_GPU_streaming;
			} else {
				cout << "Using single transfer method\n";
				buildMultisliceOutput = buildMultisliceOutput_GPU_singlexfer;
			}

#else
			buildMultisliceOutput = buildMultisliceOutput_CPUOnly;
#endif //PRISM_ENABLE_GPU
		}
	}
}
