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

#include "parseInput.h"
#include <iostream>
#include <map>
#include <string>
#include <stdlib.h>

namespace Prismatic {
    using namespace std;

    void printHelp() {
        Metadata<PRISMATIC_FLOAT_PRECISION> defaults;
        std::cout << "Basic usage is prismatic -i filename [other options]" << std::endl;
        std::cout << "The following options are available with prismatic, each documented as long form (short form) *parameters* : description\n"
                "\n"
                "* --input-file (-i) filename : the filename containing the atomic coordinates, see www.prism-em.com/about for details (default: " << defaults.filenameAtoms << ")\n" <<
                "* --output-file(-o) filename : output filename (default: " << defaults.filenameOutput << ")\n" <<
                "* --interp-factor (-f) number : PRISM interpolation factor, used for both X and Y (default: " << defaults.interpolationFactorX << ")\n" <<
		        "* --interp-factor-x (-fx) number : PRISM interpolation factor in X (default: " << defaults.interpolationFactorX << ")\n" <<
		        "* --interp-factor-y (-fy) number : PRISM interpolation factor in Y (default: " << defaults.interpolationFactorY << ")\n" <<
                "* --num-threads (-j) value : number of CPU threads to use (default: " << defaults.numThreads << ")\n" <<
                "* --num-streams (-S) value : number of CUDA streams to create per GPU (default: " << defaults.numStreamsPerGPU << ")\n" <<
		        "* --num-gpus (-g) value : number of GPUs to use. A runtime check is performed to check how many are actually available, and the minimum of these two numbers is used. (default: " << defaults.numGPUs << ")\n" <<
                "* --slice-thickness (-s) thickness : thickness of each slice of projected potential (in Angstroms) (default: " << defaults.sliceThickness << ")\n" <<
		        "* --batch-size (-b) value : number of probes/beams to propagate simultaneously for both CPU and GPU workers. (default: " << defaults.batchSizeCPU << ")\n" <<
		        "* --batch-size-cpu (-bc) value : number of probes/beams to propagate simultaneously for CPU workers. (default: " << defaults.batchSizeCPU << ")\n" <<
		        "* --batch-size-gpu (-bg) value : number of probes/beams to propagate simultaneously for GPU workers. (default: " << defaults.batchSizeGPU << ")\n" <<
                "* --help(-h) : print information about the available options\n"
                "* --pixel-size (-p) pixel_size : size of simulated potential/probe X/Y pixel size (default: " << defaults.realspacePixelSize[0] << "). Note this is different from the size of a pixel in the output, which is determined by probe_stepX(Y)\n" <<
		        "* --pixel-size-x (-px) pixel_size : size of simulated potential/probe X pixel size (default: " << defaults.realspacePixelSize[1] << "). Note this is different from the size of a pixel in the output, which is determined by probe_stepX(Y)\n" <<
		        "* --pixel-size-y (-py) pixel_size : size of simulated potential/probe Y pixel size (default: " << defaults.realspacePixelSize[0] << "). Note this is different from the size of a pixel in the output, which is determined by probe_stepX(Y)\n" <<
		        "* --detector-angle-step (-d) step_size : angular step size for detector integration bins (in mrad) (default: " << (1000 * defaults.detectorAngleStep) << ")\n" <<
                "* --cell-dimension (-c) x y z : size of sample in x, y, z directions (in Angstroms) (default: " << defaults.cellDim[2] << " " << defaults.cellDim[1] << " " << defaults.cellDim[0] << ")\n" <<
		        "* --tile-uc (-t) x y z : tile the unit cell x, y, z number of times in x, y, z directions, respectively (default: " << defaults.tileX << " " << defaults.tileY << " " << defaults.tileZ << ")\n" <<
                "* --algorithm (-a) p/m : the simulation algorithm to use, either (p)rism or (m)ultislice (default: PRISM)\n" <<
                "* --energy (-E) value : the energy of the electron beam (in keV) (default: " << defaults.E0/1000 << ")\n" <<
                "* --alpha-max (-A) angle : the maximum probe angle to consider (in mrad) (default: " << 1000*defaults.alphaBeamMax << ")\n" <<
                "* --potential-bound (-P) value : the maximum radius from the center of each atom to compute the potental (in Angstroms) (default: " << defaults.potBound << ")\n" <<
                "* --also-do-cpu-work (-C) bool=true : boolean value used to determine whether or not to also create CPU workers in addition to GPU ones (default: 1)\n" <<
                "* --streaming-mode 0/1 : boolean value to force code to use (true) or not use (false) streaming versions of GPU codes. The default behavior is to estimate the needed memory from input parameters and choose automatically. (default: Auto)\n" <<
                "* --probe-step (-r) step_size : step size of the probe for both X and Y directions (in Angstroms) (default: " << defaults.probeStepX << ")\n" <<
                "* --probe-step-x (-rx) step_size : step size of the probe in X direction (in Angstroms) (default: " << defaults.probeStepX << ")\n" <<
                "* --probe-step-y (-ry) step_size : step size of the probe in Y direction (in Angstroms) (default: " << defaults.probeStepY << ")\n" <<
		        "* --random-seed (-rs) step_size : random integer number seed\n"
	            "* --probe-xtilt (-tx) value : probe X tilt (in mrad) (default: " << defaults.probeXtilt << ")\n" <<
                "* --probe-ytilt (-ty) value : probe X tilt (in mrad) (default: " << defaults.probeYtilt << ")\n" <<
                "* --probe-defocus (-df) value : probe defocus (in mrad) (default: " << defaults.probeDefocus << ")\n" <<
                "* -C3 value : microscope C3 aberration constant (in Angstrom) (default: " << defaults.C3 << ")\n" <<
                "* -C5 value : microscope C5 aberration constant (in Angstrom) (default: " << defaults.C5 << ")\n" <<
                "* --probe-semiangle (-sa) value : maximum probe semiangle (in mrad) (default: " << 1000*defaults.probeSemiangle << ")\n" <<
                "* --scan-window-x (-wx) min max : size of the window to scan the probe in X (in fractional coordinates between 0 and 1) (default: " << defaults.scanWindowXMin << " " << defaults.scanWindowXMax << ")\n" <<
                "* --scan-window-y (-wy) min max : size of the window to scan the probe in Y (in fractional coordinates between 0 and 1) (default: " << defaults.scanWindowYMin << " " << defaults.scanWindowYMax << ")\n" <<
                "* --num-FP (-F) value : number of frozen phonon configurations to calculate (default: " << defaults.numFP << ")\n" <<
		        "* --thermal-effects (-te) bool : whether or not to include Debye-Waller factors (thermal effects) (default: True)\n" <<
                "* --occupancy (-oc) bool : whether or not to consider occupancy values for likelihood of atoms existing at each site (default: True)\n" <<
		        "* --save-2D-output (-2D) ang_min ang_max : save the 2D STEM image integrated between ang_min and ang_max (in mrads) (default: Off)\n" <<
	            "* --save-3D-output (-3D) bool=true : Also save the 3D output at the detector for each probe (3D output mode) (default: On)\n" <<
                "* --save-4D-output (-4D) bool=false : Also save the 4D output at the detector for each probe (4D output mode) (default: Off)\n";
    }

    bool parse_a(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No algorithm provided for -a (syntax is -a algorithm). Choices are (m)ultislice or (p)rism\n";
            return false;
        }
        std::string algo = std::string((*argv)[1]);
        if (algo == "m" | algo == "multislice"){
            meta.algorithm = Prismatic::Algorithm::Multislice;
        } else if (algo == "p" | algo == "prism"){
            meta.algorithm = Prismatic::Algorithm::PRISM;
        } else {
            cout << "Unrecognized algorithm \"" << (*argv)[1] << "\"\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_A(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No maximum probe angle provided for -A (syntax is -A angle (in mrad))\n";
            return false;
        }
        if ( (meta.alphaBeamMax =  ( (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) / 1000 ) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for maximum probe angle (syntax is -A angle (in mrad))\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };


	bool parse_b(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No batch size provided for -b (syntax is -b batch_size)\n";
			return false;
		}
		if ( (meta.batchSizeTargetCPU = atoi((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for batch size (syntax is -b batch_size)\n";
			return false;
		}
		meta.batchSizeTargetGPU = meta.batchSizeTargetCPU;
		meta.batchSizeGPU = meta.batchSizeTargetGPU;
		meta.batchSizeCPU = meta.batchSizeTargetCPU;
		argc-=2;
		argv[0]+=2;
		return true;
	};

	bool parse_bc(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No batch size provided for -bc (syntax is -bc batch_size)\n";
			return false;
		}
		if ( (meta.batchSizeTargetCPU = atoi((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for CPU batch size (syntax is -bc batch_size)\n";
			return false;
		}
		meta.batchSizeCPU = meta.batchSizeTargetCPU;
		argc-=2;
		argv[0]+=2;
		return true;
	};

	bool parse_bg(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No batch size provided for -bg (syntax is -bg batch_size)\n";
			return false;
		}
		if ( (meta.batchSizeTargetGPU = atoi((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for GPU batch size (syntax is -bg batch_size)\n";
			return false;
		}
		meta.batchSizeGPU = meta.batchSizeTargetGPU;
		argc-=2;
		argv[0]+=2;
		return true;
	};

    bool parse_c(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 4){
            cout << "Insufficient cell dimensions provided (syntax is -c x y z)\n";
            return false;
        }

        // the indexing in PRISM stores the cell dimensions as Z, Y, X so we must rearrange the
        // order of the inputs which are X, Y, Z
        if ( (meta.cellDim[2] = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for X cell dimension (syntax is -c x, y, z)\n";
            return false;
        }
        if ( (meta.cellDim[1] = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[2])) == 0){
            cout << "Invalid value \"" << (*argv)[2] << "\" provided for Y cell dimension (syntax is -c x, y, z)\n";
            return false;
        }
        if ( (meta.cellDim[0] = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[3])) == 0){
            cout << "Invalid value \"" << (*argv)[3] << "\" provided for Z cell dimension (syntax is -c x, y, z)\n";
            return false;
        }
        meta.userSpecifiedCelldims = true;
        argc-=4;
        argv[0]+=4;
        return true;
    };

    bool parse_C(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No state provided for -C (syntax is -f 0/1)\n";
            return false;
        }
        meta.alsoDoCPUWork = std::string((*argv)[1]) == "0" ? false : true;
        argc-=2;
        argv[0]+=2;
        return true;
    };

	bool parse_d(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No detector angle step provided for -d (syntax is -d detector_step (in mrad))\n";
			return false;
		}
		if ( (meta.detectorAngleStep = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1]) / 1000) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for potential bound (syntax is -d detector_step (in mrad)\n";
			return false;
		}
		argc-=2;
		argv[0]+=2;
		return true;
	};

    bool parse_streaming_mode(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No state provided for -C (syntax is -f 0/1)\n";
            return false;
        }
        meta.transferMode = std::string((*argv)[1]) == "0" ? Prismatic::StreamingMode::SingleXfer :  Prismatic::StreamingMode::Stream;
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_h(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        printHelp();
        return false;
    };

    bool parse_i(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No filename provided for -i (syntax is -i filename)\n";
            return false;
        }
        meta.filenameAtoms = std::string((*argv)[1]);
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_f(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No interpolation factor provided for -f (syntax is -f interpolation_factor)\n";
            return false;
        }
        if ( (meta.interpolationFactorX = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for PRISM interpolation factors (syntax is -f interpolation_factor)\n";
            return false;
        }
	    meta.interpolationFactorY = meta.interpolationFactorX;
        argc-=2;
        argv[0]+=2;
        return true;
    };

	bool parse_fx(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No interpolation factor provided for -fx (syntax is -fx interpolation_factor_x)\n";
			return false;
		}
		if ( (meta.interpolationFactorX = atoi((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for PRISM interpolation factor (syntax is -fx interpolation_factor_x)\n";
			return false;
		}
		argc-=2;
		argv[0]+=2;
		return true;
	};

	bool parse_fy(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No interpolation factor provided for -fy (syntax is -fy interpolation_factor_y)\n";
			return false;
		}
		if ( (meta.interpolationFactorY = atoi((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for PRISM interpolation factor (syntax is -fy interpolation_factor_y)\n";
			return false;
		}
		argc-=2;
		argv[0]+=2;
		return true;
	};

    bool parse_j(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No number of threads provided (syntax is -j numThreads)\n";
            return false;
        }
        if ( (meta.numThreads = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of threads  (syntax is -j numThreads)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_E(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No energy provided for -E (syntax is -E energy (in keV))\n";
            return false;
        }
        if ( (meta.E0 = 1000 * (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for energy  (syntax is -E energy (in keV))\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_F(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No number of frozen phonon configurations provided for -F (syntax is -F #)\n";
            return false;
        }
        if ( (meta.numFP = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of frozen phonon configurations (syntax is -F #)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_g(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No number of GPUs provided (syntax is -g numGPUs)\n";
            return false;
        }
        if ( (string((*argv)[1]) == "0") ){
            meta.numGPUs = 0;
            argc-=2;
            argv[0]+=2;
            return true;
        }
        if ( (meta.numGPUs = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of GPUs (syntax is -g numGPUs)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_s(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No slice thickness provided (syntax is -s slice_thickness (in Angstroms))\n";
            return false;
        }
        if ( (meta.sliceThickness = atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for slice_thickness (syntax is -s slice_thickness (in Angstroms))\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_S(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No number of CUDA streams per GPU provided (syntax is -S num_streams)\n";
            return false;
        }
        if ( (meta.numStreamsPerGPU = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of streams (syntax is -S num_streams)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_t(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                   int& argc, const char*** argv){
        if (argc < 4){
            cout << "Insufficient arguments provided for unit cell tiling (syntax is --tile-uc x y z)\n";
            return false;
        }

        // the indexing in PRISM stores the cell dimensions as Z, Y, X so we must rearrange the
        // order of the inputs which are X, Y, Z
        if ( (meta.tileX = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for unit cell tiling in X (syntax is --tile-uc x y z)\n";
            return false;
        }
        if ( (meta.tileY = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[2])) == 0){
            cout << "Invalid value \"" << (*argv)[2] << "\" provided for unit cell tiling in Y (syntax is --tile-uc x y z)\n";
            return false;
        }
        if ( (meta.tileZ = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[3])) == 0){
            cout << "Invalid value \"" << (*argv)[3] << "\" provided for unit cell tiling in Z (syntax is --tile-uc x y z)\n";
            return false;
        }

        argc-=4;
        argv[0]+=4;
        return true;


    };

    bool parse_o(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No filename provided for -o (syntax is -o filename)\n";
            return false;
        }
        meta.filenameOutput = std::string((*argv)[1]);
        //cout <<"meta.filenameAtoms = " << meta.filenameAtoms << endl;
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_p(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                              int& argc, const char*** argv){
        if (argc < 2){
            cout << "No pixel size provided for -p (syntax is -p pixel_size)\n";
            return false;
        }
        if ( (meta.realspacePixelSize[0] = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for pixel size  (syntax is -p pixel_size)\n";
            return false;
        }
	    meta.realspacePixelSize[1] = meta.realspacePixelSize[0];
        argc-=2;
        argv[0]+=2;
        return true;
    };

	bool parse_px(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No pixel size provided for -px (syntax is -px pixel_size)\n";
			return false;
		}
		if ( (meta.realspacePixelSize[1] = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for X pixel size  (syntax is -px pixel_size)\n";
			return false;
		}
		argc-=2;
		argv[0]+=2;
		return true;
	};

	bool parse_py(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No pixel size provided for -py (syntax is -py pixel_size)\n";
			return false;
		}
		if ( (meta.realspacePixelSize[0] = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for Y pixel size  (syntax is -py pixel_size)\n";
			return false;
		}
		argc-=2;
		argv[0]+=2;
		return true;
	};

    bool parse_P(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No bounding potential radius provided for -P (syntax is -P potential_bound (in Angstroms))\n";
            return false;
        }
        if ( (meta.potBound = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for potential bound (syntax is -P potential_bound (in Angstroms))\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_r(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                       int& argc, const char*** argv){
        if (argc < 2){
            cout << "No probe step provided for -r (syntax is -r probe_step (in Angstroms))\n";
            return false;
        }
        if ( (meta.probeStepX = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for probe_step (syntax is -r probe_step (in Angstroms))\n";
            return false;
        }
	    meta.probeStepY = meta.probeStepX;
        argc-=2;
        argv[0]+=2;
        return true;
    };


	bool parse_rx(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No probe step provided for -rx (syntax is -rx probe_step (in Angstroms))\n";
			return false;
		}
		if ( (meta.probeStepX = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for probe_step (syntax is -rx probe_step (in Angstroms))\n";
			return false;
		}
		meta.probeStepY = meta.probeStepX;
		argc-=2;
		argv[0]+=2;
		return true;
	};



	bool parse_ry(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No probe step provided for -ry (syntax is -ry probe_step (in Angstroms))\n";
			return false;
		}
		if ( (meta.probeStepY = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for probe_step (syntax is -ry probe_step (in Angstroms))\n";
			return false;
		}
		argc-=2;
		argv[0]+=2;
		return true;
	};


	bool parse_rs(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 2){
			cout << "No random seed provided for -rs (syntax is -rs integer)\n";
			return false;
		}
		if ( ((meta.randomSeed = atoi((*argv)[1])) == 0) & std::string(((*argv)[1]))!="0"){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for random seed (syntax is -rs integer)\n";
			return false;
		}
		argc-=2;
		argv[0]+=2;
		return true;
	};

	bool parse_tx(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                 int& argc, const char*** argv){
        if (argc < 2){
            cout << "No probe tilt provided for -tx (syntax is -tx probe_tilt)\n";
            return false;
        }
        if ( (meta.probeXtilt = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1]) / 1000) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for -tx (syntax is -tx probe_tilt\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_ty(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                  int& argc, const char*** argv){
        if (argc < 2){
            cout << "No probe tilt provided for -ty (syntax is -ty probe_tilt)\n";
            return false;
        }
        if ( (meta.probeYtilt = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1]) / 1000) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for -ty (syntax is -ty probe_tilt\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_df(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                  int& argc, const char*** argv){
        if (argc < 2){
            cout << "No defocus value provided for -df (syntax is -df defocus_value (in Angstroms))\n";
            return false;
        }
        if ( (meta.probeDefocus = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for -df (syntax is -df defocus_value (in Angstroms)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

	bool parse_C3(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	              int& argc, const char*** argv){
		if (argc < 2){
			cout << "No C3 value provided for -C3 (syntax is -C3 value (in Angstroms))\n";
			return false;
		}
		if ( (meta.C3 = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for -C3 (syntax is -C3 value (in Angstroms)\n";
			return false;
		}
		argc-=2;
		argv[0]+=2;
		return true;
	};

	bool parse_C5(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	              int& argc, const char*** argv){
		if (argc < 2){
			cout << "No C5 value provided for -C5 (syntax is -C5 value (in Angstroms))\n";
			return false;
		}
		if ( (meta.C5 = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1])) == 0){
			cout << "Invalid value \"" << (*argv)[1] << "\" provided for -C5 (syntax is -C5 value (in Angstroms)\n";
			return false;
		}
		argc-=2;
		argv[0]+=2;
		return true;
	};


	bool parse_sa(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                  int& argc, const char*** argv){
        if (argc < 2){
            cout << "No probe semiangle provided for -sa (syntax is -sa probe_semiangle in mrads)\n";
            return false;
        }
        if ( (meta.probeSemiangle = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1]) / 1000) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for -sa (syntax is -sa probe_semiangle in mrads)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

	bool parse_wx(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	             int& argc, const char*** argv){
		if (argc < 3){
			cout << "Invalid window provided for -wx (syntax is -wx min max (in fractional coordinates))\n";
			return false;
		}
		PRISMATIC_FLOAT_PRECISION minval, maxval;
		minval = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1]);
        maxval = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[2]);

		if ( (minval == 0) & (std::string((*argv)[1]) != "0")){
			cout << "Invalid lower bound \"" << (*argv)[1] << "\" provided for scan window X (syntax is -wx min max (in fractional coordinates))\n";
			return false;
		}
		if ( (maxval == 0) & (std::string((*argv)[2]) != "0")){
			cout << "Invalid upper bound \"" << (*argv)[2] << "\" provided for scan window X (syntax is -wx min max (in fractional coordinates))\n";
			return false;
		}
        if (maxval < minval){
            cout << "The provided lower bound(" << minval << ") for the X scan is greater than the maximum(" << maxval <<")." << endl;
         return false;
        }
        meta.scanWindowXMin = minval;
        meta.scanWindowXMax = maxval;
		argc-=3;
		argv[0]+=3;
		return true;
	};

    bool parse_wy(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                  int& argc, const char*** argv){
        if (argc < 3){
            cout << "Invalid window provided for -wy (syntax is -wy min max (in fractional coordinates))\n";
            return false;
        }
        PRISMATIC_FLOAT_PRECISION minval, maxval;
        minval = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1]);
        maxval = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[2]);
        if ( (minval == 0) & (std::string((*argv)[1]) != "0")){
            cout << "Invalid lower bound \"" << (*argv)[1] << "\" provided for scan window y (syntax is -wx min max (in fractional coordinates))\n";
            return false;
        }
        if ( (maxval == 0) & (std::string((*argv)[2]) != "0")){
            cout << "Invalid upper bound \"" << (*argv)[2] << "\" provided for scan window y (syntax is -wy min max (in fractional coordinates))\n";
            return false;
        }
        if (maxval < minval){
            cout << "The provided lower bound(" << minval << ") for the X scan is greater than the maximum(" << maxval <<")." << endl;
            return false;
        }
        meta.scanWindowYMin = minval;
        meta.scanWindowYMax = maxval;
        argc-=3;
        argv[0]+=3;
        return true;
    };

	bool parse_te(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
	              int& argc, const char*** argv){
		if (argc < 2){
			cout << "No value provided for -te (syntax is -te bool)\n";
			return false;
		}
		meta.includeThermalEffects = std::string((*argv)[1]) == "0" ? false : true;
		argc-=2;
		argv[0]+=2;
		return true;
	};

    bool parse_oc(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                  int& argc, const char*** argv){
        if (argc < 2){
            cout << "No value provided for -oc (syntax is -oc bool)\n";
            return false;
        }
        meta.includeOccupancy = std::string((*argv)[1]) == "0" ? false : true;
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_2D(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                  int& argc, const char*** argv){
        if (argc < 3){
            cout << "Not enough arguments for -2D (syntax is -2D ang_min ang_max)\n";
            return false;
        }
        meta.save2DOutput = true;
        if ( (string((*argv)[1]) != "0") & ((meta.integrationAngleMin = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[1]) / 1000) == 0)){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for minimum integration angle (syntax is -2D ang_min ang_max (in mrad)\n";
            return false;
        }
        if ( (meta.integrationAngleMax = (PRISMATIC_FLOAT_PRECISION)atof((*argv)[2]) / 1000) == 0){
            cout << "Invalid value \"" << (*argv)[2] << "\" provided for maximum integration angle (syntax is -2D ang_min ang_max (in mrad))\n";
            return false;
        }
        argc-=3;
        argv[0]+=3;
        return true;
    };

    bool parse_3D(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                  int& argc, const char*** argv){
        if (argc < 2){
            cout << "No value provided for -3D (syntax is -3D bool)\n";
            return false;
        }
        meta.save3DOutput = std::string((*argv)[1]) == "0" ? false : true;
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_4D(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                 int& argc, const char*** argv){
        if (argc < 2){
            cout << "No value provided for -4D (syntax is -4D bool)\n";
            return false;
        }
        meta.save4DOutput = std::string((*argv)[1]) == "0" ? false : true;
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parseInputs(Metadata<PRISMATIC_FLOAT_PRECISION> &meta,
                     int &argc, const char ***argv) {
        if (argc==1)return true; // case of no inputs to parse
        --argc;++(argv[0]);
        do {
            if (argc==0) return true; // successfully parsed all inputs
        } while (parseInput(meta, argc, argv));
        return false;
    }

    // use a lookup table to map the option switches to the corresponding function that validates/handles the arguments
    using parseFunction = bool (*)(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                                          int& argc, const char*** argv);
    static std::map<std::string, parseFunction> parser{
            {"--input-file", parse_i}, {"-i", parse_i},
            {"--interp-factor", parse_f}, {"-f", parse_f},
            {"--interp-factor-x", parse_fx}, {"-fx", parse_fx},
            {"--interp-factor-y", parse_fy}, {"-fy", parse_fy},
            {"--output-file", parse_o}, {"-o", parse_o},
            {"--num-threads", parse_j}, {"-j", parse_j},
            {"--num-streams", parse_S}, {"-S", parse_S},
            {"--slice-thickness", parse_s}, {"-s", parse_s},
            {"--num-gpus", parse_g}, {"-g", parse_g},
            {"--batch-size", parse_b}, {"-b", parse_b},
            {"--batch-size-cpu", parse_bc}, {"-bc", parse_bc},
            {"--batch-size-gpu", parse_bg}, {"-bg", parse_bg},
            {"--help", parse_h}, {"-h", parse_h},
            {"--pixel-size", parse_p}, {"-p", parse_p},
            {"--pixel-size-x", parse_px}, {"-px", parse_px},
            {"--pixel-size-y", parse_py}, {"-py", parse_py},
            {"--detector-angle-step", parse_d}, {"-d", parse_d},
            {"--cell-dimension", parse_c}, {"-c", parse_c},
            {"--algorithm", parse_a}, {"-a", parse_a},
            {"--energy", parse_E}, {"-E", parse_E},
            {"--alpha-max", parse_A}, {"-A", parse_A},
            {"--potential-bound", parse_P}, {"-P", parse_P},
            {"--also-do-cpu-work", parse_C}, {"-C", parse_C},
            {"--streaming-mode", parse_streaming_mode},
            {"--probe-step", parse_r}, {"-r", parse_r},
            {"--probe-step-x", parse_rx}, {"-rx", parse_rx},
            {"--probe-step-y", parse_ry}, {"-ry", parse_ry},
            {"--random-seed", parse_rs}, {"-rs", parse_rs},
            {"--probe-xtilt", parse_tx}, {"-tx", parse_tx},
            {"--probe-ytilt", parse_ty}, {"-ty", parse_ty},
            {"--probe-defocus", parse_df}, {"-df", parse_df},
            {"-C3", parse_C3},
            {"-C5", parse_C5},
            {"--probe-semiangle", parse_sa}, {"-sa", parse_sa},
            {"--scan-window-y", parse_wy}, {"-wy", parse_wy},
            {"--scan-window-x", parse_wx}, {"-wx", parse_wx},
            {"--tile-uc", parse_t}, {"-t", parse_t},
            {"--num-FP", parse_F}, {"-F", parse_F},
            {"--thermal-effects", parse_te}, {"-te", parse_te},
            {"--occupancy", parse_oc}, {"-oc", parse_oc},
            {"--save-2D-output", parse_2D}, {"-2D", parse_2D},
            {"--save-3D-output", parse_3D}, {"-3D", parse_3D},
            {"--save-4D-output", parse_4D}, {"-4D", parse_4D}
    };
    bool parseInput(Metadata<PRISMATIC_FLOAT_PRECISION>& meta,
                           int& argc, const char*** argv){
        parseFunction f = parser[std::string((*argv)[0])];
        if (f != NULL){
            return f(meta, argc, argv);
        } else{
            cout << "Invalid option \"" << (*argv)[0] << "\" provided\n";
            return false;
        }
    }
}
