// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include "parseInput.h"
#include <iostream>
#include <map>
#include <string>
#include <stdlib.h>

namespace PRISM {
    using namespace std;

    void printHelp() {
        std::cout << "The following options are available with `prism`, each documented as long form (short form) *parameters* : description\n"
                "\n"
                "* --input-file (-i) filename : the filename containing the atomic coordinates, which should be a plain text file with comma-separated values in the format x, y, z, Z \n"
                "* --output-file(-o) filename : output filename\n"
                "* --interp-factor (-f) number : PRISM interpolation factor\n"
                "* --num-threads (-j) value : number of CPU threads to use\n"
                "* --num-streams (-S) value : number of CUDA streams to create per GPU\n"
                "* --slice-thickness (-s) thickness : thickness of each slice of projected potential (in Angstroms)\n"
                "* --num-gpus (-g) value : number of GPUs to use. A runtime check is performed to check how many are actually available, and the minimum of these two numbers is used.\n"
                "* --help(-h) : print information about the available options\n"
                "* --pixel-size (-p) pixel_size : size of simulation pixel size\n"
                "* --cell-dimension (-c) x y z : size of sample in x, y, z directions (in Angstroms)\n"
                "* --algorithm (-a) p/m : the simulation algorithm to use, either (p)rism or (m)ultislice\n"
                "* --energy (-E) value : the energy of the electron beam (in keV)\n"
                "* --alpha-max (-A) angle : the maximum probe angle to consider (in mrad)\n"
                "* --potential-bound (-P) value : the maximum radius from the center of each atom to compute the potental (in Angstroms)\n"
                "* --also-do-cpu-work (-C) 0/1 : boolean value used to determine whether or not to also create CPU workers in addition to GPU ones\n"
                "* --streaming-mode 0/1 : boolean value to force code to use (true) or not use (false) streaming versions of GPU codes. The default behavior is to estimate the needed memory from input parameters and choose automatically.\n"
                "* --probe-step (-r) step_size : step size of the probe (in Angstroms)\n"
                "* --num-FP (-F) value : number of frozen phonon configurations to calculate\n"
                "* --save-4D-output (-4D) bool : Also save the 2D output at the detector for each probe (4D output mode)\n";
    }

    bool parse_a(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No algorithm provided for -a (syntax is -a algorithm). Choices are (m)ultislice or (p)rism\n";
            return false;
        }
        std::string algo = std::string((*argv)[1]);
        if (algo == "m" | algo == "multislice"){
            meta.algorithm = PRISM::Algorithm::Multislice;
        } else if (algo == "p" | algo == "prism"){
            meta.algorithm = PRISM::Algorithm::PRISM;
        } else {
            cout << "Unrecognized algorithm \"" << (*argv)[1] << "\"\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_A(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No maximum probe angle provided for -A (syntax is -A angle (in mrad))\n";
            return false;
        }
        if ( (meta.alphaBeamMax =  ( (PRISM_FLOAT_PRECISION)atof((*argv)[1])) / 1000 ) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for maximum probe angle (syntax is -A angle (in mrad))\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_c(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 4){
            cout << "Insufficient cell dimensions provided (syntax is -c x y z)\n";
            return false;
        }

        // the indexing in PRISM stores the cell dimensions as Z, Y, X so we must rearrange the
        // order of the inputs which are X, Y, Z
        if ( (meta.cellDim[2] = (PRISM_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for X cell dimension (syntax is -c x, y, z)\n";
            return false;
        }
        if ( (meta.cellDim[1] = (PRISM_FLOAT_PRECISION)atof((*argv)[2])) == 0){
            cout << "Invalid value \"" << (*argv)[2] << "\" provided for Y cell dimension (syntax is -c x, y, z)\n";
            return false;
        }
        if ( (meta.cellDim[0] = (PRISM_FLOAT_PRECISION)atof((*argv)[3])) == 0){
            cout << "Invalid value \"" << (*argv)[3] << "\" provided for Z cell dimension (syntax is -c x, y, z)\n";
            return false;
        }
        argc-=4;
        argv[0]+=4;
        return true;
    };

    bool parse_C(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No state provided for -C (syntax is -f 0/1)\n";
            return false;
        }
        meta.also_do_CPU_work = std::string((*argv)[1]) == "0" ? false : true;
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_streaming_mode(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No state provided for -C (syntax is -f 0/1)\n";
            return false;
        }
        meta.transfer_mode = std::string((*argv)[1]) == "0" ? PRISM::StreamingMode::SingleXfer :  PRISM::StreamingMode::Stream;
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_h(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        printHelp();
        return false;
    };

    bool parse_i(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No filename provided for -i (syntax is -i filename)\n";
            return false;
        }
        meta.filename_atoms = std::string((*argv)[1]);
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_f(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No interpolation factor provided for -f (syntax is -f interpolation_factor)\n";
            return false;
        }
        if ( (meta.interpolationFactor = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for PRISM interpolation factor (syntax is -f interpolation_factor)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_j(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No number of threads provided (syntax is -j num_threads)\n";
            return false;
        }
        if ( (meta.NUM_THREADS = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of threads  (syntax is -j num_threads)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_E(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No energy provided for -E (syntax is -E energy (in keV))\n";
            return false;
        }
        if ( (meta.E0 = 1000 * (PRISM_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for energy  (syntax is -E energy (in keV))\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_F(Metadata<PRISM_FLOAT_PRECISION>& meta,
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

    bool parse_g(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No number of GPUs provided (syntax is -g num_GPUs)\n";
            return false;
        }
        if ( (string((*argv)[1]) == "0") ){
            meta.NUM_GPUS = 0;
            argc-=2;
            argv[0]+=2;
            return true;
        }
        if ( (meta.NUM_GPUS = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of GPUs (syntax is -g num_GPUs)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_s(Metadata<PRISM_FLOAT_PRECISION>& meta,
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

    bool parse_S(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No number of CUDA streams per GPU provided (syntax is -S num_streams)\n";
            return false;
        }
        if ( (meta.NUM_STREAMS_PER_GPU = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of streams (syntax is -S num_streams)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_o(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No filename provided for -o (syntax is -o filename)\n";
            return false;
        }
        meta.filename_output = std::string((*argv)[1]);
        //cout <<"meta.filename_atoms = " << meta.filename_atoms << endl;
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_p(Metadata<PRISM_FLOAT_PRECISION>& meta,
                              int& argc, const char*** argv){
        if (argc < 2){
            cout << "No pixel size provided for -p (syntax is -p pixel_size)\n";
            return false;
        }
        if ( (meta.realspace_pixelSize = (PRISM_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for pixel size  (syntax is -p pixel_size)\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_P(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No bounding potential radius provided for -P (syntax is -P potential_bound (in Angstroms))\n";
            return false;
        }
        if ( (meta.potBound = (PRISM_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for potential bound (syntax is -P potential_bound (in Angstroms))\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_r(Metadata<PRISM_FLOAT_PRECISION>& meta,
                       int& argc, const char*** argv){
        if (argc < 2){
            cout << "No probe step provided for -r (syntax is -r probe_step (in Angstroms))\n";
            return false;
        }
        if ( (meta.probe_step = (PRISM_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for probe_step (syntax is -r probe_step (in Angstroms))\n";
            return false;
        }
        argc-=2;
        argv[0]+=2;
        return true;
    };

    bool parse_4D(Metadata<PRISM_FLOAT_PRECISION>& meta,
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

    bool parseInputs(Metadata<PRISM_FLOAT_PRECISION> &meta,
                     int &argc, const char ***argv) {
        if (argc==1)return true; // case of no inputs to parse
        --argc;++(argv[0]);
        do {
            if (argc==0) return true; // successfully parsed all inputs
        } while (parseInput(meta, argc, argv));
        return false;
    }

    // use a lookup table to map the option switches to the corresponding function that validates/handles the arguments
    using parseFunction = bool (*)(Metadata<PRISM_FLOAT_PRECISION>& meta,
                                          int& argc, const char*** argv);
    static std::map<std::string, parseFunction> parser{
            {"--input-file", parse_i}, {"-i", parse_i},
            {"--interp-factor", parse_f}, {"-f", parse_f},
            {"--output-file", parse_o}, {"-o", parse_o},
            {"--num-threads", parse_j}, {"-j", parse_j},
            {"--num-streams", parse_S}, {"-S", parse_S},
            {"--slice-thickness", parse_s}, {"-s", parse_s},
            {"--num-gpus", parse_g}, {"-g", parse_g},
            {"--help", parse_h}, {"-h", parse_h},
            {"--pixel-size", parse_p}, {"-p", parse_p},
            {"--cell-dimension", parse_c}, {"-c", parse_c},
            {"--algorithm", parse_a}, {"-a", parse_a},
            {"--energy", parse_E}, {"-E", parse_E},
            {"--alpha-max", parse_A}, {"-A", parse_A},
            {"--potential-bound", parse_P}, {"-P", parse_P},
            {"--also-do-cpu-work", parse_C}, {"-C", parse_C},
            {"--streaming-mode", parse_streaming_mode},
            {"--probe-step", parse_r}, {"-r", parse_r},
            {"--num-FP", parse_F}, {"-F", parse_F},
            {"--save-4D-output", parse_4D}, {"-4D", parse_4D}
    };
    bool parseInput(Metadata<PRISM_FLOAT_PRECISION>& meta,
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
