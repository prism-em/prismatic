//
// Created by Alan Pryor on 4/12/17.
//

#include "parseInput.h"
#include <iostream>
#include <map>
#include <string>
#include <stdlib.h>

namespace PRISM {
    using namespace std;


    void printSyntax() {
        std::cout << "Syntax is ./prism [options]" << std::endl;
    }
    void printHelp() {
        std::cout << "Options:\n";
        std::cout << "-f interpolation_factor : PRISM interpolation factor\n";
        std::cout << "-g num_GPUs : number of GPUs to use\n";
        std::cout << "-i filename : input file containing atomic coordinates and species. Should be a \
                  comma-separated text file with one row per atom of the form x, y, z, Z where Z is the atomic number.\n";
        std::cout << "-j num_threads : number of CPU threads to use\n";
        std::cout << "-o filename : output filename\n";
        std::cout << "-s num_streams : number of CUDA streams to create per GPU\n";
    }

    ParseResult parse_a(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No algorithm provided for -a (syntax is -a algorithm). Choices are (m)ultislice or (p)rism\n";
            return ParseResult::Failure;
        }
        std::string algo = std::string((*argv)[1]);
        if (algo == "m" | algo == "multislice"){
            meta.algorithm = PRISM::Algorithm::Multislice;
        } else if (algo == "p" | algo == "prism"){
            meta.algorithm = PRISM::Algorithm::PRISM;
        } else {
            cout << "Unrecognized algorithm \"" << (*argv)[1] << "\"\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_A(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No maximum probe angle provided for -A (syntax is -A angle (in mrad))\n";
            return ParseResult::Failure;
        }
        if ( (meta.alphaBeamMax =  ( (PRISM_FLOAT_PRECISION)atof((*argv)[1])) / 1000 ) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for maximum probe angle (syntax is -A angle (in mrad))\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_c(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 4){
            cout << "Insufficient cell dimensions provided (syntax is -c x y z)\n";
            return ParseResult::Failure;
        }

        // the indexing in PRISM stores the cell dimensions as Z, Y, X so we must rearrange the
        // order of the inputs which are X, Y, Z
        if ( (meta.cellDim[2] = (PRISM_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for X cell dimension (syntax is -c x, y, z)\n";
            return ParseResult::Failure;
        }
        if ( (meta.cellDim[1] = (PRISM_FLOAT_PRECISION)atof((*argv)[2])) == 0){
            cout << "Invalid value \"" << (*argv)[2] << "\" provided for Y cell dimension (syntax is -c x, y, z)\n";
            return ParseResult::Failure;
        }
        if ( (meta.cellDim[0] = (PRISM_FLOAT_PRECISION)atof((*argv)[3])) == 0){
            cout << "Invalid value \"" << (*argv)[3] << "\" provided for Z cell dimension (syntax is -c x, y, z)\n";
            return ParseResult::Failure;
        }
        argc-=4;
        argv[0]+=4;
        return ParseResult::Success;
    };

    ParseResult parse_C(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No state provided for -C (syntax is -f 0/1)\n";
            return ParseResult::Failure;
        }
        meta.also_do_CPU_work = std::string((*argv)[1]) == "0" ? false : true;
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_force_streaming_mode(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        cout << "--force_streaming_mode not implemented\n";
        return ParseResult::Failure;
//        if (argc < 2){
//            cout << "No state provided for -C (syntax is -f 0/1)\n";
//            return ParseResult::Failure;
//        }
//        meta.stream_data = std::string((*argv)[1]) == "0" ? false : true;
//        argc-=2;
//        argv[0]+=2;
//        return ParseResult::Success;
    };

    ParseResult parse_h(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        printHelp();
        return ParseResult::Failure;
    };

    ParseResult parse_i(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No filename provided for -i (syntax is -i filename)\n";
            return ParseResult::Failure;
        }
        meta.filename_atoms = std::string((*argv)[1]);
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_f(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No interpolation factor provided for -f (syntax is -f interpolation_factor)\n";
            return ParseResult::Failure;
        }
        if ( (meta.interpolationFactor = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for PRISM interpolation factor (syntax is -f interpolation_factor)\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_j(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No number of threads provided (syntax is -j num_threads)\n";
            return ParseResult::Failure;
        }
        if ( (meta.NUM_THREADS = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of threads  (syntax is -j num_threads)\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_E(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No energy provided for -E (syntax is -E energy (in keV))\n";
            return ParseResult::Failure;
        }
        if ( (meta.E0 = 1000 * (PRISM_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for energy  (syntax is -E energy (in keV))\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_F(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No number of frozen phonon configurations provided for -F (syntax is -F #)\n";
            return ParseResult::Failure;
        }
        if ( (meta.numFP = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of frozen phonon configurations (syntax is -F #)\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_g(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No number of GPUs provided (syntax is -g num_GPUs)\n";
            return ParseResult::Failure;
        }
        if ( (meta.NUM_GPUS = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of GPUs (syntax is -g num_GPUs)\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_s(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No slice thickness provided (syntax is -s slice_thickness (in Angstroms))\n";
            return ParseResult::Failure;
        }
        if ( (meta.sliceThickness = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for slice_thickness (syntax is -s slice_thickness (in Angstroms))\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_S(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No number of CUDA streams per GPU provided (syntax is -S num_streams)\n";
            return ParseResult::Failure;
        }
        if ( (meta.NUM_STREAMS_PER_GPU = atoi((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for number of streams (syntax is -S num_streams)\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_o(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){

        if (argc < 2){
            cout << "No filename provided for -o (syntax is -o filename)\n";
            return ParseResult::Failure;
        }
        meta.filename_output = std::string((*argv)[1]);
        //cout <<"meta.filename_atoms = " << meta.filename_atoms << endl;
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_p(Metadata<PRISM_FLOAT_PRECISION>& meta,
                              int& argc, const char*** argv){
        if (argc < 2){
            cout << "No pixel size provided for -p (syntax is -p pixel_size)\n";
            return ParseResult::Failure;
        }
        if ( (meta.realspace_pixelSize = (PRISM_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for pixel size  (syntax is -p pixel_size)\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_P(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No bounding potential radius provided for -P (syntax is -P potential_bound (in Angstroms))\n";
            return ParseResult::Failure;
        }
        if ( (meta.potBound = (PRISM_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for potential bound (syntax is -P potential_bound (in Angstroms))\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    ParseResult parse_r(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
        if (argc < 2){
            cout << "No probe step provided for -r (syntax is -r probe_step (in Angstroms))\n";
            return ParseResult::Failure;
        }
        if ( (meta.probe_step = (PRISM_FLOAT_PRECISION)atof((*argv)[1])) == 0){
            cout << "Invalid value \"" << (*argv)[1] << "\" provided for probe_step (syntax is -r probe_step (in Angstroms))\n";
            return ParseResult::Failure;
        }
        argc-=2;
        argv[0]+=2;
        return ParseResult::Success;
    };

    bool parseInputs(Metadata<PRISM_FLOAT_PRECISION> &meta,
                     int &argc, const char ***argv) {
        if (argc==1)return true; // case of no inputs to parse
        --argc;++(argv[0]);
        ParseResult result;
        do {
            if (argc==0) return true; // successfully parsed all inputs
            result = parseInput(meta, argc, argv);
        } while (result == ParseResult::Success);
        if (result == ParseResult::Help)printHelp();
        return false;
    }

    // use a lookup table to map the option switches to the corresponding function that validates/handles the arguments
    using parseFunction = ParseResult (*)(Metadata<PRISM_FLOAT_PRECISION>& meta,
                                          int& argc, const char*** argv);
    static std::map<std::string, parseFunction> parser{
            {"--input_file", parse_i}, {"-i", parse_i},
            {"--interp_factor", parse_f}, {"-f", parse_f},
            {"--output_file", parse_o}, {"-o", parse_o},
            {"--num_threads", parse_j}, {"-j", parse_j},
            {"--num_streams", parse_S}, {"-S", parse_S},
            {"--slice_thickness", parse_s}, {"-s", parse_s},
            {"--num_gpus", parse_g}, {"-g", parse_g},
            {"--help", parse_h}, {"-h", parse_h},
            {"--pixel_size", parse_p}, {"-p", parse_p},
            {"--cell_dimension", parse_c}, {"-c", parse_c},
            {"--algorithm", parse_a}, {"-a", parse_a},
            {"--energy", parse_E}, {"-E", parse_E},
            {"--alpha_max", parse_A}, {"-A", parse_A},
            {"--potential_bound", parse_P}, {"-P", parse_P},
            {"--also_do_cpu_work", parse_C}, {"-C", parse_C},
            {"--force_streaming_mode", parse_force_streaming_mode},
            {"--probe_step", parse_r}, {"-r", parse_r},
            {"--num-FP", parse_F}, {"-F", parse_F}
    };
    ParseResult parseInput(Metadata<PRISM_FLOAT_PRECISION>& meta,
                           int& argc, const char*** argv){
        parseFunction f = parser[std::string((*argv)[0])];
        if (f != NULL){
            return f(meta, argc, argv);
        } else{
            cout << "Invalid option \"" << (*argv)[0] << "\" provided\n";
            return ParseResult::Failure;
        }
    }
}
