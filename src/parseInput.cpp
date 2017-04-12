//
// Created by Alan Pryor on 4/12/17.
//

#include "parseInput.h"
#include <iostream>
#include <map>
#include <string>


namespace PRISM {
    using namespace std;
    using parseFunction = ParseResult (*)(Metadata<PRISM_FLOAT_PRECISION>& meta,
                                          int& argc, const char*** argv);


    void printSyntax() {
        std::cout << "test" << std::endl;
    }
    void printHelp() {
        std::cout << "help" << std::endl;
    }
    ParseResult parse_f(Metadata<PRISM_FLOAT_PRECISION>& meta,
                        int& argc, const char*** argv){
       // cout << "f parser" << endl;

        if (argc < 2){
            cout << "Invalid filename provided with -f (syntax is -f filename)\n";
            return ParseResult::Failure;
        }
        meta.filename_atoms = std::string((*argv)[1]);
        //cout <<"meta.filename_atoms = " << meta.filename_atoms << endl;
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
    static std::map<std::string, parseFunction> parser{{"-f", parse_f}};
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
