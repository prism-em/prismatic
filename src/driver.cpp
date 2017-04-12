#include "params.h"
#include "configure.h"
#include "parseInput.h"

using namespace std;
int main(int argc, const char** argv) {
	PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;
	/*
	if (argc>1) {
		meta.filename_atoms = std::string(argv[1]);
	} else{
		cout << "PRISM: Correct syntax is prism filename [interpolation factor] [algorithm] [probe spacing]  [pixel size] [slice thickness] [num threads] [num gpus]" << endl;
		return 0;
	}
	*/
	if (!PRISM::parseInputs(meta, argc, &argv))return 1;
    cout << "Successful parsing" << endl;
    meta.toString();
    return 0;
//	std::string filename;
//	meta.filename_output = "stack.mrc";
//	meta.toString();
//	PRISM::configure(meta);
//	PRISM::execute_plan(meta);
//	return 0;
}
