#include "params.h"
#include "configure.h"
#include "parseInput.h"

using namespace std;
int main(int argc, const char** argv) {
	PRISM::Metadata<PRISM_FLOAT_PRECISION> meta;
	if (!PRISM::parseInputs(meta, argc, &argv))return 1;
    cout << "Successful parsing" << endl;
    meta.toString();
	meta.toString();
	PRISM::configure(meta);
	PRISM::execute_plan(meta);
	return 0;
}
