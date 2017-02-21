#include "include/fparams.h"
#include "include/ArrayND.h"
#include <iostream>
#include <stdlib.h>
using namespace std;
int main(){
    const double pixelSize = 100.0/1000.0;
    const double potBound = 1.0;
    const double numFP = 8.0/8.0;
    const size_t sliceThickness = 2;
    const size_t interpolationFactor = 10;
    
    cout << "test" << endl;
    cout << "fparams[0]" <<  fparams[0] << endl;
    cout << "fparams[10]" <<  fparams[10] << endl;
    cout << "fparams[20]" <<  fparams[20] << endl;

    using vec_d = std::vector<double>;
    using Array2D = PRISM::ArrayND<2, vec_d>;



    Array2D u = PRISM::ones_ND<2, double>({118,1}) * 0.8;
    for (auto& i : u) std::cout << i << std::endl;
    return 0;
}