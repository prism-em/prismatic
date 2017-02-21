//
// Created by AJ Pryor on 1/27/17.
//

//Test construction of a PRISM:Array2D with std::vector
#include <iostream>
#include "ArrayND.h"
//#include <array>
using namespace std;
int main(){
    cout << "Creating a host-side 2x3 PRISM::Array2D with std::vector" << endl;
    vector<double> test{1,2,3,4,5,6};
    PRISM::ArrayND<2, std::vector<double> > arr(test,{2,3});
    cout << "nrows = " << arr.get_nrows() << endl;
    cout << "ncols = " << arr.get_ncols() << endl;
    cout << "Test Successful!" << endl;
}