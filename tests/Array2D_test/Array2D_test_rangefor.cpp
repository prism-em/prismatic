//
// Created by AJ Pryor on 1/30/17.
//
//Test range-based for loops with PRISM::Array2D

#include <iostream>
#include "Array2D.h"

using namespace std;
int main(){
    cout << "Creating a host-side 2x3 PRISM::Array2D with std::vector" << endl;
    vector<int> test{1,2,3,4,5,6};
    PRISM::Array2D< std::vector<int> > arr(test,2,3);
    cout << "Printing elements of data held in PRISM::Array2D" << endl;
    for (auto i:arr) std::cout << i << std::endl;
    cout << "Test Successful!" << endl;

}