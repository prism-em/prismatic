//
// Created by AJ Pryor on 1/27/17.
//

#include <iostream>
#include "Array2D.h"

using namespace std;
int main(){
    cout << "Creating a host-side PRISM::Array2D with std::vector" << endl;
    vector<int> test{1,2,3,4,5,6};
    PRISM::Array2D< std::vector<int> > arr(test,2,3);
    cout << "nrows = " << arr.get_nrows() << endl;
    cout << "ncols = " << arr.get_ncols() << endl;
}