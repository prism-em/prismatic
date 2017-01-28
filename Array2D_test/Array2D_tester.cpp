//
// Created by aj on 1/27/17.
//

#include <iostream>
#include "Array2D.h"

using namespace std;
int main(){
    vector<int> test {1,2,3,4,5,6};
    PRISM::Array2D<int> arr(test,2,3);
    cout << "arr nrows = " << arr.get_nrows() << endl;
    cout << "arr ncols = " << arr.get_ncols() << endl;
}