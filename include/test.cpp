//
// Created by AJ Pryor on 2/21/17.
//
#include "ArrayND.h"
#include <iostream>
#include <vector>

using namespace std;
int main(){
    //PRISM::ArrayND<3, std::vector<double> > arr(std::vector<double>(12, 0), std::vector<size_t>{2,3,2,4});
//    PRISM::ArrayND<3, std::vector<double> > arr(std::vector<double>(12, 0), std::array<size_t, 3>{2,3,2,4});
    //PRISM::ArrayND<3, std::vector<double> > arr(std::vector<double>(12, 0), {2,3,2});
    using vec_f = std::vector<double>;
    std::array<size_t, 4> dims{2,3,2,1};
    dims[3] = 2;
    PRISM::ArrayND<4, vec_f > arr4(vec_f(24, 0), dims);
    PRISM::ArrayND<3, vec_f > arr3(vec_f(12, 0), {2,3,2});
    PRISM::ArrayND<2, vec_f > arr2(vec_f(6, 0), {2,3});
    PRISM::ArrayND<1, vec_f > arr1(vec_f(2, 0), {2});
    for (auto &i:arr2)cout << i << '\n';
    cout << "arr.nrows = " << arr3.get_nrows() << endl;
    cout << "arr.ncols = " << arr3.get_ncols() << endl;
    cout << "arr.nlayers = " << arr3.get_nlayers() << endl;
   // cout << "arr.ndim4 = " << arr3.get_ndim4() << endl;
    return 0;
}
