//
// Created by aj on 1/27/17.
//

#ifndef PRISM_ARRAY2D_H
#define PRISM_ARRAY2D_H
#include <vector>
#include <cstddef>

namespace PRISM {

    template <class T>
    class Array2D {
    public:
        Array2D(std::vector<T> _data, size_t _nrows, size_t _ncols):data(_data), nrows(_nrows), ncols(_ncols){};
        size_t get_nrows() const {return this->nrows;}
        size_t get_ncols() const {return this->ncols;}
    private:
        std::vector<T> data;
        size_t nrows;
        size_t ncols;
    };
}


#endif //PRISM_ARRAY2D_H
