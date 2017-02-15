//
// Created by AJ Pryor on 2/10/17.
//

#ifndef PRISM_ARRAY3D_H
#define PRISM_ARRAY3D_H
#include <vector>
#include <cstddef>

namespace PRISM {

    template <class T>
    class Array3D {
        // 3D array class for nrows x ncolums x nlayers data indexed as C-style
        // T is expected to be a std::vector
    public:
        Array3D(T _data,
                const size_t& _nrows,
                const size_t& _ncols,
                const size_t& _nlayers);
        Array3D(){};
        size_t get_nrows()   const {return this->nrows;}
        size_t get_ncols()   const {return this->ncols;}
        size_t get_nlayers() const {return this->nlayers;}
        size_t size()        const {return this->N;}
        typename T::iterator begin();
        typename T::iterator end();
        typename T::value_type& at(const size_t& i, const size_t& j,const size_t& k);
        typename T::value_type& operator[](const size_t& i);

    private:
        T data;
        size_t nrows;
        size_t ncols;
        size_t nlayers;
        size_t N;
    };

    template <class T>
    Array3D<T>::Array3D(T _data,
                        const size_t& _nrows,
                        const size_t& _ncols,
                        const size_t& _nlayers):
            data(_data), nrows(_nrows), ncols(_ncols), nlayers(_nlayers){
        if (_data.size() != (_nrows * _ncols * _nlayers)) throw "PRISM: Size mismatch. Array size does not equal nrows * ncols * nlayers";
        this->N = _nrows * _ncols * _nlayers;
    };

    template <class T>
    typename T::iterator Array3D<T>::begin(){return this->data.begin();}

    template <class T>
    typename T::iterator Array3D<T>::end(){return this->data.end();}

    template <class T>
    typename T::value_type& Array3D<T>::at(const size_t& i, const size_t& j,const size_t& k){
        return data[i*ncols*nlayers + j*nlayers + k];
    }

    template <class T>
    typename T::value_type& Array3D<T>::operator[](const size_t& i){return data[i];}


    template <typename T>
    PRISM::Array3D< std::vector<T> > ones_3D(const size_t& nrows, const size_t& ncols, const size_t& nlayers){
        return PRISM::Array3D< std::vector<T> >(std::vector<T>(nrows*ncols,1),nrows, ncols, nlayers);
    }

    template <typename T>
    PRISM::Array3D< std::vector<T> > zeros_3D(const size_t& nrows, const size_t& ncols, const size_t& nlayers){
        return PRISM::Array3D< std::vector<T> >(std::vector<T>(nrows*ncols,0),nrows, ncols, nlayers);
    }
}


#endif //PRISM_ARRAY3D_H
