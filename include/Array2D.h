//
// Created by AJ Pryor on 1/27/17
//

#ifndef PRISM_ARRAY2D_H
#define PRISM_ARRAY2D_H
#include <vector>
#include <cstddef>

namespace PRISM {

    template <class T>
    class Array2D {
    public:
        Array2D(T _data,
                const size_t& _nrows,
                const size_t& _ncols);
        Array2D(){};
        Array2D(const Array2D<T>& other);
        size_t get_nrows() const {return this->nrows;}
        size_t get_ncols() const {return this->ncols;}
        size_t size()      const {return this->N;}
        typename T::iterator begin();
        typename T::iterator end();
        typename T::iterator begin() const;
        typename T::iterator end() const;
        typename T::value_type& at(const size_t& i, const size_t& j);
        typename T::value_type& operator[](const size_t& i);

        Array2D<T> operator-(const Array2D<T>& other);
        Array2D<T> operator+(const Array2D<T>& other);
        Array2D<T> operator*(const Array2D<T>& other) ;
        Array2D<T> operator/(const Array2D<T>& other);
        Array2D<T> operator-(const typename T::value_type& val);
        Array2D<T> operator+(const typename T::value_type& val);
        Array2D<T> operator*(const typename T::value_type& val);
        Array2D<T> operator/(const typename T::value_type& val);

    private:
        T data;
        size_t nrows;
        size_t ncols;
        size_t N;
    };

    template <class T>
    Array2D<T>::Array2D(T _data,
                        const size_t& _nrows,
                        const size_t& _ncols):data(_data), nrows(_nrows), ncols(_ncols){
        if (_data.size() != (_nrows * _ncols)) throw "PRISM: Size mismatch. Array size does not equal nrows * ncols";
        this->N = _nrows * _ncols;
    };

    template <class T>
    Array2D<T>::Array2D(const Array2D<T>& other)
            : data(other.data), nrows(other.nrows), ncols(other.ncols),N(other.N){};

    template <class T>
    typename T::iterator Array2D<T>::begin(){return this->data.begin();}

    template <class T>
    typename T::iterator Array2D<T>::end(){return this->data.end();}

    template <class T>
    typename T::iterator Array2D<T>::begin()const{return this->data.begin();}

    template <class T>
    typename T::iterator Array2D<T>::end()const{return this->data.end();}

    template <class T>
    typename T::value_type& Array2D<T>::at(const size_t& i, const size_t& j){
        return this->data[i*ncols + j];
    }

    template <class T>
    typename T::value_type& Array2D<T>::operator[](const size_t& i){return data[i];}

    template <class T>
    Array2D<T> Array2D<T>::operator-(const Array2D<T>& other){
        Array2D<T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i-=*o++;
    }

    template <class T>
    Array2D<T> Array2D<T>::operator+(const Array2D<T>& other){
        Array2D<T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i+=*o++;
    }

    template <class T>
    Array2D<T> Array2D<T>::operator*(const Array2D<T>& other){
        Array2D<T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i*=*o++;
    }

    template <class T>
    Array2D<T> Array2D<T>::operator/(const Array2D<T>& other){
        Array2D<T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i/=*o++;
    }

    template <class T>
    Array2D<T> Array2D<T>::operator-(const typename T::value_type& val){
        Array2D<T> result(*this);
        for (auto& i:result)i-=val;
        return result;
    }

    template <class T>
    Array2D<T> Array2D<T>::operator+(const typename T::value_type& val){
        Array2D<T> result(*this);
        for (auto& i:result)i+=val;
        return result;
    }

    template <class T>
    Array2D<T> Array2D<T>::operator*(const typename T::value_type& val){
        Array2D<T> result(*this);
        for (auto& i:result)i*=val;
        return result;
    }

    template <class T>
    Array2D<T> Array2D<T>::operator/(const typename T::value_type& val){
        Array2D<T> result(*this);
        for (auto& i:result)i/=val;
        return result;
    }
}


#endif //PRISM_ARRAY2D_H
