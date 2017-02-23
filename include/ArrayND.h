//
// Created by AJ Pryor on 2/21/17.
//

#ifndef PRISM_ARRAYND_H
#define PRISM_ARRAYND_H
#include <vector>
#include <cstddef>
#include <iostream>
#include <array>
#include <utility>
namespace PRISM {
    template <size_t N, class T>
    class ArrayND {
            // ND array class for data indexed as C-style
            // T is expected to be a std::vector
        public:
            ArrayND(T _data,
                    std::array<size_t, N> _dims);
            ArrayND(){};
            size_t get_nrows()   const {return this->dims[0];}
            size_t get_ncols()   const {return this->dims[1];}
            size_t get_nlayers() const {return this->dims[2];}
            size_t get_ndim4()   const {return this->dims[3];}
            size_t size()        const {return this->arr_size;}
            typename T::iterator begin();
            typename T::iterator end();
            typename T::iterator begin() const;
            typename T::iterator end()   const;
            typename T::value_type& at(const size_t& i);
            typename T::value_type& at(const size_t& i, const size_t& j);
            typename T::value_type& at(const size_t& i, const size_t& j,const size_t& k);
            typename T::value_type& at(const size_t& i, const size_t& j,const size_t& k, const size_t& l);
	        typename T::value_type at(const size_t& i)const;
	        typename T::value_type at(const size_t& i, const size_t& j)const;
	        typename T::value_type at(const size_t& i, const size_t& j,const size_t& k)const;
	        typename T::value_type at(const size_t& i, const size_t& j,const size_t& k, const size_t& l)const;

	        typename T::value_type& operator[](const size_t& i);
            typename T::value_type operator[](const size_t& i)const;
            ArrayND<N, T> operator-(const ArrayND<N, T>& other);
            ArrayND<N, T> operator+(const ArrayND<N, T>& other);
            ArrayND<N, T> operator*(const ArrayND<N, T>& other) ;
            ArrayND<N, T> operator/(const ArrayND<N, T>& other);
            ArrayND<N, T> operator-(const typename T::value_type& val);
            ArrayND<N, T> operator+(const typename T::value_type& val);
            ArrayND<N, T> operator*(const typename T::value_type& val);
            ArrayND<N, T> operator/(const typename T::value_type& val);
            ArrayND<N, T>& operator-=(const typename T::value_type& val);
            ArrayND<N, T>& operator+=(const typename T::value_type& val);
            ArrayND<N, T>& operator*=(const typename T::value_type& val);
            ArrayND<N, T>& operator/=(const typename T::value_type& val);

        private:
            std::array<size_t, N> dims;
            std::array<size_t, N-1> strides;
            size_t arr_size;
            T data;
            
        };

        template <size_t N, class T>
        ArrayND<N, T>::ArrayND(T _data,
                               std::array<size_t, N> _dims):data(_data){
            //static_assert(N==_dims.size(),"exit");
           // if (_dims.size() != N) {
             //   std::cout << "PRISM: Dimension mismatch! Provided number of dimensions does not equal N for ND array\n";
               // throw "PRISM: Dimension mismatch! Provided number of dimensions does not equal N for ND array";
            //}
            size_t _size = 1;
            for (auto& i:_dims)_size*=i;
            if (_data.size() != _size){
                throw std::domain_error("PRISM: Size mismatch! Desired array size does not match size of input data\n");
            }
            this->arr_size = _size;
            this->dims     = _dims;

            size_t stride = 1;
            for (auto i = (N-1); i > 0; --i) {
                stride *= this->dims[i];
                this->strides[i-1] = stride;
            }
        };

    template <size_t N, class T>
    typename T::iterator ArrayND<N, T>::begin(){return this->data.begin();}

    template <size_t N, class T>
    typename T::iterator ArrayND<N, T>::end(){return this->data.end();}

    template <size_t N, class T>
    typename T::iterator ArrayND<N, T>::begin() const {return this->data.begin();}

    template <size_t N, class T>
    typename T::iterator ArrayND<N, T>::end() const {return this->data.end();}

    template <size_t N, class T>
    typename T::value_type& ArrayND<N, T>::at(const size_t& i){
        return data[i];
    }

    template <size_t N, class T>
    typename T::value_type& ArrayND<N, T>::at(const size_t& j, const size_t& i){
        return data[j*strides[0] + i];
    }

    template <size_t N, class T>
    typename T::value_type& ArrayND<N, T>::at(const size_t& k, const size_t& j,const size_t& i){
        return data[k*strides[0] + j*strides[1] + i];
    }

    template <size_t N, class T>
    typename T::value_type& ArrayND<N, T>::at(const size_t& l, const size_t& k,const size_t& j, const size_t& i){
        return data[l*strides[0] + k*strides[1] + j*strides[2] + i];
    }

	template <size_t N, class T>
	typename T::value_type ArrayND<N, T>::at(const size_t& i)const{
		return data[i];
	}

	template <size_t N, class T>
	typename T::value_type ArrayND<N, T>::at(const size_t& j, const size_t& i)const{
		return data[j*strides[0] + i];
	}

	template <size_t N, class T>
	typename T::value_type ArrayND<N, T>::at(const size_t& k, const size_t& j,const size_t& i)const{
		return data[k*strides[0] + j*strides[1] + i];
	}

	template <size_t N, class T>
	typename T::value_type ArrayND<N, T>::at(const size_t& l, const size_t& k,const size_t& j, const size_t& i)const{
		return data[l*strides[0] + k*strides[1] + j*strides[2] + i];
	}



    template <size_t N, class T>
    typename T::value_type& ArrayND<N, T>::operator[](const size_t& i){return data[i];}

    template <size_t N, class T>
    typename T::value_type ArrayND<N, T>::operator[](const size_t& i)const{return data[i];}

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator-(const ArrayND<N, T>& other){
        ArrayND<N, T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i-=*o++;
    }

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator+(const ArrayND<N, T>& other){
        ArrayND<N, T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i+=*o++;
    }

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator*(const ArrayND<N, T>& other){
        ArrayND<N, T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i*=*o++;
    }

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator/(const ArrayND<N, T>& other){
        ArrayND<N, T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i/=*o++;
    }

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator-(const typename T::value_type& val){
        ArrayND<N, T> result(*this);
        for (auto& i:result)i-=val;
        return result;
    }

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator+(const typename T::value_type& val){
        ArrayND<N, T> result(*this);
        for (auto& i:result)i+=val;
        return result;
    }

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator*(const typename T::value_type& val){
        ArrayND<N, T> result(*this);
        for (auto& i:result)i*=val;
        return result;
    }

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator/(const typename T::value_type& val){
        ArrayND<N, T> result(*this);
        for (auto& i:result)i/=val;
        return result;
    }


    template <size_t N, class T>
    ArrayND<N, T>& ArrayND<N, T>::operator-=(const typename T::value_type& val){
        for (auto& i:(*this))i-=val;
        return *this;
    }

	template <size_t N, class T>
	ArrayND<N, T>& ArrayND<N, T>::operator+=(const typename T::value_type& val){
		for (auto& i:(*this))i+=val;
		return *this;
	}

	template <size_t N, class T>
	ArrayND<N, T>& ArrayND<N, T>::operator*=(const typename T::value_type& val){
		for (auto& i:(*this))i*=val;
		return *this;
	}

	template <size_t N, class T>
	ArrayND<N, T>& ArrayND<N, T>::operator/=(const typename T::value_type& val){
		for (auto& i:(*this))i/=val;
		return *this;
	}


    template <>
    size_t ArrayND<3, std::vector<double> >::get_ndim4()   const = delete;

    template <size_t N, class T>
    PRISM::ArrayND<N, std::vector<T> > ones_ND(const std::array<size_t, N>& dims){
        size_t size = 1;
        for (auto& i:dims)size*=i;
        return PRISM::ArrayND<N, std::vector<T> >(std::vector<T>(size,1), dims);
    }

    template <size_t N, class T>
    PRISM::ArrayND<N, std::vector<T> > zeros_ND(const std::array<size_t, N>& dims){
        size_t size = 1;
        for (auto& i:dims)size*=i;
        return PRISM::ArrayND<N, std::vector<T> >(std::vector<T>(size,0), dims);
    }

	template <class T>
	using Array2D = PRISM::ArrayND<2, std::vector<T> >;
	template <class T>
	using Array1D = PRISM::ArrayND<1, std::vector<T> >;
	template <class T>
	std::pair<Array2D<T>, Array2D<T>> meshgrid(const Array1D<T>& X, const Array1D<T>& Y){
		Array2D<T> xx = zeros_ND<2, T>({Y.size(), X.size()});
		Array2D<T> yy = zeros_ND<2, T>({Y.size(), X.size()});
		for (auto i = 0; i < xx.get_nrows(); ++i){
			for (auto j = 0; j < xx.get_ncols(); ++j){
				xx.at(i,j) = X[j];
				yy.at(i,j) = Y[i];
			}
		}
		return std::pair<Array2D<T>, Array2D<T> >(xx,yy);
	}
}

#endif //PRISM_ARRAYND_H
