// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// Prismatic is distributed under the GNU General Public License (GPL)
// If you use Prismatic, we kindly ask that you cite the following papers:

// 1. Ophus, C.: A fast image simulation algorithm for scanning
//    transmission electron microscopy. Advanced Structural and
//    Chemical Imaging 3(1), 13 (2017)

// 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
//    Implementation of Image Simulation Algorithms for Scanning
//	  Transmission Electron Microscopy. arXiv:1706.08563 (2017)

#ifndef PRISM_ARRAYND_H
#define PRISM_ARRAYND_H

#include <vector>
#include <cstddef>
#include <iostream>
#include <array>
#include <utility>
#include <fstream>
#include <cstring>
#include <complex>
namespace Prismatic {
    template <size_t N, class T>
    class ArrayND {
            // ND array class for data indexed as C-style, i.e. arr.at(k,j,i) where i is the fastest varying index
            // and k is the slowest

            // T is expected to be a std::vector
        public:
            ArrayND(T _data,
                    std::array<size_t, N> _dims);
            ArrayND(){};
			size_t get_dimi() const {return this->dims[N-1];}
			size_t get_dimj() const {return this->dims[N-2];}
			size_t get_dimk() const {return this->dims[N-3];}
			size_t get_diml() const {return this->dims[N-4]; }
			size_t get_dimm() const {return this->dims[N-5]; }
            size_t size()     const {return this->arr_size;}
            typename T::iterator begin();
            typename T::iterator end();
            typename T::const_iterator begin() const;
            typename T::const_iterator end()   const;
            typename T::value_type& at(const size_t& i);
            typename T::value_type& at(const size_t& j, const size_t& i);
            typename T::value_type& at(const size_t& k, const size_t& j,const size_t& i);
            typename T::value_type& at(const size_t& l, const size_t& k,const size_t& j, const size_t& i);
	        typename T::value_type at(const size_t& i)const;
	        typename T::value_type at(const size_t& j, const size_t& i)const;
	        typename T::value_type at(const size_t& k, const size_t& j,const size_t& i)const;
	        typename T::value_type at(const size_t& l, const size_t& k,const size_t& j, const size_t& i)const;

	        typename T::value_type& operator[](const size_t& i);
            typename T::value_type operator[](const size_t& i)const;
            ArrayND<N, T>  operator-( const ArrayND<N, T>& other);
            ArrayND<N, T>  operator+( const ArrayND<N, T>& other);
            ArrayND<N, T>  operator*( const ArrayND<N, T>& other);
            ArrayND<N, T>  operator/( const ArrayND<N, T>& other);
            ArrayND<N, T>  operator-( const typename T::value_type& val);
            ArrayND<N, T>  operator+( const typename T::value_type& val);
            ArrayND<N, T>  operator*( const typename T::value_type& val);
            ArrayND<N, T>  operator/( const typename T::value_type& val);
	        ArrayND<N, T>  operator-( const ArrayND<N, T>& other)const;
	        ArrayND<N, T>  operator+( const ArrayND<N, T>& other)const;
	        ArrayND<N, T>  operator*( const ArrayND<N, T>& other)const;
	        ArrayND<N, T>  operator/( const ArrayND<N, T>& other)const;
			ArrayND<N, T>& operator-=(const ArrayND<N, T>& other);
			ArrayND<N, T>& operator+=(const ArrayND<N, T>& other);
			ArrayND<N, T>& operator*=(const ArrayND<N, T>& other);
			ArrayND<N, T>& operator/=(const ArrayND<N, T>& other);
	        ArrayND<N, T>  operator-( const typename T::value_type& val)const;
	        ArrayND<N, T>  operator+( const typename T::value_type& val)const;
	        ArrayND<N, T>  operator*( const typename T::value_type& val)const;
	        ArrayND<N, T>  operator/( const typename T::value_type& val)const;
            ArrayND<N, T>& operator-=(const typename T::value_type& val);
            ArrayND<N, T>& operator+=(const typename T::value_type& val);
            ArrayND<N, T>& operator*=(const typename T::value_type& val);
            ArrayND<N, T>& operator/=(const typename T::value_type& val);

	    inline void toMRC_f(const char* filename)const;
        private:
            std::array<size_t, N> dims;
            std::array<size_t, N-1> strides;
            size_t arr_size;
            T data;
            
        };

        template <size_t N, class T>
        ArrayND<N, T>::ArrayND(T _data,
                               std::array<size_t, N> _dims):data(_data){
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
    typename T::const_iterator ArrayND<N, T>::begin() const {return this->data.begin();}

    template <size_t N, class T>
    typename T::const_iterator ArrayND<N, T>::end() const {return this->data.end();}

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
	    return result;
    }

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator+(const ArrayND<N, T>& other){
        ArrayND<N, T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i+=*o++;
	    return result;
    }

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator*(const ArrayND<N, T>& other){
        ArrayND<N, T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i*=*o++;
	    return result;
    }

    template <size_t N, class T>
    ArrayND<N, T> ArrayND<N, T>::operator/(const ArrayND<N, T>& other){
        ArrayND<N, T> result(*this);
        typename T::value_type* o = other.begin();
        for (auto& i:result)i/=*o++;
	    return result;
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
	ArrayND<N, T> ArrayND<N, T>::operator-(const ArrayND<N, T>& other)const{
		ArrayND<N, T> result(*this);
		typename T::value_type* o = other.begin();
		for (auto& i:result)i-=*o++;
		return result;
	}

	template <size_t N, class T>
	ArrayND<N, T> ArrayND<N, T>::operator+(const ArrayND<N, T>& other)const{
		ArrayND<N, T> result(*this);
		typename T::value_type* o = other.begin();
		for (auto& i:result)i+=*o++;
		return result;
	}

	template <size_t N, class T>
	ArrayND<N, T> ArrayND<N, T>::operator*(const ArrayND<N, T>& other)const{
		ArrayND<N, T> result(*this);
		typename T::value_type* o = other.begin();
		for (auto& i:result)i*=*o++;
		return result;
	}

	template <size_t N, class T>
	ArrayND<N, T> ArrayND<N, T>::operator/(const ArrayND<N, T>& other)const{
		ArrayND<N, T> result(*this);
		typename T::value_type* o = other.begin();
		for (auto& i:result)i/=*o++;
		return result;
	}

	template <size_t N, class T>
	ArrayND<N, T> ArrayND<N, T>::operator-(const typename T::value_type& val)const{
		ArrayND<N, T> result(*this);
		for (auto& i:result)i-=val;
		return result;
	}

	template <size_t N, class T>
	ArrayND<N, T> ArrayND<N, T>::operator+(const typename T::value_type& val)const{
		ArrayND<N, T> result(*this);
		for (auto& i:result)i+=val;
		return result;
	}

	template <size_t N, class T>
	ArrayND<N, T> ArrayND<N, T>::operator*(const typename T::value_type& val)const{
		ArrayND<N, T> result(*this);
		for (auto& i:result)i*=val;
		return result;
	}

	template <size_t N, class T>
	ArrayND<N, T> ArrayND<N, T>::operator/(const typename T::value_type& val)const{
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

	template <size_t N, class T>
	ArrayND<N, T>& ArrayND<N, T>::operator-=(const ArrayND<N, T>& other){
        auto o = other.begin();
		for (auto& i:*this)i-=*o++;
		return *this;
	}

	template <size_t N, class T>
	ArrayND<N, T>& ArrayND<N, T>::operator+=(const ArrayND<N, T>& other){
		auto o = other.begin();
		for (auto& i:*this)i+=*o++;
		return *this;
	}

	template <size_t N, class T>
	ArrayND<N, T>& ArrayND<N, T>::operator*=(const ArrayND<N, T>& other){
		auto o = other.begin();
		for (auto& i:*this)i*=*o++;
		return *this;
	}

	template <size_t N, class T>
	ArrayND<N, T>& ArrayND<N, T>::operator/=(const ArrayND<N, T>& other){
		auto o = other.begin();
		for (auto& i:*this)i/=*o++;
		return *this;
	}

    template <size_t N, class T>
    Prismatic::ArrayND<N, std::vector<T> > ones_ND(const std::array<size_t, N> dims){
        size_t size = 1;
        for (auto& i:dims)size*=i;
        return Prismatic::ArrayND<N, std::vector<T> >(std::vector<T>(size,1), dims);
    }

    template <size_t N, class T>
    Prismatic::ArrayND<N, std::vector<T> > zeros_ND(const std::array<size_t, N> dims){
        size_t size = 1;
        for (auto& i:dims)size*=i;
        return Prismatic::ArrayND<N, std::vector<T> >(std::vector<T>(size,0), dims);
    }

	template <class T>
	using Array2D_T = Prismatic::ArrayND<2, std::vector<T> >;
	template <class T>
	using Array1D_T = Prismatic::ArrayND<1, std::vector<T> >;
	template <class T>
	std::pair<Array2D_T<T>, Array2D_T<T>> meshgrid(const Array1D_T<T>& Y, const Array1D_T<T>& X){
		Array2D_T<T> yy = zeros_ND<2, T>({{Y.size(), X.size()}});
		Array2D_T<T> xx = zeros_ND<2, T>({{Y.size(), X.size()}});
		for (auto j = 0; j < xx.get_dimj(); ++j){
			for (auto i = 0; i < xx.get_dimi(); ++i){
				yy.at(j,i) = Y[j];
				xx.at(j,i) = X[i];
			}
		}
		return std::pair<Array2D_T<T>, Array2D_T<T> >(yy,xx);
	}

	template <>
	inline void ArrayND<3, std::vector<double> >::toMRC_f(const char* filename) const{
		// output to an MRC file in float format
		// see http://bio3d.colorado.edu/imod/doc/mrc_format.txt for details
		std::ofstream f(filename, std::ios::binary |std::ios::out);
		if (f) {
			int int_header[56];
			char char_header[800];
			std::memset((void *) char_header, 0, 800);
			std::memset((void *) int_header, 0, 56 * 4);
			int_header[0] = (int) dims[2]; //nx
			int_header[1] = (int) dims[1]; //ny
			int_header[2] = (int) dims[0]; //nz
			int_header[3] = 2; //mode, float
			f.write((char*)int_header,56*4); //use 4 instead of sizeof(int) because architecture may change but file format won't
			f.write(char_header,800);
			float* data_buffer = new float[this->size()];
			for (auto i = 0; i < this->size(); ++i)data_buffer[i] = (float)data[i];
			f.write((char*)data_buffer,this->size()*sizeof(float));
			delete[] data_buffer;
		}
	}


	template <>
	inline void ArrayND<3, std::vector<float> >::toMRC_f(const char* filename) const{
		// output to an MRC file in float format
		// see http://bio3d.colorado.edu/imod/doc/mrc_format.txt for details
		std::ofstream f(filename, std::ios::binary |std::ios::out);
		if (f) {
			int int_header[56];
			char char_header[800];
			std::memset((void *) char_header, 0, 800);
			std::memset((void *) int_header, 0, 56 * 4);
			int_header[0] = (int) dims[2]; //nx
			int_header[1] = (int) dims[1]; //ny
			int_header[2] = (int) dims[0]; //nz
			int_header[3] = 2; //mode, float
			f.write((char*)int_header,56*4); //use 4 instead of sizeof(int) because architecture may change but file format won't
			f.write(char_header,800);

			float* data_buffer = new float[this->size()];
			for (auto i = 0; i < this->size(); ++i)data_buffer[i] = data[i];
			f.write((char*)data_buffer,this->size()*sizeof(float));
			delete[] data_buffer;
		} else {
			std::cout << "error opening file " << filename << std::endl;
		}
	}

	template <>
	inline void ArrayND<2, std::vector<double> >::toMRC_f(const char* filename) const{
		// output to an MRC file in float format
		// see http://bio3d.colorado.edu/imod/doc/mrc_format.txt for details
		std::ofstream f(filename, std::ios::binary |std::ios::out);
		if (f) {
			int int_header[56];
			char char_header[800];
			std::memset((void *) char_header, 0, 800);
			std::memset((void *) int_header, 0, 56 * sizeof(int));
			int_header[0] = (int) dims[1]; //nx
			int_header[1] = (int) dims[0]; //ny
			int_header[2] = (int) 1; //nz
			int_header[3] = 2; //mode, float
			f.write((char*)int_header,56*4); //use 4 instead of sizeof(int) because architecture may change but file format won't
			f.write(char_header,800);
			float* data_buffer = new float[this->size()];
			for (auto i = 0; i < this->size(); ++i)data_buffer[i] = (float)data[i];
			f.write((char*)data_buffer,this->size()*sizeof(float));
			delete[] data_buffer;
		}
	}

	template <>
	inline void ArrayND<2, std::vector<float> >::toMRC_f(const char* filename) const{
		// output to an MRC file in float format
		// see http://bio3d.colorado.edu/imod/doc/mrc_format.txt for details
		std::ofstream f(filename, std::ios::binary |std::ios::out);
		if (f) {
			int int_header[56];
			char char_header[800];
			std::memset((void *) char_header, 0, 800);
			std::memset((void *) int_header, 0, 56 * sizeof(int));
			int_header[0] = (int) dims[1]; //nx
			int_header[1] = (int) dims[0]; //ny
			int_header[2] = (int) 1; //nz
			int_header[3] = 2; //mode, float
			f.write((char*)int_header,56*4); //use 4 instead of sizeof(int) because architecture may change but file format won't
			f.write(char_header,800);
			float* data_buffer = new float[this->size()];
			for (auto i = 0; i < this->size(); ++i)data_buffer[i] = (float)data[i];
			f.write((char*)data_buffer,this->size()*sizeof(float));
			delete[] data_buffer;
		}
	}

	template <>
	inline void ArrayND<2, std::vector<unsigned int> >::toMRC_f(const char* filename) const{
		// output to an MRC file in float format
		// see http://bio3d.colorado.edu/imod/doc/mrc_format.txt for details
		std::ofstream f(filename, std::ios::binary |std::ios::out);
		if (f) {
			int int_header[56];
			char char_header[800];
			std::memset((void *) char_header, 0, 800);
			std::memset((void *) int_header, 0, 56 * sizeof(int));
			int_header[0] = (int) dims[1]; //nx
			int_header[1] = (int) dims[0]; //ny
			int_header[2] = (int) 1; //nz
			int_header[3] = 2; //mode, float
			f.write((char*)int_header,56*4); //use 4 instead of sizeof(int) because architecture may change but file format won't
			f.write(char_header,800);
			float* data_buffer = new float[this->size()];
			for (auto i = 0; i < this->size(); ++i)data_buffer[i] = (float)data[i];
			f.write((char*)data_buffer,this->size()*sizeof(float));
			delete[] data_buffer;
		}
	}


    // We want to prevents programming errors like querying the size of the 3rd dimension of
    // a 2D object. Rather than introducing a check at runtime, I'll just delete them. Unfortunately,
    // this does require full template instantiation, so each type has to be done separately.
	template <>
	size_t ArrayND<1, std::vector<double> >::get_dimj()   const = delete;
	template <>
	size_t ArrayND<1, std::vector<double> >::get_dimk()   const = delete;
	template <>
	size_t ArrayND<1, std::vector<double> >::get_diml()   const = delete;
	template <>
	size_t ArrayND<1, std::vector<double> >::get_dimm()   const = delete;
	template <>
	size_t ArrayND<2, std::vector<double> >::get_dimk()   const = delete;
	template <>
	size_t ArrayND<2, std::vector<double> >::get_diml()   const = delete;
	template <>
	size_t ArrayND<2, std::vector<double> >::get_dimm()   const = delete;
	template <>
	size_t ArrayND<3, std::vector<double> >::get_diml()   const = delete;
	template <>
	size_t ArrayND<3, std::vector<double> >::get_dimm()   const = delete;
	template <>
	size_t ArrayND<4, std::vector<double> >::get_dimm()   const = delete;

	template <>
	size_t ArrayND<1, std::vector< std::complex<double> > >::get_dimj()   const = delete;
	template <>
	size_t ArrayND<1, std::vector< std::complex<double> > >::get_dimk()   const = delete;
	template <>
	size_t ArrayND<1, std::vector< std::complex<double> > >::get_diml()   const = delete;
	template <>
	size_t ArrayND<1, std::vector< std::complex<double> > >::get_dimm()   const = delete;
	template <>
	size_t ArrayND<2, std::vector< std::complex<double> > >::get_dimk()   const = delete;
	template <>
	size_t ArrayND<2, std::vector< std::complex<double> > >::get_diml()   const = delete;
	template <>
	size_t ArrayND<2, std::vector< std::complex<double> > >::get_dimm()   const = delete;
	template <>
	size_t ArrayND<3, std::vector< std::complex<double> > >::get_diml()   const = delete;
	template <>
	size_t ArrayND<3, std::vector< std::complex<double> > >::get_dimm()   const = delete;
	template <>
	size_t ArrayND<4, std::vector< std::complex<double> > >::get_dimm()   const = delete;

    template <>
    size_t ArrayND<1, std::vector<float> >::get_dimj()   const = delete;
    template <>
    size_t ArrayND<1, std::vector<float> >::get_dimk()   const = delete;
    template <>
    size_t ArrayND<1, std::vector<float> >::get_diml()   const = delete;
    template <>
    size_t ArrayND<1, std::vector<float> >::get_dimm()   const = delete;
    template <>
    size_t ArrayND<2, std::vector<float> >::get_dimk()   const = delete;
    template <>
    size_t ArrayND<2, std::vector<float> >::get_diml()   const = delete;
    template <>
    size_t ArrayND<2, std::vector<float> >::get_dimm()   const = delete;
    template <>
    size_t ArrayND<3, std::vector<float> >::get_diml()   const = delete;
    template <>
    size_t ArrayND<3, std::vector<float> >::get_dimm()   const = delete;
    template <>
    size_t ArrayND<4, std::vector<float> >::get_dimm()   const = delete;

    template <>
    size_t ArrayND<1, std::vector< std::complex<float> > >::get_dimj()   const = delete;
    template <>
    size_t ArrayND<1, std::vector< std::complex<float> > >::get_dimk()   const = delete;
    template <>
    size_t ArrayND<1, std::vector< std::complex<float> > >::get_diml()   const = delete;
    template <>
    size_t ArrayND<1, std::vector< std::complex<float> > >::get_dimm()   const = delete;
    template <>
    size_t ArrayND<2, std::vector< std::complex<float> > >::get_dimk()   const = delete;
    template <>
    size_t ArrayND<2, std::vector< std::complex<float> > >::get_diml()   const = delete;
    template <>
    size_t ArrayND<2, std::vector< std::complex<float> > >::get_dimm()   const = delete;
    template <>
    size_t ArrayND<3, std::vector< std::complex<float> > >::get_diml()   const = delete;
    template <>
    size_t ArrayND<3, std::vector< std::complex<float> > >::get_dimm()   const = delete;
    template <>
    size_t ArrayND<4, std::vector< std::complex<float> > >::get_dimm()   const = delete;

    template <>
    size_t ArrayND<1, std::vector<size_t> >::get_dimj()   const = delete;
    template <>
    size_t ArrayND<1, std::vector<size_t> >::get_dimk()   const = delete;
    template <>
    size_t ArrayND<1, std::vector<size_t> >::get_diml()   const = delete;
    template <>
    size_t ArrayND<1, std::vector<size_t> >::get_dimm()   const = delete;
    template <>
    size_t ArrayND<2, std::vector<size_t> >::get_dimk()   const = delete;
    template <>
    size_t ArrayND<2, std::vector<size_t> >::get_diml()   const = delete;
    template <>
    size_t ArrayND<2, std::vector<size_t> >::get_dimm()   const = delete;
    template <>
    size_t ArrayND<3, std::vector<size_t> >::get_diml()   const = delete;
    template <>
    size_t ArrayND<3, std::vector<size_t> >::get_dimm()   const = delete;
    template <>
    size_t ArrayND<4, std::vector<size_t> >::get_dimm()   const = delete;

}


#endif //PRISM_ARRAYND_H
