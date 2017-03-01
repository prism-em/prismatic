# PRISM
C++/CUDA implementation of compact S-matrix formulism for fast Multislice simulation of electron micrographs. PRISM is capable of producing simulated micrographs with tolerable error compared with the original implementation by [Kirkland](http://people.ccmr.cornell.edu/~kirkland/) while providing computational acceleration of over ###x.


## LINUX/OS X Installation

*\AJ: Doing the above should build all binary executables
within the source tree. Currently, this includes a
number of tests for creating PRISM::Array2D objects
with std/Thrust backends for CPU and GPU arrays in 
./Array2D_test/. There is also a basic CUDA example
for doubling an array in ./CUDA_test/ that shows how
to do so using both low level CUDA and the higher 
level Thrust API calls. The main executable is
PRISM at the top of the source tree -- one day
this will be either an entry point to a CLI or 
a GUI, but for now is a placeholder and computes
an FFT*

PRISM is built on top of [CMake](https://cmake.org/), a cross-platform compilation utility that allows a single source tree to be compiled into a variety of formats including UNIX Makefiles, Microsoft Visual Studio projects, Mac OS XCode projects, etc. Only having to maintain one project means PRISM developers spend less time managing multiple code-bases and more time optimizing, debugging, and extending PRISM, resulting in a better end-user experience.  

To install PRISM, you must first [install Cmake](https://cmake.org/install/) and [FFTW](http://www.fftw.org/fftw2_doc/fftw_6.html). Once that is finished, open a terminal and get the PRISM source from Github using `git clone`:

```
git clone git@github.com:apryor6/PRISM.git
cd PRISM/
```

Conventional CMake practice is to use out-of-source builds, which means we will compile the source code into a separate directory. This has a number of advantages including providing the flexibility to build multiple version of PRISM (such as compiling with/without GPU support), and allowing for easier cleanup. First, make a build directory (the name doesn't matter) at the top of the PRISM source tree.

```
mkdir build
cd build
```
Then invoke CMake

```
cmake ../
```

This will generate a Makefile with the necessary dependencies and paths to compile PRISM. Finally, compile and install PRISM with:

```
make
```

If this completes successfully, you can then run PRISM with the following syntax

```
./PRISM /path/to/atoms.txt interpolation_factor
```

where `/path/to/atoms.txt` file is a txt file containing csv values x,y,z,Z for each atom (1 per row), and `interpolation_factor` is PRISM's *f* parameter.

## Building Test Programs
PRISM contains a number of test programs that can be useful
for debugging. They can be built by setting the 
CMake variable `PRISM_BUILD_TESTS=1`. This can be done on the command line with the -D flag during the cmake
call from within the build directory like so
```
cmake -DPRISM_BUILD_TESTS=1 ../
```
Alternatively, if you have already run `cmake` there will
be a CMakeCache.txt file that contains all of the option 
settings. You can edit options directly from this file or interactively edit options with `ccmake` like so
```$xslt
ccmake .
```
After any changes, rerun `make` to rebuild with the current options.
## Implementation

*Note: we following the standard NVIDIA naming convention and interchangeably refer to the CPU as the "host" and to the GPU as the "device".*

PRISM is implemented in C++/CUDA using mainly the [C++ standard library](http://en.cppreference.com/w/) and the [CUDA Thrust library](https://github.com/thrust/thrust). Host FFTs are computed using the [FFTW library](http://www.fftw.org/). 

Multi-dimensional arrays in PRISM such as [PRISM::Array2D](Array2D.h) are implemented as a thin wrapper containing the size of each dimension and a vector-like array (that is, a [std::vector](http://en.cppreference.com/w/cpp/container/vector), [thrust::host_vector](https://thrust.github.io/doc/classthrust_1_1host__vector.html), or [thrust::device_vector](https://thrust.github.io/doc/classthrust_1_1device__vector.html)). These vector objects must support at least `.begin()` and `.end()` iterators as PRISM uses modern C++11 constructs such as [range-based for loops](http://en.cppreference.com/w/cpp/language/range-for). The reason for using a custom array class is to maximize performance by ensuring the underlying data is contiguous in memory while also leveraging the host-device memory transfer overhead convenience provided by `thrust::device_vector` and `thrust::host_vector`.


