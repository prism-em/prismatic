# Installing `Prismatic`  
---
Table of Contents    
	- [Dependencies](#dependencies)  
	- [Installers](#installers)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  - [Mac OS X](#binary-installers-mac)  
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  - [Windows](#binary-installers-win)  
	- [Building `Prismatic` from the source code](#from-source)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [Setting environmental variables](#environmental-setup)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [Linux](#environmental-setup-linux)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [Mac OS X](#environmental-setup-mac)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [Windows](#environmental-setup-win)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [Installing on Linux](#linux-install)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  - [Compiling with `CMake` from the command line](#compiling-linux)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [Installing on Mac OS X](#mac-install)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  - [Compiling with `CMake` from the command line](#compiling-mac)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; - [Installing on Windows](#win-install)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  - [Compiling with `CMake` from the command line](#compiling-win)  
	- [Python: Installing `PyPrismatic`](#python-installing-pyprismatic)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  - [Building the `cuPrismatic` library](#cuprismatic)  
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  - [Installing `PyPrismatic` with Pip](#installing-with-pip-cuda)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  - [Installing with setup.py](#installing-with-setup)  
	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  - [Python: Testing `PyPrismatic`](#testing-PyPrismatic)  
	- [Setting `CMake` Options](#setting-cmake-options)  
	- [List of `Prismatic` `CMake` options](#make-options)  
	- [Command line options](#cli-options) 



## Dependencies

The following dependencies are needed by `Prismatic`:

*Required*

* [CMake](https://CMake.org/) (*For compiling the source code*)  
* [Boost](http://www.boost.org)  
* [FFTW](http://www.fftw.org) (compiled with `--enable-float`, `--enable-shared`, and `--enable-threads`)  

*Optional*

* [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit) (*For GPU support*)  
* Python 3, a good choice is [the Anaconda distribution](https://www.continuum.io/downloads) (*For the python package*, `PyPrismatic`)  
* [Qt 5](https://www.qt.io/) (*For building the GUI*)  

`Prismatic` was developed using CUDA 8.0, but likely works with older versions as well and we welcome feedback from any user who attempts to do so (CUDA 7.0, 7.5 also have been reported to work).
*Note: Even if you download a binary version of the GPU codes, you will still need to have the CUDA toolkit installed so that the `cuFFT` libraries can be found at runtime.*


## Installers
<a name="binary-installers-mac"></a>
### Binary Installers (Mac OS X) 


<a name="binary-installers-win"></a>
### Binary Installers (Windows) 

*Links to binary installers should go here in the future*

*Links to binary installers should go here in the future*

<a name="from-source"></a>
## Building `Prismatic` from the source code

`Prismatic` is built using [`CMake`](https://`CMake`.org/), a cross-platform compilation utility that 
allows a single source tree to be compiled into a variety of formats including UNIX Makefiles, 
Microsoft Visual Studio projects, Mac OS XCode projects, etc.


<a name="get-source-code"></a>
### Getting the source code 

Once the [dependencies](#dependencies) are installed get the `Prismatic` source either from [compressed source files](https://github.com/prism-em/`Prismatic`/archive/master.zip) or directly 
from [Github](https://github.com/prism-em/`Prismatic`) using `git clone`. Next, follow the instructions appropriate for your operating system.

<a name="environmental-setup"></a>
## Environmental setup

`CMake` and/or `setuptools` needs to know where to find the [dependencies](#dependencies) in order to build `Prismatic` and `PyPrismatic`, respectively. In my opinion the easiest way to do this is by setting environmental variables. You can also manually provide paths to compilation commands, but if you take the time to setup your environment then everything should be automatically found. Follow the appropriate instructions below for your operating system.

<a name="environmental-setup-linux"></a>

#### Linux

These are the relevant environmental variables on Linux

* `CPLUS_INCLUDE_PATH` - default search location for C++ header files  
* `LIBRARY_PATH` - default search location for C/C++ libraries  
* `LD_LIBRARY_PATH` - default search location for shared libraries to be loaded at runtime 

These variables can be set from the terminal with the `export` command. For example, to prepend an example header search path for boost "/usr/local/boost_1_60_0", header search for CUDA located at "/usr/local/cuda-8.0/include", and CUDA library path at "/usr/local/cuda-8.0/lib64" one could invoke

~~~
export CPLUS_INCLUDE_PATH=/usr/local/boost_1_60_0:/usr/local/cuda-8.0/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=/usr/local/cuda-8.0/lib64:$LIBRARY_PATH
~~~

You can make these changes persistent by adding the same lines to your `~/.bashrc` file so that they are executed at startup every time you open a terminal. 

Depending how your system is configured and what portions of `Prismatic` you are building you may need to add additional paths. For example, if you are building the GUI, you will also need to provide the paths to Qt5 headers and libraries. See the [dependencies](#dependencies) for details about what is required.


<a name="environmental-setup-mac"></a>

#### Mac OS X

These are the relevant environmental variables on Mac OS X

* `CPLUS_INCLUDE_PATH` - default search location for C++ header files  
* `LIBRARY_PATH` - default search location for C/C++ libraries  
* `DYLD_LIBRARY_PATH` - default search location for shared libraries to be loaded at runtime 

These variables can be set from the terminal with the `export` command. For example, to prepend an example header search path for boost "/usr/local/boost_1_60_0", header search for CUDA located at "/usr/local/cuda-8.0/include", and CUDA library path at "/usr/local/cuda-8.0/lib64" one could invoke

~~~
export CPLUS_INCLUDE_PATH=/usr/local/boost_1_60_0:/usr/local/cuda-8.0/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=/usr/local/cuda-8.0/lib64:$LIBRARY_PATH
~~~

You can make these changes persistent by adding the same lines to your `~/.bash_profile` file so that they are executed at startup every time you open a terminal. 

Depending how your system is configured and what portions of `Prismatic` you are building you may need to add additional paths. For example, if you are building the GUI, you will also need to provide the paths to Qt5 headers and libraries. See the [dependencies](#dependencies) for details about what is required.

<a name="environmental-setup-win"></a>

#### Windows

These are the relevant environmental variables on Windows

* `PATH` - default search location for executables (.exe.) and shared libraries (.dll)
* `INCLUDE` - default search location for C++ headers  
* `LIB` - default search location for C++ libraries

These environmental variables can be set graphically through system settings. The specfic details of how to this will vary depending on which version of Windows you are using, but a quick Google search should be able to provide you step-by-step instructions. For example, on Windows 10, typing "variable" into the search feature on the taskbar reveals "Edit the system environmental variables".

Depending how your system is configured and what portions of `Prismatic` you are building you may need to add additional paths. For example, if you are building the GUI, you will also need to provide the paths to Qt5 headers and libraries. See the [dependencies](#dependencies) for details about what is required.

<a name="linux-install"></a>
## Installing on Linux  

<a name="compiling-linux"></a>
### Compiling with `CMake` from the command line on Ubuntu Linux

To build `Prismatic` from the command line with `CMake`, open a terminal and navigate to the top of the source directory 

```
cd /path/to/`Prismatic`/
```

Conventional `CMake` practice is to use out-of-source builds, which means we will compile the source code into
a separate directory. This has a number of advantages including providing the flexibility to build multiple
version of `Prismatic` (such as compiling with/without GPU support), and allowing for easier cleanup. First, 
make a build directory (the name doesn't matter) at the top of the `Prismatic` source tree.

```
mkdir build
cd build
```
Then invoke `CMake`

```
`CMake` ../
```

This will generate a Makefile  project with the necessary dependencies and paths to compile `Prismatic`. The default
behavior is to build only the CLI without GPU support. These options can be enabled as described in later sections.
Finally, compile and install `Prismatic` with:

```
make
```

For faster compilation, add the `-j` switch to `make` to use multiple threads, for example to compile with 8 threads

```
make -j 8
```

If this succeeds, the executable file `Prismatic` has been built and can be run from within the build directory. To permanently 
install them, invoke

``` 
make install
```

which may require `sudo` privileges. This will place the files in `/usr/local/bin`. 

`CMake` will attempt to locate the various dependencies needed by `Prismatic` (see the section on [setting up your environment](#environmental-setup)), but if it cannot then it will produce an error and set the variable to NOTFOUND. For example, if the `Boost_INCLUDE_DIR` (the location of the Boost libraries), is not found, it will be set to `Boost_INCLUDE_DIR-NOTFOUND`. You will need to manually set the path to boost (see [below](#cmake-options) for how to set options), and then rerun `CMake`.


<a name="mac-install"></a>
## Installing on Mac OS X  


<a name="compiling-mac"></a>
### Compiling with `CMake` from the command line on OS X

*If you prefer a graphical approach, you can use the ``CMake`-gui` and follow analagous steps as [in Windows](#compiling-win)*

To build `Prismatic` from the command line with `CMake`, open a terminal and navigate to the top of the source directory 

```
cd /path/to/`Prismatic`/
```

Conventional `CMake` practice is to use out-of-source builds, which means we will compile the source code into
a separate directory. This has a number of advantages including providing the flexibility to build multiple
version of `Prismatic` (such as compiling with/without GPU support), and allowing for easier cleanup. First, 
make a build directory (the name doesn't matter) at the top of the `Prismatic` source tree.

```
mkdir build
cd build
```
Then invoke `CMake`

```
`CMake` ../
```

This will generate a Makefile  project with the necessary dependencies and paths to compile `Prismatic`. The default
behavior is to build only the CLI without GPU support. These options can be enabled as described in later sections.
Finally, compile and install `Prismatic` with:

```
make
```

For faster compilation, add the `-j` switch to `make` to use multiple threads, for example to compile with 8 threads

```
make -j 8
```

If this succeeds, the executable file `Prismatic` has been built and can be run from within the build directory. To permanently 
install them, invoke

``` 
make install
```

which may require `sudo` privileges. This will place the files in `/usr/local/bin`. 

`CMake` will attempt to locate the various dependencies needed by `Prismatic` (see the section on [setting up your environment](#environmental-setup)), but if it cannot then it will produce an error and set the variable to NOTFOUND. For example, if the `Boost_INCLUDE_DIR` (the location of the Boost libraries), is not found, it will be set to `Boost_INCLUDE_DIR-NOTFOUND`. You will need to manually set the path to boost (see [below](#cmake-options) for how to set options), and then rerun `CMake`.


<a name="win-install"></a>
## Installing on Windows  


<a name="compiling-win"></a>
### Compiling with `CMake`'s GUI on Windows

To build `Prismatic` on Windows with the `CMake` GUI, first open `CMake` and set the location of the source code to the top level directory of `Prismatic` (this is the directory  containing "CMakeLists.txt"). Next, choose the location to build the binaries. It is recommended to create a separate folder, perhaps called "build", for this purpose. Click `Configure`, and choose the C++ compiler you would like to use. I have successfully tested Microsoft Visual Studio 2015 (64-bit) and would recommend this version if possible, particularly if you are compiling with GPU-support due to the fact that both NVIDIA and Microsoft are proprietary vendors, there is often some conflict between the newest versions of `nvcc` and MSVS. 

Based on the [option settings](cmake-options), `CMake` will then attempt to find the necessary dependencies. If you [have fully setup your environment](#environmental-setup-win), then configuration should succeed. If it fails, then variables it cannot satisfy will be set to NOTFOUND. For example, if the `Boost_INCLUDE_DIR` (the location of the Boost libraries), is not found, it will be set to `Boost_INCLUDE_DIR-NOTFOUND`. You will need to manually set the path to boost (see [below](#cmake-options) for how to set options), and then rerun `Configure`. Sometimes it is not necessary to specify a value for every "NOTFOUND" variable if there are multiple ones from the same dependency as CMake may be able to find them relative to one of the variables you do specify. For example, if you specify the path to one of the Qt libraries and reconfigure the other ones will likely be found automatically. Other `CMake` variables, like `FFTW_MPI_LIBRARY`, are not used by `Prismatic` and if they can be left "NOTFOUND". In my opinion the easiest thing to do is to specify the value of whatever variables are causing error messages one-by-one until configuration succeeds.

Once configuration is complete, click `Generate` and a MSVS .sln file will be created. Open this file, set the build mode to "Release" (*make sure you done forget this! If the build mode is not Release when you build `cuPrismatic.dll` then bad things will happen later down the road*), and then run "Build Solution" to compile the code. You may find that the build process randomly fails and then works after trying a second (or third, or fourth) time. You can then find the executable `prismatic-gui.exe` inside of the "Release" folder within the build directory you selected in the `CMake` GUI.


<a name="python-installing-pyprismatic"></a>
## Python: Installing `PyPrismatic`
`PyPrismatic` is a Python package for invoking the C++/CUDA code in `Prismatic`. It can be installed easily with `pip` provided the following dependencies are installed:  

*Required*

* [Boost](http://www.boost.org/)  
* [FFTW](http://www.fftw.org) (compiled with `--enable-float`, `--enable-shared`, and `--enable-threads`)    

*Optional*
* [The CUDA Toolkit](https://developer.nvidia.com/cuda-downloads) (*For GPU support*)    
* [The `cuPrismatic` library](#`cuPrismatic`) (*For GPU support*)  


*The optional dependencies are only necessary if you wish to use the GPU code. You will also need an NVIDIA GPU with compute capability 3.0 or greater and will add the `--enable-gpu` to the installation command. More details can be found below*  



<a name="installing-with-pip-cuda"></a>
### Installing `PyPrismatic` with Pip

If you have installed the above dependencies and [setup your environmental variables](#environmental-setup), `PyPrismatic` can be installed easily with `pip` using either 

~~~
pip install `PyPrismatic`
~~~

for CPU-only mode or for GPU mode:

~~~
pip install `PyPrismatic` --install-option="--enable-gpu"
~~~


Alternatively, you can tell `pip` where the dependencies are using `--global-option` like so (*you should change the actual names of the paths to match your machine, this is just an example*):

~~~
pip install `PyPrismatic` --global-option=build_ext --global-option="-I/usr/local/boost_1_60_0"
~~~

for CPU-only mode or for GPU mode:

~~~
pip install `PyPrismatic` --global-option=build_ext --global-option="-I/usr/local/boost_1_60_0:/usr/local/cuda-8.0/include" --global-option="-L/usr/local/cuda-8.0/lib64" --install-option="--enable-gpu"
~~~

*On Windows, a list containing multiple path names must be separated by `;` instead of `:`*

<a name="installing-with-setup"></a>

### Installing with setup.py

To install the python package from the source code with `setup.py`, first [get the source code](#get-source-code). Then navigate to the top directory (the one with `setup.py`) and invoke 

~~~
pip3 install -r requirements.txt
~~~ 

To install the Python dependencies. Next enter either (*you should change the actual names of the paths to match your machine, this is just an example*)

~~~
python3 setup.py build_ext --include-dirs=/usr/local/boost_1_60_0 install
~~~

to compile in CPU-only mode, or the following to compile the GPU version

~~~
python3 setup.py build_ext --include-dirs=/usr/local/boost_1_60_0:/usr/local/cuda-8.0/include --library-dirs=/usr/local/cuda-8.0/lib64 install --enable-gpu
~~~

*On Windows, a list containing multiple path names must be separated by `;` instead of `:`*


If you have [setup your environmental variables](#environmental-setup), you can ignore the extra arguments and just use `python3 setup.py install` or `python3 setup.py install --enable-gpu` 


<a name="cuprismatic"></a>

### Building the `cuPrismatic` library

*This step is only required for GPU support in `PyPrismatic`*

One of the great features about `CMake` is that it can easily tolerate the complicated task of compiling a single program with multiple compilers, which is necessary when mixing CUDA and Qt (discussed more [here](http://prism-em.com/source-code/#combining)).  For a CPU-only version of `PyPrismatic`, the necessary C++ source code can be distributed with the Python package, and `setuptools` can easily compile the necessary C++ extension module. Unfortunately, to my knowledge `setuptools` does not play nicely with `nvcc`, which makes distributing a Python package that utilizes custom GPU code more challenging. My solution was to compile the CUDA code into a a single shared library called `cuPrismatic`, which can then be linked against in the Python package as if it were any other C++ library. There is no "new" code in `cuPrismatic`, it simply serves as an intermediate step to help Python out by bundling the GPU code together into something it can work with. As an aside, this type of step is all `CMake` is doing under-the-hood to make CUDA and Qt play nicely in the first place -- it's just compiling the various formats of source code into a commonly understood form.

With that being said, in order to install the GPU-enabled version of `PyPrismatic`, you must first build `cuPrismatic`. To do so, you will need to [get the source code](#get-source-code), set the ``Prismatic`_ENABLE_PYTHON=1` variable in `CMake`, then configure and compile the project. More detail on this process of using `CMake` is [described above](#from-source). Once the library is installed, then proceed below.

<a name="testing-`PyPrismatic`"></a>
### Testing `PyPrismatic`

You can test your installation of `PyPrismatic` with the following commands from within `python3`

~~~ python
import `PyPrismatic` as pr
pr.demo()
~~~

*Note: If you receive an error like `ImportError: cannot import name core` but the installation process appeared to work, make sure to change directories out of the top level of the source code and try again. This occurs because the `PyPrismatic` package was built and installed globally, but if you are currently in the top level directory that contains the source code folder `PyPrismatic` then python will attempt to use that incorrectly as the package.*

<a name ="setting-cmake-options"></a>
## Setting `CMake` options

All aspects of how `Prismatic` is compiled, such as whether or not to include GUI or GPU support, are controlled through `CMake` variables.
There are at least four different ways to adjust these:

If you are using the `CMake` GUI, then options are turned on/off with check boxes
and paths are set by clicking browse and navigating to the file/folder of interest.   

If you are using the command line tool, `CMake`, then options are set with the `-D` (Define) flag. For example, 
to set `My_Variable` to 0 one would add `-DMY_VARIABLE=0` to the call to `CMake` (see the sections on enabling GUI or GPU support for more examples).

There is also the hybrid "command-line GUI" option, `ccmake`, which provides an interactive way to adjust `CMake` options from the command line.

Finally, you can also directly edit a file called "CMakeCache.txt". The first time you run `CMake` for a given project
this special file will be created containing all of the option settings. Whenever you generate a project or Makefile, the options
will be read for the "CMakeCache.txt" file, so options can be changed here directly.


**_Note_**: Any time you change `CMake` options for a particular project you must regenerate the build files and recompile
before the changes will actually take effect

### Enabling Double Precision

*Currently, this double precision setting is not supported by `PyPrismatic`*

The default behavior for `Prismatic` is to use single precision (type float). You can use double precision instead by setting ``Prismatic`_ENABLE_DOUBLE_PRECISION=1`. Note that as of this writing double precision operations are ~4x slower on the GPU, and by every test I have done the precision difference is entirely unnoticeable. However, I leave it as an option . If you find a case where using double precision is impactful, I would be very interested to hear about it.


<a name ="cmake-options"></a>
## List of `Prismatic` `CMake` options

Here's a list of the various custom options you can set and `CMake` and what they represent

* `Prismatic_ENABLE_CLI` - Build the command line interface `prismatic`"
* `Prismatic_ENABLE_GPU` - Enable GPU supprt. Requires locating CUDA headers/libraries
* `Prismatic_ENABLE_GUI` - Build the GUI, `prismatic-gui`. Requires locating Qt5 headers/libraries. A good way to do this is by setting `CMAKE_PREFIX_PATH` to point to the root of your Qt installation. For example, with Qt5 installed on a Mac through Homebrew one might add `-DCMAKE_PREFIX_PATH=/usr/local/Cellar/qt/5.9.0/` to the call to `cmake`
* `Prismatic_ENABLE_PYTHON_GPU` - Build the `cuPrismatic` shared library, which is used by the GPU version
* `Prismatic_ENABLE_DOUBLE_PRECISION` - Use type `double` for float precision. This requires locating the double precision `FFTW` libraries.


<a name ="cli-options"></a>
## Using `Prismatic` from the command line

`Prismatic` contains a command line tool, `prismatic`, that can be used to run simulations from within a terminal, bash script, etc. Building it requires the `CMake` variable `PRISM_ENABLE_CLI=1` at compilation time, which is the default behavior.

The following options are available with `prismatic` (you can also print available options with `prismatic --help`), each documented as **_long form_** **_(short form)_** *parameters* : description

* --input-file (-i) filename : the filename containing the atomic coordinates, which should be a plain text file with comma-separated values in the format x, y, z, Z 
* --output-file(-o) filename : output filename
* --interp-factor (-f) number : PRISM interpolation factor, used for both X and Y
* --interp-factor-x (-fx) number : PRISM interpolation factor in X
* --interp-factor-y (-fy) number : PRISM interpolation factor in Y
* --num-threads (-j) value : number of CPU threads to use
* --num-streams (-S) value : number of CUDA streams to create per GPU
* --num-gpus (-g) value : number of GPUs to use. A runtime check is performed to check how many are actually available, and the minimum of these two numbers is used.
* --batch-size (-b) value : number of probes/beams to propagate simultaneously for both CPU and GPU workers.
* --batch-size-cpu (-bc) value : number of probes/beams to propagate simultaneously for CPU workers.
* --batch-size-gpu (-bg) value : number of probes/beams to propagate simultaneously for GPU workers.
* --slice-thickness (-s) thickness : thickness of each slice of projected potential (in Angstroms)
* --help(-h) : print information about the available options
* --pixel-size (-p) pixel_size : size of simulation pixel size
* --detector-angle-step (-d) step_size : angular step size for detector integration bins (in mrad)
* --cell-dimension (-c) x y z : size of sample in x, y, z directions (in Angstroms)
* --tile-uc (-t) x y z : tile the unit cell x, y, z number of times in x, y, z directions, respectively
* --algorithm (-a) p/m : the simulation algorithm to use, either (p)rism or (m)ultislice
* --energy (-E) value : the energy of the electron beam (in keV)
* --alpha-max (-A) angle : the maximum probe angle to consider (in mrad)
* --potential-bound (-P) value : the maximum radius from the center of each atom to compute the potental (in Angstroms)
* --also-do-cpu-work (-C) bool=true : boolean value used to determine whether or not to also create CPU workers in addition to GPU ones
* --streaming-mode 0/1 : boolean value to force code to use (true) or not use (false) streaming versions of GPU codes. The default behavior is to estimate the needed memory from input parameters and choose automatically.
* --probe-step (-r) step_size : step size of the probe for both X and Y directions (in Angstroms)
* --probe-step-x (-rx) step_size : step size of the probe in X direction (in Angstroms)
* --probe-step-y (-ry) step_size : step size of the probe in Y direction (in Angstroms)
* --random-seed (-rs) step_size : random number seed
* --probe-xtilt (-tx) value : probe X tilt
* --probe-ytilt (-ty) value : probe X tilt
* --probe-defocus (-df) value : probe defocus (in Angstroms)
* --probe-semiangle (-sa) value : maximum probe semiangle
* --scan-window-x (-wx) min max : size of the window to scan the probe in X (in fractional coordinates between 0 and 1)
* --scan-window-y (-wy) min max : size of the window to scan the probe in Y (in fractional coordinates between 0 and 1)
* --num-FP (-F) value : number of frozen phonon configurations to calculate
* --thermal-effects (-te) bool : whether or not to include Debye-Waller factors (thermal effects)
* --save-2D-output (-2D) ang_min ang_max : save the 2D STEM image integrated between ang_min and ang_max (in mrads)
* --save-3D-output (-3D) bool=true : Also save the 3D output at the detector for each probe (3D output mode)
* --save-4D-output (-4D) bool=false : Also save the 4D output at the detector for each probe (4D output mode)
