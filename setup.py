# Copyright Alan (AJ) Pryor, Jr. 2017
# Transcribed from MATLAB code by Colin Ophus
# Prismatic is distributed under the GNU General Public License (GPL)
# If you use Prismatic, we kindly ask that you cite the following papers:

# 1. Ophus, C.: A fast image simulation algorithm for scanning
#    transmission electron microscopy. Advanced Structural and
#    Chemical Imaging 3(1), 13 (2017)

# 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
#    Implementation of Image Simulation Algorithms for Scanning
#    Transmission Electron Microscopy. arXiv:1706.08563 (2017)

from setuptools import setup, Extension
from setuptools.command.install import install
import sys
import os

prismatic_extra_definitions = []
prismatic_libs  		    = []

# In CPU-only mode, all Prismatic C++ source files will be compiled into the Python package. For GPU support, a shared library is
# built from the CUDA/C++ sources and the Python package links against it. The C++ files that are potentially compiled into
# the shared library (along with CUDA files) are the "extra" sources and in CPU-mode they are added to the Python package sources.
# The reason for doing it this way is for user convenience so that compiling/linking to an extra Prismatic library is not required until
# enabling for GPU, at which point it is already necessary for the other CUDA libraries.
prismatic_sources_base  = ['pyprismatic/core.cpp'
						  ]
prismatic_sources_extra = [
	'src/go.cpp',
	'src/atom.cpp',
	'src/configure.cpp',
	'src/Multislice_calcOutput.cpp',
	'src/Multislice_entry.cpp',#
	'src/parseInput.cpp',#
	'src/PRISM01_calcPotential.cpp',
	'src/PRISM02_calcSMatrix.cpp',#
	'src/PRISM03_calcOutput.cpp',#
	'src/PRISM_entry.cpp',
	'src/projectedPotential.cpp',#
	'src/utility.cpp',#
	'src/WorkDispatcher.cpp'
							]
prismatic_sources 		    = prismatic_sources_base
prismatic_include_dirs 		= ["./include"]
prismatic_library_dirs 		= []

if os.name == "nt": #check for Windows OS
	prismatic_libs 		   		= ["libfftw3f-3"]
else:
	prismatic_libs 		   		= ["fftw3f", "fftw3f_threads"]
prismatic_extra_compile_defs= ['-std=c++11']

class InstallCommand(install):
	user_options = install.user_options + [
		('enable-gpu', None, None),
	]

	def initialize_options(self):
		install.initialize_options(self)
		self.enable_gpu = False

	def finalize_options(self):
		install.finalize_options(self)

	def run(self):
		install.run(self)


if ("--enable-gpu" in sys.argv): # GPU-mode, add some macro definitions and link with the CUDA libraries
	prismatic_extra_definitions.extend([("PRISMATIC_ENABLE_GPU",1), ("BUILDING_CUPRISMATIC",1)])
	prismatic_libs.extend(["cuprismatic", "cufft", "cudart"])
else: # CPU-only mode, add the extra C++ sources to the python package
	prismatic_sources.extend(prismatic_sources_extra)

pyprimsatic_core = Extension('pyprismatic.core',
	sources=prismatic_sources,
	include_dirs=prismatic_include_dirs,
	extra_compile_args=prismatic_extra_compile_defs,
	define_macros=prismatic_extra_definitions,
	library_dirs=prismatic_library_dirs,
	libraries=prismatic_libs)

setup(name = 'PyPrismatic',
	author = 'Alan (AJ) Pryor, Jr.', 
	author_email='apryor6@gmail.com',
	version = '1.1.11',
	description="Python wrapper for Prismatic package for\
	 fast image simulation using the PRISM and multislice algorithms in\
	  Scanning Transmission Electron Microscopy (STEM)",
	ext_modules=[pyprimsatic_core],
	packages=['pyprismatic'],
    install_requires=['numpy>=1.13.0','matplotlib>=2.0.2'],
    cmdclass={'install': InstallCommand})
