from distutils.core import setup, Extension

pyprimsatic = Extension('pyprismatic',
	sources=[
	'pyprismatic.cpp',
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
	'src/WorkDispatcher.cpp'],
	include_dirs=["./include", "/usr/local/boost_1_60_0/"],
	extra_compile_args=['-std=c++11'],
	define_macros=[],
	library_dirs=["/usr/local/lib/"],
	libraries=["fftw3f", "fftw3f_threads"])

setup(name = 'PyPrismatic', 
	  version = '1.0.01',
	  description="Python wrapper for Prismatic\
	 			   package for fast image simulation using the PRISM and multislice\
	 			   algorithms in Scanning Transmission Electron Microscopy (STEM)",
	  ext_modules=[pyprimsatic])
