from setuptools import setup, Extension
from setuptools.command.install import install
import sys

prismatic_extra_definitions = []
prismatic_libs  		    = []
prismatic_sources 		    = [
							'pyprismatic/core.cpp',
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
prismatic_include_dirs 		= ["./include"]
prismatic_library_dirs 		= []
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


if ("--enable-gpu" in sys.argv):
	# print("GPU ENABLED")
	prismatic_extra_definitions.extend([("PRISMATIC_ENABLE_GPU",1)])
	prismatic_libs.extend(["cuprismatic", "cufft", "cudart"])
	# prismatic_include_dirs.extend(["/usr/local/cuda-8.0/include/"])
	# prismatic_library_dirs.extend(["/usr/local/cuda-8.0/lib64/"])

print("prismatic_libs = " , prismatic_libs)
pyprimsatic_core = Extension('pyprismatic.core',
	sources=prismatic_sources,
	include_dirs=prismatic_include_dirs,
	extra_compile_args=prismatic_extra_compile_defs,
	define_macros=prismatic_extra_definitions,
	library_dirs=prismatic_library_dirs,
	libraries=prismatic_libs)
	# cmdclass={'install': InstallCommand})


setup(name = 'PyPrismatic',
	author = 'Alan (AJ) Pryor, Jr.', 
	author_email='apryor6@gmail.com',
	version = '1.1.1',
	description="Python wrapper for Prismatic package for\
	 fast image simulation using the PRISM and multislice algorithms in\
	  Scanning Transmission Electron Microscopy (STEM)",
	ext_modules=[pyprimsatic_core],
	packages=['pyprismatic'],
    install_requires=['numpy>=1.13.0','matplotlib>=2.0.2'],
    cmdclass={'install': InstallCommand})
