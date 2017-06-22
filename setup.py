from setuptools import setup, Extension
from setuptools.command.install import install
import sys

global libs
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
	print("GPU ENABLED")
	# libs.extend("bla")
else:
	print("GPU disabled")

pyprimsatic_core = Extension('pyprismatic.core',
	sources=[
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
	'src/WorkDispatcher.cpp'],
	include_dirs=["./include"],
	extra_compile_args=['-std=c++11'],
	define_macros=[],
	library_dirs=["/usr/local/lib/"],
	libraries=["fftw3f", "fftw3f_threads"],
	cmdclass={'install': InstallCommand})


setup(name = 'PyPrismatic',
	author = 'Alan (AJ) Pryor, Jr.', 
	author_email='apryor6@gmail.com',
	version = '1.0.5',
	description="Python wrapper for Prismatic package for\
	 fast image simulation using the PRISM and multislice algorithms in\
	  Scanning Transmission Electron Microscopy (STEM)",
	ext_modules=[pyprimsatic_core],
	packages=['pyprismatic'],
    install_requires=['numpy>=1.13.0','matplotlib>=2.0.2'],
    cmdclass={'install': InstallCommand})
