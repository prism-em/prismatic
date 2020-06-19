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
from setuptools.command.develop import develop
from setuptools.command.build_ext import build_ext
import sys
import os
import platform
import subprocess

# Filename for the C extension module library
c_module_name = 'pyprismatic.core'

# Parse command line flags
options = {k: 'OFF' for k in ['--opt', '--debug', '--enable-gpu']}
for flag in options.keys():
    if flag in sys.argv:
        options[flag] = 'ON'
        sys.argv.remove(flag)

# Command line flags forwarded to CMake
cmake_cmd_args = []
for f in sys.argv:
    if f.startswith('-D'):
        cmake_cmd_args.append(f)

for f in cmake_cmd_args:
    sys.argv.remove(f)

# In CPU-only mode, all Prismatic C++ source files will be compiled into the Python package. For GPU support, a shared library is
# built from the CUDA/C++ sources and the Python package links against it. The C++ files that are potentially compiled into
# the shared library (along with CUDA files) are the "extra" sources and in CPU-mode they are added to the Python package sources.
# The reason for doing it this way is for user convenience so that compiling/linking to an extra Prismatic library is not required until
# enabling for GPU, at which point it is already necessary for the other CUDA libraries.

class InstallCommand(install):
    user_options = install.user_options + [("enable-gpu", None, None)]

    def initialize_options(self):
        install.initialize_options(self)
        self.enable_gpu = False

    def finalize_options(self):
        install.finalize_options(self)

    def run(self):
        install.run(self)

class DevelopCommand(develop):
    user_options = develop.user_options + [("enable-gpu", None, None)]

    def initialize_options(self):
        develop.initialize_options(self)
        self.enable_gpu = False

    def finalize_options(self):
        develop.finalize_options(self)

    def run(self):
        develop.run(self)

class CMakeExtension(Extension):
    def __init__(self, name, cmake_lists_dir='.', **kwa):
        Extension.__init__(self, name, sources=[], **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)

class cmake_build_ext(build_ext):
    def build_extensions(self):
        # Ensure that CMake is present and working
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError('Cannot find CMake executable')

        for ext in self.extensions:

            extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
            cfg = 'Debug' if options['--debug'] == 'ON' else 'Release'
            if options['--enable-gpu'] == 'ON':
                gpu_enable = '1'
            else:
                gpu_enable = '0'

            cmake_args = [
                '-DCMAKE_BUILD_TYPE=%s' % cfg,
                # Ask CMake to place the resulting library in the directory
                # containing the extension
                '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                # Other intermediate static libraries are placed in a
                # temporary build directory instead
                '-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), self.build_temp),
                # Hint CMake to use the same Python executable that
                # is launching the build, prevents possible mismatching if
                # multiple versions of Python are installed
                '-DPYTHON_EXECUTABLE={}'.format(sys.executable),
                '-DPRISMATIC_ENABLE_PYPRISMATIC=1',
                '-DPRISMATIC_ENABLE_CLI=0',
                '-DPRISMATIC_ENABLE_GPU={}'.format(gpu_enable)
            ]

            # We can handle some platform-specific settings at our discretion
            if platform.system() == 'Windows':
                plat = ('x64' if platform.architecture()[0] == '64bit' else 'Win32')
                cmake_args += [
                    # These options are likely to be needed under Windows
                    '-DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=TRUE',
                    '-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir),
                ]
                # Assuming that Visual Studio and MinGW are supported compilers
                if self.compiler.compiler_type == 'msvc':
                    # When using Visual Studio generator, set the platform
                    # otherwise don't since other generator, such as NMake,
                    # doesn't support platform specification.
                    if 'Visual Studio' in os.environ.get('CMAKE_GENERATOR', ''):
                        cmake_args += [
                            '-DCMAKE_GENERATOR_PLATFORM=%s' % plat,
                        ]
                else:
                    cmake_args += [
                        '-G', 'MinGW Makefiles',
                    ]

            cmake_args += cmake_cmd_args

            if not os.path.exists(self.build_temp):
                os.makedirs(self.build_temp)

            # Config
            subprocess.check_call(['cmake', ext.cmake_lists_dir] + cmake_args,
                                  cwd=self.build_temp)

            # Build
            subprocess.check_call(['cmake', '--build', '.', '--config', cfg],
                                  cwd=self.build_temp)

setup(
    name="PyPrismatic",
    author="Alan (AJ) Pryor, Jr.",
    author_email="apryor6@gmail.com",
    version="1.2.0",
    description="Python wrapper for Prismatic package for fast image simulation using the PRISM and multislice algorithms in Scanning Transmission Electron Microscopy (STEM)",
    ext_modules=[CMakeExtension(c_module_name)],
    packages=["pyprismatic"],
    install_requires=["numpy>=1.13.0", "matplotlib>=2.0.2","h5py>=2.9.0"],
    cmdclass={"install": InstallCommand,
    "develop":DevelopCommand,
    "build_ext":cmake_build_ext},
)
