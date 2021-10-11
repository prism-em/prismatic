# Copyright Alan (AJ) Pryor, Jr. 2017
# Transcribed from MATLAB code by Colin Ophus
# Prismatic is distributed under the GNU General Public License (GPL)
# If you use Prismatic, we kindly ask that you cite the following papers:

# 1. Ophus, C.: A fast image simulation algorithm for scanning
#    transmission electron microscopy. Advanced Structural and
#    Chemical Imaging 3(1), 13 (2017)

# 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
#    Implementation of Image Simulation Algorithms for Scanning
# 	 Transmission Electron Microscopy. arXiv:1706.08563 (2017)

from . import core  # noqa
from . import fileio  # noqa
from . import process
from pyprismatic.params import Metadata

def keySearch(dictionary,layer):
    layer+=1
    try:
        keys = dictionary.keys()
        for key in keys:
            printStr = ""
            for i in range(0,layer):
                printStr += "--"
            printStr += key
            print(printStr)
            keySearch(dictionary[key],layer)
    except:
        a = 1

def demo():
    import os

    with open("temp.XYZ", "w") as fid:
        fid.write(
            "one unit cell of 100 silicon\n\
  5.43    5.43    5.43\n\
14  0.0000  0.0000  0.0000  1.0  0.076\n\
14  2.7150  2.7150  0.0000  1.0  0.076\n\
14  1.3575  4.0725  1.3575  1.0  0.076\n\
14  4.0725  1.3575  1.3575  1.0  0.076\n\
14  2.7150  0.0000  2.7150  1.0  0.076\n\
14  0.0000  2.7150  2.7150  1.0  0.076\n\
14  1.3575  1.3575  4.0725  1.0  0.076\n\
14  4.0725  4.0725  4.0725  1.0  0.076\n\
-1"
        )
    meta = Metadata(filenameAtoms="temp.XYZ", filenameOutput="demo.h5")
    meta.algorithm = "multislice"
    meta.go()

    try:
        import h5py

        demoFile = h5py.File('demo.h5','r')
        print('demo.h5 filestructure:')
        keySearch(demoFile, 0)
    except ImportError:
        pass

    os.remove("temp.XYZ")
