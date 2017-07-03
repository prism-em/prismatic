# Copyright Alan (AJ) Pryor, Jr. 2017
# Transcribed from MATLAB code by Colin Ophus
# Prismatic is distributed under the GNU General Public License (GPL)
# If you use Prismatic, we kindly ask that you cite the following papers:

# 1. Ophus, C.: A fast image simulation algorithm for scanning
#    transmission electron microscopy. Advanced Structural and
#    Chemical Imaging 3(1), 13 (2017)

# 2. Pryor, Jr., A., Ophus, C., and Miao, J.: A Streaming Multi-GPU
#    Implementation of Image Simulation Algorithms for Scanning
#	 Transmission Electron Microscopy. arXiv:1706.08563 (2017)

from . import core
from . import fileio
from pyprismatic.params import Metadata

def demo():
	import os 
	with open('temp.XYZ', 'w') as fid:
		fid.write("one unit cell of 100 silicon\n\
  5.43    5.43    5.43\n\
14  0.0000  0.0000  0.0000  1.0  0.076\n\
14  2.7150  2.7150  0.0000  1.0  0.076\n\
14  1.3575  4.0725  1.3575  1.0  0.076\n\
14  4.0725  1.3575  1.3575  1.0  0.076\n\
14  2.7150  0.0000  2.7150  1.0  0.076\n\
14  0.0000  2.7150  2.7150  1.0  0.076\n\
14  1.3575  1.3575  4.0725  1.0  0.076\n\
14  4.0725  4.0725  4.0725  1.0  0.076\n\
-1")
	meta = Metadata(filenameAtoms='temp.XYZ',
				    filenameOutput='output.mrc');
	meta.go()
	import numpy as np
	from pyprismatic.fileio import readMRC
	import matplotlib.pyplot as plt
	result = readMRC("output.mrc")
	plt.figure()
	plt.imshow(np.squeeze(np.sum(result,axis=2)))
	plt.show()
	os.remove("temp.XYZ")