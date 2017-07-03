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

def readMRC(filename, dtype=float, order="C"):
    """
    * readMRC *

    Read in a volume in .mrc file format. See http://bio3d.colorado.edu/imod/doc/mrc_format.txt

    :param filename: Filename of .mrc
    :return: NumPy array containing the .mrc data

    Author: Alan (AJ) Pryor, Jr.
    Jianwei (John) Miao Coherent Imaging Group
    University of California, Los Angeles
    Copyright 2015-2016. All rights reserved.

    """
    import numpy as np
    import struct
    headerIntNumber = 56
    sizeof_int = 4
    headerCharNumber = 800
    sizeof_char = 1
    with open(filename,'rb') as fid:
        int_header = struct.unpack('=' + 'i'*headerIntNumber, fid.read(headerIntNumber * sizeof_int))
        char_header = struct.unpack('=' + 'c'*headerCharNumber, fid.read(headerCharNumber * sizeof_char))
        dimz, dimy, dimx, data_flag= int_header[:4]
        if (data_flag == 0):
            datatype='u1'
        elif (data_flag ==1):
            datatype='i1'
        elif (data_flag ==2):
            datatype='f4'
        elif (data_flag ==3):
            datatype='c'
        elif (data_flag ==4):
            datatype='f4'
        elif (data_flag ==6):
            datatype='u2'
        else:
            raise ValueError("No supported datatype found!\n")

        return np.fromfile(file=fid, dtype=datatype,count=dimx*dimy*dimz).reshape((dimx,dimy,dimz),order=order).astype(dtype)
