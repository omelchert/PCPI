""" FILE: dataIO.py

Module implementing functions to parse input data contained in file 

Author: O. Melchert
Date:   25.06.2016 
"""
import sys
import numpy as np
import scipy 
import scipy.io

def fetchData_conv(inFileName):
        """retrieve data from file

        \param inFileName path to file where data is stored in
        """
        rArray   = []
        zArray   = []
        ArzArray = []
        with open(inFileName,'r') as f:
              # skip header
              h = f.readline()
              # sift through remainder
              for line in f:
                 c = line.split()
                 r,z,Arz = float(c[0]), float(c[1]), float(c[2])
                 rArray.append(r)
                 zArray.append(z)
                 ArzArray.append(Arz)
        return np.array(rArray,dtype=np.float),\
               np.array(zArray,dtype=np.float),\
               np.array(ArzArray,dtype=np.float)


# EOF: dataIO.py
