""" FILE: mciBin.py

Module implenting function that computes volumetric energy density [J/m^3] 
from input and output of mcxyz.c C-code available under Ref. [1].

Refs:

    [1] Monte Carlo simulation of light transport in 3D heterogeneous
        tissues (mcxyz.c)
        Jacues S L, Li T
        http://omlc.org/software/mc/mcxyz/index.html 
        (acesssed 07.08.2016)

"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import numpy as np
import ctypes as ct
import os
import sys


def getEnergyDensity(fNameBase):
        """compute energy deposition for voxels.

        Reads voxel fluence rate`Fxyz` [J/m^2 per J delivered] of from 
        XXX_F.bin file, voxel tissue types `Txyz` from XXX_T.bin file as well 
        as absorption coefficients `mua` from XXX_H.mci file and computes 
        volumetric energy density `Hxyz` [J/m^3 per J delivered] for each 
        voxel.

        Args:
            fNameBase (str): Basename of file 3-tuple consisting of 
                `fNameBase_H.mci`, `fNameBase_F.bin`, and, `fNameBase_T.bin`.

        Returns:
            (x,y,z) (3-tuple, numpy arrays): Equi-spaced x, y, z grids [cm]. 
            Hxyz (numpy array, ndim=3): Volumetric energy density [J/m^3 per 
                J delivered] per voxel. 

        Notes:
            The 3-tuple of files represent input (`fBaseName_H.mci`, 
            `fBaseName_T.bin`) and output (`fBaseName_F.bin`) of the 
            mcxyz.c C-code available under Ref. [1].
        
            Ref. [2] writes "The C-cpde mcxyz.c program reads the input files,
            executes the simulation, and saves a file that holds the spatial
            distribution of deposited energy, W(x,y,z) [J/cm^3 per J delivered]
            or [1/cm^3]". Note that this is not entierly true! The C-code 
            computes W(x,y,z), but saves F(x,y,z) = W(x,y,z) / mua(x,y,z) 
            having units [J/m^2 per J delivered] as _F.bin file. Hence, the 
            auxiliary information in _T.bin and _H.mci is needed to compute
            the volumetric energy density per voxel.

        Refs:
            [1] Monte Carlo simulation of light transport in 3D heterogeneous
                tissues (mcxyz.c)
                Jacques S L, Li T
                http://omlc.org/software/mc/mcxyz/index.html 
                http://omlc.org/software/mc/mcxyz/mcxyzOct10_2014.zip
                (acesssed 07.08.2016)

            [2] Coupling 3D Monte Carlo light transport in optically 
                heterogeneous tissues to photoacoustic signal generation
                Jacques S L
                Photoacoustics 2 (2014) 137

        """

        def fetchData(file):
            with open(file) as f:
                return f.readlines()

        def fetchFluenceRate(file,roiSize):
            """read fluence rate  [J/m^2 per J delivered] from .bin file."""
            return np.fromfile(file,dtype = ct.c_float).reshape(roiSize)

        def fetchTissueType(file,roiSize):
            """read tissue type from .bin file."""
            return np.fromfile(file,dtype = ct.c_int8).reshape(roiSize)

        def fetchParameters(fName):
            """read simulation parameters from .mci file."""
            #data = open(fName,'r').readlines()
            data = fetchData(fName) 
            cast = lambda x,myType: myType(x.strip())
            roiShape = (cast(data[1],int), cast(data[2],int), cast(data[3],int))
            d3x = (cast(data[4],float), cast(data[5],float), cast(data[6],float))
            Nt = cast(data[21],int)
            muaDict = {i:cast(data[22+i*3],float) for i in range(Nt)}
            muaFunc = np.vectorize(lambda i: muaDict[i], otypes=[np.float]) 
            return roiShape, d3x, muaFunc

        # READ SIMULATION PARAMETERS FROM .mci FILE ---------------------------
        ((Nx,Ny,Nz), (dx,dy,dz), mua) = fetchParameters(fNameBase+'_H.mci')

        # READ FLUENCE RATE [J/m^2 per J delivered] FROM _F.bin FILE ----------
        Fxyz = fetchFluenceRate(fNameBase+'_F.bin',(Nz,Ny,Nx))
        
        # READ TISSUE TYPE FROM _T.bin FILE =----------------------------------
        Txyz = fetchTissueType(fNameBase+'_T.bin',(Nz,Ny,Nx))

        # COMPUTE ABSORBED ENERGY DENSITY [J/m^3 per J delivered ] ------------
        Hxyz = np.multiply(Fxyz,mua(Txyz))

        # SET COORDINATE AXES -------------------------------------------------
        x = np.linspace(0.5*dx, (Nx+0.5)*dx, Nx, endpoint=False)
        y = np.linspace(0.5*dy, (Ny+0.5)*dy, Ny, endpoint=False)
        z = np.linspace(0.5*dz, (Nz+0.5)*dz, Nz, endpoint=False)
        
        return (x, y, z), Hxyz


# EOF: mciBin.py
