''' FILE: pureAbs_useCase2_main.py

Script to model Source volume of the laboratory setup underlying Figs. 8(a,b)
in Refs. [1] and to compute distribution of absorbed energy and OA signal at
detector location.

Refs:
    [1] Pulsed optoacoustic characterization of layered media
        Paltauf, G and Schmidt-Kloiber, H 
        J. Appl. Phys. 88 (2000) 1624

AUTHOR: O. Melchert
DATE: 01.09.2016
'''

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import sys; sys.path.append('../../../')
import numpy as np
import PyPCPI.voxelizedMedia.dataIO as io 
import PyPCPI.voxelizedMedia.pureAbsorber.sourceVolume as sv 
import PyPCPI.voxelizedMedia.pureAbsorber.irradiationSourceProfile as isp 
import PyPCPI.voxelizedMedia.poissonIntegralSolver as ps



def modelSourceVolume(ma,ftd):
        """Model source volume

        Source volume model of the laboratory situation underlying Figs. 8(a,b) 
        in Ref. [1], Section IV.B on far field measurements.

        Args:
            ma (float): absorption coefficient in (cm^-1).
            ftd (float): diameter of flat-top ISP in (cm)

        Returns:
            (x,y,z) (3-tuple, numpy array, ndim=1): 1D coordinate grids.
            (x0,y0,z0) (2-tuple, floats): x and y center position of ISP 
                symmetry axis and z position of absorbing layer.
            roi (numpy array, ndim=3): computational region of interest 
                containing absorbed energy density in units (J/m^3).

        Note: 
            Used units are: mu (cm^-1) and F0 (J/cm^2), hence after beam
            propagation one yields roi (J/cm^3 = 10^6 J/m^3). Output is
            adjusted to roi in units (J/m^3).

            The parameters used for the flat-top beam profile and the
            irradiated fluence are adjusted to reproduce the experimental 
            data displayed in Fig. 8 of Ref. [1].

        Refs:
            [1] Pulsed optoacoustic characterization of layered media
                Paltauf, G and Schmidt-Kloiber, H 
                J. Appl. Phys. 88 (2000) 1624
        """
        # SOURCE VOLUME PARS --------------------------------------------------
        xMax, Nx = 0.4, 600    # bdry, meshpts: x-axis
        yMax, Ny = 0.4, 600    # bdry, meshpts: y-axis
        zMax, Nz = 0.1, 100     # bdry, meshpts: z-axis

        # ABSORBING LAYER PARS ------------------------------------------------
        z0, dz = 0.0, 0.03      # start, width of absorbing layer        

        # FLAT TOP BEAM PROFILE  PARS -----------------------------------------
        x0, y0 = xMax/2, yMax/2 # x,y pos of symmetry axis
        fta, ftR = ftd/2, 6.  # radius, edge width
        f0 = 1.0 # incident fluence (J/cm^2)

        # SET OPTICAL PROPERTIES OF SOURCE VOLUME -----------------------------
        (x,y,z), roi = sv.setROI((xMax,yMax,zMax), (Nx,Ny,Nz))
        sv.addAbsorbingLayer(((x,y,z),roi),z0,dz,ma)

        # CHOOSE ISP AND PROPAGATE BEAM ---------------------------------------
        iProf = f0*isp.flatTop((x, y), (x0, y0), fta, ftR)
        sv.propagateBeam((z, roi), iProf) 
        # roi in units (J/cm^3 = 10^6 J/m^3). Scale output to (J/m^3). 
        return (x,y,z), (x0,y0,z0), roi*10**6


def main():
        # SIMULATION PARAMETERS -----------------------------------------------
        G = 1.0    # Grueneisen parameter
        c0 = 150000. # speed of sound (cm/s)
        zD = -0.12    # detector location 
        ftd = float(sys.argv[1]) # flat-top diameter (cm)
        ma = float(sys.argv[2]) # absorption coefficient (cm^-1)


        # OBTAIN ABSORBED ENERGY WITHIN SOURCE VOLUME -------------------------
        (x,y,z), (x0,y0,z0), Wxyz = modelSourceVolume(ma,ftd)
        
        # COMPUTE OA EXCESS PRESSURE SIGNAL -----------------------------------
        t, p = ps.pressure(((x,y,z),G*Wxyz),(x0,y0,z0+zD),c0)
        norm = 1./max(p)

        # LIST RESULTS --------------------------------------------------------
        tau = t + zD/c0 # retarded time
        print "# (t in s) (c0 tau in cm) (p(tau) in MPa) (p(tau) in arb. units)"
        for i in range(tau.size):
               print t[i], c0*tau[i], p[i], p[i]*norm


main()
# EOF: pureAbs_useCase2_main.py 
