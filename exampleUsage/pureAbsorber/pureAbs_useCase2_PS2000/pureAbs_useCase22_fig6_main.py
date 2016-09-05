''' FILE: pureAbs_useCase22_fig6_main.py

Script to model Source volume of the laboratory setup underlying Fig. 6 in
Refs. [1] and to compute distribution of absorbed energy and OA signal at
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
import PyPCPI.voxelizedMedia.dataIO as io 
import PyPCPI.voxelizedMedia.pureAbsorber.sourceVolume as sv 
import PyPCPI.voxelizedMedia.pureAbsorber.irradiationSourceProfile as isp 
import PyPCPI.voxelizedMedia.poissonIntegralSolver as ps

def modelSourceVolume():
        """Model source volume

        Source volume model of the laboratory situation described in 
        Fig. 6 in Ref. [1], Section IV.A on near field measurements.

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
            data displayed in Fig. 6 of Ref. [1].

        Refs:
            [1] Pulsed optoacoustic characterization of layered media
                Paltauf, G and Schmidt-Kloiber, H 
                J. Appl. Phys. 88 (2000) 1624
        """
        # SOURCE VOLUME PARS --------------------------------------------------
        xMax, Nx = 0.3, 1000    # bdry, meshpts: x-axis
        yMax, Ny = 0.3, 1000    # bdry, meshpts: y-axis
        zMax, Nz = 0.07, 100     # bdry, meshpts: z-axis

        # ABSORBING LAYER PARS ------------------------------------------------
        z01, dz1, ma1 = 0.0, 0.017, 25.      # absorbing layer 1 
        z02, dz2, ma2 = 0.017, 0.2-0.017, 100.      # absorbing layer 2 

        # FLAT TOP BEAM PROFILE  PARS -----------------------------------------
        x0, y0 = xMax/2, yMax/2 # x,y pos of symmetry axis
        fta, ftR = 0.45/2, 8.  # radius, edge width
        f0 = 1.0 # incident fluence (J/cm^2)

        # SET OPTICAL PROPERTIES OF SOURCE VOLUME -----------------------------
        (x,y,z), roi = sv.setROI((xMax,yMax,zMax), (Nx,Ny,Nz))
        sv.addAbsorbingLayer(((x,y,z),roi),z01,dz1,ma1)
        sv.addAbsorbingLayer(((x,y,z),roi),z02,dz2,ma2)

        # CHOOSE ISP AND PROPAGATE BEAM ---------------------------------------
        iProf = f0*isp.flatTop((x, y), (x0, y0), fta, ftR)
        sv.propagateBeam((z, roi), iProf) 
        # roi in units (J/cm^3 = 10^6 J/m^3). Scale output to (J/m^3). 
        return (x,y,z), (x0,y0,z01), roi*10**6


def main():
        # SIMULATION PARAMETERS -----------------------------------------------
        G = 1.0    # Grueneisen parameter
        c0 = 150000. # speed of sound (cm/s)
        zD = -0.1    # detector location 

        # OBTAIN ABSORBED ENERGY WITHIN SOURCE VOLUME -------------------------
        (x,y,z), (x0,y0,z0), Wxyz = modelSourceVolume()
        
        # COMPUTE OA EXCESS PRESSURE SIGNAL -----------------------------------
        t, p = ps.pressure(((x,y,z),G*Wxyz),(x0,y0,z0+zD),c0)

        # LIST RESULTS --------------------------------------------------------
        tau = t + zD/c0 # retarded time
        print "# (t in s) (c0 tau in cm) (p(tau) in au)"
        for i in range(tau.size):
               print t[i], c0*tau[i], p[i]


main()
# EOF: pureAbs_useCase22_fig6_main.py 
