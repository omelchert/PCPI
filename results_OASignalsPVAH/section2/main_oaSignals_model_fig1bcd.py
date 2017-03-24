''' FILE: main_oaSignals_model_fig1a.py 


AUTHOR: O. Melchert
DATE: 20.01.2017
'''

import sys; sys.path.append('../../')
import numpy as np
import PyPCPI.voxelizedMedia.dataIO as io 
import PyPCPI.voxelizedMedia.pureAbsorber.sourceVolume as sv 
import PyPCPI.voxelizedMedia.pureAbsorber.irradiationSourceProfile as isp 
import PyPCPI.voxelizedMedia.poissonIntegralSolver as ps



def modelSourceVolume():
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
        #xMax, Nx = 1.0, 900    # bdry, meshpts: x-axis
        #yMax, Ny = 1.0, 900    # bdry, meshpts: y-axis
        #zMax, Nz = 1.0, 300     # bdry, meshpts: z-axis
        xMax, Nx = 1.0, 900    # bdry, meshpts: x-axis
        yMax, Ny = 1.0, 900    # bdry, meshpts: y-axis
        zMax, Nz = 1.0, 263     # bdry, meshpts: z-axis

        # ABSORBING LAYER PARS ------------------------------------------------
        z0, dz = 0.0, 0.1      # start, width of absorbing layer        

        # FLAT TOP BEAM PROFILE  PARS -----------------------------------------
        x0, y0 = xMax/2, yMax/2 # x,y pos of symmetry axis
        fta, ftR = 0.1, 2.  # radius, edge width
        f0 = 1.0 # incident fluence (J/cm^2)

        # SET OPTICAL PROPERTIES OF SOURCE VOLUME -----------------------------
        ma = 24.0 # absorption coeficcient
        (x,y,z), roi = sv.setROI((xMax,yMax,zMax), (Nx,Ny,Nz))
        sv.addAbsorbingLayer(((x,y,z),roi),z0,dz,ma)

        # CHOOSE ISP AND PROPAGATE BEAM ---------------------------------------
        iProf = f0*isp.flatTop((x, y), (x0, y0), fta, ftR)
        sv.propagateBeam((z, roi), iProf) 
        # roi in units (J/cm^3 = 10^6 J/m^3). Scale output to (J/m^3). 
        return (x,y,z), (x0,y0,z0), roi*10**6


def detectorPositions((x0,z0),RD):
        pos = [(np.degrees(phi),x0+RD*np.cos(phi),z0-RD*np.sin(phi))   for phi in np.radians([90.0,180.0,270.0])]
        return pos


def sinogram():
        # SIMULATION PARAMETERS -----------------------------------------------
        G = 1.0    # Grueneisen parameter
        c0 = 150000. # speed of sound (cm/s)
        R = 0.5001

        # OBTAIN ABSORBED ENERGY WITHIN SOURCE VOLUME -------------------------
        (x,y,z), (x0,y0,z0), Wxyz = modelSourceVolume()
        detPos = detectorPositions((x0,z0+0.2),R)

        for (phi,xD,zD) in detPos:
       
                sys.stderr.write("# phi=%lf \n"%(phi))
                print "# phi, xD, zD = ", phi, xD, zD
                print "# D = ", abs(zD-0.1)*2/24/0.15/0.15
                # COMPUTE OA EXCESS PRESSURE SIGNAL -----------------------------------
                t, p = ps.pressure(((x,y,z),G*Wxyz),(xD,y0,z0+zD),c0)

                # LIST RESULTS --------------------------------------------------------
                tau = t + zD/c0 # retarded time
                print "# (t in s) (c0 tau in cm) (p(tau) in au)"
                norm = 1./max(p)
                for i in range(tau.size):
                       print phi, c0*t[i], c0*tau[i], p[i], p[i]*norm
                print

def main():
        # SIMULATION PARAMETERS -----------------------------------------------
        G = 1.0    # Grueneisen parameter
        c0 = 150000. # speed of sound (cm/s)
        zD = -0.5    # detector location 


        # OBTAIN ABSORBED ENERGY WITHIN SOURCE VOLUME -------------------------
        (x,y,z), (x0,y0,z0), Wxyz = modelSourceVolume()
        
        # COMPUTE OA EXCESS PRESSURE SIGNAL -----------------------------------
        t, p = ps.pressure(((x,y,z),G*Wxyz),(x0,y0,z0+zD),c0)

        # LIST RESULTS --------------------------------------------------------
        tau = t + zD/c0 # retarded time
        print "# (t in s) (c0 tau in cm) (p(tau) in au)"
        for i in range(tau.size):
               print t[i], c0*tau[i], p[i]


#main()
sinogram()
# EOF: main_oaSignals_model_fig1a.py 
