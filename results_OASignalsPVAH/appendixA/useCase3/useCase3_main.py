''' FILE: useCase2_main.py

Script to model Source volume of the laboratory setup described in Refs. [1,2] 
and to compute distribution of absorbed energy and OA signal at detector 
location.

Refs:
    [1] Light distribution measurement sin absorbing materials by 
        optical detection of laser-induced stress waves
        Paltauf, G and Schmidt-Kloiber, H and Guss, H
        Appl. Phys. Lett. 69 (1996) 1526


    [2] Measurement of laser-induced acoustic waves with a calibrated
        optical transducer
        Paltauf, G and Schmidt-Kloiber, H
        J. Appl. Phys. 82 (1997) 1525 

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
import PyPCPI.signalPostProcessing.finitePulse as pp

def modelSourceVolume():
        """Model source volume

        Source volume model of the laboratory situation described in 
        Refs. [1,2].

        Args:
            None

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
            irradiated fluence are taken from Ref. [2].

        Refs:
            [1] Light distribution measurement sin absorbing materials by 
                optical detection of laser-induced stress waves
                Paltauf, G and Schmidt-Kloiber, H and Guss, H
                Appl. Phys. Lett. 69 (1996) 1526

            [2] Measurement of laser-induced acoustic waves with a calibrated
                optical transducer
                Paltauf, G and Schmidt-Kloiber, H
                J. Appl. Phys. 82 (1997) 1525 
        """
        # SOURCE VOLUME PARS --------------------------------------------------
        xMax, Nx = 10., 750 #800     # bdry, meshpts: x-axis
        yMax, Ny = 10., 750 #800     # bdry, meshpts: y-axis
        zMax, Nz = 1., 50 #200     # bdry, meshpts: z-axis

        # ABSORBING LAYER PARS ------------------------------------------------
        z0, dz = 0.0, 0.9       # start, width of absorbing layer        
        mu = 1.65               # absorption coeffiecient

        # FLAT TOP BEAM PROFILE  PARS -----------------------------------------
        x0, y0 = xMax/2, yMax/2 # x,y pos of symmetry axis
        ftr = 4.
        fta = float(sys.argv[1])/(1.+1./ftr) # radius, radius/edge width ratio
        f0 = 1. # incident fluence (J/cm^2)

        # SET OPTICAL PROPERTIES OF SOURCE VOLUME -----------------------------
        (x,y,z), roi = sv.setROI((xMax,yMax,zMax), (Nx,Ny,Nz))
        sv.addAbsorbingLayer(((x,y,z),roi),z0,dz,mu)

        # CHOOSE ISP AND PROPAGATE BEAM ---------------------------------------
        #iProf = f0*isp.topHat((x, y), (x0, y0), fta)
        iProf = f0*isp.flatTop((x, y), (x0, y0), fta, ftr)
        sv.propagateBeam((z, roi), iProf) 
        # roi in units (J/cm^3 = 10^6 J/m^3). Scale output to (J/m^3). 
        return (x,y,z), (x0,y0,z0), roi*10**6


def main():
        # SIMULATION PARAMETERS -----------------------------------------------
        G = 1.0    # Grueneisen parameter
        c0 = 1.5 # speed of sound (cm/s)
        zD = 3.150    # detector location 
        R0 = float(sys.argv[2])

        # OBTAIN ABSORBED ENERGY WITHIN SOURCE VOLUME -------------------------
        (x,y,z), (x0,y0,z0), Wxyz = modelSourceVolume()

        rDList = np.linspace(0.,R0,30,endpoint=True)

        pList = {}
        nMin = None
        for idx,rD in enumerate(rDList):
                # COMPUTE OA EXCESS PRESSURE SIGNAL -----------------------------------
                sys.stderr.write("# %d \n"%(idx))
                #t, p = ps.pressure(((x,y,z),G*Wxyz),(x0+rD,y0,zD),c0)
                t, p = ps.pressure_LambertianDetector(((x,y,z),G*Wxyz),(x0+rD,y0,zD),c0)
                if len(p)<nMin or nMin==None:
                        nMin=len(p)
                pList[idx] = (t[:nMin],p[:nMin])

        pSum = np.zeros(nMin)
        dr = rDList[1]-rDList[0]
        t0,p0 = pList[0]
        pSum += 0.5*rDList[0]*p0*2*np.pi 
        for i in range(1,len(pList)-1):
                ti,pi = pList[i]
                pSum += rDList[i]*pi*2*np.pi 
        ti,pi = pList[len(pList)-1]
        pSum += 0.5*rDList[-1]*pi*2*np.pi 

        #for i in range(1,len(pList)):
        #        ti,pi = pList[i]
        #        tim,pim = pList[i-1]
        #        ri = rDList[i]
        #        rim = rDList[i-1]
        #        pSum += 0.5*(ri-rim)*(ri*pi+rim*pim)*2*np.pi 

        t,p = pList[0]
        # LIST RESULTS --------------------------------------------------------
        tau = t + zD/c0 # retarded time
        dz = 0.1
        pAv = pp.convolveFlat(pSum,(c0*t,dz))
        norm = 1./max(pSum)
        print "# (t in s) (c0 tau in cm) (p(tau) in bar)"
        for i in range(tau.size):
               print t[i], c0*tau[i], pSum[i]*norm,pAv[i]*norm

        print "# sum(p) = ", sum(pSum)

def main2():
        # SIMULATION PARAMETERS -----------------------------------------------
        G = 1.0    # Grueneisen parameter
        c0 = 1.5 # speed of sound (cm/s)
        zD = 3.250    # detector location 

        rDList = np.linspace(0.,1.5,20,endpoint=True)

        # OBTAIN ABSORBED ENERGY WITHIN SOURCE VOLUME -------------------------
        (x,y,z), (x0,y0,z0), Wxyz = modelSourceVolume()

        for idx,rD in enumerate(rDList):
                # COMPUTE OA EXCESS PRESSURE SIGNAL -----------------------------------
                t, p = ps.pressure(((x,y,z),G*Wxyz),(x0+rD,y0,zD),c0)

                tau = t + zD/c0 # retarded time
                print "# (t in s) (c0 tau in cm) (p(tau) in bar)"
                for i in range(tau.size):
                       print rD, t[i], c0*tau[i], p[i]
                print

main()
# EOF: useCase2_main.py 
