""" FILE: finitePulse.py

Module implementing functions to transform the solution obtained for a 
delta-pulse to a Laser pulse having a Gaussian profile with finite duration.
See, e.g. 

"Photoacoustic waves excited in liquids by fiber-transmitted laser pulses",
Paltauf, G. and Schmidt-Kloiber, H. and Franz, M.
J. Acoust. Soc. Am. 104 (1998) 890

and 

"A One-Dimensional Solution of the Photoacoustic Wave Equation and its
Relationship with Optical Absorption",
Cywiak, D. et al.
Int. J. Thermophys. 34 (2013) 1473

Author: O. Melchert
Date:   22.12.2015
        25.06.2016 (last changes)
"""
import numpy as np


def pulseKernelFlat(tAxis,tp):
        """Flat temporal irradiation profile

        Implements Laser pulse kernel with flat-top temporal shape
        
        \param tAxis time domain
        \param tp width of flat-top temporal profile
        """
        dt = tAxis[1]-tAxis[0]
        n  = int(tp/dt)/2
        k  = np.zeros(tAxis.size) 
        k[tAxis.size/2 - n: tAxis.size/2 + n] = 1
        return k


def pulseKernelGauss(tAxis,tp):
        """Gaussian temporal irradiation profile

        Implements Laser pulse kernel with gaussian temporal shape. Therein
        the standard deviation of the irradiation profile translates to 
        standard discrete gaussian via sigma = tp / sqrt(8)
        
        \param tAxis time domain
        \param tp width of Gaussian temporal profile
        """
        return np.exp(-4.*(tAxis-tAxis[int(0.5*tAxis.shape[0])])**2/tp**2)


def convolve(pDelta,kernel):
        """Finite-time Laser pulse

        Convolves solution obtained for delta-pulse to solution for finite
        time irradiation profile (specified by kernel)

        \param pDelta solution for delta-pulse
        \parem kernel finite-time irradiation profile
        """
        return np.convolve(kernel,pDelta,mode='same')/kernel.sum()


def convolveGauss(pDelta,(tAxis,tp)):
        return convolve(pDelta, pulseKernelGauss(tAxis,tp))


def convolveFlat(pDelta,(tAxis,tp)):
        return convolve(pDelta, pulseKernelFlat(tAxis,tp))


# FINITE-PULSE TEST SUITE: 
if __name__=="__main__":

        # SIMULATION PARAMETER ------------------------------------------------
        zMin   = 0.             # medium initial point
        zMax   = 1.             # medium final point
        N      = 1000
        
        # SET BEAM AXIS -------------------------------------------------------
        (tAxis,dt) = np.linspace(zMin,zMax,N, retstep=True)

        s = np.zeros(N) 
        s[N/2]=1
        s[N/4]=0.4
        tp = 100.*dt

        r  = convolve(s, pulseKernelGauss(tAxis,tp))

        for i in range(N):
                print i,r[i]

# EOF: finitePulse.py
