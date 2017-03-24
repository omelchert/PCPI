import sys
import numpy as np


def OAVolterra2ndType_Karabutov(p0,wD,dt):
        """optoacoustic Volterra equation of 2nd type 

        compute optoacoustic (OA) waveform using the linear Volterra integral 
        equation of 2nd type discussed by Ref. [1]

        Args:
            p0 [numpy array, ndim=1]: initial pressure profile
            wD [float]: characteristic frequency of diffraction process
            dt [float]: time increment
        
        Returns:
            pD [numpy array, ndim=1]: excess pressure at detection point

        Refs:
            [1] Time-resolved laser optoacoustic tomography of inhomogeneous
                media
                Karabutov, A. A. and Podymova, N. B. and Letokhov, V. S.
                Appl. Phys. B, 63 (1996) 545
        """
        # INITIALIZATION 
        pD        = np.zeros(p0.size)           
        K0        = 1.0                         
        K1        = np.exp(-wD*dt)              
        K1_K0     = np.exp(-wD*dt)              

        # RECURRENCE RELATION implementing EQ. (40) of REF. [1] 
        I     = 0               
        pD[0] = p0[0]           
        for i in range(1,p0.size):
            I     = I*K1_K0 + 0.5*wD*dt*(K1*p0[i-1] + K0*p0[i])
            pD[i] = p0[i] - I
        return pD


def OAVolterra1stType_Sigrist(p0,wD,wa,dt):
        """optoacoustic Volterra equation of 1st type 

        compute optoacoustic (OA) waveform using the linear Volterra integral 
        equation of 1st type discussed by Refs. [1,2]

        Args:
            p0 [numpy array, ndim=1]: initial pressure profile
            wD [float]: characteristic frequency of diffraction process
            wa [float]: characteristic frequency of OA signal spectrum
            dt [float]: time increment
        
        Returns:
            pDf [numpy array, ndim=1]: excess pressure for free surface

        Refs:
            [1] Diffraction characteristics of laser-induced acoustic waves 
                in liquids
                Terzic, M. and Sigrist, M. W.
                J. Appl. Phys. 56 (1984) 93

            [2] Laser generation of acoustic waves in liquids and gases
                Sigrist, M. W.
                J. Appl. Phys. 60 (1986) R83
        """
        # INITIALIZE PROPAGATION KERNEL 
        K = np.exp(-wD*np.arange(p0.size)*dt)
        # EXCESS PRESSURE FOR RIGID BOUNDARY, EQ. (3) OF REF. [1]
        pDr = wa*np.convolve(p0,K,mode='full')*dt
        # EXCESS PRESSURE FOR FREE BOUNDARY, EQ. (4) of REF. [1]
        pDf = np.gradient(pDr)/wa/dt
        return pDf


def main():
        cmd = sys.argv
        # DOMAIN PARAMETERS ###################################################
        zMax = 10.              # DOMAIN: max z-Range  
        N = 6000                # DOMAIN: number of meshpoints
        c0 = 1.5                # DOMAIN: sound velocity [mm/mus]
        # ABSORBING LAYER PARAMETERS ##########################################
        w = 0.3                 # LAYER: absorbing layer width [mm]
        mu = float(cmd[1])      # LAYER: absorption coeff [mm^-1]
        # LASER AND DETECTOR PARAMETERS #######################################
        a0 = float(cmd[2])      # ISP: Karabutov spot radius
        d = a0/np.sqrt(2)       # ISP: diameter of laser spot
        zD = 1.2                # DET: detection distance
        # DERIVED PARAMETERS ##################################################
        wa = c0*mu              # frequency of OA signal spectrum 
        lac = 2.*np.pi/mu       # acoustic wavelength [due to Karabutov]  
        D = zD*lac/d/d          # diffraction parameter
        wD = c0*mu*D/2/np.pi    # related parameter  

        (z,dz) = np.linspace(0,zMax,N,endpoint=False,retstep=True)
        layerMask = np.logical_and(z>zD,z<zD+w)
        tau, dtau = (z - zD)/c0, dz/c0

        mua = np.zeros(N); mua[layerMask] = mu
        p0 = mua*np.exp(-np.cumsum(mua)*dz)/mu


        pD_K = OAVolterra2ndType_Karabutov(p0,wD,dtau)

        pD_S = OAVolterra1stType_Sigrist(p0,wD,wa,dtau)

        print "# (tau) (mu) (p0) (pD Karabutov) (pD Sigrist)"
        for i in range(p0.size):
            print c0*tau[i], mua[i], p0[i], pD_K[i], pD_S[i]


main()
