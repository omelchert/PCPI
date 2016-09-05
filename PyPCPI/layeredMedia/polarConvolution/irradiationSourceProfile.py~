""" FILE: irradiationSourceProfile.py

Module implementing different common irradiation source profiles 

Refs:
    [1] Pulsed optoacoustic characterization of layered media
        Paltauf, G. and Schmidt-Kloiber, H.
        Journal of Applied Physics, 88 (2000) 1624

    [2] Simulation study of melanoma detection in human skin tissues by 
        laser-generated surface acoustic waves 
        Chen, K. et al. 
        Journal of Biomedical Optics, 19 (2014) 077007 

    [3] Detection, numerical simulation and approximate inversion of 
        optoacoustic signals generated in multi-layered PVA hydrogel based 
        tissue phantoms
        Blumenr\"other, E. and Melchert, O. and Wollweber, M. and Roth, B. 
        (not published) arXiv:1605.05657.

    [4] 3-D Volume Reconstruction of Skin Lesions for Melanin and Blood Volume 
        Estimation and Lesion Severity Analysis. 
        D'Alessandro, B. and Dhawan, A. P. 
        IEEE Transactions on Medical Imaging, 31 (2012) 2083

"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import numpy as np 


def beamProfNormalization(r,f,P=1.0):
        """Beam profile normalization.

        Normalization of radially symmetric beam profile through solving 
        N = 2 pi \int_0_rMax r f(r) dr 
        and returning the normalization factor P/N.

        Args:
            r (numpy array, ndim=1): Equispaced 1D grid for radial coordinate.
            f (numpy array, ndim=1): Radial beam profile.
            P (float): Desired laser beam intensity.

        Returns:
            norm (float): Beam profile normalization factor.
        """
        return P/(2*np.pi*np.trapz(r*f,r))

def topHat(r,A0=0.2):
        """Top-hat beam profile as used in Ref. [1].
        
        Args:
            r (numpy array, ndim=1): Equispaced 1D grid for radial coordinate.
            A0 (float): Top hat width (default: 0.2).

        Returns: 
            f (numpy array, ndim=1): Top-hat beam profile.
        """
        return r<A0

def Gaussian(r,(R0,A0)=(0.,0.2)):
        """Gaussian beam profile as used in Ref. [1,2].
        
        Args:
            r (numpy array, ndim=1): Equispaced 1D grid for radial coordinate.
            A0 (float): 1/e2-width of Gaussian (default: 0.2).

        Returns:
            f (numpy array, ndim=1): Gaussian beam profile.
        """
        return np.exp(-2*(r-R0)**2/A0/A0)

def flatTop(r, (R0,D0)=(0.5,0.05)):
        """Flat-top beam profile as used in Ref. [1,3].
        
        Args:
            r (numpy array, ndim=1): Equispaced 1D grid for radial coordinate.
            R0 (float): Top-hat width (default: 0.5).
            D0 (float): 1/e2-width of Gaussian (default: 0.05).

        Returns:
            f (numpy array, ndim=1): Flat-top beam profile.
        """
        condList = [r<R0, r>=R0]
        funcList = [lambda r: 1, 
                    lambda r: np.exp(-2*(r-R0)**2/D0**2)]
        return np.piecewise(r,condList,funcList)

def flatTopDonut(r, (R0,R1,D0,D1)=(0.5,0.7,0.1,0.1)):
        """Flat-top donut beam profile as used in Ref. [4].
        
        Args:
            r (numpy array, ndim=1): Equispaced 1D grid for radial coordinate.
            R0 (float): Inner radius of radial profile (default: 0.5).
            R1 (float): Outer radius of radial profile (default: 0.7).
            D0 (float): 1/e2-width of inner Gaussian (default: 0.1).
            D1 (float): 1/e2-width of outer Gaussian (default: 0.1).

        Returns:
            f (numpy array, ndim=1): Flat-top donut beam profile.
        """
        condList = [r<=R0, (r>R0) & (r<R1), r>=R1]
        funcList = [lambda r: np.exp(-2*(r-R0)**2/D0**2), 
                    lambda r: 1., 
                    lambda r: np.exp(-2*(r-R1)**2/D1**2)]
        return np.piecewise(r,condList,funcList)

# EOF: irradiationSouceProfile.py
