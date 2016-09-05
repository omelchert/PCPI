""" FILE: irradiationSourceProfile.py

Module implementing transverse irradiation source profiles for beam propagation  
through source volume.

"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import numpy as np

def beamProfNormalization((x,y),isp,P=1.0):
        """Beam profile normalization.

        Normalization of radially symmetric beam profile through solving 
        N = 2 pi \int_0_rMax r f(r) dr 
        and returning the normalization factor P/N.

        Args:
            x (numpy array, ndim=1): Equispaced 1D grid for x coordinate.
            y (numpy array, ndim=1): Equispaced 1D grid for x coordinate.
            isp (numpy array, ndim=2): 2D irradiation source profile. 
            P (float): Desired laser beam intensity.

        Returns:
            norm (float): Beam profile normalization factor.
        """
        I = np.zeros(y.size)
        for i in range(y.size):
            I[i] = np.trapz(isp[i,:],x)
        norm = np.trapz(I,y)
        return P/norm

def planeWave((x,y)):
        """plane wave irradiation source profile.

        Args:
            x (numpy array, ndim=1): Equispaced 1D grid for x coordinate.
            y (numpy array, ndim=1): Equispaced 1D grid for x coordinate.

        Returns:
            isp (numpy array, ndim=2) 2D plane wave beam profile. 
        """
        return np.ones((len(y),len(x))) 


def Gaussian((x, y), (x0,y0)=(0.5,0.5), sigma=0.1):
        """Gaussian irradiation source profile.

        Args:
            x (numpy array, ndim=1): Equispaced 1D grid for x coordinate.
            y (numpy array, ndim=1): Equispaced 1D grid for x coordinate.
            (x0,y0) (2-tuple, floats): x and y center location of Gaussian
                profile (default: (x0,y0)=(0.5,0.5)).
            sigma (float): exp(-2) width of Gaussian profile (default: 0.1).

        Returns:
            isp (numpy array, ndim=2) 2D Gaussian beam profile. 
        """
        yy,xx = np.meshgrid(y, x, indexing='ij')
        return np.exp(-( (xx-x0)**2 + (yy-y0)**2)/(2*sigma*sigma)) 


def topHat((x, y), (x0,y0)=(0.5,0.5), a=0.1):
        """Top-hat irradiation source profile

        Args:
            x (numpy array, ndim=1): Equispaced 1D grid for x coordinate.
            y (numpy array, ndim=1): Equispaced 1D grid for x coordinate.
            (x0,y0) (2-tuple floats): x and y center of top hat location
                (default: (x0,y0)=(0.5,0.5)).
            a (float): width of top hat profile (default: 0.1).

        Returns:
            isp (numpy array, ndim=2) 2D top hat beam profile. 
        """
        yy,xx = np.meshgrid(y, x, indexing='ij')
        return np.sqrt((xx-x0)**2 + (yy-y0)**2) <= a 

def flatTop((x,y),(x0,y0)=(0.5,0.5),a=0.1,adRatio=2.0):
        """Smoothed top-hat irradiation source profile

        Args:
            x (numpy array, ndim=1): Equispaced 1D grid for x coordinate.
            y (numpy array, ndim=1): Equispaced 1D grid for x coordinate.
            (x0,y0) (2-tuple, floats): x and y center location of Gaussian 
                profile (default: (x0,y0)=(0.5,0.5)).
            a (float): width of top hat profile (default: 0.1).
            adRatio (float): ratio of top-hat width to 1/e-depth of smoothed
                region (default: adRatio=2. [Elias old Laser profile]).

        Returns:
            isp (numpy array, ndim=2) 2D flat top beam profile. 
        """
        yy,xx = np.meshgrid(y, x, indexing='ij')
        mask = (xx-x0)**2 + (yy-y0)**2 <= a*a
        isp  = np.exp(-((np.sqrt((xx-x0)**2+(yy-y0)**2) - a)**2)/(a/adRatio)**2) 
        isp[mask] = 1.
        return isp 

#EOF: irradiationSourceProfile.py
