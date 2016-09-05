""" FILE: hankelTransform.py

Module implementing functions to compute zero order hankel transform. 

Refs:
    [1] Operational and convolution properties of two-dimensional Fourier 
        transforms in polar coordinates 
        Baddour, N.
        J. Opt. Soc. Am. A 26 (2009) 1767-1777 
  
    [2] Algorithms to Numerically Evaluate the Hankel Transform
        Cree, M. J. and Bones, P. J.
        Computers Math. Applic., 26 (1993) 1-12
  
    [3] Fast computation of zero order Hankel transform
        Gopalan, K. and Chen, C. S. 
        Journal of the Franklin Institute, 316 (1983) 317-326

"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import scipy
import scipy.special as scs
import numpy as np 

def zeroOrderHankelTrafo(r,fr):
        """Compute zero order Hankel transform.

        Compute zeroth order hankel transform using trapezoidal quadrature 
        rule for equally spaced samples as explained in section 4 of Ref. [2].

        Args:
            r (numpy array, ndim=1): equispaced 1D grid for radial coordinate.
            fr (numpy array, ndim=1): objective function.

        Returns:
            rho (numpy array, ndim=1): equi-spaced complementary grid.
            F0 (numpy array, ndim=1): zeroth order Hankel transform of 
                objective function.

        Notes:
            Note that the notation follows Refs. [1], [3]. The computational 
            efficiency is O(N^2).

        Refs:
            [1] Operational and convolution properties of two-dimensional 
                Fourier transforms in polar coordinates 
                Baddour, N.
                J. Opt. Soc. Am. A 26 (2009) 1767-1777 
          
            [2] Algorithms to Numerically Evaluate the Hankel Transform
                Cree, M. J. and Bones, P. J.
                Computers Math. Applic., 26 (1993) 1-12
          
            [3] Fast computation of zero order Hankel transform
                Gopalan, K. and Chen, C. S. 
                Journal of the Franklin Institute, 316 (1983) 317-326
        """
        dr    = r[1]-r[0]
        rho   = np.linspace(0,1./(2.*dr),r.size,endpoint=False)
        rfrJ0 = r*fr*scs.j0(r*rho[:,np.newaxis])
        return rho,np.trapz(rfrJ0,x=r,axis=1) 

# EOF: hankelTransform.py
