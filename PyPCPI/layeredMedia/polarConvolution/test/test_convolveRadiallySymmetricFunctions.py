""" FILE: test_radiallySymmetricFunctions.py

Unittest module for radiallySymmetricFunctions.py
"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import sys; sys.path.append("../")
import unittest
import scipy
import scipy.special as scs
import numpy as np
import hankelTransform as ht 
import convolveRadiallySymmetricFunctions as crs


def eRMS(Fn,Fx):
    """Compute root mean square error for input arrays.

    Implements root mean square error (eRMS) according to Eq (25) in Ref. [1].

    Args:
        Fn: numpy array containing values of Hankel transform of test function.
        Fx: numpy array with same length as Fn containing values of exact
            Hankel transform of test function.

    Refs:
        [1] Algorithms to Numerically Evaluate the Hankel Transform
            Cree, M. J. and Bones P. J.
            Computers Math. Applic. 26 (1993) 1-12 
    """
    return np.sqrt(((Fn-Fx)**2).mean()/(Fx*Fx).mean())


class BeamProfilePolarConvolutionTestCase(unittest.TestCase):
        """Unit test for radiallySymmetricFunctions.py.
        
        Inherits from unittest.TestCase
        """

        def setUp(self):
            """Sets up unit test prerequisits
            
            Attributes:
                r: Equidistant mesh points for test function.
                rho: Equidistant mesh points of complementary coordinate.
                h: function in r domain. 
                G0: factor 1 of function in rho domain. 
                F0: factor 2 of function in rho domain. 
                H0: 2D fourier trafo of h.
            """
            rMin    = 0.01
            rMax    = 10.
            N       = 2000
            self.r  = np.linspace(rMin,rMax,N,retstep=False,endpoint=False)
            self.h  = lambda r: r*r*np.exp(-np.pi*r*r) 
            self.G0 = lambda rho: np.exp(-rho*rho/4/np.pi)/2/np.pi**2
            self.F0 = lambda rho: (1-rho*rho/4/np.pi)/2/np.pi
            self.rho, self.H0 = ht.zeroOrderHankelTrafo(self.r,self.h(self.r))


        def tearDown(self):
            """Deletes all attributes set by method setUp()."""
            del self.r
            del self.h
            del self.F0
            del self.G0
            del self.rho
            del self.H0


        def test_convolve(self):
            """Perform 2D polar convolution unit test."""
            hr = crs.convolve(self.r, self.rho, 
                              self.G0(self.rho), self.F0(self.rho)) 
            self.assertAlmostEqual(eRMS(hr,self.h(self.r)), 0.,3)


if __name__=="__main__":
        unittest.main()

# EOF: test_radiallySymmetricFunctions.py
