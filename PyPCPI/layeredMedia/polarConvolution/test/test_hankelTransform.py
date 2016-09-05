""" FILE: test_hankelTransform.py

Unittest module for hankelTransform.py
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
from hankelTransform import zeroOrderHankelTrafo 


def sombrero():
    a0 = 2.
    S  = lambda r: a0*scs.j1(a0*r)/r 
    H  = lambda r: r<a0 
    return S,H

def exp():
    a0  = 1.5
    S   = lambda r: np.exp(-a0*r)
    H   = lambda r: a0/(a0*a0+r*r)**1.5
    return S,H

def Gauss():
    a0 = 1./np.sqrt(np.pi)
    S  = lambda r: np.exp(-r*r/a0/a0/(2*np.pi)**2) 
    H  = lambda r: np.exp(-r*r/a0/a0)*2*np.pi     
    return S,H

def Bessel():
    a0 = 20.
    S  = lambda r: scs.j0(a0*r)
    H  = lambda r: np.pi*2./r if r==a0 else 0
    return S,H

def eRMS(Fn,Fx):
    """Compute root mean square error for input arrays.

    Implements root mean square error (eRMS) according to Eq (25) in Ref

    Algorithms to Numerically Evaluate the Hankel Transform
    Cree, M. J. and Bones P. J.
    Computers Math. Applic. 26 (1993) 1-12 

    Args:
        Fn: numpy array containing values of Hankel transform of test function.
        Fx: numpy array with same length as Fn containing values of exact
            Hankel transform of test function.

    """
    return np.sqrt(((Fn-Fx)**2).mean()/(Fx*Fx).mean())


class HankelTransformTestCase(unittest.TestCase):
        """Unit test for hankelTransform.py.
        
        Implements two unit tests for the hankel transform module. 
        Inherits from unittest.TestCase.
        """

        def setUp(self):
            """Sets up unit test prerequisits.
            
            Attributes:
                r: numpy array containing mesh points for test function.
                f: test function.
                Fx: exact Hankel transform of test function.
            """
            rMin   = 0.0
            rMax   = 10.
            N      = 1000
            self.r = np.linspace(rMin,rMax,N,retstep=False,endpoint=False)
            self.f, self.Fx = Gauss()

        def tearDown(self):
            """Deletes all attributes set by method setUp()."""
            del self.r
            del self.f
            del self.Fx

        def test_hankelTransform_FourierPair(self):
            """Perform unit test on Fourier pair."""
            rho,Fn = zeroOrderHankelTrafo(self.r,self.f(self.r))
            self.assertLessEqual( eRMS(Fn,self.Fx(rho)), 1.0)

        def test_hankelTransform_selfReciprocal(self): 
            """Perform unit test checking self reciprocality."""
            rho, Fn = zeroOrderHankelTrafo(self.r,self.f(self.r))
            r, finv = zeroOrderHankelTrafo(rho,Fn)
            self.assertLessEqual( eRMS(finv,self.f(r)), 1.0)

        def test_hankelTransform_Parseval(self):
            """Perform unit test based on Parseval relationship.

            Uses Parseval relation as detailed in Eq. (99) of Ref. [1] for 
            a radially symmetric function, i.e. a function for which the 
            Fourier series only contains the zeroth order n=0.
            """
            rho, F0 = zeroOrderHankelTrafo(self.r,self.f(self.r))
            # NOTE:
            # Hankel Transform of f: F0 = \int f(r) J0(r rho) rho drho
            # Fourier Transform of f: F = 2 \pi F0
            # 
            # LHS of Eq. (99) in Ref. [1]
            # \int |f|^2 d^3r = 2 \pi \int r f^2 dr 
            #     \equiv 2 \pi IrDomain
            IrDomain = np.trapz(self.r*self.f(self.r)*self.f(self.r),self.r)
            # RHS of Eq. (99) in Ref. [1]
            # (2 \pi)^-1 \int |F|^2 rho drho = 2 \pi \int |F0|^2 rho drho   
            #     \equiv 2 \pi IrhoDomain          
            IrhoDomain = np.trapz(rho*F0*F0,rho)
            self.assertAlmostEqual(IrDomain, IrhoDomain,1)
            

if __name__ == "__main__":
        unittest.main()

# EOF: test_hankelTransform.py
