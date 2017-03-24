""" FILE: test_irradiationSource.py

Unittest module for irradiationSource.py
"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import sys; sys.path.append("../")
import unittest
import numpy as np
import irradiationSourceProfile as isp


class BeamProfileNormalizationTestCase(unittest.TestCase):
        """Unit test for irradiationSource.py.
        
        Implements normalization unit tests for the irradiation source module. 
        Unittest considers Gaussian beam profile.
        Inherits from unittest.TestCase.
        """

        def setUp(self):
            """Sets up unit test prerequisits.
            
            Attributes:
                r: Numpy array containing mesh points for test function.
                f: Irradiation source profile.
                R0: 1/e2-width of gaussian beam profile.
                P: Desired beam intensity.
            """
            rMin    = 0.0
            rMax    = 1.
            N       = 100000
            self.R0 = 0.2
            self.P  = 2.0
            self.r  = np.linspace(rMin,rMax,N,retstep=False,endpoint=False)
            self.f  = isp.Gaussian(self.r,(0.,self.R0))

        def tearDown(self):
            """Deletes all attributes set by method setUp()."""
            del self.r
            del self.f
            del self.R0
            del self.P

        def test_normalization(self):
            """Perform beam profile normalization unit test."""
            Inum = isp.beamProfNormalization(self.r,self.f,self.P)
            Iex = 2.*self.P/(np.pi*self.R0**2)
            self.assertAlmostEqual(Inum,Iex,places=6)


if __name__=="__main__":
        unittest.main()

# EOF: test_irradiationSource.py
