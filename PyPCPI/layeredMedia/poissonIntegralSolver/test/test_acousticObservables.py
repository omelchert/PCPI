""" FILE: test_acousticObservables.py

Unittest module for acousticObservables.py
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
from acousticObservables import pressure, velocityPotential

class AcousticObservablesTestCase(unittest.TestCase):
        """Unit test for acousticObservables.py.
        
        Implements unit test for region of interest (ROI) with constant
        absorbed volumetric energy density.
        Inherits from unittest.TestCase.
        """

        def setUp(self):
            """Sets up unit test prerequisits.
            
            Attributes:
                r: numpy array containing equidistant radial gridpoints. 
                z: numpy array containing equidistant z-axis gridpoints.
                Wrz: numpy array containing absorbed volumetric energy density.
            """
            N      = 200
            rMin   = 0.0
            rMax   = N*20e-6
            zMin   = 0.0
            zMax   = N*20e-6

            self.r = np.linspace(rMin,rMax,N,retstep=False,endpoint=False)
            self.z = np.linspace(zMin,zMax,N,retstep=False,endpoint=False)
            self.Wrz = np.ones((self.r.size,self.z.size))*10**6


        def tearDown(self):
            """Deletes all attributes set by method setUp()."""
            del self.r
            del self.z
            del self.Wrz


        def test_pressureSignalGeneration(self):
            """Perform unit test for medium with uniform energy deposition.
            
            The pressure signal for a pointlike detector in the middle of 
            a medium with constant absorbed volumetric energy density is 
            computed. The integration over the region of interest is carried
            out in cylindrical polar coordinates. The pressure signal in the 
            range [0,zMax), where the signal is not affected by the spatial 
            limits of the ROI should be equal to 0.138 MPa. The unit test 
            averages the signal over the range [0,zMax), normalizes by 
            0.138 MPa and checks for unity up to the second decimal.

            Note that this unit test is based on an exemplary simulation 
            summarized in Fig 4(a,b) of 

            "Coupling 3D Monte carlo light transport in optically 
            heterogeneous tissues to photoacoustic signal generation",
            Jacques, S. L.,
            Photoacoustics 2 (2014) 137-142.
            """
            c0  = 1500.
            rho = 1000.
            G   = 0.138
            rD  = 0.
            zD  = 0.5*(self.z[self.z.size/2]+self.z[1+self.z.size/2])

            t,p = pressure(
                (self.r,self.z,G*self.Wrz),
                (rD,zD),
                c0,
                1)

            t,I = velocityPotential(
                (self.r,self.z,self.Wrz),
                (rD,zD),
                (G,rho,c0),
                1)

            iMax = 50
            ratio = np.sum(p[:iMax])/iMax/G/10**6

            self.assertAlmostEqual(ratio, 1.,1)


if __name__ == "__main__":
        unittest.main()

# EOF: test_acousticObservables.py
