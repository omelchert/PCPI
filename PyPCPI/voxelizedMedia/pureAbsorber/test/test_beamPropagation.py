""" FILE: test_beamPropagation.py

Unittest module for sourceVolume.py and irradiationSourceprofile.py.
"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import sys; sys.path.append("../")
import unittest
import numpy as np
import sourceVolume as sv
import irradiationSourceProfile as isp

class BeamPropagationTestCase(unittest.TestCase):
        """Unit test for sourceVolume.py and irradiationSource.py.
        
        Implements unit test for region of interest (ROI) with known 
        absorption profile.
        Inherits from unittest.TestCase.
        """

        def setUp(self):
            """Sets up unit test prerequisits.
            
            Attributes:
                rMax (3-tuple, floats): x, y, z coordinates upper boundaries. 
                rN (3-tuple, ints): Nx, Ny, Nz number of equispaced samples.
                par (3-tuple, floats): z0 (starting location), w (width), and 
                    ma (constant absorption coefficient) of absorbing layer. 
                f (function): function implementinc expected exact absorption
                    profile for comparison.
            """
            N      = 200
            xMax   = 1.0 
            yMax   = 1.0
            zMax   = 1.0
            z0 = 0.2
            w = 0.6
            mu = 10.0

            self.rMax = (xMax,yMax,zMax)
            self.rN = (N,N,N)
            self.par = (z0,w,mu)
       
            def func(r):
                condList = [r<=z0, (r>z0) & (r<z0+w), r>=z0+w]
                funcList = [lambda r: 0.,
                            lambda r: mu*np.exp(-mu*(r-z0)),
                            lambda r: 0.]
                return np.piecewise(r,condList,funcList)

            self.f = lambda x: func(x)


        def tearDown(self):
            """Deletes all attributes set by method setUp()."""
            del self.rMax
            del self.rN
            del self.par
            del self.f


        def test_beamPropagation(self):
            """Perform unit test for beam propagation.

            Check computed depth profile of absorption profile vs expected
            (exact) profile for single layered source volume.
            
            """
            (z0, w, mu) = self.par
            (x, y, z), roi = sv.setROI(self.rMax, self.rN)
            sv.addAbsorbingLayer(((x, y, z), roi), z0, w, mu)
            sv.propagateBeam((z, roi), isp.planeWave((x, y)))

            # 1D slice of roi
            roi_z = roi[:,0,0]
            # expected profile
            f_z = self.f(z)
            # sum of squared residuals
            err = np.sum((roi_z-f_z)**2)

            self.assertAlmostEqual(err, 0.,6)


if __name__ == "__main__":
        unittest.main()

# EOF: test_beamPropagation.py
