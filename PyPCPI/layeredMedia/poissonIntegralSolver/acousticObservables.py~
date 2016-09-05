""" FILE: acousticObservables.py

Module implementing functions to compute acoustic observables. 

Refs:
    [1] Landau, L. D. and Lifshitz, E. M.,
        Hydrodynamik (4th Ed.),
        Akademie-Verlag (1981, Berlin)

"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import sys
import scipy
import scipy.special as scs
import numpy as np 
from poissonIntegral_cython import polarPoissonIntegralSolver


def pressure((r,z,p0rz),(rD,zD),c0=1.,Nphi=360):
        """Compute acoustic observables.

        Compute on-axis and off-axis excess pressure signals for given 
        material response to extended radially symmetric transverse 
        irradiation source profile as observed by a pointlike detector.

        Args:
            r (numpy array, ndim=1): equispaced 1D grid for radial coordinate.
            z (numpy array, ndim=1): equispaced 1D grid for z-coordinate.
            p0rz (numpy array, ndim=2): region of interest (ROI) containing 
                initial acoustic stress response to extended irradiation 
                source profile.
            rD (float): axial deflection of detector position from beam axis.
            zD (float): z-position of detector relative to first layer.
            Nphi (int): number of mesh-points for azimuthal coordinate, only
                important if rD is different from zero, i.e. for off axis 
                signals. 

        Returns:
            t (numpy array, ndim=1): equi-spaced complementary grid.
            p (numpy array, ndim=1): excess pressure at detector position as
                function of time.

        Notes:
            The excess pressure integral at a given field point over time sould
            yield zero, see chapter 68 of Ref. [1]. The integral is available
            on output as `Ip`.

            If on input Hrz = initial acoustic stress, then on output 
            p = optoacoustic excess pressure, if on input Hrz = absorbed 
            volumetric energy density, then on output p = optoacoustic 
            excess pressure/Grueneisenparameter

        Refs:
            [1] Landau, L. D. and Lifshitz, E. M.,
                Hydrodynamik (4th Ed.),
                Akademie-Verlag (1981, Berlin)

        """
        c0t, I = polarPoissonIntegralSolver(r, z, p0rz, rD, zD, Nphi)

        I = I/(4.*np.pi*c0)/(c0t[1]-c0t[0])
        p = np.gradient(I)*c0/(c0t[1]-c0t[0])

        return c0t/c0, p 


def velocityPotential((r,z,p0rz),(rD,zD),(Gamma,rho,c0),Nphi=360):
        """Compute velocity potential at detector field point.

        Compute on-axis and off-axis velocity potentials for given 
        material response to extended radially symmetric transverse 
        irradiation source profile as observed by a pointlike detector.

        Args:
            r (numpy array, ndim=1): equispaced 1D grid for radial coordinate.
            z (numpy array, ndim=1): equispaced 1D grid for z-coordinate.
            p0rz (numpy array, ndim=2): region of interest (ROI) containing 
                initial acoustic stress response to extended irradiation 
                source profile.
            rD (float): axial deflection of detector position from beam axis.
            zD (float): z-position of detector relative to first layer.
            Nphi (int): number of mesh-points for azimuthal coordinate, only
                important if rD is different from zero, i.e. for off axis 
                signals. 

        Returns:
            t (numpy array, ndim=1): equi-spaced time grid of detector signal.
            I (numpy array, ndim=1): velocity potential at detector position as
                function of time.

        """
        c0t, I = polarPoissonIntegralSolver(r, z, p0rz, rD, zD, Nphi)

        I = I/(4.*np.pi*c0)/(c0t[1]-c0t[0])
        phi = -Gamma/rho * I

        return c0t/c0, phi

# EOF: acousticObservables.py 
