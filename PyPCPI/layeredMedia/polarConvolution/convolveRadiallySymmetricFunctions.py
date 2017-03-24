""" FILE: convolveRadiallySymmetricFunctions.py

Module implementing 2D convolution of two radially symmetric functions.

Refs:
    [1] Operational and convolution properties of two-dimensional Fourier 
        transforms in polar coordinates 
        Baddour, N.
        J. Opt. Soc. Am. A 26 (2009) 1767-1777 

"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import scipy
import scipy.special as scs
import numpy as np 
import irradiationSourceProfile
import hankelTransform


def convolve(r,rho,G0,F0):
        """Compute 2D convolution of two radially symmetric functions.

        Implements full 2D convolution of two radial functions based on 
        their zeroth order Fourier coefficients `G0` and `F0` as detailed
        in Eq. (84) of Ref. [1].

        Args:
            r (numpy array, ndim=1): equispaced 1D grid for radial coordinate.
            rho (numpy array, ndim=1): equispaced 1D grid for complementary 
                coordinate relative to `r`.
            G0 (numpy array, ndim=1): Zeroth order Fourier coefficient of 
                radially symmetric function specifying material response to 
                pencil beam in Fourier space.
            F0 (numpy array, ndim=1): Zeroth order Fourier coefficient of 
                radially symmetric function specifying beam profile in Fourier 
                space.
        
        Returns:
            h (numpy array, ndim=1): 2D convolution of two radially symmetric 
                functions based on their Fourier representation.

        Notes:
            The 2D convolution integral given by Eq. (84) of see Ref. [1] 
            relies on the zeroth order coefficients G0(rho)=2\pi H0(g(r)) and 
            F0(rho)=2\pi H0(f(r)) of the Fourier Transforms of f and g, NOT 
            their plain Hankel transforms. For radially symmetric functions 
            only this zeroth order Fourier coefficient is nonzero. Thus, the 
            integrahd formulated in terms of their Hankel transforms introduces 
            an additional factor of (2\pi)^2 in Eq. (84).

        Refs:
            [1] Operational and convolution properties of two-dimensional 
                Fourier transforms in polar coordinates 
                Baddour, N.
                J. Opt. Soc. Am. A 26 (2009) 1767-1777 
        """
        rhoG0F0J0 = rho*G0*F0*scs.j0(rho*r[:,np.newaxis])
        return 2*np.pi*np.trapz(rhoG0F0J0, rho, axis=1)

def convolveROI(r,z,grz,fr,P=1.0):
        """Compute response to spatially extended irradiation source profile.

        Computes full 2D convolution of pencil beam response with extended
        beam profile at each z-coordinate of region of interest (ROI).

        Args:
            r (numpy array, ndim=1): equispaced 1D grid for radial coordinate.
            z (numpy array, ndim=1): equispaced 1D grid for depth coordinate.
            grz (numpy array, ndim=2): response to infinitely narrow beam.
            fr (numpy array, ndim=1): radial beam profile.
            P (float): beam intensity (default: 1.0).

        Returns:
            Wrz (numpy array, ndim=2): response to spatially extended 
                irradiation source profile.

        See Also:
            convolve: 2D convolution of two radially symmetric functions.

        Notes: 
            Due to the spatial invariance of the region of interest, the 
            material response to an extended photon beam with given radial 
            profile is obtained as convolution of two radially symmetric 
            functions for given z-coordinate, see Eq. (84) of Ref. [1]. For the
            integration, a numpy native trapezoidal rule is used.

        Refs:
            [1] Operational and convolution properties of two-dimensional 
                Fourier transforms in polar coordinates 
                Baddour, N.
                J. Opt. Soc. Am. A 26 (2009) 1767-1777 
        """
        Wrz  = np.zeros(grz.shape)

        # Normalize beam profile using function from isp module
        f0 = irradiationSourceProfile.beamProfNormalization(r,fr,P)
        # Hankel transform of beam profile will be used repeatedly so 
        # precompute it here for time-efficiency 
        rho,F0 = hankelTransform.zeroOrderHankelTrafo(r,fr)

        for iz in range(z.size):
            # Hankel transform of response to pencil beam
            rho,G0 = hankelTransform.zeroOrderHankelTrafo(r,grz[:,iz])
            # Gibbs-phenomenon workaround: trimm off all negative valued 
            # contributions to volumetric energy density at given depth
            hr = convolve(r,rho,G0,F0)
            hr = np.abs(hr)
            #hr = np.maximum(hr,np.zeros(hr.size))
            # Normalize contribution to volumetric energy density
            Wrz[:,iz] = f0*hr
        return Wrz

# EOF: convolveRadiallySymmetricFunctions.py
