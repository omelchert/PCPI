"""FILE: customPolarSolver.pyx

Numpy based cython implementation of a Poisson-type integral equation that
might be used to obtain a solution of the inhomogeneous wave equation for a
fixed point in coordinate space as function of time. The solution is derived
from the corresponding Greens function in free space. See e.g. Sects II., III.
of [DeanBen:2012]: 


Refs:
    [1] Landau, L. D. and Lifshitz, E. M.,
        Hydrodynamik (4th Ed.),
        Akademie-Verlag (1981, Berlin)

    [2] Accurate Model-Based Reconstruction Algorithm for Three-Dimensional
        Optoacoustic Tomography, 
        Dean-Ben, X.L. and Buehler, A. and Ntziachristos, V. and Razansky, D. 
        IEEE Transaction on Medical Imaging, 31 (2012) 11922

author: O. Melchert
date: 2016/06/23
"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import  numpy as np
cimport numpy as np
cimport cython

# prefer C function with small overhead 
cdef extern from "math.h":
        double sqrt(double value) 

# WARNING: the following features were turned off to yield speed-up
# (i) bound-checking for array indices
# (ii) support for negative array indices 
@cython.boundscheck(False)
@cython.wraparound(False)
def polarPoissonIntegralSolver(
        np.ndarray[double, ndim=1] r,\
        np.ndarray[double, ndim=1] z,\
        np.ndarray[double, ndim=2] p0rz,\
        double rD,\
        double zD,\
        int Nphi
        ):
        """
        Polar coordinate based poisson integral solver 
        
        fast solver for the optoacoustic Poisson integral based on a 
        representation of the volumetric energy density within the source 
        volume in cylindrical polar coordinates

        \param[in]  rD    axial deflection of detector position from beam axis
        \param[in]  zD    z-coordinate of detector (zD<0: backward mode) 
        \param[in]  r     r-axis (ndim=1)
        \param[in]  z     z-axis (ndim=1)
        \param[in]  p0rz  initial acoustic stress profile (ndim=2)
        \param[in]  Nphi  interpolation points for azimuthal angle
        \param[out] tau   retarded signal depth
        \param[out] p     exess pressure at detection point 
        """
        # DECLARATION ---------------------------------------------------------
        cdef int iTauMax, i, j, k, binId
        cdef double dz, RR, d, dV
        cdef np.ndarray[double, ndim=1] tau, cosPhi, I
        
        # INITIALIZATION ------------------------------------------------------
        z = z - zD
        dz = z[1]-z[0]
        iTauMax = int(1.5*sqrt(np.max(np.abs(z))**2+(np.max(r)+rD)**2)/dz + 1)
        tau = np.linspace(0.,iTauMax*dz,iTauMax,endpoint=False)
        cosPhi = np.cos(np.linspace(0,np.pi,Nphi,endpoint=False))
        I = np.zeros(tau.size)
        dV = (z[1]-z[0])*(r[1]-r[0])*np.pi/Nphi

        # INTEGRATION OVER COMPUTATIONAL DOMAIN -------------------------------
        for i in range(z.size):
            for k in range(r.size):
                for j in range(Nphi):
                    RR = rD*rD + r[k]*r[k] + 2.0*rD*r[k]*cosPhi[j]
                    d = sqrt(RR + z[i]*z[i])
                    binId = int(d/dz)
                    I[binId] += 2.0*r[k]*p0rz[k,i]*dV/d
                
        return tau, I

# EOF: customPolarSolver.pyx 
