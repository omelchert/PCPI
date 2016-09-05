"""FILE: customCartesianSolverMcxyz.pyx

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
def cartPoissonIntegralSolver(
        np.ndarray[double, ndim=1] x,\
        np.ndarray[double, ndim=1] y,\
        np.ndarray[double, ndim=1] z,\
        np.ndarray[double, ndim=3] p0xyz,\
        double xD,\
        double yD,\
        double zD
        ):
        """cartesian coordinate based poisson integral solver. 
        
        fast solver for the optoacoustic Poisson integral based on a 
        representation of the volumetric energy density within the source 
        volume in cartesian coordinates

        Args:
            x (numpy array, ndim=1) x-coordinate of detector 
            y (numpy array, ndim=1) y-coordinate of detector 
            z (numpy array, ndim=1) z-coordinate of detector
            p0xyz (numpy array, ndim=3) initial acoustic stress profile 
            xD (double) x-coordinate of detector position 
            yD (double) y-coordinate of detector position 
            zD (double) z-coordinate of detector position (zD<0: backward mode) 

        Returns:
            tau (numpy array, ndim=1) measurement depth
            I (numpy array, ndim=1) Poisson Integral base of oa pressure
        """
        # DECLARATION ---------------------------------------------------------
        cdef int iTauMax, i, j, k, binId
        cdef double dz, d, dV, r, rr, dI
        cdef np.ndarray[double, ndim=1] tau, I
        
        # INITIALIZATION ------------------------------------------------------
        x = x - xD
        y = y - yD
        z = z - zD
        dz = z[1]-z[0]
        dV = (x[1]-x[0])*(y[1]-y[0])*dz 
        iTauMax = int((np.max(np.abs(z))+np.max(np.abs(x))+np.max(np.abs(y)))/dz)
        tau = np.linspace(0.,iTauMax*dz,iTauMax,endpoint=False)
        I = np.zeros(tau.size)
        
        # INTEGRATION OVER COMPUTATIONAL DOMAIN -------------------------------
        for i in range(x.size):
            for j in range(y.size):
                rr = x[i]*x[i] + y[j]*y[j]
                for k in range(z.size):
                    d = sqrt(rr + z[k]*z[k])
                    binId = int(d/dz)
                    dI = p0xyz[k,j,i]*dV/d
                    I[binId] += dI 
                
        return tau, I

# EOF: customCartesianSolverMcxyz.pyx 
