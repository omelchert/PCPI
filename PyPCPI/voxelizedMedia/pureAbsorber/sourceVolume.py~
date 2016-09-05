""" FILE: sourceVolume.py

Module implementing functions to set up absorption coefficient values within 
source volume voxels and to propagate given transverse irradiation source
profile.

"""

__authors__   = "O. Melchert"
__copyright__ = "(c) 2016, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import numpy as np


def setROI((xMax,yMax,zMax)=(0.1,0.1,0.1),(Nx,Ny,Nz)=(100,100,100)):
        """Set region of region of interest.

        Args:
            (xMax, yMax, zMax) (3-tuple floats): upper boundary coordinates
                for x, y, z coordinate axes.
            (Nx, Ny, Nz) (3-tuple ints): number of samples for x, y, z axes. 

        Returns:
            (x,y,z) (3-tuple numpy arrays, ndim=1): Equispaced x, y, z 
                coordinate meshes.
            roi (numpy array, ndim=3): Empty source volume. 
        """
        x = np.linspace(0,xMax,Nx, endpoint=False)
        y = np.linspace(0,yMax,Ny, endpoint=False)
        z = np.linspace(0,zMax,Nz, endpoint=False)
        roi = np.zeros((Nz,Ny,Nx))
        return (x,y,z), roi

def addAbsorbingSphere(((x,y,z),roi),(x0,y0,z0),r,ma):
        """Add purely absorbing shpere to ROI.

        Args:
            (x,y,z) (3-tuple numpy arrays, ndim=1): Equispaced x, y, z 
                coordinate meshes.
            roi (numpy array, ndim=3): Source volume. 
            (x0, y0, z0) (3-tuple floats): center location of absorbing sphere. 
            r (float): radius of sphere.
            ma (float): absorption coefficient of sphere.

        Returns:
            none
        """
        zz,yy,xx = np.meshgrid(z,y,x, indexing='ij')
        mask = (xx-x0)**2 + (yy-y0)**2  + (zz-z0)**2 <= r*r 
        roi[mask] = ma 

def addAbsorbingCylinder(((x,y,z),roi),(x0,z0),r,ma):
        """Add purely absorbing cylinder to ROI.

        Cylinder extends along y direction.

        Args:
            (x,y,z) (3-tuple numpy arrays, ndim=1): Equispaced x, y, z 
                coordinate meshes.
            roi (numpy array, ndim=3): Source volume. 
            (x0, z0) (2-tuple floats): center location of absorbing cylinder. 
            r (float): radius of cylinder.
            ma (float): absorption coefficient of cylinder.

        Returns:
            none
        """
        zz,yy,xx = np.meshgrid(z,y,x, indexing='ij')
        mask = (xx-x0)**2 + (zz-z0)**2 <= r*r 
        roi[mask] = ma

def addAbsorbingLayer(((x,y,z),roi),z0,w,ma):
        """Add purely absorbing layer to ROI.

        Layer extends along x,y plane.

        Args:
            (x,y,z) (3-tuple numpy arrays, ndim=1): Equispaced x, y, z 
                coordinate meshes.
            roi (numpy array, ndim=3): Source volume. 
            z0 (float): initial boundary location of absorbing layer. 
            w (float): width of absorbing layer. 
            ma (float): absorption coefficient within layer.

        Returns:
            none
        """
        zz,yy,xx = np.meshgrid(z,y,x, indexing='ij')
        mask = (zz > z0) * (zz<z0+w) 
        roi[mask] = ma

def propagateBeam((z,roi),isp):
        """Propagate beam along z-directions through ROI.

        Args:
            z (numpy array, ndim=1): Equispaced z coordinate axis.
            roi (numpy array, ndim=3): Source volume. 
            isp (numpy array, ndim=2): irradiation source profile. 

        Returns:
            none
        """
        dz=z[1]-z[0]
        expFac = np.exp(-roi[0,:,:]*dz)
        roi[0,:,:] = isp[:,:]*roi[0,:,:]*expFac 
        for i in range(1,z.size):
                expFac *= np.exp(-roi[i,:,:]*dz)
                roi[i,:,:] = isp[:,:]*roi[i,:,:]*expFac


# EOF: sourceVolume.py
