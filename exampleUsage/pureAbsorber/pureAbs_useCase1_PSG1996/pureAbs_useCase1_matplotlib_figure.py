''' FILE: pureAbs_useCase1_matplotlib_figure.py

Script to model Source volume of the laboratory setup described in Refs. [1,2] 
and to compute distribution of absorbed energy and OA signal at detector 
location for different values of the ISP top-hat width. The resulting data
curves are visualized using matplotlib.

Refs:
    [1] Light distribution measurement sin absorbing materials by 
        optical detection of laser-induced stress waves
        Paltauf, G and Schmidt-Kloiber, H and Guss, H
        Appl. Phys. Lett. 69 (1996) 1526


    [2] Measurement of laser-induced acoustic waves with a calibrated
        optical transducer
        Paltauf, G and Schmidt-Kloiber, H
        J. Appl. Phys. 82 (1997) 1525 

AUTHOR: O. Melchert
DATE: 10.01.2017
'''

__authors__   = "O. Melchert"
__copyright__ = "(c) 2017, Hannover Centre for Optical Technologies"
__license__   = "3-clause BSD License"
__contact__   = "oliver.melchert@hot.uni-hannover.de"

import sys; sys.path.append('../../../')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import PyPCPI.voxelizedMedia.dataIO as io 
import PyPCPI.voxelizedMedia.pureAbsorber.sourceVolume as sv 
import PyPCPI.voxelizedMedia.pureAbsorber.irradiationSourceProfile as isp 
import PyPCPI.voxelizedMedia.poissonIntegralSolver as ps

def modelSourceVolume(fta):
        """Model source volume

        Source volume model of the laboratory situation described in 
        Refs. [1,2].

        Args:
            None

        Returns:
            (x,y,z) (3-tuple, numpy array, ndim=1): 1D coordinate grids.
            (x0,y0,z0) (2-tuple, floats): x and y center position of ISP 
                symmetry axis and z position of absorbing layer.
            roi (numpy array, ndim=3): computational region of interest 
                containing absorbed energy density in units (J/m^3).

        Note: 
            Used units are: mu (cm^-1) and F0 (J/cm^2), hence after beam
            propagation one yields roi (J/cm^3 = 10^6 J/m^3). Output is
            adjusted to roi in units (J/m^3).

            The parameters used for the flat-top beam profile and the
            irradiated fluence are taken from Ref. [2].

        Refs:
            [1] Light distribution measurement sin absorbing materials by 
                optical detection of laser-induced stress waves
                Paltauf, G and Schmidt-Kloiber, H and Guss, H
                Appl. Phys. Lett. 69 (1996) 1526

            [2] Measurement of laser-induced acoustic waves with a calibrated
                optical transducer
                Paltauf, G and Schmidt-Kloiber, H
                J. Appl. Phys. 82 (1997) 1525 
        """
        # SOURCE VOLUME PARS --------------------------------------------------
        xMax, Nx = 0.6, 600     # bdry, meshpts: x-axis
        yMax, Ny = 0.6, 600     # bdry, meshpts: y-axis
        zMax, Nz = 0.6, 200     # bdry, meshpts: z-axis

        # ABSORBING LAYER PARS ------------------------------------------------
        z0, dz = 0.0, 0.6       # start, width of absorbing layer        
        mu = 24.0               # absorption coeffiecient

        # FLAT TOP BEAM PROFILE  PARS -----------------------------------------
        x0, y0 = xMax/2, yMax/2 # x,y pos of symmetry axis
        fta, ftr = fta, 4.5  # radius, radius/edge width ratio
        f0 = 0.18 # incident fluence (J/cm^2)

        # SET OPTICAL PROPERTIES OF SOURCE VOLUME -----------------------------
        (x,y,z), roi = sv.setROI((xMax,yMax,zMax), (Nx,Ny,Nz))
        sv.addAbsorbingLayer(((x,y,z),roi),z0,dz,mu)

        # CHOOSE ISP AND PROPAGATE BEAM ---------------------------------------
        iProf = f0*isp.flatTop((x, y), (x0, y0), fta, ftr)
        sv.propagateBeam((z, roi), iProf) 
        # roi in units (J/cm^3 = 10^6 J/m^3). Scale output to (J/m^3). 
        return (x,y,z), (x0,y0,z0), roi*10**6


def figure():
        # SIMULATION PARAMETERS -----------------------------------------------
        G = 0.11    # Grueneisen parameter
        c0 = 150000. # speed of sound (cm/s)
        zD = -0.1    # detector location 

        # GENERAL FIGURE SETTINGS ---------------------------------------------
        plt.rc('font',family='serif',size=12)
        plt.figure(figsize=(7.,7.))

        #io.gpl.writeROI(x,z,Wxyz[:,y.size/2,:],'./dataUseCase1/Wxy0z.prof')

        # FIG (a): FLAT-TOP ISP
        ax12= plt.subplot2grid((2,2), (0,0))
        ax12.set_aspect('equal')
        (x,y,z), (x0,y0,z0), Wxyz = modelSourceVolume(0.025)
        ax12.contour(x,z,Wxyz[:,y.size/2,:]/10**6,100,lw=0.1)
        c2 = ax12.contourf(x,z,Wxyz[:,y.size/2,:]/10**6,100)
        divider2 = make_axes_locatable(ax12)
        cax12 = divider2.append_axes("right", "5%", pad="3%")
        c2Bar = plt.colorbar(c2,cax=cax12)
        c2Bar.set_label('$W(x,y_0,z)$ (MJ/m$^2$)',fontsize=10)
        c2Bar.ax.tick_params(labelsize=10)
        ax12.set_ylabel('$z$  (cm)')
        ax12.set_xlabel('$x$  (cm)')
        ax12.text(0.,1.05,'(a)',transform=ax12.transAxes,fontsize=12)
        ax12.text(0.05,0.9,'$a = 0.025$ cm',transform=ax12.transAxes,fontsize=12,color='white')

        # FIG (b): FLAT-TOP ISP 
        ax13= plt.subplot2grid((2,2), (0,1))
        ax13.set_aspect('equal')
        (x,y,z), (x0,y0,z0), Wxyz = modelSourceVolume(0.194)
        ax13.contour(x,z,Wxyz[:,y.size/2,:]/10**6,100,lw=0.1)
        c3 = ax13.contourf(x,z,Wxyz[:,y.size/2,:]/10**6,100)
        divider3 = make_axes_locatable(ax13)
        cax13 = divider3.append_axes("right", "5%", pad="3%")
        c3Bar = plt.colorbar(c3,cax=cax13)
        c3Bar.set_label('$W(x,y_0,z)$ (MJ/m$^2$)',fontsize=10)
        c3Bar.ax.tick_params(labelsize=10)
        ax13.set_ylabel('$z$  (cm)')
        ax13.set_xlabel('$x$  (cm)')
        ax13.text(0.,1.05,'(b)',transform=ax13.transAxes,fontsize=12)
        ax13.text(0.05,0.9,'$a = 0.194$ cm',transform=ax13.transAxes,fontsize=12,color='white')

        ax21 = plt.subplot2grid((2,2), (1,0), colspan=2)
        for fta in [0.025, 0.05, 0.075, 0.1, 0.194]:
        #for fta in [0.025, 0.194]:
            # OBTAIN ABSORBED ENERGY WITHIN SOURCE VOLUME ---------------------
            (x,y,z), (x0,y0,z0), Wxyz = modelSourceVolume(fta)
            # COMPUTE OA EXCESS PRESSURE SIGNAL -------------------------------
            t, p = ps.pressure(((x,y,z),G*Wxyz),(x0,y0,z0+zD),c0)
            # ADD OA SIGNAL TO PLOT -------------------------------------------
            ax21.plot((t+zD/c0)*10**6,p/10**5, label='$a$=%3.2lf cm'%(fta))
        ax21.set_xlim([-0.01*10**6/c0,0.23*10**6/c0])
        ax21.set_ylabel('$p$ (bar)')
        ax21.set_xlabel('$\\tau$ ($\\mu$s)')
        ax21.axhline(y=0,color='black', ls='--')
        ax21.legend(loc='upper right', prop={'size':10}, frameon=False)
        ax21.text(0.,1.05,'(c)',transform=ax21.transAxes,fontsize=12)
        ax21.text(0.5,0.92,'$z_D = -0.1$ cm',transform=ax21.transAxes)

        plt.tight_layout()
        plt.savefig('useCase_pureAbs_abc.eps',dpi=600,format='eps')


figure()
# EOF: pureAbs_useCase1_matplotlib_figure.py 
