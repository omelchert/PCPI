""" FILE: useCase13_matplotlib_fullFigure.py

Exemplary script, reproducing the results shown in Fig. 4 and
Fig. 5 using matplotlib.

Figure is setup to be 7 inches wide (600 dpi) and have fontsize 
of 12 pt.

AUTHOR: O. Melchert
DATE: 11.01.2017
"""
import sys; sys.path.append('../../../')
import matplotlib.pyplot as plt
import PyPCPI.voxelizedMedia.dataIO as io 
import PyPCPI.voxelizedMedia.poissonIntegralSolver as ps
from mpl_toolkits.axes_grid1 import make_axes_locatable

def main():
        # SIMUMATION PARAMETERS -----------------------------------------------
        iPath = "../inputData/skinvessel" 
        G = 0.138 # Grueneisen parameter
        c0 = 150000. # speed of sound (cm/s)

        # VOLUMETRIC ENERGY DENSITY PER VOXEL --------------------------------- 
        ((x,y,z), Wxyz) = io.mcxyzio.getEnergyDensity(iPath)


        # GENERAL FIGURE SETTINGS ---------------------------------------------
        plt.rc('font',family='serif',size=12)
        plt.figure(figsize=(7.,9.0))

        # FIG (a): x-scan -----------------------------------------------------
        ax11= plt.subplot2grid((3,2), (0,0))
        ax11.set_aspect('equal')
        ax11.contourf(y,z,Wxyz[:,:,x.size/2]/1000,100,lw=0.1)
        c1 = ax11.contourf(y,z,Wxyz[:,:,x.size/2]/1000,100)
        divider1 = make_axes_locatable(ax11)
        cax11 = divider1.append_axes("right", "5%", pad="3%")
        c1Bar = plt.colorbar(c1,cax=cax11)
        c1Bar.set_label('$W(x_0,y,z)$ (kJ/m$^2$)',fontsize=10)
        c1Bar.ax.tick_params(labelsize=10)
        ax11.set_ylabel('$z$  (cm)')
        ax11.set_xlabel('$y$  (cm)')
        ax11.text(0.,1.05,'(a)',transform=ax11.transAxes,fontsize=12)


        # FIG (b): y-scan -----------------------------------------------------
        ax12= plt.subplot2grid((3,2), (0,1))
        ax12.set_aspect('equal')
        ax12.contour(x,z,Wxyz[:,y.size/2,:]/1000,100,lw=0.1)
        c2 = ax12.contourf(x,z,Wxyz[:,y.size/2,:]/1000,100)
        divider2 = make_axes_locatable(ax12)
        cax12 = divider2.append_axes("right", "5%", pad="3%")
        c2Bar = plt.colorbar(c2,cax=cax12)
        c2Bar.set_label('$W(x,y_0,z)$ (kJ/m$^2$)',fontsize=10)
        c2Bar.ax.tick_params(labelsize=10)
        ax12.set_ylabel('$z$  (cm)')
        ax12.set_xlabel('$x$  (cm)')
        ax12.text(0.,1.05,'(b)',transform=ax12.transAxes,fontsize=12)

        # FIG (c): OFF-AXIS OA SIGNALS ----------------------------------------
        ax21 = plt.subplot2grid((3,2), (1,0), colspan=2)
        ax21.text(0.5,0.75,'$x_D$=0.05 cm\n$z_D$=-0.20 cm',transform=ax21.transAxes,fontsize=10)
        xD=0.05
        zD=-0.20
        for yD in [0.01, 0.03, 0.05]:
            t, p = ps.pressure(((x,y,z),G*Wxyz/1000),(xD,yD,zD),c0)
            ax21.plot(c0*t+zD,p, label='$y_D$=%3.2lf cm'%(yD))
        ax21.set_xlim([-0.002,0.06])
        ax21.set_ylabel('$p$ (kPa per J)')
        ax21.set_xlabel('$c \\tau$ (cm)')
        ax21.axhline(y=0,color='black', ls='--')
        ax21.legend(loc='upper right', prop={'size':10}, frameon=False)
        ax21.text(0.,1.05,'(c)',transform=ax21.transAxes,fontsize=12)

        # FIG (d): ON-AXIS OA SIGNALS -----------------------------------------
        ax31 = plt.subplot2grid((3,2), (2,0), colspan=2)
        ax31.text(0.5,0.75,'$x_D$=0.05 cm\n$y_D$=0.05 cm',transform=ax31.transAxes,fontsize=10)
        xD=0.05
        yD=0.05
        for zD in [-0.1, -0.4, -1.6]:
            t, p = ps.pressure(((x,y,z),G*Wxyz/1000),(xD,yD,zD),c0)
            ax31.plot(c0*t+zD,p, label='$z_D$=%3.2lf cm'%(zD))
        ax31.set_xlim([-0.002,0.06])
        ax31.set_ylabel('$p$ (kPa per J)')
        ax31.set_xlabel('$c \\tau$ (cm)')
        ax31.axhline(y=0,color='black', ls='--')
        ax31.legend(loc='upper right', prop={'size':10}, frameon=False)
        ax31.text(0.,1.05,'(d)',transform=ax31.transAxes,fontsize=12)

        plt.tight_layout()
        plt.savefig('useCase_mcxyz_abcd.eps',dpi=600,format='eps')

main()
# EOF: useCase13_matplotlib_fullFigure.py 
