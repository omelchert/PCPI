""" FILE: useCase24_matplotlib_fullFigure.py

Exemplary script, reproducing the results shown in Fig. 2 and
Fig. 3 using matplotlib

AUTHOR: O. Melchert
DATE: 09.01.2017
"""
import sys; sys.path.append('../../../')
import numpy as np
import matplotlib.pyplot as plt
import PyPCPI.layeredMedia.dataIO as io 
import PyPCPI.layeredMedia.poissonIntegralSolver as ps
import PyPCPI.signalPostProcessing.finitePulse as pp
from mpl_toolkits.axes_grid1 import make_axes_locatable


def figure():

        # PARAMETERS FOR OA SIGNAL GENERATION
        c0 = 150000.0           # sonic velocity
        G = 1.0                 # Grueneissen parameter
        rD = 0.0                # detector location relative to beam axis
        tp = 2.*10**(-8)        # gaussian pulse width (s) 

        # GENERAL FIGURE SETTINGS ---------------------------------------------
        plt.rc('font',family='serif',size=12)
        plt.figure(figsize=(7.,6.0))

        # FIG (a): GAUSSIAN ISP
        ax11= plt.subplot2grid((2,3), (0,0))
        r,z,W_g  = io.npz.readROI("./dataUseCase21/Wrz_Gaussian.npz")
        ax11.set_aspect('equal')
        c1 = ax11.contourf(z,r,W_g,100)
        divider1 = make_axes_locatable(ax11)
        cax11 = divider1.append_axes("right", "5%", pad="3%")
        c1Bar = plt.colorbar(c1,cax=cax11)
        c1Bar.set_label('$W(r,z)$ (J/m$^2$)',fontsize=10)
        c1Bar.ax.tick_params(labelsize=10)
        ax11.set_ylabel('$r$  (cm)')
        ax11.set_xlabel('$z$  (cm)')


        # FIG (b): FLAT-TOP ISP
        ax12= plt.subplot2grid((2,3), (0,1))
        r,z,W_ft = io.npz.readROI("./dataUseCase21/Wrz_flatTop.npz")
        c2 = ax12.contourf(z,r,W_ft,100)
        divider2 = make_axes_locatable(ax12)
        cax12 = divider2.append_axes("right", "5%", pad="3%")
        c2Bar = plt.colorbar(c2,cax=cax12)
        c2Bar.set_label('$W(r,z)$ (J/m$^2$)',fontsize=10)
        c2Bar.ax.tick_params(labelsize=10)
        ax12.set_ylabel('$r$  (cm)')
        ax12.set_xlabel('$z$  (cm)')

        # FIG (c): DONUT ISP
        ax13= plt.subplot2grid((2,3), (0,2))
        r,z,W_d  = io.npz.readROI("./dataUseCase21/Wrz_Donut.npz")
        c3 = ax13.contourf(z,r,W_d,100)
        divider3 = make_axes_locatable(ax13)
        cax13 = divider3.append_axes("right", "5%", pad="3%")
        c3Bar = plt.colorbar(c3,cax=cax13)
        c3Bar.set_label('$W(r,z)$ (J/m$^2$)', fontsize=10)
        c3Bar.ax.tick_params(labelsize=10)
        ax13.set_ylabel('$r$  (cm)')
        ax13.set_xlabel('$z$  (cm)')

        # FIG (d): ON-AXIS OA SIGNALS FOR DONUT ISP 
        ax21 = plt.subplot2grid((2,3), (1,0), colspan=3)
        for zD in [-0.5,-1.0,-2.0,-4.0,-8.0]:
            t,p_d = ps.pressure((r,z,G*W_d),(rD,zD),c0,1) 
            p_d = pp.convolveGauss(p_d,(t,tp))
            ax21.plot(c0*t+zD,p_d, label='$z_D$=%3.2lf cm'%(zD))
        ax21.set_xlim([-0.005,0.18])
        ax21.set_ylabel('$p$ (Pa per J)')
        ax21.set_xlabel('$c \\tau$ (cm)')
        ax21.axhline(y=0,color='black', ls='--')
        ax21.legend(loc='upper right', prop={'size':10}, frameon=False)


        plt.tight_layout()
        plt.show()


figure()
# EOF: useCase24_matplotlib_fullFigure.py
