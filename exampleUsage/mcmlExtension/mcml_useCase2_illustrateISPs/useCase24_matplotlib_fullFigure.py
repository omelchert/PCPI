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
        plt.figure(figsize=(7.,7.))

        # FIG (b): FLAT-TOP ISP
        ax12= plt.subplot2grid((2,2), (0,0))
        ax12.set_aspect('equal')
        r,z,W_ft = io.npz.readROI("./dataUseCase21/Wrz_flatTop.npz")
        ax12.contourf(z,r,W_ft,100,lw=0.1)
        c2 = ax12.contourf(z,r,W_ft,100)
        divider2 = make_axes_locatable(ax12)
        cax12 = divider2.append_axes("right", "5%", pad="3%")
        c2Bar = plt.colorbar(c2,cax=cax12)
        c2Bar.set_label('$W(r,z)$ (J/m$^2$)',fontsize=10)
        c2Bar.ax.tick_params(labelsize=10)
        ax12.set_ylabel('$r$  (cm)')
        ax12.set_xlabel('$z$  (cm)')
        ax12.text(0.,1.05,'(a)',transform=ax12.transAxes,fontsize=12)

        # FIG (c): DONUT ISP
        ax13= plt.subplot2grid((2,2), (0,1))
        ax13.set_aspect('equal')
        r,z,W_d  = io.npz.readROI("./dataUseCase21/Wrz_Donut.npz")
        ax13.contour(z,r,W_d,100,lw=0.1)
        c3 = ax13.contourf(z,r,W_d,100)
        divider3 = make_axes_locatable(ax13)
        cax13 = divider3.append_axes("right", "5%", pad="3%")
        c3Bar = plt.colorbar(c3,cax=cax13)
        c3Bar.set_label('$W(r,z)$ (J/m$^2$)', fontsize=10)
        c3Bar.ax.tick_params(labelsize=10)
        ax13.set_ylabel('$r$  (cm)')
        ax13.set_xlabel('$z$  (cm)')
        ax13.text(0.,1.05,'(b)',transform=ax13.transAxes,fontsize=12)

        # FIG (d): ON-AXIS OA SIGNALS FOR DONUT ISP 
        ax21 = plt.subplot2grid((2,2), (1,0), colspan=2)
        for zD in [-0.5,-1.0,-2.0,-4.0,-8.0]:
            t,p_d = ps.pressure((r,z,G*W_d),(rD,zD),c0,1) 
            p_d = pp.convolveGauss(p_d,(t,tp))
            ax21.plot(c0*t+zD,p_d, label='$z_D$=%3.2lf cm'%(zD))
        ax21.set_xlim([-0.005,0.18])
        ax21.set_ylabel('$p$ (Pa per J)')
        ax21.set_xlabel('$c \\tau$ (cm)')
        ax21.axhline(y=0,color='black', ls='--')
        ax21.legend(loc='upper right', prop={'size':10}, frameon=False)
        ax21.text(0.,1.05,'(c)',transform=ax21.transAxes,fontsize=12)

        plt.tight_layout(pad=2.)
        #plt.show()
        plt.savefig('useCase_mcml_abc.eps',dpi=600,format='eps')


figure()
# EOF: useCase24_matplotlib_fullFigure.py
