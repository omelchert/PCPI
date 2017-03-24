""" FILE: useCase15_matplotlib_oaSignal.py

AUTHOR: O. Melchert
DATE: 10.01.2017
"""
import sys; sys.path.append('../../../')
import matplotlib.pyplot as plt
import PyPCPI.layeredMedia.dataIO as io 
import PyPCPI.layeredMedia.poissonIntegralSolver as ps

def figure():
        fName = './dataUseCase11/Wrz_2Layer_mcml_ISPGauss_R00.15.npz'

        c0 = 1.0 # speed of sound in medium 
        G = 1.0 # efficiency of transforming absorbed energy to acoustic stress
        zD = -2.0 # z-location of pointlike detector 
        Nphi = 1 # number of azimutal mesh points considered for integration 
                 # (in case of on-axis signals Nphi=1 suffices)

        # READ VOLUMETRIC ENERGY DENSITY FROM COMPRESSED .npz FILE
        (r,z,Wrz) = io.npz.readROI(fName)

        # COMPUTE OA SIGNAL AT DETECTOR POSITION
        t,p = ps.pressure((r,z,G*Wrz),(0.0,zD),c0,Nphi) 
        
        # VISUALIZE USING MATPLOTLIB
        plt.rc('font',family='serif',size=14)
        plt.figure(figsize=(8.,6.), dpi=80)
        plt.plot(c0*t+zD,p, label='$z_D$=%3.2lf cm'%(zD))
        plt.xlim(-0.005,0.18)
        plt.ylabel('$p$ (Pa)', size=14)
        plt.xlabel('$c \\tau$ (cm)', size=14)
        plt.axhline(y=0,color='black', ls='--')
        plt.legend(loc='upper right', prop={'size':14}, frameon=False)
        plt.show()


figure()
# EOF: useCase15_matplotlib_oaSignal.py
