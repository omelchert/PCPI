""" FILE: useCase13_main_oaSignal.py

AUTHOR: O. Melchert
DATE: 29.08.2016
"""
import sys; sys.path.append('../../../')
import PyPCPI.layeredMedia.dataIO as io 
import PyPCPI.layeredMedia.poissonIntegralSolver as ps

def main():
        fName = sys.argv[1] # E.g. 'Wrz_SkinOA_1G_ISPGauss_R00.80.npz'

        c0 = 1.0 # speed of sound in medium 
        G = 1.0 # efficiency of transforming absorbed energy to acoustic stress
        zD = -2.0 # z-location of pointlike detector 
        Nphi = 1 # number of azimutal mesh points considered for integration 
                 # (in case of on-axis signals Nphi=1 suffices)

        # READ VOLUMETRIC ENERGY DENSITY FROM COMPRESSED .npz FILE
        (r,z,Wrz) = io.npz.readROI(fName)

        t,p = ps.pressure((r,z,G*Wrz),(0.0,zD),c0,Nphi) 
        
        tau = t + zD/c0
        print "# (c0 tau) (p(tau))"
        for i in range(tau.size):
                print c0*tau[i],p[i]


main()
# EOF: useCase13_main_oaSignal.py
