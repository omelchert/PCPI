""" FILE: useCase22_main_axialOASignal.py

AUTHOR: O. Melchert
DATE: 30.08.2016
"""
import sys; sys.path.append('../../../')
import PyPCPI.layeredMedia.dataIO as io 
import PyPCPI.layeredMedia.poissonIntegralSolver as ps

def main():
        iPathNpz = sys.argv[1] 

        (c0, G) = (1.0, 1.0)
        (rD, zD) = (0.0, -3.0)

        (r,z,Wrz) = io.npz.readROI(iPathNpz)
        t,p = ps.pressure((r,z,G*Wrz),(rD,zD),c0,1) 

        print "# (c0 tau) (p)"
        for i in range(t.size):
                print c0*t[i] + zD, p[i]

main()
# EOF: useCase22_main_axialOASignal.py
