""" FILE: useCase22_main_axialOASignal.py

AUTHOR: O. Melchert
DATE: 30.08.2016
"""
import sys; sys.path.append('../../../')
import PyPCPI.layeredMedia.dataIO as io 
import PyPCPI.layeredMedia.poissonIntegralSolver as ps
import PyPCPI.signalPostProcessing.finitePulse as pp

def main():
        iPathNpz = sys.argv[2] 

        c0, G = 150000.0, 1.0
        rD, zD = 0.0, float(sys.argv[1])
        tp = 2.*10**(-8) # gaussian pulse width (s) 

        (r,z,Wrz) = io.npz.readROI(iPathNpz)
        t,p = ps.pressure((r,z,G*Wrz),(rD,zD),c0,1) 

        ptp = pp.convolveGauss(p,(t,tp))

        print "# tp = %g s"%(tp)
        print "# (c0 tau) (p)"
        for i in range(t.size):
                print c0*t[i] + zD, p[i], ptp[i]

main()
# EOF: useCase22_main_axialOASignal.py
