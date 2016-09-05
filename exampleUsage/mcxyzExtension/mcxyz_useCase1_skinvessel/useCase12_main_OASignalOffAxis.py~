import sys; sys.path.append('../../../')
import numpy as np
import PyPCPI.voxelizedMedia.dataIO as io 
import PyPCPI.voxelizedMedia.poissonIntegralSolver as ps


def main():

        # SIMULATION PARAMETERS -----------------------------------------------
        basePath = "../inputData/skinvessel"
        G = 0.138 # Grueneisen parameter
        c0 = 150000. # speed of sound (cm/s)
        xD = float(sys.argv[1]) # detector x-position (cm)
        yD = float(sys.argv[2]) # detector y-position (cm)
        zD = float(sys.argv[3]) # detector z-position (cm)

        # VOLUMETRIC ENERGY DENSITY PER VOXEL --------------------------------- 
        ((x,y,z), Wxyz) = io.mcxyzio.getEnergyDensity(basePath)
        # OA EXCESS PRESSURE SIGNAL -------------------------------------------
        t, p = ps.pressure(((x,y,z),G*Wxyz),(xD,yD,zD),c0)

        tau = t + zD/c0
        print "# (c0 tau) (p(tau))"
        for i in range(tau.size):
                print c0*tau[i],p[i]


main()
# EOF: mcxyz_useCase12_main_offYAxis.py 
