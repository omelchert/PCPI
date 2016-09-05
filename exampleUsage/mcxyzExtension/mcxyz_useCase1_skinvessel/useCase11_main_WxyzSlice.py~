import sys; sys.path.append('../../../')
import PyPCPI.voxelizedMedia.dataIO as io 

def main():
        # PATH TO INPUT FILES -------------------------------------------------
        iPath = "../inputData/skinvessel" 

        # VOLUMETRIC ENERGY DENSITY PER VOXEL --------------------------------- 
        ((x,y,z), Wxyz) = io.mcxyzio.getEnergyDensity(iPath)

        io.gpl.writeROI(y,z,Wxyz[:,:,x.size/2],'./dataUseCase11/Wxyz_yScan.dat')
        io.gpl.writeROI(x,z,Wxyz[:,y.size/2,:],'./dataUseCase11/Wxyz_xScan.dat')

main()
