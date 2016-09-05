""" FILE: useCase12_main_yieldGnuplotOutput.py

AUTHOR: O. Melchert
DATE: 29.08.2016
"""
import sys; sys.path.append('../../../')
import PyPCPI.layeredMedia.dataIO as io


def main():
        fName = sys.argv[1] # E.g. 'Wrz_SkinOA_1G_ISPGauss_R00.80.npz'

        # READ VOLUMETRIC ENERGY DENSITY FROM COMPRESSED .npz FILE
        r,z,Wrz = io.npz.readROI(fName)

        # WRITE VOLUMETRIC ENERGY DENSITY IN GNUPLOT IMAGE FORMAT 
        io.gpl.writeROI(r,z,Wrz)


main()
# EOF: useCase12_main_yieldGnuplotOutput.py
