""" FILE: useCase011_main_convolve.py

AUTHOR: O. Melchert
DATE: 29.08.2016
"""
import sys; sys.path.append('../../../')
import PyPCPI.layeredMedia.dataIO as io
import PyPCPI.layeredMedia.polarConvolution as conv 
                                    

def main():
        R0 = 0.15 # 1/e2 WIDTH OF BEAM PROFILE
        inFileName = sys.argv[1]
        outFileName = "./dataUseCase11/Wrz_%s_ISPGauss_R0%3.2lf.npz"%\
                      (inFileName.split('/')[-1].split('.')[0],R0)

        # FETCH IMPULSE RESPONSE FROM mco. FILE
        r,z,Arz = io.mcmlio.fetchRawData(inFileName)

        # SET CUSTOM IRRADIATION SOURCE PROFILE FROM OSG LIBRARY
        myIsp = conv.isp.Gaussian(r, R0)
        # CONVOLVE IMPULSE RESPONSE USING CUSTOM IRRADIATION SOURCE PROFILE
        Wrz = conv.convolveROI(r, z, Arz, myIsp)

        # SAVE DATA IN NUMPY FORMAT
        io.npz.writeROI(r, z, Wrz, outFileName)
        

main()
# EOF: useCase011_main_convolve.py
