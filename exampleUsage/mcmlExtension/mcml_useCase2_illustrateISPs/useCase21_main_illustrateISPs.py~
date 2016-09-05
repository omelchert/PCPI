""" FILE: useCase21_main_illustrateISPs.py

Example script illustrating the convolution of the impulse response obtained
for a infinitely narrow photon beam (stored in a given .mco file) to the 
material response for a photon beam of finite extend and given transversal
axially symmetric irradiation source profile.

AUTHOR: O. Melchert
DATE: 29.08.2016
"""
import sys; sys.path.append('../../../')
import numpy as np
import PyPCPI.layeredMedia.dataIO as io 
import PyPCPI.layeredMedia.polarConvolution as conv


def main():
        # INITIALIZATION AND DECLARATION --------------------------------------
        iPathMco = sys.argv[1] 

        gParam = (0.0,0.15)  # GAUSSIAN ISP 
        dParam = (0.10, 0.15, 0.03, 0.08) # DONUT ISP
        fParam = (0.10, 0.05)  # FLAT-TOP ISP

        def convHelper(isp,s):
            Wrz = conv.convolveROI(r, z, Arz, isp )
            io.gpl.writeROI(r, z, Wrz, "./dataUseCase21/Wrz_%s.prof"%(s))
            io.npz.writeROI(r, z, Wrz, "./dataUseCase21/Wrz_%s.npz"%(s))

        def customFlatTop(r, (R0, D0)):
            condList = [r<R0, r>=R0]
            funcList = [lambda r: 1, lambda r: np.exp(-2*(r-R0)**2/D0**2)]
            return np.piecewise(r,condList,funcList)

        # FETCH IMPULSE RESPONSE FROM .mco FILE -------------------------------
        r,z,Arz = io.mcmlio.fetchRawData(iPathMco)

        # YIELD RESPONSE FOR CUSTOM IRRADIATION SOURCE PROFILES ---------------
        convHelper(conv.isp.Gaussian(r, gParam), "Gaussian")
        convHelper(conv.isp.flatTopDonut(r,dParam), "Donut")
        convHelper(customFlatTop(r, fParam), "flatTop")

        
main()
# EOF: useCase21_main_illustrateISPs.py
