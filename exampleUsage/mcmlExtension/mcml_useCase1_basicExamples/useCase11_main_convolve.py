""" FILE: useCase011_main_convolve.py

NOTE: 
  10.01.2017 - added functions figure(), test_figure() for plotting results 
      using matplotlib tools

AUTHOR: O. Melchert
DATE: 29.08.2016
"""
import sys; sys.path.append('../../../')
import matplotlib.pyplot as plt
import PyPCPI.layeredMedia.dataIO as io
import PyPCPI.layeredMedia.polarConvolution as conv 
                                    

def figure((r,z,W)):
        c = plt.contourf(z,r,W,100)
        cBar = plt.colorbar(c)
        cBar.set_label('$W(r,z)$ (J/m$^2$)')
        plt.title('Absorbed energy density')
        plt.xlabel('$z$ (cm)')
        plt.ylabel('$r$ (cm)')
        plt.show()


def test_figure():
        figure(io.npz.readROI(
            './dataUseCase11/Wrz_2Layer_mcml_ISPGauss_R00.15.npz'))


def main():
        R0 = 0.15 # 1/e2 WIDTH OF BEAM PROFILE
        inFileName = sys.argv[1]
        outFileName = "./dataUseCase11/Wrz_%s_ISPGauss_R0%3.2lf.npz"%\
                      (inFileName.split('/')[-1].split('.')[0],R0)

        # FETCH IMPULSE RESPONSE FROM mco. FILE
        r,z,Arz = io.mcmlio.fetchRawData(inFileName)

        # SET CUSTOM IRRADIATION SOURCE PROFILE FROM OSG LIBRARY
        myIsp = conv.isp.Gaussian(r, (0.,R0))
        # CONVOLVE IMPULSE RESPONSE USING CUSTOM IRRADIATION SOURCE PROFILE
        Wrz = conv.convolveROI(r, z, Arz, myIsp)

        # SAVE DATA IN NUMPY FORMAT
        io.npz.writeROI(r, z, Wrz, outFileName)

        # SHOW ISP RESPONSE USING MATPLOTLIB
        figure((r,z,Wrz))

        
        
#test_figure()
main()
# EOF: useCase011_main_convolve.py
