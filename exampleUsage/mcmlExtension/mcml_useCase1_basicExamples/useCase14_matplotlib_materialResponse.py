""" FILE: useCase14_figure_materialResponse.py

stand-alone python script for plotting full material response (stored in .npz
format) using matplotlib tools

AUTHOR: O. Melchert
DATE: 10.01.2017
"""
import sys; sys.path.append('../../../')
import matplotlib.pyplot as plt
import PyPCPI.layeredMedia.dataIO as io

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

test_figure()
