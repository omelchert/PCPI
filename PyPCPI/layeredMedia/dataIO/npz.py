import numpy as np


def writeROI(r,z,Wrz,fName='Wrz_dataROI.npz'):
        np.savez_compressed(fName, r=r, z=z, Wrz=Wrz)

def readROI(fName):
        npzDict = np.load(fName)
        return npzDict['r'], npzDict['z'], npzDict['Wrz']
