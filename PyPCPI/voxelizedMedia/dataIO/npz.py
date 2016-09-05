import numpy as np


def writeROI(x,y,z,Wxyz,fName='Wxyz_dataROI.npz'):
        np.savez_compressed(fName, x=x, y=y, z=z, Wxyz=Wxyz)

def readROI(file):
        npzDict = np.load(file)
        return npzDict['x'], npzDict['y'], npzDict['z'], npzDict['Wxyz']
