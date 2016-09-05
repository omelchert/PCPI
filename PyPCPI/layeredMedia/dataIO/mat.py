
import scipy 
import scipy.io

def writeROI(r,z,Wrz,fName='Wrz_dataROI.mat'):
        scipy.io.savemat(fName,dict(r=r, z=z, Wrz=Wrz))
                          
def readROI(fName):
        matDict = scipy.io.loadmat(fName)
        return matDict['r'], matDict['z'], matDict['Wrz']
