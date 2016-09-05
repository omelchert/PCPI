
import scipy 
import scipy.io

def writeROI(x,y,z,Wxyz,fName='Wxyz_dataROI.mat'):
        scipy.io.savemat(fName,dict(x=x, y=y, z=z, Wxyz=Wxyz))
                          
def readROI(file):
        matDict = scipy.io.loadmat(file)
        return matDict['x'], matDict['y'], matDict['z'], matDict['Wxyz']
