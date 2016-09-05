import sys; sys.path.append('../'); sys.path.append('../../')
import sourceVolume 
import irradiationSourceProfile
import poissonIntegralSolver.acousticObservables as solver

def writeROI(r,z,Wrz,f=None):
        f = open(f,'w') if f else sys.stdout 
        f.write("# Nr Nz %d %d\n"%(r.size, z.size))
        f.write("%d "%(r.size))
        for ir in range(r.size):
            f.write("%lf "%(r[ir]))
        f.write("\n")
        for iz in range(z.size):
            f.write("%lf "%(z[iz]))
            for ir in range(r.size):
                f.write("%lf "%(Wrz[ir,iz]))
            f.write("\n")


def main():
        (x,y,z),roi = sourceVolume.setROI((0.1,0.1,0.1),(200,200,200))
        
        sourceVolume.addAbsorbingLayer(((x,y,z),roi),0.02,0.01,10)
        sourceVolume.addAbsorbingCylinder(((x,y,z),roi),(0.051,0.051),0.01,20)
        
        #isp = irradiationSourceProfile.Gaussian((x,y),(0.05,0.05),0.01)
        isp = irradiationSourceProfile.flatTop((x,y),(0.05,0.05),0.01,2.)

        norm = irradiationSourceProfile.beamProfNormalization((x,y),isp)
        #exit()

        sourceVolume.propagateBeam((z,roi),isp) 

        t,p = solver.pressure(
                ((x,y,z),roi),
                (0.05,0.05,-0.2),
                1.)

        #writeROI(x,z,roi[:,y.size/2,:])
        writeROI(x,y,norm*isp)
        sys.exit()

        for i in range(t.size):
                print t[i],p[i]



main()
