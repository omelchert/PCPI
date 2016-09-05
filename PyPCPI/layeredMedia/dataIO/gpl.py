
import sys

def writeROI_sp(r,z,Wrz,f=None):
        f = open(f,'w') if f else sys.stdout 
        f.write("# Nr Nz %d %d\n"%(r.size, z.size))
        f.write("# (r) (z) (Wrz)\n")
        for ir in range(r.size):
            for iz in range(z.size):
                f.write("%lf %lf %lf\n"%(r[ir],z[iz],Wrz[ir,iz]))
            f.write("\n")

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

