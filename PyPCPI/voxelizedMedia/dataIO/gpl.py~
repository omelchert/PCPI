import sys

def writeROI(x,y,Wxy,f=None):
        f = open(f,'w') if f else sys.stdout 
        f.write("# Nx Ny %d %d\n"%(x.size, y.size))
        f.write("%d "%(x.size))
        for ix in range(x.size):
            f.write("%lf "%(x[ix]))
        f.write("\n")
        for iy in range(y.size):
            f.write("%lf "%(y[iy]))
            for ix in range(x.size):
                f.write("%lf "%(Wxy[iy,ix]))
            f.write("\n")

