import numpy as np

def fetchRawData(inFileName):
        rArray   = []
        zArray   = []
        ArzArray = []
        ctr = 0
        with open(inFileName,'r') as f:
            for line in f:
                    c = line.split()
                    # FILTER FOR r,z MESH PROPERTIES --------------------------
                    if len(c)>1 and c[0]!='#':
                        ctr+=1
                        if ctr == 5: 
                            dz, dr = float(c[0]), float(c[1])
                        if ctr == 6:
                            Nz, Nr = int(c[0]), int(c[1])

                    # FILTER FOR INTERNAL PHOTON DISTRIBUTION -----------------
                    if len(c)==1 and c[0]=='A_rz':
                            c = f.next().split()
                            while(len(c)>=1):
                                ArzArray.append(c)
                                c = f.next().split()

        #print Nz,Nr
        #print Nz*Nr
        #print np.asarray(ArzArray,dtype=float).size
        rAxis = np.linspace(0,dr*Nr,Nr,endpoint=False)
        zAxis = np.linspace(0,dz*Nz,Nz,endpoint=False)
        ArzArray = np.asarray(ArzArray,dtype=float).reshape((Nr,Nz))

        return rAxis,zAxis,ArzArray
