"""
Create random catalogs.
"""
import os,sys
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")

def create_halo_random_catalog(outpath,edges,Nh,ndivs,do_DM):
    """
    Create random catalogs for the full region
    and a JK subregion.
    """
    #Create the randoms output directory
    os.system("mkdir -p %s"%outpath+"/randoms")
    print "\tRandom directories created."

    M = 1 #Multiplicative factor
    Nr = int(M*Nh) #Randoms number
    Njk = ndivs**3
    Nrjk = int(Nr/Njk)
    Nrdm = int(Nr*1) #Arbitrary amount
    Nrjkdm = int(Nrdm/Njk)
    print "Creating random catalogs with:"
    print "\tN_halo_randoms full = %d"%Nr
    print "\tN_halo_randoms/JK   = %d"%Nrjk
    dl,dr = edges
    width = dr-dl
    x = np.random.rand(Nr)*width
    y = np.random.rand(Nr)*width
    z = np.random.rand(Nr)*width
    pos = np.array([x,y,z]).T
    np.savetxt(outpath+"/randoms/full_halo_random.txt",pos)

    widthjk = width/ndivs
    xjk = np.random.rand(Nrjk)*widthjk
    yjk = np.random.rand(Nrjk)*widthjk
    zjk = np.random.rand(Nrjk)*widthjk
    posjk = np.array([xjk,yjk,zjk]).T
    np.savetxt(outpath+"/randoms/jk_halo_random.txt",posjk)

    if do_DM:
        print "\tN_DM_randoms full   = %d"%Nrdm
        print "\tN_DM_randoms/JK     = %d"%Nrjkdm
        x = np.random.rand(Nrdm)*width
        y = np.random.rand(Nrdm)*width
        z = np.random.rand(Nrdm)*width
        pos = np.array([x,y,z]).T
        np.savetxt(outpath+"/randoms/full_dm_random.txt",pos)

        xjk = np.random.rand(Nrjkdm)*widthjk
        yjk = np.random.rand(Nrjkdm)*widthjk
        zjk = np.random.rand(Nrjkdm)*widthjk
        posjk = np.array([xjk,yjk,zjk]).T
        np.savetxt(outpath+"/randoms/jk_dm_random.txt",posjk)

    print "\tRandom catalogs created."
    return
