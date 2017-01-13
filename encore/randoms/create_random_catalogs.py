"""
Create random catalogs.
"""
import os,sys
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")

def create_halo_random_catalog(outpath,edges,Nh,ndivs):
    """
    Create random catalogs for the full region
    and a JK subregion.
    """
    M = 1 #Multiplicative factor
    Nr = int(M*Nh) #Randoms number
    Njk = ndivs**3
    Nrjk = int(Nr/Njk)
    Nrjkdm = Nrjk*2#00 #Halos have to have at LEAST 200 particles
    print "Creating random catalogs with:"
    print "\tN_halo_randoms full = %d"%Nr
    print "\tN_halo_randoms/JK   = %d"%Nrjk
    print "\tN_DM_randoms/JK     = %d"%Nrjkdm
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

    xjkdm = np.random.rand(Nrjkdm)*widthjk
    yjkdm = np.random.rand(Nrjkdm)*widthjk
    zjkdm = np.random.rand(Nrjkdm)*widthjk
    posjk = np.array([xjkdm,yjkdm,zjkdm]).T
    np.savetxt(outpath+"/randoms/jk_dm_random.txt",posjk)

    print "Creation of random catalogs complete."
    return
