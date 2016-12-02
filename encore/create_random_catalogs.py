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
    M = 10 #Multiplicative factor
    Nr = int(M*Nh) #Randoms number
    Njk = ndivs**3
    Nrjk = int(Nr/Njk)
    print "Creating random catalogs with:"
    print "\tN_randoms full = %d"%Nr
    print "\tN_randoms JK   = %d"%Nrjk
    dl,dr = edges
    width = dr-dl
    x = np.random.rand(Nr)*width
    y = np.random.rand(Nr)*width
    z = np.random.rand(Nr)*width
    pos = np.array([x,y,z]).T
    np.savetxt(outpath+"/halohalo_correlation_function/randoms/full_random.txt",pos)

    widthjk = width/ndivs
    xjk = np.random.rand(Nrjk)*widthjk
    yjk = np.random.rand(Nrjk)*widthjk
    zjk = np.random.rand(Nrjk)*widthjk
    posjk = np.array([xjk,yjk,zjk]).T
    np.savetxt(outpath+"/halohalo_correlation_function/randoms/jk_random.txt",posjk)

    print "Creation of random catalogs complete."
    return
