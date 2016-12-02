"""
Compute the JK halo-halo correlation function
and the covariance matrix.
"""
import os,sys
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")
try: import treecorr
except ImportError: raise Exception("Must install treecorr.")

def calculate_JK_hhcf(outpath,nbins,limits,edges,Nh,randoms,ndivs):
    """
    Calculate the halo-halo correlation function for the JK subregions.
    """
    #Interface with treecorr
    config = {'nbins':nbins,'min_sep':limits[0],'max_sep':limits[1]}
    random_cat = treecorr.Catalog(x=randoms[:,0],y=randoms[:,1],z=randoms[:,2],config=config)
    RR = treecorr.NNCorrelation(config)
    RR.process(random_cat)
    #Jackknife subregion step size
    step = (edges[1]-edges[0])/ndivs
    Njk = int(ndivs**3)
    TC_output = [] #Will be of length Njk*(Njk-1)/2
    jkpath = outpath+"/JK_halo_cats/jk_halo_cat_%d.txt"
    for index in range(1):#Njk):
        infile = open(jkpath%index,"r")
        halos = [] #Will be Nhjk X 3
        for line in infile:
            if line[0] is "#": continue
            parts = line.split()
            halos.append([float(parts[8]),float(parts[9]),float(parts[10])])
        halos = np.array(halos)
        infile.close()
        i = index%ndivs
        j = (index/ndivs)%ndivs
        k = index/ndivs**2
        random_cat = treecorr.Catalog(x=randoms[:,0]+i*step,\
                                      y=randoms[:,1]+j*step,\
                                      z=randoms[:,2]+k*step,\
                                      config=config)
        halo_cat = treecorr.Catalog(x=halos[:,0],y=halos[:,1],z=halos[:,2],config=config)
        DD = treecorr.NNCorrelation(config)
        DR = treecorr.NNCorrelation(config)
        DD.process(halo_cat)
        DR.process(halo_cat,random_cat)
        xi,varxi = DD.calculateXi(RR,DR)
        print xi
        print "weight",DD.weight
        DD.write("JK_test.txt",RR,DR)
        
    print "HHCF JK singles not implemented yet!"
    return
