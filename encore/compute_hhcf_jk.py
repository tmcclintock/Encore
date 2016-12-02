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

    Inputs:
       outpath: the base directory where the output directories exist
       nbins: the number of radial bins to calculate xi_hh for
       limits: an array with two entries with the min/max separation
       edges: the spatial edges of the simulation, i.e. xmin and xmax
       do_JK: boolean for wheather to calculate jackknife values
       ndivs: number of JK subregions
    """
    #Treecorr interface
    config = {'nbins':nbins,'min_sep':limits[0],'max_sep':limits[1]}
    
    #Calculate RR autocorrelation once
    random_cat = treecorr.Catalog(x=randoms[:,0],y=randoms[:,1],z=randoms[:,2],config=config)
    RRa = treecorr.NNCorrelation(config)
    RRa.process(random_cat)

    #Jackknife subregion step size
    step = (edges[1]-edges[0])/ndivs
    Njk = int(ndivs**3)

    #Get all autocorrelations
    DDa_all, DRa_all = calculate_autos(outpath,config,randoms,step,ndivs,Njk)
    print len(DDa_all),len(DRa_all)
    print "HHCF JK not implemented yet!"
    return

def calculate_autos(outpath,config,randoms,step,ndivs,Njk):
    """
    Calculate the DD and DR autocorrelations
    """
    DDa_all = []
    DRa_all = []
    jkpath = outpath+"/JK_halo_cats/jk_halo_cat_%d.txt"
    for index in range(Njk):
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
        DDa_all.append(DD)
        DRa_all.append(DR)
    print "For HHCF all DD and DR autocorrelations computed."
    return DDa_all,DRa_all
