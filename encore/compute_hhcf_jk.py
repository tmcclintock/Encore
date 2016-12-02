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
    #Jackknife subregion step size
    step = (edges[1]-edges[0])/ndivs
    Njk = int(ndivs**3)

    #Read in all halos
    all_halos = read_halos(outpath,Njk)

    #Treecorr interface
    config = {'nbins':nbins,'min_sep':limits[0],'max_sep':limits[1]}
    
    #Calculate RR autocorrelation once
    random_cat = treecorr.Catalog(x=randoms[:,0],y=randoms[:,1],z=randoms[:,2],config=config)
    RRa = treecorr.NNCorrelation(config)
    RRa.process(random_cat)

    #Get all autocorrelations
    #These are all of length Njks
    DDa_all, DRa_all = calculate_autos(outpath,config,all_halos,randoms,step,ndivs,Njk)

    #Get all cross correlations
    DDc_all,DRc_all,RRc_all = calculate_cross(outpath,config,all_halos,randoms,step,ndivs,Njk)

    #Allocate the totals
    DDt = treecorr.NNCorrelation(config)
    DRt = treecorr.NNCorrelation(config)
    RRt = treecorr.NNCorrelation(config)
    for i in range(Njk):
        DDt+=DDa_all[i]
        DDt.meanr[:]+=0.5*DDc_all[i].meanr[:]
        DDt.meanlogr[:]+=0.5*DDc_all[i].meanlogr[:]
        DDt.weight[:]+=0.5*DDc_all[i].weight[:]
        DDt.npairs[:]+=0.5*DDc_all[i].npairs[:]
        DDt.tot+=0.5*DDc_all[i].tot

        DRt+=DRa_all[i]
        DRt.meanr[:]+=0.5*DRc_all[i].meanr[:]
        DRt.meanlogr[:]+=0.5*DRc_all[i].meanlogr[:]
        DRt.weight[:]+=0.5*DRc_all[i].weight[:]
        DRt.npairs[:]+=0.5*DRc_all[i].npairs[:]
        DRt.tot+=0.5*DRc_all[i].tot

        RRt+=RRa
        RRt.meanr[:]+=0.5*RRc_all[i].meanr[:]
        RRt.meanlogr[:]+=0.5*RRc_all[i].meanlogr[:]
        RRt.weight[:]+=0.5*RRc_all[i].weight[:]
        RRt.npairs[:]+=0.5*RRc_all[i].npairs[:]
        RRt.tot+=0.5*RRc_all[i].tot
    DDt.write("test_jkhhcf.txt",RRt,DRt)

    print "HHCF JK not implemented yet!"
    return

def calculate_cross(outpath,config,all_halos,randoms,step,ndivs,Njk):
    """
    Calcualte the DD, DR and RR cross correlations
    """
    #Inititalize these objects
    DDc_all = []
    DRc_all = []
    RRc_all = []
    for i in range(Njk):
        DDc_all.append(treecorr.NNCorrelation(config))
        DRc_all.append(treecorr.NNCorrelation(config))
        RRc_all.append(treecorr.NNCorrelation(config))
        
    for index1 in range(Njk):
        halos1 = all_halos[index1]
        halo_cat1 = treecorr.Catalog(x=halos1[:,0],y=halos1[:,1],z=halos1[:,2],config=config)
        i = index1%ndivs
        j = (index1/ndivs)%ndivs
        k = index1/ndivs**2
        random_cat1 = treecorr.Catalog(x=randoms[:,0]+i*step,\
                                      y=randoms[:,1]+j*step,\
                                      z=randoms[:,2]+k*step,\
                                      config=config)
        for index2 in xrange(index1+1,Njk):
            halos2 = all_halos[index2]
            halo_cat2 = treecorr.Catalog(x=halos2[:,0],y=halos2[:,1],z=halos2[:,2],config=config)
            i = index2%ndivs
            j = (index2/ndivs)%ndivs
            k = index2/ndivs**2
            random_cat2 = treecorr.Catalog(x=randoms[:,0]+i*step,\
                                           y=randoms[:,1]+j*step,\
                                           z=randoms[:,2]+k*step,\
                                           config=config)
            DD = treecorr.NNCorrelation(config)
            DR = treecorr.NNCorrelation(config)
            RR = treecorr.NNCorrelation(config)
            DD.process(halo_cat1,halo_cat2)
            DR.process(halo_cat1,random_cat2)
            RR.process(random_cat1,random_cat2)
            DDc_all[index1]+=DD
            DDc_all[index2]+=DD
            DRc_all[index1]+=DR
            DRc_all[index2]+=DR
            RRc_all[index1]+=RR
            RRc_all[index2]+=RR
    return DDc_all,DRc_all,RRc_all

def calculate_autos(outpath,config,all_halos,randoms,step,ndivs,Njk):
    """
    Calculate the DD and DR autocorrelations
    """
    DDa_all = []
    DRa_all = []
    for index in range(Njk):
        halos = all_halos[index]
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

def read_halos(outpath,Njk):
    all_halos = []
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
        all_halos.append(halos)
    return np.array(all_halos)
