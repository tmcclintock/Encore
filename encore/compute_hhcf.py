"""
Compute the halo-halo correlation function.
"""
import os,sys
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")
try: import treecorr
except ImportError: raise Exception("Must install treecorr.")


def compute_hhcf(outpath,nbins,limits,edges,do_randoms,do_JK,ndivs):
    #Step 1: calculate the mean masses
    if os.path.exists(outpath+"/info_files/halo_mass_info.txt"):
        halo_mass_info = np.loadtxt(outpath+"/info_files/halo_mass_info.txt")
    else:
        halo_mass_info = calculate_mean_mass(outpath,do_JK,ndivs)
        np.savetxt(outpath+"/info_files/halo_mass_info.txt",halo_mass_info)
    print "Halo masses averaged, results:"
    print "\tMmean  = %e"%halo_mass_info[0]
    print "\tMtotal = %e"%halo_mass_info[1]
    print "\tNhalos = %d"%halo_mass_info[2]

    #Step 2: create random catalogs.
    randpath = outpath+"/halohalo_correlation_function/randoms/full_random.txt"
    jkrandpath = outpath+"/halohalo_correlation_function/randoms/jk_random.txt"
    if os.path.exists(randpath): randoms = np.loadtxt(randpath)
    else: raise Exception("Must create random catalog first.")
    if do_JK: randomsjk = np.loadtxt(jkrandpath)
    else: randomsjk = None
    print "Random catalogs loaded.\nUsing random catalogs with:"
    print "\tN_randoms full = %d"%len(randoms)
    if do_JK: print "\tN_randoms JK   = %d"%len(randomsjk)

    #Step 3: calculate the full HH correlation function
    calcalate_hhcf_full(outpath,nbins,limits,halo_mass_info[2],randoms)

    #Step 4: calculate the HH correlation function within JK subregions
    if do_JK: 
        calculate_JK_singles_hhcf(outpath,nbins,limits,edges,halo_mass_info[2],randomsjk,ndivs)

    print "Jackknifed halo-halo correlation function not implemented yet!"
    return

def calculate_JK_singles_hhcf(outpath,nbins,limits,edges,Nh,randoms,ndivs):
    """
    Calculate the halo-halo correlation function
    for the JK subregions.
    """
    #Interface with treecorr
    config = {'nbins':nbins,'min_sep':limits[0],'max_sep':limits[1]}
    random_cat = treecorr.Catalog(x=randoms[:,0],y=randoms[:,1],z=randoms[:,2],config=config)
    RR = treecorr.NNCorrelation(config)
    RR.process(random_cat)
    #RR.write("RR_test.txt")
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

def calcalate_hhcf_full(outpath,nbins,limits,Nh,randoms):
    """
    Calculate the halo-halo correlation function
    for the full volume.
    """
    redpath = outpath+"/reduced_halo_cats/reduced_halo_cat.txt"
    infile = open(redpath,"r")
    halos = np.zeros((int(Nh),3))
    i = 0
    for line in infile:
        if line[0] is "#": continue
        parts = line.split()
        halos[i,:] = float(parts[8]),float(parts[9]),float(parts[10])
        i+=1
    infile.close()
    #Interface with treecorr
    config = {'nbins':nbins,'min_sep':limits[0],'max_sep':limits[1]}
    halo_cat = treecorr.Catalog(x=halos[:,0],y=halos[:,1],z=halos[:,2],config=config)
    random_cat = treecorr.Catalog(x=randoms[:,0],y=randoms[:,1],z=randoms[:,2],config=config)
    DD = treecorr.NNCorrelation(config)
    DR = treecorr.NNCorrelation(config)
    RR = treecorr.NNCorrelation(config)
    DD.process(halo_cat)
    DR.process(halo_cat,random_cat)
    RR.process(random_cat)
    DD.write(outpath+"/halohalo_correlation_function/full_hhcf/full_hhcf.txt",RR,DR)

    print "Calculating the FULL hh correlation function complete."
    return

def calculate_mean_mass(outpath,do_JK,ndivs):
    """
    Calculate the mean masses of the halos.
    If JK is on then calculate mean masses
    for the halos in the JK files.
    """
    Mtotal = 0.0
    N = 0
    redpath = outpath+"/reduced_halo_cats/reduced_halo_cat.txt"
    infile = open(redpath,"r")
    for line in infile:
        if line[0] is "#": continue
        parts = line.split()
        m = float(parts[2])
        Mtotal += m
        N += 1
    Mmean = Mtotal/N
    print "Mean masses for JK subregions not implemented yet!"
    return [Mmean,Mtotal,N]
