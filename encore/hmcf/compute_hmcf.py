"""
Compute the halo-matter correlation function.
"""
import os,sys
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")
try: import treecorr
except ImportError: raise Exception("Must install treecorr.")

#Pull out the indices
indices = {}
with open("rockstar_config") as myfile:
    for line in myfile:
        name, var = line.partition("=")[::2]
        indices[name.strip()] = int(var)
x_index = indices['x']
y_index = indices['y']
z_index = indices['z']
m_index = indices['m']

def compute_hmcf(outpath,nbins,limits,edges,do_JK,ndivs):
    """
    Compute the halo-matter correlation function.

    Inputs:
       outpath: the base directory where the output directories exist
       nbins: the number of radial bins for the xi_hm
       limits: an array with two entries with the min/max separation
       edges: the spatial edges of the simulation, i.e. xmin and xmax
       do_JK: boolean for wheather to calculate jackknife values
       ndivs: number of JK subregions
    """
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
    Mmean,Mtotal,Nh = halo_mass_info

    #Step 2: get random catalogs.
    halorandpath = outpath+"/randoms/jk_halo_random.txt"
    dmrandpath   = outpath+"/randoms/jk_dm_random.txt"
    if os.path.exists(halorandpath) and os.path.exists(dmrandpath): 
        halorandoms = np.loadtxt(halorandpath)
        dmrandoms = np.loadtxt(dmrandpath)
    else: raise Exception("Must create random catalog first.")
    print "Random catalogs loaded.\nUsing random catalogs with:"
    print "\tN_halo_randoms JK = %d"%len(halorandoms)
    print "\tN_dm_randoms JK   = %d"%len(dmrandoms)

    #Step 3: calculate the full HH correlation function
    calcalate_hmcf(outpath,nbins,limits,edges,halorandoms,dmrandoms,ndivs)

    print "Halo-matter correlation function not implemented yet!"
    return

def calcalate_hmcf(outpath,nbins,limits,edges,halorandoms,dmrandoms,ndivs):
    """
    Calcualte the halo-matter correlation function.

    Note: both set of randoms are already jackknifed.
    """
    #Jackknife subregion step size
    step = (edges[1]-edges[0])/ndivs
    Njk = int(ndivs**3)

    #Read in all halos; Note: the DM particles will be read in one at a time
    all_halos = read_halos(outpath,Njk)

    #Treecorr interface
    config = {'nbins':nbins,'min_sep':limits[0],'max_sep':limits[1]}

    #Calculate RR once
    halorandom_cat = treecorr.Catalog(x=halorandoms[:,0],y=halorandoms[:,1],z=halorandoms[:,2],config=config)
    dmrandom_cat   = treecorr.Catalog(x=dmrandoms[:,0],y=dmrandoms[:,1],z=dmrandoms[:,2],config=config)
    RRa = treecorr.NNCorrelation(config)
    RRa.process(halorandom_cat,dmrandom_cat)
    print "HMCF RR autocorrelation calculated."

    #Calculate all DR/RD pairs
    #NOTE: DR is for halos and RD is for dark matter
    #DR_all = compute_DR(outpath,config,step,halorandoms,ndivs)
    #RD_all,mapping = compute_RD(config,step,dmrandoms,ndivs) #NEED TO ADD A DM PATH!!

    #Calculate all DD pairs
    #DD_all = compute_DD(outpath,dmpath,config,step,ndivs)
    
    #Jackknife everything and output results
    #compute_hmcf_jk(outpath,DD_all,DR_all,RD_all,RRa,ndivs)

    print "HMCF JK not implemented yet!"
    return

def read_halos(outpath,Njk):
    all_halos = []
    jkpath = outpath+"/JK_halo_cats/jk_halo_cat_%d.txt"
    for index in range(Njk):
        infile = open(jkpath%index,"r")
        halos = [] #Will be Nhjk X 3
        for line in infile:
            if line[0] is "#": continue
            parts = line.split()
            halos.append([float(parts[x_index]),float(parts[y_index]),float(parts[z_index])])
        halos = np.array(halos)
        infile.close()
        all_halos.append(halos)
    return np.array(all_halos)
