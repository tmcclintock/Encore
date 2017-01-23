"""
Compute the halo-matter correlation function.
"""
import os,sys
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")
try: import treecorr
except ImportError: raise Exception("Must install treecorr.")
try: import pygadgetreader as pgr
except ImportError: raise Exception("Must install pygadgetreader.")

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

def compute_hmcf(outpath,nbins,limits,edges,do_JK,ndivs,DSF):
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
    halorandpath = outpath+"/randoms/full_halo_random.txt"
    dmrandpath   = outpath+"/randoms/full_dm_random.txt"
    if os.path.exists(halorandpath) and os.path.exists(dmrandpath): 
        halorandoms = np.loadtxt(halorandpath)
        dmrandoms = np.loadtxt(dmrandpath)
    else: raise Exception("Must create random catalog first.")
    print "Random catalogs loaded.\nUsing random catalogs with:"
    print "\tN_halo_randoms full = %d"%len(halorandoms)
    print "\tN_dm_randoms full   = %d"%len(dmrandoms)

    #Step 3: calculate the full HM correlation function
    calcalate_hmcf_full(outpath,nbins,limits,Nh,halorandoms,dmrandoms,DSF)

    #Step 4: calculate the HH correlation function within JK subregions
    if do_JK:
        halorandpath = outpath+"/randoms/jk_halo_random.txt"
        dmrandpath   = outpath+"/randoms/jk_dm_random.txt"
        if os.path.exists(halorandpath) and os.path.exists(dmrandpath): 
            halorandoms = np.loadtxt(halorandpath)
            dmrandoms = np.loadtxt(dmrandpath)
        import compute_hmcf_jk
        compute_hmcf_jk.calculate_JK_hmcf(outpath,nbins,limits,edges,Nh,halorandoms,dmrandoms,ndivs,DSF)

    print "HMCF JK not implemented yet!"
    print "Halo-matter correlation function not implemented yet!"
    return

def calcalate_hmcf_full(outpath,nbins,limits,Nh,halorandoms,dmrandoms,DSF):
    """
    Calcualte the halo-matter correlation function.
    """
    print "Calculating full halo-matter correlation function."
    #Read in the halos
    redpath = outpath+"/reduced_halo_cats/reduced_halo_cat.txt"
    infile = open(redpath,"r")
    halos = np.zeros((int(Nh),3))
    i = 0
    for line in infile:
        if line[0] is "#": continue
        parts = line.split()
        halos[i,:] = float(parts[x_index]),float(parts[y_index]),float(parts[z_index])
        i+=1
    infile.close()

    #Read in the dm particles
    dmpath = outpath+"/down_sampled_dm/down_sampled_dm_DSF%d"%DSF
    dm = pgr.readsnap(dmpath,"pos","dm")

    #Treecorr interface
    config = {'nbins':nbins,'min_sep':limits[0],'max_sep':limits[1]}
    halorandom_cat = treecorr.Catalog(x=halorandoms[:,0],y=halorandoms[:,1],z=halorandoms[:,2],config=config)
    dmrandom_cat   = treecorr.Catalog(x=dmrandoms[:,0],y=dmrandoms[:,1],z=dmrandoms[:,2],config=config)
    halo_cat = treecorr.Catalog(x=halos[:,0],y=halos[:,1],z=halos[:,2],config=config)
    dm_cat = treecorr.Catalog(x=dm[:,0],y=dm[:,1],z=dm[:,2],config=config)

    DD = treecorr.NNCorrelation(config)
    DR = treecorr.NNCorrelation(config)
    RD = treecorr.NNCorrelation(config)
    RR = treecorr.NNCorrelation(config)

    DD.process(halo_cat,dm_cat)
    DR.process(halo_cat,halorandom_cat)
    RD.process(dm_cat,dmrandom_cat)
    RR.process(halorandom_cat,dmrandom_cat)
    DD.write(outpath+"/halomatter_correlation_function/full_hmcf/full_hmcf.txt",RR,DR,RD)

    print "\tHalo-matter correlation function full complete."
    return

#MOVE THE BELOW TO A JK FILE
def read_halos(outpath,Njk):
    all_halos = []
    jkpath = outpath+"/JK_halo_cats/jk_halo_cat_%d.txt"
    for index in xrange(0,Njk):
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
