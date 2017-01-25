"""
Compute the halo-halo correlation function.
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

def compute_hhcf(outpath,randompath,nbins,limits,edges,do_JK,ndivs):
    """
    Compute the halo-halo correlation function.

    Inputs:
       outpath: the base directory where the output directories exist
       nbins: the number of radial bins for the xi_hh
       limits: an array with two entries with the min/max separation
       edges: the spatial edges of the simulation, i.e. xmin and xmax
       do_JK: boolean for wheather to calculate jackknife values
       ndivs: number of JK subregions
    """
    print "Calculating halo-halo mass function."

    #Step 0: create the output directories for HHCF
    create_hhcf_directories(outpath)

    #Step 1: calculate the mean masses
    if os.path.exists(outpath+"/info_files/halo_mass_info.txt"):
        halo_mass_info = np.loadtxt(outpath+"/info_files/halo_mass_info.txt")
    else:
        halo_mass_info = calculate_mean_mass(outpath,do_JK,ndivs)
        np.savetxt(outpath+"/info_files/halo_mass_info.txt",halo_mass_info)
    print "\tHalo masses averaged, results:"
    print "\t\tMmean  = %e"%halo_mass_info[0]
    print "\t\tMtotal = %e"%halo_mass_info[1]
    print "\t\tNhalos = %d"%halo_mass_info[2]
    Mmean,Mtotal,Nh = halo_mass_info

    #Step 2: get random catalogs.
    randpath = randompath+"/randoms/full_halo_random.txt"
    jkrandpath = randompath+"/randoms/jk_halo_random.txt"
    if os.path.exists(randpath): randoms = np.loadtxt(randpath)
    else: raise Exception("Must create random catalog first.")
    if do_JK: randomsjk = np.loadtxt(jkrandpath)
    else: randomsjk = None
    print "\tRandom catalogs loaded. Using random catalogs with:"
    print "\t\tN_randoms full = %d"%len(randoms)
    if do_JK: print "\t\tN_randoms JK   = %d"%len(randomsjk)

    #Step 3: calculate the full HH correlation function
    calcalate_hhcf_full(outpath,nbins,limits,Nh,randoms)

    #Step 4: calculate the HH correlation function within JK subregions
    if do_JK:
        import compute_hhcf_jk
        compute_hhcf_jk.calculate_JK_hhcf(outpath,nbins,limits,edges,Nh,randomsjk,ndivs)

    print "Halo-halo correlation function complete."
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
        halos[i,:] = float(parts[x_index]),float(parts[y_index]),float(parts[z_index])
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

    print "\tFull halo-halo correlation function complete."
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
        m = float(parts[m_index])
        Mtotal += m
        N += 1
    Mmean = Mtotal/N
    print "\tMean masses for JK subregions not implemented yet!"
    return [Mmean,Mtotal,N]

def create_hhcf_directories(outpath):
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/full_hhcf")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/JK_single")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/JK_combined")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/cov_matrix")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/final_hhcf")
    print "\tHalo-halo correlation function directories created."
