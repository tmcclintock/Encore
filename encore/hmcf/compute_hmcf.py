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

def compute_hmcf(outpath, catalog, halorandompath,
                 dmpath, dmrandompath, edges,
                 nbins, Rlimits, do_JK, 
                 jkcatalog, jkdms, jkrands, jkdmrands,
                 ndivs):
    """Compute the halo-matter correlation function.
    """
    print "Calculating halo-matter correlation function."

    #Step 1: create the HMCF output directories
    create_hmcf_directories(outpath)

    #Step 2: get random catalogs.
    if os.path.exists(halorandompath) and os.path.exists(dmrandompath): 
        halorandoms = np.loadtxt(halorandompath)
        dmrandoms = np.loadtxt(dmrandompath)
    else: raise Exception("Must create random catalog first.")
    print "\tRandom catalogs loaded.\n\tUsing random catalogs with:"
    print "\t\tN_halo_randoms full = %d"%len(halorandoms)
    print "\t\tN_dm_randoms full   = %d"%len(dmrandoms)

    #Step 3: calculate the full HM correlation function
    Nhalos = calcalate_hmcf_full(outpath, catalog, dmpath, nbins, Rlimits, halorandoms, dmrandoms)

    #Step 4: calculate the HM correlation function within JK subregions
    if do_JK:
        halorandoms = np.loadtxt(jkrands)
        dmrandoms = np.loadtxt(jkdmrands)
        import compute_hmcf_jk
        compute_hmcf_jk.calculate_JK_hmcf(outpath, jkcatalog, jkdms, edges, nbins, Rlimits, 
                                          halorandoms, dmrandoms, ndivs)
        print "still working on JK HMCF"
        sys.exit()
    print "Halo-matter correlation function calculated."
    return

def calcalate_hmcf_full(outpath, catalog, dmpath, nbins, Rlimits, halorandoms, dmrandoms):
    """Calcualte the halo-matter correlation function.
    """
    print "Calculating full halo-matter correlation function."
    #Read in the halos
    Nh = 0 #Number of halos
    halos = []
    with open(catalog,"r") as infile:
        for line in infile:
            if line[0] is "#": continue
            parts = line.split()
            halos.append([float(parts[x_index]),float(parts[y_index]),float(parts[z_index])])
            Nh+=1
    halos = np.array(halos)

    #Read in the dm particles
    dm = pgr.readsnap(dmpath,"pos","dm")

    #Treecorr interface
    config = {'nbins':nbins,'min_sep':Rlimits[0],'max_sep':Rlimits[1]}
    halorandom_cat = treecorr.Catalog(x=halorandoms[:,0],
                                      y=halorandoms[:,1],
                                      z=halorandoms[:,2],
                                      config=config)
    dmrandom_cat   = treecorr.Catalog(x=dmrandoms[:,0],
                                      y=dmrandoms[:,1],
                                      z=dmrandoms[:,2],
                                      config=config)
    halo_cat       = treecorr.Catalog(x=halos[:,0],
                                      y=halos[:,1],
                                      z=halos[:,2],
                                      config=config)
    dm_cat         = treecorr.Catalog(x=dm[:,0],
                                      y=dm[:,1],
                                      z=dm[:,2],
                                      config=config)
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
    return Nh

def create_hmcf_directories(outpath):
    os.system("mkdir -p %s"%outpath+"/halomatter_correlation_function/full_hmcf")
    os.system("mkdir -p %s"%outpath+"/halomatter_correlation_function/JK_single")
    os.system("mkdir -p %s"%outpath+"/halomatter_correlation_function/JK_combined")
    os.system("mkdir -p %s"%outpath+"/halomatter_correlation_function/cov_matrix")
    os.system("mkdir -p %s"%outpath+"/halomatter_correlation_function/final_hmcf")
    print "\tHalo-matter correlation function directories created."
    return
