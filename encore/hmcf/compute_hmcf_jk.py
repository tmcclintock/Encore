"""
Compute the JK halo-matter correlation function
and the covariance matrix.
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

def calculate_JK_hmcf(outpath,nbins,limits,edges,Nh,halorandoms,dmrandoms,ndivs,DSF):
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
    print "\tComputing JK HMCF."
    #Jackknife subregion step size
    step = (edges[1]-edges[0])/ndivs
    Njk = int(ndivs**3)

    #Read in all halos
    all_halos = read_halos(outpath,Njk)

    #Read in all dm
    all_dms = read_dm(outpath,Njk,DSF)

    #Treecorr interface
    config = {'nbins':nbins,'min_sep':limits[0],'max_sep':limits[1]}

    print "\tHMCF JK NOT IMPLEMENTED YET!!"
    return

    
def read_halos(outpath,Njk):
    """
    Read in halos from the jackknife files.
    Returns an array of Njk X N_halos_i X 3 where
    there are N_halos_i in the i'th JK file.
    This is not a constant number.
    """
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

def read_dm(outpath,Njk,DFS):
    """
    Read in the dark matter particles,
    and then relable into jackknife subregions.
    
    """
    dmpath = outpath+"/down_sampled_dm/down_sampled_dm_DSF%d"%DSF
    dm = pgr.readsnap(dmpath,"pos","dm")

    all_dms = []
    for index in xrange(0,Njk):
        
