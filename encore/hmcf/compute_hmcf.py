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
    #GOTTA REWORK THE RANDOMS BECAUSE THEY AREN'T IN THE RIGHT PLACE

    print "Halo-matter correlation function not implemented yet!"
    return
