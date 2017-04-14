"""
Compute the halo-halo correlation function.
"""
import os,sys
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")
try: 
    from Corrfunc.theory.DD import DD
    from Corrfunc.theory.xi import xi
except ImportError: raise Exception("Must install Corrfunc.")


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

def compute_hhcf(outpath, halopath, jkhalopath, jkrandompath, nbins, limits, edges, do_JK, ndivs):
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

    #Step 1: calculate the full HH correlation function
    #calcalate_hhcf_full(outpath, halopath, nbins, limits, edges)

    #Step 2: calculate the HH correlation function within JK subregions
    if do_JK:
        print "HHCF JK Not implemented yet"
        import compute_hhcf_jk
        compute_hhcf_jk.calculate_JK_hhcf(outpath, jkhalopath, jkrandompath, nbins, limits, edges, ndivs)
        return

    print "Halo-halo correlation function complete."
    return

def calcalate_hhcf_full(outpath, halopath, nbins, limits, edges):
    """
    Calculate the halo-halo correlation function
    for the full volume.
    """
    X, Y, Z, M = np.genfromtxt(halopath, unpack=True)[[x_index,y_index,z_index, m_index]]
    N = len(M)
    Mmean = np.mean(M)
    print "\tHalo masses averaged, results:"
    print "\t\tMmean  = %e"%Mmean
    print "\t\tNhalos = %d"%N
    boxsize = max(edges) - min(edges)
    nthreads = 4 #arbitrary
    bins = np.logspace(np.log10(min(limits)), np.log10(max(limits)), nbins+1)
    results = xi(boxsize, nthreads, bins, X, Y, Z, output_ravg=True)
    header = "rmin rmax ravg xi npairs weightavg"
    np.savetxt(outpath+"/halohalo_correlation_function/full_hhcf/full_hhcf.txt", results, header=header)
    print "\tFull halo-halo correlation function complete."
    return

def create_hhcf_directories(outpath):
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/full_hhcf")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/JK_single")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/JK_combined")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/cov_matrix")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/final_hhcf")
    print "\tHalo-halo correlation function directories created."
