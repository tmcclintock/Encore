"""
Compute the halo mass function.
"""
import os
import numpy as np

def compute_mass_function(outpath,nbins,do_JK,ndivs,limits=None):
    #Step 1: figure out the Min/Max masses
    if limits is None: 
        if os.path.exists(outpath+"/info_files/mass_limits.txt"):
            limits = np.loadtxt(outpath+"/info_files/mass_limits.txt")
        else: 
            limits = find_mass_limits(outpath)
            np.savetxt(outpath+"/info_files/mass_limits.txt",limits)

    #Step 2: calcalate the full mass function
    calculate_full_mass_function(outpath,limits,nbins)

    #Step 3: do JK calculation
    if do_JK:
        calculate_JK_mass_function(outpath,limits,nbins,ndivs)
    
    print "Mass function not implemented yet."
    return

def find_mass_limits(outpath):
    mmin, mmax = 1e999, 0.0
    redpath = outpath+"/reduced_halo_cats/reduced_halo_cat.txt"
    infile = open(redpath,"r")
    for line in infile:
        if line[0] is "#": continue
        parts = line.split()
        m = float(parts[2])
        if m < mmin: mmin = m
        if m > mmax: mmax = m*1.00001 #Go slightly high for completeness
    infile.close()
    return np.array([mmin,mmax])

def calculate_full_mass_function(outpath,limits,nbins):
    redpath = outpath+"/reduced_halo_cats/reduced_halo_cat.txt"
    fullpath = outpath+"/mass_function/full_N/full_N.txt"
    N = np.zeros((nbins)) #Holds the mass function
    mmin,mmax = limits
    lmmin,lmmax = np.log(limits)
    edges = np.exp(np.linspace(lmmin,lmmax,nbins+1))
    bins = np.array([edges[:-1],edges[1:]]).T
    infile = open(redpath,"r")
    for line in infile:
        if line[0] is "#": continue
        parts = line.split()
        m = float(parts[2])
        for b,i in zip(bins,range(len(bins))):
            if m >= b[0] and m < b[1]: 
                N[i]+=1
                break
            else: continue #Halo not in any bins
    outfile = open(fullpath,"w")
    outfile.write("#Bin_left\tBin_right\tN_halos\n")
    for i in range(len(bins)):
        outfile.write("%.4e\t%.4e\t%d\n"%(bins[i,0],bins[i,1],N[i]))
    outfile.close()
    print "Successfully created full N(M)."
    return

def calculate_JK_mass_function(outpath,limits,nbins,ndivs):
    mmin,mmax = limits
    lmmin,lmmax = np.log(limits)
    edges = np.exp(np.linspace(lmmin,lmmax,nbins+1))
    bins = np.array([edges[:-1],edges[1:]]).T

    Njks = ndivs**3
    jkoutbase = outpath+"/mass_function/JK_single_N/jk_single_N_%d.txt"
    for i in range(Njks): 
        outfile = open(jkoutbase%i,"w")
        outfile.write("#Bin_left\tBin_right\tN_halos\n")
    
        N = np.zeros((nbins)) #Holds the mass function
        jkredpath = outpath+"/JK_halo_cats/jk_halo_cat_%d.txt"%i
        infile = open(jkredpath,"r")
        for line in infile:
            if line[0] is "#": continue
            parts = line.split()
            m = float(parts[2])
            for b,j in zip(bins,range(len(bins))):
                if m >= b[0] and m < b[1]: 
                    N[j]+=1
                    break
                else: continue #Halo not in any bins
        for j in range(len(bins)):
            outfile.write("%.4e\t%.4e\t%d\n"%(bins[j,0],bins[j,1],N[j]))
        outfile.close()
    print "Successfully created JK N(M) singles."
