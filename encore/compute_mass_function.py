"""
Compute the halo mass function.
"""
import os
import numpy as np

def compute_mass_function(outpath,nbins,limits=None):
    #Step 1: figure out the Min/Max masses
    if limits is None: limits = find_mass_limits(outpath)

    #Step 2: calcalate the full mass function
    calculate_full_mass_function(outpath,limits,nbins)
    
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
    return [mmin,mmax]

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
