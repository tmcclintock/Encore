"""
Compute the halo mass function.
"""
import os
import numpy as np

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

def compute_mass_function(catalog,outpath,jkcatalog,nbins,do_JK,ndivs):
    print "Computing mass function."

    #Step 0: create the paths
    create_mass_function_directories(outpath)

    #Step 1: figure out the Min/Max masses
    limits = find_mass_limits(catalog)
    print "\tUsing M_min = %.2e and M_max = %.2e"%(limits[0],limits[1])

    #Step 2: calcalate the full mass function
    calculate_full_mass_function(catalog, outpath, limits, nbins)

    #Step 3: do JK calculation
    if do_JK and jkcatalog is None: 
        raise Exception("Cannot jackknife the mass function without a jkcatalog. Please see the docstring for encore.compute_mass_function().")
    if do_JK: calculate_JK_mass_function(jkcatalog, outpath, limits, nbins, ndivs)
    
    print "\tMass function successfully computed."
    return

def find_mass_limits(catalog):
    mmin, mmax = 1e999, 0.0
    infile = open(catalog,"r")
    for line in infile:
        if line[0] is "#": continue
        parts = line.split()
        m = float(parts[m_index])
        if m < mmin: mmin = m
        if m > mmax: mmax = m*1.00001 #Go slightly high for completeness
    infile.close()
    return np.array([mmin,mmax])

def calculate_full_mass_function(catalog, outpath, limits, nbins):
    print "\tFinding full mass function."
    fullpath = outpath+"/mass_function/full_N/full_N.txt"
    N = np.zeros((nbins)) #Holds the mass function
    mmin,mmax = limits
    lmmin,lmmax = np.log(limits)
    edges = np.exp(np.linspace(lmmin ,lmmax, nbins+1))
    bins = np.array([edges[:-1], edges[1:]]).T
    infile = open(catalog, "r")
    for line in infile:
        if line[0] is "#": continue
        parts = line.split()
        m = float(parts[m_index])
        for b,i in zip(bins, range(len(bins))):
            if m >= b[0] and m < b[1]: 
                N[i]+=1
                break
            else: continue #Halo not in any bins
    outfile = open(fullpath, "w")
    outfile.write("#Bin_left\tBin_right\tN_halos\n")
    for i in range(len(bins)):
        outfile.write("%.4e\t%.4e\t%d\n"%(bins[i,0], bins[i,1], N[i]))
    outfile.close()
    return

def calculate_JK_mass_function(jkcatalog, outpath, limits, nbins, ndivs):
    print "\tFinding JK mass function."
    calculate_JK_mass_function_singles(jkcatalog, outpath, limits, nbins, ndivs)
    N = combine_JK_mass_function(outpath,nbins,ndivs)
    make_data_and_cov(outpath,nbins,ndivs,N)
    return

def make_data_and_cov(outpath ,nbins, ndivs, N):
    """
    Make the final data and the covariance matrix.
    """
    Njks = ndivs**3
    jkinbase  = outpath+"/mass_function/JK_combined_N/jk_combo_N_%d.txt"
    finalpath = outpath+"/mass_function/final_mass_function/final_mass_function.txt"
    covpath   = outpath+"/mass_function/cov_matrix/cov_matrix.txt"
    Nall = np.zeros((Njks,nbins))
    for i in range(Njks):
        jkcombo = np.genfromtxt(jkinbase%i)
        Nall[i] = jkcombo[:,2]
        bins = jkcombo[:,:2]
    Nmean = np.mean(Nall,0)
    cov = np.zeros((nbins,nbins))
    covfile = open(covpath,"w")
    for i in range(nbins):
        for j in range(nbins):
            cov[i,j] = (Njks-1.)/Njks*np.sum((Nall[:,i]-Nmean[i])*(Nall[:,i]-Nmean[i]))
            covfile.write("%e\t"%cov[i,j])
        covfile.write("\n")
    outfile = open(finalpath,"w")
    outfile.write("#Bin_left\tBin_right\tN_halos\tN_halos_err\n")
    for i in range(nbins):
        outfile.write("%.4e\t%.4e\t%d\t%e\n"%(bins[i,0],bins[i,1],N[i],np.sqrt(cov[i,i])))
    outfile.close()
    print "\tFinal mass function JK data and covariance matrix created."
    return

def combine_JK_mass_function(outpath, nbins, ndivs):
    """
    Combine the JKs.
    """
    Njks = ndivs**3
    jkinbase  = outpath+"/mass_function/JK_single_N/jk_single_N_%d.txt"
    jkoutbase = outpath+"/mass_function/JK_combined_N/jk_combo_N_%d.txt"
    N = np.zeros((nbins))
    Nerr = np.zeros((nbins))
    Nall = np.zeros((Njks,nbins))
    for i in range(Njks):
        jksingle = np.genfromtxt(jkinbase%i)
        Nall[i] = jksingle[:,2]
        bins = jksingle[:,:2]
    N = np.sum(Nall,0)
    for i in range(Njks):
        Nout = N - Nall[i]
        outfile = open(jkoutbase%i,"w")
        outfile.write("#Bin_left\tBin_right\tN_halos\n")
        for j in range(nbins):
            outfile.write("%.4e\t%.4e\t%d\n"%(bins[j,0],bins[j,1],Nout[j]))
        outfile.close()
    print "\tSuccessfully combined JK N(M) files."
    return N

def calculate_JK_mass_function_singles(jkcatalog, outpath, limits, nbins, ndivs):
    """
    Get N(M) for individual JK regions.
    """
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
        infile = open(jkcatalog%i,"r")
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
    print "\tJK mass function singles created."
    return

def create_mass_function_directories(outpath):
    os.system("mkdir -p %s"%outpath+"/mass_function/full_N")
    os.system("mkdir -p %s"%outpath+"/mass_function/JK_single_N")
    os.system("mkdir -p %s"%outpath+"/mass_function/JK_combined_N")
    os.system("mkdir -p %s"%outpath+"/mass_function/cov_matrix")
    os.system("mkdir -p %s"%outpath+"/mass_function/final_mass_function")
    print "\tMass function directories created."
    return
