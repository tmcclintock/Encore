"""
Compute the JK halo-halo correlation function
and the covariance matrix.
"""
import os,sys
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")
#try: import treecorr
#except ImportError: raise Exception("Must install treecorr.")
try: 
    from Corrfunc.theory.DD import DD as DDcf
    from Corrfunc.theory.xi import xi as xicf
    nthreads = 8
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
#A constant we use
sqrt3 = np.sqrt(3)

def calculate_JK_hhcf(outpath, halopath, randompath, nbins, limits, edges, ndivs):
    """Calculate the halo-halo correlation function for the JK subregions.

    Inputs:
       outpath: the base directory where the output directories exist
       halopath: the path to the halo catalog
       nbins: the number of radial bins to calculate xi_hh for
       limits: an array with two entries with the min/max separation
       edges: the spatial edges of the simulation, i.e. xmin and xmax
       Nh: number of halos
       randoms: catalog of randoms within a JK region
       ndivs: number of JK subregions
    """
    print "\tComputing JK HHCF."
    print "This is broken right now."
    return 0

    #Jackknife subregion step size
    step = (edges[1]-edges[0])/ndivs
    Njk = int(ndivs**3)

    #Figure out the bins
    bins = np.logspace(np.log10(limits[0]), np.log10(limits[1]), nbins+1)

    #Read in all halos
    all_halos = read_jk_cats(halopath,Njk)
    all_rands = read_jk_cats(randompath,Njk)

    #Get all autocorrelations
    #These are all of length Njks
    DDa_all, DRa_all, RRa_all = calculate_autos(outpath, all_halos, all_rands, bins, Njk)
    print "Got autos"
    sys.exit()

    #Get all cross correlations
    DDc_all, DRc_all, RRc_all = calculate_cross(outpath,config,all_halos,all_rands,step,ndivs,Njk)

    #Find the totals
    DDt, DRt, RRt = calculate_total(outpath,DDa_all,DRa_all,RRa_all,DDc_all,DRc_all,RRc_all,Njk,config)

    #Find each of the xi curves
    xi_all = calculate_xi_LOO(DDt,DRt,RRt,DDa_all,DRa_all,RRa_all,
                              DDc_all,DRc_all,RRc_all,Njk)

    #Find the covariance matrix
    cov = calculate_cov_matrix(outpath,xi_all,nbins,Njk)
    err = np.sqrt(np.diagonal(cov))

    #Output the final data vector
    fullpath = outpath+"halohalo_correlation_function/full_hhcf/full_hhcf.txt"
    data = np.loadtxt(fullpath)
    R = data[:,0]
    xi_true = data[:,3]
    make_final_hhcf_data(outpath,R,xi_true,err)

    print "\tHHCF JK complete."
    return

def make_final_hhcf_data(outpath,R,xi,err):
    """Create the finalized data file with the errors.
    """
    finalpath = outpath+"/halohalo_correlation_function/final_hhcf/final_hhcf.txt"
    outfile = open(finalpath,"w")
    outfile.write("#R Mpc/h; Xi; Xi_JK_err\n")
    for i in range(len(R)):
        outfile.write("%.4e\t%.4e\t%.4e\n"%(R[i],xi[i],err[i]))
    outfile.close()
    print "\tFinal HHCF JK data created."
    return
    
def calculate_cov_matrix(outpath,xi_all,nbins,Njk):
    """
    Calculate the JK covariance matrix.
    """
    fullpath = outpath+"halohalo_correlation_function/full_hhcf/full_hhcf.txt"
    data = np.loadtxt(fullpath)
    xi_true = data[:,3]
    cov = np.zeros((nbins,nbins))
    C = (Njk-1.0)/Njk #Hartlap correction
    for i in range(nbins):
        for j in range(nbins):
            cov[i,j] = C*np.sum((xi_true[i]-xi_all[:,i])*(xi_true[j]-xi_all[:,j]))
        #if i < 5:
        #    print xi_true[i], xi_all[:10,i]
    covpath = outpath+"/halohalo_correlation_function/cov_matrix/cov_matrix.txt"
    np.savetxt(covpath,cov)
    return cov

def add_with_mult(NN1, NN2, x):
    NN1.meanr  += x * NN2.meanr
    NN1.meanlogr += x * NN2.meanlogr
    NN1.weight += x * NN2.weight
    NN1.npairs += x * NN2.npairs
    NN1.tot    += x * NN2.tot
    return

def subtract_with_mult(NN1, NN2, x):
    NN1.meanr  -= x * NN2.meanr
    NN1.meanlogr -= x * NN2.meanlogr
    NN1.weight -= x * NN2.weight
    NN1.npairs -= x * NN2.npairs
    NN1.tot    -= x * NN2.tot
    return

def calculate_xi_LOO(DDt,DRt,RRt,DDa_all,DRa_all,RRa_all,DDc_all,DRc_all,RRc_all,Njk):
    """Calculate the leave-one-out HHCFs.
    """
    xi_all = []
    for i in range(Njk):
        DDl = DDt.copy() #Copy the totals
        DRl = DRt.copy()
        RRl = RRt.copy()
        #Remove the auto-correlations
        subtract_with_mult(DDl, DDa_all[i], 1.0)
        subtract_with_mult(DRl, DRa_all[i], 1.0)
        subtract_with_mult(RRl, RRa_all[i], 1.0)

        #Remove cross correlations
        subtract_with_mult(DDl, DDc_all[i], 1.0)
        subtract_with_mult(DRl, DRc_all[i], 1.0)
        subtract_with_mult(RRl, RRc_all[i], 1.0)
            
        #Calculate xi_LOO
        xi,varxi = DDl.calculateXi(RRl,DRl)
        xi_all.append(xi)
    print "\t\tLeave-one-out HHCFs calculated."
    return np.array(xi_all)

def calculate_total(outpath,DDa_all,DRa_all,RRa_all,DDc_all,DRc_all,RRc_all,Njk,config):
    """Resum the auto and cross correlations to get total quantities.
    """
    DDt = treecorr.NNCorrelation(config)
    DRt = treecorr.NNCorrelation(config)
    RRt = treecorr.NNCorrelation(config)
    for i in range(Njk):
        DDt+=DDa_all[i]
        DRt+=DRa_all[i]
        RRt+=RRa_all[i]
        add_with_mult(DDt, DDc_all[i], 0.5)
        add_with_mult(DRt, DRc_all[i], 0.5)
        add_with_mult(RRt, RRc_all[i], 0.5)
    DDt.write(outpath+"/halohalo_correlation_function/JK_combined/hhcf_resum.txt",RRt,DRt)
    print "\t\tResumming HHCF JK total complete."
    return DDt,DRt,RRt

def calculate_cross(outpath,config,all_halos,all_rands,step,ndivs,Njk):
    """Calcualte the DD, DR and RR cross correlations.
    """
    #Inititalize these objects
    DDc_all = []
    DRc_all = []
    RRc_all = []
    for i in xrange(0,Njk):
        DDc_all.append(treecorr.NNCorrelation(config))
        DRc_all.append(treecorr.NNCorrelation(config))
        RRc_all.append(treecorr.NNCorrelation(config))
        
    for index1 in xrange(0,Njk):
        i1 = index1%ndivs
        j1 = (index1/ndivs)%ndivs
        k1 = index1/ndivs**2
        halos1 = all_halos[index1]
        rands1 = all_rands[index1]
        halo_cat1 = treecorr.Catalog(x=halos1[:,0],
                                     y=halos1[:,1],
                                     z=halos1[:,2],
                                     config=config)
        rand_cat1 = treecorr.Catalog(x=rands1[:,0],
                                     y=rands1[:,1],
                                     z=rands1[:,2],
                                     config=config)
        for index2 in xrange(index1+1,Njk):
            i2 = index2%ndivs
            j2 = (index2/ndivs)%ndivs
            k2 = index2/ndivs**2
            #Figure out if the separation is large enough to skip
            max_sep = config['max_sep']
            if max_sep < step*(np.sqrt((i1-i2)**2+(j1-j2)**2+(k1-k2)**2)-sqrt3):
                continue #Skip this cross correlation
            #print "Cross correlating %d vs %d"%(index1, index2)
            halos2 = all_halos[index2]
            rands2 = all_rands[index2]
            halo_cat2 = treecorr.Catalog(x=halos2[:,0],
                                         y=halos2[:,1],
                                         z=halos2[:,2],
                                         config=config)
            rand_cat2 = treecorr.Catalog(x=rands2[:,0],
                                         y=rands2[:,1],
                                         z=rands2[:,2],
                                         config=config)
            DD = treecorr.NNCorrelation(config)
            DR = treecorr.NNCorrelation(config)
            RD = treecorr.NNCorrelation(config)
            RR = treecorr.NNCorrelation(config)
            #DD.process(halo_cat1,halo_cat2)
            #DR.process(halo_cat1,rand_cat2)
            #RD.process(rand_cat1,halo_cat2)
            #RR.process(rand_cat1,rand_cat2)
            DD.process_cross(halo_cat1, halo_cat2)
            DR.process_cross(halo_cat1, rand_cat2)
            RD.process_cross(rand_cat1, halo_cat2)
            RR.process_cross(rand_cat1, rand_cat2)
            DDc_all[index1]+=DD
            DDc_all[index2]+=DD
            DRc_all[index1]+=DR
            DRc_all[index1]+=RD
            DRc_all[index2]+=DR
            DRc_all[index2]+=RD
            RRc_all[index1]+=RR
            RRc_all[index2]+=RR
    #Write them all out here
    for i in xrange(0, Njk):
        DDc_all[i].write(outpath+"/halohalo_correlation_function/JK_single/CROSS%d.txt"%i,RRc_all[i],DRc_all[i])
    return DDc_all, DRc_all, RRc_all

def calculate_autos(outpath, all_halos, all_rands, bins, Njk):
    """Calculate the DD and DR autocorrelations.
    This calculates the correlations within the same spatial regions.
    """
    DDa_all = []
    DRa_all = []
    RRa_all = []
    for index in xrange(0,Njk):
        xh, yh, zh = all_halos[index].T
        xr, yr, zr = all_rands[index].T
        print xh.shape, yh.shape, zh.shape
        print bins
        DDresult = DDcf(1, 4, bins, xh, yh, zh, periodic=False, verbose=True)
        DRresult = DDcf(0, nthreads, bins, xh, yh, zh,
                        X2=xr, Y2=yr, Z2=zr, 
                        periodic=True, verbose=False)
        RRresult = DDcf(1, nthreads, bins, xr, yr, zr,
                        periodic=True, verbose=False)
        
        DD = []
        DR = []
        RR = []
        for i in range(len(DDresult)):
            DD.append(DDresult[i]['npairs'])
            DR.append(DRresult[i]['npairs'])
            RR.append(RRresult[i]['npairs'])
        DD = np.array(DD)
        DR = np.array(DR)
        RR = np.array(RR)
        DDa_all.append(DD)
        DRa_all.append(DR)
        RRa_all.append(RR)
        DD.write(outpath+"/halohalo_correlation_function/JK_single/AUTO%d.txt"%index,RR,DR)
    print "\t\tHHCF DD, DR and RR autocorrelations computed."
    return DDa_all, DRa_all, RRa_all

def read_jk_cats(halopath, Njk):
    """Read in halos from the jackknife files.
    Returns an array of Njk X N_halos_i X 3 where
    there are N_halos_i in the i'th JK file.
    This is not a constant number.
    """
    all_halos = []
    for index in range(Njk):
        infile = open(halopath%index,"r")
        halos = [] #Will be Nhjk X 3
        for line in infile:
            if line[0] is "#": continue
            parts = line.split()
            halos.append([float(parts[x_index]),float(parts[y_index]),float(parts[z_index])])
        halos = np.array(halos)
        infile.close()
        all_halos.append(halos)
    return np.array(all_halos)
