"""
Compute the JK halo-matter correlation function
and the covariance matrix.
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

def calculate_JK_hmcf(outpath, jkcatalog, jkdms, edges, nbins, Rlimits, halorandoms, dmrandoms, ndivs):
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
    print "\t\tReading in halos."
    all_halos = read_halos(jkcatalog,Njk)

    #Read in all dm
    print "\t\tReading in DM particles."
    all_dms = read_dm(jkdms,Njk)

    #Treecorr interface
    config = {'nbins':nbins,'min_sep':Rlimits[0],'max_sep':Rlimits[1]}

    #Calculate RR autocorrelation once
    print "\t\tPerforming RR autocorrelation."
    halorandom_cat = treecorr.Catalog(x=halorandoms[:,0],y=halorandoms[:,1],z=halorandoms[:,2],config=config)
    dmrandom_cat   = treecorr.Catalog(x=dmrandoms[:,0],y=dmrandoms[:,1],z=dmrandoms[:,2],config=config)
    RRa = treecorr.NNCorrelation(config)
    RRa.process(halorandom_cat,dmrandom_cat)

    #Get all autocorrelations
    #These are all of length Njks
    DDa_all, DRa_all, RDa_all = calculate_autos(config,all_halos,all_dms,halorandoms,dmrandoms,step,ndivs,Njk)
    sys.exit()

    #Get all cross correlations
    DDc_all,DRc_all,RDc_all,RRc_all = calculate_cross(config,all_halos,all_dms,halorandoms,dmrandoms,step,ndivs,Njk)

    #Find the totals
    DDt,DRt,RDt,RRt = calculate_total(RRa,DDa_all,DRa_all,RDa_all,DDc_all,DRc_all,RDc_all,RRc_all,Njk,config)

    #Find each of the xi curves
    xi_all = calculate_xi_LOO(RRa,DDt,DRt,RDt,RRt,DDa_all,DRa_all,RDa_all,
                              DDc_all,DRc_all,RDc_all,RRc_all,Njk)

    #Find the covariance matrix
    cov = calculate_cov_matrix(outpath,xi_all,nbins,Njk)
    err = np.sqrt(np.diagonal(cov))

    #Output the final data vector
    fullpath = outpath+"halomatter_correlation_function/full_hmcf/full_hmcf.txt"
    data = np.loadtxt(fullpath)
    R = data[:,0]
    xi_true = data[:,3]
    make_final_hmcf_data(outpath,R,xi_true,err)

    print "\tHMCF JK NOT IMPLEMENTED YET!!"
    return

def make_final_hmcf_data(outpath,R,xi,err):
    """
    Create the finalized data file with the errors.
    """
    finalpath = outpath+"/halomatter_correlation_function/final_hmcf/final_hmcf.txt"
    outfile = open(finalpath,"w")
    outfile.write("#R Mpc/h; Xi; Xi_JK_err\n")
    for i in range(len(R)):
        outfile.write("%.4e\t%.4e\t%.4e\n"%(R[i],xi[i],err[i]))
    outfile.close()
    print "\tFinal HMCF JK data created."
    return

def calculate_cov_matrix(outpath,xi_all,nbins,Njk):
    """
    Calculate the JK covariance matrix.
    """
    fullpath = outpath+"halomatter_correlation_function/full_hmcf/full_hmcf.txt"
    data = np.loadtxt(fullpath)
    xi_true = data[:,3]
    cov = np.zeros((nbins,nbins))
    C = (Njk-1.0)/Njk
    for i in range(nbins):
        for j in range(nbins):
            cov[i,j] = C*np.sum((xi_true[i]-xi_all[:,i])*(xi_true[j]-xi_all[:,j]))
    covpath = outpath+"halomatter_correlation_function/cov_matrix/cov_matrix.txt"
    np.savetxt(covpath,cov)
    return cov


def calculate_xi_LOO(RRa,DDt,DRt,RDt,RRt,DDa_all,DRa_all,RDa_all,DDc_all,DRc_all,RDc_all,RRc_all,Njk):
    """
    Calculate the leave-one-out HHCFs.
    """
    xi_all = []
    for i in range(Njk):
        DDl = DDt.copy() #Copy the totals
        DRl = DRt.copy()
        RDl = RDt.copy()
        RRl = RRt.copy()
        #Remove the auto-correlations
        DDl.weight[:]-=DDa_all[i].weight[:]
        DDl.npairs[:]-=DDa_all[i].npairs[:]
        DDl.tot-=DDa_all[i].tot
        DRl.weight[:]-=DRa_all[i].weight[:]
        DRl.npairs[:]-=DRa_all[i].npairs[:]
        DRl.tot-=DRa_all[i].tot
        RDl.weight[:]-=RDa_all[i].weight[:]
        RDl.npairs[:]-=RDa_all[i].npairs[:]
        RDl.tot-=RDa_all[i].tot
        RRl.weight[:]-=RRa.weight[:]
        RRl.npairs[:]-=RRa.npairs[:]
        RRl.tot-=RRa.tot
        for j in range(Njk): #Remove the cross-correlations
            DDl.weight[:]-=0.5*DDc_all[i].weight[:]
            DDl.npairs[:]-=0.5*DDc_all[i].npairs[:]
            DDl.tot-=0.5*DDc_all[i].tot
            DRl.weight[:]-=DRc_all[i].weight[:]
            DRl.npairs[:]-=DRc_all[i].npairs[:]
            DRl.tot-=DRc_all[i].tot
            RDl.weight[:]-=RDc_all[i].weight[:]
            RDl.npairs[:]-=RDc_all[i].npairs[:]
            RDl.tot-=RDc_all[i].tot
            RRl.weight[:]-=0.5*RRc_all[i].weight[:]
            RRl.npairs[:]-=0.5*RRc_all[i].npairs[:]
            RRl.tot-=0.5*RRc_all[i].tot
        #Calculate xi_LOO
        xi,varxi = DDl.calculateXi(RRl,DRl,RDl)
        xi_all.append(xi)
    print "\t\tLeave-one-out HMCFs calculated."
    return np.array(xi_all)

def calculate_total(RRa,DDa_all,DRa_all,RDa_all,DDc_all,DRc_all,RDc_all,RRc_all,Njk,config):
    """
    Resum the auto and cross correlations
    to get total quantities.
    """
    DDt = treecorr.NNCorrelation(config)
    DRt = treecorr.NNCorrelation(config)
    RDt = treecorr.NNCorrelation(config)
    RRt = treecorr.NNCorrelation(config)
    for i in range(Njk):
        DDt+=DDa_all[i]
        DRt+=DRa_all[i]
        DRt+=DRc_all[i]
        RDt+=RDa_all[i]
        RDt+=RDc_all[i]
        RRt+=RRa
        DDt.weight[:]+=0.5*DDc_all[i].weight[:]
        DDt.npairs[:]+=0.5*DDc_all[i].npairs[:]
        DDt.tot+=0.5*DDc_all[i].tot
        RRt.weight[:]+=0.5*RRc_all[i].weight[:]
        RRt.npairs[:]+=0.5*RRc_all[i].npairs[:]
        RRt.tot+=0.5*RRc_all[i].tot
    print "\t\tResumming HMCF JK total complete."
    return DDt,DRt,RDt,RRt

def calculate_cross(config,all_halos,all_dms,halorandoms,dmrandoms,step,ndivs,Njk):
    """
    Calcualte the DD, DR, RD and RR cross correlations
    """
    #Inititalize these objects
    DDc_all = []
    DRc_all = []
    RDc_all = []
    RRc_all = []
    for i in xrange(0,Njk):
        DDc_all.append(treecorr.NNCorrelation(config))
        DRc_all.append(treecorr.NNCorrelation(config))
        RDc_all.append(treecorr.NNCorrelation(config))
        RRc_all.append(treecorr.NNCorrelation(config))

    for index1 in xrange(0,Njk):
        halos = all_halos[index1]
        halo_cat = treecorr.Catalog(x=halos[:,0],
                                    y=halos[:,1],
                                    z=halos[:,2],
                                    config=config)
        i1 = index1%ndivs
        j1 = (index1/ndivs)%ndivs
        k1 = index1/ndivs**2
        dmrandom_cat = treecorr.Catalog(x=dmrandoms[:,0]+i1*step,
                                        y=dmrandoms[:,1]+j1*step,
                                        z=dmrandoms[:,2]+k1*step,
                                        config=config)
        for index2 in xrange(index1+1,Njk):
            dm = all_dms[index2]
            dm_cat = treecorr.Catalog(x=dm[:,0],
                                      y=dm[:,1],
                                      z=dm[:,2],config=config)
            i2 = index2%ndivs
            j2 = (index2/ndivs)%ndivs
            k2 = index2/ndivs**2
            halorandom_cat = treecorr.Catalog(x=halorandoms[:,0]+i2*step,
                                              y=halorandoms[:,1]+j2*step,
                                              z=halorandoms[:,2]+k2*step,
                                              config=config)
            DD = treecorr.NNCorrelation(config)
            DR = treecorr.NNCorrelation(config)
            RD = treecorr.NNCorrelation(config)
            RR = treecorr.NNCorrelation(config)
            DD.process(halo_cat,dm_cat)
            DR.process(halo_cat,halorandom_cat)
            RD.process(dm_cat,dmrandom_cat)
            RR.process(dmrandom_cat,halorandom_cat)
            DDc_all[index1]+=DD
            DDc_all[index2]+=DD
            DRc_all[index1]+=DR
            RDc_all[index2]+=RD
            RRc_all[index1]+=RR
            RRc_all[index2]+=RR
    return DDc_all,DRc_all,RDc_all,RRc_all

def calculate_autos(config,all_halos,all_dms,halorandoms,dmrandoms,step,ndivs,Njk):
    """
    Calculate the DD, DR, RD autocorrelations
    """
    DDa_all = []
    DRa_all = []
    RDa_all = []
    for index in range(Njk):
        i = index%ndivs
        j = (index/ndivs)%ndivs
        k = index/ndivs**2
        halos = all_halos[index]
        dm = all_dms[index]

        halorandom_cat = treecorr.Catalog(x=halorandoms[:,0]+i*step,
                                          y=halorandoms[:,1]+j*step,
                                          z=halorandoms[:,2]+k*step,
                                          config=config)
        halo_cat = treecorr.Catalog(x=halos[:,0],
                                    y=halos[:,1],
                                    z=halos[:,2],
                                    config=config)
        dmrandom_cat = treecorr.Catalog(x=dmrandoms[:,0]+i*step,
                                        y=dmrandoms[:,1]+j*step,
                                        z=dmrandoms[:,2]+k*step,
                                        config=config)
        dm_cat = treecorr.Catalog(x=dm[:,0],
                                  y=dm[:,1],
                                  z=dm[:,2],
                                  config=config)
        DD = treecorr.NNCorrelation(config)
        DR = treecorr.NNCorrelation(config)
        RD = treecorr.NNCorrelation(config)
        DD.process(halo_cat,dm_cat)
        DR.process(halo_cat,halorandom_cat)
        RD.process(dm_cat,dmrandom_cat)
        DDa_all.append(DD)
        DRa_all.append(DR)
        RDa_all.append(RD)
    print "\t\tHMCF DD, DR, RD autocorrelations computed."
    return DDa_all,DRa_all,RDa_all
    
def read_halos(jkcatalog,Njk):
    """
    Read in halos from the jackknife files.
    Returns an array of Njk X N_halos_i X 3 where
    there are N_halos_i in the i'th JK file.
    This is not a constant number.
    """
    all_halos = []
    for index in xrange(0,Njk):
        infile = open(jkcatalog%index,"r")
        halos = [] #Will be Nhjk X 3
        for line in infile:
            if line[0] is "#": continue
            parts = line.split()
            halos.append([float(parts[x_index]),float(parts[y_index]),float(parts[z_index])])
        halos = np.array(halos)
        infile.close()
        all_halos.append(halos)
    return np.array(all_halos)

def read_dm(jkdms,Njk):
    """
    Read in dms from the jackknife files.
    Returns an array of Njk X N_dm_i X 3 where
    there are N_dm_i in the i'th JK file.
    This is not a constant number.
    """

    all_dms = []
    for index in xrange(0,Njk):
        dms = np.genfromtxt(jkdms%index)
        all_dms.append(dms)
    all_dms = np.array(all_dms)
    return all_dms

