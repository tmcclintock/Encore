"""
Compute the JK halo-halo correlation function
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

def calculate_JK_hhcf(outpath,nbins,limits,edges,Nh,randoms,ndivs):
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
    print "\tComputing JK HHCF."
    #Jackknife subregion step size
    step = (edges[1]-edges[0])/ndivs
    Njk = int(ndivs**3)

    #Read in all halos
    all_halos = read_halos(outpath,Njk)

    #Treecorr interface
    config = {'nbins':nbins,'min_sep':limits[0],'max_sep':limits[1]}
    
    #Calculate RR autocorrelation once
    random_cat = treecorr.Catalog(x=randoms[:,0],y=randoms[:,1],z=randoms[:,2],config=config)
    RRa = treecorr.NNCorrelation(config)
    RRa.process(random_cat)

    #Get all autocorrelations
    #These are all of length Njks
    DDa_all, DRa_all = calculate_autos(config,all_halos,randoms,step,ndivs,Njk)

    #Get all cross correlations
    DDc_all,DRc_all,RRc_all = calculate_cross(config,all_halos,randoms,step,ndivs,Njk)

    #Find the totals
    DDt,DRt,RRt = calculate_total(RRa,DDa_all,DRa_all,DDc_all,DRc_all,RRc_all,Njk,config)

    #Find each of the xi curves
    xi_all = calculate_xi_LOO(RRa,DDt,DRt,RRt,DDa_all,DRa_all,
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
    """
    Create the finalized data file with the errors.
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
    C = (Njk-1.0)/Njk
    for i in range(nbins):
        for j in range(nbins):
            cov[i,j] = C*np.sum((xi_true[i]-xi_all[:,i])*(xi_true[j]-xi_all[:,j]))
    covpath = outpath+"halohalo_correlation_function/cov_matrix/cov_matrix.txt"
    np.savetxt(covpath,cov)
    return cov

def calculate_xi_LOO(RRa,DDt,DRt,RRt,DDa_all,DRa_all,DDc_all,DRc_all,RRc_all,Njk):
    """
    Calculate the leave-one-out HHCFs.
    """
    xi_all = []
    for i in range(Njk):
        DDl = DDt.copy() #Copy the totals
        DRl = DRt.copy()
        RRl = RRt.copy()
        #Remove the auto-correlations
        DDl.weight[:]-=DDa_all[i].weight[:]
        DDl.npairs[:]-=DDa_all[i].npairs[:]
        DDl.tot-=DDa_all[i].tot
        DRl.weight[:]-=DRa_all[i].weight[:]
        DRl.npairs[:]-=DRa_all[i].npairs[:]
        DRl.tot-=DRa_all[i].tot
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
            RRl.weight[:]-=0.5*RRc_all[i].weight[:]
            RRl.npairs[:]-=0.5*RRc_all[i].npairs[:]
            RRl.tot-=0.5*RRc_all[i].tot
        #Calculate xi_LOO
        xi,varxi = DDl.calculateXi(RRl,DRl)
        xi_all.append(xi)
    print "\t\tLeave-one-out HHCFs calculated."
    return np.array(xi_all)

def calculate_total(RRa,DDa_all,DRa_all,DDc_all,DRc_all,RRc_all,Njk,config):
    """
    Resum the auto and cross correlations
    to get total quantities.
    """
    DDt = treecorr.NNCorrelation(config)
    DRt = treecorr.NNCorrelation(config)
    RRt = treecorr.NNCorrelation(config)
    for i in range(Njk):
        DDt+=DDa_all[i]
        DRt+=DRa_all[i]
        DRt+=DRc_all[i]
        RRt+=RRa
        DDt.weight[:]+=0.5*DDc_all[i].weight[:]
        DDt.npairs[:]+=0.5*DDc_all[i].npairs[:]
        DDt.tot+=0.5*DDc_all[i].tot
        RRt.weight[:]+=0.5*RRc_all[i].weight[:]
        RRt.npairs[:]+=0.5*RRc_all[i].npairs[:]
        RRt.tot+=0.5*RRc_all[i].tot
    print "\t\tResumming HHCF JK total complete."
    return DDt,DRt,RRt

def calculate_cross(config,all_halos,randoms,step,ndivs,Njk):
    """
    Calcualte the DD, DR and RR cross correlations
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
        halos1 = all_halos[index1]
        halo_cat1 = treecorr.Catalog(x=halos1[:,0],
                                     y=halos1[:,1],
                                     z=halos1[:,2],config=config)
        i1 = index1%ndivs
        j1 = (index1/ndivs)%ndivs
        k1 = index1/ndivs**2
        random_cat1 = treecorr.Catalog(x=randoms[:,0]+i1*step,
                                       y=randoms[:,1]+j1*step,
                                       z=randoms[:,2]+k1*step,
                                       config=config)
        for index2 in xrange(index1+1,Njk):
            halos2 = all_halos[index2]
            halo_cat2 = treecorr.Catalog(x=halos2[:,0],
                                         y=halos2[:,1],
                                         z=halos2[:,2],config=config)
            i2 = index2%ndivs
            j2 = (index2/ndivs)%ndivs
            k2 = index2/ndivs**2
            random_cat2 = treecorr.Catalog(x=randoms[:,0]+i2*step,
                                           y=randoms[:,1]+j2*step,
                                           z=randoms[:,2]+k2*step,
                                           config=config)
            DD = treecorr.NNCorrelation(config)
            DR = treecorr.NNCorrelation(config)
            RD = treecorr.NNCorrelation(config)
            RR = treecorr.NNCorrelation(config)
            DD.process(halo_cat1,halo_cat2)
            DR.process(halo_cat1,random_cat2)
            RD.process(halo_cat2,random_cat1)
            RR.process(random_cat1,random_cat2)
            DDc_all[index1]+=DD
            DDc_all[index2]+=DD
            DRc_all[index1]+=DR
            DRc_all[index2]+=RD
            RRc_all[index1]+=RR
            RRc_all[index2]+=RR
    return DDc_all,DRc_all,RRc_all

def calculate_autos(config,all_halos,randoms,step,ndivs,Njk):
    """
    Calculate the DD and DR autocorrelations
    """
    DDa_all = []
    DRa_all = []
    for index in range(Njk):
        halos = all_halos[index]
        i = index%ndivs
        j = (index/ndivs)%ndivs
        k = index/ndivs**2
        random_cat = treecorr.Catalog(x=randoms[:,0]+i*step,
                                      y=randoms[:,1]+j*step,
                                      z=randoms[:,2]+k*step,config=config)
        halo_cat = treecorr.Catalog(x=halos[:,0],
                                    y=halos[:,1],
                                    z=halos[:,2],config=config)
        DD = treecorr.NNCorrelation(config)
        DR = treecorr.NNCorrelation(config)
        DD.process(halo_cat)
        DR.process(halo_cat,random_cat)
        DDa_all.append(DD)
        DRa_all.append(DR)
    print "\t\tHHCF DD and DR autocorrelations computed."
    return DDa_all,DRa_all

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
