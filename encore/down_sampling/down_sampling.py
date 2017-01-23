"""
Perform the dark matter down sampling.

Note: DSF stands for "down sampling factor"
and is the multiplicative factor that the number
of dark matter particles will be reduced by. So if
DSF is 10 then only one tenth of the particles
will remain.
"""
import os, inspect
dirname = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")
try: import pygadgetreader as pgr
except ImportError: raise Exception("Must install pygadgetreader.")


def down_sample(outpath,dmpath,DSF):
    print "Down sampling on file: %s"%dmpath
    print "\tDSF = %d"%DSF

    Ndm = pgr.readheader(dmpath,"dmcount")
    Nds = Ndm/DSF

    DSpath = outpath+"/down_sampled_dm/down_sampled_dm_DSF%d"%DSF
    if os.path.exists(DSpath): print "Down sampled DM catalog already exists."
    else:
        command = dirname+"/subsample_code/subsamp_parts LGADGET %s %d %s"%(dmpath,Nds,DSpath)
        os.system(command)

    print "\tDown sampling complete."
    return

def jackknife_dm(outpath,DSF,ndivs):
    """
    Jackknife the dark matter particles from the down-sampled
    dark matter catalog.
    """
    print "Jackknifing DM particles."
    if os.path.exists(outpath+"/info_files/spatial_limits.txt"):
        limits = np.loadtxt(outpath+"/info_files/spatial_limits.txt")
    else:
        raise Exception("Must find spatial limits before jackknifing dark matter particles.")

    dx = (limits[0,1]-limits[0,0])/ndivs
    dy = (limits[1,1]-limits[1,0])/ndivs
    dz = (limits[2,1]-limits[2,0])/ndivs

    jkoutbase = outpath+"/down_sampled_dm/JK_dm_cats/jk_halo_cat_%d.txt"
    jkarray = []
    Njks = ndivs**3
    for i in xrange(0,Njks): jkarray.append(open(jkoutbase%i,"w"))

    #Read in the dm particles
    dmpath = outpath+"/down_sampled_dm/down_sampled_dm_DSF%d"%DSF
    dm = pgr.readsnap(dmpath,"pos","dm")
    
    for i in xrange(0,len(dm)):
        x,y,z = dm[i]
        i = min(np.floor(x/dx),ndivs-1)
        j = min(np.floor(y/dy),ndivs-1)
        k = min(np.floor(z/dz),ndivs-1)
        index = int(k*ndivs*ndivs + j*ndivs + i)
        jkarray[index].write("%e\t%e\t%e\n"%(x,y,z))

    for i in range(Njks): jkarray[i].close()
    print "\tJackknifing DM particles complete."
    return
