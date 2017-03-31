"""
Perform the dark matter down sampling.

Note: DSF stands for "down sampling fraction"
and is the multiplicative factor that the number
of dark matter particles will be reduced by. So if
DSF is 0.1 then only one tenth of the particles
will remain.
"""
import os, inspect
dirname = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")
try: import pygadgetreader as pgr
except ImportError: raise Exception("Must install pygadgetreader.")


def down_sample(outpath, dmpath, DSF):
    #Create the directories to put the down-sampling output
    os.system("mkdir -p %s"%outpath+"/down_sampled_dm/")

    print "Down sampling on file: %s"%dmpath
    print "\tDSF = %f"%DSF

    Ndm = pgr.readheader(dmpath,"dmcount")
    Nds = int(Ndm*DSF) #Number of particles to keep

    DSpath = outpath+"/down_sampled_dm/down_sampled_dm_DSF%.2f"%DSF
    if os.path.exists(DSpath): print "Down sampled DM catalog already exists."
    else:
        command = dirname+"/subsample_code/subsamp_parts LGADGET %s %d %s"%(dmpath,Nds,DSpath)
        os.system(command)

    print "\tDown sampling complete."
    return

def jackknife_dm(outpath, dmpath, edges, ndivs):
    """
    Jackknife the dark matter particles from the down-sampled
    dark matter catalog.
    """
    #Create the output directory for this function
    os.system("mkdir -p %s"%outpath+"/down_sampled_dm/JK_dm_cats")

    print "Jackknifing DM particles."
    dx = dy = dz = (edges[1] - edges[0])/ndivs

    jkoutbase = outpath+"/down_sampled_dm/JK_dm_cats/jk_dm_cat_%d.txt"
    jkarray = []
    Njks = ndivs**3
    for i in xrange(0,Njks): jkarray.append(open(jkoutbase%i,"w"))

    #Read in the dm particles
    dmpath = dmpath
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
