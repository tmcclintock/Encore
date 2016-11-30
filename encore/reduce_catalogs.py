"""
Reduce the halo catalog.
"""
import os

def reduce_halo_catalog(outpath,pmass):
    print "Reducing halo catalog."
    redpath = outpath+"/reduced_halo_cats/reduced_halo_cat.txt"
    if os.path.exists(redpath): print "Reduced halo catalog already exists."
    else: 
        print "Not implemented yet."
    return

def jackknife_halo_catalog(outpath):
    return
