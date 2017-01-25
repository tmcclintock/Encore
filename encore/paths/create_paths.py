"""
Create the paths to the output directories.
"""
import os

def create_paths(outpath):
    print "Creating directories."
    #Create the info files directories
    os.system("mkdir -p %s"%outpath+"/info_files")
    
    #Create the randoms directories
    os.system("mkdir -p %s"%outpath+"/randoms")
    print "\tRandom directories created."

    #Create the halo-matter correlation function output directories

    print "MM Correlation function directory not implemented yet!"
    
    print "DeltaSigma directories not implemented yet!"
    return
