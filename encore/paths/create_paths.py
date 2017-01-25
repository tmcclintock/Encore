"""
Create the paths to the output directories.
"""
import os

def create_paths(outpath):
    print "Creating directories."
    #Create the info files directories
    os.system("mkdir -p %s"%outpath+"/info_files")
    
    print "MM Correlation function directory not implemented yet!"
    
    print "DeltaSigma directories not implemented yet!"
    return
