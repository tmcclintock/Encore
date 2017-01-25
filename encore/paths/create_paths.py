"""
Create the paths to the output directories.
"""
import os

def create_paths(outpath):
    print "Creating directories."
    #Create the info files directories
    os.system("mkdir -p %s"%outpath+"/info_files")

    #Create the down sampled DM particle directory
    os.system("mkdir -p %s"%outpath+"/down_sampled_dm/")
    os.system("mkdir -p %s"%outpath+"/down_sampled_dm/JK_dm_cats")
    
    #Create the randoms directories
    os.system("mkdir -p %s"%outpath+"/randoms")
    print "\tRandom directories created."

    #Create the halo-halo correlation function output directories
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/full_hhcf")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/JK_single")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/JK_combined")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/cov_matrix")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/final_hhcf")
    print "\tHalo-halo correlation function directories created."

    #Create the halo-halo correlation function output directories
    os.system("mkdir -p %s"%outpath+"/halomatter_correlation_function/full_hmcf")
    os.system("mkdir -p %s"%outpath+"/halomatter_correlation_function/JK_single")
    os.system("mkdir -p %s"%outpath+"/halomatter_correlation_function/JK_combined")
    os.system("mkdir -p %s"%outpath+"/halomatter_correlation_function/cov_matrix")
    os.system("mkdir -p %s"%outpath+"/halomatter_correlation_function/final_hmcf")
    print "\tHalo-matter correlation function directories created."

    print "MM Correlation function directory not implemented yet!"
    
    print "DeltaSigma directories not implemented yet!"
    return
