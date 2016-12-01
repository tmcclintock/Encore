"""
Create the paths to the output directories.
"""
import os

def create_paths(outpath):
    #Create the info files directories
    os.system("mkdir -p %s"%outpath+"/info_files")

    #Create general halo catalog directories
    os.system("mkdir -p %s"%outpath+"/reduced_halo_cats")
    os.system("mkdir -p %s"%outpath+"/JK_halo_cats")
    
    #Create the mass function output directories
    os.system("mkdir -p %s"%outpath+"/mass_function/full_N")
    os.system("mkdir -p %s"%outpath+"/mass_function/JK_single_N")
    os.system("mkdir -p %s"%outpath+"/mass_function/JK_combined_N")
    os.system("mkdir -p %s"%outpath+"/mass_function/cov_matrix")
    os.system("mkdir -p %s"%outpath+"/mass_function/final_mass_function")
    print "Mass function directories created."

    #Create the halo-halo correlation function output directories
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/randoms")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/full_hhcf")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/JK_single")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/JK_combined")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/cov_martrix")
    os.system("mkdir -p %s"%outpath+"/halohalo_correlation_function/final_hhcf")
    print "Halo-halo correlation function directories created."

    print "HM and MM Correlation function directories not implemented yet!"
    
    print "DeltaSigma directories not implemented yet!"
    return
