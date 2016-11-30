"""
The encore class.
"""
import os

class encore(object):
    def __init__(self,outpath="./"):
        self.outpath = outpath
        self.make_paths()

    def make_paths(self):
        """
        Create paths to the output directories.
        """
        outpath = self.outpath
        #Create general halo catalog directories
        os.system("mkdir -p %s"%outpath+"/reduced_halo_cats")
        os.system("mkdir -p %s"%outpath+"/JK_halo_cats")

        #Create the mass function output directories
        os.system("mkdir -p %s"%outpath+"/mass_function/full_N")
        os.system("mkdir -p %s"%outpath+"/mass_function/JK_single_N")
        os.system("mkdir -p %s"%outpath+"/mass_function/JK_combined_N")
        os.system("mkdir -p %s"%outpath+"/mass_function/cov_matrices")
        os.system("mkdir -p %s"%outpath+"/mass_function/final_mass_functions")
        print "Mass function directories created."
        
        print "Correlation function directories not implemented yet!"

        print "DeltaSigma directories not implemented yet!"
        return

if __name__=="__main__":
    my_encore = encore("./output/")
    print "Unit test complete"
