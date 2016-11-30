"""
The encore class.
"""

class encore(object):
    def __init__(self,particle_mass,outpath="./"):
        self.particle_mass = particle_mass #Msun/h
        self.outpath = outpath
        self.create_paths()

    def create_paths(self):
        """
        Create paths to the output directories.
        """
        import create_paths
        create_paths.create_paths(self.outpath)
        return

    def reduce_halo_catalogs(self):
        """
        Reduce the halo catalog.
        """
        import reduce_catalogs
        reduce_catalogs.reduce_halo_catalog(self.outpath,self.particle_mass)
        #if do_JK: reduce_catalogs.jackknife_halo_catalog()
        return

    def compute_mass_function(self,nbins=10):
        """
        Compute the halo mass function.
        """
        import compute_mass_function
        compute_mass_function.compute_mass_function(self.outpath,nbins)
        return

if __name__=="__main__":
    particle_mass = 3.98769e10 #Msun/h
    my_encore = encore(particle_mass,"./output/")
    my_encore.reduce_halo_catalogs()
    my_encore.compute_mass_function()
    print "Unit test complete"
