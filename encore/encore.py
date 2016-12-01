"""
The encore class.
"""

class encore(object):
    def __init__(self,particle_mass,outpath="./",do_JK=False):
        self.particle_mass = particle_mass #Msun/h
        self.outpath = outpath
        self.do_JK = do_JK
        self.create_paths()

    def create_paths(self):
        """
        Create paths to the output directories.
        """
        import create_paths
        create_paths.create_paths(self.outpath)
        return

    def reduce_halo_catalogs(self,ndivs=4):
        """
        Reduce the halo catalog.
        """
        import reduce_catalogs
        reduce_catalogs.reduce_halo_catalog(self.outpath,\
                                            self.particle_mass,\
                                            self.do_JK,ndivs)
        #if do_JK: reduce_catalogs.jackknife_halo_catalog()
        return

    def compute_mass_function(self,nbins=10,do_JK=None,ndivs=4):
        """
        Compute the halo mass function.
        """
        import compute_mass_function
        if do_JK is None: do_JK = self.do_JK
        compute_mass_function.compute_mass_function(self.outpath,nbins,do_JK,ndivs)
        return

if __name__=="__main__":
    particle_mass = 3.98769e10 #Msun/h
    my_encore = encore(particle_mass,"./output/",do_JK=False)
    my_encore.reduce_halo_catalogs()
    my_encore.compute_mass_function(do_JK=True)
    print "Unit test complete"
