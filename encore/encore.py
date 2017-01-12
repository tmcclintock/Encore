"""
The encore class.

There is a unit test at the bottom of this file.
"""

class encore(object):
    def __init__(self,outpath="./",particle_mass=3e10,do_JK=False,ndivs=2):
        self.particle_mass = particle_mass #Msun/h
        self.outpath = outpath
        self.do_JK = do_JK
        self.ndivs = ndivs
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
        reduce_catalogs.reduce_halo_catalog(self.outpath,self.particle_mass,self.do_JK,self.ndivs)
        return

    def create_random_catalogs(self,edges,N,do_JK=False,do_DM=False):
        """
        Create random catalogs.
        """
        import create_random_catalogs as crc
        crc.create_halo_random_catalog(self.outpath,edges,N,self.ndivs)
        print "Only halo randoms implemented right now!"
        return

    def compute_mass_function(self,nbins=10,do_JK=None):
        """
        Compute the halo mass function.
        """
        import mass_function
        if do_JK is None: do_JK = self.do_JK
        mass_function.compute_mass_function(self.outpath,nbins,do_JK,self.ndivs)
        return

    def compute_hhcf(self,edges,nbins=10,limits=[1.0,50.0],do_JK=None):
        """
        Compute the halo-halo correlation function.
        
        The only required input is the spatial edges
        of the snapshot (e.g. 0 to 1000 Mpc/h).
        """
        import hhcf
        if do_JK is None: do_JK = self.do_JK
        hhcf.compute_hhcf(self.outpath,nbins,limits,edges,do_JK,self.ndivs)
        return

if __name__=="__main__":
    particle_mass = 3.98769e10 #Msun/h
    my_encore = encore("../output/",particle_mass,do_JK=True)
    my_encore.reduce_halo_catalogs()
    my_encore.compute_mass_function(do_JK=True)
    edges = [0.0,1050.0] #Mpc/h; spatial edges of the snapshot
    my_encore.create_random_catalogs(edges,N=100000)
    my_encore.compute_hhcf(edges,do_JK=True)
    print "Unit test complete"
