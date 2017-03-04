"""
The encore class.

There is a unit test at the bottom of this file.
"""
import os

class encore(object):
    def __init__(self,halopath='NOT INITIALIZED',dmpath="NOT INITIALIZED",
                 randompath=None,DSdmpath=None,outpath="./",
                 particle_mass=3e10,do_JK=False,ndivs=2,DSF=1000):
        self.particle_mass = particle_mass #Msun/h
        self.halopath = halopath
        self.dmpath = dmpath
        self.randompath = randompath
        self.DSdmpath = DSdmpath
        self.outpath = outpath
        self.do_JK = do_JK
        self.ndivs = ndivs
        self.DSF = DSF
        #Make a path for the info files
        os.system("mkdir -p %s"%outpath+"/info_files")

    def reduce_halo_catalogs(self):
        """
        Reduce the halo catalog.
        """
        import reduce_catalogs
        reduce_catalogs.reduce_halo_catalogs.reduce_halo_catalog(self.halopath,self.outpath,self.particle_mass,self.do_JK,self.ndivs)
        return

    def create_random_catalogs(self,edges,N,do_JK=False,do_DM=False,recreate=False):
        """Create catalogs of random points in a volume specified by the user.

        Args:
            edges (array_like): Spatial edges that contain the random points. Assumes a cube.
            N (int): Number of random points to be created in the entire volume.
            do_JK (bool): Flag to turn on partitioning into jackknife regions; default False.
            do_DM (bool): Flag to also create randoms for dark matter particles, in addition to halos; Default False.
            recreate (bool): Recreate random catalogs if they are found to exist already; default False.

        """
        import randoms
        if self.randompath is not None:
            print "Randompath specified as: %s"%self.randompath
            if not recreate: print "Using random_catalog file found there."
            else:
                print "Creating randoms there."
                randoms.create_random_catalogs.create_halo_random_catalog(self.randompath,edges,N,self.ndivs,do_DM)
            return
        else: self.randompath = self.outpath
        randoms.create_random_catalogs.create_halo_random_catalog(self.outpath,edges,N,self.ndivs,do_DM)
        return

    def compute_mass_function(self,nbins=10,do_JK=None):
        """Computes the halo mass function.

        Args:
            nbins (int): Number of mass bins to put halos in; default is 10.
            do_JK (bool): Flag to turn on partitioning into jackknife regions; default uses the value passed at initialization.
        
        """
        import mass_function
        if do_JK is None: do_JK = self.do_JK
        mass_function.compute_mass_function(self.outpath,nbins,do_JK,self.ndivs)
        return

    def compute_hhcf(self,edges,nbins=10,limits=[1.0,50.0],do_JK=None):
        """Compute the halo-halo correlation function.
        
        Args:
            edges (array_like): Spatial edges that contain the random points. Assumes a cube.
            nbins (int): Number of mass bins to put halos in; default is 10.
            limits (double): Radial limits of the bins of the correlation function; default is [1.0,50.0].
            do_JK (bool): Flag to turn on partitioning into jackknife regions; default uses the value passed at initialization.
        
        """
        import hhcf
        if do_JK is None: do_JK = self.do_JK
        hhcf.compute_hhcf(self.outpath,self.randompath,nbins,limits,edges,do_JK,self.ndivs)
        return

    def down_sample_dm(self,DSF=None):
        """Down sample the dark matter particles by a factor of DSF, which is short for "down sampling factor".

        Args:
            DSF (int): Factor by which the dark matter catalog is reduced by. E.g. 10 means only one out of ten particles is kept; default uses the value passed at initialization.

        """
        import down_sampling
        if DSF is None: DSF = self.DSF
        else: self.DSF = DSF
        if self.DSdmpath is not None:
            if os.path.exists(self.DSdmpath+"/down_sampled_dm/down_sampled_dm_DSF%d"%DSF):
                print "Down sampled DM catalog found in %s"%self.DSdmpath
                print "\tUsing that instead.\n\tDown sampling complete."
                return
            else:
                print "Down sampled DM path specified but no down sampled catalog found.\n\tCreating a down sampled catalog in %s"%outpath
                down_sampling.down_sampling.down_sample(self.outpath,
                                                        self.dmpath,DSF)
                return
        self.DSdmpath = self.outpath
        down_sampling.down_sampling.down_sample(self.DSdmpath,
                                                self.dmpath,DSF)
        return

    def jackknife_dm(self):
        """Jackknife the dark matter particles.

        """
        import down_sampling
        down_sampling.down_sampling.jackknife_dm(self.outpath,self.DSF,self.ndivs)
        return

    def compute_hmcf(self,edges,nbins=10,limits=[1.0,50.0],do_JK=None):
        """Compute the halo-matter correlation function.
        
        Args:
            edges (array_like): Spatial edges that contain the random points. Assumes a cube.
            nbins (int): Number of mass bins to put halos in; default is 10.
            limits (double): Radial limits of the bins of the correlation function; default is [1.0,50.0].
            do_JK (bool): Flag to turn on partitioning into jackknife regions; default uses the value passed at initialization.

        """
        import hmcf
        if do_JK is None: do_JK = self.do_JK
        hmcf.compute_hmcf(self.outpath,nbins,limits,edges,do_JK,self.ndivs,self.DSF,self.DSdmpath)
        return

if __name__=="__main__":
    particle_mass = 3.98769e10 #Msun/h
    my_encore = encore(outpath="../output/",particle_mass=particle_mass,do_JK=True)
    my_encore.reduce_halo_catalogs()
    my_encore.compute_mass_function(do_JK=True)
    edges = [0.0,1050.0] #Mpc/h; spatial edges of the snapshot
    my_encore.create_random_catalogs(edges,N=100000,do_DM=True)
    my_encore.compute_hhcf(edges,do_JK=True)
    my_encore.compute_hmcf(edges,do_JK=True)
    print "Unit test complete"
