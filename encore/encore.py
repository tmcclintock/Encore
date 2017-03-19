"""
The encore class.

There is a unit test at the bottom of this file.
"""
import os

class encore(object):
    """The encore class. Used to analyze rockstar halo catalogs along with accompanying dark matter particle catalogs.
    
    """

    def __init__(self,halopath='NOT INITIALIZED',dmpath="NOT INITIALIZED",
                 randompath=None,DSdmpath=None,reducedhalopath=None,
                 outpath="./",particle_mass=3e10,do_JK=False,ndivs=2,DSF=1000):
        """Create the encore object.

        Args:
            halopath (string; optional): TODO
            dmpath (string; optional): TODO
            randompath (string; optional): TODO
            DSdmpath (string; optional): TODO
            reducedhalopath (string; optional): TODO
            outpath (string; optional): Path for all output files; default is ./.
            particle_mass (float; optional): Mass of dark matter particles; default is 3e10 Msun/h; note that units are arbitrary.
            do_JK (boolean; optional): Compute jackknifes for all quantities, including jackknifing the DM and halo catalogs; default is False.
            ndivs (int; optional): TODO
            DSF (int; optional): TODO

        """
        self.outpath = outpath #Path to the output
        self.particle_mass = particle_mass #Msun/h
        self.halopath = halopath #Unreduced halos
        self.dmpath = dmpath #Undown sampled dark matter
        self.do_JK = do_JK #Boolean to do jackknifing
        self.ndivs = ndivs #Number of JK divisions in a dimension
        self.DSF = DSF #Down sampling factor
        if DSdmpath is not None: self.DSdmpath = DSdmpath
        else: self.DSdmpath = self.outpath
        if os.path.exists(self.DSdmpath+"/down_sampled_dm/down_sampled_dm_DSF%d"%DSF):
            print "Down sampled DM catalog found in %s"%self.DSdmpath
            print "\tDown sampling complete."
            self.down_sampled = True
        else: self.down_sampled = False
        if randompath is not None: self.randompath = randompath
        else: self.randompath = self.outpath
        if os.path.exists(self.randompath+"/randoms/"):
            print "Random catalogs found in %s"%self.randompath
            print "\tFor new cats: create_random_catalogs(...,recreate=True)."
            self.have_randoms = True
        else: self.have_randoms = False
        if reducedhalopath is not None: self.reducedhalopath = reducedhalopath
        else: self.reducedhalopath = self.outpath
        if os.path.exists(self.reducedhalopath+"/reduced_halo_cats/reduced_halo_cat.txt"):
            print "Reduced halo catalog found in %s"%self.reducedhalopath
            print "\tFor new cats: reduce_halo_catalogs(...,recreate=True)."
            self.have_reducedhalos = True
        else: self.have_reducedhalos = False
        #Make a path for the info files
        os.system("mkdir -p %s"%outpath+"/info_files/")

    def reduce_halo_catalogs(self,recreate=False):
        """
        Reduce the halo catalog.
        """
        import reduce_catalogs
        if not self.have_reducedhalos or recreate:
            reduce_catalogs.reduce_halo_catalogs.reduce_halo_catalog(self.halopath,self.reducedhalopath,self.particle_mass,self.do_JK,self.ndivs)
            self.have_reducedhalos = True
        else: print "Halos already reduced."
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
        if not self.have_randoms or recreate:
            randoms.create_random_catalogs.create_halo_random_catalog(self.randompath,edges,N,self.ndivs)
            self.have_randoms = True
        else: print "Random catalogs already created."
        return

    def down_sample_dm(self,DSF=None):
        """Down sample the dark matter particles by a factor of DSF, which is short for "down sampling factor".

        Args:
            DSF (int): Factor by which the dark matter particles is reduced by. E.g. 10 means only one tenth of the particles are kept; default is value passed at initialization.
        """
        import down_sampling
        if DSF is None: DSF = self.DSF
        else: self.DSF = DSF
        if not self.down_sampled:
            down_sampling.down_sampling.down_sample(self.DSdmpath,
                                                    self.dmpath,DSF)
            self.down_sampled = True
        else: print "Already down sampled."
        return

    def jackknife_dm(self):
        """Jackknife the down sampled dark matter particle catalog.

        """
        import down_sampling
        down_sampling.down_sampling.jackknife_dm(self.DSdmpath,
                                                 self.DSF,self.ndivs)
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
        hhcf.compute_hhcf(self.outpath,self.reducedhalopath,self.randompath,nbins,limits,edges,do_JK,self.ndivs)
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
        hmcf.compute_hmcf(self.outpath,self.DSdmpath,self.randompath,
                          nbins,limits,edges,do_JK,self.ndivs,self.DSF)
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
