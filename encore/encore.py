"""
The encore class.
"""
import os

class encore(object):
    """The encore class. Used to analyze rockstar halo catalogs along with accompanying dark matter particle catalogs.
    """
    def __init__(self,**kwargs):
        """Create a rockstar object.
        """
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self):
        outstr = "Rockstar encore with:"
        for name in self.__dict__:
            outstr += "\n\t%s: %s"%(name,getattr(self,name))
        return outstr

    def reduce_halo_catalogs(self,recreate=False):
        """
        Reduce the halo catalog by removing halos that are too small (<200 particles).

        Args:
            recreate (bool): Flag to re-reduce the rockstar catalog, even if it's already reduced; default False.
        """
        import reduce_catalogs
        cat = getattr(self,"catalog")
        args = {"outpath":"./", "particle_mass":None, "do_JK":False, "ndivs":2}
        for key in args.keys():
            try:
                args[key] = getattr(self,key)
            except AttributeError: pass
        reduce_catalogs.reduce_halo_catalogs.reduce_halo_catalog(cat,args['outpath'],args['particle_mass'],args['do_JK'],args['ndivs'],recreate)
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

    def compute_mass_function(self):
        """Computes the halo mass function.

        Note: jkcatalog must be a formatted string that takes a single integer representing which jackknife region it is.
        """
        import mass_function
        cat = getattr(self,"catalog")
        args = {"outpath":"./", "jkcatalog":None, "nbins":10, "do_JK":False, "ndivs":2}
        for key in args.keys():
            try:
                args[key] = getattr(self,key)
            except AttributeError: pass
        mass_function.compute_mass_function(cat,args['outpath'],args['jkcatalog'],args['nbins'],args['do_JK'],args['ndivs'])
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
    my_encore = encore(outpath="../output/",particle_mass=3e10,do_JK=True)
    print my_encore
