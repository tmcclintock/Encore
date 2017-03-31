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

    def set(self,**kwargs):
        """Set new attributes or replace old ones.
        """
        for key, value in kwargs.items():
            setattr(self, key, value)
        return

    def reduce_halo_catalogs(self,recreate=False):
        """
        Reduce the halo catalog by removing halos that are too small (<200 particles).

        Note: This assumes that the snapshot is a cubic region.

        Args:
            recreate (bool): Flag to re-reduce the rockstar catalog, even if it's already reduced; default False.
        """
        import reduce_catalogs
        cat   = getattr(self,"catalog")
        edges = getattr(self,"edges")
        args  = {"outpath":"./", "particle_mass":None, "do_JK":False, "ndivs":4}
        for key in args.keys():
            try:
                args[key] = getattr(self,key)
            except AttributeError: pass
        reduce_catalogs.reduce_halo_catalogs.reduce_halo_catalog(cat, args['outpath'], edges, args['particle_mass'], args['do_JK'], args['ndivs'], recreate)
        return

    def create_random_catalogs(self):
        """Create catalogs of random points in a volume specified by the user.
        """
        import randoms
        edges     = getattr(self,"edges")
        N_randoms = getattr(self,"N_randoms")
        args      = {"outpath":"./", "ndivs":2}
        for key in args.keys():
            try:
                args[key] = getattr(self,key)
            except AttributeError: pass
        randoms.create_random_catalogs.create_halo_random_catalog(args['outpath'], N_randoms, edges, args['ndivs'])
        return

    def down_sample_dm(self):
        """Down sample the dark matter particles by a factor of DSF, which is short for "down sampling fraction".

        Note: default down sample fraction is 1/1000, since anything larger than that has some chance of filling up the memory.
        """
        import down_sampling
        dmpath = getattr(self,"dmpath")
        args      = {"outpath":"./", "dm_down_sample_fraction":0.001}
        for key in args.keys():
            try:
                args[key] = getattr(self,key)
            except AttributeError: pass
        down_sampling.down_sampling.down_sample(args['outpath'], dmpath, args['dm_down_sample_fraction'])
        return

    def jackknife_dm(self):
        """Jackknife the down sampled dark matter particle catalog.
        """
        import down_sampling
        dmpath = getattr(self,"dmpath")
        edges = getattr(self,"edges")
        args      = {"outpath":"./", "ndivs":4}
        for key in args.keys():
            try:
                args[key] = getattr(self,key)
            except AttributeError: pass
        down_sampling.down_sampling.jackknife_dm(args['outpath'], dmpath, edges, args['ndivs'])
        return

    def compute_mass_function(self):
        """Computes the halo mass function.

        Note: jkcatalog must be a formatted string that takes a single integer representing which jackknife region it is.
        E.g.: jk_halo_catalog%d.txt
        """
        import mass_function
        cat  = getattr(self,"catalog")
        args = {"outpath":"./", "jkcatalog":None, "nbins":10, "do_JK":False, "ndivs":2}
        for key in args.keys():
            try:
                args[key] = getattr(self,key)
            except AttributeError: pass
        mass_function.compute_mass_function(cat,args['outpath'],args['jkcatalog'],args['nbins'],args['do_JK'],args['ndivs'])
        return

    def compute_hhcf(self):
        """Compute the halo-halo correlation function.
        """
        import hhcf
        cat   = getattr(self,"catalog")
        rands = getattr(self,"randoms")
        edges = getattr(self,"edges")
        args  = {"outpath":"./", "nbins":10, "Rlimits":[1.0,50.0], "do_JK":False, "jkcatalog":None, "jkrands":None, "ndivs":2}
        for key in args.keys():
            try:
                args[key] = getattr(self,key)
            except AttributeError: pass
        hhcf.compute_hhcf(args['outpath'], cat, rands, args['jkcatalog'], args['jkrands'], args['nbins'], args['Rlimits'], edges, args['do_JK'], args['ndivs'])
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
