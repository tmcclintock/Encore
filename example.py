"""
An example of how to use encore.

Steps:
0) Add a path to encore to your python path.

1) Define the DM particle mass and the spatial edges of the snapshot.

2) Create an Encore object with a path to the output, the particle mass,
and (optional) specify whether you want jackknifing.

3) Reduce the halo catalog to halos you want, downsample dark matter
particles, jackknife the dark matter particles.

4) Compute things you want, such as the mass function or
correlation functions. To accomplish the latter you
need to create random catalogs.
"""

#Step 0
import sys
sys.path.insert(0,"../Encore/")
import encore

#Step 1
particle_mass = 3.98769e10 #Msun/h
edges = [0.0,1050.0] #Mpc/h; spatial edges of the snapshot

#Step 2
my_encore = encore.encore(outpath="./output/", 
                          particle_mass=particle_mass,do_JK=True)

#Step 3
#my_encore.down_sample_dm() #Already done
#my_encore.jackknife_dm()

#Step 4
#my_encore.compute_mass_function(do_JK=True)
#my_encore.create_random_catalogs(edges,N=800000) #Comment this out once it is run one time
my_encore.compute_hhcf(edges,do_JK=True)
#You can also define your own radial bins
#Note: more bins and smaller scales means a longer run time
#limits = [0.1,50.0]
#nbins = 20
#my_encore.compute_hmcf(edges,nbins=nbins,limits=limits,do_JK=False)
