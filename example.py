"""
An example of how to use encore.

Steps:
0) Add a path to encore to your python path.

1) Define the DM particle mass and the spatial edges of the snapshot.

2) Create an Encore object with a path to the output, the particle mass,
and (optional) specify whether you want jackknifing.

3) Reduce the halo catalog to halos you want.

4) Compute things you want, such as the mass function or
correlation functions. To accomplish the latter you
need to create random catalogs.
"""
import sys
sys.path.insert(0,"../Encore/")

import encore

#Step 1
particle_mass = 3.98769e10 #Msun/h
edges = [0.0,1050.0] #Mpc/h; spatial edges of the snapshot

#Step 2
my_encore = encore.encore("./output/",particle_mass,do_JK=True)

#Step 3
my_encore.reduce_halo_catalogs()

#Step 4
my_encore.compute_mass_function(do_JK=True)
my_encore.create_random_catalogs(edges,N=100000)
my_encore.compute_hhcf(edges,do_JK=True)
