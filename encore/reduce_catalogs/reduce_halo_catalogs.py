"""
Reduce the halo catalog.
"""
import os
import numpy as np

#Pull out the indices
indices = {}
with open("rockstar_config") as myfile:
    for line in myfile:
        name, var = line.partition("=")[::2]
        indices[name.strip()] = int(var)
x_index = indices['x']
y_index = indices['y']
z_index = indices['z']
m_index = indices['m']
if 'n' in indices:
    n_index = indices['n']
else:
    n_index = None
    print "\tWARNING: Particle number index missing."

def reduce_halo_catalog(halopath,outpath,pmass,do_JK,ndivs):
    print "Reducing halo catalog."
    redpath = outpath+"/reduced_halo_cats/reduced_halo_cat.txt"
    if os.path.exists(redpath): print "Reduced halo catalog already exists."
    else: 
        outfile = open(redpath,"w")
        with open(halopath) as infile:
            if n_index is not None:
                for line in infile:
                    if line[0] is "#": 
                        outfile.write(line)
                        continue
                    parts = line.split()
                    Np = int(parts[n_index])
                    if Np >= 200:outfile.write(line)
            else:
                for line in infile:
                    if line[0] is "#": 
                        outfile.write(line)
                        continue
                    parts = line.split()
                    M = float(parts[m_index])
                    Np = int(M/pmass+0.01) #extra added for rounding
                    if Np >= 200:outfile.write(line)
        outfile.close()
    if do_JK: jackknife_halo_catalog(outpath,ndivs)
    return

def jackknife_halo_catalog(outpath,ndivs):
    if os.path.exists(outpath+"/info_files/spatial_limits.txt"):
        limits = np.loadtxt(outpath+"/info_files/spatial_limits.txt")
    else:
        limits = find_spatial_limits(outpath)
        np.savetxt(outpath+"/info_files/spatial_limits.txt",limits)
    dx = (limits[0,1]-limits[0,0])/ndivs
    dy = (limits[1,1]-limits[1,0])/ndivs
    dz = (limits[2,1]-limits[2,0])/ndivs

    redpath = outpath+"/reduced_halo_cats/reduced_halo_cat.txt"
    jkoutbase = outpath+"/JK_halo_cats/jk_halo_cat_%d.txt"
    jkarray = []
    Njks = ndivs**3
    for i in range(Njks): jkarray.append(open(jkoutbase%i,"w"))
    
    infile = open(redpath,"r")
    for line in infile:
        if line[0] is "#": #Write headers
            for i in range(Njks): 
                jkarray[i].write(line)
            continue
        parts = line.split()
        x,y,z = float(parts[x_index]),float(parts[y_index]),float(parts[z_index])
        i = min(np.floor(x/dx),ndivs-1)
        j = min(np.floor(y/dy),ndivs-1)
        k = min(np.floor(z/dz),ndivs-1)
        index = int(k*ndivs*ndivs + j*ndivs + i)
        jkarray[index].write(line)
            
    for i in range(Njks): jkarray[i].close()
    print "JK halo catalogs complete."
    return

def find_spatial_limits(outpath):
    xmin,xmax = 1e99,-1e99
    ymin,ymax = zmin,zmax = xmin,xmax
    redpath = outpath+"/reduced_halo_cats/reduced_halo_cat.txt"
    infile = open(redpath,"r")
    for line in infile:
        if line[0] is "#": continue
        parts = line.split()
        x,y,z = float(parts[x_index]),float(parts[y_index]),float(parts[z_index])
        if x < xmin: xmin = np.floor(x)
        if x > xmax: xmax = np.ceil(x)
        if y < ymin: ymin = np.floor(y)
        if y > ymax: ymax = np.ceil(y)
        if z < zmin: zmin = np.floor(z)
        if z > zmax: zmax = np.ceil(z)
    infile.close()
    return np.array([[xmin,xmax],[ymin,ymax],[zmin,zmax]])
