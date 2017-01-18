"""
Perform the dark matter down sampling.

Note: DSF stands for "down sampling factor"
and is the multiplicative factor that the number
of dark matter particles will be reduced by. So if
DSF is 10 then only one tenth of the particles
will remain.
"""
try: import numpy as np
except ImportError: raise Exception("Must install numpy.")
try: import pygadgetreader as pgr
except ImportError: raise Exception("Must install pygadgetreader.")

def down_sample(dmpath,DSF):
    print "Down sampling on file: %s"%dmpath
    print "\tDSF = %d"%DSF

    Ndm = pgr.readheader(dmpath,"dmcount")
    Nds = Ndm/DSF

    

    return
