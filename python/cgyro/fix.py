# Convert large dlux datafiles to new compact form (no radial dependence)

import sys
import numpy as np
from gacodeplotdefs import *
from gacodefuncs import *
from cgyro.data import cgyrodata

sim = cgyrodata('./')
sim.getflux()

if hasattr(sim,'ky_flux'):
    print "INFO: (fix.py) Detected out.cgyro.ky_flux.  You're up to date."
    sys.exit()
    
ns = sim.n_species

data = 0

# Old/large flux format
sim.getbigflux()          
ys = np.zeros([ns,3,sim.n_field,sim.n_n,sim.n_time])
if hasattr(sim,'kxky_flux_n'):
    print 'INFO: (fix.py) Found out.cgyro.kxk_flux_n' 
    ys[:,0,0,:,:] = np.sum(sim.kxky_flux_n,axis=0)
    data=1
if hasattr(sim,'kxky_flux_e'):       
    print 'INFO: (fix.py) Found out.cgyro.kxk_flux_e' 
    ys[:,1,0,:,:] = np.sum(sim.kxky_flux_e,axis=0)
    data=1

if data==0:
    print 'Error: (fix.py) No flux data found.'   

file = open('out.cgyro.ky_flux','w') 
for n in range(sim.n_time):
    for m in range(sim.n_n):
        for k in range(sim.n_field):
            for j in range(3):
                for i in range(ns):
                    file.write(str(ys[i,j,k,m,n])+'\n')

print 'INFO: (fix.py) Created out.cgyro.ky_flux' 
