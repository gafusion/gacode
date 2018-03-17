# Tool to repair old output formats

import os
import sys
import numpy as np
from gacodefuncs import *
from cgyro.data import cgyrodata


# Remove old parameters
os.system('cp input.cgyro tmp')
os.system('grep -v -e IMPLICIT_FLAG tmp > input.cgyro')

# Run cgyro -t
os.system('cgyro -t')

# Update time vector (old 2-column versus new 3-column)
data = np.loadtxt('out.cgyro.time')
nt = data.shape[0]
if (data.shape[1] == 2):
   print 'INFO: (fix) Found old 2-column time array'
   f=open('out.cgyro.time','w')
   for i in range(data.shape[0]):
      f.write('{0:.5e} {1:.5e} {1:.5e}\n'.format(data[i,0],data[i,1],data[i,1]))
   f.close()

# Read data in simdir
sim = cgyrodata('./')
    
ns = sim.n_species

# Old/large flux format
data = 0
sim.getbigflux()          
ys = np.zeros([ns,3,sim.n_field,sim.n_n,sim.n_time])
if hasattr(sim,'kxky_flux_n'):
    print 'INFO: (fix.py) Found out.cgyro.kxky_flux_n' 
    ys[:,0,0,:,:] = np.sum(sim.kxky_flux_n,axis=0)
    data=1
if hasattr(sim,'kxky_flux_e'):       
    print 'INFO: (fix.py) Found out.cgyro.kxky_flux_e' 
    ys[:,1,0,:,:] = np.sum(sim.kxky_flux_e,axis=0)
    data=1

if data==0:
    print 'Error: (fix.py) No flux data found.'   

f = open('out.cgyro.ky_flux','w') 
for n in range(sim.n_time):
    for m in range(sim.n_n):
        for k in range(sim.n_field):
            for j in range(3):
                for i in range(ns):
                    f.write(str(ys[i,j,k,m,n])+'\n')

f.close()

# Update to binary format
os.system('cgyro -bin')
# Remove ASCII
os.system('cgyro -cbin')
