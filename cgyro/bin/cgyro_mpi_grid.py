import os
import sys
import numpy as np
from gacodefuncs import *

def gcd(m,n):

   a = m ; b = n

   if a < b:
      c = a
      a = b
      b = c

   while True:
     c = np.mod(a,b)
     if c == 0:
        break
     a = b
     b = c

   d = b

   return d

n_xi       = int(parameter_val('input.cgyro.gen','N_XI'))
n_energy   = int(parameter_val('input.cgyro.gen','N_ENERGY'))
n_species  = int(parameter_val('input.cgyro.gen','N_SPECIES'))
n_radial   = int(parameter_val('input.cgyro.gen','N_RADIAL'))
n_theta    = int(parameter_val('input.cgyro.gen','N_THETA'))
n_toroidal = int(parameter_val('input.cgyro.gen','N_TOROIDAL'))

# Velocity-space (v) and configuration-space (c) dimensions
nv = n_energy*n_xi*n_species
nc = n_radial*n_theta

d = gcd(nv,nc)

print('          [coll]     [str]      [NL]  ')
print(' n_MPI    nc_loc    nv_loc   n_split  v_ord')
print('------    ------    ------   -------  -----')
for it in range(d*n_toroidal):
   j = it+1
   if np.mod(d*n_toroidal,j) == 0 and np.mod(j,n_toroidal) == 0:
      if np.mod(j,n_species*n_toroidal) == 0:
         extra = '   1,2'
      else:
         extra = '   1'

      n_proc_1 = j//n_toroidal
      nc_loc = nc//n_proc_1          
      nv_loc = nv//n_proc_1           
      nsplit = 1+(nv_loc*n_theta-1)//n_toroidal
      print(' {:5d}     {:5d}     {:5d}     {:5d}'.format(j,nc_loc,nv_loc,nsplit)+extra)


