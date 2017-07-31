# Dump text report of flux

import sys
import numpy as np
from gacodefuncs import *
from cgyro.data import cgyrodata

w = float(sys.argv[1])
ext = sys.argv[2]

sim = cgyrodata('./')
sim.getflux()

nt = sim.n_time

# Determine imin
imin=iwindow(sim.t,w)

if imin == nt-1:
    print "Averaging Window too small." 
    sys.exit()

print 'INFO: Average Window:',str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])

if ext == 'dump':

   # Datafile output

   b = np.zeros([sim.n_n])

   tag = ['fluxn','fluxe','fluxv']

   for i in range(3):
      f=open('out.cgyroplot.'+tag[i],'w')
      for ispec in range(sim.n_species):
         f.write(str(ispec)+' ')
         for i_n in range(sim.n_n):
            y0 = np.sum(sim.ky_flux,axis=2)
            y  = y0[:,i,:,:]
            b[i_n] = average(y[ispec,i_n,:],sim.t,w)
            f.write("{:.3e}".format(b[i_n])+' ')
            f.write('\n')
      print 'Wrote output to out.cgyroplot.'+tag[i]
else:
    # Screen output

    # Set print precision
    np.set_printoptions(precision=3,suppress=False,threshold=100000)

    b = np.zeros([sim.n_species])
 
    tag = [
       'GAMMA [GB]',
       'Q     [GB]',
       'PI    [GB]']

    for i in range(3):
       try:
          for ispec in range(sim.n_species):
             y = np.sum(sim.ky_flux,axis=(2,3))
             b[ispec] = average(y[ispec,i,:],sim.t,w)
             print tag[i],b
       except:
             pass

    
