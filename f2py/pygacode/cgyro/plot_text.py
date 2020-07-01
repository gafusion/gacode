# Dump text report of flux

import sys
import numpy as np
from ..gacodefuncs import *
from .data import cgyrodata

def print_freq():
   # Set print precision
   np.set_printoptions(precision=5,suppress=True)

   print('   omega    gamma')
   print(sim.freq[:,0,-1])
   
def print_flux(ascii):
   sim.getflux('auto')

   b = np.zeros([sim.n_species])
 
   tag = [
      'GAMMA [GB]',
      'Q     [GB]',
      'PI    [GB]']


   # Determine imin
   imin,imax=iwindow(sim.t,w,wmax)

   print('INFO: (text.py) Average Window:'+str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[imax]))

   title = '        '
   for ispec in range(sim.n_species):
      title = title+'       '+specmap(sim.mass[ispec],sim.z[ispec])
   print(title)

   y = np.sum(sim.ky_flux,axis=(2,3))
   for i in range(3):
      bstr=''
      for ispec in range(sim.n_species):
         b[ispec] = average(y[ispec,i,:],sim.t,w,wmax)
         bstr = bstr+"{:7.3f}".format(b[ispec])+' '
      print(tag[i]+' '+bstr)

   if ext == 'ascii':
      try:
         import gnuplotlib as gp
         x = sim.t
         y0 = y[0,1,:]
         gp.plot(x,y0,_with='lines',terminal='dumb 120,40')
      except:
         print('HINT: pip install --user gnuplotlib') 
         
#-------------------------------------------------------------------
        
w    = float(sys.argv[1])
wmax = float(sys.argv[2])
ext  = sys.argv[3]

sim = cgyrodata('./')

nt = sim.n_time

if sim.n_n > 1:
   # Print flux (assuming nonlinear case)
   print_flux(ext)
else:
   print_flux(ext)
   # Pring frequency (assuming linear)
   print_freq()
        
