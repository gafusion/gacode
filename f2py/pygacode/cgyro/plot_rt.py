import struct
import sys
import numpy as np
import os
from matplotlib import rc
import matplotlib.pyplot as plt
from ..gacodefuncs import *
from .data import cgyrodata

PREC='f' ; BIT=4

# Use first 3 args to define plot and font size
rc('text',usetex=True)
rc('font',size=int(sys.argv[12]))

ftype = sys.argv[1]
moment = sys.argv[2]
species = int(sys.argv[3])
ymin = sys.argv[4]
ymax = sys.argv[5]
nx = int(sys.argv[6])
nd = int(sys.argv[7])
istr = sys.argv[8]
fmin = sys.argv[9]
fmax = sys.argv[10]
colormap = sys.argv[11]
font = int(sys.argv[12])
theta = int(sys.argv[13])

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species
nth = sim.theta_plot

if nx < 0:
   nx = nr+1
   usefft = True
else:
   usefft = False

epx = np.zeros([nx,nr],dtype=complex)
x = np.zeros([nx])

itheta,thetapi = indx_theta(theta,nth)

#------------------------------------------------------------------------
# Some setup
#
if usefft:
   for i in range(nx):
      x[i] = i*2*np.pi/nx
else:
   # Fourier arrays
   for i in range(nx):
      x[i] = i*2*np.pi/(nx-1)
      for p in range(nr):
         epx[i,p]=np.exp(1j*(p-nr//2)*x[i])

#------------------------------------------------------------------------
# Real-space field reconstruction
def maptoreal(nr,nx,c):

    import numpy as np
    import time

    start = time.time()

    f = np.zeros([nx])
    for p in range(nr):
       f[:] = f[:]+np.real(c[p]*epx[:,p])

    end = time.time()

    return f,end-start
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# FFT version
def maptoreal_fft(nr,nx,c):

   import numpy as np
   import time

   d = np.zeros([nx],dtype=complex)

   # Mapping
   # d[ ix,0 ] = c[ ix,0]
   # d[ ix,0 ] = c[-ix,0]^*

   for ix in range(-nr//2+1,nr//2):
      i = ix
      if ix < 0:
         i = ix+nx
      d[i] = c[-ix+nr//2]

   start = time.time()

   # Sign convention negative exponent exp(-inx)
   f = np.real(np.fft.fft(d))*0.5

   end = time.time()

   return f,end-start
#------------------------------------------------------------------------

# Get filename and tags
fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],moment)

# Check to see if data exists (try binary data first)
if os.path.isfile('bin'+fdata):
    fdata = 'bin'+fdata
    print('INFO: (plot_rt) Found binary data in '+fdata)
    hasbin = True
else:
   raise ValueError('(plot_rt) No data for -moment '+moment+' exists.  Try -moment phi')

if usefft:
    print('INFO: (plot_rt) Using FFT (fast)')

if isfield:
    n_chunk = 2*nr*nth*nn
else:
    n_chunk = 2*nr*nth*ns*nn
 
# Open binary file
with open(fdata,'rb') as fbin:

   f = np.zeros([nt,nx])

   for i in range(nt):
      aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
      print('INFO: (plot_fluct) Time index '+str(i)) 

      if isfield:
         a = np.reshape(aa,(2,nr,nth,nn),order='F')
         c = a[0,:,itheta,0]+1j*a[1,:,itheta,0]
      else:
         a = np.reshape(aa,(2,nr,nth,ns,nn),order='F')
         c = a[0,:,itheta,species,0]+1j*a[1,:,itheta,species,0]

      # Derivative (d/dx)**nd
      for p in range(-nr//2,nr//2):
         c[p+nr//2] = c[p+nr//2]*(1j*p)**nd

      ff = np.zeros([nx],order='F')
      if usefft:
         ff,t = maptoreal_fft(nr,nx,c)
      else:
         ff,t = maptoreal(nr,nx,c)

      f[i,:] = ff[:]

xp = x/(2*np.pi)*sim.length
tp = np.arange(nt)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
ax.set_title(title)
ax.set_xlabel(r'$(c_s/a) t$')
ax.set_ylabel(r'$x/\rho_s$')
if fmin == 'auto':
   f0=np.min(f)
   f1=np.max(f)
else:
   f0=float(fmin)
   f1=float(fmax)

levels = np.arange(f0,f1,(f1-f0)/256)
ax.contourf(tp,xp,np.transpose(f),levels,cmap=plt.get_cmap(colormap))
fig.tight_layout(pad=0.3)
plt.subplots_adjust(top=0.94)
plt.show()

