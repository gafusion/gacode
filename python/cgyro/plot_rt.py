import struct
import sys
import numpy as np
import os
from matplotlib import rc
import matplotlib.pyplot as plt
from gacodefuncs import *
from cgyro.data import cgyrodata

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

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species

if nx < 0:
   nx = nr+1
   usefft = True
else:
   usefft = False
   
epx = np.zeros([nx,nr],dtype=np.complex)
x = np.zeros([nx])

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
         epx[i,p]=np.exp(1j*(p-nr/2)*x[i])
      
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
   
   d = np.zeros([nx],dtype=np.complex)

   # Mapping
   # d[ ix,0 ] = c[ ix,0] 
   # d[ ix,0 ] = c[-ix,0]^* 

   for ix in range(-nr/2+1,nr/2):
      i = ix
      if ix < 0:
         i = ix+nx
      d[i] = c[-ix+nr/2]

   start = time.time()

   # Sign convention negative exponent exp(-inx)
   f = np.real(np.fft.fft(d))*0.5
    
   end = time.time()
  
   return f,end-start
#------------------------------------------------------------------------

u=specmap(sim.mass[species],sim.z[species])

# Set filename root and title
isfield = True
if (moment == 'n'):
    fdata = '.cgyro.kxky_n'
    title = r'${\delta \mathrm{n}}_'+u+'$'
    isfield = False
elif (moment == 'e'):
    fdata = '.cgyro.kxky_e'
    title = r'${\delta \mathrm{E}}_'+u+'$'
    isfield = False
elif (moment == 'phi'):
    fdata = '.cgyro.kxky_phi'
    title = r'$\delta\phi$'
elif (moment == 'apar'):
    fdata = '.cgyro.kxky_apar'
    title = r'$\delta A_\parallel$'
elif (moment == 'bpar'):
    fdata = '.cgyro.kxky_bpar'
    title = r'$\delta B_\parallel$'

# ERROR CHECKS

# Check to see if data exists (try binary data first)
if os.path.isfile('bin'+fdata):
    fdata = 'bin'+fdata
    print 'INFO: (plot_fluct) Found binary data in '+fdata 
    hasbin = True
elif os.path.isfile('out'+fdata):
    fdata = 'out'+fdata
    print ' BAD: (plot_fluct) Using inefficient ASCII data in '+fdata 
    hasbin = False
else:
    print 'ERROR: (plot_fluct) No data for -moment '+moment+' exists.  Try -moment phi'
    sys.exit()

if usefft:
    print 'INFO: (plot_fluct) Using FFT (fast)'

# **WARNING** Assumes theta_plot=1 
if isfield:
    n_chunk = 2*nr*nn
else:
    n_chunk = 2*nr*ns*nn
     
# Open binary file
fbin = open(fdata,'rb') 

f = np.zeros([nt,nx])

for i in range(nt):
   aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
   print 'INFO: (plot_fluct) Time index '+str(i) 

   if isfield:
      a = np.reshape(aa,(2,nr,nn),order='F')
      c = a[0,:,0]+1j*a[1,:,0]
   else:
      a = np.reshape(aa,(2,nr,ns,nn),order='F')
      c = a[0,:,species,0]+1j*a[1,:,species,0]

   # Derivative (d/dx)**nd   
   for p in range(-nr/2,nr/2):
      c[p+nr/2] = c[p+nr/2]*(1j*p)**nd
      
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

