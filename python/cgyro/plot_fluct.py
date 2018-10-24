import struct
import sys
import numpy as np
import os
from matplotlib import rc
import matplotlib.pyplot as plt
from gacodefuncs import *
from cgyro.data import cgyrodata

try:
   import gapy
   hasgapy = True
except:
   print ' BAD: (plot_fluct) Please build gapy.so library!'
   hasgapy = False

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
ny = int(sys.argv[7])
istr = sys.argv[8]
fmin = sys.argv[9]
fmax = sys.argv[10]
colormap = sys.argv[11]

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species

if nx < 0 or ny < 0:
   nx = nr+1
   ny = 2*nn-1
   usefft = True
else:
   usefft = False
   
epx = np.zeros([nx,nr],dtype=np.complex)
eny = np.zeros([ny,nn],dtype=np.complex)
x = np.zeros([nx])
y = np.zeros([ny])

#------------------------------------------------------------------------
# Some setup 
#
if usefft:
   for i in range(nx):
      x[i] = i*2*np.pi/nx
   for j in range(ny):
      y[j] = j*2*np.pi/ny
elif hasgapy:
   for i in range(nx):
      x[i] = i*2*np.pi/(nx-1)
   for j in range(ny):
      y[j] = j*2*np.pi/(ny-1)
else:
   # Fourier arrays
   for i in range(nx):
      x[i] = i*2*np.pi/(nx-1)
      for p in range(nr):    
         epx[i,p]=np.exp(1j*(p-nr/2)*x[i])

   for j in range(ny):
      y[j] = j*2*np.pi/(ny-1)
      for n in range(nn):    
         eny[j,n]=np.exp(-1j*n*y[j])

   # factor of 1/2 for n=0
   eny[:,0] = 0.5*eny[:,0] 
      
#------------------------------------------------------------------------
# Real-space field reconstruction (if no gapy)
def maptoreal(nr,nn,nx,ny,c):

    import numpy as np
    import time

    start = time.time()

    # This needs to be fast, so we use numpy.outer
    f = np.zeros([nx,ny])
    for p in range(nr):
        for n in range(nn):
            f[:,:] = f[:,:]+np.real(c[p,n]*np.outer(epx[:,p],eny[:,n]))
    
    end = time.time()
  
    return f,end-start
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# FFT version
def maptoreal_fft(nr,nn,nx,ny,c):

   import numpy as np
   import time
   
   d = np.zeros([nx,ny],dtype=np.complex)

   # Mapping
   # d[ ix, iy] = c[ ix,iy] 
   # d[ ix,-iy] = c[-ix,iy]^* 

   for ix in range(-nr/2+1,nr/2):
      i = ix
      if ix < 0:
         i = ix+nx
      d[i,0:nn] = c[-ix+nr/2,0:nn]

   for ix in range(-nr/2,nr/2-1):
      i = ix
      if ix < 0:
         i = ix+nx
      for iy in range(1,nn):
         d[i,ny-iy] = np.conj(c[ix+nr/2,iy])
          
   start = time.time()

   # Sign convention negative exponent exp(-inx)
   f = np.real(np.fft.fft2(d))*0.5
    
   end = time.time()
  
   return f,end-start
#------------------------------------------------------------------------


# Generate vector of time frames 
if istr == '-1':
    ivec = range(nt)
else:
    ivec = str2list(istr)

# Get filename and tags 
fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],moment)

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

# This is the logic to generate a frame
def frame():

   if i in ivec:
      if isfield:
         a = np.reshape(aa,(2,nr,nn),order='F')
         c = a[0,:,:]+1j*a[1,:,:]
      else:
         a = np.reshape(aa,(2,nr,ns,nn),order='F')
         c = a[0,:,species,:]+1j*a[1,:,species,:]
                
      f = np.zeros([nx,ny],order='F')
      if hasgapy:
         gapy.realfluct(c,f)
         t = 0.0
      else:
         if usefft:
            f,t = maptoreal_fft(nr,nn,nx,ny,c)
         else:
            f,t = maptoreal(nr,nn,nx,ny,c)
         
      if fmin == 'auto':
         f0=np.min(f)
         f1=np.max(f)
      else:
         f0=float(fmin)
         f1=float(fmax)

      xp = x/(2*np.pi)*sim.length
      # ky[1] < 0 is possible
      yp = y/np.abs(sim.ky[1])
      aspect = max(abs(yp))/max(abs(xp))

      fig = plt.figure(figsize=(8,8*aspect))
      ax = fig.add_subplot(111)
      ax.set_title(title)
      ax.set_xlabel(r'$x/\rho_s$')
      ax.set_ylabel(r'$y/\rho_s$')
      ax.set_aspect('equal')
        
      levels = np.arange(f0,f1,(f1-f0)/256)
      ax.contourf(xp,yp,np.transpose(f),levels,cmap=plt.get_cmap(colormap))
      print 'INFO: (plot_fluct) min=%e , max=%e  (t=%e)' % (f0,f1,t)

      fig.tight_layout(pad=0.3)
      plt.subplots_adjust(top=0.94)
      if ftype == 'screen':
         plt.show()
      else:
         fname = fdata+str(i)
         # Filename uses frame number 
         plt.savefig(str(i)+'.'+ftype)
         # Close each time to prevent memory accumulation
         plt.close()
            
      if i == max(ivec):
       sys.exit()

    
i = 0

if hasbin:

   # Binary data
   work = True
   # Open binary file
   fbin = open(fdata,'rb') 

   while work:
      try:
         aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
      except:
         sys.exit()
      
      i = i+1
      print 'INFO: (plot_fluct) Time index '+str(i) 
      frame()

else:

   m = 0

   aa = np.zeros([n_chunk])
   for line in open(fdata):
      aa[m] = float(line)
      m = m+1
      if m == n_chunk:
         i = i+1
         m = 0
         print 'INFO: (plot_fluct) Time index '+str(i) 
         frame()

