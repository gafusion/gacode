import struct
import sys
import numpy as np
import os
from matplotlib import rc
import matplotlib.pyplot as plt
from gacodefuncs import *
from cgyro.data import cgyrodata
from mpi4py import MPI

iproc = MPI.COMM_WORLD.Get_rank()
nproc = MPI.COMM_WORLD.Get_size()

print nproc

PREC='f' ; BIT=4  

ext = sys.argv[1]
moment = sys.argv[2]
species = int(sys.argv[3])
px = int(sys.argv[4])
py = int(sys.argv[5])
nx = int(sys.argv[6])
ny = int(sys.argv[7])
istr = sys.argv[8]
fmin = sys.argv[9]
fmax = sys.argv[10]
colormap = sys.argv[11]
font = int(sys.argv[12])
land = int(sys.argv[13])
theta = float(sys.argv[14])

# Use first 3 args to define plot and font size 
rc('text',usetex=True)
rc('font',size=font)

# Extension handling
pre,ftype=mkfile(ext)
   
sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species
nth = sim.theta_plot

ivec = time_vector(istr,nt)
itheta = theta_indx(theta,nth)

# FFT only
nx = nr+1
ny = 2*nn-1
   
epx = np.zeros([nx,nr],dtype=np.complex)
eny = np.zeros([ny,nn],dtype=np.complex)
x = np.zeros([nx])
y = np.zeros([ny])

#------------------------------------------------------------------------
# Some FFT setup (only option)
#
mode = 'FFT'
for i in range(nx):
   x[i] = i*2*np.pi/nx
for j in range(ny):
   y[j] = j*2*np.pi/ny
 
def maptoreal_fft(nr,nn,nx,ny,c):

   import numpy as np
   import time
   
   d = np.zeros([nx,ny],dtype=np.complex)

   # Mapping
   # d[ ix, iy] = c[ ix,iy] 
   # d[ ix,-iy] = c[-ix,iy]^* 

   start = time.time()

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
          

   # Sign convention negative exponent exp(-inx)
   f = np.real(np.fft.fft2(d))*0.5
    
   end = time.time()
  
   return f,end-start
#------------------------------------------------------------------------

# Get filename and tags 
fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],moment)

# Check to see if data exists (try binary data first)
if os.path.isfile('bin'+fdata):
    fdata = 'bin'+fdata
    print 'INFO: (plot_fluct_mpi) Found binary data in '+fdata 
    hasbin = True
else:
    print 'ERROR: (plot_fluct_mpi) No data for -moment '+moment+' exists.  Try -moment phi'
    sys.exit()

if isfield:
    n_chunk = 2*nr*nth*nn
else:
    n_chunk = 2*nr*nth*ns*nn

# This is the logic to generate a frame
def frame():

   if isfield:
      a = np.reshape(aa,(2,nr,nth,nn),order='F')
      c = a[0,:,itheta,:]+1j*a[1,:,itheta,:]
   else:
      a = np.reshape(aa,(2,nr,nth,ns,nn),order='F')
      c = a[0,:,itheta,species,:]+1j*a[1,:,itheta,species,:]
                
   f = np.zeros([nx,ny],order='F')
      
   f,t = maptoreal_fft(nr,nn,nx,ny,c)
      
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
      
   levels = np.arange(f0,f1,(f1-f0)/256)
   if land == 0:
      fig = plt.figure(figsize=(px/100.0,py/100.0))
      ax = fig.add_subplot(111)
      ax.set_xlabel(r'$x/\rho_s$')
      ax.set_ylabel(r'$y/\rho_s$')
      ax.contourf(xp,yp,np.transpose(f),levels,cmap=plt.get_cmap(colormap),extend='both')
      plt.subplots_adjust(top=0.94)
   else:
      fig = plt.figure(figsize=(px/100.0,py/100.0))
      ax = fig.add_subplot(111)
      ax.set_xlabel(r'$y/\rho_s$')
      ax.set_ylabel(r'$x/\rho_s$')
      ax.contourf(yp,xp,f,levels,cmap=plt.get_cmap(colormap),extend='both')
 
   print 'INFO: (plot_fluct_mpi '+mode+') min=%2.4f , max=%2.4f  (t=%3.3f)' % (f0,f1,t)

   ax.set_title(title)
   ax.set_aspect('equal')
   fig.tight_layout()

   if ftype == 'screen':
      plt.show()
   else:
      # Filename uses frame number 
      plt.savefig(pre+str(i)+'.'+ftype)
      # Close each time to prevent memory accumulation
      plt.close()
            
i = 0

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
   if np.mod(i,nproc) == iproc:
      print 'INFO: (plot_fluct_mpi) Time index '+str(i),'[',iproc,']' 
      if i in ivec:
         frame() 
   if i == max(ivec):
      sys.exit()
