import struct
import sys
import numpy as np
import os
import time
from matplotlib import rc
import matplotlib.pyplot as plt
from ..gacodefuncs import *
from .data import cgyrodata

try:
   from pygacode import vis
   haspygacode = True
except:
   print("BAD: (plot_fluct) Please 'pip install pygacode'")
   haspygacode = False

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
theta = int(sys.argv[14])

# Use first 3 args to define plot and font size
rc('text',usetex=False)
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
itheta,thetapi = indx_theta(theta,nth)

if nx < 0 or ny < 0:
   nx = nr+1
   ny = 2*nn-1
   usefft = True
else:
   usefft = False

epx = np.zeros([nx,nr],dtype=complex)
eny = np.zeros([ny,nn],dtype=complex)

# FFT-consistent real-space meshpoints
x = np.arange(nx)*2*np.pi/nx
y = np.arange(ny)*2*np.pi/ny

#------------------------------------------------------------------------
# Some setup
#

if usefft:
   mode = 'FFT'
elif haspygacode:
   mode = 'pygacode'
else:
   mode = 'slow'
   # Fourier arrays
   for i in range(nx):
      for p in range(nr):
         epx[i,p]=np.exp(1j*(p-nr//2)*x[i])

   for j in range(ny):
      for n in range(nn):
         eny[j,n]=np.exp(-1j*n*y[j])

   # Only computing half sum 
   eny[:,0] = 0.5*eny[:,0]

print('INFO: (plot_fluct) Computation method: '+mode)
   
#------------------------------------------------------------------------
# Real-space field reconstruction (if no pygacode)
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

   # Storage for numpy inverse real transform (irfft2)
   d = np.zeros([nx,nn],dtype=complex)

   start = time.time()

   for i in range(nr):
      p = i-nr//2
      # k is the "standard FFT index"
      if -p < 0:
         k = -p+nx
      else:
         k = -p
      # Use identity f(p,-n) = f(-p,n)* 
      d[k,0:nn] = np.conj(c[i,0:nn])

   # 2D inverse real Hermitian transform
   # NOTE: using inverse FFT with convention exp(ipx+iny), so need n -> -n 
   # NOTE: need factor of 0.5 to match half-sum method of slow maptoreal()
   f = np.fft.irfft2(d,s=[nx,ny],norm='forward')*0.5
   
   end = time.time()

   return f,end-start
#------------------------------------------------------------------------

# Get filename and tags
print('INFO: (plot_fluct) '+moment+' fluctuations')
if moment == 't':
	fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],'n')
	fdata_n = 'bin' + fdata
	fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],'e')	
	fdata_e = 'bin' + fdata
	u = specmap(sim.mass[species],sim.z[species])
	title = r'${\delta \mathrm{T}}_'+u+'$'
else:
	fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],moment)

# Check to see if data exists (try binary data first)
if os.path.isfile('bin'+fdata):
    fdata = 'bin'+fdata
    print('INFO: (plot_fluct) Found binary data in '+fdata)
    hasbin = True
elif os.path.isfile('out'+fdata):
    fdata = 'out'+fdata
    print('BAD: (plot_fluct) Using inefficient ASCII data in '+fdata)
    hasbin = False
else:
    print('ERROR: (plot_fluct) No data for -moment '+moment+' exists.  Try -moment phi')
    sys.exit()

if isfield:
    n_chunk = 2*nr*nth*nn
else:
    n_chunk = 2*nr*nth*ns*nn

# This is the logic to generate a frame
def frame():

   if i in ivec:
      if isfield:
         a = np.reshape(aa,(2,nr,nth,nn),order='F')
         c = a[0,:,itheta,:]+1j*a[1,:,itheta,:]
      elif moment == 't':
         a = np.reshape(aa,(2,nr,nth,ns,nn),order='F')
         b = np.reshape(bb,(2,nr,nth,ns,nn),order='F')
         a = ((2./3)*b - sim.temp[species]*a)/sim.dens[species] #dE = (3/2)*(T*delta n + n*deltaT)
         c = a[0,:,itheta,species,:]+1j*a[1,:,itheta,species,:]
      else:
         a = np.reshape(aa,(2,nr,nth,ns,nn),order='F')
         c = a[0,:,itheta,species,:]+1j*a[1,:,itheta,species,:]

      f = np.zeros([nx,ny],order='F')
      if usefft:
         f,t = maptoreal_fft(nr,nn,nx,ny,c)
      elif haspygacode:
         start = time.time()
         # see f2py/vis/vis.f90
         vis.realfluct(c,f)
         end = time.time()
         t = end-start
      else:
         f,t = maptoreal(nr,nn,nx,ny,c)

      # Correct for half-sum
      f = 2*f
            
      if fmin == 'auto':
         f0=np.min(f)
         f1=np.max(f)

         #CH alternate suggestion ensure even range around 0 for fluctuations
         flim = max(abs(f0),abs(f1))
         f0=-flim
         f1=flim

      else:
         f0=float(fmin)
         f1=float(fmax)

      # Physical maxima
      xmax = sim.length
      ymax = (2*np.pi)/np.abs(sim.ky[1])
      xp = x/(2*np.pi)*xmax
      yp = y/(2*np.pi)*ymax

      # Periodic extensions
      xp = np.append(xp,xmax)
      yp = np.append(yp,ymax)
      fp = np.zeros([nx+1,ny+1])
      fp[0:nx,0:ny] = f[:,:]
      fp[-1,:] = fp[0,:]
      fp[:,-1] = fp[:,0]
      
      levels = np.arange(f0,f1,(f1-f0)/256)
      if land == 0:
         fig = plt.figure(figsize=(px/100.0,py/100.0))
         ax = fig.add_subplot(111)
         ax.set_xlabel(r'$x/\rho_s$')
         ax.set_ylabel(r'$y/\rho_s$')
         ax.set_xlim([0,xmax])
         ax.set_ylim([0,ymax])
         ax.contourf(xp,yp,np.transpose(fp),levels,cmap=plt.get_cmap(colormap))
      else:
         fig = plt.figure(figsize=(px/100.0,py/100.0))
         ax = fig.add_subplot(111)
         ax.set_xlabel(r'$y/\rho_s$')
         ax.set_ylabel(r'$x/\rho_s$')
         ax.set_xlim([0,ymax])
         ax.set_ylim([0,xmax])
         ax.contourf(yp,xp,fp,levels,cmap=plt.get_cmap(colormap))

      print('INFO: (plot_fluct '+mode+') min=%e , max=%e  (t=%e)' % (f0,f1,t))
      print('INFO: (plot_fluct) Shape = '+str(f.shape))

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
            
      if i == max(ivec):
       sys.exit()

    
i = 0

if moment == 't':
   # Binary data
   work = True
   # Open binary file
   with open(fdata_n,'rb') as fbin_n, open(fdata_e,'rb') as fbin_e:
       while work:
          try:
             aa = struct.unpack(PREC*n_chunk,fbin_n.read(BIT*n_chunk))
             bb = struct.unpack(PREC*n_chunk,fbin_e.read(BIT*n_chunk))
          except:
             sys.exit()

          i = i+1
          print('INFO: (plot_fluct) Time index '+str(i)) 
          frame()

elif hasbin:
   # Binary data
   work = True
   # Open binary file
   with open(fdata,'rb') as fbin:
       while work:
          try:
             aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
          except:
             sys.exit()

          i = i+1
          print('INFO: (plot_fluct) Time index '+str(i)) 
          frame()

else:

   m = 0

   aa = np.zeros([n_chunk])
   with open(fdata) as f:
      for line in f:
         aa[m] = float(line)
         m = m+1
         if m == n_chunk:
            i = i+1
            m = 0
            print('INFO: (plot_fluct) Time index '+str(i)) 
            frame()

