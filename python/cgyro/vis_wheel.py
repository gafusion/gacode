import struct
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rc
from gacodefuncs import *
from cgyro.data import cgyrodata
from mayavi import mlab 
try:
   import gapy
except:
   print 'ERROR: (vis_torcut) Please build gapy.so library!'
   sys.exit()
   
ext      = sys.argv[1]
moment   = sys.argv[2]
species  = int(sys.argv[3])
nx       = int(sys.argv[4])
ny       = int(sys.argv[5])
nz       = int(sys.argv[6])
istr     = sys.argv[7]
fmin     = sys.argv[8]
fmax     = sys.argv[9]
colormap = sys.argv[10]
font     = int(sys.argv[11])

# Define plot and font size 
rc('text',usetex=True)
rc('font',size=font)

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species
nth = sim.theta_plot

ivec = time_vector(istr,nt)

xp = sim.length
yp = 2*np.pi/np.abs(sim.ky[1])
aspect = yp/xp

pre,ftype = mkfile(ext)

#------------------------------------------------------------------------
# (r,theta)=(x,y) mesh setup 
#
if nth == 1:
   print 'ERROR: (vis_wheel) This visualization requires THETA_PLOT > 1.'
   sys.exit()
   
if nx < 0:
   nx = 128
if ny < 0:
   ny = 128
if nz < 0:
   nz = 128

x = np.zeros([nx]) # x
y = np.zeros([ny]) # alpha
z = np.zeros([nz]) # theta

for i in range(nx):
   x[i] = i*2*np.pi/(nx-1)
for j in range(ny):
   y[j] = j*2*np.pi/(ny-1)
for k in range(nz):
   z[k] = k*np.pi/(nz-1)-np.pi

# 1. 
xp1 = np.zeros([nx,ny])
yp1 = np.zeros([nx,ny])
zp1 = np.zeros([nx,ny])

for i in range(nx):
   for j in range(ny):
      dx = x[i]/(4*np.pi)
      xp1[i,j] = 0.2+dx
      yp1[i,j] = 0.0
      zp1[i,j] = y[j]/(4*np.pi)*aspect

# 2. 
xp2 = np.zeros([nx,nz])
yp2 = np.zeros([nx,nz])
zp2 = np.zeros([nx,nz])

for i in range(nx):
   for k in range(nz):
      dx = x[i]/(4*np.pi)
      xp2[i,k] = (0.2+dx)*np.cos(z[k])
      yp2[i,k] = (0.2+dx)*np.sin(z[k])
      zp2[i,k] = 0.5*aspect

# 3. 
xp3 = np.zeros([ny,nz])
yp3 = np.zeros([ny,nz])
zp3 = np.zeros([ny,nz])

for j in range(ny):
   for k in range(nz):
      dx = x[0]/(4*np.pi)
      xp3[j,k] = (0.2+dx)*np.cos(z[k])
      yp3[j,k] = (0.2+dx)*np.sin(z[k])
      zp3[j,k] = y[j]/(4*np.pi)*aspect
# 4. 
xp4 = np.zeros([ny,nz])
yp4 = np.zeros([ny,nz])
zp4 = np.zeros([ny,nz])

for j in range(ny):
   for k in range(nz):
      dx = x[-1]/(4*np.pi)
      xp4[j,k] = (0.2+dx)*np.cos(z[k])
      yp4[j,k] = (0.2+dx)*np.sin(z[k])
      zp4[j,k] = y[j]/(4*np.pi)*aspect

#------------------------------------------------------------------------

# Get filename and tags 
fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],moment)

# Check to see if data exists 
if os.path.isfile('bin'+fdata):
    fdata = 'bin'+fdata
    print 'INFO: (vis_torcut) Found binary data in '+fdata 
else:
    print 'ERROR: (vis_torcut) No data for -moment '+moment+' exists.  Try -moment phi'
    sys.exit()

if isfield:
    n_chunk = 2*nr*nn*nth
else:
    n_chunk = 2*nr*ns*nn*nth
   
# This is the logic to generate a frame
def frame():
   
   if isfield:
      a = np.reshape(aa,(2,nr,nth,nn),order='F')
   else:
      a = np.reshape(aa,(2,nr,nth,ns,nn),order='F')

   mlab.figure(size=(1000,800))
   if isfield:
      c = a[0,:,:,:]+1j*a[1,:,:,:]
   else:
      c = a[0,:,:,species,:]+1j*a[1,:,:,species,:]

   
   # 1a
   f = np.zeros([nx,ny],order='F')
   gapy.realfluct(c[:,nth/2,:],f)

   if fmin == 'auto':
      f0=np.min(f)
      f1=np.max(f)
   else:
      f0=float(fmin)
      f1=float(fmax)

   mlab.mesh(xp1,yp1,zp1,scalars=f,colormap=colormap,vmin=f0,vmax=f1,opacity=1.0)

   # 1b
   f = np.zeros([nx,ny],order='F')
   gapy.realfluct(c[:,0,:],f)
   mlab.mesh(-xp1,yp1,zp1,scalars=f,colormap=colormap,vmin=f0,vmax=f1,opacity=1.0)

   # 2a,b
   f = np.zeros([nx,nz],order='F')
   gapy.wheel1(c,f)
   mlab.mesh(xp2,yp2,zp2,scalars=f,colormap=colormap,vmin=f0,vmax=f1,opacity=1.0)
   mlab.mesh(xp2,yp2,0*zp2,scalars=f,colormap=colormap,vmin=f0,vmax=f1,opacity=1.0)

   # 3a,b
   f = np.zeros([ny,nz],order='F')
   gapy.wheel2(c,f)
   mlab.mesh(xp3,yp3,zp3,scalars=f,colormap=colormap,vmin=f0,vmax=f1,opacity=1.0)
   mlab.mesh(xp4,yp4,zp4,scalars=f,colormap=colormap,vmin=f0,vmax=f1,opacity=1.0)

   # View from positive z-axis
   mlab.view(azimuth=75, elevation=65,distance=3.2)
   print 'INFO: (vis_wheel) min=%e , max=%e' % (f0,f1)

   if ftype == 'screen':
      mlab.show()
   else:
      # Filename uses frame number 
      mlab.savefig(pre+str(i)+'.'+ftype)
      # Close each time to prevent memory accumulation
      mlab.close()
                

PREC='f' ; BIT=4

# Open binary file
fbin = open(fdata,'rb') 

i = 0
while True:
   try:
      aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
   except:
      sys.exit()

   i += 1
   print 'INFO: (vis_torcut) Time index '+str(i) 
   if i in ivec:
      frame()
   if i == max(ivec):
      sys.exit()
