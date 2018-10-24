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
   hasgapy = True
except:
   print 'ERROR: (vis_torcut) Please build gapy.so library!'
   sys.exit()
   
PREC='f' ; BIT=4

ftype    = sys.argv[1]
moment   = sys.argv[2]
species  = int(sys.argv[3])
nx       = int(sys.argv[4])
ny       = int(sys.argv[5])
istr     = sys.argv[6]
fmin     = sys.argv[7]
fmax     = sys.argv[8]
colormap = sys.argv[9]
font     = int(sys.argv[10])

# Define plot and font size 
rc('text',usetex=True)
rc('font',size=font)

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species
nth = sim.theta_plot

# Generate vector of time frames 
if istr == '-1':
   ivec = [nt]
elif istr == 'all':
   ivec = range(nt)
else:
   ivec = str2list(istr)

#------------------------------------------------------------------------
# (r,theta)=(x,y) mesh setup 
#
if nth == 1:
   print 'WARNING: (vis_torcut) Should use THETA_PLOT > 1 in CGYRO.'

if nx < 0 or ny < 0:
   nx = nr+1
   ny = 2*nn-1

x = np.zeros([nx])
y = np.zeros([ny])

for i in range(nx):
   x[i] = i*2*np.pi/(nx-1)
for j in range(ny):
   y[j] = j*2*np.pi/(ny-1)-np.pi

xp = np.zeros([nx,ny])
yp = np.zeros([nx,ny])
zp = np.zeros([nx,ny])

for i in range(nx):
   for j in range(ny):
      r = 0.1+x[i]/(4*np.pi)
      xp[i,j] = 1.0+r*np.cos(y[j])
      yp[i,j] = sim.kappa*r*np.sin(y[j])
      zp[i,j] = 0.0

# Shape functions
gapy.geo.geo_rmin_in=sim.rmin
gapy.geo.geo_rmaj_in=sim.rmaj
gapy.geo.geo_q_in=sim.q
gapy.geo.geo_s_in=sim.shear
gapy.geo.geo_kappa_in=sim.kappa
gapy.geo.geo_delta_in=sim.delta
gapy.geo.geo_drmaj_in=sim.shift

gapy.geo.geo_interp(y,True)
# g1 -> q*theta
g1 = -gapy.geo.geo_nu
# g2 -> theta 
g2 = gapy.geo.geo_b*gapy.geo.geo_captheta/gapy.geo.geo_s_in/gapy.geo.geo_grad_r**2

print g1
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

   mlab.figure(size=(900,900))
   if isfield:
      c = a[0,:,:,:]+1j*a[1,:,:,:]
   else:
      c = a[0,:,:,species,:]+1j*a[1,:,:,species,:]
                
   f = np.zeros([nx,ny],order='F')
   gapy.torcut(sim.m_box,sim.q,g1,g2,c,f)

   if fmin == 'auto':
      f0=np.min(f)
      f1=np.max(f)
   else:
      f0=float(fmin)
      f1=float(fmax)

   mlab.mesh(xp,yp,zp,scalars=f,colormap=colormap,vmin=f0,vmax=f1,opacity=1.0)
   # View from positive z-axis
   mlab.view(azimuth=0, elevation=0)
   print 'INFO: (vis_torcut) min=%e , max=%e' % (f0,f1)

   if ftype == 'screen':
      mlab.show()
   else:
      fname = fdata+str(i)
      # Filename uses frame number 
      mlab.savefig(str(i)+'.'+ftype)
      # Close each time to prevent memory accumulation
      mlab.close()
                
i = 0

work = True
# Open binary file
fbin = open(fdata,'rb') 

while work:
   try:
      aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
   except:
      sys.exit()
      
   i = i+1
   print 'INFO: (vis_torcut) Time index '+str(i) 
   if i in ivec:
      frame()
   if i == max(ivec):
      sys.exit()
