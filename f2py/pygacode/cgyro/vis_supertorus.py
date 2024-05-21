import struct
import sys
import time
import numpy as np
import os
from matplotlib import cm
from ..gacodefuncs import *
from .data import cgyrodata
from .vis_mesh import *

ext      = sys.argv[1]
moment   = sys.argv[2]
species  = int(sys.argv[3])
nx       = int(sys.argv[4])
nz       = int(sys.argv[5])
nphi     = int(sys.argv[6])
phimax   = float(sys.argv[7])
istr     = sys.argv[8]
fmin_in  = sys.argv[9]
fmax_in  = sys.argv[10]
colormap = sys.argv[11]
font     = int(sys.argv[12])
legacy   = bool(int(sys.argv[13]))
dn       = int(sys.argv[14])
mag      = float(sys.argv[15])
nozonal  = bool(int(sys.argv[16]))
onlyzonal = bool(int(sys.argv[17]))
px        = int(sys.argv[18])
py        = int(sys.argv[19])

try:
   from vis.vis import vis
   print("INFO: (vis_torcut) Successfully imported vis")
except:
   print("ERROR: (vis_torcut) Please 'make vismod' in gacode/f2py")
   sys.exit()
   
try:
   from pygacode import geo
   print("INFO: (vis_torcut) Successfully imported geo")
except:
   print("ERROR: (vis_torcut) Please 'pip install pygacode'")
   sys.exit()

try:
   from mayavi import mlab
   print("INFO: (vis_torcut) Successfully imported mayavi")
except:
   print("ERROR: (vis_torcut) Mayavi module not found.")
   sys.exit()

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species
nth = sim.theta_plot

lovera = sim.length*sim.rho/dn*mag

print('HINT: adjust -dn to match experimental dn (rho/a and Lx/a will shrink)')
print('Lx/rho = {:.2f}'.format(sim.length))
print('rho/a  = {:.4f}'.format(sim.rho/dn))
print('Lx/a   = {:.4f}'.format(lovera))
print('px = {}  py = {}'.format(px,py))
print('nx = {}  nz = {}  nphi = {}'.format(nx,nz,nphi))

ivec = time_vector(istr,nt)

s=ext.split('.')
if len(s) == 2:
   pre   = s[0]
   ftype = s[1]
else:
   pre = ''
   ftype = s[0]

#------------------------------------------------------------------------
# Mesh setup 
#
start = time.time()

if nth == 1:
   print('WARNING: (vis_supertorus) Should use THETA_PLOT > 1 in CGYRO.')

if nx < 0:
   nx = 128
if nz < 0:
   nz = 128

x,z,xp,yp,zp = vis_mesh(sim,nx,nz,dn,lovera)

print('MESH TIME = '+'{:.3e}'.format(time.time()-start)+' s.')

# Define shape functions for GEO evaluation
geo = vis_geo(sim,geo)
geo.geo_interp(z,True)

if legacy:
   # s-alpha approximate (apparently used in legacy GYRO movies)
   # g1 -> q*theta
   # g2 -> theta 
   g1 = sim.q*z 
   g2 = z
else:
   # Correct form of Clebsch angle expansion nu(r,theta)
   g1 = -geo.geo_nu 
   g2 = geo.geo_b*geo.geo_captheta/geo.geo_s_in/geo.geo_grad_r**2
   
#------------------------------------------------------------------------

# Get filename and tags 
fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],moment)

# Check to see if data exists 
if os.path.isfile('bin'+fdata):
    fdata = 'bin'+fdata
    print('INFO: (vis_supertorus) Found binary data in '+fdata) 
else:
    print('ERROR: (vis_supertorus) No data for -moment '+moment+' exists.  Try -moment phi')
    sys.exit()

if isfield:
    n_chunk = 2*nr*nn*nth
else:
    n_chunk = 2*nr*ns*nn*nth
    
def frame():

   global c
   
   if isfield:
      a = np.reshape(aa,(2,nr,nth,nn),order='F')
   else:
      a = np.reshape(aa,(2,nr,nth,ns,nn),order='F')

   mlab.figure(bgcolor=(1,1,1))
   if isfield:
      c = a[0,:,:,:]+1j*a[1,:,:,:]
   else:
      c = a[0,:,:,species,:]+1j*a[1,:,:,species,:]
                
   if nozonal and nn > 1:
      c[:,:,0] = 0.0
   if onlyzonal and nn > 1:
      c[:,:,1:] = 0.0
   
# This is the logic to generate an (r,theta) frame
def subframe_torcut(iframe):

   global c,fmin,fmax
   
   phi = -2*np.pi*phimax*iframe
   xpp = xp*np.cos(phi)
   ypp = yp
   zpp = xp*np.sin(phi)

   f = np.zeros([nx,nz],order='F')
   
   vis.torcut(dn,sim.m_box,sim.q,sim.thetap,g1-phi/dn,g2,c,f)

   if iframe == 0:
      if fmin_in == 'auto':
         fmin = np.min(f)
         fmax = np.max(f)
      else:
         fmin = float(fmin_in)
         fmax = float(fmax_in)
      print('INFO: (vis_supertorus) min={:.3f} | max={:.3f}'.format(fmin,fmax))
         
   image = mlab.mesh(xpp,ypp,zpp,scalars=f,
                     colormap=colormap,vmin=fmin,vmax=fmax,opacity=1.0)

# This is the logic to generate an (theta,phi) frame
def subframe_torside():

   global c,fmin,fmax
   
   phi = -np.linspace(0,2*np.pi*phimax,nphi)
   xpp = np.zeros([nz,nphi])
   ypp = np.zeros([nz,nphi])
   zpp = np.zeros([nz,nphi])
   
   f = np.zeros([2,nz,nphi],order='F')
   
   vis.torside(dn,sim.m_box,sim.q,phimax,sim.thetap,g1,g2,c,f)

   for i in [0,1]:
      xpp[:,:] = xp[-i,:,None]*np.cos(phi[None,:])
      ypp[:,:] = yp[-i,:,None]
      zpp[:,:] = xp[-i,:,None]*np.sin(phi[None,:])

      image = mlab.mesh(xpp[:,:],ypp[:,:],zpp[:,:],scalars=f[i,:,:],
                        colormap=colormap,vmin=fmin,vmax=fmax,opacity=1.0) 


PREC='f' ; BIT=4

print('ivec {}'.format(ivec))

# Open binary file
with open(fdata,'rb') as fbin:
   i = 0
   while True:
      try:
         aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
      except:
         sys.exit()

      i += 1
      print('INFO: (vis_supertorus) Time index {:d} '.format(i))
      if i in ivec:
         frame()
         start = time.time()
         # cap at phi=0
         subframe_torcut(0)
         # cap at phi=-phimax
         subframe_torcut(1)
         # inner and outer body of torus
         subframe_torside()
         # View from positive z-axis
         mlab.view(azimuth=0,elevation=0)
         print('MAP TIME = '+'{:.3e}'.format(time.time()-start)+' s.')
         if ftype == 'screen':
            mlab.show()
         else:
            # Filename uses frame number 
            mlab.savefig(pre+str(i)+'.'+ftype,size=(px,py))
            # Close each time to prevent memory accumulation
            mlab.close()
         
      if i == max(ivec):
         sys.exit()
