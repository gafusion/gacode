import struct
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import cm
from ..gacodefuncs import *
from .data import cgyrodata
try:
    from mayavi import mlab
except ModuleNotFoundError:
   print("WARNING: (vis_torcut) Mayavi module not found - will use matplotlib")
   no_mayavi = True 
   import matplotlib.tri as tri
   from matplotlib.colors import Normalize
try:
   from pygacode import vis
   from pygacode import geo
   print("INFO: (vis_torcut) Successfully imported vis")
except:
   print("ERROR: (vis_torcut) Please type 'pip install pygacode'")
   sys.exit()

ext      = sys.argv[1]
moment   = sys.argv[2]
species  = int(sys.argv[3])
nx       = int(sys.argv[4])
nz       = int(sys.argv[5])
istr     = sys.argv[6]
fmin     = sys.argv[7]
fmax     = sys.argv[8]
colormap = sys.argv[9]
font     = int(sys.argv[10])
legacy   = bool(int(sys.argv[11]))
dn       = int(sys.argv[12])
mag      = float(sys.argv[13])
nozonal  = bool(int(sys.argv[14]))
onlyzonal = bool(int(sys.argv[15]))

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species
nth = sim.theta_plot

print(f'INFO: (vis_torcut) dn = {dn} (cf. dn in out.cgyro.info')
print('HINT: adjust -dn to match experimental dn (rho/a and Lx/a will shrink)')
print('Lx/rho = {:.2f}'.format(sim.length))
print('rho/a  = {:.4f}'.format(sim.rho/dn))
lovera = sim.length*sim.rho/dn*mag
print('Lx/a   = {:.4f}'.format(lovera))

# To calculate dn we could use the following formula given rhosunit,
# However, rhosunit is not an attribute of sim.
#dn = sim.ky[1]/sim.q * sim.rmin/sim.rhosunit

ivec = time_vector(istr,nt)
print(f"INFO: (vis_torcut) will produce {len(ivec)} frames")

s=ext.split('.')
if len(s) == 2:
   pre   = s[0]
   ftype = s[1]
else:
   pre = ''
   ftype = s[0]

#------------------------------------------------------------------------
# (r,theta)=(x,y) mesh setup 
#
if nth == 1:
   print('WARNING: (vis_torcut) Should use THETA_PLOT > 1 in CGYRO.')

if nx < 0:
   nx = 128
if nz < 0:
   nz = 128

x = np.zeros([nx])
z = np.zeros([nz])

for i in range(nx):
   x[i] = i*2*np.pi/(nx-1.0)/dn
for k in range(nz):
   z[k] = k*2*np.pi/(nz-1.0)-np.pi
   
xp = np.zeros([nx,nz])
yp = np.zeros([nx,nz])
zp = np.zeros([nx,nz])

for i in range(nx):
   r = sim.rmin+(dn*x[i]/(2*np.pi)-0.5)*lovera
   for k in range(nz):
      a = z[k] + (sim.shape_cos[0]
         + sim.shape_cos[1]*np.cos(z[k]) 
         + sim.shape_cos[2]*np.cos(2*z[k]) 
         + sim.shape_cos[3]*np.cos(3*z[k])
         + np.arcsin(sim.delta)*np.sin(z[k]) 
         - sim.zeta*np.sin(2*z[k]) 
         + sim.shape_sin[3]*np.sin(3*z[k]))
      xp[i,k] = sim.rmaj+r*np.cos(a)
      yp[i,k] = sim.zmag+sim.kappa*r*np.sin(z[k])
      zp[i,k] = 0.0

# Shape functions 
geo.signb_in=1 # fix
geo.geo_rmin_in=sim.rmin
geo.geo_rmaj_in=sim.rmaj
geo.geo_drmaj_in=sim.shift
geo.geo_zmag_in=sim.zmag
geo.geo_dzmag_in=sim.dzmag
geo.geo_q_in=sim.q
geo.geo_s_in=sim.shear
geo.geo_kappa_in=sim.kappa
geo.geo_delta_in=sim.delta
geo.geo_zeta_in=sim.zeta
geo.geo_s_kappa_in=sim.s_kappa
geo.geo_s_delta_in=sim.s_delta
geo.geo_s_zeta_in=sim.s_zeta
geo.geo_beta_star_in=sim.beta_star
geo.geo_shape_cos0_in = sim.shape_cos[0]
geo.geo_shape_cos1_in = sim.shape_cos[1]
geo.geo_shape_cos2_in = sim.shape_cos[2]
geo.geo_shape_cos3_in = sim.shape_cos[3]
geo.geo_shape_cos4_in = 0.0
geo.geo_shape_cos5_in = 0.0
geo.geo_shape_cos6_in = 0.0
geo.geo_shape_sin3_in = sim.shape_sin[3]
geo.geo_shape_sin4_in = 0.0
geo.geo_shape_sin5_in = 0.0
geo.geo_shape_sin6_in = 0.0
geo.geo_shape_s_cos0_in = sim.shape_s_cos[0]
geo.geo_shape_s_cos1_in = sim.shape_s_cos[1]
geo.geo_shape_s_cos2_in = sim.shape_s_cos[2]
geo.geo_shape_s_cos3_in = sim.shape_s_cos[3]
geo.geo_shape_s_cos4_in = 0.0
geo.geo_shape_s_cos5_in = 0.0
geo.geo_shape_s_cos6_in = 0.0
geo.geo_shape_s_sin3_in = sim.shape_s_sin[3]
geo.geo_shape_s_sin4_in = 0.0
geo.geo_shape_s_sin5_in = 0.0
geo.geo_shape_s_sin6_in = 0.0

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
   
if int(mag) == 0:
   showco=True
else:
   showco=False
if showco:
   # (Optional) plot of the geometry functions, then exit
   rc('text',usetex=True) ; rc('font',size=font)
   fig = plt.figure(figsize=(10,8))

   ax = fig.add_subplot(111)
   ax.plot(z/np.pi,g1/sim.q,label='g1')
   ax.plot(z/np.pi,g2,label='g2')
   ax.plot(z/np.pi,z,'--k')
   ax.set_xlim([-1,1])
   ax.legend()
   plt.show()
   sys.exit()
   
#------------------------------------------------------------------------

# Get filename and tags 
fdata,title,isfield = tag_helper(sim.mass[species],sim.z[species],moment)

# Check to see if data exists 
if os.path.isfile('bin'+fdata):
    fdata = 'bin'+fdata
    print('INFO: (vis_torcut) Found binary data in '+fdata) 
else:
    print('ERROR: (vis_torcut) No data for -moment '+moment+' exists.  Try -moment phi')
    sys.exit()

if isfield:
    n_chunk = 2*nr*nn*nth
else:
    n_chunk = 2*nr*ns*nn*nth

# If we use matplotlib instead of mayavi, triangulate outside of frame() routine for efficiency
if no_mayavi:
    x = xp.flatten()
    y = yp.flatten()
    triang = tri.Triangulation(x,y)
    triangles = triang.triangles
    # Remove large triangles that reach across the flux tube.
    dR = np.max(np.diff(xp[:, nz//2])) # theta = 0 mesh spacing in dR
    dZ = np.max(np.diff(yp[-1,:])) # outer flux surface dZ around theta.
    max_length = 2*max([dR, dZ])
    print(f"INFO: (vis_torcut) triangulation max side length [/a] = {max_length:.3f}")
    # Mask off unwanted triangles.
    xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
    ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
    maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
    mask = maxi < max_length
    triang.set_mask(~mask)
    print(f"INFO: (vis_torcut) mesh contains {len(triangles)} triangles")
   
# This is the logic to generate a frame
def frame():
   
   if isfield:
      a = np.reshape(aa,(2,nr,nth,nn),order='F')
   else:
      a = np.reshape(aa,(2,nr,nth,ns,nn),order='F')
   if no_mayavi:
       fig, ax = plt.subplots(1,1,figsize=(8,8),clear=True)
       ax.plot(xp[0,:],yp[0,:],'k-')
       ax.plot(xp[-1,:],yp[-1,:],'k-')
       ax.set_aspect("equal")
       ax.set_ylabel("Z/a")
       ax.set_xlabel("R/a")
   else:
       mlab.figure(size=(900,900),bgcolor=(1,1,1))
       
   if isfield:
      c = a[0,:,:,:]+1j*a[1,:,:,:]
   else:
      c = a[0,:,:,species,:]+1j*a[1,:,:,species,:]
                
   if nozonal and nn > 1:
      c[:,:,0] = 0.0
   if onlyzonal and nn > 1:
      c[:,:,1:] = 0.0

   f = np.zeros([nx,nz],order='F')
   vis.torcut(dn,sim.m_box,sim.q,sim.thetap,g1,g2,c,f)

   if fmin == 'auto':
      f0=np.min(f)
      f1=np.max(f)
   else:
      f0=float(fmin)
      f1=float(fmax)
      
   if no_mayavi:
       norm = Normalize(f0, f1)
       cmap = plt.get_cmap(colormap)
       ax.tripcolor(triang, f.flatten(), shading='gouraud',norm=norm,cmap=cmap)
   else:
       image = mlab.mesh(xp,yp,zp,scalars=f,colormap=colormap,vmin=f0,vmax=f1,opacity=1.0)
       # View from positive z-axis
       mlab.view(azimuth=0, elevation=0)
   
   print('INFO: (vis_torcut) min={:.3f} | max={:.3f}'.format(f0,f1))
   
   if ftype == 'screen':
      if no_mayavi:
          plt.show()
      else:
          mlab.show()
   else:
       if no_mayavi:
           plt.savefig(pre+str(i)+'.'+ftype)
           plt.close("all")
       else:
          # Filename uses frame number 
          mlab.savefig(pre+str(i)+'.'+ftype)
          # Close each time to prevent memory accumulation
          mlab.close()
                
PREC='f' ; BIT=4

# Open binary file
with open(fdata,'rb') as fbin:
   i = 0
   while True:
      try:
         aa = struct.unpack(PREC*n_chunk,fbin.read(BIT*n_chunk))
      except:
         sys.exit()

      i += 1
      print('INFO: (vis_torcut) Time index {:d} '.format(i))
      if i in ivec:
         frame()
      if i == max(ivec):
         sys.exit()
