import sys
import string
import numpy as np

from ..gacodefuncs import *
from .prgen_shape_util import *

def prgen_shape(r,z,narc,nf,xplot):

   # Number of theta-points for plotting
   dx = np.zeros(narc)
   ur = np.zeros(narc) ; uz = np.zeros(narc)
   vr = np.zeros(narc) ; vz = np.zeros(narc)

   # Extrema ; definitions of rmin, rmaj, etc.
   n1 = np.argmax(z) ; m1 = np.argmin(z)
   n2 = np.argmax(r) ; m2 = np.argmin(r)
   za = z[n1] ; zb = z[m1]
   ra = r[n2] ; rb = r[m2]
   zmaj = 0.5*(za+zb) ; zmin = 0.5*(za-zb)
   rmaj = 0.5*(ra+rb) ; rmin = 0.5*(ra-rb)

   if z[1] < z[0]:
      # Reverse order (may be needed)
      r = np.flip(r,0) ; z = np.flip(z,0)

   # Compute generalized angles
   eps = 1.0-1e-9
   for i in range(narc):
      # (ur,uz): principle angles (discontinuous)
      uz[i] = np.arcsin(eps*(z[i]-zmaj)/zmin)
      ur[i] = np.arccos(eps*(r[i]-rmaj)/rmin)
      # (vr,vz): proper (continuous) branches
      if i > 4:
         if uz[i] > uz[i-1] and ur[i] > ur[i-1]:
            vz[i] = uz[i] ; vr[i] = ur[i]
         elif uz[i] < uz[i-1] and ur[i] > ur[i-1]:
            vz[i] = np.pi-uz[i] ; vr[i] = ur[i]
         elif uz[i] < uz[i-1] and ur[i] < ur[i-1]:
            vz[i] = np.pi-uz[i] ; vr[i] = 2*np.pi-ur[i]
         elif uz[i] > uz[i-1] and ur[i] < ur[i-1]:
            vz[i] = 2*np.pi+uz[i] ; vr[i] = 2*np.pi-ur[i]
      else:
         vz[i] = uz[i] ; vr[i] = ur[i]

   # Define vr as deviation from vz
   vr = vr-vz ; vr[-1] = vr[0]

   x = vz
   # Now:
   #
   # R = rmaj + rmin*cos(x+vr)
   # Z = zmaj + zmin*sin(x)
   #
   # NOTE: kappa=zmin/rmin

   xr = np.zeros(4)
   xr[0] = rmin
   xr[1] = rmaj
   xr[2] = zmin/rmin
   xr[3] = zmaj

   for i in range(narc-1):
      dx[i] = x[i+1]-x[i]

   dx[-1] = dx[0]

   # Compute expansion coefficients (cr):
   #  vr = sum cr*cos(nx)+sr*sin(nx)
   cr = np.zeros(nf+1)
   sr = np.zeros(nf+1)

   cr[0] = moment(narc,vr,np.ones(narc),dx)
   for p in range(1,nf+1):
      cr[p] = moment(narc,vr,np.cos(p*x),dx)
      sr[p] = moment(narc,vr,np.sin(p*x),dx)

   if xplot > 0.0:
      outfile = '{:.3f}'.format(xplot)+'.png'
      plot_ang(r,z,x,vr,xr,cr,sr,outfile)

   return cr,sr,xr

def prgen_fshape(rd,zd,nf):

    nd = len(rd)

    # Construct equally-spaced poloidal angle
    theta  = np.linspace(0,1,nd)*2*np.pi
    dtheta = theta[1]-theta[0]

    ar = np.zeros(nf+1) ; br = np.zeros(nf+1)
    az = np.zeros(nf+1) ; bz = np.zeros(nf+1)

    ds = dtheta/np.pi

    # Trapezoidal integration spectrally-accurate
    for n in range(nf+1):
       y = np.cos(n*theta[:])*rd[:] ; ar[n] = ds*np.sum(y[:-1])
       y = np.sin(n*theta[:])*rd[:] ; br[n] = ds*np.sum(y[:-1])
       y = np.cos(n*theta[:])*zd[:] ; az[n] = ds*np.sum(y[:-1])
       y = np.sin(n*theta[:])*zd[:] ; bz[n] = ds*np.sum(y[:-1])

    return ar,br,az,bz

#-----------------------------------------------------------------------------
# Old Fourier expansion
#-----------------------------------------------------------------------------
def oldfourier(ri,zi,nf,rnorm):

   print("INFO: (oldfourier) Generating legacy Fourier coefficients")
   npsi = len(rnorm)

   ari = np.zeros([nf+1,npsi]) ; bri = np.zeros([nf+1,npsi])
   azi = np.zeros([nf+1,npsi]) ; bzi = np.zeros([nf+1,npsi])

   for i in range(npsi-1):
      r=ri[:,i+1] ; z=zi[:,i+1]
      ar,br,az,bz = prgen_fshape(r,z,nf)
      ari[:,i+1] = ar[:] ; bri[:,i+1] = br[:]
      azi[:,i+1] = az[:] ; bzi[:,i+1] = bz[:]

   # Repair origin
   ari[0,:] = extrap(rnorm,ari[0,:])
   azi[0,:] = extrap(rnorm,azi[0,:])
   for i in range(1,nf+1):
      ari[i,:] = zero(rnorm,ari[i,:]) ; bri[i,:] = zero(rnorm,bri[i,:])
      azi[i,:] = zero(rnorm,azi[i,:]) ; bzi[i,:] = zero(rnorm,bzi[i,:])

   u = ari
   u = np.append(u,bri)
   u = np.append(u,azi)
   u = np.append(u,bzi)
   u.tofile('fluxfit.geo')

   return
