import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from .prgen_geqdsk import *
from .prgen_contour import *
from .prgen_shape import *

if len(sys.argv) > 1:
   gfile   = sys.argv[1]
   nrz     = int(sys.argv[2])
   narc    = int(sys.argv[3])
   npsi    = int(sys.argv[4])
   nfourier = int(sys.argv[5])
else:
   print('Usage: python prgen_shapeprofile.py <gfile> <nrz> <narc> <npsi> <nfourier>')
   sys.exit()

efit = prgen_geqdsk(gfile)
nf = 3

# Call OMFIT mapper
ri,zi,psi,q,p,fpol = prgen_contour(efit,nrz=nrz,levels=npsi,psinorm=0.996,narc=narc,quiet=False)

pnorm = ((psi[:]-psi[0])/(psi[-1]-psi[0]))
rnorm = np.sqrt(pnorm)

# ci -> cos terms
# si -> sin terms
# xi -> rmin,rmaj,kappa,zmaj

# Standard contour mode
ci = np.zeros([nf+1,npsi])
si = np.zeros([nf+1,npsi])
xi = np.zeros([4,npsi])
for i in range(npsi-1):
   r=ri[:,i+1] ; z=zi[:,i+1]
   if (np.mod(i,4) == 0) and (nfourier > 0):
      xplot = pnorm[i+1]
   else:
      xplot = 0.0
   cr,sr,xr = prgen_shape(r,z,narc,nf,xplot)
   ci[:,i+1] = cr[:]
   si[:,i+1] = sr[:]
   xi[:,i+1] = xr[:]

# Repair functions near origin

xi[0,:] = zero(rnorm,xi[0,:]) # rmin
xi[1,:] = extrap(pnorm,xi[1,:]) # rmaj
xi[2,:] = extrap(rnorm,xi[2,:]) # kappa
xi[3,:] = extrap(rnorm,xi[3,:]) # zmaj

si[1,:] = zero(rnorm,si[1,:]) # delta
si[2,:] = zero(pnorm,si[2,:]) # zeta
si[3,:] = zero(pnorm,si[3,:]) # s3

ci[0,:] = extrap(pnorm,ci[0,:]) # tilt
ci[1,:] = zero(rnorm,ci[1,:]) # c1
ci[2,:] = zero(pnorm,ci[2,:]) # c2
ci[3,:] = zero(pnorm,ci[3,:]) # c3

f=open('out.dim','w')
f.write(str(npsi)+'\n')
f.write(str(nf+1)+'\n')
f.write(str(efit['RCENTR'])+'\n')
f.write(str(efit['BCENTR'])+'\n')
f.write(str(efit['CURRENT']*1e-6)+'\n')
f.close()
u = psi
u = np.append(u,q)
u = np.append(u,p)
u = np.append(u,si[:,:])
u = np.append(u,ci[:,:])
u = np.append(u,xi[:,:])
u.tofile('out.data')

# Finished standard mode
# If nfourier > 0, generate Fourier coefficients and diagnostics

if nfourier > 0:
   # Generate pure Fourier expansion coefficients, write to fluxfit.geo
   oldfourier(ri,zi,nfourier,rnorm)
   # Plot radial profiles of cos,sin,etc
   plot_coef(pnorm,ci,si,xi)
