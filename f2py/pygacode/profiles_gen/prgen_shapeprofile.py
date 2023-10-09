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
   mag     = sys.argv[2]
   narc    = int(sys.argv[3])
   npsi    = int(sys.argv[4])
   nharm   = int(sys.argv[5])
   nfourier = int(sys.argv[6])
   plotpng = bool(int(sys.argv[7]))
   psinorm = float(sys.argv[8])
else:
   print('Usage: python prgen_shapeprofile.py <gfile> <mag> <narc> <npsi> <nfourier>')
   sys.exit()

efit = prgen_geqdsk(gfile)
psi0 = efit['SIMAG']
psi1 = efit['SIBRY']

# Call GACODE mapper (uses n=npsi-1 and excludes magnetic axis)
ri,zi,psi,p,fpol,q,psi_sep = prgen_contour(efit,mag=mag,nc=npsi,psinorm=psinorm,narc=narc)

pnorm = (psi[:]-psi0)/(psi1-psi0)
rnorm = np.sqrt(pnorm)

# ci -> cos terms
# si -> sin terms
# xi -> rmin,rmaj,kappa,zmaj

nf = nharm

# Standard contour mode
ci = np.zeros([nf+1,npsi])
si = np.zeros([nf+1,npsi])
xi = np.zeros([4,npsi])
for i in range(1,npsi):
   r=ri[:,i] ; z=zi[:,i]
   if i == npsi-1:
      xplot = pnorm[i]
   else:
      xplot = 0.0
   cr,sr,xr = prgen_shape(r,z,narc,nf,xplot)
   # Use i+1 to make room for magnetic axis (i=0)
   ci[:,i] = cr[:]
   si[:,i] = sr[:]
   xi[:,i] = xr[:]

# Repair functions near origin

xi[0,0] = 0.0                   # rmin
xi[1,:] = extrap(pnorm,xi[1,:]) # rmaj
xi[2,:] = extrap(rnorm,xi[2,:]) # kappa
xi[3,:] = extrap(rnorm,xi[3,:]) # zmaj

ci[0,:] = extrap(pnorm,ci[0,:]) # tilt

si[1:,0] = 0.0 # delta
ci[1:,0] = 0.0 # c1

with open('out.dim','w') as f:
   f.write(str(npsi)+'\n')
   f.write(str(nf+1)+'\n')
   f.write(str(efit['RCENTR'])+'\n')
   f.write(str(efit['BCENTR'])+'\n')
   f.write(str(efit['CURRENT']*1e-6)+'\n') # Convert Amps [EFIT] to MA
   f.write(str(psi0)+'\n')
   f.write(str(psi1)+'\n')
   f.write(str(psi_sep)+'\n')

u = psi
u = np.append(u,q)
u = np.append(u,p)
u = np.append(u,fpol)
u = np.append(u,si[:,:])
u = np.append(u,ci[:,:])
u = np.append(u,xi[:,:])
u.tofile('out.data')

# Finished standard mode
# If nfourier > 0, generate Fourier coefficients and diagnostics

if nfourier > 0:
   # Generate pure Fourier expansion coefficients, write to fluxfit.geo
   oldfourier(ri,zi,nfourier,rnorm)
   if plotpng:
      # Plot radial profiles of cos,sin,etc
      plot_coef(pnorm,ci,si,xi)
