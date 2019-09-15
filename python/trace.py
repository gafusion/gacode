from omfit.classes.omfit_eqdsk import OMFITeqdsk
import sys
import numpy as np
import gapy
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text',usetex=True)
rc('font',size=18)

gfile=''          # list of g-files to trace the flux-surfaces of
resolution=0.002  # change resolution with bicubic-spline interpolation (default 0)
uniform=200       # uniform number of points in theta angle (default 0)
levels=16         # number of flux surfaces (default 0, same as gfile psi grid)
maxPSI=0.9999     # normalized PSI/RHO at which to place the separatrix
rhoVSpsi='psi'    # levels are specified in `rho` or `psi` (default `psi`)
filetype='none'   # dump flux surfaces traces to file `flux_xxx` (default 1)

gfile = sys.argv[1]
ix = int(sys.argv[2])

EQDSK = OMFITeqdsk(gfile)

if levels > 0:
   #read the gfile and calculate flux surfaces at given resolution
   EQDSK.addFluxSurfaces(resolution=resolution,
                         map=rhoVSpsi,
                         forceFindSeparatrix=True,
                         quiet=False,
                         calculateAvgGeo=True,
                         maxPSI=float(maxPSI),
                         levels=levels)
else:
   #read the gfile and calculate flux surfaces at given resolution
   EQDSK.addFluxSurfaces(resolution=resolution,
                         map=rhoVSpsi,
                         forceFindSeparatrix=True,
                         quiet=False,
                         calculateAvgGeo=True,
                         maxPSI=float(maxPSI))

if uniform > 0:
   print('INFO (fluxSurfaceTracer) Resampling around contour')
   EQDSK['fluxSurfaces'].resample(uniform)

if filetype == 'ascii':
   print('INFO (fluxSurfaceTracer) output = '+outfile)
   EQDSK['fluxSurfaces']['flux'].deploy(outfile)

r = np.zeros([levels,uniform])
z = np.zeros([levels,uniform])
for i in range(levels):
   r[i,:] = EQDSK['fluxSurfaces']['flux'][i]['R']
   z[i,:] = EQDSK['fluxSurfaces']['flux'][i]['Z']

gapy.surfpar.nf=4
gapy.surfpar.iselect=ix+1

gapy.surfpar.surfpar_do(r,z)

#--------------------------------------------------------------------
# PLOTTING

fig = plt.figure(figsize=(13,8))

# PLOT 1
ax = fig.add_subplot(121,aspect='equal')
ax.set_xlabel(r'$R$')
ax.set_ylabel(r'$Z$')

ax.plot(gapy.surfpar.r,gapy.surfpar.z,'-k',linewidth=2)
ax.plot(gapy.surfpar.rp,gapy.surfpar.zp,'-m',linewidth=2)
ax.plot(gapy.surfpar.r[0],gapy.surfpar.z[0],'o',linewidth=1)

# PLOT 2
ax = fig.add_subplot(122)
ax.set_xlabel('$\ell$')

y = gapy.surfpar.x/(2*np.pi)

ax.set_xlim([0,1])
ax.plot(y,gapy.surfpar.vr,'-r',linewidth=1,label=r'$\theta_R$')
ax.plot(y,gapy.surfpar.vz,'-b',linewidth=1,label=r'$\theta_Z$')
ax.plot(y,gapy.surfpar.pr,'--r',linewidth=2)
ax.plot(y,gapy.surfpar.pz,'--b',linewidth=2)

ax.legend(ncol=2,loc=1)

plt.tight_layout()

plt.show()
