import sys
import numpy as np
from gacodeplotdefs import *
from gacodefuncs import *
from gyro.data import GYROData

sim      = GYROData(sys.argv[1])
field    = sys.argv[2]
i_moment = int(sys.argv[3])
window   = float(sys.argv[4])
ftype    = sys.argv[5]
ymin     = sys.argv[6]
ymax     = sys.argv[7]

n_field   = int(sim.profile['n_field'])
n_kinetic = int(sim.profile['n_kinetic'])

t = sim.t['(c_s/a)t']

# Read data in gbflux_i
sim.read_gbflux_i()

flux = sim.gbflux_i

# Manage field
if field == 's':
    flux0 = np.sum(flux,axis=1)
    ftag  = sim.tagfield[3]
else:
    i_field = int(field)
    flux0 = flux[:,i_field,:,:,:]
    ftag  = sim.tagfield[i_field]

# Manage moment
mtag = sim.tagmom[i_moment]

# Determine tmin
for i in range(len(t)):
    if t[i] < (1.0-window)*t[len(t)-1]:
        imin = i

#======================================
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$r/a$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$'+mtag+' \;('+ftag+')$',color='k',fontsize=GFONTSIZE)
ax.set_title(r'$'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$',fontsize=GFONTSIZE)
#=====================================

color = ['k','m','b','c']

n_x = sim.profile['n_x']
ave = np.zeros(n_x)

# Loop over species
for i in range(n_kinetic):
    stag = sim.tagspec[i]
    ave[:] = average_n(flux0[i,i_moment,:,:],t,window,n_x)
    ax.plot(sim.profile['r'],ave[:],label=stag,color=color[i])

if ymax != 'auto':
    ax.set_ylim([float(ymin),float(ymax)])

ax.legend()

if ftype == 'screen':
    plt.show()
else:
    outfile = 'gbflux_i.'+ftype
    plt.savefig(outfile)
