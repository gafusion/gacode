import sys
import numpy as np
from gacodeplotdefs import *
from gacodefuncs import *
from gyro.data import GYROData

sim       = GYROData(sys.argv[1])
field     = sys.argv[2]
i_moment  = int(sys.argv[3])
w         = float(sys.argv[4])
plotfile  = sys.argv[5]
datafile  = sys.argv[6]
lx        = float(sys.argv[7])
ly        = float(sys.argv[8])
title     = sys.argv[9]
ymin      = sys.argv[10]
ymax      = sys.argv[11]
xspan1    = float(sys.argv[12])
xspan2    = float(sys.argv[13])

n_field   = int(sim.profile['n_field'])
n_kinetic = int(sim.profile['n_kinetic'])

t = sim.t['(c_s/a)t']

# Read data in gbflux_i and make gbflux
sim.read_gbflux_i()
sim.make_gbflux()

flux = sim.gbflux

# Manage field
if field == 's':
    flux0 = np.sum(flux,axis=1)
    ftag  = sim.tagfield[3]
else:
    i_field = int(field)
    flux0   = flux[:,i_field,:,:]
    ftag    = sim.tagfield[i_field]

# Manage moment
mtag = sim.tagmom[i_moment]

#======================================
fig = plt.figure(figsize=(lx,ly))
fig.subplots_adjust(left=0.1,right=0.95,top=0.92,bottom=0.12)
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$(c_s/a) t$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$'+mtag+' \;('+ftag+')$',color='k',fontsize=GFONTSIZE)
#=====================================

cvec = ['k','m','b','c','g']

# Determine tmin
for i in range(len(t)):
    if t[i] < (1.0-w)*t[len(t)-1]:
        imin = i

if title=='null':
    ax.set_title(r'$'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$')
else:
    ax.set_title(r'$\mathrm{'+title+'}$')
    
# Loop over species
if datafile == 'none':
    # Plot data to screen or image file.
    for i in range(n_kinetic):
        ave   = average(flux0[i,i_moment,:],t,w)
        if i > n_kinetic-2 and sim.profile['electron_method'] > 1:
            stag = r'$e'
            color = 'r'
        else:
            stag = r'$i_'+str(i+1)
            color = cvec[i]

        label = stag+' : '+str(round(ave,3))+'$'
        y     = ave*np.ones(len(t))
        ax.plot(t[imin:],y[imin:],'--',color=color)
        ax.plot(t,flux0[i,i_moment,:],label=label,color=color)
else:
    # Write data to datafile
    print 'INFO: (gyro_plot) Output to datafile not supported.  Use raw out.gyro.gbflux.'

if xspan1 > 0.0:
    ax.axvspan(xspan1,xspan2,facecolor='g',alpha=0.1)

ax.set_xlim([0,t[-1]])

if ymax != 'auto':
    ax.set_ylim([float(ymin),float(ymax)])
        
ax.legend(loc=1)

if plotfile == 'screen':
    plt.show()
else:
    plt.savefig(plotfile)
    print "INFO: (gyro_plot) Wrote plot to " + plotfile + "."
