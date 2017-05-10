import sys
import numpy as np
from gacodeplotdefs import *
from gacodefuncs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
field = sys.argv[2]
ymin  = sys.argv[3]
ymax  = sys.argv[4]
nstr  = sys.argv[5]

sim = cgyrodata('./')
sim.getbigfield()

#============================================================
# Set figure size and axes
fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.08,right=0.96,top=0.92,bottom=0.12)
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(TIME)
ax.set_ylabel(r'$\left| \delta\phi_n \right|$')
ax.set_yscale('log')
ax.set_title(r'$\mathrm{Fluctuation~intensity} \quad k_\theta = nq/r$')

p2 = np.sum(sim.phisq[:,0,:,:],axis=0)/sim.rho**2

if nstr == 'null':
    nvec = range(sim.n_n)
else:
    nvec = str2list(nstr)
    
for n in nvec:
    num = r'$n='+str(n)+'$'
    if n==0:
        ax.plot(sim.t,np.sqrt(p2[n,:]),linewidth=2,label=num)
    else:
        ax.plot(sim.t,np.sqrt(p2[n,:]),label=num)

ax.set_xlim([0,max(sim.t)])

if sim.n_n > 16:
    ax.legend(loc=4, ncol=5, prop={'size':12})
else:
    ax.legend(loc=4, ncol=6, prop={'size':12})

ymin,ymax=setlimits(ax.get_ylim(),ymin,ymax)

ax.set_ylim([ymin,ymax])

fname = 'out.cgyro.ky_phi'

if ftype == 'screen':
    plt.show()
else:
    fname=fname+'.'+ftype
    print 'INFO: (plot_phi) Created '+fname
    plt.savefig(fname)
