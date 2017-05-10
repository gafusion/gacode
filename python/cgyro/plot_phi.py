import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
field = sys.argv[2]

sim = cgyrodata('./')
sim.getbigfield()

#======================================
# Set figure size and axes
fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.09,right=0.96,top=0.92,bottom=0.12)
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(TIME)
ax.set_ylabel(r'$\left| \Phi \right|$')
ax.set_yscale('log')
ax.set_title(r'$\mathrm{Fluctuation~intensity} \quad k_\theta = nq/r$')
#======================================

p2 = np.sum(sim.phisq,axis=0)/sim.rho**2
itheta=0

# n=0 intensity
y0 = p2[itheta,0,:]
ax.plot(sim.t,np.sqrt(y0),label=r'$n=0$',linewidth=2)

# finite-n intensity
yn = p2[itheta,1,:]
for n in range(2,sim.n_n):
    yn = yn+p2[itheta,n,:]
ax.plot(sim.t,np.sqrt(yn),label=r'$n>0$')
        
ax.set_xlim([0,max(sim.t)])

ax.legend(loc=4)

fname = 'out.cgyro.phi'

if ftype == 'screen':
    plt.show()
elif ftype == 'dump':
    data = np.column_stack((sim.t,np.sqrt(y0)))
    data = np.column_stack((data,np.sqrt(yn)))
    head = '(cs/a) t     Phi_0/rho_*    Phi_n/rho_*'
    np.savetxt(fname,data,fmt='%.8e',header=head)
    print 'INFO: (plot_phi) Created '+fname 
else:
    fname=fname+'.'+ftype
    print 'INFO: (plot_phi) Created '+fname
    plt.savefig(fname)
