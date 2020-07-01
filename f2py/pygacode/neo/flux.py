"""flux.py is executed by the bash script neo_plot."""

import sys
import numpy as np
import matplotlib.pyplot as plt
from .data import NEOData

GFONTSIZE=16

sim    = NEOData(sys.argv[1])
moment = int(sys.argv[2])
ftype  = sys.argv[3]

n_species = sim.grid['n_species']
n_radial  = sim.grid['n_radial']

x = sim.equil['r_over_a']

y = np.zeros([n_radial,n_species])

if moment == 0:
    y[:,:] = sim.transport['Gamma'][:,:]
    label='$\Gamma$'
elif moment == 1:
    y[:,:] = sim.transport['Pi'][:,:]
    label='$\Pi$'
elif moment == 2:
    y[:,:] = sim.transport['Q'][:,:]
    label='$Q$'

#======================================
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax.grid(which="both",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel('$r/a$',fontsize=GFONTSIZE)
ax.set_ylabel(label,fontsize=GFONTSIZE)
ax.set_title('Title Goes Here',fontsize=GFONTSIZE)
#=====================================

# Loop over species
for i in range(n_species):
    ax.plot(x,y[:,i],label='Species '+str(i))

ax.set_xlim([min(x),max(x)])
ax.legend()

if ftype == 'screen':
    plt.show()
else:
    outfile = 'flux.'+ftype
    plt.savefig(outfile)
