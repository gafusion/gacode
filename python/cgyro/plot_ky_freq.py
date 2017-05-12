import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]

sim = cgyrodata('./')

fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.085,right=0.97,top=0.92,bottom=0.12)

#======================================
# Omega
ax = fig.add_subplot(121)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$k_y \rho_s$')
ax.set_ylabel(r'$(a/c_s)\, \omega$')

ax.plot(sim.ky,sim.freq[0,:,-1],color='blue')
ax.plot(sim.ky,sim.freq[0,:,-1],"o",color='k')
if len(sim.ky) > 1:
   ax.set_xlim([0,max(sim.ky)])
#======================================

#======================================
# Gamma
ax = fig.add_subplot(122)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$k_y \rho_s$')
ax.set_ylabel(r'$(a/c_s)\, \gamma$')

ax.plot(sim.ky,sim.freq[1,:,-1],color='red')
ax.plot(sim.ky,sim.freq[1,:,-1],"o",color='k')
if len(sim.ky) > 1:
   ax.set_xlim([0,max(sim.ky)])
#======================================

if ftype == 'screen':
    plt.show()
else:
    outfile = 'out.cgyro.ky_freq.'+ftype
    plt.savefig(outfile)
