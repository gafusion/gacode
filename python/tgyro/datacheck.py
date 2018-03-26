# Text-based checker for TGYRO data class

import sys
import numpy as np
from tgyro.data import tgyrodata

sim = tgyrodata(sys.argv[1])

print 'Number of ions  : ',sim.n_ion
print 'Number of radii : ',sim.n_r
print 'Evolution eqns  : ',sim.n_evolve
print 'Completed iter  : ',sim.n_iterations

res =sim.data['E(eflux_i)'][-1]+sim.data['E(eflux_e)'][-1]+sim.data['E(mflux)'][-1]+sim.data['E(pflux_e)'][-1]+sim.data['E(pflux_i1)'][-1]+sim.data['E(pflux_i2)'][-1]

print np.sum(res)/(sim.n_evolve*(sim.n_r-1))
p=0
for key in sorted(sim.data.keys()):
    p=p+1
    print p,' : ',key
