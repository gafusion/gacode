# 2025 version:
# CHANGES: module load omfit, use 5e13 for density, add Li and B.
# Thanks to Tomas

import aurora
import numpy as np
import numpy.polynomial.chebyshev as pc

#------------------------------------------------------------------------
# Inputs

ion_vec = ['H','He','Be','C','O','N','F','Ne','Al','Si','Ar','Ca','Fe','Ni','Kr','Mo','Xe','W','Li','B']

# Temperature limits in keV
te_max = 50.0
te_min = 0.05

# Chebyshev nodes
nc = 64
#-----------------------------------------------------------

pue_base = '/fusion/projects/toolbox/sciortinof/atomlib/atomdat_master/pue2021_data/'

# log-temperature bounds (convert to eV)
lte1 = np.log(te_max*1e3)
lte0 = np.log(te_min*1e3)
dt   = lte1-lte0

# Chebyshev nodes
x = pc.chebpts1(nc)

te_ev = np.exp(0.5*dt*(x+1)+lte0)

# Updated for 2025: 5e13/cm^3
ne_cm3 = 5e13*np.ones(nc)

for i in ion_vec:
    line,cont = aurora.get_cooling_factors(
        i, ne_cm3, te_ev,
        line_rad_file=pue_base+f'plt_caic_mix_{i}.dat',
        cont_rad_file=pue_base+f'prb_{i}.dat',
        plot=False
        )    
    print('Writing '+i)
    # Total cooling is line+cont
    np.savetxt(i+'.txt',line+cont)

np.savetxt('te.txt',te_ev)
np.savetxt('x.txt',x)

