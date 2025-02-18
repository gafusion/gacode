'''
Cooling curves for multiple ions from Aurora
No charge exchange included. Simple ionization equilibrium.
sciortino, 2021
'''

import aurora
import numpy as np
import numpy.polynomial.chebyshev as pc

#------------------------------------------------------------------------
# Inputs

ions_list = ['H','He','Be','C','O','N','F','Ne','Al','Si','Ar','Ca','Fe','Ni','Kr','Mo','Xe','W','Li']

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

x = pc.chebpts1(nc)

te_ev = np.exp(0.5*dt*(x+1)+lte0)
ne_cm3 = np.ones_like(te_ev)

for imp in ions_list:
    rfile = pue_base+f'plt_caic_mix_{imp}.dat'
    line,cont = aurora.get_cooling_factors(
        imp, ne_cm3, te_ev,
        line_rad_file=pue_base+f'plt_caic_mix_{imp}.dat',
        cont_rad_file=pue_base+f'prb_{imp}.dat'
        )    
    ofile = imp+'.txt'
    print('Writing '+imp)
    # NOTE: total cooling curve is cont+line
    np.savetxt(ofile,line+cont)

np.savetxt('te.txt',te_ev)
np.savetxt('x.txt',x)

