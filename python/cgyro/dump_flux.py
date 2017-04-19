import sys
import numpy as np
from gacodeplotdefs import *
from cgyro.data import cgyrodata

sim = cgyrodata('./')
sim.getflux()

ns = sim.n_species
t = sim.t

ys = np.sum(sim.ky_flux,axis=(2,3))
for moment in ['n','e','v']:
    if moment == 'n':
        ntag = 'Density~flux'
        mtag = '\Gamma'
        ttag = 'G'
        ftag = 'flux_n'
        y = ys[:,0,:]
    elif moment == 'e':
        ntag = 'Energy~flux'
        mtag = 'Q'
        ttag = 'Q'
        ftag = 'flux_e'
        y = ys[:,1,:]
    elif moment == 'v':
        ntag = 'Momentum~flux'
        mtag = '\Pi'
        ttag = 'Pi'
        ftag = 'flux_v'
        y = ys[:,2,:]

    data  = np.column_stack((sim.t,y[0,:]))
    head  = '(cs/a) t     '+ttag+'_1/'+ttag+'_GB'
    fname = 'out.cgyro.'+ftag
    for ispec in range(1,ns,1):
        head = head+'       '+ttag+'_'+str(ispec+1)+'/'+ttag+'_GB'
        data = np.column_stack((data,y[ispec,:]))
    np.savetxt(fname,data,fmt='%.8e',header=head)
    print 'INFO: (dump_flux.py) Created '+fname
