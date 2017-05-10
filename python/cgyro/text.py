# Dump text report of flux

import sys
import numpy as np
from gacodefuncs import *
from cgyro.data import cgyrodata

w = float(sys.argv[1])
ext = sys.argv[2]

sim = cgyrodata('./')
sim.getflux()

nt = sim.n_time

# Determine imin
imin=iwindow(sim.t,w)

if imin == nt-1:
    print "Averaging Window too small." 
    sys.exit()

print 'INFO: Average Window:',str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])

if ext == 'dump':
    # Datafile output

    b = np.zeros([2,sim.n_n])

    fn=open('out.cgyroplot.fluxn','w')
    fe=open('out.cgyroplot.fluxe','w')
    for ispec in range(sim.n_species):
        fn.write(str(ispec)+' ')
        fe.write(str(ispec)+' ')
        for i_n in range(sim.n_n):
            y = np.sum(sim.ky_flux,axis=2)
            yn = y[:,0,:,:]
            ye = y[:,1,:,:]
            b[0,i_n] = average(yn[ispec,i_n,:],sim.t,w)
            b[1,i_n] = average(ye[ispec,i_n,:],sim.t,w)
            fn.write("{:.3e}".format(b[0,i_n])+' ')
            fe.write("{:.3e}".format(b[1,i_n])+' ')
        fn.write('\n')
        fe.write('\n')            
        print 'Wrote output to out.cgyroplot.flux*'
else:
    # Screen output

    # Set print precision
    np.set_printoptions(precision=3,suppress=False,threshold=100000)

    b = np.zeros([sim.n_species])

    try:
        for ispec in range(sim.n_species):
	    y = np.sum(sim.ky_flux,axis=(2,3))
            b[ispec] = average(y[ispec,0,:],sim.t,w)
        print 'GAMMA [GB]  ',b
    except:
        pass
            

    try:
        for ispec in range(sim.n_species):
            y = np.sum(sim.ky_flux,axis=(2,3))
            b[ispec] = average(y[ispec,1,:],sim.t,w)
        print 'Q     [GB]  ',b
    except:
        pass



    
