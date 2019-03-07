# Dump text report of flux

import sys
import numpy as np
from gacodefuncs import *
from cgyro.data import cgyrodata

def print_freq():
    # Set print precision
    np.set_printoptions(precision=5,suppress=True)

    print('   omega    gamma')
    for i in range(nt):
        print(sim.freq[:,0,i])
    
def print_flux():
    # Set print precision
    np.set_printoptions(precision=4,suppress=True)

    b = np.zeros([sim.n_species])
 
    tag = [
       'GAMMA [GB]',
       'Q     [GB]',
       'PI    [GB]']

    sim.getflux('auto')

    # Determine imin
    imin,imax=iwindow(sim.t,w,wmax)

    print('INFO: (text.py) Average Window:',str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[imax]))

    title = '        '
    for ispec in range(sim.n_species):
        title = title+'       '+specmap(sim.mass[ispec],sim.z[ispec])
    print(title)

    for i in range(3):
        try:
            for ispec in range(sim.n_species):
                y = np.sum(sim.ky_flux,axis=(2,3))
                b[ispec] = average(y[ispec,i,:],sim.t,w,wmax)
            print(tag[i],b)
        except:
            pass

#-------------------------------------------------------------------
        
w    = float(sys.argv[1])
wmax = float(sys.argv[2])
ext  = sys.argv[3]

sim = cgyrodata('./')

nt = sim.n_time

if sim.n_n > 1:
    # Print flux (assuming nonlinear case)
    print_flux()
else:
    # Pring frequency (assuming linear)
    print_freq()
        
