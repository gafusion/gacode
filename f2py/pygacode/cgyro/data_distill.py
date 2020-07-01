import sys
import os
import numpy as np
from gacodefuncs import *
from .data import cgyrodata

if len(sys.argv) < 2:
    print("Usage: cgyro_distill.py <window> <theta_index>")
    sys.exit()

if not os.path.isfile('bin.cgyro.kxky_phi'):
    print("Error: no data found")
    sys.exit()

w = float(sys.argv[1])
theta = float(sys.argv[2])
                 
# read basic simulation data
sim = cgyrodata('./')

# read larger files
sim.getbigfield()

nt  = sim.n_time
nr  = sim.n_radial
nn  = sim.n_n
ns  = sim.n_species
nth = sim.theta_plot

t = sim.t

# Theta-slide index
i0 = theta_index(theta,nth)

# Get index for average window
imin,imax=iwindow(t,w,0.0)

header = 'Averaging over '+str(t[imin])+' < (c_s/a) t < '+str(t[imax])
header = header+' | nkx = '+str(nr)+' | nky = '+str(nn)+'\n'

# kxky_phi -> [2,self.n_radial,self.theta_plot,self.n_n,nt]

# 1. Average phi^2/rho_{*D}
phi2 = (sim.kxky_phi[0,:,i0,:,:]**2+sim.kxky_phi[1,:,i0,:,:]**2)/sim.rho**2

ave_phi2 = np.zeros([nr,nn])
for j in range(nn):
    ave_phi2[:,j] = average_n(phi2[:,j,:],t,w,0.0,nr)

np.savetxt('out.cgyro.phi2kxky',ave_phi2,fmt='%1.6e',header=header)

# 2. Average Re(phi_0)
phi0 = sim.kxky_phi[:,:,i0,0,:]/sim.rho
ave_phi0 = np.zeros([nr,2])
ave_phi0[:,0] = average_n(phi0[0,:,:],t,w,0.0,nr)
ave_phi0[:,1] = average_n(phi0[1,:,:],t,w,0.0,nr)

np.savetxt('out.cgyro.phi0kx',ave_phi0,fmt='%1.6e',header=header)

# Diagnostic print

print('INFO: (data_distill.py) Wrote out.cgyro.phi2kxky = <|phi_k|^2>/rho^2_*D')
print('INFO: (data_distill.py) Wrote out.cgyro.phi0kx   = <Re(phi_0)>/rho_*D, Im(phi_0)>/rho_{D*}')
print('INFO: (data_distill.py) All field normalized to rho_{*D}')
print('INFO: (data_distill.py) Selecting theta-index {} of {}'.format(i0,nth))
print('INFO: (data_distill.py) '+header)
