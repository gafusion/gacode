import sys
import os
import numpy as np
from gacodefuncs import *
from .data import cgyrodata

if not os.path.isfile('bin.cgyro.kxky_phi'):
    print("ERROR: (data_distill) no data found")
    sys.exit()

w = sys.argv[1]
theta = int(sys.argv[2])

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

# Get index for average window
imin,imax=time_index(t,w)

i0,thetapi = indx_theta(theta,nth)

header = 'Averaging over '+str(t[imin])+' < (c_s/a) t < '+str(t[imax])
header = header+' | nkx = '+str(nr)+' | nky = '+str(nn)

# kxky_phi -> [2,self.n_radial,self.theta_plot,self.n_n,nt]

# 1. Average phi^2/rho_{*D}
phi2 = (sim.kxky_phi[0,:,i0,:,:]**2+sim.kxky_phi[1,:,i0,:,:]**2)/sim.rho**2

ave_phi2 = np.zeros([nr,nn])
ave_phi2[:,:] = time_average(phi2[:,:,:],t,imin,imax)

np.savetxt('out.cgyro.phi2kxky',ave_phi2,fmt='%1.6e',header=header)

# 2. Average Re(phi_0)
phi0 = sim.kxky_phi[:,:,i0,0,:]/sim.rho
ave_phi0 = np.zeros([nr,2])
ave_phi0[:,0] = time_average(phi0[0,:,:],t,imin,imax)
ave_phi0[:,1] = time_average(phi0[1,:,:],t,imin,imax)

np.savetxt('out.cgyro.phi0kx',ave_phi0,fmt='%1.6e',header=header)

# Diagnostic print

print('INFO: (data_distill.py) Fields normalized to rho_{*D}')
print('INFO: (data_distill.py) '+header)
print('INFO: (data_distill.py) Wrote out.cgyro.phi2kxky = <|phi_k|^2>/rho^2_*D')
print('INFO: (data_distill.py) Wrote out.cgyro.phi0kx   = <Re(phi_0)>/rho_*D, Im(phi_0)>/rho_{D*}')
