import os
import numpy as np
#from gacodefuncs import *
from cgyro.data import cgyrodata
import matplotlib.pyplot as plt
from matplotlib import rc
from pydmd import DMD
import matplotlib.gridspec as gridspec


#---------------------------------------------------------------------------
# INPUTS
mydir='/home/candy/reg03/'

# time downsample
k = 2
# theta downsample
l = 6

# SVD rank to perform DMD 
svd_rank = 0

#---------------------------------------------------------------------------
# COLLECT DATA
sim = cgyrodata(mydir)
sim.getbigfield()

t = sim.t
n_radial  = sim.n_radial  
n_theta   = sim.theta_plot
n_species = sim.n_species   

ovec = {}
ovec['phi']  = sim.kxky_phi[0,:,:,0,:] +1j*sim.kxky_phi[1,:,:,0,:]
ovec['apar'] = sim.kxky_apar[0,:,:,0,:]+1j*sim.kxky_apar[1,:,:,0,:]
ovec['bpar'] = sim.kxky_bpar[0,:,:,0,:]+1j*sim.kxky_bpar[1,:,:,0,:]

# step for DMD
dt = k*(t[1]-t[0])

#---------------------------------------------------------------------------
# RUN DMD
dmd = DMD(svd_rank=svd_rank,exact=True)

evec = {}
for x in ['phi','apar','bpar']:
    dmd.fit(ovec[x][:,::l,::k])
    evec[x] = 1j*np.log(dmd.eigs)/dt

#---------------------------------------------------------------------------
# PLOTTING

rc('font',size=25)
rc('text',usetex=True)

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)

ax.set_xlim([-1.2,1.2])
ax.set_ylim([-0.25,0.5])
ax.set_xlabel(r"$\omega$")
ax.set_ylabel(r"$\gamma$")
ax.grid(which="both",ls=":")
ax.grid(which="major",ls=":")

# symbols
args = {}
args['phi']  = {'color':'r','marker':'o','alpha':0.4}
args['apar'] = {'color':'b','marker':'s','facecolors':'none'}
args['bpar'] = {'color':'k','marker':'+'}

for x in ['phi','apar','bpar']:
    ax.scatter(evec[x].real,evec[x].imag,s=50,**args[x])

plt.tight_layout()
plt.show()
