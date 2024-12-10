import os
import numpy as np
#from gacodefuncs import *
from cgyro.data import cgyrodata
import matplotlib.pyplot as plt
from pydmd import DMD
import matplotlib.gridspec as gridspec

DTYPE='float32'

def read_binary_file(file_path, dtype=np.complex64):

    with open(file_path, 'rb') as file:
        file_content = file.read()
        
    type_size = np.dtype(dtype).itemsize
    
    num_elements = len(file_content) // type_size
    
    data = np.frombuffer(file_content, dtype=np.complex64, count=num_elements)
    
    return data

mydir='/home/candy/reg03/'

# main
filename      = mydir+"bin.cgyro.kxky_phi"
dens_filename = mydir+"bin.cgyro.kxky_n"
apar_filename = mydir+"bin.cgyro.kxky_apar"

# optional 
bpar_filename = mydir+"bin.cgyro.kxky_bpar"
v_filename    = mydir+"bin.cgyro.kxky_v"

sim = cgyrodata(mydir)
sim.getbigfield()

print('ky*rho = ', sim.ky0)
#print('omega = ', sim.freq[0,0,-1])
#print('gamma = ', sim.freq[1,0,-1])

t = sim.t
N_radial = sim.n_radial   # CGYRO N_RADIAL
N_theta = sim.theta_plot   # CGYRO N_THETA
Nspecies = sim.n_species   # number of gyrokinetic species 

sim_fields = np.zeros([len(t), N_theta, N_radial], dtype=DTYPE)
sim_ptcls = np.zeros([len(t), N_theta, N_radial, Nspecies], dtype=DTYPE)

phi_kxky  = sim.kxky_phi[0,:,:,0,:] +1j*sim.kxky_phi[1,:,:,0,:]
apar_kxky = sim.kxky_apar[0,:,:,0,:]+1j*sim.kxky_apar[1,:,:,0,:]
bpar_kxky = sim.kxky_bpar[0,:,:,0,:]+1j*sim.kxky_bpar[1,:,:,0,:]

data = read_binary_file(dens_filename)
dens_kxky = data.reshape(sim_ptcls.shape)

data = read_binary_file(v_filename)
v_kxky = data.reshape(sim_ptcls.shape)

# SVD rank to perform DMD 
svd_rank = 0

# Switch to DMD time step
k=1

# step for DMD   
delt = t[::k][1]-t[::k][0]

#-----------------------------------------------------------------------------------------
# fluctuating potential 

dmd = DMD(svd_rank = svd_rank, exact = True)

fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = fig.add_subplot(gs[0,0])

ax1.axvline(0., linestyle="dashed", color="k", linewidth=0.5)
ax1.axhline(0., linestyle="dashed", color="k", linewidth=0.5)

for j in range(N_radial):

    dmd.fit(phi_kxky[j,:,::k])
    realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)

    ax1.plot(realEigs.real, realEigs.imag, 'o', color="tab:blue")

#-----------------------------------------------------------------------------------------
# fluctuating density 1

for j in range(N_radial):

    dmd.fit(dens_kxky[::k,:,j,0].T)
    realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)

    ax1.plot(realEigs.real, realEigs.imag, 's', color="tab:green", mfc="none", markersize=9.5)

#-----------------------------------------------------------------------------------------
# fluctuating vector potential 

for j in range(N_radial):

    dmd.fit(apar_kxky[j,:,::k])
    realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)

    ax1.plot(realEigs.real, realEigs.imag, 'o', color="red", mfc="none")

#-----------------------------------------------------------------------------------------
# fluctuating magnetic field  

for j in range(N_radial):

    dmd.fit(bpar_kxky[j,::k,:])
    realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)

    #ax1.plot(realEigs.real, realEigs.imag, 's', color="red", mfc="none", markersize=7.5)

ax1.set_xlabel(r"realEigs.real, i.e. $\omega$", fontsize=15)
ax1.set_ylabel(r"realEigs.imag, i.e. $\gamma$", fontsize=15)
ax1.tick_params(labelsize=15)

plt.show()

#-----------------------------------------------------------------------------------------
# fluctuating v 

#dmd = DMD(svd_rank = svd_rank, exact = True)
#dmd.fit(v_kxky[:,:,j,1].T)

#realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)

#ax1.axvline(0., linestyle="dashed", color="k", linewidth=0.5)
#ax1.axhline(0., linestyle="dashed", color="k", linewidth=0.5)

#ax1.plot(realEigs.real, realEigs.imag, '^', color="cyan", mfc="none")

#ax1.set_xlabel(r"realEigs.real, i.e. $\omega$", fontsize=15)
#ax1.set_ylabel(r"realEigs.imag, i.e. $\gamma$", fontsize=15)
#ax1.tick_params(labelsize=15)

#plt.show()






