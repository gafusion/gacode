import os
import numpy as np
from gacodefuncs import *
from cgyro.data import cgyrodata
import matplotlib.pyplot as plt
from pydmd import DMD
import matplotlib.gridspec as gridspec

DTYPE='float64'


def read_binary_file(file_path, dtype=np.complex64):

    with open(file_path, 'rb') as file:
        file_content = file.read()
        
    type_size = np.dtype(dtype).itemsize
    
    num_elements = len(file_content) // type_size
    
    data = np.frombuffer(file_content, dtype=np.complex64, count=num_elements)
    
    return data


# main
filename = "/path/to/bin.cgyro.kxky_phi"
dens_filename = "/path/to/bin.cgyro.kxky_n"
apar_filename = "/path/to/bin.cgyro.kxky_apar"

# optional 
bpar_filename = "/path/to/bin.cgyro.kxky_bpar"
v_filename = "/path/to/bin.cgyro.kxky_v"

time = np.loadtxt("/path/to/out.cgyro.time")
sim = cgyrodata('/path/to/all_CGYRO_outputs/')



print('ky*rho = ', sim.ky0)
print('omega = ', sim.freq[0,0,-1])
print('gamma = ', sim.freq[1,0,-1])



t = time[:,0]
N_radial = 4   # CGYRO N_RADIAL
N_theta = 32   # CGYRO N_THETA
Nspecies = 2   # number of gyrokinetic species 
theta = np.linspace(-np.pi, np.pi, N_theta * N_radial) * N_radial
theta0 = np.linspace(-np.pi, np.pi, N_theta)



data = read_binary_file(filename)
sim0 = np.zeros([len(t), N_theta, N_radial], dtype=DTYPE)
sim0 = data.reshape(sim0.shape)
phi_kxky = sim0

data = read_binary_file(dens_filename)
sim0 = np.zeros([len(t), N_theta, N_radial, Nspecies], dtype=DTYPE)
sim0 = data.reshape(sim0.shape)
dens_kxky = sim0

data = read_binary_file(apar_filename)
sim0 = np.zeros([len(t), N_theta, N_radial], dtype=DTYPE)
sim0 = data.reshape(sim0.shape)
apar_kxky = sim0

data = read_binary_file(bpar_filename)
sim0 = np.zeros([len(t), N_theta, N_radial], dtype=DTYPE)
sim0 = data.reshape(sim0.shape)
bpar_kxky = sim0

data = read_binary_file(v_filename)
sim0 = np.zeros([len(t), N_theta, N_radial, Nspecies], dtype=DTYPE)
sim0 = data.reshape(sim0.shape)
v_kxky = sim0


# SVD rank to perform DMD 
svd_rank = 0
# j to label radial grid 
j=2


# step for DMD   
delt = t[1]-t[0]

#-----------------------------------------------------------------------------------------
# fluctuating potential 

dmd = DMD(svd_rank = svd_rank, exact = True)
dmd.fit(phi_kxky[:,:,j].T)

realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)


fig = plt.figure()
gs = gridspec.GridSpec(1,1)
ax1 = fig.add_subplot(gs[0,0])

ax1.axvline(0., linestyle="dashed", color="k", linewidth=0.5)
ax1.axhline(0., linestyle="dashed", color="k", linewidth=0.5)

ax1.plot(realEigs.real, realEigs.imag, 'o')

ax1.set_xlabel(r"realEigs.real, i.e. $\omega$", fontsize=15)
ax1.set_ylabel(r"realEigs.imag, i.e. $\gamma$", fontsize=15)
ax1.tick_params(labelsize=15)



print ('DMD predictions based on Phi:')
print (realEigs)
#plt.show()


#-----------------------------------------------------------------------------------------
# fluctuating density 1

dmd = DMD(svd_rank = svd_rank, exact = True)
dmd.fit(dens_kxky[:,:,j,0].T)

realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)


ax1.axvline(0., linestyle="dashed", color="k", linewidth=0.5)
ax1.axhline(0., linestyle="dashed", color="k", linewidth=0.5)

ax1.plot(realEigs.real, realEigs.imag, 's', color="tab:green", mfc="none", markersize=9.5)

ax1.set_xlabel(r"realEigs.real, i.e. $\omega$", fontsize=15)
ax1.set_ylabel(r"realEigs.imag, i.e. $\gamma$", fontsize=15)
ax1.tick_params(labelsize=15)



print ('DMD predictions based on density (species 1):')
print (realEigs)
#plt.show()


#-----------------------------------------------------------------------------------------
# fluctuating density 2

dmd = DMD(svd_rank = svd_rank, exact = True)
dmd.fit(dens_kxky[:,:,j,1].T)

realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)


ax1.axvline(0., linestyle="dashed", color="k", linewidth=0.5)
ax1.axhline(0., linestyle="dashed", color="k", linewidth=0.5)

ax1.plot(realEigs.real, realEigs.imag, 's', color="green", mfc="none", markersize=8.5)

ax1.set_xlabel(r"realEigs.real, i.e. $\omega$", fontsize=15)
ax1.set_ylabel(r"realEigs.imag, i.e. $\gamma$", fontsize=15)
ax1.tick_params(labelsize=15)



print ('DMD predictions based on density (species 2):')
print (realEigs)
#plt.show()


#-----------------------------------------------------------------------------------------
# fluctuating vector potential 

dmd = DMD(svd_rank = svd_rank, exact = True)
dmd.fit(apar_kxky[:,:,j].T)

realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)


ax1.axvline(0., linestyle="dashed", color="k", linewidth=0.5)
ax1.axhline(0., linestyle="dashed", color="k", linewidth=0.5)

ax1.plot(realEigs.real, realEigs.imag, 'o', color="red", mfc="none")

ax1.set_xlabel(r"realEigs.real, i.e. $\omega$", fontsize=15)
ax1.set_ylabel(r"realEigs.imag, i.e. $\gamma$", fontsize=15)
ax1.tick_params(labelsize=15)



print ('DMD predictions based on Apar:')
print (realEigs)
#plt.show()


#-----------------------------------------------------------------------------------------
# fluctuating magnetic field  

dmd = DMD(svd_rank = svd_rank, exact = True)
dmd.fit(bpar_kxky[:,:,j].T)

realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)


ax1.axvline(0., linestyle="dashed", color="k", linewidth=0.5)
ax1.axhline(0., linestyle="dashed", color="k", linewidth=0.5)

ax1.plot(realEigs.real, realEigs.imag, 'o', color="deeppink", mfc="none", markersize=7.5)

ax1.set_xlabel(r"realEigs.real, i.e. $\omega$", fontsize=15)
ax1.set_ylabel(r"realEigs.imag, i.e. $\gamma$", fontsize=15)
ax1.tick_params(labelsize=15)



print ('DMD predictions based on Bpar:')
print (realEigs)
plt.show()


#-----------------------------------------------------------------------------------------
# fluctuating v 

dmd = DMD(svd_rank = svd_rank, exact = True)
dmd.fit(v_kxky[:,:,j,1].T)

realEigs = np.log(dmd.eigs)/(-complex(0,1)*delt)


ax1.axvline(0., linestyle="dashed", color="k", linewidth=0.5)
ax1.axhline(0., linestyle="dashed", color="k", linewidth=0.5)

ax1.plot(realEigs.real, realEigs.imag, '^', color="cyan", mfc="none")

ax1.set_xlabel(r"realEigs.real, i.e. $\omega$", fontsize=15)
ax1.set_ylabel(r"realEigs.imag, i.e. $\gamma$", fontsize=15)
ax1.tick_params(labelsize=15)



print (realEigs)
plt.show()






