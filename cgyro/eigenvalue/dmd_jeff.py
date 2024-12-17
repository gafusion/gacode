import os,sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from pydmd import DMD
from cgyro.data import cgyrodata
from dmd_util import *

#---------------------------------------------------------------------------
# INPUTS
mydir = sys.argv[1]+'/'
nmode = 16

tol = 0.01

# time downsample
k = int(sys.argv[2])

# max time
tmax = float(sys.argv[3])

# theta downsample
l = 1

# SVD rank to perform DMD 
svd_rank = 0

# observables
ostr = ['phi','apar','bpar']

#---------------------------------------------------------------------------
# COLLECT DATA
sim = cgyrodata(mydir,silent=True)
sim.getbigfield()

t = sim.t
n_radial  = sim.n_radial  
n_theta   = sim.theta_plot
n_species = sim.n_species
#n_time = len(t)

# gacode utility for time indices
imin = 0
imax = int(len(t)*tmax)

# balloonings-space potentials
ovec = {}
y = sim.kxky_phi[0,:,:,0,imin:imax]+1j*sim.kxky_phi[1,:,:,0,imin:imax]
ovec['phi'],a0 = map1d(y,sim.q)
y = sim.kxky_apar[0,:,:,0,imin:imax]+1j*sim.kxky_apar[1,:,:,0,imin:imax]
ovec['apar'],a0 = map1d(y,sim.q)
y = sim.kxky_bpar[0,:,:,0,imin:imax]+1j*sim.kxky_bpar[1,:,:,0,imin:imax]
ovec['bpar'],a0 = map1d(y,sim.q)

# step for DMD
dt = k*(t[1]-t[0])

#---------------------------------------------------------------------------
# RUN DMD
dmd = DMD(svd_rank=0,exact=True,sorted_eigs='abs')

evec = {}
modes = {}
for x in ostr:
    down = downsample(ovec[x],k) 
    dmd.fit(down)
    evec[x] = 1j*np.log(dmd.eigs)/dt
    modes[x] = dmd.modes

# determine most unstable modes (zvec)

zvec = np.zeros([3,nmode],dtype=complex)
mvec = np.zeros([3,n_radial*n_theta,nmode],dtype=complex)
for i,x in enumerate(ostr):
    z  = evec[x] 
    zi = z.imag
    k = np.flip(np.argsort(zi))
    zvec[i,:] = z[k[:nmode]]
    mvec[i,:,:] = modes[x][:,k[:nmode]]

# find true eigenmodes and compute errors

m01 = []
m02 = []
for i in range(nmode):
    for j in range(nmode):
        em = abs(zvec[0,i]-zvec[1,j])
        en = abs(zvec[0,i]-zvec[2,j])
        if em < tol:
            m01.append([i,j])
        if en < tol:
            m02.append([i,j])

#---------------------------------------------------------------------------
# PLOTTING

rc('font',size=25)
rc('text',usetex=True)

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)

ax.set_xlabel(r"$\omega$")
ax.set_ylabel(r"$\gamma$")
ax.grid(which="both",ls=":")
ax.grid(which="major",ls=":")

# symbols
args = {}
args['phi']  = {'color':'r','marker':'o'}
args['apar'] = {'color':'b','marker':'s','facecolors':'none'}
args['bpar'] = {'color':'k','marker':'+'}

#for x in ostr:
#    ax.scatter(evec[x].real,evec[x].imag,s=60,**args[x])
for x in ostr:
    ax.scatter(evec[x].real,evec[x].imag,s=60,**args[x],alpha=0.2)

print('-----------------------------------------------')
print('tmax = {}'.format(t[imax-1]))
print('dt   = {}'.format(dt))

freqs = []
field = []
for v in m01:
    for w in m02:
        if v[0] == w[0]:
            i = v[0] ; j = v[1] ; k = w[1]
            e0 = (zvec[0,i]+zvec[1,j]+zvec[2,k])/3 
            err = abs(zvec[0,i]-e0)+abs(zvec[1,j]-e0)+abs(zvec[2,k]-e0)
            print("gamma = {:.3f} omega = {:+.3f} | err = {:.3e}".format(e0.imag,e0.real,err))
            ax.scatter(e0.real,e0.imag,s=60)
            freqs.append(e0)
            field.append(mvec[:,:,i]) 

n = len(freqs)
# theta_* = 0
p0 = n_radial//2*n_theta+n_theta//2
a0 = np.zeros(n,dtype=complex)

#ax.set_xlim([-1.2,1.2])
#ax.set_ylim([-0.2,0.5])

plt.tight_layout()
plt.show()

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)
for i in range(n):
    a0[i] = field[i][0,p0]
    ax.plot(np.real(field[i][0,:]/a0[i]))
plt.show()

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)
for i in range(n):
    ax.plot(np.real(field[i][1,:]/a0[i]))
plt.show()

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111)
for i in range(n):
    ax.plot(np.real(field[i][2,:]/a0[i]))
plt.show()
