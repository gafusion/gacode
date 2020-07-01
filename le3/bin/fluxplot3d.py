# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import sys
import numpy as np
from gacodeplotdefs import *

# Number of theta-points for plotting
narc = 128
ntarc = 11

infile = 'out.m3d.geo'
ftype  = 'screen'

data=np.loadtxt(infile)
n_t = int(data[0])+1
n_p = int(data[1])+1
n = len(data)-2
vec = data[2:len(data)]

print('Sanity check: ',n-n_t*n_p*4*4)

vec = np.reshape(vec,(4,4,n_p,n_t),order='F')
ar = vec[0,0,:,:]
br = vec[0,1,:,:]
cr = vec[0,2,:,:]
dr = vec[0,3,:,:]

az = vec[1,0,:,:]
bz = vec[1,1,:,:]
cz = vec[1,2,:,:]
dz = vec[1,3,:,:]

fig = plt.figure(figsize=(8,12))
ax = fig.add_subplot(111,aspect='equal')
ax.set_xlabel(r'$R$')
ax.set_ylabel(r'$Z$')

t = 2*np.pi*np.arange(0,narc)/float(narc-1)
tt = 2*np.pi*np.arange(0,ntarc)/float(ntarc-1)

# Fourier geometry flux-surfaces 

for p in range(ntarc):
    x = np.zeros([narc])
    y = np.zeros([narc])
    for m in range(n_t):
        for n in range(n_p):
            x = x+\
                ar[n,m]*np.sin(m*t)*np.cos(n*tt[p])+\
                br[n,m]*np.sin(m*t)*np.sin(n*tt[p])+\
                cr[n,m]*np.cos(m*t)*np.cos(n*tt[p])+\
                dr[n,m]*np.cos(m*t)*np.sin(n*tt[p]) 
            y = y+\
                az[n,m]*np.sin(m*t)*np.cos(n*tt[p])+\
                bz[n,m]*np.sin(m*t)*np.sin(n*tt[p])+\
                cz[n,m]*np.cos(m*t)*np.cos(n*tt[p])+\
                dz[n,m]*np.cos(m*t)*np.sin(n*tt[p]) 

    ax.plot(x,y,linewidth=1)

for i in range(narc):
    print(t[i],x[i],y[i])

ax.legend()
if ftype == 'screen':
    plt.show()
else:
    outfile = ftype
    plt.savefig(outfile)
