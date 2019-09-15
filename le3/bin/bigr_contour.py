import sys
import os
import numpy as np
import matplotlib.cm as cm
from gacodeplotdefs import *

infile = 'out.bigr'
ftype  = 'screen'

data=np.loadtxt(infile)

nphi=32
ntheta=201

x=np.arange(nphi)*1.0/(nphi-2)
y=np.arange(ntheta)*1.0/(ntheta-1)

f=data[:]
print(len(f))

z=f.reshape((ntheta,nphi),order='F')
print(z)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
#ax.set_xlabel(r'$n$')
#ax.set_ylabel(r'$\delta$')

# Plotting
#fig = plt.figure(figsize=(6,6))
#fig.subplots_adjust(left=0.15, right=0.97, top=1.05, bottom=0.0)
#fig = plt.figure(figsize=(6*0.7,6*0.7))
#fig.subplots_adjust(left=0.16, right=0.97, top=1.05, bottom=0.0)

ax.set_aspect('equal')

#dz=(np.amax(f)-np.amin(f))/256
#levels = np.arange(np.amin(f)-dz,np.amax(f)+dz,dz)
ax.contourf(x,y,z)
#ax.set_ylabel(r'$\theta/\pi$')
#ax.set_xlabel(r'$\varphi/\pi$')
#ax.set_xlim([0,2])
#ax.set_ylim([0,2])
plt.show()

