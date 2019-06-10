import sys
import os
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
#from gacodeplotdefs import *
from gacodefuncs import *

simdir = sys.argv[1]
imgfile = sys.argv[2]

x = np.loadtxt(simdir+'/out.neo.gxi_x')
y = np.loadtxt(simdir+'/out.neo.gxi_t')/np.pi
z_in = np.loadtxt(simdir+'/out.neo.gxi')

grid = np.loadtxt(simdir+'/out.neo.grid')

nx = len(x)
ny = len(y)
ne = int(grid[1]+1)
ns = int(grid[0])

z_in = np.reshape(z_in,(nx,ny,ne,ns),order='F')

# Plotting
fig = plt.figure(figsize=(20,6))
#fig.subplots_adjust(left=0.15, right=0.97, top=1.05, bottom=0.0)
#fig = plt.figure(figsize=(6*0.7,6*0.7))
#fig.subplots_adjust(left=0.16, right=0.97, top=1.05, bottom=0.0)

for j in range(ns):
    for ie in range(ne):
        ax = fig.add_subplot(ns,ne,ie+j*ne+1)
        ax.set_aspect('equal')
        ax.set_title(ie,size=18)

        z = z_in[:,:,ie,j]
        dz=(np.amax(z)-np.amin(z))/128
        levels = np.arange(np.amin(z)-dz,np.amax(z)+dz,dz)
        ax.contourf(y,x,z,levels,cmap=cm.jet,origin='lower')
        ax.set_ylabel(r'$\theta$')
        ax.set_xlabel(r'$\xi$')
        ax.set_xlim([-1,1])
        ax.set_ylim([-1,1])

if imgfile == 'screen':
    plt.show()
else:
    plt.savefig(imgfile)
    print "INFO: () Wrote plot to "+imgfile+"."
