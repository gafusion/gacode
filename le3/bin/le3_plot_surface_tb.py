# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import sys
import os
import numpy as np
import matplotlib.cm as cm
from gacodeplotdefs import *
from scipy.interpolate import griddata

simdir  = sys.argv[1]
imgfile = sys.argv[2]
index   = int(sys.argv[3])

rc('lines',linewidth=1.5)

#-----------------------------------------------
# Spectral dimensions
dim = np.loadtxt(simdir+'/out.le3.geoscalar')
nts = int(dim[0])
nps = int(dim[1])
ns  = int(dim[2])

dx = 2*np.ones(nts+1)
dy = 2*np.ones(nps+1)
dx[0] = 1
dy[0] = 1

z = np.zeros(ns)
#------------------------------------------------

#------------------------------------------------
# Grid variables
nx = 128
ny = 128

t = 2*np.pi*np.arange(nx)/float(nx-1)
p = 2*np.pi*np.arange(ny)/float(ny-1)

f  = np.zeros([nx,ny])
fb = np.zeros([nx,ny])
tb = np.zeros([nx,ny])

vec = np.loadtxt(simdir+'/out.le3.geovector')

# Theta_bar
z=vec[:,0]

i=0
for ips in range(nps+1):
    for its in range(nts+1):
        if its > 0:
            # A
            f = f+np.outer(np.sin(its*t),np.cos(ips*p))*z[i]*dx[its]*dy[ips]
            i = i+1
        if ips > 0 and its > 0:
            # B
            f = f+np.outer(np.sin(its*t),np.sin(ips*p))*z[i]*dx[its]*dy[ips]
            i = i+1
        if ips+its > 0:
            # C
            f = f+np.outer(np.cos(its*t),np.cos(ips*p))*z[i]*dx[its]*dy[ips]
            i = i+1
        if ips > 0:
            # D
            f = f+np.outer(np.cos(its*t),np.sin(ips*p))*z[i]*dx[its]*dy[ips]
            i = i+1
        if ips+its == 0:
            f = f+z[i]
            i = i+1

tb = f

if index < 7:
    z = vec[:,index]
    if index == 0:
        symbol = '\\bar{\\theta}-\\theta'
    if index == 1:
        symbol = 'vx'
    if index == 2:
        symbol = 'flux'
    if index == 3:
        symbol = 'B'
    if index == 4:
        symbol = '1'
    if index == 5:
        symbol = 'Ave'
    if index == 6:
        symbol = 'B'
else:
    vec = np.loadtxt(simdir+'/out.le3.rho')
    z = vec[:,index-7]
    if index == 7:
        symbol = '\\delta\\theta'
    if index == 8:
        symbol = '\\delta\chi'

i=0
for ips in range(nps+1):
    for its in range(nts+1):
        if its > 0:
            # A
            f = f+np.outer(np.sin(its*t),np.cos(ips*p))*z[i]*dx[its]*dy[ips]
            i = i+1
        if ips > 0 and its > 0:
            # B
            f = f+np.outer(np.sin(its*t),np.sin(ips*p))*z[i]*dx[its]*dy[ips]
            i = i+1
        if ips+its > 0:
            # C
            f = f+np.outer(np.cos(its*t),np.cos(ips*p))*z[i]*dx[its]*dy[ips]
            i = i+1
        if ips > 0:
            # D
            f = f+np.outer(np.cos(its*t),np.sin(ips*p))*z[i]*dx[its]*dy[ips]
            i = i+1
        if ips+its == 0:
            f = f+z[i]
            i = i+1

for j in range(ny):
    fb[:,j] = griddata(tb[:,j],f[:,j],t)

# Plotting
fig = plt.figure(figsize=(6,6))
fig.subplots_adjust(left=0.15, right=0.97, top=1.05, bottom=0.0)

ax = fig.add_subplot(111)
ax.set_aspect('equal')
ax.set_title(r'$'+symbol+'$',size=18)

dz=(np.amax(f)-np.amin(f))/256
levels = np.arange(np.amin(f)-dz,np.amax(f)+dz,dz)
ax.contourf(p/np.pi,t/np.pi,fb,levels,cmap=cm.jet,origin='lower')
ax.set_ylabel(r'$\bar{\theta}/\pi$')
ax.set_xlabel(r'$\varphi/\pi$')
ax.set_xlim([0,2])
ax.set_ylim([0,2])

if imgfile == 'screen':
    plt.show()
else:
    plt.savefig(imgfile)
    print("INFO: (le3_plot_surface) Wrote plot to "+imgfile+".")
