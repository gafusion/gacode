#!/usr/bin/env python
# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range

# USAGE: 
#
#    python contour.py <file> <nx> <ny>

import sys
#import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from gacodeplotdefs import *
from pylab import *
from scipy.interpolate import griddata
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

from matplotlib import rc
from matplotlib import axes

#os.environ['PATH'] = os.environ['PATH'] + ':/usr/bin/latex'
ioff()
#matplotlib.use('Agg')

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 30}

rc('font', **font)

rc('text',usetex=True)

mpl.rcParams['axes.linewidth'] = 4
mpl.rcParams['xtick.minor.size'] = 4
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['lines.linewidth'] = 4
#mpl.rcParams['xtick.major.width'] = 4
#mpl.rcParams['xtick.major.width'] = 4

#----------------------------------------------------------------------

nx=64
ny=64

file=sys.argv[1]
nx=int(sys.argv[2])
ny=int(sys.argv[3])

data = np.loadtxt(file)

xm = data[:,0]
ym = data[:,1]
zm1 = data[:,2]
zm2 = data[:,3]

xmin = min(xm)#-0.001
xmax = max(xm)#+0.001

ymin = min(ym)
ymax = max(ym)
zmin1 = min(zm1)
zmax1 = max(zm1)
zmin2 = min(zm2)
zmax2 = max(zm2)

#zmin1 = -(round(max(fabs(zmin1), fabs(zmax1))+0.05, 1))
zmax1 = round(max(fabs(zmin1), fabs(zmax1))+0.05, 1)
zmin1 = -zmax1

#zmin1 = -0.5
#zmax1 = 0.5

zmin2 = 1.0e-2
zmax2 = round(zmax2+0.05, 1)

grid_x, grid_y = np.mgrid[xmin:xmax:nx*1j,ymin:ymax:ny*1j]

grid1 = griddata((xm,ym),zm1,(grid_x, grid_y),fill_value=100)
grid2 = griddata((xm,ym),zm2,(grid_x, grid_y),fill_value=100)

#------------------------------------------------------
# omega
fig = plt.figure(figsize=(13.8,10))
ax = fig.add_subplot(111)

#Set number of levels in contour plot
levels = np.arange(zmin1,zmax1,(zmax1-zmin1)/100)
#Make contour plot
#The standard map is cmap=plt.cm.jet 
plot = ax.contourf(grid_x,grid_y,grid1,levels, cmap=get_cmap('RdBu'))#, extend='both')#cmap=plt.cm.get_cmap('coolwarm'))

#Set contour lines
#contour_ticks = np.logspace(floor(log10(1.0e-1)), ceil(log10(zmax1)),num=(1 + ceil(log10(zmax1)) - floor(log10(1.0e-1)))) 
#contour_ticks = append(contour_ticks,(-1)*contour_ticks)
#contour_ticks = append(contour_ticks,5*contour_ticks)
contour_ticks = np.linspace(round(zmin1,1),round(zmax1,1),9)
#Add contour lines with labels
plot2 = ax.contour(grid_x,grid_y,grid1, levels=contour_ticks, colors = 'k', hold='on', linestyles='dashed', linewidths=3)
ax.clabel(plot2, fmt = '%.1e', colors = 'k', fontsize=14)

#Set ticks in color bar
bar_ticks = np.linspace(zmin1,zmax1,10)
#Add color bar
bar = plt.colorbar(plot, ticks=bar_ticks, format='%.1e', fraction=0.20)
bar.set_label(r'$\omega~\left[c_s/a\right]$', fontsize=40, rotation=270)

#Set labels
ax.set_xlabel(r'$r/a$', fontsize=40)
ax.set_ylabel(r'$k_\theta \rho_s$', fontsize=40)

ax.set_xticklabels(ax.get_xticks(), font)
ax.set_yticklabels(ax.get_yticks(), font)
#ax.tick_param(width=5)
#ax.set_xtickparam(width=5)
#ax.yaxis.Axis.set_tick_params(width=5)


#Plot limits
xlim_min = xmin
xlim_max = xmax
ylim_min = ymin
ylim_max = ymax
#margin = 0.005
#xlim_min = round(xmin-margin, 2)
#xlim_max = round(xmax+margin, 2)
#ylim_min = round(ymin-margin, 2)
#ylim_max = round(ymax+margin, 2)

ax.set_xlim([xlim_min,xlim_max])
ax.set_ylim([ylim_min,ylim_max])

plt.title('Mode frequency')
plt.savefig('omega.pdf')

#------------------------------------------------------
# gamma
fig = plt.figure(figsize=(13.8,10))
ax = fig.add_subplot(111)

#Set number of levels in contour plot
#levels = np.arange(zmin2,zmax2,(zmax2-zmin2)/100)
levels = np.logspace(log10(zmin2),log10(zmax2),num=50) #logarithmic

#Make contour plot
#The standard map is cmap=plt.cm.jet 
plot = ax.contourf(grid_x,grid_y,grid2,levels, cmap=get_cmap('Reds'), norm=LogNorm())#cmap=plt.cm.jet)

#Set contour lines
contour_ticks = np.logspace(floor(log10(zmin2)),ceil(log10(zmax2)),num=(1 + ceil(log10(zmax2)) - floor(log10(zmin2)))) 
contour_ticks = append(contour_ticks,append(2*contour_ticks,5*contour_ticks))
#Add contour lines with labels
plot2 = ax.contour(grid_x,grid_y,grid2, levels=contour_ticks, norm=LogNorm() , colors = 'k', hold='on', linestyles='dashed', linewidths=3)
ax.clabel(plot2, fmt = '%.1e', colors = 'k', fontsize=14)

#Set ticks in color bar
bar_ticks = np.logspace(log10(zmin2),log10(zmax2),num=2*(1 + ceil(log10(zmax2)) - floor(log10(zmin2)))) 
#Add color bar
bar = plt.colorbar(plot, ticks=bar_ticks, format='%.1e', fraction=0.20)
bar.set_label(r'$\gamma~\left[c_s/a\right]$', fontsize=40, rotation=270)

#Set labels
ax.set_xlabel(r'$r/a$', fontsize=40)
ax.set_ylabel(r'$k_\theta \rho_s$', fontsize=40)

#Plot limits
ax.set_xlim([xmin,xmax])
ax.set_ylim([ymin,ymax])

plt.title('Growth rate')

plt.savefig('gamma.pdf')
