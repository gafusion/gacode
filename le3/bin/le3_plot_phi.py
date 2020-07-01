# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import sys
import os
import numpy as np
from gacodeplotdefs import *

simdir  = sys.argv[1]
imgfile = sys.argv[2]

rc('lines',linewidth=2)
rc('mathtext',default='bf')

t  = np.loadtxt(simdir+'/out.le3.t')
p  = np.loadtxt(simdir+'/out.le3.p')
tb = np.loadtxt(simdir+'/out.le3.tb')

#======================================
fig = plt.figure(figsize=(7,6))
fig.subplots_adjust(left=0.17, right=0.97, top=0.95, bottom=0.12)
ax = fig.add_subplot(111)
ax.grid(which="both",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'\boldmath{$\varphi$}')
ax.set_ylabel(r'\boldmath{$\bar\theta$}',color='k')
#=====================================

mt = len(t)
mp = len(p)

for j in range(mt):
    ax.plot(p[:]/np.pi,tb[j,:]/np.pi,'-')

TICKS=[0,1,2]
LABELS=[r'\boldmath{$0$}',r'\boldmath{$\pi$}',r'\boldmath{$2\pi$}']
ax.set_xticks(TICKS)
ax.set_xticklabels(LABELS)
ax.set_xlim(0,2)

TICKS=[0,1,2]
LABELS=[r'\boldmath{$0$}',r'\boldmath{$\pi$}',r'\boldmath{$2\pi$}']
ax.set_yticks(TICKS)
ax.set_yticklabels(LABELS)
ax.set_ylim(0,2)

if imgfile == 'screen':
    plt.show()
else:
    plt.savefig(imgfile)
    print("INFO: (le3_plot_phi) Wrote plot to "+imgfile+".")
