#-------------------------------------------------------------
# plot_template.py
#
# PURPOSE:
#  Simple template for plotting profiles
#-------------------------------------------------------------

import os
import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from gacodefuncs import *
from pygacode import expro
from pygacode import gapystr_get as gapystr

matplotlib.rc('text',usetex=True)
matplotlib.rc('font',size=18)

rmin = 0.1
rmax = 0.98
ext = 'screen'
loc = 2
dot = 1
therm = 1

m1 = 0 ; m2 = 0

alls = not bool(therm)

#=====================================
fig = plt.figure(figsize=(6.5,5))
ax = fig.add_subplot(111)
ax.grid(which="both",ls=":",alpha=0.3)
ax.grid(which="major",ls=":",alpha=0.3)
ax.set_xlabel(r'$r/a$')
#=====================================

def plotit(ax,x,y,ystr):

   global m1,m2,dot
   
   ax.plot(x[m1:m2],y[m1:m2],label=r'$'+ystr+'$')
   if dot:
      ax.plot(x[m1:m2],y[m1:m2],'o',color='k',alpha=0.3,ms=4)     
   return

expro.expro_read('./input.gacode')

x = expro.expro_rmin ; x = x/max(x)
n = expro.expro_n_ion
   
# normalization
a = expro.expro_rmin[-1]
csa = expro.expro_cs/a

sname = gapystr(expro.expro_name) 
stype = gapystr(expro.expro_type)
   
# Set x-range
m1 = 0 ; m2 = len(x)
if rmax != 'auto':
   m2 = np.argmin(np.abs(x-np.float(rmax)))
if rmin != 'auto':
   m1 = np.argmin(np.abs(x-np.float(rmin)))
      
ax.set_xlim([x[m1],x[m2-1]])
ax.set_xlabel(r'$r/a$')

m=m2

y = a*expro.expro_dlntedr ; ystr = 'a/L_{Te}' ; plotit(ax,x,y,ystr)
for p in range(n):
   if alls or stype[p] == '[therm]':
      y = a*expro.expro_dlntidr[p,:]
      ystr = 'a/L_{Ti}~['+sname[p]+']'
      plotit(ax,x,y,ystr)

plt.show()
