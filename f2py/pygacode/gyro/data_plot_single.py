# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from .data_plot import gyrodata_plot

data_in = gyrodata_plot('./')

# Use first 3 args to define plot and font size 
rc('text',usetex=True)
rc('font',size=int(sys.argv[1]))
data_in.lx = int(sys.argv[2])
data_in.ly = int(sys.argv[3])

# Shift list by 3
sys.argv = sys.argv[3:]

plot_type = sys.argv[1]

if plot_type == 'freq':

   w   = float(sys.argv[2])
   ext = sys.argv[3]
   wmax = 0.0

   data_in.plot_freq(w=w,wmax=wmax)

   outfile = 'freq.'+ext

elif plot_type == 'ball':

   index = sys.argv[2]
   ext   = sys.argv[3]
   tmax  = float(sys.argv[4])

   key = data_in.plot_balloon(index=index,tmax=tmax)

   outfile = key+'.'+ext

elif plot_type == 'zf':
   
   w   = float(sys.argv[2])
   ext = sys.argv[3]
   
   data_in.plot_zf(w=w)

   outfile = 'zf.'+ext

elif plot_type == 'phi':

   ymax = sys.argv[2]
   ext  = sys.argv[3]

   ax = data_in.plot_phi(ymax=ymax)
   
   outfile = 'phi.'+ext

elif plot_type == 'flux':

   w      = float(sys.argv[2])
   field  = sys.argv[3]
   moment = int(sys.argv[4])
   ymin   = sys.argv[5]
   ymax   = sys.argv[6]
   ext    = sys.argv[7]
   
   data_in.plot_gbflux(w=w,field=field,moment=moment,ymin=ymin,ymax=ymax)

   outfile = 'flux.'+ext

elif plot_type == 'flux_i':

   w      = float(sys.argv[2])
   aw     = int(sys.argv[3])
   field  = sys.argv[4]
   moment = int(sys.argv[5])
   ymin   = sys.argv[6]
   ymax   = sys.argv[7]
   ext    = sys.argv[8]
   
   data_in.plot_gbflux_i(w=w,aw=aw,field=field,moment=moment,ymin=ymin,ymax=ymax)

   outfile = 'flux_i.'+ext

elif plot_type == 'flux_n':

   w      = float(sys.argv[2])
   field  = sys.argv[3]
   moment = int(sys.argv[4])
   ymin   = sys.argv[5]
   ymax   = sys.argv[6]
   ext    = sys.argv[7]
   
   data_in.plot_gbflux_n(w=w,field=field,moment=moment,ymin=ymin,ymax=ymax)

   outfile = 'flux_n.'+ext

elif plot_type == 'flux_exc':

   w      = float(sys.argv[2])
   ext    = sys.argv[3]
   
   data_in.plot_gbflux_exc(w=w)

   outfile = 'flux_exc.'+ext

elif plot_type == 'flux_rt':

   w      = float(sys.argv[2])
   field  = sys.argv[3]
   moment = int(sys.argv[4])
   ext    = sys.argv[5]
   
   data_in.plot_gbflux_rt(w=w,field=field,moment=moment)

   outfile = 'flux_rt.'+ext
  
elif plot_type == 'hmom':

   w       = float(sys.argv[2])
   species = int(sys.argv[3])
   moment  = int(sys.argv[4])
   ext     = sys.argv[5]
   wmax = 0.0

   ax = data_in.plot_moment_zero(w=w,wmax=wmax,species=species,moment=moment)
   
   outfile = 'moment_zero.'+ext
  
elif plot_type == 'dprof':

   w       = float(sys.argv[2])
   species = int(sys.argv[3])
   moment  = int(sys.argv[4])
   ymin    = sys.argv[5]
   ymax    = sys.argv[6]
   ext     = sys.argv[7]

   ax = data_in.plot_profile_tot(w=w,species=species,moment=moment,ymin=ymin,ymax=ymax)
   
   outfile = 'profile_tot.'+ext
  
#---------------------------------------------------------------
# Plot to screen or to image file
if ext == 'screen':
    plt.show()
else:
    plt.savefig(outfile)

if ext != 'screen':
       print('INFO: (data_plot_single) Created '+outfile)
#---------------------------------------------------------------
