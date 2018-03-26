import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from gyro.data_plot import gyrodata_plot

rc('text',usetex=True)
rc('font',size=18)

plot_type = sys.argv[1]

if plot_type == 'freq':

   w      = float(sys.argv[2])
   ftype  = sys.argv[3]

   gyrodata_plot('./').plot_freq(w=w)

   outfile = 'freq.'+ftype

elif plot_type == 'balloon':

   simdir = sys.argv[2]
   index  = sys.argv[3]
   ftype  = sys.argv[4]
   tmax   = float(sys.argv[5])

   key = gyrodata_plot(simdir).plot_balloon(index=index,tmax=tmax)

   outfile = key+'.'+ftype

elif plot_type == 'zf':

   simdir = sys.argv[2]
   w      = float(sys.argv[3])
   ftype  = sys.argv[4]
   
   gyrodata_plot(simdir).plot_zf(w=w)

   outfile = 'zf.'+ftype

elif plot_type == 'gbflux':

   simdir   = sys.argv[2]
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   w        = float(sys.argv[5])
   lx       = float(sys.argv[6])
   ly       = float(sys.argv[7])
   title    = sys.argv[8]
   ymin     = sys.argv[9]
   ymax     = sys.argv[10]
   span1   = float(sys.argv[11])
   span2   = float(sys.argv[12])
   ftype    = sys.argv[13]
   
   gyrodata_plot(simdir).plot_gbflux(field=field,i_moment=i_moment,w=w,
                                     lx=lx,ly=ly,title=title,ymin=ymin,ymax=ymax,
                                     span1=span1)

   outfile = 'gbflux.'+ftype

elif plot_type == 'gbflux_i':

   simdir   = sys.argv[2]
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   w        = float(sys.argv[5])
   ymin     = sys.argv[6]
   ymax     = sys.argv[7]
   ftype    = sys.argv[8]
   
   gyrodata_plot(simdir).plot_gbflux_i(field=field,i_moment=i_moment,w=w,
                                     ymin=ymin,ymax=ymax)

   outfile = 'gbflux_i.'+ftype

elif plot_type == 'gbflux_n':

   simdir   = sys.argv[2]
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   w        = float(sys.argv[5])
   datafile = sys.argv[6]
   ftype    = sys.argv[7]
   
   gyrodata_plot(simdir).plot_gbflux_n(field=field,i_moment=i_moment,w=w,datafile=datafile)

   outfile = 'gbflux_n.'+ftype

elif plot_type == 'gbflux_rt':

   simdir   = sys.argv[2]
   field    = sys.argv[3]
   i_moment = int(sys.argv[4])
   w        = float(sys.argv[5])
   ftype    = sys.argv[6]
   
   gyrodata_plot(simdir).plot_gbflux_rt(field=field,i_moment=i_moment,w=w)

   outfile = 'gbflux_rt.'+ftype

elif plot_type == 'gbflux_exc':

   simdir   = sys.argv[2]
   w        = float(sys.argv[3])
   ftype    = sys.argv[4]
   
   gyrodata_plot(simdir).plot_gbflux_exc(w=w)

   outfile = 'gbflux_exc.'+ftype

elif plot_type == 'phi_n0':

   simdir = sys.argv[2]
   lx     = float(sys.argv[3])
   ly     = float(sys.argv[4])
   ymax   = sys.argv[5]
   span1  = float(sys.argv[6])
   span2  = float(sys.argv[7])
   ftype  = sys.argv[8]

   ax = gyrodata_plot(simdir).plot_phi_n0(lx=lx,ly=ly,ymax=ymax,span1=span1,span2=span2)
   
   outfile = 'phi_n0.'+ftype
  
elif plot_type == 'moment_zero':

   simdir   = sys.argv[2]
   i_moment = int(sys.argv[3])
   w        = float(sys.argv[4])
   ftype    = sys.argv[5]

   ax = gyrodata_plot(simdir).plot_moment_zero(i_moment=i_moment,w=w)
   
   outfile = 'moment_zero.'+ftype
  

#---------------------------------------------------------------
# Plot to screen or to image file
if ftype == 'screen':
    plt.show()
else:
    plt.savefig(outfile)
#---------------------------------------------------------------
