import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from .data_plot import cgyrodata_plot
from .data_dump import cgyrodata_dump

data_in = cgyrodata_plot('./')

# Use first 4 args to define plot and font size

# Use latex fonts if set
if int(sys.argv[1]) == 0:
   rc('text',usetex=False)
else:
   rc('text',usetex=True)

rc('font',size=int(sys.argv[2]))
data_in.lx = int(sys.argv[3])
data_in.ly = int(sys.argv[4])

# Shift list by 4
sys.argv = sys.argv[4:]

plot_type = sys.argv[1]

doplot=True

if plot_type == 'freq':

   w     = float(sys.argv[2])
   wmax  = float(sys.argv[3])
   norm  = sys.argv[4]
   ftype = sys.argv[5]

   data_in.plot_freq(w=w,wmax=wmax,norm=norm)

   outfile = 'out.cgyro.freq.'+ftype

elif plot_type == 'ky_freq':

   w     = float(sys.argv[2])
   wmax  = float(sys.argv[3])
   norm  = sys.argv[4]
   ftype = sys.argv[5]

   data_in.plot_ky_freq(w=w,wmax=wmax,norm=norm)

   outfile = 'out.cgyro.ky_freq.'+ftype

elif plot_type == 'ky_phi':

   field = int(sys.argv[2])
   theta = float(sys.argv[3])
   ymin  = sys.argv[4]
   ymax  = sys.argv[5]
   nstr  = sys.argv[6]
   norm  = sys.argv[7]
   ftype = sys.argv[8]

   data_in.plot_ky_phi(field=field,theta=theta,ymin=ymin,ymax=ymax,nstr=nstr,norm=norm)

   outfile = 'out.cgyro.ky_phi.'+ftype

elif plot_type == 'rcorr_phi':

   field = int(sys.argv[2])
   theta = float(sys.argv[3])
   w     = float(sys.argv[4])
   wmax  = float(sys.argv[5])
   ftype = sys.argv[6]

   data_in.plot_rcorr_phi(field=field,theta=theta,w=w,wmax=wmax)

   outfile = 'out.cgyro.rcorr_phi.'+ftype

elif plot_type == 'geo':

   ftype = sys.argv[2]

   data_in.plot_geo()

   outfile = 'out.cgyro.geo.'+ftype

elif plot_type == 'error':

   ftype = sys.argv[2]

   data_in.plot_error()

   outfile = 'out.cgyro.error.'+ftype

elif plot_type == 'ball':

   itime = int(sys.argv[2])
   field = int(sys.argv[3])
   tmax  = float(sys.argv[4])
   ftype = sys.argv[5]

   head,x,y1,y2 = data_in.plot_ball(itime=itime,field=field,tmax=tmax)

   outfile = 'out.cgyro.ball.'+ftype

elif plot_type == 'zf':

   w     = float(sys.argv[2])
   wmax  = float(sys.argv[3])
   field = int(sys.argv[4])
   ftype = sys.argv[5]

   data_in.plot_zf(w=w,wmax=wmax,field=field)

   outfile = 'out.cgyro.zf.'+ftype

elif plot_type == 'phi':

   w     = float(sys.argv[2])
   wmax  = float(sys.argv[3])
   field = int(sys.argv[4])
   theta = float(sys.argv[5])
   ymin  = sys.argv[6]
   ymax  = sys.argv[7]
   rms   = int(sys.argv[8])
   ftype = sys.argv[9]

   head,x,y1,y2 = data_in.plot_phi(w=w,wmax=wmax,field=field,theta=theta,ymin=ymin,ymax=ymax,rms=rms)

   outfile = 'out.cgyro.phi.'+ftype

elif plot_type == 'low':

   w     = float(sys.argv[2])
   wmax  = float(sys.argv[3])
   spec  = int(sys.argv[4])
   moment = sys.argv[5]
   ftype = sys.argv[6]
   theta = 0.0
   ymin  = sys.argv[7]
   ymax  = sys.argv[8]

   data_in.plot_low(w=w,wmax=wmax,spec=spec,moment=moment,theta=theta,ymin=ymin,ymax=ymax)

   outfile = 'out.cgyro.low.'+ftype

elif plot_type == 'corrug':

   w     = float(sys.argv[2])
   wmax  = float(sys.argv[3])
   spec  = int(sys.argv[4])
   moment = sys.argv[5]
   ftype = sys.argv[6]
   theta = 0.0
   ymin  = sys.argv[7]
   ymax  = sys.argv[8]

   data_in.plot_corrug(w=w,wmax=wmax,spec=spec,moment=moment,theta=theta,ymin=ymin,ymax=ymax)

   outfile = 'out.cgyro.corrug.'+ftype

elif plot_type == 'shift':

   w     = float(sys.argv[2])
   wmax  = float(sys.argv[3])
   ftype = sys.argv[4]
   theta = 0.0
   ymin  = sys.argv[5]
   ymax  = sys.argv[6]

   data_in.plot_shift(w=w,wmax=wmax,theta=theta,ymin=ymin,ymax=ymax)

   outfile = 'out.cgyro.shift.'+ftype

elif plot_type == 'flux':

   w      = float(sys.argv[2])
   wmax   = float(sys.argv[3])
   field  = int(sys.argv[4])
   moment = sys.argv[5]
   ymin   = sys.argv[6]
   ymax   = sys.argv[7]
   fc     = int(sys.argv[8])
   ftype  = sys.argv[9]
   loc    = int(sys.argv[10])
   nscale = int(sys.argv[11])
   cflux  = sys.argv[12]
   norm   = sys.argv[13]

   if ftype == 'dump':
      cgyrodata_dump('./').dump_flux(fc=fc)
      doplot=False
   else:
      data_in.plot_flux(w=w,wmax=wmax,field=field,moment=moment,
                        ymin=ymin,ymax=ymax,fc=fc,loc=loc,nscale=nscale,cflux=cflux,norm=norm)

   outfile = 'out.cgyro.flux.'+ftype

elif plot_type == 'ky_flux':

   w      = float(sys.argv[2])
   wmax   = float(sys.argv[3])
   field  = int(sys.argv[4])
   moment = sys.argv[5]
   ymin   = sys.argv[6]
   ymax   = sys.argv[7]
   fc     = int(sys.argv[8])
   ftype  = sys.argv[9]
   diss   = int(sys.argv[10])
   cflux  = sys.argv[11]

   if ftype == 'dump':
      cgyrodata_dump('./').dump_ky_flux(w=w,wmax=wmax,field=field,moment=moment,fc=fc)
      doplot=False
   else:
      data_in.plot_ky_flux(w=w,field=field,moment=moment,
                           ymin=ymin,ymax=ymax,fc=fc,diss=diss,cflux=cflux)

   outfile = 'out.cgyro.ky_flux.'+ftype

elif plot_type == 'xflux':

   w      = float(sys.argv[2])
   wmax   = float(sys.argv[3])
   moment = sys.argv[4]
   ymin   = sys.argv[5]
   ymax   = sys.argv[6]
   ftype  = sys.argv[7]
   nscale = int(sys.argv[8])

   data_in.plot_xflux(w=w,wmax=wmax,moment=moment,ymin=ymin,ymax=ymax,nscale=nscale)

   outfile = 'out.cgyro.xflux.'+ftype

elif plot_type == 'kxky_phi':

   field = int(sys.argv[2])
   theta = float(sys.argv[3])
   w     = float(sys.argv[4])
   wmax  = float(sys.argv[5])
   ftype = sys.argv[6]

   data_in.plot_kxky_phi(field=field,theta=theta,w=w,wmax=wmax)

   outfile = 'out.cgyro.kxky_phi.'+ftype

elif plot_type == 'kx_phi':

   field = int(sys.argv[2])
   theta = float(sys.argv[3])
   w     = float(sys.argv[4])
   wmax  = float(sys.argv[5])
   ymin  = sys.argv[6]
   ymax  = sys.argv[7]
   nstr  = sys.argv[8]
   ftype = sys.argv[9]
   diss = int(sys.argv[10])

   data_in.plot_kx_phi(field=field,theta=theta,w=w,wmax=wmax,ymin=ymin,ymax=ymax,nstr=nstr,diss=diss)

   outfile = 'out.cgyro.kx_phi.'+ftype

elif plot_type == 'hb':

   itime = int(sys.argv[2])
   spec  = int(sys.argv[3])
   tmax  = float(sys.argv[4])
   mesh  = int(sys.argv[5])
   ftype = sys.argv[6]

   data_in.plot_hb(itime=itime,spec=spec,tmax=tmax,mesh=mesh)

   outfile = 'out.cgyro.hb.'+ftype

elif plot_type == 'hbcut':

   itime = int(sys.argv[2])
   spec  = int(sys.argv[3])
   tmax  = float(sys.argv[4])
   theta = float(sys.argv[5])
   ftype = sys.argv[6]

   data_in.plot_hbcut(itime=itime,spec=spec,tmax=tmax,theta=theta)

   outfile = 'out.cgyro.hbcut.'+ftype

elif plot_type == 'hball':

   itime = int(sys.argv[2])
   spec  = int(sys.argv[3])
   ymin  = sys.argv[4]
   ymax  = sys.argv[5]
   tmax  = float(sys.argv[6])
   nstr  = sys.argv[7]
   ie    = int(sys.argv[8])
   ftype = sys.argv[9]

   data_in.plot_hball(itime=itime,spec=spec,ymin=ymin,ymax=ymax,tmax=tmax,nstr=nstr,ie=ie)

   outfile = 'out.cgyro.hball.'+ftype

else:

   print('ERROR: (data_plot_single) Plot type not found')

#---------------------------------------------------------------
# Plot to screen or to image file

if doplot:
   if ftype == 'screen':
      plt.show()
   elif ftype == 'dump':
      data = np.column_stack((x,y1,y2))
      np.savetxt(outfile,data,fmt='%.8e',header=head)
   else:
      plt.savefig(outfile)

   if ftype != 'screen':
      print('INFO: (data_plot_single) Created '+outfile)
#---------------------------------------------------------------

