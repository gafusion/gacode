import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from cgyro.data_plot import cgyrodata_plot
from cgyro.data_dump import cgyrodata_dump

data_in = cgyrodata_plot('./')

# Use first 3 args to define plot and font size 
rc('text',usetex=True)
rc('font',size=int(sys.argv[1]))
data_in.lx = int(sys.argv[2])
data_in.ly = int(sys.argv[3])

# Shift list by 3
sys.argv = sys.argv[3:]

plot_type = sys.argv[1]

if plot_type == 'freq':
   
   w     = float(sys.argv[2])
   ftype = sys.argv[3]

   data_in.plot_freq(w=w)

   outfile = 'out.cgyro.freq.'+ftype

elif plot_type == 'ky_freq':
   
   w     = float(sys.argv[2])
   ftype = sys.argv[3]

   data_in.plot_ky_freq(w=w)

   outfile = 'out.cgyro.ky_freq.'+ftype

elif plot_type == 'ky_phi':
   
   field = int(sys.argv[2])
   ymin  = sys.argv[3]
   ymax  = sys.argv[4]
   nstr  = sys.argv[5]
   ftype = sys.argv[6]

   data_in.plot_ky_phi(field=field,ymin=ymin,ymax=ymax,nstr=nstr)

   outfile = 'out.cgyro.ky_phi.'+ftype

elif plot_type == 'rcorr_phi':
   
   field = int(sys.argv[2])
   w     = float(sys.argv[3])
   ftype = sys.argv[4]

   data_in.plot_rcorr_phi(field=field,w=w)

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
   field = int(sys.argv[3])
   ftype = sys.argv[4]

   data_in.plot_zf(w=w,field=field)
   
   outfile = 'out.cgyro.zf.'+ftype

elif plot_type == 'phi':

   field = int(sys.argv[2])
   ftype = sys.argv[3]

   head,x,y1,y2 = data_in.plot_phi(field=field)
   
   outfile = 'out.cgyro.phi.'+ftype

elif plot_type == 'flux':

   w      = float(sys.argv[2])
   field  = int(sys.argv[3])
   moment = sys.argv[4]
   ymin   = sys.argv[5]
   ymax   = sys.argv[6]
   fc     = int(sys.argv[7])
   ftype  = sys.argv[8]
   loc    = int(sys.argv[9])
   nscale = int(sys.argv[10])

   if ftype == 'dump':
      cgyrodata_dump('./').dump_flux(fc=fc)
      sys.exit()
   else:
      data_in.plot_flux(w=w,field=field,moment=moment,ymin=ymin,ymax=ymax,fc=fc,loc=loc,nscale=nscale)

   outfile = 'out.cgyro.flux.'+ftype

elif plot_type == 'ky_flux':

   w      = float(sys.argv[2])
   field  = int(sys.argv[3])
   moment = sys.argv[4]
   ymin   = sys.argv[5]
   ymax   = sys.argv[6]
   fc     = int(sys.argv[7])
   ftype  = sys.argv[8]

   if ftype == 'dump':
      cgyrodata_dump('./').dump_ky_flux(w=w,field=field,moment=moment,fc=fc)
      sys.exit()
   else:
      data_in.plot_ky_flux(w=w,field=field,moment=moment,ymin=ymin,ymax=ymax,fc=fc)

   outfile = 'out.cgyro.ky_flux.'+ftype

elif plot_type == 'xflux':

   w      = float(sys.argv[2])
   moment = sys.argv[3]
   ymin   = sys.argv[4]
   ymax   = sys.argv[5]
   ftype  = sys.argv[6]
   nscale = int(sys.argv[7])

   data_in.plot_xflux(w=w,moment=moment,ymin=ymin,ymax=ymax,nscale=nscale)

   outfile = 'out.cgyro.xflux.'+ftype

elif plot_type == 'kxky_phi':

   field = int(sys.argv[2])
   w     = float(sys.argv[3])
   ftype = sys.argv[4]

   data_in.plot_kxky_phi(field=field,w=w)
   
   outfile = 'out.cgyro.kxky_phi.'+ftype

elif plot_type == 'kx_phi':

   field = int(sys.argv[2])
   w     = float(sys.argv[3])
   ymin  = sys.argv[4]
   ymax  = sys.argv[5]
   nstr  = sys.argv[6]
   ftype = sys.argv[7]

   data_in.plot_kx_phi(field=field,w=w,ymin=ymin,ymax=ymax,nstr=nstr)
   
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
   theta = sys.argv[5]
   ftype = sys.argv[6]

   data_in.plot_hbcut(itime=itime,spec=spec,tmax=tmax,theta=theta)
   
   outfile = 'out.cgyro.hbcut.'+ftype


else:

   print 'ERROR: (data_plot_single) Plot type not found'
   
#---------------------------------------------------------------
# Plot to screen or to image file

if ftype == 'screen':
    plt.show()
elif ftype == 'dump':
    data = np.column_stack((x,y1,y2))
    np.savetxt(outfile,data,fmt='%.8e',header=head)
else:
    plt.savefig(outfile)

if ftype != 'screen':
       print 'INFO: (data_plot_single) Created '+outfile
#---------------------------------------------------------------

