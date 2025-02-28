import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from .data_plot import cgyrodata_plot
from .data_dump import cgyrodata_dump

data_in = cgyrodata_plot('./')

# Use latex fonts if set
if int(sys.argv[1]) == 0:
   rc('text',usetex=False)
else:
   rc('text',usetex=True)

rc('font',size=int(sys.argv[2]))

plot_type = sys.argv[3]

# Rule out most plots for ZF test
if data_in.zf_test:
   zf_test = ['zf','geo','error']
   if plot_type not in zf_test:
      plot_type = 'zf_not'

xin = data_in.plot_dictinit()

xin['fig']    = None
xin['lx']     = int(sys.argv[4])
xin['ly']     = int(sys.argv[5])
xin['w']      = sys.argv[6]
xin['norm']   = sys.argv[7]
xin['ftype']  = sys.argv[8]
xin['itime']  = int(sys.argv[9])
xin['field']  = int(sys.argv[10])
xin['fnorm']  = int(sys.argv[11])
xin['moment'] = sys.argv[12]
xin['tmax']   = float(sys.argv[13])
xin['theta']  = int(sys.argv[14])
xin['ymin']   = sys.argv[15]
xin['ymax']   = sys.argv[16]
xin['kxmin']  = sys.argv[17]
xin['kxmax']  = sys.argv[18]
xin['nstr']   = sys.argv[19]
xin['abs']    = int(sys.argv[20])
xin['fc']     = int(sys.argv[21])
xin['loc']    = int(sys.argv[22])
xin['nscale'] = int(sys.argv[23])
xin['cflux']  = sys.argv[24]
xin['spec']   = int(sys.argv[25])
xin['bar']    = bool(int(sys.argv[26]))
xin['ie']     = int(sys.argv[27])
xin['mesh']   = int(sys.argv[28])

# plotting control
ftype   = xin['ftype']
outfile = 'out.cgyro.'+plot_type+'.'+ftype
doplot  = True

if plot_type == 'freq':
   
   head = data_in.plot_freq(xin)

elif plot_type == 'ky_freq':
   
   head,x,y1,y2 = data_in.plot_ky_freq(xin)

elif plot_type == 'error':

   head = data_in.plot_error(xin)

elif plot_type == 'geo':

   head = data_in.plot_geo(xin)

elif plot_type == 'ball':

   head,x,y1,y2 = data_in.plot_ball(xin)

elif plot_type == 'ky_phi':

   head = data_in.plot_ky_phi(xin)

elif plot_type == 'phi':

   head,x,y1,y2 = data_in.plot_phi(xin)

elif plot_type == 'triad':

   head = data_in.plot_triad(xin)

elif plot_type == 'triad_v2':

   head = data_in.plot_triad_v2(xin)

elif plot_type == 'triad_v3':

   head = data_in.plot_triad_v3(xin)

elif plot_type == 'flux':

   if ftype == 'nox' or ftype == 'dump':
       doplot = False
  
   if ftype == 'dump':
      cgyrodata_dump('./').dump_flux(xin)
   else:
      data_in.plot_flux(xin)

elif plot_type == 'rcorr_phi':

   head = data_in.plot_rcorr_phi(xin)
                                 
elif plot_type == 'low':

   head = data_in.plot_low(xin)

elif plot_type == 'zf':

   head = data_in.plot_zf(xin)

elif plot_type == 'corrug':

   head = data_in.plot_corrug(xin)

elif plot_type == 'shift':

   head,x,y1,y2 = data_in.plot_shift(xin)

elif plot_type == 'ky_flux':
  
   if ftype == 'nox' or ftype == 'dump':
       doplot = False

   if ftype == 'dump':
      cgyrodata_dump('./').dump_ky_flux(xin)
   else:
      data_in.plot_ky_flux(xin)

elif plot_type == 'xflux':

   head = data_in.plot_xflux(xin)

elif plot_type == 'kxky_phi':

   head = data_in.plot_kxky_phi(xin)

elif plot_type == 'kx_phi':
  
   head = data_in.plot_kx_phi(xin)

elif plot_type == 'kx_shift':
   
   head = data_in.plot_kx_shift(xin)

elif plot_type == 'poly_phi':
  
   head = data_in.plot_poly_phi(xin)

elif plot_type == 'hb':

   head = data_in.plot_hb(xin)

elif plot_type == 'hbcut':

   head = data_in.plot_hbcut(xin)

elif plot_type == 'hball':

   head = data_in.plot_hball(xin)

elif plot_type == 'zf_not':

   print('ERROR: (data_plot_single) Not available with ZF test.')
   
else:

   print('ERROR: (data_plot_single) Plot type not found')

#---------------------------------------------------------------
# Plot to screen or to image file

if doplot:
   if ftype == 'screen':
      plt.show()
   elif ftype == 'dump':
      if head == None:
         print('WARNING: (data_plot_single) dump not available')
      else:
         if y2 is None:
            data = np.column_stack((x,y1))
         else:
            data = np.column_stack((x,y1,y2))
         np.savetxt(outfile,data,fmt='%.8e',header=head)
         print('INFO: (data_plot_single) Created '+outfile)
   else:
      plt.savefig(outfile)
#---------------------------------------------------------------

