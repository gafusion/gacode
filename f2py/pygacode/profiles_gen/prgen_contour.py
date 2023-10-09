import sys
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from scipy import interpolate

def efit_rzmap(g,x,y,mx,my):

   # Map from pixels to physical dimensions
   r = g['RLEFT']+g['RDIM']/(mx-1)*x
   z = g['ZMID']-g['ZDIM']/2+g['ZDIM']/(my-1)*y

   return r,z

def efit_rzmapi(g,r,z,mx,my):

   # Map from physical units to (float) pixels
   x = (r-g['RLEFT'])/g['RDIM']*(mx-1)
   y = (z-g['ZMID']+g['ZDIM']/2)/g['ZDIM']*(my-1)

   return x,y   

def prgen_contour(g,mag,nc,psinorm,narc):

   plot=True
   
   #-------------------------------------------------------------------------
   # Extract data from gfile
   #
   # Del* psi = -mu R J_phi = -mu R^2 p' 
   #
   # J_phi = R p' + II'/R
   
   # phi [Wb/rad]
   
   nx    = g['NW']
   ny    = g['NH']
   psi0  = g['SIMAG'] # Poloidal flux at EFIT magnetic axis
   psi1  = g['SIBRY'] # Poloidal flux at EFIT separatrix
   efitp = g['PRES']  # Plasma pressure [Pa]
   efitq = g['QPSI']  # Safety factor
   efitf = g['FPOL']  # Poloidal current function I = R Bt [T-m]
   lx    = g['RDIM']
   ly    = g['ZDIM']

   npsi = len(efitp)
   psi_efit = g['PSIRZ']
   #-------------------------------------------------------------------------

   #-------------------------------------------------------------------------
   # Pixel representation
   xv = np.linspace(0,1,nx) ; yv = np.linspace(0,1,ny)

   if mag == 'auto':
      mag = int(np.rint(512/np.sqrt(nx*ny)))
      if mag < 1:
         mag = 1
   else:
      mag = int(mag)

   print('INFO: (prgen_contour) Contour magnification (mag) = {:d}'.format(mag))
  
   if mag > 1:
      # spline interpolation (quadratic, not cubic, to avoid oscillation)
      spl = interpolate.RectBivariateSpline(yv,xv,psi_efit,kx=2,ky=2)
      mx = mag*nx ; my = mag*ny
      xv = np.linspace(0,1,mx) ; yv = np.linspace(0,1,my)
      psi_efit = spl(yv,xv)
   else:
      mx = nx ; my = ny

   # Pixel values of magnetic axis
   ixa,iya = efit_rzmapi(g,g['RMAXIS'],g['ZMAXIS'],mx,my)
   #-------------------------------------------------------------------------

   #-------------------------------------------------------------------------
   # Find separatrix using bisection
   #
   # Notation:
   #    (ixc,iyc)=contour, (ixa,iya)=axis, (ixs,iys)=separatrix
   dpsi = psi1-psi0
   z0 = psi1-0.9*dpsi
   dz = dpsi/10
   tol = 1e-15
   iycv = []
   
   while abs(dz) > tol: 
      z0 = z0+dz
      frac = (z0-psi0)/(psi1-psi0)
      contours = measure.find_contours(psi_efit,z0)
      found = 0
      for contour in contours:
         # Test for closed contour
         if contour[0,0] == contour[-1,0] and contour[0,1] == contour[-1,1]:
            ixc = contour[:,1] ; iyc = contour[:,0]
            # Test for circling the origin
            if min(ixc) < ixa < max(ixc) and min(iyc) < iya < max(iyc):
               iycv.append(min(iyc))
               if len(iycv) > 3:
                  d1 = iycv[-1]-iycv[-2] ; d2 = iycv[-2]-iycv[-3]
                  # Test for large jump in contour
                  if abs(d1) > 2*abs(d2):
                     del iycv[-1]
                  else:
                     # Success: this is a legitimate flux-surface
                     found = 1
                     ixs = ixc ; iys = iyc

      if found == 0:
         z0 = z0-dz
         dz = dz/2

   psi_sep = z0
   #-------------------------------------------------------------------------

   psi_efit_min = np.min(psi_efit[int(iya)-8:int(iya)+8,int(ixa)-8:int(ixa)+8])
   print('INFO: (prgen_contour) psi_sep (% efit) = {:.3f}'.format(100*(psi_sep-psi0)/(psi1-psi0)))

   #-------------------------------------------------------------------------
   # Compute interior contours
   #
   rv = np.zeros([narc,nc])
   zv = np.zeros([narc,nc])
   psic = np.zeros(nc)

   # Move inside of separatrix by amount psinorm
   dpsi = psinorm*(psi1-psi0)/(nc-1)

   for i in range(nc):
      z0 = psi0+i*dpsi
      psic[i] = z0
      if i == 0:
         # Leave rv=zv=0
         continue
      contours = measure.find_contours(psi_efit,z0)
      for contour in contours:
         if contour[0,0] == contour[-1,0] and contour[0,1] == contour[-1,1]:
            ixc = contour[:,1] ; iyc = contour[:,0]
            if min(ixc) < ixa < max(ixc) and min(iyc) < iya < max(iyc):
               # Arc length
               dl = np.sqrt(np.diff(ixc)**2+np.diff(iyc)**2)
               larc = np.zeros([len(ixc)]) ; larc[1:] = np.cumsum(dl)
               # Cubic interpolation from fine contour mesh to coarse t-mesh
               t = np.linspace(0,1,narc)*larc[-1]
               cs = interpolate.splrep(larc,ixc,per=True) ; r=interpolate.splev(t,cs) 
               cs = interpolate.splrep(larc,iyc,per=True) ; z=interpolate.splev(t,cs)

               # Shift elements so that first index at max(R).
               s = np.argmax(r)
               r[0:-1] = np.roll(r[0:-1],-s) ; r[-1] = r[0] 
               z[0:-1] = np.roll(z[0:-1],-s) ; z[-1] = z[0]

               rv[:,i],zv[:,i] = efit_rzmap(g,r,z,mx,my)

   #-------------------------------------------------------------------------

   #-------------------------------------------------------------------------
   # Radial profiles of p,f,q
   #
   efitpsi = np.linspace(psi0,psi1,npsi)
   cs = interpolate.interp1d(efitpsi,efitp,kind='quadratic') ; out_p = cs(psic)
   cs = interpolate.interp1d(efitpsi,efitf,kind='quadratic') ; out_f = cs(psic)
   cs = interpolate.interp1d(efitpsi,efitq,kind='quadratic') ; out_q = cs(psic)

   # Recalculate q (as a check) based on definition (and some identities)
   loopint = np.zeros(nc-1)
   for k in range(1,narc-1):
      loopint[:] = loopint[:]+(rv[k+1,1:]-rv[k,1:])*(zv[k+1,1:]+zv[k,1:])/(rv[k+1,1:]+rv[k,1:])

   cs,tck,ier,msg = interpolate.splrep(psic[1:],loopint,full_output=True)
   if ier <= 0:
      qi = interpolate.splev(psic[1:],cs,der=1)
      qi = out_f[1:]*qi/(2*np.pi)
   else:
      print('WARNING: (prgen_contour) Problem in splrep')
      qi = 0
            
   #-------------------------------------------------------------------------

   if plot:      
   
      # 1. Flux contours
      
      asp = ly/lx*(nx/ny)
      fig,ax = plt.subplots(figsize=(8*lx/ly,6))
      ax.imshow(g['PSIRZ'],cmap=plt.cm.hsv,aspect=asp,origin='lower')
      ax.set_xlabel('EFIT cell')
      ax.set_ylabel('EFIT cell')

      # Rescaling factor for resampled contours
      xscale = (nx-1)/(mx-1)
      yscale = (ny-1)/(my-1)
      
      # A few contours
      for i in [nc//2,nc-1]:
         x,y = efit_rzmapi(g,rv[:,i],zv[:,i],mx,my)
         label = r'$\psi = {:.4f}$'.format((psic[i]-psi0)/(psi1-psi0))
         ax.plot(x*xscale,y*yscale,'--k',linewidth=1,label=label)
         
      contours = measure.find_contours(psi_efit,psi1)
      for contour in contours:
         ixc = contour[:,1] ; iyc = contour[:,0]
         label = r'$\psi_\mathrm{separatrix}^\mathrm{EFIT}$' 
         ax.plot(ixc*xscale,iyc*yscale,color='b',linewidth=1,label=label)
               
      # GACODE-computed Separatrix
      label = r'$\psi_\mathrm{separatrix}^\mathrm{GACODE}$' 
      ax.plot(ixs*xscale,iys*yscale,'w',linewidth=1,label=label)

      ax.legend()
      plt.tight_layout()
      
      ofile = 'prgen_efit.pdf'
      print('INFO: (prgen_contour) Writing '+ofile)
      plt.savefig(ofile)

      # 2. q-profiles

      fig,ax = plt.subplots(figsize=(8,5))
      x = np.sqrt((psic-psi0)/(psi1-psi0))
      if ier <= 0:
         ax.plot(x[1:],np.sqrt(abs(qi)),color='magenta',linewidth=1,
                 label=r'$\mathbf{GACODE~integral}$')
      ax.plot(x,np.sqrt(abs(out_q)),color='black',marker='o',markersize=2,linestyle='--',linewidth=1,
              label=r'$\mathbf{EFIT}$')
      ax.set_xlabel(r'$\sqrt{\psi/\psi_\mathrm{sep}}$')
      ax.set_ylabel(r'$\sqrt{q}$')
      ax.legend()
      
      plt.tight_layout()

      ofile = 'prgen_q.pdf'
      print('INFO: (prgen_contour) Writing '+ofile)
      plt.savefig(ofile)
      
   return rv,zv,psic,out_p,out_f,out_q,psi_sep
