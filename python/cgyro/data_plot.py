import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from gacodefuncs import *
import data

MYDIR=os.path.basename(os.getcwd())

class cgyrodata_plot(data.cgyrodata):

   TEXPHI  = r'{\delta\phi}'
   TEXAPAR = r'\delta {A_\parallel}'
   TEXBPAR = r'\delta {B_\parallel}'

   def kxky_select(self,theta,field,moment,species):

      # Select theta index
      if self.theta_plot == 1:
         itheta = 0
      else:
         # theta=0 check just to be safe
         if theta == 0.0:
            itheta = self.theta_plot/2
         else:
            itheta = int((theta+1.0)/2.0*self.theta_plot)
      
      if moment == 'phi':
         if field == 0:
            f = self.kxky_phi_abs[:,itheta,:,:]
            ft = self.TEXPHI
         elif field == 1:
            f = self.kxky_apar_abs[:,itheta,:,:]
            ft = self.TEXAPAR
         else:
            f = self.kxky_bpar_abs[:,itheta,:,:]
            ft = self.TEXBPAR
      elif moment == 'n':
         f = self.kxky_n_abs[:,itheta,species,:,:]
         ft = ''
      elif moment == 'e':
         f = self.kxky_e_abs[:,itheta,species,:,:]
         ft = ''

      print('INFO: (kxky_select) Selected theta index {:d} of {:d} '.
            format(itheta+1,self.theta_plot))
      return f,ft
       
         
   def plot_freq(self,w=0.5,wmax=0.0,fig=None):
      '''
      Plot gamma and omega vs time

      ARGUMENTS:
      w: fractional width of time window
      '''

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      #======================================
      # Omega
      ax = fig.add_subplot(121)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$(a/c_s)\, \omega$')

      for i in range(self.n_n):
         ax.plot(self.t,self.freq[0,i,:])

      ax.set_xlim([0,self.t[-1]])
      #======================================

      #======================================
      # Gamma
      ax = fig.add_subplot(122)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$(a/c_s)\, \gamma$')

      for i in range(self.n_n):
         ax.plot(self.t,self.freq[1,i,:])

      ax.set_xlim([0,self.t[-1]])
      #======================================

      fig.tight_layout(pad=0.3)
      
   def plot_ky_freq(self,w=0.5,wmax=0.0,fig=None):
      '''
      Plot mode frequency versus ky

      ARGUMENTS:
      w: fractional width of time window
      '''

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      #======================================
      # Omega
      ax = fig.add_subplot(121)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$k_y \rho_s$')
      ax.set_ylabel(r'$(a/c_s)\, \omega$')

      ax.plot(self.ky,self.freq[0,:,-1],color='blue')
      ax.plot(self.ky,self.freq[0,:,-1],"o",color='k')
      if len(self.ky) > 1:
         ax.set_xlim([0,self.ky[-1]])
      #======================================

      #======================================
      # Gamma
      ax = fig.add_subplot(122)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$k_y \rho_s$')
      ax.set_ylabel(r'$(a/c_s)\, \gamma$')
         
      ax.plot(self.ky,self.freq[1,:,-1],color='red')
      ax.plot(self.ky,self.freq[1,:,-1],"o",color='k')
      if len(self.ky) > 1:
         ax.set_xlim([0,self.ky[-1]])
      #======================================

      fig.tight_layout(pad=0.3)

   def plot_ky_phi(self,field=0,theta=0.0,ymin='auto',ymax='auto',nstr='null',fig=None):

      # Plot fields versus time for particular values of ky

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      self.getbigfield()

      f,ft = self.kxky_select(theta,field,'phi',0)
      p = np.sum(f[:,:,:],axis=0)/self.rho
      
      ax = fig.add_subplot(111)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$\left| '+ft+'_n \\right|$')
      ax.set_yscale('log')
      ax.set_title(r'$\mathrm{Fluctuation~intensity} \quad k_\theta = nq/r$')
        
      if nstr == 'null':
         nvec = range(self.n_n)
      else:
         nvec = str2list(nstr)
    
      for n in nvec:
         num = r'$n='+str(n)+'$'
         if n==0:
            ax.plot(self.t,p[n,:],linewidth=2,label=num)
         else:
            ax.plot(self.t,p[n,:],label=num)

      ax.set_xlim([0,max(self.t)])

      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))
      if ymin != 'auto':
         ax.set_ylim(bottom=float(ymin))

      if self.n_n > 16:
         ax.legend(loc=4, ncol=5, prop={'size':12})
      else:
         ax.legend(loc=4, ncol=6, prop={'size':12})

      fig.tight_layout(pad=0.3)

      
   def plot_rcorr_phi(self,field=0,theta=0.0,w=0.5,wmax=0.0,fig=None):
      '''
      Plot radial correlation 

      ARGUMENTS:
      w: fractional width of time window
      '''

      import scipy as scipy
      import scipy.signal as signal
      from scipy.optimize import curve_fit

      def absexp(x,tau):
         return np.exp(-np.abs(x)/tau)

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      self.getbigfield()

      t   = self.t
      kx  = self.kx
      ave = np.zeros(self.n_radial)

      imin,imax=iwindow(self.t,w,wmax)
    
      dk = kx[1]-kx[0]
      x0 = kx[-1]+dk

      color = ['m','k','b','c']
      xlabel=r'$r / \rho_s$'
      windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[imax])+']$'

      ax = fig.add_subplot(1,1,1)
      ax.set_title(r'$\mathrm{Average~radial~correlation} \quad $'+windowtxt)
      ax.set_xlabel(xlabel)

      f,ft = self.kxky_select(theta,field,'phi',0)
      y = np.sum(f[:,1:,:],axis=1)
      
      for j in range(self.n_radial):
         ave[j] = average(y[j,:],self.t,w,wmax)

      ave = np.roll(ave,-self.n_radial/2)
      ave[0] = 0.0
      corr = np.fft.fft(ave,self.n_radial)
      corr = np.fft.fftshift(corr)
      corr /= np.max(np.abs(corr))
      corr = corr.real
      delta_r = np.fft.fftfreq(self.n_radial)
      delta_r = np.fft.fftshift(delta_r)
      Lx = 2*np.pi/dk
      delta_r *= Lx

      # calculate envelope
      corr_hilbert = signal.hilbert(corr)
      corr_env = np.abs(corr_hilbert)
      ax.set_ylabel(r'$C_{'+ft+'}(\Delta r)$',color='k')
      ax.plot(delta_r, 0*delta_r, color='k', ls='--')
      ax.plot(delta_r, corr, color=color[0])

      l_corr, pcov = curve_fit(absexp, delta_r, corr_env, p0=10.0)
      ax.plot(delta_r,absexp(delta_r,l_corr),color=color[1],ls='-.')

      ax.set_xlim([np.min(delta_r),np.max(delta_r)])
      ax.set_ylim(-1,1)

      fig.tight_layout(pad=0.3)

      print('INFO: (data_plot.py) l_corr = {:d}'.format(l_corr))

   def plot_phi(self,field=0,theta=0.0,fig=None):

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      self.getbigfield()

      f,ft = self.kxky_select(theta,field,'phi',0)

      #======================================
      # Set figure size and axes
      ax = fig.add_subplot(111)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$\left|'+ft+r'\right|$')
      ax.set_yscale('log')
      ax.set_title(r'$\mathrm{Fluctuation~intensity} \quad k_\theta = nq/r$')
      #======================================

      y0 = np.sum(f[:,0,:],axis=0)/self.rho      

      # n=0 intensity
      ax.plot(self.t,y0,label=r'$n=0$',linewidth=2)

      # finite-n intensity
      yn = np.sum(f[:,1,:],axis=0)/self.rho 
      for n in range(2,self.n_n):
         yn = yn+np.sum(f[:,n,:],axis=0)/self.rho 

      ax.plot(self.t,yn,label=r'$n>0$')
        
      ax.set_xlim([0,max(self.t)])

      ax.legend(loc=4)

      head = '(cs/a) t     Phi_0/rho_*    Phi_n/rho_*'

      fig.tight_layout(pad=0.3)

      return head,self.t,y0,yn

   def plot_zf(self,w=0.5,wmax=0.0,field=0,fig=None):

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      if self.n_n > 1:
         raise ValueError('(plot_zf.py) This plot option valid for ZF test only.')

      t  = self.t
      k0 = self.kx[0]

      print('INFO: (plot_zf.py) Using index theta index n_theta/3+1')
      if field == 0:
         f = self.phib[0,self.n_theta/3,:]
      elif field == 1:
         f = self.aparb[0,self.n_theta/3,:]
      else:
         f = self.bparb[0,self.n_theta/3,:]

      # Initialization in CGYRO is with 1e-6*besselj0 # phic[0]
      gfactor = 1e6*(1-np.i0(k0**2)*np.exp(-k0**2))/(np.i0(k0**2)*np.exp(-k0**2))

      y = f*gfactor
      
      #----------------------------------------------------
      # Average calculations
      imin,imax = iwindow(t,w,wmax)
      ave  = average(y[:],t,w,wmax)
      print('INFO: (plot_zf) Integral time-average = %.6f' % ave)

      ave_vec = ave*np.ones(len(t))
      #----------------------------------------------------

      ax = fig.add_subplot(111)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$\mathrm{Re}\left( \delta\phi/\delta\phi_0 \right)$')

      ax.plot(t,y,label=r'$k_x=%.3g$' % k0)
    
      ax.plot(t[imin:],ave_vec[imin:],color='b',
              label=r'$\mathrm{Average}$',linewidth=1)

      theory = 1.0/(1.0+1.6*self.q**2/np.sqrt(self.rmin/self.rmaj))
      ax.plot([0,max(t)],[theory,theory],color='grey',
              label=r'$\mathrm{RH \; theory}$',alpha=0.3,linewidth=4)

      theory2 = 1./(1.0+2.*self.q**2)
      ax.plot([0,max(t)],[theory2,theory2],color='m',
              label=r'$\mathrm{fluid \; theory}$',alpha=0.3,linewidth=4)

      ax.legend(loc=1,prop={'size':14})

      fig.tight_layout(pad=0.3)

   def plot_geo(self,fig=None):

      self.getgeo()

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      theta = self.geo[:,0]/np.pi

      for p in range(8):
         p1 = p+1
         y = self.geo[:,p1]

         ax = fig.add_subplot(2,4,p1)
         ax.grid(which="majorminor",ls=":")
         ax.grid(which="major",ls=":")
         ax.set_xlabel(r'$\theta/\pi$')
         ax.set_title(r'$'+self.geotag[p1]+'$')
         ax.plot(theta,y)
         ax.set_xlim([-1,1])

      fig.tight_layout(pad=0.3)

   def plot_error(self,fig=None):

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      ax = fig.add_subplot(111)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$\mathrm{Integration~Error}$')
      ax.set_yscale('log')

      ax.plot(self.t[2:],self.err1[2:],label=r'$\mathrm{Total~error}$')
      ax.plot(self.t[2:],self.err2[2:],label=r'$\mathrm{RK4~error}$')
      ax.set_xlim([0,self.t[-1]])

      ax.legend()

      fig.tight_layout(pad=0.3)

   def plot_ball(self,itime=-1,field=0,tmax=-1.0,fig=None):

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      if itime > self.n_time-1:
         itime = self.n_time-1

      # Construct complex eigenfunction at selected time
      if field == 0:
         f = self.phib[0,:,itime]+1j*self.phib[1,:,itime]
         ytag = self.TEXPHI
      elif field == 1:
         f = self.aparb[0,:,itime]+1j*self.aparb[1,:,itime]
         ytag = self.TEXAPAR
      elif field == 2:
         f = self.bparb[0,:,itime]+1j*self.bparb[1,:,itime]
         ytag = self.TEXBPAR

      ax = fig.add_subplot(111)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$\theta_*/\pi$')
      ax.set_ylabel(r'$'+ytag+'$')

      if self.n_radial == 1:
         # Manage n=0 (ZF) case
         x = self.theta/np.pi
         ax.set_xlim([-1,1])
      else:
         # Assume n > 0 (ballooning mode) if n_radial > 1
         x = self.thetab/np.pi
         if tmax < 0.0:
            ax.set_xlim([1-self.n_radial,-1+self.n_radial])
         else:
            ax.set_xlim([-tmax,tmax])

      y1 = np.real(f)
      y2 = np.imag(f)

      ax.plot(x,y1,'-o',color='black',markersize=2,label=r'$\mathrm{Re}$')
      ax.plot(x,y2,'-o',color='red',markersize=2,label=r'$\mathrm{Im}$')
      
      ax.legend()

      fig.tight_layout(pad=0.3)

      return 'ang  Re(f)  Im(f)',x,y1,y2
         
   def plot_flux(self,w=0.5,wmax=0.0,field=0,moment='e',ymin='auto',ymax='auto',
                 fc=0,fig=None,loc=2,nscale=0,cflux='auto'):
         
      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      usec = self.getflux(cflux)

      ns = self.n_species
      t  = self.t

      field_tag = '\mathrm{Total}'

      # Total flux or components
      if fc == 0:
         ys = np.sum(self.ky_flux,axis=(2,3))
      else:
         ys = np.sum(self.ky_flux[:,:,field,:,:],axis=2)
         if field == 0:
            field_tag = '\phi'
         elif field == 1:
            field_tag = 'A_\parallel'
         else:
            field_tag = 'B_\parallel'

      # Now, ys -> {n_species,3,nt}

      if moment == 'n':
         ntag = 'Density~flux'
         mtag = '\Gamma'
         ttag = 'G'
         ftag = 'flux_n'
         y = ys[:,0,:]
      elif moment == 'e':
         ntag = 'Energy~flux'
         mtag = 'Q'
         ttag = 'Q'
         ftag = 'flux_e'
         y = ys[:,1,:]
      elif moment == 'v':
         ntag = 'Momentum~flux'
         mtag = '\Pi'
         ttag = 'Pi'
         ftag = 'flux_v'
         y = ys[:,2,:]
      else:
         raise ValueError('(plot_flux.py) Invalid moment.')

      if usec:
         ntag = ntag+'~(central)'

      # Normalizations
      if nscale == 0:
         norm_vec = np.ones(ns)
         mnorm = ''
      else:
         norm_vec = 1.0/self.dens
         mnorm = '^\mathrm{norm}'

      # Get index for average window
      imin,imax=iwindow(t,w,wmax)

      # Otherwise plot
      ax = fig.add_subplot(111)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)

      color = ['k','m','b','c','g','r']

      windowtxt = '['+str(t[imin])+' < (c_s/a) t < '+str(t[imax])+']'

      ax.set_title(r'$\mathrm{'+ntag+'} \quad '+windowtxt+'\quad ['+field_tag+']$')

      for ispec in range(ns):
         y_norm = y[ispec,:]*norm_vec[ispec]
         ave    = average(y_norm,t,w,wmax)
         y_ave  = ave*np.ones(len(t))
         u = specmap(self.mass[ispec],self.z[ispec])
         label = r'$'+mtag+mnorm+'_'+u+'/'+mtag+'_\mathrm{GB}: '+str(round(ave,3))+'$'
         # Average
         ax.plot(t[imin:imax+1],y_ave[imin:imax+1],'--',color=color[ispec])
         # Time trace
         ax.plot(self.t,y_norm,label=label,color=color[ispec])

      ax.legend(loc=loc)

      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))
      if ymin != 'auto':
         ax.set_ylim(bottom=float(ymin))

      fig.tight_layout(pad=0.3)

   def plot_xflux(self,w=0.5,wmax=0.0,moment='e',ymin='auto',ymax='auto',fig=None,nscale=0):

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      self.getxflux()
      
      ns = self.n_species
      nl = self.n_global+1
      t  = self.t

      ky  = self.ky
      ave = np.zeros((self.n_n,ns))

      # NOTE: lky_flux_* -> [ 2, nl , ns , n_n , nt ]
      #                       0  1    2     3    4 

      if moment == 'n':
         ntag = 'Density~flux'
         mtag = '\Gamma'
         z = np.sum(self.lky_flux_n,axis=3)
         ftag = 'xflux_n'
      elif moment == 'e':
         ntag = 'Energy~flux'
         mtag = 'Q'
         z = np.sum(self.lky_flux_e,axis=3)
         ftag = 'xflux_e'
      elif moment == 'v':
         ntag = 'Momentum~flux'
         mtag = '\Pi'
         z = np.sum(self.lky_flux_v,axis=3)
         ftag = 'xflux_v'
      else:
         raise ValueError('(plot_xflux.py) Invalid moment.')


      # Call routine for domain average
      e = 0.2
      self.xfluxave(w,moment,e=e,nscale=nscale)

      # Rescale with density ratio
      if nscale == 1:
         mnorm = '^\mathrm{norm}'
      else:
         mnorm = ''

      # Determine tmin
      imin,imax=iwindow(t,w,wmax)

      #============================================================
      # Otherwise plot
      ax = fig.add_subplot(111)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$r/L_x$')

      color = ['k','m','b','c','g','r']

      windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[imax])+']$'

      ax.set_title(r'$\mathrm{'+ntag+'} \quad $'+windowtxt)
    
      t = -np.pi+2*np.pi*np.arange(0.0,1.0,0.001)

      for ispec in range(ns):

         u = specmap(self.mass[ispec],self.z[ispec])

         # Flux curve
         g = np.zeros(len(t))
         g = self.lky_xr[ispec,0] 
         for l in range(1,nl):
            g = g+2*(np.cos(l*t)*self.lky_xr[ispec,l]-np.sin(l*t)*self.lky_xi[ispec,l])
         ax.plot(t/(2*np.pi),g,color=color[ispec])

         #---------------------------------
         # Flux partial average over [-e,e]
         g0 = self.lky_flux_ave[ispec,0]
         label = r'$'+mtag+mnorm+'_'+u+'/'+mtag+'_\mathrm{GB}: '+str(round(g0,3))+'$'
         ax.plot([-e,e],[g0,g0],'o-',color=color[ispec],alpha=0.2,linewidth=3,label=label)
         #---------------------------------

         #---------------------------------
         # Flux partial average over "negative" interval
         g1 = self.lky_flux_ave[ispec,1]
         ax.plot([0.5-e,0.5],[g1,g1],'o--',color=color[ispec],alpha=0.2,linewidth=3)
         ax.plot([-0.5,-0.5+e],[g1,g1],'o--',color=color[ispec],alpha=0.2,linewidth=3)
         #---------------------------------

         #---------------------------------
         # Flux spectral average
         gs = self.lky_xr[ispec,0]+2*np.pi/4*self.lky_xr[ispec,1]
         #---------------------------------
         
         #---------------------------------
         # Flux domain average
         ga = self.lky_xr[ispec,0]
         ax.plot([-0.5,0.5],[ga,ga],color=color[ispec],alpha=0.5)
         #---------------------------------

         print('INFO: (plot_xflux) Ave [inner/inner_spec, outer, domain] = '
               '{:.2f}/{:.2f}, {:.2f}, {:.2f}'.format(g0,gs,g1,ga)) 

         if ymax != 'auto':
            ax.set_ylim(top=float(ymax))
         if ymin != 'auto':
            ax.set_ylim(bottom=float(ymin))

         ax.axvspan(-0.25,0.25,facecolor='g',alpha=0.1)
         ax.set_xlim([-0.5,0.5])

         ax.legend(loc=2)

      fig.tight_layout(pad=0.3)

   def plot_ky_flux(self,w=0.5,wmax=0.0,field=0,moment='e',ymin='auto',ymax='auto',
                    fc=0,diss=0,fig=None,cflux='auto'):

      if self.n_n == 1:
         raise ValueError('(plot_ky_flux.py) Plot not available with a single mode.')

      ns = self.n_species
      t  = self.t

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.ly*ns,self.ly))

      usec = self.getflux(cflux)

      ky  = self.ky
      ave = np.zeros((self.n_n,ns))

      field_tag = '\mathrm{Total}'

      if fc == 0:
         ys = np.sum(self.ky_flux,axis=(2))
      else:
         ys = self.ky_flux[:,:,field,:,:]
         if field == 0:
            field_tag = '\phi'
         elif field == 1:
            field_tag = 'A_\parallel'
         else:
            field_tag = 'B_\parallel'

      if moment == 'n':
         ntag = 'Density~flux'
         mtag = '\Gamma'
         ttag = 'G'
         ftag = 'flux_n'
         y = ys[:,0,:,:]
      elif moment == 'e':
         ntag = 'Energy~flux'
         mtag = 'Q'
         ttag = 'Q'
         ftag = 'flux_e'
         y = ys[:,1,:,:]
      elif moment == 'v':
         ntag = 'Momentum~flux'
         mtag = '\Pi'
         ttag = 'Pi'
         ftag = 'flux_v'
         y = ys[:,2,:,:]
      else:
         raise ValueError('(plot_ky_flux.py) Invalid moment.')
      
      # Determine tmin
      imin,imax=iwindow(t,w,wmax)

      color = ['m','k','b','c']

      if usec:
         cstr = '~\mathrm{(central)}'
      else:
         cstr = ''
         
      windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[imax])+']'+cstr+'$'


      if ky[-1] < 0.0:
         ky = -ky
         xlabel=r'$-k_\theta \rho_s$'
      else:
         xlabel=r'$k_\theta \rho_s$'
    
      dk = ky[1]-ky[0]
    
      for ispec in range(ns):
         for j in range(self.n_n):
            ave[j,ispec] = average(y[ispec,j,:],self.t,w,wmax)

      # One plot per species
      for ispec in range(ns):
         ax = fig.add_subplot(1,ns,ispec+1)
         ax.set_xlabel(xlabel)
         u = specmap(self.mass[ispec],self.z[ispec])
         ax.set_ylabel(r'$'+mtag+'_'+u+'$',color='k')
         ax.set_title(windowtxt)
         ax.bar(ky,ave[:,ispec],width=dk/1.1,color=color[ispec],
                alpha=0.5,edgecolor='black',align='center')
         
         # Dissipation curve             
         if diss == 1:
            ax.plot(ky,self.alphadiss*ax.get_ylim()[1]*0.5,linewidth=2,color='k',alpha=0.2)

         # Maximum
         j = np.argmax(ave[:,ispec])
         if j < len(ky)-1:
            xs,ys = quadratic_max(ky[j-1:j+2],ave[j-1:j+2,ispec])
         else:
            xs = ky[-1]
         print('INFO: (data_plot.py) Max(flux) occurs at ky*rho = {:.3f}'.format(xs))
         
         # Set axis ranges
         ax.set_xlim([0,ky[-1]+dk])
         if ymax != 'auto':
            ax.set_ylim(top=float(ymax))
         if ymin != 'auto':
            ax.set_ylim(bottom=float(ymin))

      fig.tight_layout(pad=0.3)

   def plot_kxky_phi(self,field=0,theta=0.0,w=0.5,wmax=0.0,fig=None):

      from mpl_toolkits.mplot3d import Axes3D

      x0 = max(abs(self.kx))*0.25
      y0 = max(abs(self.ky))

      asp=y0/(2*x0)

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.lx*asp))
 
      self.getbigfield()

      #-----------------------------------------------------------------
      # Note array structure
      # self.phi = np.reshape(data,(2,self.n_radial,self.n_n,nt),'F')

      t  = self.t  
      nx = self.n_radial
      ny = self.n_n

      f = np.zeros([nx-1,ny])

      # Field data selector
      fx,ft = self.kxky_select(theta,field,'phi',0)

      imin,imax=iwindow(t,w,wmax)
      for i in np.arange(imin,self.n_time):
         f = f+fx[1:,:,i]
      
      # Fix (0,0)
      i0 = nx/2-1
      f[i0,0] = 1e-6

      # Reverse y order for image plotting
      f = f[:,::-1]

      # Scale data
      f = np.log(f)

      ax = fig.add_subplot(111)

      windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[imax])+']$'

      ax.set_xlabel(r'$k_x \rho_s/4$')
      ax.set_ylabel(r'$k_y \rho_s$')
      ax.set_title(r'$\mathrm{Time}$-$\mathrm{averaged~'+ft+'~intensity} \quad $'+windowtxt)

      ax.imshow(np.transpose(f),extent=[-x0,x0,0,y0],interpolation='none')

      fig.tight_layout(pad=0.5)

   def plot_kx_phi(self,field=0,theta=0.0,w=0.5,wmax=0.0,ymin='auto',ymax='auto',nstr='null',diss=0,fig=None):

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      self.getbigfield()

      t   = self.t
      kx  = self.kx
      ave = np.zeros(self.n_radial)

      imin,imax=iwindow(self.t,w,wmax)
    
      dk = kx[1]-kx[0]
      x0 = kx[-1]+dk

      ax = fig.add_subplot(1,1,1)

      color = ['m','k','b','c']
      xlabel=r'$k_x \rho_s$'
      windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[imax])+']$'

      ax.set_title(r'$\mathrm{Average~fluctuation~intensity} \quad $'+windowtxt)
      ax.set_xlabel(xlabel)

      f,ft = self.kxky_select(theta,field,'phi',0)
  
      if nstr == 'null':
         y = np.sum(f[:,:,:],axis=1)
         for j in range(self.n_radial):
            ave[j] = average(y[j,:],self.t,w,wmax)
         ax.set_ylabel(r'$\overline{'+ft+'_\mathrm{total}}$',color='k')
         ax.step(kx+dk/2,np.sqrt(ave[:]),color=color[0])
      else:
         y = np.zeros([self.n_radial,self.n_time])
         nvec = str2list(nstr)
         print('INFO: (plot_kx_phi) n = '+str(nvec))
         ax.set_ylabel(r'$\overline{'+ft+'_n}$',color='k')
         for n in nvec:
            num = r'$n='+str(n)+'$'
            y[:] = self.kxky_phi_abs[:,0,n,:]
            for j in range(self.n_radial):
               ave[j] = average(f[j,n,:],self.t,w,wmax)
            ax.plot(kx+dk/2,np.sqrt(ave[:]),ls='steps',label=num)
            if self.n_n > 16:
               ax.legend(loc=4, ncol=5, prop={'size':12})
            else:
               ax.legend(loc=4, ncol=6, prop={'size':12})

      ax.set_xlim([-x0,x0])
      ax.set_yscale('log')
      if ymin == 'auto':
         ax.set_ylim(bottom=0.5*np.sqrt(ave[-1]))
      else:
         ax.set_ylim(bottom=float(ymin))
      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))
         
         
      # Dissipation curve             
      if diss == 1:
         ax.plot(kx,self.radialdiss*ax.get_ylim()[1]*0.5,linewidth=2,color='k',alpha=0.2)

      fig.tight_layout(pad=0.3)

   def plot_hb(self,itime=-1,spec=0,tmax=-1.0,mesh=0,fig=None):
            
      import matplotlib.cm as cm
      
      u = specmap(self.mass[spec],self.z[spec])

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.lx))

      theta=0.0

      if itime > self.n_time-1:
         itime = self.n_time-1

      # Compute index for theta value in pitch angle and energy plots
      i0 = int(round((1.0+float(theta))*self.n_theta/2.0))
      if i0 > self.n_theta-1:
         i0 = self.n_theta-1

      if self.n_radial > 1:
         n0 = (self.n_radial/2)*self.n_theta+i0
         x = self.thetab/np.pi
      else:
         n0 = self.n_theta/3
         x = self.theta/np.pi
        
      if tmax < 0.0:
         if self.n_radial == 1:
            tmax = 1.0
         else:
            tmax = self.n_radial-1

      p = 0
      for row in range(3):

         p = p+1

         if row == 0:
            ie = 0
         if row == 1:
            ie = self.n_energy/2
         if row == 2:
            ie = self.n_energy-1

         #======================================
         ax = fig.add_subplot(3,2,p)

         ax.set_title(r'${\rm Re} \, h_'+u+' \quad \mathrm{ie}='+str(ie)+'$')
         ax.set_xlabel(r'$\theta/\pi$')
         ax.set_ylabel(r'$\xi = v_\parallel/v$')

         hp = np.transpose(np.array(self.hb[0,:,spec,:,ie,itime]))
         h_norm = 0.5*(hp[self.n_xi/2-1,n0]+hp[self.n_xi/2,n0])
         hp = hp/h_norm
         hmin = hp.min()
         hmax = hp.max()
         dh = (hmax-hmin)/100.0

         levels = np.arange(hmin-dh,hmax+dh,dh)

         ax.contourf(x,self.xi,hp,levels,cmap=cm.jet,origin='lower')
         ax.set_xlim([-tmax,tmax])

         # Plot dots for mesh points
         if row == 1 and mesh == 1:
            for i in range(self.n_theta*self.n_radial):
               for j in range(self.n_xi):
                  ax.plot([self.thetab[i]/np.pi],[self.xi[j]],marker='.',color='k',markersize=4)

         #======================================
         p = p+1
         #======================================
         ax = fig.add_subplot(3,2,p)

         ax.set_title(r'${\rm Im} \, h_'+u+' \quad \mathrm{ie}='+str(ie)+'$')
         ax.set_xlabel(r'$\theta/\pi$')
         ax.set_ylabel(r'$\xi = v_\parallel/v$')

         hp = np.transpose(np.array(self.hb[1,:,spec,:,ie,itime]))
         hmin = hp.min()
         hmax = hp.max()
         dh = (hmax-hmin)/100.0

         levels = np.arange(hmin-dh,hmax+dh,dh)

         ax.contourf(x,self.xi,hp,levels,cmap=cm.jet,origin='lower')
         ax.set_xlim([-tmax,tmax])
         
         #======================================

      fig.tight_layout(pad=0.3)

   def plot_hbcut(self,itime=-1,spec=0,tmax=-1.0,theta=0.0,fig=None):

      u = specmap(self.mass[spec],self.z[spec])

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.lx))
       
      if itime > self.n_time-1:
         itime = self.n_time-1

      func = self.hb

      # Compute index for theta value in pitch angle and energy plots
      i0 = int(round((1.0+theta)*self.n_theta/2.0))
      if i0 > self.n_theta-1:
         i0 = self.n_theta-1

      if self.n_radial > 1:
         x = self.thetab/np.pi
      else:
         x = self.theta/np.pi

      if tmax < 0.0:
         if self.n_radial == 1:
            tmax = 1.0
         else:
            tmax = self.n_radial-1

      p = 0
      for row in range(3):

         p = p+1

         if row == 0:
            ie = 0
            ix = 0
         if row == 1:
            ie = self.n_energy/2
            ix = self.n_xi/2
         if row == 2:
            ie = self.n_energy-1
            ix = self.n_xi-1

         #========================================================
         ax = fig.add_subplot(3,3,p)
         ax.grid(which="majorminor",ls=":")
         ax.grid(which="major",ls=":")

         ax.set_title(r'$'+u+': \\xi=0 \quad \mathrm{ie}='+str(ie)+'$')
         ax.set_xlabel(r'$\theta/\pi$')

         if self.n_xi%2 == 0:
            hp = np.array(func[:,:,spec,self.n_xi/2,ie,itime]+
                          func[:,:,spec,self.n_xi/2-1,ie,itime])*0.5
         else:
            hp = np.array(func[:,:,spec,self.n_xi/2,ie,itime])
            
         ax.plot(x,hp[0,:],'-o',color='black',markersize=2)
         ax.plot(x,hp[1,:],'-o',color='blue',markersize=2)
         ax.set_xlim([-tmax,tmax])

         #========================================================

         p = p+1

         #========================================================
         ax = fig.add_subplot(3,3,p)
         ax.grid(which="majorminor",ls=":")
         ax.grid(which="major",ls=":")

         ax.set_title(r'$'+u+': \\theta/\pi='+str(theta)+' \quad \mathrm{ie}='+str(ie)+'$')
         ax.set_xlabel(r'$\xi = v_\parallel/v$')

         n0 = (self.n_radial/2)*self.n_theta+i0
         
         hp = np.array(func[0,:,spec,:,ie,itime])
         ax.plot(self.xi,hp[n0,:],'-o',color='black',markersize=2)
         hp = np.array(func[1,:,spec,:,ie,itime])
         ax.plot(self.xi,hp[n0,:],'-o',color='blue',markersize=2)
         ax.set_xlim([-1,1])
         #========================================================

         p = p+1

         #========================================================
         ax = fig.add_subplot(3,3,p)
         ax.grid(which="majorminor",ls=":")
         ax.grid(which="major",ls=":")

         ax.set_title(r'$'+u+': \\theta/\pi='+str(theta)+' \quad \mathrm{ix}='+str(ix)+'$')
         ax.set_xlabel(r'$x=\sqrt{\varepsilon}$')

         n0 = (self.n_radial/2)*self.n_theta+i0

         hpa = np.array(func[0,n0,spec,:,:,itime])
         hpb = np.array(func[1,n0,spec,:,:,itime])
         p0 = np.zeros(self.n_energy)
         p1 = np.zeros(self.n_energy)
         if row == 0:
            for ix in range(self.n_xi):
               p0 = p0+hpa[ix,:]*0.5
               p1 = p1+hpb[ix,:]*0.5
         elif row == 1:
            for ix in range(self.n_xi):
               p0 = p0+hpa[ix,:]*self.xi[ix]*1.5
               p1 = p1+hpb[ix,:]*self.xi[ix]*1.5
         elif row == 2:
            for ix in range(self.n_xi):
               p0 = p0+hpa[ix,:]*(3*self.xi[ix]**2-1)/2*2.5
               p1 = p1+hpb[ix,:]*(3*self.xi[ix]**2-1)/2*2.5
            
         ax.plot(np.sqrt(self.energy),p0,'-o',color='black',markersize=2)
         ax.plot(np.sqrt(self.energy),p1,'-o',color='blue',markersize=2)
         #========================================================

      fig.tight_layout(pad=0.3)

   def plot_hball(self,itime=-1,spec=0,tmax=-1.0,ymin='auto',ymax='auto',nstr='null',ie=0,fig=None):

      if nstr == 'null':
         nvec = range(self.n_n)
      else:
         nvec = str2list(nstr)

      u = specmap(self.mass[spec],self.z[spec])

      # Diagnostics
      print('l    = '+nstr)
      print('e    = '+str(ie))
      print('spec = '+u)

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))
       
      if itime > self.n_time-1:
         itime = self.n_time-1

      ax = fig.add_subplot(111)
      ax.grid(which="majorminor",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$\theta/\pi$')

      x = self.thetab/np.pi
      if tmax < 0.0:
         ax.set_xlim([1-self.n_radial,-1+self.n_radial])
      else:
         ax.set_xlim([-tmax,tmax])

      # y = y[re/im,theta,xi]
      try:
         y = np.array(self.hb[:,:,spec,:,ie,itime])
      except:
         raise ValueError('(plot_hball.py) Need to run with H_PRINT_FLAG=1')
         
      xp,wp = np.polynomial.legendre.leggauss(self.n_xi)
      c = np.zeros(self.n_xi)
      alr = np.zeros([len(x)])
      ali = np.zeros([len(x)])
      cvec = ['black','red','blue','green','purple','magenta']
      for l in nvec:
         c[:] = 0.0 ; c[l] = 1.0
         pl = np.polynomial.legendre.legval(xp,c)
         for j in range(len(x)):
            alr[j] = (l+0.5)*np.sum(pl[:]*y[0,j,:])
            ali[j] = (l+0.5)*np.sum(pl[:]*y[1,j,:])
      
         ax.plot(x,alr,'-',color=cvec[l%6],label=r'$\ell='+str(l)+'$')
         ax.plot(x,ali,'--',color=cvec[l%6])

      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))
      if ymin != 'auto':
         ax.set_ylim(bottom=float(ymin))

      ax.legend(loc=1)
      fig.tight_layout(pad=0.3)
