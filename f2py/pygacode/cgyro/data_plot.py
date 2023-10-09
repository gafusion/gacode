import sys
import os
import numpy as np
import scipy.signal as signal
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rc
from ..gacodefuncs import *
from . import data

MYDIR=os.path.basename(os.getcwd())

class cgyrodata_plot(data.cgyrodata):

   def plot_dictinit(self):

      # Init for OMFIT

      xin = {}

      xin['fig']    = None
      xin['lx']     = 12
      xin['ly']     = 6
      xin['w']      = '0.5'
      xin['norm']   = 'elec'
      xin['ftype']  = 'screen'
      xin['itime']  = -1
      xin['field']  = 0
      xin['moment'] = 'phi'
      xin['tmax']   = -1.0
      xin['theta']  = -1
      xin['ymin']   = 'auto'
      xin['ymax']   = 'auto'
      xin['kxmin']  = 'auto'
      xin['kxmax']  = 'auto'
      xin['nstr']   = 'null'
      xin['abs']    = 0
      xin['fc']     = 0
      xin['loc']    = 2
      xin['nscale'] = 0
      xin['cflux']  = 'auto'
      xin['spec']   = 0
      xin['bar']    = 0
      xin['ie']     = 0
      xin['mesh']   = 0

      return xin

   def plot_freq(self,xin):

      # Function: plot gamma and omega vs time

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      # set normalizations
      t = self.getnorm(xin['norm'])

      #======================================
      # Omega
      ax = fig.add_subplot(121)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(self.tstr)
      ax.set_ylabel(self.fstr[0])

      for i in range(self.n_n):
         ax.plot(t,self.fnorm[0,i,:])

      ax.set_xlim([0,t[-1]])
      #======================================

      #======================================
      # Gamma
      ax = fig.add_subplot(122)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(self.tstr)
      ax.set_ylabel(self.fstr[1])

      for i in range(self.n_n):
         ax.plot(t,self.fnorm[1,i,:])

      ax.set_xlim([0,t[-1]])
      #======================================

      fig.tight_layout(pad=0.3)


   def plot_ky_freq(self,xin):

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      t = self.getnorm(xin['norm']) ; ky = self.kynorm

      #======================================
      # Omega
      ax = fig.add_subplot(121)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(self.kystr)
      ax.set_ylabel(self.fstr[0])

      y1 = self.fnorm[0,:,-1] 
      ax.plot(ky,y1,color='blue')
      ax.plot(ky,y1,"o",color='k')
      if len(ky) > 1:
         ax.set_xlim([0,ky[-1]])
      #======================================

      #======================================
      # Gamma
      ax = fig.add_subplot(122)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(self.kystr)
      ax.set_ylabel(self.fstr[1])

      y2 = self.fnorm[1,:,-1]
      ax.plot(ky,y2,color='red')
      ax.plot(ky,y2,"o",color='k')
      if len(ky) > 1:
         ax.set_xlim([0,ky[-1]])
      #======================================

      fig.tight_layout(pad=0.3)

      return 'ky            omega            gamma',ky,y1,y2


   def plot_error(self,xin):

      t = self.getnorm(xin['norm'])

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(self.tstr)
      ax.set_ylabel(r'$\mathrm{Integration~Error}$')
      ax.set_yscale('log')

      ax.plot(t[2:],self.err1[2:],label=r'$\mathrm{Total~error}$')
      ax.plot(t[2:],self.err2[2:],label=r'$\mathrm{RK~error}$')
      ax.set_xlim([0,t[-1]])

      ax.legend()

      fig.tight_layout(pad=0.3)

      return


   def plot_geo(self,xin):

      self.getgeo()

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(1.2*xin['lx'],1.2*xin['ly']))

      # Decrease font size a bit for this plot
      rc('font',size=12)
      # Create 3x4 subplot grid
      gs = GridSpec(3,4)

      theta = self.geo[:,0]/np.pi

      # CGYRO geometry functions
      for p in range(9):
         p1 = p+1
         if p < 4:
            a = 1.0
         else:
            a = 1.0/self.rho
         y = a*self.geo[:,p1]

         ax = fig.add_subplot(gs[p//3,np.mod(p,3)])
         ax.grid(which="both",ls=":")
         ax.grid(which="major",ls=":")
         ax.set_xlabel(r'$\theta/\pi$')
         ax.set_title(r'$'+self.geotag[p1]+'$')
         ax.plot(theta,y,'m')
         ax.plot(theta,y,'o',color='k',markersize=2)
         ax.set_xlim([-1,1])

      # captheta 
      if len(self.geo[0,:]) > 12:
         y = self.geo[:,12]
         ax = fig.add_subplot(gs[0,3])
         ax.grid(which="both",ls=":")
         ax.grid(which="major",ls=":")
         ax.set_xlabel(r'$\theta/\pi$')
         ax.set_title(r'$'+self.geotag[12]+'$')
         ax.plot(theta,y,'m')
         ax.plot(theta,y,'o',color='k',markersize=2)
         ax.set_xlim([-1,1])
         ax.set_ylim([y[0],-y[0]])

      # Flux surface
      ax = fig.add_subplot(gs[1:,3],aspect='equal')
      ax.set_title(r'$r/a='+str(self.rmin)+'$')
      ax.set_facecolor('lightcyan')
      ax.set_xlabel(r'$R$')
      ax.set_ylabel(r'$Z$')
      t = 2*np.pi*np.linspace(0,1,200)
      rmaj = self.rmin
      zmaj = self.zmag
      r = self.rmin
      k = self.kappa

      # Map legacy shape parameters
      self.shape_sin[1] = np.arcsin(self.delta)
      self.shape_sin[2] = -self.zeta

      arg = t+self.shape_cos[0]
      for p in range(1,7):
         arg = arg+self.shape_sin[p]*np.sin(p*t)+self.shape_cos[p]*np.cos(p*t)

      x = rmaj+r*np.cos(arg)
      y = zmaj+k*r*np.sin(t)
      
      ax.plot(x,y,'k')
   
      fig.tight_layout(pad=0.3)

      return

   
   def plot_ball(self,xin):

      itime = xin['itime']
      field = xin['field']
      tmax  = xin['tmax']

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      if itime > self.n_time-1:
         itime = self.n_time-1

      # Construct complex eigenfunction at selected time
      if field == 0:
         f = self.phib[0,:,itime]+1j*self.phib[1,:,itime]
         ytag = r'$\delta\phi$'
      elif field == 1:
         f = self.aparb[0,:,itime]+1j*self.aparb[1,:,itime]
         ytag = r'$A_\parallel$'
      elif field == 2:
         f = self.bparb[0,:,itime]+1j*self.bparb[1,:,itime]
         ytag = r'$B_\parallel$'

      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$\theta_*/\pi$')
      ax.set_ylabel(ytag)

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

   
   def plot_ky_phi(self,xin):

      # Plot fields versus time for each ky

      theta  = xin['theta']
      field  = xin['field']
      nstr   = xin['nstr']
      ymin   = xin['ymin']
      ymax   = xin['ymax']
      moment = xin['moment']
      spec   = xin['spec']

      t = self.getnorm(xin['norm'])  

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      self.getbigfield()

      f,ft = self.kxky_select(theta,field,moment,spec,gbnorm=True)

      p = np.sum(abs(f[:,:,:]),axis=0)

      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(self.tstr)
      ax.set_ylabel(self.ylabeler('n',ft))
      ax.set_yscale('log')
      ax.set_title(r'$\mathrm{Fluctuation~intensity}$')

      if nstr == 'null':
         nvec = list(range(self.n_n))
      else:
         nvec = str2list(nstr)

      for n in nvec:
         num = '$n='+str(n)+'$'
         if n==0:
            ax.plot(t,p[n,:],linewidth=2,label=num)
         else:
            ax.plot(t,p[n,:],label=num)

      ax.set_xlim([0,max(t)])

      if ymin != 'auto':
         ax.set_ylim(bottom=float(ymin))
      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))

      if self.n_n > 16:
         ax.legend(loc=4, ncol=5, prop={'size':12})
      else:
         ax.legend(loc=4, ncol=6, prop={'size':12})

      fig.tight_layout(pad=0.3)


   def plot_phi(self,xin):

      w       = xin['w']
      theta   = xin['theta']
      field   = xin['field']
      ymin    = xin['ymin']
      ymax    = xin['ymax']
      absnorm = xin['abs']
      moment  = xin['moment']
      spec    = xin['spec']

      t = self.getnorm(xin['norm'])  

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      self.getbigfield()
      
      f,ft = self.kxky_select(theta,field,moment,spec,gbnorm=True)

      #======================================
      # Set figure size and axes
      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(self.tstr)
      ax.set_yscale('log')
      ax.set_title(r'$\mathrm{Fluctuation~intensity}$')
      #======================================

      # Get index for average window
      imin,imax=time_index(t,w)
    
      # n=0 intensity (gyroBohm)
      y0 = np.sum(abs(f[:,0,:]),axis=0)
      s0 = np.sum(abs(f[:,0,:])**2,axis=0)  

      # finite-n intensity (gyroBohm)
      yn = np.sum(abs(f[:,1:,:]),axis=(0,1))
      sn = np.sum(abs(f[:,1:,:])**2,axis=(0,1))
      
      s = np.ones(imax-imin+1)
      rstr = self.rhoi
      if absnorm == 0:
         s0_ave = time_average(s0,t,imin,imax)
         sn_ave = time_average(sn,t,imin,imax)
         lab0 = self.ylabeler('0',ft,sq=True,tave=True,sqrt=True)
         labn = self.ylabeler('+',ft,sq=True,tave=True,sqrt=True)
         print('INFO: (plot_phi) [RMS NORM n=0] = {:.4f}'.format(np.sqrt(s0_ave)))
         print('INFO: (plot_phi) [RMS NORM n>0] = {:.4f}'.format(np.sqrt(sn_ave)))
         ax.plot(t,np.sqrt(s0),label=lab0,linewidth=2)
         ax.plot(t,np.sqrt(sn),label=labn)
         ax.plot(t[imin:imax+1],np.sqrt(s0_ave)*s,'--k')
         ax.plot(t[imin:imax+1],np.sqrt(sn_ave)*s,'--k')
      else:
         y0_ave = time_average(y0,t,imin,imax)
         yn_ave = time_average(yn,t,imin,imax)
         lab0 = self.ylabeler('0',ft,tave=True)
         labn = self.ylabeler('+',ft,tave=True)
         print('INFO: (plot_phi) [ABS NORM n=0] = {:.4f}'.format(y0_ave))
         print('INFO: (plot_phi) [ABS NORM n>0] = {:.4f}'.format(yn_ave))
         ax.plot(t,y0,label=lab0,linewidth=2)
         ax.plot(t,yn,label=labn)
         ax.plot(t[imin:imax+1],y0_ave*s,'--k')
         ax.plot(t[imin:imax+1],yn_ave*s,'--k')

      ax.set_xlim(0,t[-1])
      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))
      if ymin != 'auto':
         ax.set_ylim(bottom=float(ymin))
      ax.legend(loc=3)

      # JC: labels only correct for default norm
      head = '(cs/a) t     Phi_0/rho_*    Phi_n/rho_*'

      fig.tight_layout(pad=0.3)

      return head,self.t,y0,yn


   def plot_flux(self,xin):
      
      w      = xin['w']
      field  = xin['field']
      moment = xin['moment']
      ymin   = xin['ymin']
      ymax   = xin['ymax']
      fc     = xin['fc']
      ftype  = xin['ftype']
      loc    = xin['loc']
      nscale = xin['nscale']
      cflux  = xin['cflux']
      norm   = xin['norm']
      theta  = xin['theta']

      
      if xin['fig'] is None and ftype != 'nox':
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      if moment == 'phi':
         moment = 'e'
  
      t    = self.getnorm(xin['norm'])
      usec = self.getflux(cflux)
      ns   = self.n_species

      # Total flux or components
      if fc == 0:
         ys = np.sum(self.ky_flux,axis=(2,3))
      else:
         ys = np.sum(self.ky_flux[:,:,field,:,:],axis=2)

      # Now, ys -> {n_species,3,nt}

      if moment == 'n':
         mtag = '\Gamma'
         ttag = 'G'
         ftag = 'flux_n'
         y = ys[:,0,:]
      elif moment == 'e':
         mtag = 'Q'
         ttag = 'Q'
         ftag = 'flux_e'
         y = ys[:,1,:]/self.qc
      elif moment == 'v':
         # JC: correct for other norms?
         mtag = '\Pi'
         ttag = 'Pi'
         ftag = 'flux_v'
         y = ys[:,2,:]
      elif moment == 's':
         # JC: correct for other norms?
         mtag = 'S'
         ttag = 'S'
         ftag = 'exch'
         y = ys[:,3,:]
      else:
         raise ValueError('(plot_flux.py) Unrecognized moment. Try: n,e,v,s')

      # Normalizations
      if nscale == 0:
         norm_vec = np.ones(ns)
         mnorm = ''
      else:
         norm_vec = 1.0/self.dens
         mnorm = '^\mathrm{norm}'

      color = ['k','m','b','c','g','r']

      imin,imax=time_index(t,w)
      mpre,mwin = wintxt(imin,imax,t,usec=usec,fc=fc,field=field)

      # Otherwise plot
      if not ftype == 'nox':
         ax = fig.add_subplot(111)
         ax.grid(which="both",ls=":")
         ax.grid(which="major",ls=":")
         ax.set_xlabel(self.tstr)
         ax.set_title(mpre+mwin)

      for ispec in range(ns):
         y_norm = y[ispec,:]*norm_vec[ispec]
         ave   = time_average(y_norm,t,imin,imax)
         y_ave = ave*np.ones(len(t))
         u = specmap(self.mass[ispec],self.z[ispec])
         label = r'$'+mtag+mnorm+'_\mathrm{'+u+'}/'+mtag+self.gbnorm+': '+str(round(ave,3))+'$'
         if not ftype == 'nox':
            # Average
            ax.plot(t[imin:imax+1],y_ave[imin:imax+1],'--',color=color[ispec])
            # Time trace
            ax.plot(t,y_norm,label=label,color=color[ispec])

      if not ftype == 'nox':
         ax.legend(loc=loc)
         if ymax != 'auto':
            ax.set_ylim(top=float(ymax))
         if ymin != 'auto':
            ax.set_ylim(bottom=float(ymin))
         fig.tight_layout(pad=0.3)

      title = '               '
      for ispec in range(ns):
         title = title+specmap(self.mass[ispec],self.z[ispec])+'       '
      print(title)

      tag = [
         'GAMMA [GB]',
         'Q     [GB]',
         'PI    [GB]',
         'S     [GB]']
      for i in range(self.n_flux):
         bstr=''
         for ispec in range(ns):
            ave = time_average(ys[ispec,i,:],t,imin,imax) ; var=0
            bstr = bstr+"{:7.3f}".format(ave)+' '
         print(tag[i]+' '+bstr)


   def plot_ky_flux(self,xin):

      w      = xin['w']
      field  = xin['field']
      moment = xin['moment']
      ymin   = xin['ymin']
      ymax   = xin['ymax']
      fc     = xin['fc']
      ftype  = xin['ftype']
      loc    = xin['loc']
      nscale = xin['nscale']
      cflux  = xin['cflux']
      norm   = xin['norm']
      bar    = xin['bar']
      theta  = xin['theta']

      t = self.getnorm(xin['norm'])
      
      if self.n_n == 1:
         print('ERROR: (plot_ky_flux) Plot not available with a single mode.')
         return

      if moment == 'phi':
         moment = 'e'
      
      ns = self.n_species
      nn = self.n_n

      if xin['fig'] is None and ftype != 'nox':
         if ns < 4:
            nrow = 1 ; ncol = ns
         elif ns == 4:
            nrow = 2 ; ncol = 2
         elif ns > 4:
            nrow = 2 ; ncol = 3
         fig = plt.figure(MYDIR,figsize=(xin['ly']*ncol,xin['ly']*nrow))

      usec = self.getflux(cflux)
      
      ky = self.kynorm

      if fc == 0:
         ys = np.sum(self.ky_flux,axis=2)
      else:
         ys = self.ky_flux[:,:,field,:,:]

      if moment == 'n':
         mtag = '\Gamma'
         ttag = 'G'
         ftag = 'flux_n'
         y = ys[:,0,:,:]
      elif moment == 'e':
         mtag = 'Q'
         ttag = 'Q'
         ftag = 'flux_e'
         y = ys[:,1,:,:]/self.qc
      elif moment == 'v':
         mtag = '\Pi'
         ttag = 'Pi'
         ftag = 'flux_v'
         y = ys[:,2,:,:]
      elif moment == 's':
         mtag = 'S'
         ttag = 'S'
         ftag = 'exch'
         y = ys[:,3,:,:]
      elif moment == 'k':
         #  JC: needs work
         mtag = 'I'
         ttag = 'I'
         ftag = 'i'
         self.getbigfield()
         # I and K
         ns = 2
         y = np.zeros([ns,nn,self.n_time])
         f,ft = self.kxky_select(theta,field,'phi',0)
         y[0,:,:] = np.sum(abs(f[:,:,:]),axis=0)
         f,ft = self.kxky_select(theta,field,'k',0)
         y[1,:,:] = np.sum(abs(f[:,:,:]),axis=0)
      else:
         raise ValueError('(plot_ky_flux.py) Invalid moment.')
      
      color = ['magenta','k','blue','cyan','red','green']

      imin,imax=time_index(t,w)
      mpre,mwin = wintxt(imin,imax,t,usec=usec,fc=fc,field=field)

      if ky[-1] < 0.0:
         ky = -ky

      dk = ky[1]-ky[0]
    
      ave = np.zeros((self.n_n,ns))
      for ispec in range(ns):
         ave[:,ispec] = time_average(y[ispec,:,:],t,imin,imax)

      # One plot per species
      for ispec in range(ns):
         u = specmap(self.mass[ispec],self.z[ispec])
         if not ftype == 'nox':
            ax = fig.add_subplot(nrow,ncol,ispec+1)
               
            ax.set_xlabel(self.kystr)
            ax.set_ylabel(r'$'+mtag+'_\mathrm{'+u+'}/'+mtag+self.gbnorm+'$',color='k')

            if ispec < ncol:
               ax.set_title(mpre+mwin,fontsize=16)
            if bar == True:
               ax.bar(ky,ave[:,ispec],width=dk/1.1,color=color[ispec],
                      alpha=0.5,edgecolor='black',align='center')
               # Set axis ranges
               ax.set_xlim([0,ky[-1]+dk])
            else:
               ax.grid(which="both",ls=":")
               ax.grid(which="major",ls=":")
               ax.set_xscale('log')
               ax.set_yscale('log')
               ax.plot(ky[1:],ave[1:,ispec],'o',color='k',markersize=2,alpha=0.6)
               ax.plot(ky[1:],ave[1:,ispec],color=color[ispec])
               
            if ymax != 'auto':
               ax.set_ylim(top=float(ymax))
            if ymin != 'auto':
               ax.set_ylim(bottom=float(ymin))
               
         # Maximum
         j = np.argmax(ave[:,ispec])
         if j < len(ky)-1 and j > 0:
            xs,ys = quadratic_max(ky[j-1:j+2],ave[j-1:j+2,ispec])
         else:
            xs = ky[-1]
         print('INFO: (data_plot.py) Max(flux) occurs at ky*rho = {:.3f}'.format(xs))

      if not ftype == 'nox':
         fig.tight_layout(pad=0.4)
         
         
   def plot_rcorr_phi(self,xin):

      # Should be checked by Chris after updated form of kxky_select to n_radial-1 points
      
      def absexp(x,tau):
         return np.exp(-np.abs(x)/tau)
      
      w     = xin['w']
      theta = xin['theta']
      field = xin['field']

      t = self.getnorm(xin['norm'])

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      self.getbigfield()
 
      kx  = self.kxnorm
      nx  = self.n_radial
      
      imin,imax=time_index(self.t,w)

      dk = kx[1]-kx[0]
      x0 = kx[-1]+dk

      color  = ['m','k','b','c']
      xlabel = r'$r/'+self.rhoi+'$'
      mpre,mwin = wintxt(imin,imax,t)

      ax = fig.add_subplot(1,1,1)
      ax.set_title(r'$\mathrm{Average~radial~correlation} \quad $'+mwin)
      ax.set_xlabel(xlabel)

      f,ft = self.kxky_select(theta,field,'phi',0)
      # This puts the 0 element back in the radial direction
      y = np.zeros([nx,self.n_time])
      y[1:,:] = np.sum(abs(f[:,1:,:]),axis=1)

      ave = time_average(y,t,imin,imax)

      ave = np.roll(ave,-nx//2)
      ave[0] = 0.0
      corr = np.fft.fft(ave,nx)
      corr = np.fft.fftshift(corr)
      corr /= np.max(np.abs(corr))
      corr = corr.real
      delta_r = np.fft.fftfreq(nx)
      delta_r = np.fft.fftshift(delta_r)
      Lx = 2*np.pi/dk
      delta_r *= Lx

      # calculate envelope
      corr_hilbert = signal.hilbert(corr)
      corr_env = np.abs(corr_hilbert)
      ax.set_ylabel(r'$C_{'+ft+'}(\Delta r)$',color='k')
      ax.plot(delta_r,0*delta_r,color='k',ls='--')
      ax.plot(delta_r,corr,color=color[0])

      l_corr, pcov = curve_fit(absexp, delta_r, corr_env, p0=10.0)
      ax.plot(delta_r,absexp(delta_r,l_corr),color=color[1],ls='-.')

      ax.set_xlim([np.min(delta_r),np.max(delta_r)])
      ax.set_ylim(-1,1)

      fig.tight_layout(pad=0.3)

      print('INFO: (rcorr_phi) l_corr = {:.3f}'.format(l_corr[0]))


   def plot_low(self,xin):

      w      = xin['w']
      theta  = xin['theta']
      moment = xin['moment']
      field  = xin['field']
      ymin   = xin['ymin']
      ymax   = xin['ymax']
      ftype  = xin['ftype']
      spec   = xin['spec']

      t = self.getnorm(xin['norm'])
      
      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      self.getbigfield()

      #======================================
      # Set figure size and axes
      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(self.tstr)
      #======================================

      color = ['k','m','b','c','g','r']

      # Get index for average window
      imin,imax=time_index(t,w)

      mpre,mwin = wintxt(imin,imax,t)

      ax.set_title(mwin)

      nx = self.n_radial-1
      # p=1
      p1 = nx//2+1
      
      f,ft = self.kxky_select(theta,field,moment,spec,gbnorm=True)

      ax.set_ylabel(self.ylabeler('1',ft))

      yr = np.real(f[p1,0,:]) ; yi = np.imag(f[p1,0,:])
      # Re
      ave = time_average(yr,t,imin,imax) ; y_ave = ave*np.ones(len(t))
      ax.plot(self.t,yr,color=color[0],label=r'$\mathrm{Re} = '+str(round(ave,3))+'$')
      ax.plot(t[imin:imax+1],y_ave[imin:imax+1],'--',color=color[0])
      # Im
      ave = time_average(yi,t,imin,imax) ; y_ave = ave*np.ones(len(t))
      ax.plot(self.t,yi,color=color[1],label=r'$\mathrm{Im} = '+str(round(ave,3))+'$' )
      ax.plot(t[imin:imax+1],y_ave[imin:imax+1],'--',color=color[1])

      ax.legend(loc=2)

      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))
      if ymin != 'auto':
         ax.set_ylim(bottom=float(ymin))

      fig.tight_layout(pad=0.3)

      return
   

   def plot_corrug(self,xin):

      norm   = xin['norm']
      w      = xin['w']
      theta  = xin['theta']
      moment = xin['moment']
      ymin   = xin['ymin']
      ymax   = xin['ymax']
      ftype  = xin['ftype']
      spec   = xin['spec']

      t = self.getnorm(xin['norm'])

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      self.getbigfield()
      self.getnorm(norm)
      
      nx = self.n_radial-1
      
      #======================================
      # Set figure size and axes
      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$r/L$')
      #======================================

      color = ['k','m','b','c','g','r']

      # Get index for average window
      imin,imax=time_index(t,w)
      mpre,mwin = wintxt(imin,imax,t)

      ax.set_title(mwin)
      
      # f[p,n,t]
      f,ft = self.kxky_select(theta,0,moment,spec,gbnorm=True)

      # complex n=0 amplitudes
      yr = time_average(np.real(f[:,0,:]),t,imin,imax)
      yi = time_average(np.imag(f[:,0,:]),t,imin,imax)
      
      nxp = 8*nx
      fave = np.zeros(nxp)
      dave = np.zeros(nxp)
      x = np.linspace(0,2*np.pi,nxp)
      y = yr+1j*yi
      # JC: should use array notation, check definitions
      for i in range(nx):
         p = i-nx//2
         fave = fave + np.real(     np.exp(1j*p*x)*y[i])
         dave = dave + np.real(1j*p*np.exp(1j*p*x)*y[i])
      dave = dave*(2*np.pi/self.length)
         
      ax.plot(x/(2*np.pi),fave,color='k',label=self.ylabeler('0',ft,abs=False))
      ax.plot(x/(2*np.pi),dave,color='b',label=self.ylabeler('0',ft,abs=False,d=True))
      ax.legend()
      
      ax.set_xlim([0,1])

      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))
      if ymin != 'auto':
         ax.set_ylim(bottom=float(ymin))

      fig.tight_layout(pad=0.3)

      return

   def plot_shift(self,xin):

      # JC: not corrected for norm. Do we need this function?
      
      w     = xin['w']
      theta = xin['theta']
      ymin  = xin['ymin']
      ymax  = xin['ymax']

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))
            
      t = self.getnorm(xin['norm'])  
      self.getbigfield()
      
      # Get index for average window
      imin,imax=time_index(t,w)
      mpre,mwin = wintxt(imin,imax,t)

      #======================================
      # Set figure size and axes
      ax = fig.add_subplot(111)
      ax.set_title(mwin)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(self.kystr)
      ax.set_ylabel(r'$\left\langle k_x \rho_s \right\rangle$')
      #======================================

      ky = abs(self.ky)
      k0 = 2*np.pi/self.length

      f,ft = self.kxky_select(theta,0,'phi',0)

      #y1,y2 = shift_legendre(f,imin,imax)
      y1,y2 = shift_fourier(f,imin,imax)
      ax.plot(ky,k0*y1,color='k')
      ax.plot(ky,-k0*y2,linestyle='--',color='k')

      # EAB print
      #for i in range(len(ky)):
      #   print(ky[i],k0*y1[i],-k0*y2[i])

      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))
      if ymin != 'auto':
         ax.set_ylim(bottom=float(ymin))

      fig.tight_layout(pad=0.3)

      return '   ky*rho       kx*rho',ky,y1,None

   def plot_zf(self,xin):

      w     = xin['w']
      field = xin['field']
    
      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      if self.n_n > 1:
         print('ERROR: (plot_zf) This plot option valid for ZF test only.')
         return

      t  = self.t
      k0 = 2*np.pi/self.length

      print('INFO: (plot_zf) Using index theta index n_theta/3+1')
      if field == 0:
         f = self.phib[0,self.n_theta//3,:]
      elif field == 1:
         f = self.aparb[0,self.n_theta//3,:]
      else:
         f = self.bparb[0,self.n_theta/3,:]

      # Initialization in CGYRO is with 1e-6*besselj0 # phic[0]
      gfactor = 1e6*(1-np.i0(k0**2)*np.exp(-k0**2))/(np.i0(k0**2)*np.exp(-k0**2))

      y = f*gfactor
      
      #----------------------------------------------------
      # Average calculations
      imin,imax=time_index(t,w)
      ave  = time_average(y,t,imin,imax)
      print('INFO: (plot_zf) Integral time-average = {:.4f}'.format(ave))

      ave_vec = ave*np.ones(len(t))
      #----------------------------------------------------

      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
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

      return

         
   def plot_xflux(self,xin):

      w      = xin['w']
      moment = xin['moment']
      nscale = xin['nscale']
      ymin   = xin['ymin']
      ymax   = xin['ymax']

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      self.getxflux()
      
      ns = self.n_species
      nl = self.n_global+1
      t  = self.t

      ky  = self.ky
      ave = np.zeros((self.n_n,ns))

      if moment == 'phi':
         moment = 'e'

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
         print('ERROR: (plot_xflux) Invalid moment.')
         sys.exit()


      # Call routine for domain average
      e = 0.2
      self.xfluxave(w,moment,e=e,nscale=nscale)

      # Rescale with density ratio
      if nscale == 1:
         mnorm = '^\mathrm{norm}'
      else:
         mnorm = ''


      #============================================================
      # Otherwise plot
      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$r/L_x$')

      color = ['k','m','b','c','g','r']

      imin,imax=time_index(t,w)
      mpre,mwin = wintxt(imin,imax,t)

      ax.set_title(r'$\mathrm{'+ntag+'} \quad $'+mwin)
    
      a = -np.pi+2*np.pi*np.arange(0.0,1.0,0.001)

      for ispec in range(ns):

         u = specmap(self.mass[ispec],self.z[ispec])

         # Flux curve
         g = np.zeros(len(t))
         g = self.lky_xr[ispec,0] 
         for l in range(1,nl):
            g = g+2*(np.cos(l*a)*self.lky_xr[ispec,l]-np.sin(l*a)*self.lky_xi[ispec,l])
         ax.plot(a/(2*np.pi),g,color=color[ispec])

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
         ax.set_xticks([-0.5,-0.375,-0.25,-0.125,0,0.125,0.25,0.375,0.5])
         ax.set_xticklabels([r'$-0.5$',
                             r'$-0.375$',
                             r'$-0.25$',
                             r'$-0.125$',
                             r'$0$',
                             r'$0.125$',
                             r'$0.25$',
                             r'$0.375$',
                             r'$0.5$'])

         ax.legend(loc=2)

      fig.tight_layout(pad=0.3)

      
   def plot_kxky_phi(self,xin):

      norm   = xin['norm']
      w      = xin['w']
      field  = xin['field']
      moment = xin['moment']
      theta  = xin['theta']
      spec   = xin['spec']

      x0 = max(abs(self.kx))
      y0 = max(abs(self.ky))

      xl = -x0 ; xr = x0
      
      if xin['kxmin'] != 'auto':
         xl = float(xin['kxmin'])
      if xin['kxmax'] != 'auto':
         xr = float(xin['kxmax'])

      asp = y0/(xr-xl)
      
      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],asp*xin['lx']))
 
      self.getbigfield()
      t = self.getnorm(norm)

      #-----------------------------------------------------------------
      # Note array structure
      # self.phi = np.reshape(data,(2,self.n_radial,self.n_n,nt),'F')

      nx = self.n_radial
      ny = self.n_n

      fplot = np.zeros([nx-1,ny])
      
      # Field data selector
      f,ft = self.kxky_select(theta,field,moment,spec)
      
      imin,imax=time_index(t,w)
      y = np.sum(abs(f[:,:,imin:imax+1]),axis=2)
         
      # Fix (0,0)
      i0 = nx//2-1
      y[i0,0] = y[i0+1,0]

      # Reverse y order, take transpose, take log
      y = np.transpose(np.log(y[:,::-1]))

      ax = fig.add_subplot(111)

      mpre,mwin = wintxt(imin,imax,t)

      ax.set_xlabel(self.kxstr)
      ax.set_ylabel(self.kystr)
      ax.set_title(r'$\mathrm{Log} |'+ft+'| \quad $'+mwin)

      ax.imshow(y,extent=[xl,xr,0,y0],interpolation='none',cmap='plasma')
      print('INFO: (plot_kxky_phi) min={:.2e} max={:.2e}'.format(np.min(y),np.max(y)))

      fig.tight_layout(pad=0.3)

      return
   
      
   def plot_kx_phi(self,xin):

      w      = xin['w']
      field  = xin['field']
      moment = xin['moment']
      theta  = xin['theta']
      ymin   = xin['ymin']
      ymax   = xin['ymax']
      nstr   = xin['nstr']
      bar    = xin['bar']

      t = self.getnorm(xin['norm'])
      
      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      self.getbigfield()
      
      kx = self.kxnorm
      nx = self.n_radial-1 
      
      imin,imax=time_index(t,w)
    
      dk = kx[1]-kx[0]
      x0 = kx[-1]+dk

      ax = fig.add_subplot(1,1,1)

      mpre,mwin = self.wintxt(imin,imax,t,0,0,field)

      ax.set_title(r'$\mathrm{Average~fluctuation~intensity} \quad $'+mwin)
      ax.set_xlabel(self.kxstr)

      f,ft = self.kxky_select(theta,field,moment,0,gbnorm=True)

      if nstr == 'null' or nstr == '+':
         if nstr == '+':
            y = np.sum(abs(f[:,1:,:])**2,axis=1)
         else:
            y = np.sum(abs(f[:,:,:])**2,axis=1)

         ave = time_average(y,t,imin,imax)
         ax.set_ylabel(self.ylabeler(nstr,ft,sq=True,tave=True))
         if bar:
            ax.step(kx+dk/2,ave,color='m')
         else:
            ax.plot(kx,ave,color='m')
           
      else:
         nvec = str2list(nstr)
         print('INFO: (plot_kx_phi) n = '+str(nvec))
         ax.set_ylabel(self.ylabeler('n',ft,sq=True,tave=True))
         for n in nvec:
            num = r'$n='+str(n)+'$'
            ave = time_average(abs(f[:,n,:])**2,t,imin,imax)
            if bar:
               ax.step(kx+dk/2,ave,label=num)
            else:
               if n == 0:
                  ax.plot(kx[:nx//2:2],ave[:nx//2:2],color='m',alpha=0.4)
                  ax.plot(kx[nx//2+1::2],ave[nx//2+1::2],color='m',alpha=0.4)
                  ax.plot(kx[:nx//2],ave[:nx//2],label=num,color='k')
                  ax.plot(kx[nx//2+1:],ave[nx//2+1:],color='k')
               else:
                  ax.plot(kx,ave,label=num)
                  
            if self.n_n > 16:
               ax.legend(loc=4, ncol=5, prop={'size':12})
            else:
               ax.legend(loc=4, ncol=6, prop={'size':12})

      ax.set_xlim([-x0,x0])
      ax.set_yscale('log')

      if ymin == 'auto':
         ax.set_ylim(bottom=0.5*ave[-1])
      else:
         ax.set_ylim(bottom=float(ymin))
      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))

      if xin['kxmin'] != 'auto':
         ax.set_xlim(left=float(xin['kxmin']))
      if xin['kxmax'] != 'auto':
         ax.set_xlim(right=float(xin['kxmax']))

      fig.tight_layout(pad=0.3)
      
      return

   def plot_kx_shift(self,xin):

      w      = xin['w']
      field  = xin['field']
      moment = xin['moment']
      theta  = xin['theta']
      nstr   = xin['nstr']

      t = self.getnorm(xin['norm'])
      
      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      self.getbigfield()

      kx = self.kxnorm
      nx = self.n_radial
      nn = self.n_n
      nt = self.n_time
      
      ave = np.zeros(nx)

      imin,imax=time_index(t,w)
    
      dk = kx[1]-kx[0]
      x0 = kx[-1]+dk

      ax = fig.add_subplot(1,1,1)

      mpre,mwin = self.wintxt(imin,imax,t,0,0,field)

      ax.set_title(r'$\mathrm{Average~fluctuation~intensity} \quad $'+mwin)
      ax.set_xlabel(self.kxstr)

      f,ft = self.kxky_select(theta,field,moment,0,gbnorm=True)

      lst = ['-','--']
      color = ['black','magenta','blue','red','green','purple']
      
      x = kx[::2]
      if nstr == 'null' or nstr == '+':
         nstr = '+'
         nvec = np.arange(1,nn)
         ax.set_ylabel(self.ylabeler('+',ft,sq=True,tave=True))
      else:
         nvec = str2list(nstr)
         ax.set_ylabel(self.ylabeler('n',ft,sq=True,tave=True))

      y = np.zeros([2,len(x)])
      for i,n in enumerate(nvec):
         phi  = np.zeros([nx,nt],dtype=complex)
         for p in range(1,nx):
            phi[p,:] = f[p-1,n,:]

         phi_T = np.fft.ifft(np.fft.ifftshift(phi,axes=0),axis=0)

         for k in range(2):
            if k == 1:
               phi_T = np.roll(phi_T,nx//2,axis=0)
            
            fint = phi_T[nx//4:3*nx//4,:]
            phi = np.fft.fftshift(np.fft.fft(fint,axis=0))

            y[k,:] = y[k,:]+time_average(np.abs(phi[:,:])**2,self.t,imin,imax)

         if nstr != '+':
            mycol = color[np.mod(i,len(color))]
            num = r'$n='+str(n)+'$'
            ax.plot(x,y[0,:],linestyle='-',color=mycol,label=num)
            ax.plot(x,y[1,:],linestyle='--',color=mycol)
            y[:,:] = 0.0
            ax.legend(loc=4, ncol=6, prop={'size':12})
         elif nstr == '+' and n == nvec[-1]:
            ax.plot(x,y[0,:],linestyle='-',color='k')
            ax.plot(x,y[1,:],linestyle='--',color='k')

      ax.set_yscale('log')
      ax.set_xlim([-x0,x0])
      
      if xin['kxmin'] != 'auto':
         ax.set_xlim(left=float(xin['kxmin']))

      if xin['kxmax'] != 'auto':
         ax.set_xlim(right=float(xin['kxmax']))

      fig.tight_layout(pad=0.3)
      
      return
      

   def plot_poly_phi(self,xin):

      import scipy.special as sp

      w      = xin['w']
      field  = xin['field']
      moment = xin['moment']
      theta  = xin['theta']
      nstr   = xin['nstr']
      ymin   = xin['ymin']
      ymax   = xin['ymax']

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      self.getbigfield()
      nx = self.n_radial
      nt = self.n_time
      nn = self.n_n
      n0 = nx//2
      nk = 2*n0
      t = self.t

      imin,imax=time_index(self.t,w)
      mpre,mwin = self.wintxt(imin,imax,t,0,0,field)

      ax = fig.add_subplot(1,1,1)

      color = ['m','k','b','c']
      xlabel=r'$k$'
   
      ax.set_title(r'$\mathrm{Average~fluctuation~intensity~(Legendre)} \quad $'+mwin)
      ax.set_xlabel(xlabel)

      y1 = np.zeros([nn])
      y2 = np.zeros([nn])

      f,ft = self.kxky_select(theta,field,'phi',0)

      c1 = np.zeros([nk],dtype=complex)
      d1 = np.zeros([nk])
      c2 = np.zeros([nk],dtype=complex)
      d2 = np.zeros([nk])

      mat1 = np.zeros([nk,n0-1])
      mat2 = np.zeros([nk,n0-1])
      kvec = np.arange(nk)
      pvec = np.arange(1,n0)
      z = pvec*np.pi/2
      for k in kvec:
         mat1[k,:]  = sp.spherical_jn(k,z)
         mat2[k,:]  = mat1[k,:]*(-1)**pvec[:]

      ai = 1j**kvec   
      ak = 2*kvec+1   

      if nstr == 'null':
         nvec = np.arange(nn)
         nsum =  True
      else:
         nvec = str2list(nstr)
         print('INFO: (plot_kx_phi) n = '+str(nvec))
         nsum = False
         
      for n in nvec:
         for jt in np.arange(imin,imax+1):

            y = f[:,n,jt]

            phim = np.flip(y[0:n0-1])
            phi0 = y[n0-1]
            phip = y[n0:]

            mp1 = np.matmul(mat1,phip)
            mm1 = np.matmul(mat1,phim)
            mp2 = np.matmul(mat2,phip)
            mm2 = np.matmul(mat2,phim)

            c1[:] = ak[:]*(mp1[:]*ai[:]+mm1[:]*np.conj(ai[:]))
            c2[:] = ak[:]*(mp2[:]*ai[:]+mm2[:]*np.conj(ai[:]))
            
            c1[0] = c1[0]+phi0
            c2[0] = c2[0]+phi0
            
            d1[:] = d1[:]+(np.abs(c1[:]))**2/ak[:]
            d2[:] = d2[:]+(np.abs(c2[:]))**2/ak[:]

         if not nsum:
            ax.plot(kvec,d1)
            ax.plot(kvec,d2)
            # EAB print
            #for i in range(nk):
            #   print(kvec[i],d1[i],d2[i],n)
            c1[:] = 0.0
            c2[:] = 0.0
            d1[:] = 0.0
            d2[:] = 0.0

            
      if nsum:
         ax.plot(kvec,d1,color='k')
         ax.plot(kvec,d2,color='m')
         #for i in range(nk):
         #   print(kvec[i],d1[i],d2[i])

      ax.axvspan(1.5*n0,nk,alpha=0.2)
         
      ax.set_xlim([0,nk])
      ax.set_xlim([0,nk])
      ax.set_yscale('log')
 
      if ymin != 'auto':
         ax.set_ylim(bottom=float(ymin))
      if ymax != 'auto':
         ax.set_ylim(top=float(ymax))
  
      fig.tight_layout(pad=0.3)

      return

   def plot_hb(self,xin):
            
      import matplotlib.cm as cm

      spec  = xin['spec']
      itime = xin['itime']
      tmax  = xin['tmax']
      mesh  = xin['mesh']

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))

      u = specmap(self.mass[spec],self.z[spec])

      if itime > self.n_time-1:
         itime = self.n_time-1

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
         if row == 1:
            ie = self.n_energy//2
         if row == 2:
            ie = self.n_energy-1

         #======================================
         ax = fig.add_subplot(3,2,p)

         ax.set_title(r'${\rm Re} \, h_'+u+' \quad \mathrm{ie}='+str(ie)+'$')
         ax.set_xlabel(r'$\theta/\pi$')
         ax.set_ylabel(r'$\xi = v_\parallel/v$')

         hp = np.transpose(np.array(self.hb[0,:,spec,:,ie,itime]))
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

      return

   def plot_hbcut(self,xin):

      w      = xin['w']
      itime  = xin['itime']
      tmax   = xin['tmax']
      theta  = xin['theta']
      nstr   = xin['nstr']
      ymin   = xin['ymin']
      ymax   = xin['ymax']
      spec   = xin['spec']

      u = specmap(self.mass[spec],self.z[spec])

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))
       
      if itime > self.n_time-1:
         itime = self.n_time-1

      func = self.hb

      # Compute index for theta value in pitch angle and energy plots
      i0,thetapi = indx_theta(theta,self.n_theta)
       
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
            ix = 2
         if row == 1:
            ie = self.n_energy//2
            ix = 3
         if row == 2:
            ie = self.n_energy-1
            ix = 5

         #========================================================
         ax = fig.add_subplot(3,3,p)
         ax.grid(which="both",ls=":")
         ax.grid(which="major",ls=":")

         ax.set_title(r'$'+u+': \\xi=0 \quad \mathrm{ie}='+str(ie)+'$')
         ax.set_xlabel(r'$\theta/\pi$')

         if self.n_xi%2 == 0:
            hp = np.array(func[:,:,spec,self.n_xi//2,ie,itime]+
                          func[:,:,spec,self.n_xi//2-1,ie,itime])*0.5
         else:
            hp = np.array(func[:,:,spec,self.n_xi//2,ie,itime])
            
         ax.plot(x,hp[0,:],'-o',color='black',markersize=2)
         ax.plot(x,hp[1,:],'-o',color='blue',markersize=2)
         ax.set_xlim([-tmax,tmax])

         #========================================================

         p = p+1

         #========================================================
         ax = fig.add_subplot(3,3,p)
         ax.grid(which="both",ls=":")
         ax.grid(which="major",ls=":")

         ax.set_title(r'$'+u+': \\theta/\pi='+str(thetapi)+' \quad \mathrm{ie}='+str(ie)+'$')
         ax.set_xlabel(r'$\xi = v_\parallel/v$')

         n0 = (self.n_radial//2)*self.n_theta+i0
         
         hp = np.array(func[0,:,spec,:,ie,itime])
         ax.plot(self.xi,hp[n0,:],'-o',color='black',markersize=2)
         hp = np.array(func[1,:,spec,:,ie,itime])
         ax.plot(self.xi,hp[n0,:],'-o',color='blue',markersize=2)
         ax.set_xlim([-1,1])
         #========================================================

         p = p+1

         #========================================================
         ax = fig.add_subplot(3,3,p)
         ax.grid(which="both",ls=":")
         ax.grid(which="major",ls=":")

         ax.set_title(r'$'+u+': \\theta/\pi='+str(thetapi)+' \quad \mathrm{ix}='+str(ix)+'$')
         ax.set_xlabel(r'$x=\sqrt{\varepsilon}$')

         n0 = (self.n_radial//2)*self.n_theta+i0

         hpa = np.array(func[0,n0,spec,:,:,itime])
         hpb = np.array(func[1,n0,spec,:,:,itime])
         p0 = np.zeros(self.n_energy)
         p1 = np.zeros(self.n_energy)
                        
         ax.plot(np.sqrt(self.energy),hpa[ix,:],'-o',color='black',markersize=2)
         ax.plot(np.sqrt(self.energy),hpb[ix,:],'-o',color='blue',markersize=2)
         #========================================================

      fig.tight_layout(pad=0.3)

      return

   def plot_hball(self,xin):

      w      = xin['w']
      tmax   = xin['tmax']
      itime  = xin['itime']
      nstr   = xin['nstr']
      ymin   = xin['ymin']
      ymax   = xin['ymax']
      spec   = xin['spec']
      ie     = xin['ie']

      if nstr == 'null':
         nvec = list(range(self.n_n))
      else:
         nvec = str2list(nstr)

      u = specmap(self.mass[spec],self.z[spec])

      # Diagnostics
      print('l    = '+nstr)
      print('e    = '+str(ie))
      print('spec = '+u)

      if xin['fig'] is None:
         fig = plt.figure(MYDIR,figsize=(xin['lx'],xin['ly']))
       
      if itime > self.n_time-1:
         itime = self.n_time-1

      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
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

      return
