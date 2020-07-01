from . import data
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from ..gacodefuncs import *
from .data import GYROData

MYDIR=os.path.basename(os.getcwd())

class gyrodata_plot(data.GYROData):

   def plot_freq(self,w=0.5,wmax=0.0,fig=None):
      '''
      Plot gamma and omega vs time

      ARGUMENTS:
      w: fractional width of time window
      '''

      if self.profile['nonlinear_flag'] == 1:
         raise IOError("ERROR (data_plot): Not available for nonlinear run")

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx*1.2,self.ly))

      t = self.t['(c_s/a)t']

      # Determine tmin
      imin,imax = iwindow(t,w,wmax)

      color = ['k','m','b','c']
      tor_n = self.profile['n0'] + \
              self.profile['d_n']*np.arange(0,self.profile['n_n'])

      #======================================
      ax = fig.add_subplot(121)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$(a/c_s)\gamma$',color='k')
      #=====================================

      print(self.freq['(a/c_s)gamma'])

      # Gamma
      for i in range(self.profile['n_n']):
          ax.plot(t[imin:],self.freq['(a/c_s)gamma'][i,imin:],color=color[i],
                  label='gamma_n%d'%tor_n[i])

      #======================================
      ax = fig.add_subplot(122)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$(a/c_s)\omega$',color='k')
      #=====================================

      # Omega
      for i in range(self.profile['n_n']):
          ax.plot(t[imin:],self.freq['(a/c_s)w'][i,imin:],color=color[i],
                  label='omega_n%d'%tor_n[i])

      plt.tight_layout(pad=0.2)

   def plot_balloon(self,index='phi',tmax=-1,fig=None):
      '''
      Plot the ballooning eigenmode structure
       
      ARGUMENTS:
      index: one of 'phi','a','aperp' or an integer
      tmax : max value of (extended angle)/pi
      '''

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      self.read_balloon()

      if index == 'phi':
         key = 'balloon_phi'
         index = list(self.balloon.keys()).index(key)
      elif index == 'a':
         key = 'balloon_a'
         index = list(self.balloon.keys()).index(key)
      elif index == 'aperp':
         key = 'balloon_aperp'
         index = list(self.balloon.keys()).index(key)
      else:
         index = int(index)

      key = list(self.balloon)[index]

      ytitle = '\\begin{verbatim}'+key+'\\end{verbatim}'

      if key == 'balloon_a':
         ytitle = r'$\delta A_\parallel$'
      if key == 'balloon_phi':
         ytitle = r'$\delta \phi$'
      if key == 'balloon_aperp':
         ytitle = r'$\delta B_\parallel$'
      if key == 'balloon_epar':
         ytitle = r'$\delta E_\parallel$'

      #======================================
      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$\theta_*/\pi$')
      ax.set_ylabel(ytitle)
      #=====================================

      n_p   = self.profile['n_x']/self.profile['box_multiplier']
      n_ang = self.profile['n_theta_plot']*n_p

      x = -(1.0+n_p)+2.0*n_p*np.arange(n_ang)/float(n_ang)

      ax.plot(x,np.real(self.balloon[key][:,0,-1]),color='k',label='Re')
      ax.plot(x,np.imag(self.balloon[key][:,0,-1]),color='m',label='Im')

      if tmax < 0.0:
         ax.set_xlim([-(n_p-1),n_p-1])
      else:
         ax.set_xlim([-tmax,tmax])

      ax.legend()
      plt.tight_layout(pad=0.2)
        
      return key

   def plot_zf(self,w=0.5,fig=None):
      '''
      Plot the zonal (n=0) potential versus time.

      ARGUMENTS:
      w: fractional time-average window
      '''

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      # Read freq data
      self.read_moment('u')
   
      ntheta = self.profile['n_theta_plot']
      nx = self.profile['n_x']

      y = np.real(self.moment_u[ntheta/3,nx/3,0,0,:])

      if abs(y[0]) > 1e-6:
         y = y[:]/y[0]
         ylabel = r'$\delta\phi/\delta\phi(0)$'
         test_flag = 1
      else:
         ylabel = r'$\delta\phi$'
         test_flag = 0
         
         
      t = self.t['(c_s/a)t']

      #----------------------------------------------------
      # Average calculations

      imin,imax = iwindow(t,w,wmax)
      ave  = average(y[:],t,w,wmax)
      print('INFO: (plot_zf) Spatial point (nx,ntheta)=',nx,ntheta)

      ave_vec = ave*np.ones(len(t))
      #----------------------------------------------------

      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(ylabel)

      ax.plot(t,y,color='k')

      if test_flag:
         theory = 1.0/(1.0+1.6*2.0**2/np.sqrt(0.5/3.0))
         ax.plot([0,max(t)],[theory,theory],color='grey',
                 label=r'$\mathrm{RH \; theory}$',alpha=0.3,linewidth=4)
         ax.plot(t[imin:],ave_vec[imin:],color='b',
                 label=r'$\mathrm{Average}$',linewidth=1)
         print('INFO: (plot_zf) Integral time-average = %.6f' % ave)

      ax.legend()
      plt.tight_layout(pad=0.2)

   def plot_phi(self,ymax='auto',fig=None):
      '''
      Plot the n=0 AND n>0 potentials versus time.

      ARGUMENTS:        
      ymax : max vertical plot range
      '''

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx*1.2,self.ly))
         
      t = self.t['(c_s/a)t']

      #======================================
      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      ax.set_ylabel(r'$\langle e \phi/T_e \rangle/\rho_\star $',color='k')
      #=====================================

      # rho_* for normalization
      x=np.average(self.profile['rho_s'])

      ax.plot(t,self.field_rms[0,:]/x,color='k',label=r'$n=0$')
      ax.plot(t,self.field_rms[1,:]/x,color='purple',label=r'$n>0$')

      ax.set_xlim([0,t[-1]])

      if ymax != 'auto':
         ax.set_ylim([0,float(ymax)])

      ax.legend()
      plt.tight_layout(pad=0.2)

   def plot_gbflux(self,w=0.5,wmax=0.0,field='s',moment=0,ymin='0.0',ymax='auto',loc=2, fig=None):
      '''
      Plot the gyrobohm flux versus time.
      '''

      self.read_gbflux_i()

      if fig is None:
         fig = plt.figure(MYDIR, figsize=(self.lx, self.ly))


      ns = int(self.profile['n_kinetic'])

      t = self.t['(c_s/a)t']

      field_tag = '\mathrm{Total}'

      # Manage field
      if field == 's':
         flux0 = np.sum(self.gbflux,axis=1)
      else:
         i_field = int(field)
         flux0   = self.gbflux[:,i_field,:,:]
         if field == 0:
            field_tag = '\phi'
         elif field == 1:
            field_tag = 'A_\parallel'
         else:
            field_tag = 'B_\parallel'

      if moment == 0:
         ntag = 'Density~flux'
         mtag = '\Gamma'
         ttag = 'G'
         ftag = 'flux_n'
         y = flux0[:,0,:]
      elif moment == 1:
         ntag = 'Energy~flux'
         mtag = 'Q'
         ttag = 'Q'
         ftag = 'flux_e'
         y = flux0[:,1,:]
      elif moment == 2:
         ntag = 'Momentum~flux'
         mtag = '\Pi'
         ttag = 'Pi'
         ftag = 'flux_v'
         y = flux0[:,2,:]
      else:
         print('ERROR: (plot_flux.py) Invalid moment.')
         sys.exit()
       
      # Normalizations
      nscale = 0
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
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)

      color = ['k','m','b','c','g','r']

      windowtxt = '['+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+']'

      ax.set_title(r'$\mathrm{'+ntag+'} \quad '+windowtxt+'\quad ['+field_tag+']$')

      for ispec in range(ns):
         y_norm = y[ispec,:]*norm_vec[ispec]
         ave    = average(y_norm,t,w,wmax)
         y_ave  = ave*np.ones(len(t))
         u = specmap(1.0/self.profile['mu'][ispec]**2,self.profile['z'][ispec])
         label = r'$'+mtag+mnorm+'_'+u+'/'+mtag+'_\mathrm{GB}: '+str(round(ave,3))+'$'
         # Average
         ax.plot(t[imin:],y_ave[imin:],'--',color=color[ispec])
         # Time trace
         ax.plot(t,y_norm,label=label,color=color[ispec])

      ax.legend(loc=loc)

      if ymax != 'auto':
         ax.set_ylim([float(ymin),float(ymax)])

      fig.tight_layout()

   def plot_gbflux_i(self,w=0.5,aw=0,field='s',moment=0,ymin='0.0',ymax='auto',loc=2,fig=None):
      '''
      Plot flux versus radius
      '''

      # Read data in gbflux_i
      self.read_gbflux_i()

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      ns   = int(self.profile['n_kinetic'])
      n_x  = self.profile['n_x']
      nd   = self.profile['n_explicit_damp']

      r = self.profile['r']
      t = self.t['(c_s/a)t']

      field_tag = '\mathrm{Total}'

      fluxi = self.gbflux_i
      flux  = self.gbflux

      # Manage field
      if field == 's':
         fluxi0 = np.sum(self.gbflux_i,axis=1)
         flux0  = np.sum(self.gbflux,axis=1)
      else:
         i_field = int(field)
         fluxi0 = self.gbflux_i[:,i_field,:,:,:]
         flux0  = self.gbflux[:,i_field,:,:]
         if field == 0:
            field_tag = '\phi'
         elif field == 1:
            field_tag = 'A_\parallel'
         else:
            field_tag = 'B_\parallel'
         
      if moment == 0:
         ntag = 'Density~flux'
         mtag = '\Gamma'
         ttag = 'G'
         ftag = 'flux_n'
         yi = fluxi0[:,0,:,:]
         y  = flux0[:,0,:]
      elif moment == 1:
         ntag = 'Energy~flux'
         mtag = 'Q'
         ttag = 'Q'
         ftag = 'flux_e'
         yi = fluxi0[:,1,:,:]
         y  = flux0[:,1,:]
      elif moment == 2:
         ntag = 'Momentum~flux'
         mtag = '\Pi'
         ttag = 'Pi'
         ftag = 'flux_v'
         yi = fluxi0[:,2,:,:]
         y  = flux0[:,2,:]
      else:
         print('ERROR: (plot_flux.py) Invalid moment.')
         sys.exit()

      imin,imax = iwindow(t,w,wmax)

      # Normalizations
      nscale = 0
      if nscale == 0:
         norm_vec = np.ones(ns)
         mnorm = ''
      else:
         norm_vec = 1.0/self.dens
         mnorm = '^\mathrm{norm}'

      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$r/a$')

      color = ['k','m','b','c','g','r']

      windowtxt = '['+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+']'

      ax.set_title(r'$\mathrm{'+ntag+'} \quad '+windowtxt+'\quad ['+field_tag+']$')

      avei = np.zeros(n_x)
      one  = np.ones(n_x)
      
      # Loop over species
      for ispec in range(ns):
         u = specmap(1.0/self.profile['mu'][ispec]**2,self.profile['z'][ispec])

         # Time-averaged flux curve
         avei[:] = average_n(yi[ispec,:,:]*norm_vec[ispec],t,w,wmax,n_x)
         ax.plot(r,avei[:],color=color[ispec])

         # Full average
         ave = average(y[ispec,:]*norm_vec[ispec],t,w,wmax)
         print('INFO: (plot_gbflux_i) Full average = {:.2f}'.format(ave))
         # Partial-r average
         if aw > 0:
            ave = np.average(avei[aw:-(aw+1)])
         else:
            ave = np.average(avei[nd:-(nd+1)])
         label = r'$'+mtag+mnorm+'_'+u+'/'+mtag+'_\mathrm{GB}: '+str(round(ave,3))+'$'

         if aw > 0:
            ax.plot(r[aw:-(aw+1)],ave*one[aw:-(aw+1)],'--',label=label,color=color[ispec])
         else:
            if nd > 0:
               ax.plot(r[nd:-(nd+1)],ave*one[nd:-(nd+1)],'--',label=label,color=color[ispec])
            else:
               ax.plot(r,ave*one,'--',label=label,color=color[ispec])

      if nd > 0:
         ax.axvspan(r[0],r[nd],facecolor='g',alpha=0.1)
         ax.axvspan(r[-(nd+1)],r[-1],facecolor='g',alpha=0.1)
      ax.set_xlim(r[0],r[-1])

      if ymax != 'auto':
         ax.set_ylim([float(ymin),float(ymax)])

      ax.legend()
      plt.tight_layout()

   def plot_gbflux_n(self,w=0.5,wmax=0.0,field='s',moment=0,ymin='auto',ymax='auto',fig=None):
      '''
      Plot ky-dependent flux
      '''
      ns      = int(self.profile['n_kinetic'])
      n_n     = int(self.profile['n_n'])

      if n_n == 1:
         print('ERROR (plot_gbflux_n.py) Plot not available with a single mode.')
         sys.exit()

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.ly*ns,self.ly))

      # Need to read gbflux_n data
      self.read_gbflux_n()

      t    = self.t['(c_s/a)t']
      flux = self.gbflux_n

      if field == 's':
         ys = np.sum(self.gbflux_n,axis=1)
      else:
         ys = self.gbflux_n[:,field,:,:,:]
         if field == 0:
            field_tag = '\phi'
         elif field == 1:
            field_tag = 'A_\parallel'
         else:
            field_tag = 'B_\parallel'

      if moment == 0:
         ntag = 'Density~flux'
         mtag = '\Gamma'
         ttag = 'G'
         ftag = 'flux_n'
         y = ys[:,0,:,:]
      elif moment == 1:
         ntag = 'Energy~flux'
         mtag = 'Q'
         ttag = 'Q'
         ftag = 'flux_e'
         y = ys[:,1,:,:]
      elif moment == 2:
         ntag = 'Momentum~flux'
         mtag = '\Pi'
         ttag = 'Pi'
         ftag = 'flux_v'
         y = ys[:,2,:,:]
      else:
         print('ERROR (plot_ky_flux.py) Invalid moment.')
         sys.exit()

      imin,imax = iwindow(t,w,wmax)

      color = ['k','m','b','c','g','r']

      windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+']$'

      ky  = self.profile['kt_rho']
      dk  = ky[1]-ky[0]
      ave = np.zeros((n_n,ns))

      for ispec in range(ns):
         for j in range(n_n):
            ave[j,ispec] = average(y[ispec,j,:],t,w,wmax)

      for ispec in range(ns):
         ax = fig.add_subplot(1,ns,ispec+1)
         ax.set_xlabel(r'$k_\theta \rho_s$')
         u = specmap(1.0/self.profile['mu'][ispec]**2,self.profile['z'][ispec])
         ax.set_ylabel(r'$'+mtag+'_'+u+'$',color='k')
         ax.set_title(windowtxt)
         ax.bar(ky-dk/2.0,ave[:,ispec],width=dk/1.1,color=color[ispec],
                alpha=0.5,edgecolor='black')
         
      # Set axis ranges
      ax.set_xlim([0,ky[-1]+dk])
      if ymax != 'auto':
         ax.set_ylim([0,float(ymax)])

      plt.tight_layout()

   def plot_gbflux_exc(self,w=0.5,fig=None):

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(12,8))

      ns = int(self.profile['n_kinetic'])
      t = self.t['(c_s/a)t']

      # Read data in gbflux_i and make gbflux
      self.read_gbflux_exc()

      flux = self.gbflux_exc

      #======================================
      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(TIME)
      #=====================================

      imin,imax=iwindow(t,w,wmax)

      color = ['k','m','b','c']

      # Loop over species
      for i in range(ns):
         for j in range(2):
            ave   = average(flux[i,j,:],t,w,wmax)
            stag  = self.tagspec[i]
            label = stag+': '+str(round(ave,3))
            y     = ave*np.ones(len(t))
            ax.plot(t[imin:],y[imin:],'--',color=color[i])
            ax.plot(t,flux[i,j,:],label=label,color=color[i])

      ax.legend()
      plt.tight_layout()

   def plot_gbflux_rt(self,w=0.5,wmax=0.0,field='s',moment=0,fig=None):

      ns = int(self.profile['n_kinetic'])

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(7*ns,6))

      t = self.t['(c_s/a)t']

      self.read_gbflux_i()

      flux = self.gbflux_i

      # Manage field
      if field == 's':
         flux0 = np.sum(flux,axis=1)
         ftag  = self.tagfield[3]
      else:
         field = int(field)
         flux0 = flux[:,field,:,:,:]
         ftag  = self.tagfield[field]

      # Manage moment
      mtag = self.tagmom[moment]

      imin,imax = iwindow(t,w,wmax)

      color = ['k','m','b','c']

      # Loop over species
      for i in range(ns):
           stag  = self.tagspec[i]
           ax = fig.add_subplot(1,ns,i+1)
           ax.set_xlabel(r'$(c_s/a)t$')
           ax.set_ylabel(r'$r/a$')
           ax.set_title(r'$'+mtag+' \;('+ftag+'\,\mathrm{'+stag+'})$',color='k')
           ax.contourf(t[imin:],self.profile['r'],flux0[i,moment,:,imin:],
                       cmap=plt.cm.jet,nlevels=128)
           ax.set_xlim([t[imin],t[-1]])
           ax.set_ylim([self.profile['r'][0],self.profile['r'][-1]])

      plt.tight_layout()

   def plot_moment_zero(self,w=0.5,species=0,moment=0,fig=None):

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      nt  = self.n
      n_x = self.profile['n_x']
      nd  = self.profile['n_explicit_damp']
      
      t = self.t['(c_s/a)t']
      r = self.profile['r']
      
      # Get n=0 moment data
      self.read_moment_zero()

      if moment == 0:
         ntag = 'Density~moment'
         mtag = '\delta n'
         ftag = 'delta_n'
      elif moment == 1:
         ntag = 'Energy~moment'
         mtag = '\delta E'
         ftag = 'delta_e'
      elif moment == 2:
         ntag = 'vpar~moment'
         mtag = '\delta v'
         ftag = 'delta_v'

      imin,imax  = iwindow(t,w,wmax)
      title = r'$'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$'
      
      u = specmap(1.0/self.profile['mu'][species]**2,self.profile['z'][species])
      
      f = np.zeros([n_x,nt])
      g = np.zeros(n_x)

      # Compare moment versus corresponding source
      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$r/a$')
      ax.set_ylabel(r'$'+mtag+'_'+u+'/\\rho_s$')
      ax.set_title(title)

      f = self.moment_zero[:,species,moment,:]
      g = average_n(f,t,w,wmax,n_x)/self.profile['rho_s']

      ax.plot(r,g,label=r'moment',color='k')

      f = self.moment_zero[:,species,moment+3,:]
      g = average_n(f,t,w,wmax,n_x)/self.profile['rho_s']
      
      ax.plot(r,g,label=r'source',color='m',linewidth=3,alpha=0.2)

      if nd > 0:
         ax.axvspan(r[0],r[nd],facecolor='g',alpha=0.1)
         ax.axvspan(r[-(nd+1)],r[-1],facecolor='g',alpha=0.1)
      ax.set_xlim(r[0],r[-1])
      
      ax.legend()
      plt.tight_layout()

   def plot_profile_tot(self,w=0.5,species=0,moment=0,fig=None):

      if fig is None:
         fig = plt.figure(MYDIR,figsize=(self.lx,self.ly))

      nt  = self.n
      n_x = self.profile['n_x']
      nd  = self.profile['n_explicit_damp']
      
      t = self.t['(c_s/a)t']
      r = self.profile['r']
      
      # Get n=0 moment data
      self.read_moment_zero()

      n0 = self.profile['den_s'][species,:]
      T0 = self.profile['tem_s'][species,:]

      f = np.zeros([n_x,nt])
      f = self.moment_zero[:,species,0,:]
      dn = average_n(f,t,w,wmax,n_x)
       
      f  = self.moment_zero[:,species,1,:]
      dE = average_n(f,t,w,wmax,n_x)
      dT = ((2.0/3)*dE-T0*dn)/n0

      # Derivatives
      dndr = np.gradient(dn,edge_order=2)/np.gradient(r,edge_order=2)
      dTdr = np.gradient(dT,edge_order=2)/np.gradient(r,edge_order=2)

      if moment == 0:
         ntag = 'Density'
         mtag = 'a/L_n'
         ftag = 'dlnndr'
         g0 = self.profile['dlnndr_s'][species,:]
         g =  g0-dndr/n0
      elif moment == 1:
         ntag = 'Temperature'
         mtag = 'a/L_T'
         ftag = 'dlntdr'
         g0 =  self.profile['dlntdr_s'][species,:]
         g  =  g0-dTdr/T0
      else:
         sys.exit()
         
      imin,imax  = iwindow(t,w,wmax)
      title = r'$'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$'
      
      u = specmap(1.0/self.profile['mu'][species]**2,self.profile['z'][species])
      
      # Compare moment versus corresponding source
      ax = fig.add_subplot(111)
      ax.grid(which="both",ls=":")
      ax.grid(which="major",ls=":")
      ax.set_xlabel(r'$r/a$')
      ax.set_ylabel(r'$'+mtag+'_'+u+'$')
      ax.set_title(title)

      ax.plot(r,g0,color='m',linewidth=3,alpha=0.2)
      ax.plot(r,g,color='k')

      if nd > 0:
         ax.axvspan(r[0],r[nd],facecolor='g',alpha=0.1)
         ax.axvspan(r[-(nd+1)],r[-1],facecolor='g',alpha=0.1)
      ax.set_xlim(r[0],r[-1])
      
      plt.tight_layout()
    
      
