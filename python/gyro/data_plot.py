import data
import sys
import numpy as np
import matplotlib.pyplot as plt
from gacodefuncs import *
from gyro.data import GYROData

class gyrodata_plot(data.GYROData):

    def plot_freq(self,w=0.5,fig=None):
        '''
        Plot gamma and omega vs time

        ARGUMENTS:
         w: fractional width of time window
        '''

        if fig is None:
            fig = plt.figure(figsize=(12,6))
            fig.subplots_adjust(left=0.085,right=0.97,top=0.92,bottom=0.12)

        t = self.t['(c_s/a)t']

        # Read freq data
        self.read_freq()

        # Determine tmin
        imin = iwindow(t,w)

        color = ['k','m','b','c']
        tor_n = self.profile['n0'] + \
                self.profile['d_n']*np.arange(0,self.profile['n_n'])
        #======================================
        ax = fig.add_subplot(121)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(TIME)
        ax.set_ylabel(r'$(a/c_s)\gamma$',color='k')
        #=====================================

        # Gamma
        for i in range(self.profile['n_n']):
            ax.plot(t[imin:],self.freq['(a/c_s)gamma'][i,imin:],color=color[i],
                    label='gamma_n%d'%tor_n[i])

        #======================================
        ax = fig.add_subplot(122)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(TIME)
        ax.set_ylabel(r'$(a/c_s)\omega$',color='k')
        #=====================================

        # Omega
        for i in range(self.profile['n_n']):
            ax.plot(t[imin:],self.freq['(a/c_s)w'][i,imin:],color=color[i],
                    label='omega_n%d'%tor_n[i])

    def plot_balloon(self,index='phi',tmax=-1,fig=None):
        '''
        Plot the ballooning eigenmode structure
       
        ARGUMENTS:
         index: one of 'phi','a','aperp' or an integer
         tmax : max value of (extended angle)/pi
        '''

        if fig is None:
            fig = plt.figure(figsize=(10,6))

        self.read_balloon()

        if index == 'phi':
            key = 'balloon_phi'
            index = self.balloon.keys().index(key)
        elif index == 'a':
            key = 'balloon_a'
            index = self.balloon.keys().index(key)
        elif index == 'aperp':
            key = 'balloon_aperp'
            index = self.balloon.keys().index(key)
        else:
            index = int(index)

        key = self.balloon.keys()[index]

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
        ax.grid(which="majorminor",ls=":")
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
        return key

    def plot_zf(self,w=0.5,fig=None):
        '''
        Plot the zonal (n=0) potential versus time.

        ARGUMENTS:
         w: fractional time-average window
        '''

        if fig is None:
            fig = plt.figure(figsize=(10,6))
            fig.subplots_adjust(left=0.1,right=0.95,top=0.95,bottom=0.12)

        # Read freq data
        self.read_moment_u()

        ntheta = self.profile['n_theta_plot']
        nx = self.profile['n_x']

        y = np.real(self.moment_u[ntheta/3,nx/3,0,0,:])

        y = y[:]/y[0]
        t = self.t['(c_s/a)t']

        #----------------------------------------------------
        # Average calculations

        imin = iwindow(t,w)
        ave  = average(y[:],t,w)
        print 'INFO: (plot_zf) Integral time-average = %.6f' % ave
        print 'INFO: (plot_zf) (nx,ntheta)=',nx,ntheta

        ave_vec = ave*np.ones(len(t))
        #----------------------------------------------------

        ax = fig.add_subplot(111)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(TIME)
        ax.set_ylabel(r'$\Phi/\Phi_0$')

        ax.plot(t,y,color='k')
        ax.plot(t[imin:],ave_vec[imin:],color='b',
                label=r'$\mathrm{Average}$',linewidth=1)

        theory = 1.0/(1.0+1.6*2.0**2/np.sqrt(0.5/3.0))
        ax.plot([0,max(t)],[theory,theory],color='grey',
                label=r'$\mathrm{RH \; theory}$',alpha=0.3,linewidth=4)

        ax.legend()

    def plot_phi_n0(self,lx=10,ly=6,ymax='auto',span1=-1.0,span2=-1.0):
        '''
        Plot the n=0 AND n>0 potentials versus time.

        ARGUMENTS:        
         lx   : width of figure 
         ly   : height of figure 
         ymax : max vertical plot range
         span1: left end of axvspan
         span2: right end of avxspan
        '''

        fig = plt.figure(figsize=(lx,ly))
        fig.subplots_adjust(left=0.1,right=0.96,top=0.93,bottom=0.13)

        t = self.t['(c_s/a)t']

        # Read data in gbflux_i and make gbflux
        self.read_field_rms()

        #======================================
        ax = fig.add_subplot(111)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(TIME)
        ax.set_ylabel(r'$\langle e \phi/T_e \rangle/\rho_\star $',color='k')
        #=====================================

        # rho_* for normalization
        x=np.average(self.profile['rho_s'])

        ax.plot(t,self.field_rms[0,:]/x,color='k',label=r'$n=0$')
        ax.plot(t,self.field_rms[1,:]/x,color='purple',label=r'$n>0$')

        if span1 > 0.0:
           ax.axvspan(span1,span2,facecolor='g',alpha=0.1)

        ax.set_xlim([0,t[-1]])

        if ymax != 'auto':
           ax.set_ylim([0,float(ymax)])

        ax.legend()

    def plot_gbflux(self,field='s',i_moment=0,w=0.5,lx=12,ly=6,
                    title='',ymin='0.0',ymax='auto',span1=-1.0,span2=-1.0):
        '''
        Plot the n=0 AND n>0 potentials versus time.

        ARGUMENTS:        
         lx   : width of figure 
         ly   : height of figure 
         ymax : max vertical plot range
         span1: left end of axvspan
         span2: right end of avxspan
        '''

        n_field   = int(self.profile['n_field'])
        n_kinetic = int(self.profile['n_kinetic'])

        t = self.t['(c_s/a)t']

        # Read data in gbflux_i and make gbflux
        self.read_gbflux_i()
        self.make_gbflux()

        flux = self.gbflux

        # Manage field
        if field == 's':
            flux0 = np.sum(flux,axis=1)
            ftag  = self.tagfield[3]
        else:
            i_field = int(field)
            flux0   = flux[:,i_field,:,:]
            ftag    = self.tagfield[i_field]

        # Manage moment
        mtag = self.tagmom[i_moment]

        #======================================
        fig = plt.figure(figsize=(lx,ly))
        fig.subplots_adjust(left=0.1,right=0.95,top=0.92,bottom=0.12)
        ax = fig.add_subplot(111)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(TIME)
        ax.set_ylabel(r'$'+mtag+' \;('+ftag+')$',color='k')
        #=====================================

        cvec = ['k','m','b','c','g']

        # Determine tmin
        for i in range(len(t)):
            if t[i] < (1.0-w)*t[len(t)-1]:
                imin = i

        if title=='null':
            ax.set_title(r'$'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$')
        else:
            ax.set_title(r'$\mathrm{'+title+'}$')
    

        # Plot data to screen or image file.
        for i in range(n_kinetic):
            ave   = average(flux0[i,i_moment,:],t,w)
            if i > n_kinetic-2 and self.profile['electron_method'] > 1:
                stag = r'$e'
                color = 'r'
            else:
                stag = r'$i_'+str(i+1)
                color = cvec[i]

            label = stag+' : '+str(round(ave,3))+'$'
            y     = ave*np.ones(len(t))
            ax.plot(t[imin:],y[imin:],'--',color=color)
            ax.plot(t,flux0[i,i_moment,:],label=label,color=color)

        if span1 > 0.0:
            ax.axvspan(span1,span2,facecolor='g',alpha=0.1)

        ax.set_xlim([0,t[-1]])

        if ymax != 'auto':
            ax.set_ylim([float(ymin),float(ymax)])
        
            ax.legend(loc=1)

    def plot_gbflux_i(self,field='s',i_moment=0,w=0.5,ymin='0.0',ymax='auto',fig=None):
        '''
        Plot the n=0 AND n>0 potentials versus time.

        ARGUMENTS:        
         lx   : width of figure 
         ly   : height of figure 
         ymax : max vertical plot range
         span1: left end of axvspan
         span2: right end of avxspan
        '''

        if fig is None:
            fig = plt.figure(figsize=(12,8))

        n_field   = int(self.profile['n_field'])
        n_kinetic = int(self.profile['n_kinetic'])

        t = self.t['(c_s/a)t']

        # Read data in gbflux_i
        self.read_gbflux_i()

        flux = self.gbflux_i

        # Manage field
        if field == 's':
           flux0 = np.sum(flux,axis=1)
           ftag  = self.tagfield[3]
        else:
           i_field = int(field)
           flux0 = flux[:,i_field,:,:,:]
           ftag  = self.tagfield[i_field]

        # Manage moment
        mtag = self.tagmom[i_moment]

        imin = iwindow(t,w)

        #======================================
        ax = fig.add_subplot(111)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(r'$r/a$')
        ax.set_ylabel(r'$'+mtag+' \;('+ftag+')$',color='k')
        ax.set_title(r'$'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$')
        #=====================================

        color = ['k','m','b','c']

        n_x = self.profile['n_x']
        ave = np.zeros(n_x)

        # Loop over species
        for i in range(n_kinetic):
           stag = self.tagspec[i]
           ave[:] = average_n(flux0[i,i_moment,:,:],t,w,n_x)
           ax.plot(self.profile['r'],ave[:],label=stag,color=color[i])

        if ymax != 'auto':
           ax.set_ylim([float(ymin),float(ymax)])

        ax.legend()


    def plot_gbflux_n(self,field='s',i_moment=0,w=0.5,datafile='none',fig=None):
        '''
        Plot the n=0 AND n>0 potentials versus time.

        ARGUMENTS:        
         lx   : width of figure 
         ly   : height of figure 
         ymax : max vertical plot range
         span1: left end of axvspan
         span2: right end of avxspan
        '''

        n_field   = int(self.profile['n_field'])
        n_kinetic = int(self.profile['n_kinetic'])
        n_n       = int(self.profile['n_n'])

        if fig is None:
           fig = plt.figure(figsize=(7*n_kinetic,6))

        # Need to read gbflux_n data
        self.read_gbflux_n()

        t    = self.t['(c_s/a)t']
        flux = self.gbflux_n

        # Manage field
        if field == 's':
           flux0 = np.sum(flux,axis=1)
           ftag  = self.tagfield[3]
        else:
           i_field = int(field)
           flux0 = flux[:,i_field,:,:,:]
           ftag  = self.tagfield[i_field]

        # Manage moment
        mtag = self.tagmom[i_moment]

        color = ['k','m','b','c','g','r']

        k   = self.profile['kt_rho']
        dk  = k[1]-k[0]
        ave = np.zeros((n_n))

        # Determine tmin
        imin=0
        for i in range(len(t)):
           if t[i] < (1.0-w)*t[len(t)-1]:
              imin = i+1

        if datafile == 'none':

           # Plot data to screen or image file.

           for i in range(n_kinetic):
              ax = fig.add_subplot(1,n_kinetic,i+1)
              ax.set_xlabel(r'$k_\theta \rho_s$')
              ax.set_ylabel(r'$'+mtag+' \;('+ftag+')$',color='k')
              stag = self.tagspec[i]
              for j in range(n_n):
                 ave[j] = average(flux0[i,i_moment,j,:],t,w)
                 ax.set_title(stag+r': $'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$')
                 ax.bar(k-dk/2.0,ave,width=dk/1.1,color=color[i],alpha=0.4,edgecolor='black')

        else:

           # Write data to datafile

           arr = np.zeros([len(k),n_kinetic+1])
           arr[:,0] = k
           stag = '# (k_y rho_s'
           for i in range(n_kinetic):
              for j in range(n_n):
                 ave[j] = average(flux0[i,i_moment,j,:],t,w)
              arr[:,i+1] = ave
              stag = stag+' , '+sim.tagspec[i]
            
           fid = open(datafile,'w')
           fid.write('# Moment  : '+mtag+'\n')
           fid.write('# Field   : '+ftag+'\n')
           fid.write('# Time    : '+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'\n')
           fid.write(stag+')\n')
           np.savetxt(fid,arr,fmt='%.5e')
           fid.close()
        
    def plot_gbflux_exc(self,w=0.5,fig=None):
        '''
        Plot the n=0 AND n>0 potentials versus time.

        ARGUMENTS:        
         lx   : width of figure 
         ly   : height of figure 
         ymax : max vertical plot range
         span1: left end of axvspan
         span2: right end of avxspan
        '''

        if fig is None:
            fig = plt.figure(figsize=(12,8))

        n_kinetic = int(self.profile['n_kinetic'])

        t = self.t['(c_s/a)t']

        # Read data in gbflux_i and make gbflux
        self.read_gbflux_exc()

        flux = self.gbflux_exc

        #======================================
        ax = fig.add_subplot(111)
        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")
        ax.set_xlabel(TIME)
        #=====================================

        # Determine tmin
        for i in range(len(t)):
            if t[i] < (1.0-w)*t[len(t)-1]:
                imin = i

        color = ['k','m','b','c']

        # Loop over species
        for i in range(n_kinetic):
            for j in range(2):
                ave   = average(flux[i,j,:],t,w)
                stag  = self.tagspec[i]
                label = stag+': '+str(round(ave,3))
                y     = ave*np.ones(len(t))
                ax.plot(t[imin:],y[imin:],'--',color=color[i])
                ax.plot(t,flux[i,j,:],label=label,color=color[i])

        ax.legend()

        
    def plot_gbflux_rt(self,field='s',i_moment=0,w=0.5,fig=None):
        '''
        Plot the n=0 AND n>0 potentials versus time.

        ARGUMENTS:        
         lx   : width of figure 
         ly   : height of figure 
         ymax : max vertical plot range
         span1: left end of axvspan
         span2: right end of avxspan
        '''

        n_field   = int(self.profile['n_field'])
        n_kinetic = int(self.profile['n_kinetic'])

        if fig is None:
           fig = plt.figure(figsize=(7*n_kinetic,6))

        t = self.t['(c_s/a)t']

        # Read data in gbflux_i and make gbflux
        self.read_gbflux_i()
        self.make_gbflux()

        flux = self.gbflux_i

        # Manage field
        if field == 's':
           flux0 = np.sum(flux,axis=1)
           ftag  = self.tagfield[3]
        else:
           i_field = int(field)
           flux0 = flux[:,i_field,:,:,:]
           ftag  = self.tagfield[i_field]

        # Manage moment
        mtag = self.tagmom[i_moment]

        imin = iwindow(t,w)

        color = ['k','m','b','c']

        # Loop over species
        for i in range(n_kinetic):
           stag  = self.tagspec[i]
           ax = fig.add_subplot(1,n_kinetic,i+1)
           ax.set_xlabel(r'$(c_s/a)t$')
           ax.set_ylabel(r'$r/a$')
           ax.set_title(r'$'+mtag+' \;('+ftag+'\,\mathrm{'+stag+'})$',color='k')
           ax.contourf(t[imin:],self.profile['r'],flux0[i,i_moment,:,imin:],cmap=plt.cm.jet,nlevels=128)
           ax.set_xlim([t[imin],t[-1]])
           ax.set_ylim([self.profile['r'][0],self.profile['r'][-1]])


    def plot_moment_zero(self,i_moment=0,w=0.5,fig=None):
        '''
        Plot the n=0 AND n>0 potentials versus time.

        ARGUMENTS:        
         lx   : width of figure 
         ly   : height of figure 
         ymax : max vertical plot range
         span1: left end of axvspan
         span2: right end of avxspan
        '''

        if fig is None:
           fig = plt.figure(figsize=(7*3,6))

        n_field   = self.profile['n_field']
        n_kinetic = self.profile['n_kinetic']
        n_x       = self.profile['n_x']

        t = self.t['(c_s/a)t']

        # Read data in gbflux_i
        self.read_moment_zero()

        delta_n = np.empty([n_x,n_kinetic])
        delta_e = np.empty([n_x,n_kinetic])
        delta_t = np.empty([n_x,n_kinetic])

        # [n_x,n_kinetic,n_time]
        for i in range(n_kinetic):
           f   = np.array(self.moment_zero[:,i,0,:])
           delta_n[:,i] = average_n(f,t,w,n_x)
           f   = np.array(self.moment_zero[:,i,1,:])
           delta_e[:,i] = average_n(f,t,w,n_x)
           delta_t[:,i] = (delta_e[:,i]-1.5*self.profile['tem_s'][i]*delta_n[:,i])/(1.5*self.profile['den_s'][i])
        imin = iwindow(t,w)
           
        color = ['k','m','b','c']

        f = np.empty([n_x,n_kinetic])

        for j in range(3):
           ax = fig.add_subplot(1,3,j+1)
           ax.grid(which="majorminor",ls=":")
           ax.grid(which="major",ls=":")
           ax.set_xlabel(r'$r/a$')
           ax.set_title(r'$'+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+'$')

           if j == 0:
              f = delta_n
              ax.set_ylabel(r'$\delta n$',color='k')
           if j == 1:
              f = delta_e
              ax.set_ylabel(r'$\delta E$',color='k')
           if j == 2:
              f = delta_t
              ax.set_ylabel(r'$\delta T$',color='k')

           for i in range(n_kinetic):
              stag = self.tagspec[i]
              ax.plot(self.profile['r'],f[:,i],label=stag,color=color[i])

        ax.legend()

