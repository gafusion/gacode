import os
import numpy as np
import sys
import time
try:
   from gacodefuncs import *
except:
   from ..gacodefuncs import *

# class to read cgyro output data
class cgyrodata:

   # constructor reads in basic (not all) simulation data
   def __init__(self,sim_directory,silent=False,fast=False):

      self.silent = silent
      self.dir = sim_directory

      # read time vector
      hastime = self.gettime()

      # get mesh/input data
      self.getgrid()

      # set data resolution
      if self.hiprec_flag:
         # double
         self.BYTE = 'float64'
      else:
         # single
         self.BYTE = 'float32'

      if not self.silent:
         print('INFO: (data.py) Detected precision '+self.BYTE)

         
      if hastime and not fast:
         self.getdata()
      
   # standard routine to read binary or ASCII data 
   def extract(self,f):

      start = time.time()
      if os.path.isfile(self.dir+'bin'+f):
         fmt = 'bin'
         data = np.fromfile(self.dir+'bin'+f,dtype=self.BYTE)
      elif os.path.isfile(self.dir+'out'+f):
         fmt = 'out'
         data = np.fromfile(self.dir+'out'+f,dtype='float',sep=' ')
      else:
         fmt  = 'null'
         data = []

      t = 'TIME = '+'{:.3e}'.format(time.time()-start)+' s.'

      return t,fmt,data

   def gettime(self):

      if not os.path.isfile(self.dir+'out.cgyro.time'):
         print('INFO: (data.py) No time record exists.')
         return False
      
      #-----------------------------------------------------------------
      # Read time vector
      #        
      data = np.fromfile(self.dir+'out.cgyro.time',dtype='float',sep=' ')
      nt = len(data)//4
      data = np.reshape(data,(4,nt),'F')

      self.t    = data[0,:]
      self.err1 = data[1,:]
      self.err2 = data[2,:]
      self.n_time = nt
      if not self.silent:
         print('INFO: (data.py) Read time vector in out.cgyro.time.')

      #-----------------------------------------------------------------

      return True
   
   def getdata(self):

      nt = self.n_time
      
      #-----------------------------------------------------------------
      # Linear frequency
      #
      nd = 2*self.n_n*nt
      t,fmt,data = self.extract('.cgyro.freq')
      if fmt != 'null':  
         self.freq = np.reshape(data[0:nd],(2,self.n_n,nt),'F')
         if not self.silent:
            print('INFO: (data.py) Read data in '+fmt+'.cgyro.freq. '+t) 
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Ballooning potentials
      #
      nd = 2*self.n_theta*self.n_radial*nt
      f='.cgyro.phib'
      t,fmt,data = self.extract(f)
      if fmt != 'null':
         self.phib = np.reshape(data[0:nd],(2,self.n_theta*self.n_radial,nt),'F')
         if not self.silent:
            print('INFO: (data.py) Read data in '+fmt+f+'  '+t)

      f='.cgyro.aparb'
      t,fmt,data = self.extract(f)
      if fmt != 'null':
         self.aparb = np.reshape(data[0:nd],(2,self.n_theta*self.n_radial,nt),'F')
         if not self.silent:
            print('INFO: (data.py) Read data in '+fmt+f+' '+t)

      f='.cgyro.bparb'
      t,fmt,data = self.extract(f)
      if fmt != 'null':
         self.bparb = np.reshape(data[0:nd],(2,self.n_theta*self.n_radial,nt),'F')
         if not self.silent:
            print('INFO: (data.py) Read data in '+fmt+f+' '+t)
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Ballooning distribution
      #
      t,fmt,data = self.extract('.cgyro.hb')
      if fmt != 'null':
         self.hb = np.reshape(data,(2,self.n_radial*self.n_theta,
                                    self.n_species,self.n_xi,self.n_energy,nt),'F')
         if not self.silent:
            print('INFO: (data.py) Read data in '+fmt+'.cgyro.hb. '+t) 
         self.hb = self.hb/np.max(self.hb)
      #-----------------------------------------------------------------

   def getflux(self,cflux='auto'):

      if cflux == 'auto':
         if abs(self.gamma_e) > 0.0 and self.n_n > 1:
            usec = True
         else:
            usec = False
      elif cflux == 'on':
         usec = True
      else:
         usec = False

      #-----------------------------------------------------------------
      # Particle, momentum, energy fluxes, and exchange (autodetected)
      #      
      nt = self.n_time
      nd = self.n_species*self.n_field*self.n_n*nt

      if usec:
         t,fmt,data = self.extract('.cgyro.ky_cflux')
         m = data.shape[0]//nd
         if fmt != 'null':
            self.ky_flux = np.reshape(data[0:nd*m],(self.n_species,m,self.n_field,self.n_n,nt),'F')
            if not self.silent:
               print('INFO: (data.py) Read data in '+fmt+'.cgyro.ky_cflux. '+t)

      if not usec or fmt == 'null':
         t,fmt,data = self.extract('.cgyro.ky_flux')
         m = data.shape[0]//nd
         if fmt != 'null':  
            self.ky_flux = np.reshape(data[0:nd*m],(self.n_species,m,self.n_field,self.n_n,nt),'F')
            if not self.silent:
               print('INFO: (data.py) Read data in '+fmt+'.cgyro.ky_flux. '+t)

      self.n_flux = m
      #-----------------------------------------------------------------

      return usec

   def getxflux(self):

      """Global-spectral flux files (optional)"""

      #-----------------------------------------------------------------
      # Particle and energy fluxes
      #
      nt = self.n_time
      ng = self.n_global+1
      nd = 2*ng*self.n_species*self.n_field*self.n_n*nt

      t,fmt,data = self.extract('.cgyro.lky_flux_n')
      if fmt != 'null':
         self.lky_flux_n = np.reshape(data[0:nd],(2,ng,self.n_species,self.n_field,self.n_n,nt),'F')
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.lky_flux_n. '+t)

      t,fmt,data = self.extract('.cgyro.lky_flux_e')
      if fmt != 'null':
         self.lky_flux_e = np.reshape(data[0:nd],(2,ng,self.n_species,self.n_field,self.n_n,nt),'F')
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.lky_flux_e. '+t)

      t,fmt,data = self.extract('.cgyro.lky_flux_v')
      if fmt != 'null':
         self.lky_flux_v = np.reshape(data[0:nd],(2,ng,self.n_species,self.n_field,self.n_n,nt),'F')
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.lky_flux_v. '+t)
      #-----------------------------------------------------------------

   def xfluxave(self,w,moment,e=0.2,nscale=0):

      # Averaging functionfor global fluxes

      print('INFO: (xfluxave) Computing partial-domain averages')

      ns = self.n_species
      ng = self.n_global+1

      sc = np.zeros(ns)
      if nscale == 1:
         # Find ne
         for ispec in range(ns):
            if self.z[ispec] < 0.0:
               ne = self.dens[ispec]
         sc[:] = ne/self.dens[:]
      else:
         sc[:] = 1.0

      # Sum moments over fields (3) and toroidal modes (4)
      # NOTE: lky_flux_n[2,ng,ns,nfield,n_n,nt]
      if moment == 'n':
         z = np.sum(self.lky_flux_n,axis=(3,4))
      elif moment == 'e':
         z = np.sum(self.lky_flux_e,axis=(3,4))
      elif moment == 'v':
         z = np.sum(self.lky_flux_v,axis=(3,4))
         
      #--------------------------------------------
      # Arrays required outside this routine
      self.lky_xr = np.zeros([ng,ns])
      self.lky_xi = np.zeros([ng,ns])
      self.lky_flux_ave = np.zeros([ns,2])
      #--------------------------------------------

      imin,imax = time_index(self.t,w)

      # Time averages of real and imaginary parts
      self.lky_xr[:,:] = time_average(z[0,:,:,:],self.t,imin,imax)
      self.lky_xi[:,:] = time_average(z[1,:,:,:],self.t,imin,imax)

      l = np.arange(1,ng)
      u = (2*np.pi*e)*l
      for ispec in range(ns):
         # Flux partial average over [-e,e]
         g0 = self.lky_xr[0,ispec]+2*np.sum(np.sin(u)*self.lky_xr[l,ispec]/u)
         g1 = self.lky_xr[0,ispec]+2*np.sum(np.sin(u)*self.lky_xr[l,ispec]/u*(-1)**l)

         # Average over true (positive) interval
         self.lky_flux_ave[ispec,0] = g0*sc[ispec]
         # Average over negative interval
         self.lky_flux_ave[ispec,1] = g1*sc[ispec]

   def getbigfield(self):

      """Larger field files"""

      #-----------------------------------------------------------------
      # Read complex fields
      #

      nt = self.n_time
      nd = 2*self.n_radial*self.theta_plot*self.n_n*nt

      # 1a. kxky_phi
      t,fmt,data = self.extract('.cgyro.kxky_phi')
      if fmt != 'null':
         self.kxky_phi = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_n,nt),'F')
         if not self.silent:
            print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_phi. '+t)

      # 1b. kxky_apar
      t,fmt,data = self.extract('.cgyro.kxky_apar')
      if fmt != 'null':
        self.kxky_apar = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_n,nt),'F')
        if not self.silent:
           print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_apar. '+t)
 
      # 1c. kxky_bpar
      t,fmt,data = self.extract('.cgyro.kxky_bpar')
      if fmt != 'null':
         self.kxky_bpar = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_n,nt),'F')
         if not self.silent:
            print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_bpar. '+t)

      # 2. kxky_n
      nd = 2*self.n_radial*self.theta_plot*self.n_species*self.n_n*nt

      t,fmt,data = self.extract('.cgyro.kxky_n')
      if fmt != 'null':
         self.kxky_n = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_species,self.n_n,nt),'F')
         if not self.silent:
            print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_n.   '+t)

      # 3. kxky_e
      t,fmt,data = self.extract('.cgyro.kxky_e')

      if fmt != 'null':
         self.kxky_e = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_species,self.n_n,nt),'F')
         if not self.silent:
            print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_e.   '+t)

      # 4. kxky_v
      t,fmt,data = self.extract('.cgyro.kxky_v')

      if fmt != 'null':
        self.kxky_v = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_species,self.n_n,nt),'F')
        if not self.silent:
           print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_v.   '+t)

      # 5. triad
      t,fmt,data = self.extract('.cgyro.triad')
      if fmt != 'null':
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.triad.   '+t)
         nd = 2*self.n_n*self.n_radial*self.n_species*nt*8 
         self.triad = np.reshape(data[0:nd],(2,self.n_species,self.n_radial,8,self.n_n,nt),'F')

       #-----------------------------------------------------------------

   def getgeo(self):

      """Read theta-dependent geometry functions"""

      t,fmt,data = self.extract('.cgyro.geo')
      if fmt != 'null':
         nfunc = len(data)//self.n_theta
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.geo   '+t)
         self.geo = np.reshape(data,(self.n_theta,nfunc),'F')

         tags = [r'\theta',
                 r'G_\theta',
                 r'|B|',
                 r'\omega_\mathrm{stream}',
                 r'\omega_\mathrm{trap}',
                 r'\omega_\mathrm{rdrift}',
                 r'\omega_\mathrm{adrift}',
                 r'\omega_\mathrm{aprdrift}',
                 r'\omega_\mathrm{cdrift}',
                 r'\omega_\mathrm{crdrift}',
                 r'\omega_\mathrm{gammap',
                 r'k_\perp',
                 r'\Theta']

         self.geotag = []
         for i in range(nfunc):
            self.geotag.append(tags[i])
            

   # Function to pull elements out of array, checking for end of array
   def eget(self,data,p):
      
      if p >= len(data) or p == -1:
         p = -1
         d0 = 0
      else:
         d0 = data[p]
         p = p+1

      return d0,p
      
   def getgrid(self):

      #-----------------------------------------------------------------
      # Read grid data.
      #
      # NOTE: Grid data is packed, so unpack into sensible bits
      #
      data = np.fromfile(self.dir+'out.cgyro.grids',sep=' ')

      self.n_n       = int(data[0])
      self.n_species = int(data[1])
      self.n_field   = int(data[2])
      self.n_radial  = int(data[3])
      self.n_theta   = int(data[4])
      self.n_energy  = int(data[5])
      self.n_xi      = int(data[6])
      self.m_box     = int(data[7])
      self.length    = float(data[8]) # L/rho
      self.n_global  = int(data[9])
      self.theta_plot= int(data[10])
      # Set l to last data index plus one.
      l=11

      self.p = np.array(data[l:l+self.n_radial],dtype=int)

      if self.p[0] > 0:
         # zonal flow test
         self.kx = 2*np.pi*self.p[:]/self.length
         self.zf_test = True
      else:
         # kx*rho (ignore leftmost "special" element)
         self.kx = 2*np.pi*self.p[1:]/self.length
         self.zf_test = False

      mark = l+self.n_radial
      self.theta = np.array(data[mark:mark+self.n_theta])

      mark = mark+self.n_theta
      self.energy = np.array(data[mark:mark+self.n_energy])

      mark = mark+self.n_energy
      self.xi   = np.array(data[mark:mark+self.n_xi])

      mark = mark+self.n_xi
      self.thetab = np.array(data[mark:mark+self.n_theta*self.n_radial//self.m_box])

      mark = mark+self.n_theta*(self.n_radial//self.m_box)
      self.ky = np.array(data[mark:mark+self.n_n])

      mark = mark+self.n_n
      self.alphadiss = np.array(data[mark:mark+self.n_n])

      mark = mark+self.n_n
      self.radialdiss = np.array(data[mark:mark+self.n_radial])

      if not self.silent:
         print('INFO: (data.py) Read grid data in out.cgyro.grids.')
      #-----------------------------------------------------------------

      #--------------------------------------------------------
      # Construct theta_plot values (self.thetap[:])
      self.thetap = np.zeros(self.theta_plot)
      if self.theta_plot == 1:
         self.thetap[0] = 0.0
      else:
         m = self.n_theta//self.theta_plot
         for i in range(self.theta_plot):
            self.thetap[i] = self.theta[m*i]

      # Construct k_perp
      # NOTE: kx,ky and kperp are kx*rho, ky*rho, kperp*rho
      if not self.zf_test:
         nx = self.n_radial-1
         ny = self.n_n
         self.kperp = np.sqrt(np.outer(self.kx[:]**2,np.ones(ny))+
                              np.outer(np.ones(nx),self.ky[:]**2))
      
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Equil file
      #
      nshape = 7
      ns = self.n_species
      p = 0
      data = np.fromfile(self.dir+'out.cgyro.equilibrium',sep=' ')
      self.rmin,p    = self.eget(data,p)
      self.rmaj,p    = self.eget(data,p)        
      self.q,p       = self.eget(data,p)           
      self.shear,p   = self.eget(data,p)         
      self.shift,p   = self.eget(data,p)       
      self.kappa,p   = self.eget(data,p)     
      self.s_kappa,p = self.eget(data,p)       
      self.delta,p   = self.eget(data,p)       
      self.s_delta,p = self.eget(data,p)     
      self.zeta,p    = self.eget(data,p)     
      self.s_zeta,p  = self.eget(data,p)    
      self.zmag,p    = self.eget(data,p)      
      self.dzmag,p   = self.eget(data,p)    
      
      self.shape_sin   = np.zeros(nshape)
      self.shape_s_sin = np.zeros(nshape)
      self.shape_cos   = np.zeros(nshape)
      self.shape_s_cos = np.zeros(nshape)
      for i in range(3,nshape):
         self.shape_sin[i],p   = self.eget(data,p)   
         self.shape_s_sin[i],p = self.eget(data,p)  
      for i in range(nshape):
         self.shape_cos[i],p   = self.eget(data,p)
         self.shape_s_cos[i],p = self.eget(data,p)
         
      self.rho,p           = self.eget(data,p) 
      self.ky0,p           = self.eget(data,p) 
      self.betae_unit,p    = self.eget(data,p) 
      self.beta_star,p     = self.eget(data,p) 
      self.lambda_star,p   = self.eget(data,p) 
      self.gamma_e,p       = self.eget(data,p) 
      self.gamma_p,p       = self.eget(data,p) 
      self.mach,p          = self.eget(data,p) 
      self.a_meters,p      = self.eget(data,p) 
      self.b_unit,p        = self.eget(data,p) 
      self.b_gs2,p         = self.eget(data,p) 
      self.dens_norm,p     = self.eget(data,p) 
      self.temp_norm,p     = self.eget(data,p) 
      self.vth_norm,p      = self.eget(data,p) 
      self.mass_norm,p     = self.eget(data,p) 
      self.rho_star_norm,p = self.eget(data,p) 
      self.gamma_gb_norm,p = self.eget(data,p) 
      self.q_gb_norm,p     = self.eget(data,p) 
      self.pi_gb_norm,p    = self.eget(data,p) 

      # Define species vectors
      self.z      = np.zeros(ns)
      self.mass   = np.zeros(ns)
      self.dens   = np.zeros(ns)
      self.temp   = np.zeros(ns)
      self.dlnndr = np.zeros(ns)
      self.dlntdr = np.zeros(ns)
      self.nu     = np.zeros(ns)
      for i in range(ns):
         self.z[i],p      = self.eget(data,p)
         self.mass[i],p   = self.eget(data,p) 
         self.dens[i],p   = self.eget(data,p) 
         self.temp[i],p   = self.eget(data,p) 
         self.dlnndr[i],p = self.eget(data,p) 
         self.dlntdr[i],p = self.eget(data,p) 
         self.nu[i],p     = self.eget(data,p)

      # Added 3 May 2024
      self.sdlnndr = np.zeros(ns)
      self.sdlntdr = np.zeros(ns)
      self.sbeta = np.zeros(ns)
      for i in range(ns):
         self.sdlnndr[i],p = self.eget(data,p) 
         self.sdlntdr[i],p = self.eget(data,p)

      # Updated 19 Feb 2025
      self.sbeta,p = self.eget(data,p)
      self.z_eff,p = self.eget(data,p) 
      self.hiprec_flag,p = self.eget(data,p) 

      if p == -1:
         print('WARNING: (getgrid) Data format outdated. Please run cgyro -t')
      else:
         if not self.silent:
            print('INFO: (getgrid) Read {:d} entries out.cgyro.equilibrium.'.format(p))

      #-----------------------------------------------------------------

   def getnorm(self,norm):

      if norm == 'elec':

         # Standard normalizations

         self.tnorm  = self.t
         self.tstr   = TIME

         try:
            self.fnorm  = self.freq 
            self.fstr   = [r'$(a/c_s)\, \omega$',r'$(a/c_s)\, \gamma$']
         except:
            pass

         self.rhonorm = self.rho
         self.rhoi    = r'\rho_s'

         self.kxnorm = self.kx
         self.kynorm = self.ky
         self.kxstr  = r'$k_x \rho_s$'
         self.kystr  = r'$k_\theta \rho_s$'

         self.qc     = 1.0
         self.gbnorm = r'_\mathrm{GBD}'

         print('INFO: (getnorm) Using deuterium norm (rho_s = rho_sD, c_s = c_sD, etc)')

      else:

         # Species-i normalizations
         
         i = int(norm)

         te = 1.0
         for j in range(self.n_species):
            if self.z[j] < 0.0:
               te = self.temp[j]

         # Convert csD=sqrt(Te/mD) to vi=sqrt(Ti/mi)
         vc = np.sqrt(self.temp[i]/te/self.mass[i])
 
         self.tnorm = self.t*vc
         self.tstr  = r'$(v_'+str(i)+r'/a) \, t$'

         self.fnorm = self.freq/vc
         self.fstr  = [r'$(a/v_'+str(i)+r') \, \omega$',r'$(a/v_'+str(i)+r') \, \gamma$']

         # Convert rho_sD=csD/Omega_i to rho_i=vi/Omega_i
         rhoc = vc/(self.z[i]/self.mass[i])

         self.rhonorm = self.rho*rhoc
         self.rhoi = r'\rho_'+str(i)

         self.kxnorm = self.kx*rhoc
         self.kynorm = self.ky*rhoc
         self.kxstr  = r'$k_x \rho_'+str(i)+'$'
         self.kystr  = r'$k_\theta \rho_'+str(i)+'$'

          # Convert Q_GBD to Q_GBi
         self.qc = vc*(self.temp[i]/te)*rhoc**2
         self.gbnorm = r'_\mathrm{GB'+str(i)+'}'

         print('INFO: (getnorm) Using species '+norm+' norm')

      return self.tnorm

      
   def kxky_select(self,theta,field,moment,species,gbnorm=False):

      itheta,thetapi = indx_theta(theta,self.theta_plot)

      # Set norm if unset
      if 'self.rhonorm' not in locals():
         self.getnorm('elec')

      if moment == 'phi':
         if field == 0:
            f  = self.kxky_phi[0,1:,itheta,:,:]+1j*self.kxky_phi[1,1:,itheta,:,:]
            ft = TEXPHI
         elif field == 1:
            f  = self.kxky_apar[0,1:,itheta,:,:]+1j*self.kxky_apar[1,1:,itheta,:,:]
            ft = TEXAPAR
         else:
            f  = self.kxky_bpar[0,1:,itheta,:,:]+1j*self.kxky_bpar[1,1:,itheta,:,:]
            ft = TEXBPAR
      elif moment == 'k':
            f  = self.kxky_phi[0,1:,itheta,:,:]+1j*self.kxky_phi[1,1:,itheta,:,:]
            for i in range(self.n_time):
               # self.kperp -> kperp*rhos 
               f[:,:,i] = self.kperp[:,:]*f[:,:,i]*self.rhonorm/self.rho
            ft = r'k_\perp'+self.rhoi+r'\, \phi'
      elif moment == 'n':
         f  = self.kxky_n[0,1:,itheta,species,:,:]+1j*self.kxky_n[1,1:,itheta,species,:,:]
         ft = TEXDN
      elif moment == 'e':
         f  = self.kxky_e[0,1:,itheta,species,:,:]+1j*self.kxky_e[1,1:,itheta,species,:,:]
         ft = TEXDE
      elif moment == 'v':
         f  = self.kxky_v[0,1:,itheta,species,:,:]+1j*self.kxky_v[1,1:,itheta,species,:,:]
         ft = TEXDV

      # 3D structure: f[r,n,time]

      # 2024.04 moments are divided by rho in cgyro, so denormalize now
      if moment == 'n' or moment == 'e' or moment == 'v':
         f[:,:,:] = f[:,:,:]*self.rhonorm
         
      if gbnorm:
         # gyrBohm normalization
         f[:,:,:] = f[:,:,:]/self.rhonorm
      
      return f,ft

   def ylabeler(self,sub,ft,abs=True,sq=False,tave=False,sqrt=False,d=False,gb=True):

      if sub == '+' or sub == 'null':
         z = 'n'
      else:
         z = sub

      # add subscript
      u = ft+'_'+z

      if d:
         u = r'(\mathrm{Gradient~scale~length})~~\partial_r \;'+u
         
      # gyroBohm norm
      if gb:
         u = u+'/'+self.rhoi

      # Absolute value
      if abs:
         u = r'\left|'+u+r'\right|'

      if sub == '+':
         u = r'\sum_{n>0}'+u
      if sub == 'null':
         u = r'\sum_{n}'+u

      # Squared?
      if sq:
         u = u+'^2'

      if tave:
         u = r'\left\langle'+u+r'\right\rangle'

      if sqrt:
         u = r'\sqrt{'+u+'}'
         
      return '$'+u+'$'
