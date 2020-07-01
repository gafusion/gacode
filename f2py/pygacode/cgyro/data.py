import os
import numpy as np
import sys
import time
try:
   from gacodefuncs import *
except:
   from ..gacodefuncs import *
   
BYTE='float32'

# class to read cgyro output data
class cgyrodata:

   # constructor reads in basic (not all) simulation data
   def __init__(self,sim_directory,silent=False):

      self.silent = silent
      self.dir = sim_directory
      hastime = self.gettime()
      self.getgrid()
      if hastime:
         self.getdata()
      
   # standard routine to read binary or ASCII data 
   def extract(self,f):

      start = time.time()
      if os.path.isfile(self.dir+'bin'+f):
         fmt = 'bin'
         data = np.fromfile(self.dir+'bin'+f,dtype=BYTE)
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
      # Read time vector -- autodetect number of columns (ncol)
      #

      # START: 3-column output (delete this code eventually)
      tfile = self.dir+'out.cgyro.time'
      ncol = len(open(tfile,'r').readline().split())
   
      if ncol == 3:
         newlines = [] 
         for line in open(tfile,'r').readlines():
            newlines.append(line.strip()+'  0.0')
         outfile = open(tfile,'w')
         outfile.write("\n".join(newlines)) ; outfile.close()
         print('INFO: (data.py) Updated to 4-column format for out.cgyro.time')
        
      data = np.fromfile(self.dir+'out.cgyro.time',dtype='float',sep=' ')
      nt = len(data)//4
      data = np.reshape(data,(4,nt),'F')
      # END: 3-column output (delete this code eventually)

      self.t    = data[0,:]
      self.err1 = data[1,:]
      try:
         self.err2 = data[2,:]
      except:
         self.err2 = data[1,:]
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
         print('INFO: (data.py) Read data in '+fmt+f+'  '+t)

      f='.cgyro.aparb'
      t,fmt,data = self.extract(f)
      if fmt != 'null':
         self.aparb = np.reshape(data[0:nd],(2,self.n_theta*self.n_radial,nt),'F')
         print('INFO: (data.py) Read data in '+fmt+f+' '+t)

      f='.cgyro.bparb'
      t,fmt,data = self.extract(f)
      if fmt != 'null':
         self.bparb = np.reshape(data[0:nd],(2,self.n_theta*self.n_radial,nt),'F')
         print('INFO: (data.py) Read data in '+fmt+f+' '+t)
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Ballooning distribution
      #
      t,fmt,data = self.extract('.cgyro.hb')
      if fmt != 'null':
         self.hb = np.reshape(data,(2,self.n_radial*self.n_theta,
                                    self.n_species,self.n_xi,self.n_energy,nt),'F')
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.hb. '+t) 
         self.hb = self.hb/np.max(self.hb)
      #-----------------------------------------------------------------

   def getflux(self,cflux='auto'):

      if cflux == 'auto':
         if abs(self.gamma_e) > 0.0:
            usec = True
         else:
            usec = False
      elif cflux == 'on':
         usec = True
      else:
         usec = False

      #-----------------------------------------------------------------
      # Particle, momentum and energy fluxes
      #
      nt = self.n_time
      nd = self.n_species*3*self.n_field*self.n_n*nt

      if usec:
         t,fmt,data = self.extract('.cgyro.ky_cflux')
         if fmt != 'null':
            self.ky_flux = np.reshape(data[0:nd],(self.n_species,3,self.n_field,self.n_n,nt),'F')
            if not self.silent:
               print('INFO: (data.py) Read data in '+fmt+'.cgyro.ky_cflux. '+t)

      if not usec or fmt == 'null':
         t,fmt,data = self.extract('.cgyro.ky_flux')
         if fmt != 'null':  
            self.ky_flux = np.reshape(data[0:nd],(self.n_species,3,self.n_field,self.n_n,nt),'F')
            if not self.silent:
               print('INFO: (data.py) Read data in '+fmt+'.cgyro.ky_flux. '+t)
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

      """
      Do complicated spatial averages for xflux
      RESULT: self.lky_flux_ave
      """

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

      if moment == 'n':
         z = np.sum(self.lky_flux_n,axis=(3,4))
      elif moment == 'e':
         z = np.sum(self.lky_flux_e,axis=(3,4))
      elif moment == 'v':
         z = np.sum(self.lky_flux_v,axis=(3,4))
      else:
         raise ValueError('(xfluxave) Invalid moment.')

      #--------------------------------------------
      # Useful arrays required outside this routine
      self.lky_xr = np.zeros((ns,ng))
      self.lky_xi = np.zeros((ns,ng))
      self.lky_flux_ave = np.zeros((ns,2))
      #--------------------------------------------

      for ispec in range(ns):
         for l in range(ng):
            self.lky_xr[ispec,l] = average(z[0,l,ispec,:],self.t,w,0.0)*sc[ispec]
            self.lky_xi[ispec,l] = average(z[1,l,ispec,:],self.t,w,0.0)*sc[ispec]

         # Flux partial average over [-e,e]
         g0 = self.lky_xr[ispec,0]
         g1 = g0
         for l in range(1,ng):
            u = 2*np.pi*l*e
            g0 = g0+2*np.sin(u)*self.lky_xr[ispec,l]/u
            g1 = g1+2*np.sin(u)*self.lky_xr[ispec,l]/u*(-1)**l

         # Average over true (positive) interval
         self.lky_flux_ave[ispec,0] = g0
         # Average over negative interval
         self.lky_flux_ave[ispec,1] = g1

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
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_phi. '+t)
         self.kxky_phi = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_n,nt),'F')

      # 1b. kxky_apar
      t,fmt,data = self.extract('.cgyro.kxky_apar')
      if fmt != 'null':
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_apar. '+t)
         self.kxky_apar = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_n,nt),'F')

      # 1c. kxky_bpar
      t,fmt,data = self.extract('.cgyro.kxky_bpar')
      if fmt != 'null':
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_bpar. '+t)
         self.kxky_bpar = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_n,nt),'F')

      # 2. kxky_n
      nd = 2*self.n_radial*self.theta_plot*self.n_species*self.n_n*nt

      t,fmt,data = self.extract('.cgyro.kxky_n')
      if fmt != 'null':
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_n.   '+t)
         self.kxky_n = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_species,self.n_n,nt),'F')

      # 3. kxky_e
      t,fmt,data = self.extract('.cgyro.kxky_e')

      if fmt != 'null':
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_e.   '+t)
         self.kxky_e = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_species,self.n_n,nt),'F')

      # 4. kxky_e
      t,fmt,data = self.extract('.cgyro.kxky_v')

      if fmt != 'null':
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_v.   '+t)
         self.kxky_v = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_species,self.n_n,nt),'F')
      #-----------------------------------------------------------------

   def getgeo(self):

      """Read theta-dependent geometry functions"""

      t,fmt,data = self.extract('.cgyro.geo')
      if fmt != 'null':
         print('INFO: (data.py) Read data in '+fmt+'.cgyro.geo   '+t)
         self.geo = np.reshape(data,(self.n_theta,12),'F')

         self.geotag = []
         self.geotag.append('\theta')
         self.geotag.append('G_\\theta')
         self.geotag.append('|B|')
         self.geotag.append('\omega_\mathrm{stream}')
         self.geotag.append('\omega_\mathrm{trap}')
         self.geotag.append('\omega_\mathrm{rdrift}')
         self.geotag.append('\omega_\mathrm{adrift}')
         self.geotag.append('\omega_\mathrm{aprdrift}')
         self.geotag.append('\omega_\mathrm{cdrift}')
         self.geotag.append('\omega_\mathrm{crdrift}')
         self.geotag.append('\omega_\mathrm{gammap')
         self.geotag.append('k_\perp')

   def getgrid(self):

      #-----------------------------------------------------------------
      # Read grid data.
      #
      # NOTE: Grid data is packed, so unpack into sensible bits
      #
      data = np.fromfile(self.dir+'out.cgyro.grids',dtype='float32',sep=' ')

      self.n_n       = int(data[0])
      self.n_species = int(data[1])
      self.n_field   = int(data[2])
      self.n_radial  = int(data[3])
      self.n_theta   = int(data[4])
      self.n_energy  = int(data[5])
      self.n_xi      = int(data[6])
      self.m_box     = int(data[7])
      self.length    = float(data[8])
      self.n_global  = int(data[9])
      self.theta_plot= int(data[10])
      # Set l to last data index plus one.
      l=11

      self.p = np.array(data[l:l+self.n_radial],dtype=int)
      self.kx = 2*np.pi*self.p/self.length

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
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Equil file
      #
      data = np.fromfile(self.dir+'out.cgyro.equilibrium',dtype='float32',sep=' ')
      self.rmin          = data[0]
      self.rmaj          = data[1]
      self.q             = data[2]
      self.shear         = data[3]
      self.shift         = data[4]
      self.kappa         = data[5]
      self.s_kappa       = data[6]
      self.delta         = data[7]
      self.s_delta       = data[8]
      self.zeta          = data[9]
      self.s_zeta        = data[10]
      self.zmag          = data[11]
      self.dzmag         = data[12]
      self.shape_sin3    = data[13]
      self.shape_s_sin3  = data[14]
      self.shape_cos0    = data[15]
      self.shape_s_cos0  = data[16]
      self.shape_cos1    = data[17]
      self.shape_s_cos1  = data[18]
      self.shape_cos2    = data[19]
      self.shape_s_cos2  = data[20]
      self.shape_cos3    = data[21]
      self.shape_s_cos3  = data[22]
      self.rho           = data[23]
      self.ky0           = data[24]
      self.betae_unit    = data[25]
      self.beta_star     = data[26]
      self.lambda_star   = data[27]
      self.gamma_e       = data[28]
      self.gamma_p       = data[29]
      self.mach          = data[30]
      self.a_meters      = data[31]
      self.b_unit        = data[32]
      self.dens_norm     = data[33]
      self.temp_norm     = data[34]
      self.vth_norm      = data[35]
      self.mass_norm     = data[36]
      self.rho_star_norm = data[37]
      self.gamma_gb_norm = data[38]
      self.q_gb_norm     = data[39]
      self.pi_gb_norm    = data[40]
      # Define species vectors
      self.z      = np.zeros(self.n_species)
      self.mass   = np.zeros(self.n_species)
      self.dens   = np.zeros(self.n_species)
      self.temp   = np.zeros(self.n_species)
      self.dlnndr = np.zeros(self.n_species)
      self.dlntdr = np.zeros(self.n_species)
      self.nu     = np.zeros(self.n_species)
      for i in range(self.n_species):
         self.z[i]      = data[41+7*i]
         self.mass[i]   = data[42+7*i]
         self.dens[i]   = data[43+7*i]
         self.temp[i]   = data[44+7*i]
         self.dlnndr[i] = data[45+7*i]
         self.dlntdr[i] = data[46+7*i]
         self.nu[i]     = data[47+7*i]
      if not self.silent:
         print('INFO: (data.py) Read data in out.cgyro.equilibrium.')

      #-----------------------------------------------------------------

   def getnorm(self,norm):

      if norm == 'elec':

         # Standard normalizations

         self.tnorm  = self.t
         self.tstr   = TIME

         self.fnorm  = self.freq 
         self.fstr   = [r'$(a/c_{sD})\, \omega$',r'$(a/c_{sD})\, \gamma$']

         self.rhonorm = self.rho
         self.rhostr = r'$\rho_{sD}$'

         self.kynorm = self.ky
         self.kstr   = r'$k_y \rho_{sD}$'

         self.qc     = 1.0
         self.gbnorm = '_\mathrm{GBD}'
         
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
         self.tstr  = r'$(v_'+str(i)+'/a) \, t$'

         self.fnorm = self.freq/vc
         self.fstr  = [r'$(a/v_'+str(i)+') \, \omega$',r'$(a/v_'+str(i)+') \, \gamma$']

         # Convert rho_sD=csD/Omega_i to rho_i=vi/Omega_i
         rhoc = vc/(self.z[i]/self.mass[i])

         self.rhonorm = self.rho*rhoc
         self.rhostr  = r'$\rho_'+str(i)+'$'

         self.kynorm = self.ky*rhoc
         self.kstr   = r'$k_y \rho_'+str(i)+'$'

          # Convert Q_GBD to Q_GBi
         self.qc = vc*(self.temp[i]/te)*rhoc**2
         self.gbnorm = '_\mathrm{GB'+str(i)+'}'

