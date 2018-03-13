import os
import numpy as np
import sys

BYTE='float32'

class cgyrodata:

   """CGYRO output data class."""

   def __init__(self, sim_directory):

      """Constructor reads in basic (not all) simulation data."""

      self.dir = sim_directory
      self.getdata()
        
   def extract(self,f):

      import sys
      import os
      import numpy as np
      import time

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
         
      if int(sys.version_info[2]) > 6:
         t = 'TIME = '+"{:.3e}".format(time.time()-start)+' s.'
      else:        
         t = 'TIME = '+str(time.time()-start)
 

      return t,fmt,data
       
          
   def getdata(self):

      """Initialize smaller data objects (don't load larger ones)"""

      import numpy as np
      
      #-----------------------------------------------------------------
      # Read time vector.
      #
      data = np.fromfile(self.dir+'out.cgyro.time',dtype='float',sep=' ')
      nt = len(data)/3
      data = np.reshape(data,(3,nt),'F')
      self.t    = data[0,:]
      self.err1 = data[1,:]
      try:
         self.err2 = data[2,:]
      except:
         self.err2 = data[1,:]
      self.n_time = nt   
      print "INFO: (data.py) Read time vector in out.cgyro.time."
      #-----------------------------------------------------------------

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
      self.thetab = np.array(data[mark:mark+self.n_theta*(self.n_radial/self.m_box)])        
         
      mark = mark+self.n_theta*(self.n_radial/self.m_box)
      self.ky = np.array(data[mark:mark+self.n_n])

      mark = mark+self.n_n
      self.alphadiss = np.array(data[mark:mark+self.n_n])

      mark = mark+self.n_n
      self.radialdiss = np.array(data[mark:mark+self.n_radial])

      print "INFO: (data.py) Read grid data in out.cgyro.grids."
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Linear frequency
      #
      nd = 2*self.n_n*nt
      t,fmt,data = self.extract('.cgyro.freq')
      if fmt != 'null':  
         self.freq = np.reshape(data,(2,self.n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".cgyro.freq. "+t 
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Equil file
      #
      try:
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
         self.rho           = data[13]
         self.ky0           = data[14]
         self.betae_unit    = data[15]
         self.beta_star     = data[16]
         self.lambda_star   = data[17]
         self.gamma_e       = data[18]
         self.gamma_p       = data[19]
         self.mach          = data[20]
         self.a_meters      = data[21]
         self.b_unit        = data[22]
         self.dens_norm     = data[23]
         self.temp_norm     = data[24]
         self.vth_norm      = data[25]
         self.mass_norm     = data[26]
         self.rho_star_norm = data[27]
         self.gamma_gb_norm = data[28]
         self.q_gb_norm     = data[29]
         self.pi_gb_norm    = data[30]
         # Define species vectors
         self.z      = np.zeros(self.n_species)
         self.mass   = np.zeros(self.n_species)
         self.dens   = np.zeros(self.n_species)
         self.temp   = np.zeros(self.n_species)
         self.dlnndr = np.zeros(self.n_species)
         self.dlntdr = np.zeros(self.n_species)
         self.nu     = np.zeros(self.n_species)
         for i in range(self.n_species):
            self.z[i]      = data[31+7*i]
            self.mass[i]   = data[32+7*i]
            self.dens[i]   = data[33+7*i]
            self.temp[i]   = data[34+7*i]
            self.dlnndr[i] = data[35+7*i]
            self.dlntdr[i] = data[36+7*i]
            self.nu[i]     = data[37+7*i]
         print "INFO: (data.py) Read data in out.cgyro.equilibrium."
      except:
         print "WARNING: (data.py) Could not read out.cgyro.equilibrium."
         pass

      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Ballooning potentials
      #
      f='.cgyro.phib'
      t,fmt,data = self.extract(f)
      if fmt != 'null':
         self.phib = np.reshape(data,(2,self.n_theta*self.n_radial,nt),'F')
         print 'INFO: (data.py) Read data in '+fmt+f+'  '+t 

      f='.cgyro.aparb'
      t,fmt,data = self.extract(f)
      if fmt != 'null':
         self.aparb = np.reshape(data,(2,self.n_theta*self.n_radial,nt),'F')
         print 'INFO: (data.py) Read data in '+fmt+f+' '+t 

      f='.cgyro.bparb'
      t,fmt,data = self.extract(f)
      if fmt != 'null':
         self.bparb = np.reshape(data,(2,self.n_theta*self.n_radial,nt),'F')
         print 'INFO: (data.py) Read data in '+fmt+f+' '+t 
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Ballooning distribution
      #
      t,fmt,data = self.extract('.cgyro.hb')
      if fmt != 'null':
         self.hb = np.reshape(data,(2,self.n_radial*self.n_theta/self.m_box,
                                       self.n_species,self.n_xi,self.n_energy,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".cgyro.hb. "+t 
         self.hb = self.hb/np.max(self.hb)
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Compressed particle and energy fluxes
      #
      nd = self.n_species*nt
      try:
         data = np.loadtxt(self.dir+'out.cgyro.flux_n',dtype='float')
         self.flux_n = np.transpose(data[:,1:])
         print "INFO: (data.py) Read data in out.cgyro.flux_n."
      except:
         pass

      try:
         data = np.loadtxt(self.dir+'out.cgyro.flux_e',dtype='float')
         self.flux_e = np.transpose(data[:,1:])
         print "INFO: (data.py) Read data in out.cgyro.flux_e."
      except:
         pass 
      #-----------------------------------------------------------------

   def getflux(self):

      import numpy as np

      # Convenience definition
      nt = self.n_time

      #-----------------------------------------------------------------
      # Particle and energy fluxes
      #
      nd = self.n_species*3*self.n_field*self.n_n*nt
      t,fmt,data = self.extract('.cgyro.ky_flux')
      if fmt != 'null':  
         self.ky_flux = np.reshape(data[0:nd],(self.n_species,3,self.n_field,self.n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".cgyro.ky_flux. "+t 
      #-----------------------------------------------------------------

   def getbigflux(self):

      """Larger flux files (hopefully binary)"""

      import numpy as np

      # Convenience definition
      nt = self.n_time

      #-----------------------------------------------------------------
      # Particle and energy fluxes
      #
      nd = self.n_radial*self.n_species*self.n_n*nt
      t,fmt,data = self.extract('.cgyro.kxky_flux_e')
      if fmt != 'null':  
         self.kxky_flux_e = np.reshape(data[0:nd],(self.n_radial,self.n_species,self.n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".cgyro.kxky_flux_e. "+t
      #
      nd = self.n_radial*self.n_species*self.n_n*nt
      t,fmt,data = self.extract('.cgyro.kxky_flux_n')
      if fmt != 'null':  
         self.kxky_flux_n = np.reshape(data[0:nd],(self.n_radial,self.n_species,self.n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".cgyro.kxky_flux_n. "+t
      #-----------------------------------------------------------------

   def getxflux(self):

      """Global-spectral flux files"""

      import numpy as np

      # Convenience definition
      nt = self.n_time

      #-----------------------------------------------------------------
      # Particle and energy fluxes
      #
      ng = self.n_global+1
      nd = 2*ng*self.n_species*self.n_n*nt

      t,fmt,data = self.extract('.cgyro.lky_flux_n')
      if fmt != 'null':  
         self.lky_flux_n = np.reshape(data[0:nd],(2,ng,self.n_species,self.n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".cgyro.lky_flux_n. "+t

      t,fmt,data = self.extract('.cgyro.lky_flux_e')
      if fmt != 'null':  
         self.lky_flux_e = np.reshape(data[0:nd],(2,ng,self.n_species,self.n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".cgyro.lky_flux_e. "+t

      t,fmt,data = self.extract('.cgyro.lky_flux_v')
      if fmt != 'null':  
         self.lky_flux_v = np.reshape(data[0:nd],(2,ng,self.n_species,self.n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".cgyro.lky_flux_v. "+t      
      #-----------------------------------------------------------------

   def getbigfield(self):

      """Larger field files"""

      import numpy as np

      # Convenience definition
      nt = self.n_time

      #-----------------------------------------------------------------
      # Read complex fields
      #
      # 1a. kxky_phi
      nd = 2*self.n_radial*self.theta_plot*self.n_n*nt
      t,fmt,data = self.extract('.cgyro.kxky_phi')
      if fmt != 'null':
         print 'INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_phi. '+t
         self.kxky_phi = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_n,nt),'F')
         self.kxky_phi_abs = np.sqrt(self.kxky_phi[0,:,:,:,:]**2+self.kxky_phi[1,:,:,:,:]**2)

      # 1b. kxky_apar
      nd = 2*self.n_radial*self.theta_plot*self.n_n*nt
      t,fmt,data = self.extract('.cgyro.kxky_apar')
      if fmt != 'null':
         print 'INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_apar. '+t
         self.kxky_apar = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_n,nt),'F')
         self.kxky_apar_abs = np.sqrt(self.kxky_apar[0,:,:,:,:]**2+self.kxky_apar[1,:,:,:,:]**2)

      # 1c. kxky_bpar
      nd = 2*self.n_radial*self.theta_plot*self.n_n*nt
      t,fmt,data = self.extract('.cgyro.kxky_bpar')
      if fmt != 'null':
         print 'INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_bpar. '+t
         self.kxky_bpar = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_n,nt),'F')
         self.kxky_bpar_abs = np.sqrt(self.kxky_bpar[0,:,:,:,:]**2+self.kxky_bpar[1,:,:,:,:]**2)

      # 2. kxky_n
      nd = 2*self.n_radial*self.theta_plot*self.n_species*self.n_n*nt
      t,fmt,data = self.extract('.cgyro.kxky_n')

      if fmt != 'null':
         print 'INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_n.   '+t
         self.kxky_n = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_species,self.n_n,nt),'F')

      # 3. kxky_e
      nd = 2*self.n_radial*self.theta_plot*self.n_species*self.n_n*nt
      t,fmt,data = self.extract('.cgyro.kxky_e')

      if fmt != 'null':
         print 'INFO: (data.py) Read data in '+fmt+'.cgyro.kxky_e.   '+t
         self.kxky_e = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_species,self.n_n,nt),'F')
      #-----------------------------------------------------------------

   def getgeo(self):

      """Read theta-dependent geometry functions"""

      import numpy as np

      t,fmt,data = self.extract('.cgyro.geo')
      if fmt != 'null':
         print 'INFO: (data.py) Read data in '+fmt+'.cgyro.geo   '+t
         self.geo = np.reshape(data,(self.n_theta,12),'F')

         self.geotag = []
         self.geotag.append('\theta')
         self.geotag.append('w_\\theta')
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
 
