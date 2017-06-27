import os
import numpy as np
import sys

class cgyrodata:

    """CGYRO output data class."""

    def __init__(self, sim_directory):

        """Constructor reads in basic (not all) simulation data."""

        self.dir = sim_directory
        self.getdata()
   

    def getdata(self):

        """Initialize smaller data objects (don't load larger ones)"""

        import numpy as np
        import time

        #-----------------------------------------------------------------
        # Read time vector.
        #
        data = np.loadtxt(self.dir+'out.cgyro.time')
        self.t    = data[:,0]
        self.err1 = data[:,1]
        try:
            self.err2 = data[:,2]
        except:
            self.err2 = data[:,1]
        self.n_time = len(self.t)   
        print "INFO: (data.py) Read time vector in out.cgyro.time."
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Read grid data.
        # 
        # NOTE: Grid data is packed, so unpack into sensible bits
        #
        data = np.loadtxt(self.dir+'out.cgyro.grids')

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

        # Convenience definition
        nt = self.n_time

        #-----------------------------------------------------------------
        # Linear frequency
        #
        try:
            data = np.fromfile(self.dir+'out.cgyro.freq',dtype='float',sep=" ")
            self.freq = np.reshape(data,(2,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.freq."
        except:
            pass

        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Equil file
        #
        try:
            data = np.fromfile(self.dir+'out.cgyro.equilibrium',dtype='float',sep=" ")
            self.rmin        = data[0]
            self.rmaj        = data[1]
            self.q           = data[2]
            self.shear       = data[3]
            self.shift       = data[4]
            self.kappa       = data[5]
            self.s_kappa     = data[6]
            self.delta       = data[7]
            self.s_delta     = data[8]
            self.zeta        = data[9]
            self.s_zeta      = data[10]
            self.zmag        = data[11]
            self.dzmag       = data[12]
            self.rho         = data[13]
            self.ky0         = data[14]
            self.betae_unit  = data[15]
            self.beta_star   = data[16]
            self.lambda_star = data[17]
            self.gamma_e     = data[18]
            self.gamma_p     = data[19]
            self.mach        = data[20]
            self.a_meters    = data[21]
            self.b_unit      = data[22]
            self.dens_norm   = data[23]
            self.temp_norm   = data[24]
            self.vth_norm    = data[25]
            # Define species vectors
            self.z      = np.zeros(self.n_species)
            self.mass   = np.zeros(self.n_species)
            self.dens   = np.zeros(self.n_species)
            self.temp   = np.zeros(self.n_species)
            self.dlnndr = np.zeros(self.n_species)
            self.dlntdr = np.zeros(self.n_species)
            self.nu     = np.zeros(self.n_species)
            for i in range(self.n_species):
                self.z[i]      = data[26+7*i]
                self.mass[i]   = data[27+7*i]
                self.dens[i]   = data[28+7*i]
                self.temp[i]   = data[29+7*i]
                self.dlnndr[i] = data[30+7*i]
                self.dlntdr[i] = data[31+7*i]
                self.nu[i]     = data[32+7*i]
            print "INFO: (data.py) Read data in out.cgyro.equilibrium."
        except:
            print "WARNING: (data.py) Could not read out.cgyro.equilibrium."
            pass

        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Read ballooning potentials
        #
        try:
            data = np.fromfile(self.dir+'out.cgyro.phib',dtype='float',sep=" ")
            if self.n_radial == 1:
                self.phib = np.reshape(data,(2,self.n_theta,nt),'F')
            else:
                self.phib = np.reshape(data,(2,self.n_theta*self.n_radial/self.m_box,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.phib."
        except:
            pass
        try:
            data = np.fromfile(self.dir+'out.cgyro.aparb',dtype='float',sep=" ")
            if self.n_radial == 1:
                self.aparb = np.reshape(data,(2,self.n_theta,nt),'F')
            else:
                self.aparb = np.reshape(data,(2,self.n_theta*self.n_radial/self.m_box,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.aparb."
        except:
            pass
        try:
            data = np.fromfile(self.dir+'out.cgyro.bparb',dtype='float',sep=" ")
            self.bparb = np.reshape(data,(2,self.n_theta*self.n_radial/self.m_box,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.bparb."
        except:
            pass
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Ballooning distribution
        #
        try:
            data = np.fromfile(self.dir+'out.cgyro.hb',dtype='float',sep=" ")
            self.hb = np.reshape(data,(2,self.n_radial*self.n_theta/self.m_box,
                                       self.n_species,self.n_xi,self.n_energy,nt),'F')
            self.hb = self.hb/np.max(self.hb)
            print "INFO: (data.py) Read data in out.cgyro.hb."
        except:
            pass
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        # Compressed particle and energy fluxes
        #
        nd = self.n_species*nt
        try:
            start = time.time()
            data = np.loadtxt(self.dir+'out.cgyro.flux_n',dtype='float')
            end = time.time()
            self.flux_n = np.transpose(data[:,1:])
            print "INFO: (data.py) Read data in out.cgyro.flux_n."
        except:
            pass

        try:
            start = time.time()
            data = np.loadtxt(self.dir+'out.cgyro.flux_e',dtype='float')
            end = time.time()
            self.flux_e = np.transpose(data[:,1:])
            print "INFO: (data.py) Read data in out.cgyro.flux_e."
        except:
            pass 
        #-----------------------------------------------------------------

    def getflux(self):

        import numpy as np
        import time

        # Convenience definition
        nt = self.n_time

        #-----------------------------------------------------------------
        # Particle and energy fluxes
        #
        nd = self.n_species*3*self.n_field*self.n_n*nt

        try:
            start = time.time()
            data = np.fromfile(self.dir+'out.cgyro.ky_flux',dtype='float',sep=" ")
            end = time.time()
            self.ky_flux = np.reshape(data[0:nd],(self.n_species,3,self.n_field,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.ky_flux. TIME = "+str(end-start)
        except:
            pass 

        #-----------------------------------------------------------------

    def getbigflux(self):

        """Larger flux files"""

        import numpy as np
        import time

        # Convenience definition
        nt = self.n_time

        #-----------------------------------------------------------------
        # Particle and energy fluxes
        #
        nd = self.n_radial*self.n_species*self.n_n*nt

        try:
            start = time.time()
            data = np.fromfile(self.dir+'out.cgyro.kxky_flux_e',dtype='float',sep=" ")
            end = time.time()
            self.kxky_flux_e = np.reshape(data[0:nd],(self.n_radial,self.n_species,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.kxky_flux_e. TIME = "+str(end-start)
        except:
            pass

        try:
            start = time.time()
            data = np.fromfile(self.dir+'out.cgyro.kxky_flux_n',dtype='float',sep=" ")
            end = time.time()
            self.kxky_flux_n = np.reshape(data[0:nd],(self.n_radial,self.n_species,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.kxky_flux_n. TIME = "+str(end-start)
        except:
            pass 
        #-----------------------------------------------------------------

    def getxflux(self):

        """Global flux files"""

        import numpy as np
        import time

        # Convenience definition
        nt = self.n_time

        #-----------------------------------------------------------------
        # Particle and energy fluxes
        #
        ng = self.n_global+1
        nd = 2*ng*self.n_species*self.n_n*nt
        try:
            start = time.time()
            data = np.fromfile(self.dir+'out.cgyro.lky_flux_n',dtype='float',sep=" ")
            end = time.time()
            self.lky_flux_n = np.reshape(data[0:nd],(2,ng,self.n_species,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.lky_flux_n. TIME = "+str(end-start)
        except:
            pass

        try:
            start = time.time()
            data = np.fromfile(self.dir+'out.cgyro.lky_flux_e',dtype='float',sep=" ")
            end = time.time()
            self.lky_flux_e = np.reshape(data[0:nd],(2,ng,self.n_species,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.lky_flux_e. TIME = "+str(end-start)
        except:
            pass 

        try:
            start = time.time()
            data = np.fromfile(self.dir+'out.cgyro.lky_flux_v',dtype='float',sep=" ")
            end = time.time()
            self.lky_flux_v = np.reshape(data[0:nd],(2,ng,self.n_species,self.n_n,nt),'F')
            print "INFO: (data.py) Read data in out.cgyro.lky_flux_v. TIME = "+str(end-start)
        except:
            pass 
        #-----------------------------------------------------------------

    def getbigfield(self):

        """Larger field files"""

        import numpy as np
        import time

        # Convenience definition
        nt = self.n_time

        #-----------------------------------------------------------------
        # Read complex fields
        #
        try:
           nd = 2*self.n_radial*self.theta_plot*self.n_n*nt
           start = time.time()
           data = np.fromfile(self.dir+'out.cgyro.kxky_phi',dtype='float',sep=" ")
           end = time.time()
           self.kxky_phi = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_n,nt),'F')
           self.phisq = self.kxky_phi[0,:,:,:,:]**2+self.kxky_phi[1,:,:,:,:]**2
           print "INFO: (data.py) Read data in out.cgyro.kxky_phi. TIME = "+str(end-start)
        except:
           pass
        
        nd = 2*self.n_radial*self.theta_plot*self.n_species*self.n_n*nt

        try:
           start = time.time()
           data = np.fromfile(self.dir+'out.cgyro.kxky_n',dtype='float',sep=" ")
           end = time.time()
           self.n = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_species,self.n_n,nt),'F')
           self.nsq = self.n[0,:,:,:,:,:]**2+self.n[1,:,:,:,:,:]**2
           print "INFO: (data.py) Read data in out.cgyro.kxky_n. TIME = "+str(end-start)
        except:
           pass

        try:
           start = time.time()
           data = np.fromfile(self.dir+'out.cgyro.kxky_e',dtype='float',sep=" ")
           end = time.time()
           self.e = np.reshape(data[0:nd],(2,self.n_radial,self.theta_plot,self.n_species,self.n_n,nt),'F')
           self.esq = self.e[0,:,:,:,:,:]**2+self.e[1,:,:,:,:,:]**2
           print "INFO: (data.py) Read data in out.cgyro.kxky_e. TIME = "+str(end-start)
        except:
           pass
        #-----------------------------------------------------------------

    def getgeo(self):

        """Read geometry arrays"""

        import numpy as np

        # Convenience definition
        nt = self.n_time

        self.geotag = []

        #-----------------------------------------------------------------
        # Read powers
        #
        try:
            data = np.fromfile(self.dir+'out.cgyro.geo',dtype='float',sep=" ")
            self.geo = np.reshape(data,(self.n_theta,12),'F')
            print "INFO: (data.py) Read data in out.cgyro.geo."
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
        except:
            print "INFO: (data.py) Missing out.cgyro.geo."
            pass
        #-----------------------------------------------------------------
