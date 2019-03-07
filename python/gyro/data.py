import os
import numpy as np
import sys

BYTE='float32'

class GYROData:

   """GYRO output data class."""

   def __init__(self, sim_directory):

       """Constructor reads basic (not all) simulation data."""

       self.profile  = {}
       self.geometry = {}
       self.t        = {}
       self.freq     = {}
       self.balloon  = {}
       self.tagspec    = []
       self.tagmom     = []
       self.tagmomtext = []

       if not sim_directory.endswith(os.path.sep):
           sim_directory += os.path.sep
       self.dir = sim_directory
       self.getdata()
       self.read_geometry()
       self.read_units()
       self.make_tags()

       # The rest of the possible read routines should NOT be done automatically.

       """
       self.read_freq()
       self.read_balloon()
       self.read_gbflux_i()
       self.read_gbflux_n()
       self.read_gbflux_exc()
       self.read_moment_u()
       self.read_moment_n()
       self.read_moment_e()
       self.read_moment_v()
       self.read_moment_zero()
       self.read_flux_velocity()
       self.read_k_perp_squared()

       There are also routines to make "optional fluxes"

       self.make_diff()
       self.make_diff_i()
       """

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
      data = np.fromfile(self.dir+'out.gyro.t',dtype='float',sep=' ')
      nt = len(data)/2
      data = np.reshape(data,(2,nt),'F')
      self.t['n_time']    = nt
      self.t['data_step'] = data[0,:]
      self.t['(c_s/a)t']  = data[1,:]
      self.n = nt
      print "INFO: (data.py) Read time vector in out.gyro.t."
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Read grid data, then unpack
      # 
      data = np.fromfile(self.dir+'out.gyro.profile',dtype='float32',sep=' ')

      n_x    = int(data[0])
      n_spec = int(data[15])

      self.profile['n_x']             = int(data[0])
      self.profile['n_theta_section'] = int(data[1])
      self.profile['n_pass']          = int(data[2])
      self.profile['n_trap']          = int(data[3])
      self.profile['n_lambda']        = int(data[2]+data[3])
      self.profile['n_energy']        = int(data[4])
      self.profile['n_theta_plot']    = int(data[5])
      self.profile['n0']              = int(data[6])
      self.profile['n_n']             = int(data[7])
      self.profile['d_n']             = int(data[8])
      self.profile['n_explicit_damp'] = int(data[9])
      self.profile['nonlinear_flag']  = int(data[10])
      self.profile['electron_method'] = int(data[11])
      self.profile['n_field']         = int(data[12])
      self.profile['n_ion']           = int(data[13])
      self.profile['n_kinetic']       = int(data[14])
      self.profile['n_spec']          = int(data[15])
      self.profile['n_grid_exp']      = int(data[16])
      self.profile['boundary_method'] = int(data[17])
      self.profile['r']               = data[18:(18+n_x)]
      self.profile['q']               = data[(18+n_x):(18+2*n_x)]
      self.profile['r_s']             = data[(18+2*n_x):(18+3*n_x)]
      self.profile['q_s']             = data[(18+3*n_x):(18+4*n_x)]
      # The parameter "mark" is used to keep track of where in the file the
      # program is so that the indices don't get too complicated.
      mark = 18 + 4*n_x
      temp = data[mark:(mark+n_spec*n_x)]
      self.profile['dlntdr_s'] = temp.reshape((n_spec,n_x),order='F')
      mark = mark + n_spec*n_x
      temp = data[mark:(mark+n_spec*n_x)]
      self.profile['dlnndr_s'] = temp.reshape((n_spec,n_x),order='F')
      mark = mark + n_spec*n_x
      temp = data[mark:(mark+n_spec*n_x)]
      self.profile['tem_s'] = temp.reshape((n_spec,n_x),order='F')
      mark = mark + n_spec*n_x
      temp = data[mark:(mark+n_spec*n_x)]
      self.profile['den_s'] = temp.reshape((n_spec,n_x),order='F')
      mark = mark + n_spec*n_x
      self.profile['rmaj_s/r_s'] = data[mark:(mark+n_x)]
      self.profile['delta_s']    = data[(mark+n_x):(mark+2*n_x)]
      self.profile['zeta_s']     = data[(mark+2*n_x):(mark+3*n_x)]
      self.profile['kappa_s']    = data[(mark+3*n_x):(mark+4*n_x)]
      self.profile['drmaj_s']    = data[(mark+4*n_x):(mark+5*n_x)]
      self.profile['shat_s']     = data[(mark+5*n_x):(mark+6*n_x)]
      self.profile['s_delta_s']  = data[(mark+6*n_x):(mark+7*n_x)]
      self.profile['s_zeta_s']   = data[(mark+7*n_x):(mark+8*n_x)]
      self.profile['s_kappa_s']  = data[(mark+8*n_x):(mark+9*n_x)]
      self.profile['zmag_s']     = data[(mark+9*n_x):(mark+10*n_x)]
      self.profile['dzmag_s']    = data[(mark+10*n_x):(mark+11*n_x)]
      self.profile['beta_unit_s'] = data[(mark+11*n_x):(mark+12*n_x)]
      self.profile['gamma_e_s']  = data[(mark+12*n_x):(mark+13*n_x)]
      self.profile['gamma_p_s']  = data[(mark+13*n_x):(mark+14*n_x)]
      self.profile['mach_s']     = data[(mark+14*n_x):(mark+15*n_x)]
      self.profile['b_unit_s']   = data[(mark+15*n_x):(mark+16*n_x)]
      self.profile['z_eff_s']    = data[(mark+16*n_x):(mark+17*n_x)]
      self.profile['nu_s']       = data[(mark+17*n_x):(mark+18*n_x)]
      self.profile['w0_s']       = data[(mark+18*n_x):(mark+19*n_x)]
      mark = mark + 19*n_x
      self.profile['box_multiplier'] = data[mark]
      self.profile['lambda']     = data[(mark+1):(mark+1+self.profile['n_lambda'])]
      mark = mark + 1 + self.profile['n_lambda']
      self.profile['energy']     = data[mark:(mark+self.profile['n_energy'])]
      self.profile['lambda_tp']  = data[mark + self.profile['n_energy']]
      mark = mark + self.profile['n_energy'] + 1
      self.profile['kt_rho']     = data[mark:(mark+self.profile['n_n'])]
      self.profile['rho_s']      = data[mark+self.profile['n_n']]
      mark = mark + self.profile['n_n'] + 1
      self.profile['z']          = data[mark:(mark+n_spec)]
      mark = mark + n_spec 
      self.profile['mu']         = data[mark:(mark+n_spec)]
      self.profile['n_fine']     = int(data[mark+n_spec])
      # Done
      
      if self.profile['n_theta_plot'] == 1:
          self.profile['theta_plot'] = 0.0
      else:
          n_theta = self.profile['n_theta_plot']
          self.profile['theta_plot'] = -np.pi+2*np.pi*np.arange(n_theta)/float(n_theta)

      print "INFO: (data.py) Read grid data in out.gyro.profile."
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # Linear frequency
      #
      nd = 4*self.profile['n_n']*nt
      t,fmt,data = self.extract('.gyro.freq')
      if fmt != 'null':
         temp = np.reshape(data[0:nd],(4,self.profile['n_n'],nt),'F')
         self.freq['(a/c_s)w']        = temp[0,:,:]
         self.freq['(a/c_s)gamma']    = temp[1,:,:]
         self.freq['err(a/c_s)w']     = temp[2,:,:]
         self.freq['err(a/c_s)gamma'] = temp[3,:,:]
         print "INFO: (data.py) Read data in "+fmt+".gyro.freq. "+t 
      #-----------------------------------------------------------------

      #-----------------------------------------------------------------
      # RMS field
      #
      nd = 2*nt
      t,fmt,data = self.extract('.gyro.field_rms')
      if fmt != 'null':
         self.field_rms = np.reshape(data[0:nd],(2,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".gyro.field_rms. "+t 
      #-----------------------------------------------------------------
   

   def read_units(self):
      """Read out.gyro.units ."""

      try:
          f = open(self.dir+'/out.gyro.units')
      except IOError:
          raise IOError('ERROR (GYROData): Fatal error!  Missing out.gyro.units.')

      fl = f.readlines()
      f.close()
      units=[]
      for fi in fl[0:13]:
          units.append(float(fi.strip().split()[0]))
      self.units=units


   def read_geometry(self):

        try:
            geometry = np.fromfile(self.dir+'/out.gyro.geometry_arrays',dtype=float,sep=" ")
        except IOError:
            raise IOError("ERROR (GYROData): out.gyro.geometry_arrays not found.")

        temp = geometry.reshape((12, self.profile['n_fine'], self.profile['n_x']), order='F')

        self.geometry['v']       = temp[0,:,:]
        self.geometry['gsin']    = temp[1,:,:]
        self.geometry['gcos1']   = temp[2,:,:]
        self.geometry['gcos2']   = temp[3,:,:]
        self.geometry['usin']    = temp[4,:,:]
        self.geometry['ucos']    = temp[5,:,:]
        self.geometry['B']       = temp[6,:,:]
        self.geometry['G_theta'] = temp[7,:,:]
        self.geometry['grad_r']  = temp[8,:,:]
        self.geometry['G_q']     = temp[9,:,:]
        self.geometry['THETA']   = temp[10,:,:]
        self.geometry['theta_nc'] = temp[11,:,:]

    #---------------------------------------------------------------------------#

   def read_gbflux_i(self):

      #-----------------------------------------------------------------
      # Radial profiles of fluxes
      #
      nt        = self.n
      n_kinetic = self.profile['n_kinetic']
      n_field   = self.profile['n_field']
      n_x       = self.profile['n_x']
      nd        = n_kinetic*n_field*4*n_x*nt

      t,fmt,data = self.extract('.gyro.gbflux_i')
      if fmt != 'null':
         self.gbflux_i = np.reshape(data[0:nd],(n_kinetic,n_field,4,n_x,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".gyro.gbflux_i. "+t
         # NOTE: The average below used to occur in GYRO, but not it's done here.
         if self.profile['boundary_method'] == 1:
            # Periodic simulation
            self.gbflux = np.mean(self.gbflux_i,axis=3)
         else:
            # Nonperiodic simulation: don't include buffers in average
            n = self.profile['n_explicit_damp']
            self.gbflux = np.mean(self.gbflux_i[:,:,:,n:-n,:],axis=3)
      #---------------------------------------------------------------------------#

   def read_gbflux_n(self):

      #-----------------------------------------------------------------
      # ky-dependent fluxes
      #
      nt        = self.n
      n_kinetic = self.profile['n_kinetic']
      n_field   = self.profile['n_field']
      n_n       = self.profile['n_n']
      nd        = n_kinetic*n_field*4*n_n*nt

      t,fmt,data = self.extract('.gyro.gbflux_n')
      if fmt != 'null':
         self.gbflux_n = np.reshape(data[0:nd],(n_kinetic,n_field,4,n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".gyro.gbflux_n. "+t
      #---------------------------------------------------------------------------#

   def read_gbflux_exc(self):

      #-----------------------------------------------------------------
      # exchanges
      #
      nt        = self.n
      n_kinetic = self.profile['n_kinetic']
      nd        = n_kinetic*2*nt

      t,fmt,data = self.extract('.gyro.gbflux_exc')
      if fmt != 'null':
         self.gbflux_exc = np.reshape(data[0:nd],(n_kinetic,2,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".gyro.gbflux_exc. "+t
      #---------------------------------------------------------------------------#

   def read_kxkyspec(self):

        try:
            kxkyspec = np.fromfile(self.dir+'/out.gyro.kxkyspec',dtype=float,sep=" ")
        except:
            raise IOError("ERROR (GYROData): out.gyro.kxkyspec not found.")

        n_x = self.profile['n_x']
        n_n = self.profile['n_n']
        nt  = len(kxkyspec)/(n_x*n_n)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dir+
              '/out.gyro.kxkyspec too small. ')

        self.kxkyspec = kxkyspec.reshape((n_x,n_n,nt),order='F')


    #---------------------------------------------------------------------------#

   def read_moment(self,indx):

      """Read moment_u (field) data"""

      nt           = self.n
      n_theta_plot = self.profile['n_theta_plot']
      n_x          = self.profile['n_x']
      n_field      = self.profile['n_field']
      n_n          = self.profile['n_n']
      n_kinetic    = self.profile['n_kinetic']

        
      t,fmt,data = self.extract('.gyro.moment_'+indx)

      if fmt != 'null' and indx == 'u':
         nd = 2*n_theta_plot*n_x*n_field*n_n*nt
         self.moment_u = np.reshape(data[0:nd],(2,n_theta_plot,n_x,n_field,n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".gyro.moment_u. "+t
         self.moment_u = self.moment_u[0] + 1j*self.moment_u[1]

      nd = 2*n_theta_plot*n_x*n_kinetic*n_n*nt

      if fmt != 'null' and indx == 'n':
         self.moment_n = np.reshape(data[0:nd],(2,n_theta_plot,n_x,n_kinetic,n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".gyro.moment_n. "+t
         self.moment_n = self.moment_n[0] + 1j*self.moment_n[1]

      if fmt != 'null' and indx == 'e':
         self.moment_e = np.reshape(data[0:nd],(2,n_theta_plot,n_x,n_kinetic,n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".gyro.moment_e. "+t
         self.moment_e = self.moment_e[0] + 1j*self.moment_e[1]

      if fmt != 'null' and indx == 'v':
         self.moment_v = np.reshape(data[0:nd],(2,n_theta_plot,n_x,n_kinetic,n_n,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".gyro.moment_v. "+t
         self.moment_v = self.moment_v[0] + 1j*self.moment_v[1]


   def read_moment_zero(self):

      nt   = self.n
      n_x  = self.profile['n_x']
      ns   = self.profile['n_kinetic']
      nd   = n_x*ns*3*nt
      
      t,fmt,data = self.extract('.gyro.moment_zero')
      if fmt == 'bin':
         self.moment_zero = np.reshape(data[0:2*nd],(n_x,ns,6,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".gyro.moment_zero. "+t
      if fmt == 'out':
         self.moment_zero = np.reshape(data[0:nd],(n_x,ns,3,nt),'F')
         print "INFO: (data.py) Read data in "+fmt+".gyro.moment_zero. "+t
      
   def read_balloon(self):
        """Reads out.gyro.balloon*.  Data is stored in self.balloon"""

        import glob
        import string

        m     = int(self.profile['box_multiplier'])
        n_x   = int(self.profile['n_x'])
        n_ang = int(self.profile['n_theta_plot']*n_x/m)

        list = glob.glob(self.dir+'/out.gyro.balloon*')

        # If list is empty, then exit with error message.
        if len(list) == 0:
            raise IOError("ERROR (GYROData): out.gyro.balloon* not found.")

        for filename in list:
            data = np.fromfile(filename,dtype=float,sep=" ")
            u = data.reshape((2,n_ang,m,self.t['n_time']),order='F')
            ext = string.splitfields(filename,'.')[-1]
            self.balloon[ext] = u[0,...]+1j*u[1,...]
            

   def make_tags(self):
        """Generate tags for fields, moments and species"""

        n_kinetic = self.profile['n_kinetic']
        melec     = self.profile['electron_method']

        # Species tags

        for i in range(n_kinetic):

            if melec == 2 or melec == 4:
                if i == n_kinetic-1:
                    self.tagspec.append('elec')
                else:
                    self.tagspec.append('ion-'+str(i+1))

            if melec == 1:
                self.tagspec.append('ion-'+str(i+1))

            if melec == 3:
                self.tagspec.append('elec')

       # Moment tags

        self.tagmomtext = ['GAMMA [GB]','Q [GB]','PI [GB]','S [GB]']

        self.tagmom = ['\Gamma/\Gamma_\mathrm{GB}',
                       'Q/Q_\mathrm{GB}',
                       '\Pi/\Pi_\mathrm{GB}',
                       'S/S_\mathrm{GB}']

       # Field tags
        self.tagfieldtext = ['Phi',
                             'Apar',
                             'Bpar',
                             'Tot']

        self.tagfield = ['\mathrm{electrostatic}',
                         '\mathrm{flutter}',
                         '\mathrm{compression}',
                         '\mathrm{total}']

   def get_sim_exp_flux(self,field_num=None,ft1=.5,ft2=1.,moment=1,
          GB_norm=0):
        """
        Return the radially averaged, time averaged simulated flux and
        experimental flux for the radius at the box center.

        --------
        Keywords

        field_num:  0 - Electrostatic field part
                    1 - Flutter
                    2 - Compression
                    None - Sum (Default)
        ft1,ft2:    The fraction of time over which to average
                    (0.5-1.0 default)
        moment:     0 - Particle flux
                    1 - Thermal flux (Default)
                    2 - Momentum flux
                    3 - Exchange flux
        GB_norm:    0 - No normalization [unit/m^2] (Default)
                    1 - Normalized to gyrobohm flux []
                    2 - Flow [unit] # Not implemented, need input.profiles.extra
        """
        import math
        
        gbflux = self.gbflux
        if field_num==None:
          gbflux = np.sum(gbflux,axis=1)
        else:
          gbflux = gbflux[:,field_num,:,:]
        gbflux = gbflux[:,moment,:]
        ind_t1 = int(math.ceil((self.t['n_time']-1)*ft1))
        ind_t2 = int(math.floor((self.t['n_time']-1)*ft2))
        if ind_t2>len(gbflux[0,:]):
          ind_t2 = len(gbflux[0,:])
        if ind_t1<0:
          ind_t1 = 0
        gbflux_m = np.mean(gbflux[:,ind_t1:ind_t2],axis=1) #Mean
        gbflux_u = np.std(gbflux[:,ind_t1:ind_t2],axis=1)  #Standard Deviation
        ir_norm = int(self.profile['n_x']/2+1)
        r = self.profile['r'][ir_norm]
        rmin = self.input_profile.data['rmin']
        rmin = rmin/max(rmin)
        ind_r = np.argmin(abs(rmin-r))
        surface_area = self.exp_derived[22,ind_r]
        #print 'r=',r,'ir_norm=',ir_norm,'ind_r=',ind_r
        #print 'rmin=',rmin[ind_r],'rho=',self.input_profile.data['rho'][ind_r]
        #print 'surface_area=',surface_area
        fac_gyro = 1.
        fac_xp = 1.
        if GB_norm == 0:
          fac_xp = 1./surface_area #Convert MW to MW/m^2
          fac_gyro = self.units[moment+9] # Convert MW/m^2/(MW/m^2) to MW/m^2
        elif GB_norm == 1:
          fac_xp = 1./surface_area/self.units[moment+9] #Convert MW to MW/m^2/(MW/m^2)
          fac_gyro = 1.# Already in GB units
        elif GB_norm == 2:
          fac_xp = 1. # Already in MW
          fac_gyro = self.units[moment+9]*surface_area #Convert MW/m^2/(MW/m^2) to MW
        if moment==1:
          xp = {'e_exp':self.input_profile.data['pow_e'][ind_r]*fac_xp,
                'i_exp':self.input_profile.data['pow_i'][ind_r]*fac_xp}
        else:
          raise NotImplemented('Only the energy moment is currently treated')
        sim = {}
        for ti,tg in enumerate(self.tagspec):
            sim[tg] = gbflux_m[ti]*fac_gyro
            sim[tg+'_uncertainty'] = gbflux_u[ti]*fac_gyro
        return (xp,sim)

