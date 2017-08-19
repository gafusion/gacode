import os
import numpy as np
class GYROData:
    """GYRO output data class.

    Data:
    
    dirname = ""
    profile  = {}
    geometry = {}
    t        = {}
    freq     = {}
    balloon  = []
    diff     = []
    diff_i   = []
    diff_n   = []
    gbflux   = []
    gbflux_i = []
    gbflux_n = []
    moment_u = []
    moment_n = []
    moment_e = []
    moment_v = []
    moment_zero    = []
    flux_velocity  = []
    k_perp_squared = []

    Example Usage:
        >>> from gyro.data import GYROData
        >>> sim = GYROData('example_directory')
        >>> sim.make_gbflux()
    """

    #---------------------------------------------------------------------------#
    # Methods

    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object.
    def __init__(self, sim_directory, create_netcdf=1e6):
        """
        Constructor reads in data from sim_directory and creates new object.
        Keywords:
        create_netcdf - An integer, which if the number of time steps is greater
          than this number, then the method convert_to_netcdf will be run
        """

        self.init_data()
        self.set_directory(sim_directory)
        self.read_profile()
        self.read_geometry()
        self.read_units()
        self.read_input_profiles()
        self.make_tags()
        try:
            self.read_t()
        except:
            print 'Warning: Error reading %s/out.gyro.t'%(self.simdir)
        else:
            if create_netcdf and self.n>1:
                self.convert_to_netcdf()
        # The rest of the possible read routines should NOT be done automatically.

        #self.read_freq()
        #self.read_balloon()
        #self.read_gbflux_i()
        #self.read_gbflux_n()
        #self.read_gbflux_exc()
        #self.read_moment_u()
        #self.read_moment_n()
        #self.read_moment_e()
        #self.read_moment_v()
        #self.read_moment_zero()
        #self.read_flux_velocity()
        #self.read_k_perp_squared()

        # There are also routines to make "optional fluxes"

        #self.make_gbflux()
        #self.make_diff()
        #self.make_diff_i()
        from os.path import expanduser, expandvars
        import sys
        path = expanduser(expandvars('$GACODE_ROOT/shared/python/pyrats'))
        if path not in sys.path:
          sys.path.append(path)

    #---------------------------------------------------------------------------#

    def init_data(self):
        """Initialize object data."""

        self.dirname = ""

        self.profile  = {}
        self.geometry = {}
        self.t        = {}
        self.freq     = {}
        self.balloon  = {}

        self.diff           = []
        self.diff_i         = []
        self.diff_n         = []
        self.gbflux         = []
        self.gbflux_i       = []
        self.gbflux_n       = []
        self.gbflux_exc     = []
        self.moment_u       = []
        self.moment_n       = []
        self.moment_e       = []
        self.moment_v       = []
        self.moment_zero    = []
        self.flux_velocity  = []
        self.k_perp_squared = []
        self.loaded         = []

        self.tagspec    = []
        self.tagmom     = []
        self.tagmomtext = []
        self.units = []
        
        self.xp_sim_flux = {}
 
    #---------------------------------------------------------------------------#
        
    def set_directory(self, path):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        self.dirname = os.path.realpath(expanduser(expandvars(path)))

    #---------------------------------------------------------------------------#

    def read_t(self):
        """Read out.gyro.t to get time data."""

        try:
            t = np.loadtxt(self.dirname + '/out.gyro.t')
        except IOError:
            raise IOError("ERROR (GYROData): Fatal error!  Missing out.gyro.t.")

        self.t['n_time']    = len(t[:,0])
        self.t['data_step'] = t[:,0]
        self.t['(c_s/a)t']  = t[:,1]
        self.loaded.append('t')
        self.n = len(t[:,0])
    
    #---------------------------------------------------------------------------#
    
    def read_units(self):
        """Read out.gyro.units ."""
        
        try:
          f = open(self.dirname+'/out.gyro.units')
        except IOError:
          raise IOError('ERROR (GYROData): Fatal error!  Missing out.gyro.units.')
 
        fl = f.readlines()
        f.close()
        units=[]
        for fi in fl[0:13]:
            units.append(float(fi.strip().split()[0]))
        self.units=units
        self.loaded.append('units')
        
    #---------------------------------------------------------------------------#
    
    def read_input_profiles(self):
        """
        Read input.profiles, input.profiles.extra if they exist.
        Initializes self.input_profiles and self.exp_derived, respectively.
        """
        
        fn = self.dirname+'/input.profiles'
        if os.path.exists(fn):
            import pyrats.profiles_gen.data
            reload(pyrats.profiles_gen.data)
            self.input_profile = pyrats.profiles_gen.data.profiles_genData(fn)
            self.loaded.append('input.profiles')
            fn = fn+'.extra'
            if os.path.exists(fn):
                num = 25*self.input_profile.n_exp
                f = open(fn)
                fl=f.read()
                f.close()
                extra=np.array(fl.split()[0:num],dtype='Float64')
                self.exp_derived = extra.reshape(25,-1)
                self.loaded.append('input.profiles.extra')
        
    #---------------------------------------------------------------------------#

    def read_freq(self):
        """Reads frequency data.  Output is dictionary of numpy arrays with
        dimensions: n_n x n_time"""

        try:
            freq = np.loadtxt(self.dirname+'/out.gyro.freq').transpose()
        except IOError:
            raise IOError("ERROR (GYROData): Missing out.gyro.freq.")

        temp = freq.reshape((4,self.profile['n_n'],self.t['n_time']), order='F')
        
        self.freq['(a/c_s)w']        = temp[0,:,:]
        self.freq['(a/c_s)gamma']    = temp[1,:,:]
        self.freq['err(a/c_s)w']     = temp[2,:,:]
        self.freq['err(a/c_s)gamma'] = temp[3,:,:]    

        self.loaded.append('freq')

    #---------------------------------------------------------------------------#

    def read_profile(self):
        """Read out.gyro.profile to get control data.  Output is dictionary
        containing necessary information.
        """


        try:
            profile = np.loadtxt(self.dirname+'/out.gyro.profile')
        except IOError:
            raise IOError("ERROR (GYROData): Fatal error!  Missing out.gyro.profile.")
                
        n_x    = int(profile[0]) 
        n_spec = int(profile[15])

        self.profile['n_x']             = int(profile[0])
        self.profile['n_theta_section'] = int(profile[1])
        self.profile['n_pass']          = int(profile[2])
        self.profile['n_trap']          = int(profile[3])
        self.profile['n_lambda']        = int(profile[2]+profile[3])
        self.profile['n_energy']        = int(profile[4])
        self.profile['n_theta_plot']    = int(profile[5])
        self.profile['n0']              = int(profile[6])
        self.profile['n_n']             = int(profile[7])
        self.profile['d_n']             = int(profile[8])
        self.profile['n_explicit_damp'] = int(profile[9])
        self.profile['nonlinear_flag']  = int(profile[10])
        self.profile['electron_method'] = int(profile[11])
        self.profile['n_field']         = int(profile[12])
        self.profile['n_ion']           = int(profile[13])
        self.profile['n_kinetic']       = int(profile[14])
        self.profile['n_spec']          = int(profile[15])
        self.profile['n_grid_exp']      = int(profile[16])
        self.profile['boundary_method'] = int(profile[17])
        self.profile['r']               = profile[18:(18+n_x)]
        self.profile['q']               = profile[(18+n_x):(18+2*n_x)]
        self.profile['r_s']             = profile[(18+2*n_x):(18+3*n_x)]
        self.profile['q_s']             = profile[(18+3*n_x):(18+4*n_x)]
        # The parameter "mark" is used to keep track of where in the file the
        # program is so that the indices don't get too complicated.
        mark = 18 + 4*n_x
        temp = profile[mark:(mark+n_spec*n_x)]
        self.profile['dlntdr_s'] = temp.reshape((n_spec,n_x),order='F')
        mark = mark + n_spec*n_x
        temp = profile[mark:(mark+n_spec*n_x)]
        self.profile['dlnndr_s'] = temp.reshape((n_spec,n_x),order='F')
        mark = mark + n_spec*n_x
        temp = profile[mark:(mark+n_spec*n_x)]
        self.profile['tem_s'] = temp.reshape((n_spec,n_x),order='F')
        mark = mark + n_spec*n_x
        temp = profile[mark:(mark+n_spec*n_x)]
        self.profile['den_s'] = temp.reshape((n_spec,n_x),order='F')
        mark = mark + n_spec*n_x
        self.profile['rmaj_s/r_s'] = profile[mark:(mark+n_x)]
        self.profile['delta_s']    = profile[(mark+n_x):(mark+2*n_x)]
        self.profile['zeta_s']     = profile[(mark+2*n_x):(mark+3*n_x)]
        self.profile['kappa_s']    = profile[(mark+3*n_x):(mark+4*n_x)]
        self.profile['drmaj_s']    = profile[(mark+4*n_x):(mark+5*n_x)]
        self.profile['shat_s']     = profile[(mark+5*n_x):(mark+6*n_x)]
        self.profile['s_delta_s']  = profile[(mark+6*n_x):(mark+7*n_x)]
        self.profile['s_zeta_s']   = profile[(mark+7*n_x):(mark+8*n_x)]
        self.profile['s_kappa_s']  = profile[(mark+8*n_x):(mark+9*n_x)]
        self.profile['zmag_s']     = profile[(mark+9*n_x):(mark+10*n_x)]
        self.profile['dzmag_s']    = profile[(mark+10*n_x):(mark+11*n_x)]
        self.profile['beta_unit_s'] = profile[(mark+11*n_x):(mark+12*n_x)]
        self.profile['gamma_e_s']  = profile[(mark+12*n_x):(mark+13*n_x)]
        self.profile['gamma_p_s']  = profile[(mark+13*n_x):(mark+14*n_x)]
        self.profile['mach_s']     = profile[(mark+14*n_x):(mark+15*n_x)]
        self.profile['b_unit_s']   = profile[(mark+15*n_x):(mark+16*n_x)]
        self.profile['dr_eodr']    = profile[(mark+16*n_x):(mark+17*n_x)]
        self.profile['z_eff_s']    = profile[(mark+17*n_x):(mark+18*n_x)]
        self.profile['nu_s']       = profile[(mark+18*n_x):(mark+19*n_x)]
        self.profile['w0_s']       = profile[(mark+19*n_x):(mark+20*n_x)]
        mark = mark + 20*n_x
        self.profile['box_multiplier'] = profile[mark]
        self.profile['lambda']     = profile[(mark+1):(mark+1+self.profile['n_lambda'])]
        mark = mark + 1 + self.profile['n_lambda']
        self.profile['energy']     = profile[mark:(mark+self.profile['n_energy'])]
        self.profile['lambda_tp']  = profile[mark + self.profile['n_energy']]
        mark = mark + self.profile['n_energy'] + 1
        self.profile['kt_rho']     = profile[mark:(mark+self.profile['n_n'])]
        self.profile['rho_s']      = profile[mark+self.profile['n_n']]
        mark = mark + self.profile['n_n'] + 1
        self.profile['z']          = profile[mark:(mark+n_spec)]
        self.profile['n_fine']     = int(profile[mark+n_spec])
        self.profile['n_moment']   = int(profile[mark+n_spec+1])
        
        if self.profile['n_theta_plot'] == 1:
            self.profile['theta_plot'] = 0.0
        else:
            n_theta = self.profile['n_theta_plot']
            self.profile['theta_plot'] = -np.pi+2*np.pi*np.arange(n_theta)/float(n_theta)


    #---------------------------------------------------------------------------#

    def read_geometry(self):
        """Reads in geometry_array data.  Output is dictionary of numpy arrays
        with dimensions: n_fine x n_x."""        


        try:
            geometry = np.fromfile(self.dirname+'/out.gyro.geometry_arrays',dtype=float,sep=" ")
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

        self.loaded.append('geometry')

    #---------------------------------------------------------------------------#

    def read_gbflux_i(self):
        """Reads in gbflux_i data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_x x n_time"""

        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_x       = self.profile['n_x']
        self.deal_with_large_exponents('out.gyro.gbflux_i')
        
        #try:
        gbflux_i = np.fromfile(self.dirname+'/out.gyro.gbflux_i',dtype=float,sep=" ")
        #except:
            #raise #IOError("ERROR (GYROData): out.gyro.gbflux_i not found.")
            

        nt = len(gbflux_i)/(n_kinetic*n_field*4*n_x)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dirname+
              '/out.gyro.gbflux_i too small. ')            
        
        self.gbflux_i = gbflux_i.reshape((n_kinetic,n_field,4,n_x,nt),order='F')
        self.gbflux_i = self.gbflux_i[:,:,:,:,:self.n]
        self.loaded.append('gbflux_i')

    #---------------------------------------------------------------------------#

    def read_gbflux_n(self):
        """Reads gbflux_n data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_n x n_time"""
        
        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_n       = self.profile['n_n']
        
        self.deal_with_large_exponents('out.gyro.gbflux_n')

        #try:
        gbflux_n = np.fromfile(self.dirname+'/out.gyro.gbflux_n',dtype=float,sep=" ")
        #except:
            #raise IOError("ERROR (GYROData): out.gyro.gbflux_n not found.")            

        nt = len(gbflux_n)/(n_kinetic*n_field*4*n_n)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dirname+
              '/out.gyro.gbflux_n too small. ')
        
        self.gbflux_n = gbflux_n.reshape((n_kinetic,n_field,4,n_n,nt),order='F')
        self.loaded.append('gbflux_n')

     #---------------------------------------------------------------------------#

    def read_gbflux_exc(self):
        """Reads gbflux_exc data.  Output is numpy array with dimensions:
        n_kinetic x 4 x n_time"""

        n_kinetic = self.profile['n_kinetic']
        
        self.deal_with_large_exponents('out.gyro.gbflux_exc')

        #try:
        gbflux_exc = np.fromfile(self.dirname+'/out.gyro.gbflux_exc',dtype=float,sep=" ")
        #except:
        #    raise IOError("ERROR (GYROData): out.gyro.gbflux_exc not found.")
            
        nt = len(gbflux_exc)/(n_kinetic*4)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dirname+
              '/out.gyro.gbflux_exc too small. ')
            
        
        self.gbflux_exc = gbflux_exc.reshape((n_kinetic,4,nt),order='F')
        self.loaded.append('gbflux_exc')

        #---------------------------------------------------------------------------#

    def read_kxkyspec(self):

        try:
            kxkyspec = np.fromfile(self.dirname+'/out.gyro.kxkyspec',dtype=float,sep=" ")
        except:
            raise IOError("ERROR (GYROData): out.gyro.kxkyspec not found.")
            
        n_x = self.profile['n_x']
        n_n = self.profile['n_n']
        nt  = len(kxkyspec)/(n_x*n_n)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dirname+
              '/out.gyro.kxkyspec too small. ')
            
        self.kxkyspec = kxkyspec.reshape((n_x,n_n,nt),order='F')
        self.loaded.append('kxkyspec')

    #---------------------------------------------------------------------------#
    def read_field_rms(self):

        nt = self.n

        try:
            field_rms = np.fromfile(self.dirname+'/out.gyro.field_rms',dtype=float,sep=" ")
        except:
            raise IOError("ERROR (GYROData): out.gyro.field_rms not found.")
            
        self.field_rms = field_rms.reshape((2,nt),order='F')
        self.loaded.append('field_rms')

    #---------------------------------------------------------------------------#

    def read_moment_u(self):
        """Reads in moment_u data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_field x n_n x n_time"""

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_field      = self.profile['n_field']
        n_n          = self.profile['n_n']

        try:
            data = np.fromfile(self.dirname+'/out.gyro.moment_u',dtype=float,sep=" ")
        except IOError:
            raise IOError("ERROR (GYROData): out.gyro.moment_u not found.")

        nt = len(data)/(2*n_theta_plot*n_x*n_field*n_n)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dirname+
                          '/out.gyro.moment_u too small. ')

        self.moment_u = data.reshape((2,n_theta_plot,n_x,n_field,n_n,nt),order='F')
        self.moment_u = self.moment_u[0] + 1j*self.moment_u[1]

        self.loaded.append('moment_u')
        
    #---------------------------------------------------------------------------#

    def read_moment_n(self):
        """Reads in moment_n data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_n          = self.profile['n_n']

        try:
            data = np.fromfile(self.dirname+'/out.gyro.moment_n',dtype=float,sep=" ")
        except:
            raise IOError("ERROR (GYROData): out.gyro.moment_n not found.")

        nt = len(data)/(2*n_theta_plot*n_x*n_kinetic*n_n)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dirname+
                          '/out.gyro.moment_n too small. ')

        self.moment_n = data.reshape((2,n_theta_plot,n_x,n_kinetic,n_n,nt),order='F')
        self.moment_n = self.moment_n[0] + 1j*self.moment_n[1]

        self.loaded.append('moment_n')

    #---------------------------------------------------------------------------#

    def read_moment_e(self):
        """Reads in moment_e data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_field x n_n x n_time"""

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_n          = self.profile['n_n']

        try:
            data = np.fromfile(self.dirname+'/out.gyro.moment_e',dtype=float,sep=" ")
        except:
            raise IOError("ERROR (GYROData): out.gyro.moment_e not found.")

        nt = len(data)/(2*n_theta_plot*n_x*n_kinetic*n_n)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dirname+
                          '/out.gyro.moment_e too small. ')


        self.moment_e = data.reshape((2,n_theta_plot,n_x,n_kinetic,n_n,nt),order='F')
        self.moment_e = self.moment_e[0] + 1j*self.moment_e[1]

        self.loaded.append('moment_e')
        
    #---------------------------------------------------------------------------#

    def read_moment_v(self):
        """Reads in moment_v data.  Output is numpy array with dimensions:
        n_theta_plot x n_x x n_kinetic x n_n x n_time"""

        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_n          = self.profile['n_n']

        try:
            data = np.fromfile(self.dirname+'/out.gyro.moment_v',dtype=float,sep=" ")
        except:
            raise IOError("ERROR (GYROData): out.gyro.moment_v not found.")

        nt = len(data)/(2*n_theta_plot*n_x*n_kinetic*n_n)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dirname+
                          '/out.gyro.moment_e too small. ')

        self.moment_v = data.reshape((2,n_theta_plot,n_x,n_kinetic,n_n,nt),order='F')
        self.moment_v = self.moment_v[0] + 1j*self.moment_v[1]

        self.loaded.append('moment_v')        

    #---------------------------------------------------------------------------#

    def read_moment_zero(self):
        """Read data in out.gyro.moment_zero, store in self.moment_zero. 
        Dimensions: (n_x,n_kinetic,n_moment,n_time)"""

        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_moment     = self.profile['n_moment']

        try:
            data = np.fromfile(self.dirname+'/out.gyro.moment_zero',dtype=float,sep=" ")
        except:
            raise IOError("ERROR (GYROData): out.gyro.moment_zero not found.")

        t = len(data)/(n_x*n_kinetic*n_moment)
        self.moment_zero = data.reshape((n_x,n_kinetic,n_moment,t),order='F')

        self.loaded.append('moment_zero')

    #---------------------------------------------------------------------------#

    def read_flux_velocity(self):
        """Reads out.gyro.flux_velocity.  
        Output is numpy array with dimensions: 
          (n_energy,n_lambda,n_kinetic,n_field,2,n_n,n_time)"""

        n_energy  = self.profile['n_energy']
        n_lambda  = self.profile['n_lambda']
        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_n       = self.profile['n_n']

        try:
            data = np.fromfile(self.dirname+'/out.gyro.flux_velocity',dtype=float,sep=" ")
        except:
            raise IOError("ERROR (GYROData): out.gyro.flux_velocity not found.")

        t = len(data)/(n_energy*n_lambda*n_kinetic*n_field*2*n_n)
        self.flux_velocity = data.reshape(
            (n_energy,n_lambda,n_kinetic,n_field,2,n_n,t),order='F')

        self.loaded.append('flux_velocity')

    #---------------------------------------------------------------------------#

    def read_k_perp_squared(self):
        """Reads out.gyro.k_perp_squared.  
           Output is numpy array with dimensions: (n_n,n_time)"""

        try:
            data = np.fromfile(self.dirname+'/out.gyro.k_perp_squared',dtype=float,sep=" ")
        except:
            raise IOError("ERROR (GYROData): out.gyro.kperp_squared not found.")

        t = len(data)/self.profile['n_n']
        self.k_perp_squared = data.reshape((self.profile['n_n'],t),order='F')

        self.loaded.append('k_perp_squared')

    #---------------------------------------------------------------------------#

    def read_balloon(self):
        """Reads out.gyro.balloon*.  Data is stored in self.balloon"""

        import glob
        import string
        
        m     = self.profile['box_multiplier']
        n_x   = self.profile['n_x']
        n_ang = self.profile['n_theta_plot']*n_x/m

        list = glob.glob(self.dirname+'/out.gyro.balloon*')

        # If list is empty, then exit with error message.
        if len(list) == 0:
            raise IOError("ERROR (GYROData): out.gyro.balloon* not found.")

        for filename in list:
            data = np.fromfile(filename,dtype=float,sep=" ")
            u = data.reshape((2,n_ang,m,self.t['n_time']),order='F')
            ext = string.splitfields(filename,'.')[-1]
            self.balloon[ext] = u[0,...]+1j*u[1,...]
 
    #------------------------------------------------------------
    # Create data from other previously imported data

    def make_gbflux(self):
        """Makes gbflux (omitting buffers properly). Output is numpy array
           with dimensions: n_kinetic x n_field x 4 x n_time"""
        
        if len(self.gbflux_i)==0:
            self.read_gbflux_i()
        if self.profile['boundary_method'] == 1:
            # Periodic simulation
            self.gbflux = np.mean(self.gbflux_i, axis=3)
        else:
            # Nonperiodic simulation: don't include buffers in average
            n = self.profile['n_explicit_damp']
            self.gbflux = np.mean(self.gbflux_i[:,:,:,n:-n,:],axis=3)

    #---------------------------------------------------------------------------#

    def make_diff(self):
        """Makes diff.  Output is dictionary of numpy arrays with
        dimensions: n_kinetic x n_field x n_time"""

        # NOTE: deprecate n_x_offset in GYRO.

        if self.gbflux == []:
            self.make_gbflux()

        ir_norm = int(self.profile['n_x']/2+1)

        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_time    = self.t['n_time']

        self.diff = np.zeros((n_kinetic,n_field,2,n_time))
        for i in range(n_kinetic):
            # Density diffusivity
            self.diff[i,:,0,:] = self.gbflux[i,:,0,:]/(
                self.profile['dlnndr_s'][i,ir_norm]*
                self.profile['den_s'][i,ir_norm])
            # Energy diffusivity
            self.diff[i,:,1,:] = self.gbflux[i,:,1,:]/(
                self.profile['dlntdr_s'][i,ir_norm]*
                self.profile['tem_s'][i,ir_norm]*
                self.profile['den_s'][i,ir_norm])


    #---------------------------------------------------------------------------#

    def make_diff_i(self):
        """Makes diff_i.  Output is dictionary of numpy arrays with
        dimensions: n_kinetic x n_field x n_x x n_time"""

        if self.gbflux_i == []:
            self.read_gbflux_i()
 
        ir_norm = int(self.profile['n_x']/2+1)

        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_x       = self.profile['n_x']
        n_time    = self.t['n_time']

        self.diff_i = np.zeros((n_kinetic,n_field,2,n_x,n_time))
        for i in range(n_kinetic):
            # Density diffusivity
            self.diff_i[i,:,0,:,:] = self.gbflux_i[i,:,0,:,:]/(
                self.profile['dlnndr_s'][i,ir_norm]*
                self.profile['den_s'][i,ir_norm])
            # Energy diffusivity
            self.diff_i[i,:,1,:,:] = self.gbflux_i[i,:,1,:,:]/(
                self.profile['dlntdr_s'][i,ir_norm]*
                self.profile['tem_s'][i,ir_norm]*
                self.profile['den_s'][i,ir_norm])


    #---------------------------------------------------------------------------#

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

    def get_moment_t(self):
        '''
        From moment_n and moment_e get moment_t.
        
        Details:
        E = (3/2) n T
        \delta E = (3/2) [ (\delta n) T + n \delta T ]  =>
        \delta T = [ (2/3) (\delta E) - (\delta n) T ] / n
        '''
        if hasattr(self,'moment_t'):
            return self.moment_t
        nt = self.n
        n_theta_plot = self.profile['n_theta_plot']
        n_x          = self.profile['n_x']
        n_kinetic    = self.profile['n_kinetic']
        n_n          = self.profile['n_n']
        n = self.profile['den_s']
        T = self.profile['tem_s']
        #if 'moment_n' not in self.loaded:
          #self.load_nc_data('moment_n')
        #if 'moment_e' not in self.loaded:
          #self.load_nc_data('moment_e')
        #moment_n = self.moment_n
        #moment_e = self.moment_e
        moment_n = self.get_nc_data('moment_n')
        moment_e = self.get_nc_data('moment_e')
        moment_t = np.empty(moment_n.shape,order='F',dtype=moment_n.dtype)
        for i_r in range(n_x):
          for i_s in range(n_kinetic):
            moment_t[:,i_r,i_s,:,:] = ((2./3)*moment_e[:,i_r,i_s,:,:] - 
                                              T[i_s,i_r]*moment_n[:,i_r,i_s,:,:])/n[i_s,i_r]
        #self.loaded.append('moment_t')
        self.moment_t = moment_t
        return moment_t

    def get_sim_exp_flux(self,field_num=None,ft1=.5,ft2=1.,moment=1,
          GB_norm=0,n_time_bins=10):
        """
        Return the time averaged simulated flux and experimental flux for the
        radius at the box center.  Results, keyed off of the keywords, are
        cached for faster subsequent retrieval.
        
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
        n_time_bins: The average and uncertainty is obtained by splitting up
                     the region ft1*max(time) < time <ft2*max(time) into
                     n_time_bins bins, then taking the average in each of those
                     bins, and taking the standard deviation of the mean 
                     (sigma/sqrt(n_time_bins))as the uncertainty. 
                     (10 is default.)
        """
        import math
        key = 'field=%s_ft1=%g_ft2=%g_moment=%d_GBnorm=%d_n_time_bins=%d'%(
                field_num,ft1,ft2,moment,GB_norm,n_time_bins)
        if key in self.xp_sim_flux:
            return self.xp_sim_flux[key]
        #if self.gbflux_i == []:
          #read_gbflux_i()#self.make_gbflux()
        ir_norm = int(self.profile['n_x']/2+1)
        gbflux = self.get_nc_data('gbflux_i')[:,:,:,ir_norm,:]#self.gbflux_i[:,:,:,ir_norm,:]
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
        means = []
        time_inds = map(int,np.linspace(ind_t1,ind_t2,n_time_bins+1))
        for ti,i1 in enumerate(time_inds[:-1]):
          means.append(np.mean(gbflux[:,i1:time_inds[ti+1]],axis=1)) #Mean
        means = np.array(means)
        #print means.shape
        gbflux_m = np.mean(means,axis=0)
        #print gbflux_m.shape
        #print repr(means)
        gbflux_u = np.std(means,axis=0)/n_time_bins**0.5  #Standard Deviation
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
        self.xp_sim_flux[key] = (xp,sim)
        return (xp,sim)
    
    def get_ktheta_spectrum(self,field_num=None,moment=1,GB_norm=0,):
        '''
        Return the ktheta spectra as a function of time of the given flux
        moment.   The result is a dictionary whose keys correspond to the
        various species.
        
        --------
        Keywords
        field_num:  0 - Electrostatic field part
                    1 - Flutter
                    2 - Compression
                    None - Sum (Default)
        moment:     0 - Particle flux
                    1 - Thermal flux (Default)
                    2 - Momentum flux
                    3 - Exchange flux
        GB_norm:    0 - No normalization [unit/m^2] (Default)
                    1 - Normalized to gyrobohm flux []
                    2 - Flow [unit] # Not implemented, need input.profiles.extra
        
        '''
        #if self.gbflux_n == []:
          #self.read_gbflux_n()
        gbflux = self.get_nc_data('gbflux_n')[:,:,moment,:,:]
        if field_num==None:
          gbflux = np.sum(gbflux,axis=1)
        else:
          gbflux = gbflux[:,field_num,:,:]
        fac_gyro = 1.
        if GB_norm == 0:
          fac_gyro = self.units[moment+9] # Convert MW/m^2/(MW/m^2) to MW/m^2
        elif GB_norm == 1:
          fac_gyro = 1.# Already in GB units
        elif GB_norm == 2:
          fac_gyro = self.units[moment+9]*surface_area #Convert MW/m^2/(MW/m^2) to MW
        result = {}
        for ti,tg in enumerate(self.tagspec):
          result[tg] = gbflux[ti,:,:]*fac_gyro
        return result
    
    def get_flux_time_trace(self,field_num=None,moment=1,GB_norm=0):
        '''
        Return the time trace of the given flux moment for the radius at the 
        center of the box.  The result is a dictionary whose keys correspond to 
        the various species.
        
        --------
        Keywords
        
        field_num:  0 - Electrostatic field part
                    1 - Flutter
                    2 - Compression
                    None - Sum (Default)
        moment:     0 - Particle flux
                    1 - Thermal flux (Default)
                    2 - Momentum flux
                    3 - Exchange flux
        GB_norm:    0 - No normalization [unit/m^2] (Default)
                    1 - Normalized to gyrobohm flux []
                    2 - Flow [unit] # Not implemented, need input.profiles.extra
        '''
        #if self.gbflux_i == []:
          #self.read_gbflux_i()
        ir_norm = int(self.profile['n_x']/2+1)
        #gbflux = self.gbflux_i[:,:,:,ir_norm,:]
        try:
          gbflux = self.get_nc_data('gbflux_i')[:,:,:,ir_norm,:]
        except:
          if 'gbflux_i' not in self.loaded:
            self.read_gbflux_i()
          gbflux = self.gbflux_i[:,:,:,ir_norm,:]
        if field_num==None:
          gbflux = np.sum(gbflux,axis=1)
        else:
          gbflux = gbflux[:,field_num,:,:]
        gbflux = gbflux[:,moment,:]
        
        fac_gyro = 1.
        if GB_norm == 0:
          fac_gyro = self.units[moment+9] # Convert MW/m^2/(MW/m^2) to MW/m^2
        elif GB_norm == 1:
          fac_gyro = 1.# Already in GB units
        elif GB_norm == 2:
          fac_gyro = self.units[moment+9]*surface_area #Convert MW/m^2/(MW/m^2) to MW
        result = {}
        for ti,tg in enumerate(self.tagspec):
          result[tg] = gbflux[ti,:]*fac_gyro
        return result
    
    def get_flux_traces(self,field_num=None,moment=1,GB_norm=0):
        '''
        Return the radial and time trace of the given flux moment.  The result
        is a dictionary whose keys correspond to the various species.
        
        --------
        Keywords
        
        field_num:  0 - Electrostatic field part
                    1 - Flutter
                    2 - Compression
                    None - Sum (Default)
        moment:     0 - Particle flux
                    1 - Thermal flux (Default)
                    2 - Momentum flux
                    3 - Exchange flux
        GB_norm:    0 - No normalization [unit/m^2] (Default)
                    1 - Normalized to gyrobohm flux []
                    2 - Flow [unit] # Not implemented, need input.profiles.extra
        '''
        #if self.gbflux_i == []:
        #  self.read_gbflux_i()
        #gbflux = self.gbflux_i[:,:,moment,:,:]
        gbflux = self.get_nc_data('gbflux_i')[:,:,moment,:,:]
        if field_num==None:
          gbflux = np.sum(gbflux,axis=1)
        else:
          gbflux = gbflux[:,field_num,:,:]
        
        fac_gyro = 1.
        if GB_norm == 0:
          fac_gyro = self.units[moment+9] # Convert MW/m^2/(MW/m^2) to MW/m^2
        elif GB_norm == 1:
          fac_gyro = 1.# Already in GB units
        elif GB_norm == 2:
          fac_gyro = self.units[moment+9]*surface_area #Convert MW/m^2/(MW/m^2) to MW
        result = {}
        for ti,tg in enumerate(self.tagspec):
          result[tg] = gbflux[ti,:,:]*fac_gyro
        return result

    def get_rms_midplane_flucs(self,quant,ft1=.5,ft2=1,spec=-1,n1=1,n2=None,
            localnorm=True,n_time_bins=10):
        '''
        Return the rms midplane fluctuations of quant for species indexed by
        spec and field indexed by fields.
        
        Arguments:
        
        quant: 
          'temp' - Return the temperature fluctuations
          'dens' - Return the density fluctuations
          'energy' - Return the energy fluctuations
        
        Keywords:
        
        ft1: The fraction of total time to use as the starting point for the 
              time averaging [default:0.5]
        ft2: The fraction of total time to use as the final point for the time 
              averaging [default:1]
        spec: Index of which species to use (int) 0,1...(n_species-1) [default:-1]
        n1: The first mode number index [default:1]
        n2: The last mode number index [default:None, means last]
        localnorm: Boolean indicating whether to use a local equilibrium 
          quantity or the box center value for normalization [default:True]
        n_time_bins: The number of time bins to use for averaging and getting the
          standard deviation of the mean
        
        Patterned after the postGYRO plot_gyro_midplane_rmsflucamps.pro
        '''
        n_theta_plot = self.profile['n_theta_plot']
        n_x = self.profile['n_x']
        n_n = self.profile['n_n']
        nt = self.n
        tind_i = int(ft1*nt)
        tind_f = int(ft2*nt)+1
        if n2==None:
          n2 = n_n
        T = self.profile['tem_s']
        dens = self.profile['den_s']
        if quant == 'temp':
          #if 'moment_t' not in self.loaded:
          #  self.make_moment_t()
          #field = self.moment_t
          field = self.get_moment_t()
          norm = T
        elif quant == 'dens':
          #if 'moment_n' not in self.loaded:
            #self.load_nc_data('moment_n')
          #field = self.moment_n
          field = self.get_nc_data('moment_n')
          norm = dens
        elif quant == 'energy':
          #if 'moment_e' not in self.loaded:
          #  self.read_moment_e()
          #field = self.moment_e
          field = self.get_nc_data('moment_e')
          norm = 1.5*T*dens
        else:
          raise NotImplemented('quant = %s is not recognized'%quant)
        field = field[n_theta_plot/2,:,spec,:,tind_i:tind_f+1]
        if localnorm:
          norm = norm[spec,:]
        else:
          norm = norm[spec,n_x/2]
        if n1!=0:
          pwr = 2*np.sum(abs(field[:,n1:n2,:])**2,axis=1)
        else:
          pwr = 2*np.sum(abs(field[:,n1+1:n2,:])**2,axis=1)+abs(field[:,0,:])**2
        amps = []
        ntimes = tind_f-tind_i
        inds = map(int,np.linspace(0,ntimes,n_time_bins+1))
        for ti,i1 in enumerate(inds[:-1]):
          amps.append(np.mean(pwr[:,i1:inds[ti+1]],axis=1)**0.5/norm)
        amps = np.array(amps)
        amps_mean = np.mean(amps,axis=0)
        amps_err = np.std(amps,axis=0)/n_time_bins
        amp = (np.sum(pwr,axis=1)/(tind_f-tind_i+1))**0.5
        amp = amp / norm
        r = self.profile['r']
        n_bnd = self.profile['n_explicit_damp']
        fluc_mean = np.mean(amp[n_bnd:-n_bnd+1])
        fluc_err = np.std(amp[n_bnd:-n_bnd+1])
        return {'amp':amp,'r':r,'fluc_mean':fluc_mean,'fluc_err':fluc_err,
                'amps_mean':amps_mean,'amps_err':amps_err}
    
    def get_midplane_flucs(self,quant,n1=1,n2=None,localnorm=True):
        '''
        Patterned after the postGYRO plot_gyro_midplane_rmsflucamps.pro
        
        Return the midplane fluctuations of quant as a function of r and t.
        
        Arguments:
        
        quant: 
          'temp' - Return the temperature fluctuations
          'dens' - Return the density fluctuations
          'energy' - Return the energy fluctuations
        
        Keywords:
        
        n1: The first mode number index [default:1]
        n2: The last mode number index [default:None, means last]
        localnorm: Boolean indicating whether to use a local equilibrium 
          quantity or the box center value for normalization [default:True]
        '''
        n_theta_plot = self.profile['n_theta_plot']
        n_x = self.profile['n_x']
        n_n = self.profile['n_n']
        nt = self.n
        if n2==None:
          n2 = n_n
        T = self.profile['tem_s']
        dens = self.profile['den_s']
        if quant == 'temp':
          #if 'moment_t' not in self.loaded:
          #  self.make_moment_t()
          #field = self.moment_t
          field = self.get_moment_t()
          norm = T
        elif quant == 'dens':
          #if 'moment_n' not in self.loaded:
          #  self.load_nc_data('moment_n')
          #field = self.moment_n
          field = self.get_nc_data('moment_n')
          norm = dens
        elif quant == 'energy':
          #if 'moment_e' not in self.loaded:
          #  self.load_nc_data('moment_e')
          #field = self.moment_e
          field = self.get_nc_data('moment_e')
          norm = 1.5*T*dens
        else:
          raise NotImplemented('quant = %s is not recognized'%quant)
        field = field[n_theta_plot/2,:,:,:,:]
        if localnorm:
          norm = norm[:,:]
        else:
          norm = norm[:,n_x/2:n_x/2+1]
        if n1!=0:
          pwr = 2*np.sum(abs(field[:,:,n1:n2,:])**2,axis=2)
        else:
          pwr = 2*np.sum(abs(field[:,:,1:n2,:])**2,axis=2)+abs(field[:,:,0,:])**2
        pwr = pwr**0.5/norm.T[:,:,None]
        r = self.profile['r']
        n_bnd = self.profile['n_explicit_damp']
        result = {}
        for ki,k in enumerate(self.tagspec):
          result[k] = pwr[:,ki,:]
        return result

    def convert_to_netcdf(self,del_moment=True,verify=True):
        '''
        An experiment to see if converting to netCDF compacts the data, and also
        whether doing so makes it faster to read in
        
        Returns an array of accumulating times (in seconds) for accomplishing
        certain parts of the process.
        
        If verify, then verify_nc is also called at the end of this function.
        
        The keyword del_moment is used to indicate whether the moment attributes
        of this instance should be deleted to save memory
        '''
        import netCDF4
        import getpass
        import os
        import time
        t0 = time.time()
        fn = self.dirname+'/out.gyro.nc'
        if os.path.exists(fn):
            try:
                self.verify_nc(verbose=False)
            except (AssertionError,KeyError,ValueError) as e:
                print 'netcdf data is not consistent with current data; '+ \
                      'removing out.gyro.nc'
                os.system('rm %s'%fn)
            else:
                return np.array([time.time()])-t0
        print 'Creating netcdf file'
        nc = netCDF4.Dataset(fn,'w')
        nc.history = 'Created from files in %s on %s'%(self.dirname,
          os.environ['HOST'])
        nc.user = getpass.getuser()
        nc.createDimension('t',self.n)
        t = nc.createVariable('t',self.t['(c_s/a)t'].dtype,('t',))
        t[:] = self.t['(c_s/a)t']
        t.units = '(c_s/a)'
        nc.createDimension('r',self.profile['n_x'])
        r = nc.createVariable('r',self.profile['r'].dtype,('r',))
        r[:] = self.profile['r']
        r.units = 'r/a'
        nc.createDimension('units',len(self.units))
        units = nc.createVariable('units','f',('units',))
        units[:] = self.units
        nc.createDimension('species',self.profile['n_spec'])
        nc.species = ','.join(self.tagspec)
        nc.createDimension('ky',self.profile['n_n'])
        ky = nc.createVariable('ky',self.profile['kt_rho'].dtype,('ky',))
        ky[:] = self.profile['kt_rho']
        nc.createDimension('theta_plot',self.profile['n_theta_plot'])
        theta_plot = nc.createVariable('theta_plot','f',('theta_plot',))
        theta_plot[:] = self.profile['theta_plot']
        nc.createDimension('fields',self.profile['n_field'])
        nc.fields = ','.join(self.tagfieldtext[0:self.profile['n_field']])
        nc.createDimension('moments',4)
        nc.moments = ','.join(self.tagmomtext)
        times = [time.time()]
        if 'gbflux_i' not in self.loaded:
            self.read_gbflux_i()
        times.append(time.time())
        gbflux_i = nc.createVariable('gbflux_i',self.gbflux_i.dtype,
                      ('species','fields','moments', 'r','t'))
        #print gbflux_i.shape,self.gbflux_i.shape
        gbflux_i[:] = self.gbflux_i[...,:self.n]
        times.append(time.time())
        if 'gbflux_n' not in self.loaded:
            self.read_gbflux_n()
        gbflux_n = nc.createVariable('gbflux_n',self.gbflux_n.dtype,
                      ('species','fields','moments', 'ky','t'))
        #print gbflux_n.shape,self.gbflux_n.shape
        gbflux_n[:] = self.gbflux_n[...,:self.n]
        try:
            if 'gbflux_exc' not in self.loaded:
                self.read_gbflux_exc()
            gbflux_exc = nc.createVariable('gbflux_exc',self.gbflux_exc.dtype,
                          ('species','moments','t'))
            gbflux_exc[:] = self.gbflux_exc[...,:self.n]
        except IOError:
            pass
        for moment in ['u','n','v','e']:
            try:
                print 'Loading moment_'+moment
                times.append(time.time())
                if 'moment_'+moment not in self.loaded:
                    getattr(self,'read_moment_'+moment)()#.read_moment_u()
                times.append(time.time())
            except IOError:
                print 'Not loading moment_'+moment
            else:
                dim = 'species'
                if moment=='u':
                    dim = 'fields'
                moment_i = nc.createVariable('moment_%s_i'%moment,
                      self.gbflux_i.dtype, ('theta_plot','r',dim,'ky','t'))
                moment_i[:] = getattr(self,'moment_'+moment)[...,:self.n].imag
                moment_r = nc.createVariable('moment_%s_r'%moment,
                      self.gbflux_i.dtype, ('theta_plot','r',dim,'ky','t'))
                moment_r[:] = getattr(self,'moment_'+moment)[...,:self.n].real
                if del_moment and not verify:
                    delattr(self,'moment_'+moment)
                    self.loaded.pop(self.loaded.index('moment_'+moment))
        nc.close()
        times.append(time.time())
        print 'Done creating netcdf file.'
        if verify:
            self.verify_nc(del_moment=del_moment)
        times.append(time.time())
        return np.array(times)-t0
    
    def verify_nc(self,del_moment=True,verbose=True):
        '''
        Compare the netcdf version to the text version in memory
        Keywords:
        del_moment - If del_moment is True then delete the moment_{u,n,e,v} attributes after
            verifying them
        '''
        import netCDF4
        import time
        from numpy import nanmax
        if verbose:
          print 'Verifying out.gyro.nc against the version in memory'
        t0 = time.time()
        nc = netCDF4.Dataset(self.dirname+'/out.gyro.nc','r')
        vars = nc.variables
        for v in ['gbflux_i','gbflux_n','gbflux_exc']:
            if os.path.exists(self.dirname+'/out.gyro.'+v) and not v in vars:
                raise KeyError(
                  'out.gyro.%s exists, but it is not in out.gyro.nc'%v)
        for v in ['u','v','n','e']:
            mom = 'moment_'+v
            if os.path.exists(self.dirname+'/out.gyro.'+mom) and mom+'_i' not in vars:
                raise KeyError(
                  'out.gyro.%s exists, but it is not in out.gyro.nc'%mom)
        #print repr(nc)
        assert nanmax(abs(self.t['(c_s/a)t']-nc.variables['t'][:])) == 0, \
          'time not written/read correctly'
        assert nanmax(abs(self.profile['r']-nc.variables['r'][:])) == 0, \
          'rmin/a not written/read correctly'
        assert nanmax(abs(self.profile['kt_rho']-nc.variables['ky'][:])) == 0, \
          'kt_rho not written/read correctly'
        for gb in ['gbflux_i','gbflux_n','gbflux_exc']:
            if gb in self.loaded:
                gb_var = getattr(self,gb)
                assert nanmax(abs(getattr(self,gb)[...,:self.n]-nc.variables[gb][:])) == 0, \
                  '%s not written/read correctly'%gb
            elif verbose:
                print 'Skipping %s because it is not loaded in memory'%(gb)
        for moment in ['u','n','e','v']:
            if 'moment_'+moment in self.loaded:
                mom = getattr(self,'moment_'+moment)[...,:self.n]
                nc_moment = nc.variables['moment_%s_i'%moment][:]*1j+\
                            nc.variables['moment_%s_r'%moment][:]
                assert nanmax(abs(mom-nc_moment)) == 0, \
                  'moment_%s not written/read correctly'%moment
                if del_moment:
                  delattr(self,'moment_'+moment)
                  self.loaded.remove('moment_'+moment)
            elif verbose:
                print 'Skipping moment_%s because it is not loaded in memory'%(
                  moment)
        nc.close()
        if verbose:
          print 'Done verifying'
    
    def get_nc_data(self,var):
        '''
        Return the data in variable <var> from out.gyro.nc
        '''
        import netCDF4
        nc = netCDF4.Dataset(self.dirname+'/out.gyro.nc','r')
        if 'moment' not in var:
            data = nc.variables[var][:]
        else:
            data = nc.variables[var+'_i'][:]*1j+nc.variables[var+'_r'][:]
        nc.close()
        return data
    
    def load_nc_data(self,var):
        '''
        Load the data in variable <var> from out.gyro.nc into this instance
        '''
        try:
            data = self.get_nc_data(var)
            setattr(self,var,data)
            self.loaded.append(var)
        except:
            print 'Error loading from netcdf data'
            print 'Reverting to text file'
            getattr(self,'read_'+var)()
    
    def get_rho(self):
        '''
        Get the rho mapping for the minor radius variable
        '''
        r = self.profile['r']
        rmin = self.input_profile.data['rmin']
        rmin = rmin/rmin.max()
        rho = self.input_profile.data['rho']
        return np.interp(r,rmin,rho)
    
    def get_time_average(self,quant,time_axis,ft1=0.2,ft2=1,n_time_bins=10):
        '''
        Given the array quant, which has some time axis, determine the mean and
        standard deviation of the mean for n_time_bins separate samples along 
        the time axis. 
        '''
        nt = self.n
        tind_i = int(ft1*nt)
        tind_f = int(ft2*nt)
        inds = map(int,np.linspace(tind_i,tind_f,n_time_bins+1))
        means = []
        for ti,i1 in enumerate(inds[:-1]):
            means.append(np.mean(np.take(quant,range(i1,inds[ti+1]),axis=time_axis),
              axis=time_axis))
        means = np.array(means)
        #print means.shape
        mean = np.mean(means,axis=0)
        std = np.std(means,axis=0)/n_time_bins**.5
        return {'mean':mean,'std':std}
    
    def get_rho_deriv(self,quant):
        '''
        Return the rho derivative at the box center for <quant>, where <quant>
        is a valid key in the input_profile.data dict.
        '''
        from numpy import argmin
        from deriv import deriv as dydx
        q = self.input_profile.data[quant]
        rho = self.input_profile.data['rho']
        rho_ind = argmin(abs(rho-self.get_rho()[self.profile['n_x']/2]))
        return dydx(rho*self.input_profile.max_rho,q)[rho_ind]
    
    def deal_with_large_exponents(self,fn):
        '''
        When GYRO writes a number that has a three digit exponent, it dispenses with
        the E (2.5l234+234).  This function looks at the given file <fn> and rewrites
        numbers with a plus, but no E to have an E.  IDL could deal with the E-less
        form, but python cannot.
        '''
        import re
        patt=re.compile('(\d+)\+(\d+)',flags=re.MULTILINE)
        f = open(self.dirname+'/'+fn,'r')
        fr = f.read()
        f.close()
        if patt.search(fr):
            print 'Modifying ',self.dirname+'/'+fn
            f = open(self.dirname+'/'+fn,'w')
            f.write(patt.sub(r'\1E+\2',fr))
            f.close()
