import os
import numpy as np
import sys

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
        """

        self.init_data()
        self.set_directory(sim_directory)
        self.read_profile()
        self.read_geometry()
        self.read_t()
        self.read_units()
        #self.read_input_profiles()
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

        self.make_gbflux()
        self.make_diff()
        self.make_diff_i()
        """

        import sys
        from os.path import expanduser, expandvars
        import sys
        path = '$GACODE_ROOT/python'
        sys.path.append(expanduser(expandvars(path)))

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

    #---------------------------------------------------------------------------#

    def set_directory(self, path):
        """Set the simulation directory."""

        import os

        from os.path import expanduser, expandvars
        self.dirname = os.path.realpath(expanduser(expandvars(path)))

    #---------------------------------------------------------------------------#

    def read_t(self):
        """Read out.gyro.t to get time data."""

        import sys
        import numpy as np

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
            import profiles_gen.data
            self.input_profile = profiles_gen.data.profiles_genData(fn)
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

        try:
            gbflux_i = np.fromfile(self.dirname+'/out.gyro.gbflux_i',dtype=float,sep=" ")
        except IOError:
            raise IOError("ERROR (GYROData): out.gyro.gbflux_i not found.")

        nt = len(gbflux_i)/(n_kinetic*n_field*4*n_x)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dirname+
                          '/out.gyro.gbflux_i too small. ')

        self.gbflux_i = gbflux_i.reshape((n_kinetic,n_field,4,n_x,nt),order='F')

        self.loaded.append('gbflux_i')

    #---------------------------------------------------------------------------#

    def read_gbflux_n(self):
        """Reads gbflux_n data.  Output is numpy array with dimensions:
        n_kinetic x n_field x 4 x n_n x n_time"""

        n_kinetic = self.profile['n_kinetic']
        n_field   = self.profile['n_field']
        n_n       = self.profile['n_n']


        try:
            gbflux_n = np.fromfile(self.dirname+'/out.gyro.gbflux_n',dtype=float,sep=" ")
        except IOError:
            raise IOError("ERROR (GYROData): out.gyro.gbflux_n not found.")

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

        try:
            gbflux_exc = np.fromfile(self.dirname+'/out.gyro.gbflux_exc',dtype=float,sep=" ")
        except:
            raise IOError("ERROR (GYROData): out.gyro.gbflux_exc not found.")

        nt = len(gbflux_exc)/(n_kinetic*2)

        if self.n > nt:
            raise IOError('ERROR (GYROData): '+self.dirname+
              '/out.gyro.gbflux_exc too small. ')


        self.gbflux_exc = gbflux_exc.reshape((n_kinetic,2,nt),order='F')
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

        m     = int(self.profile['box_multiplier'])
        n_x   = int(self.profile['n_x'])
        n_ang = int(self.profile['n_theta_plot']*n_x/m)

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
        if self.gbflux == []:
          self.make_gbflux()
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

