class NEOData:
    """A class for management of NEO output data.

    Data:

    Example Usage:
        >>> from pyrats.neo.data import NEOData
        >>> sim1 = NEOData('example_directory')
    """

    #-------------------------------------------------------------------------#
    # Methods

    def __init__(self, sim_directory, verbose=True):
        """Constructor reads in data from sim_directory and creates new object.
        """
        self.verbose = verbose
        self.init_data()
        self.set_directory(sim_directory)
        self.read_grid()
        self.read_equil()
        self.read_theory()
        self.read_transport()
        self.read_transport_gv()
        self.read_transport_exp()
        self.read_vel()
        self.read_rotation()

    #-------------------------------------------------------------------------#

    def init_data(self):
        """Initialize object data."""

        self.dirname       = ""
        self.grid          = {}
        self.equil         = {}
        self.theory        = {}
        self.theory_nclass = {}
        self.transport     = {}
        self.transport_gv  = {}
        self.expnorm       = {}
        self.transport_exp = {}
        self.phi           = {}
        self.vel           = []
        self.veltag        = []
        self.rotation      = {}

    #-------------------------------------------------------------------------#

    def set_directory(self, path):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        self.dirname = expanduser(expandvars(path))

    #-------------------------------------------------------------------------#

    def read_grid(self):
        """Reads out.neo.grid"""

        import sys
        import numpy as np

        try:
            data = np.loadtxt(self.dirname+'/out.neo.grid')
        except:
            if self.verbose:
                print "ERROR (NEOData): Fatal error!  Missing out.neo.grid."
            return

        self.grid['n_species'] = int(data[0])
        self.grid['n_energy']  = int(data[1])
        self.grid['n_xi']      = int(data[2])
        self.grid['n_theta']   = int(data[3])
        self.grid['theta']     = data[4:4+int(data[3])]
        self.grid['n_radial']  = int(data[4+int(data[3])])
        self.grid['r_over_a']  = data[5+int(data[3]):]

    #-------------------------------------------------------------------------#

    def read_equil(self):
        """Read out.neo.equil."""

        import sys
        import numpy as np

        try:
            equil = np.loadtxt(self.dirname+'/out.neo.equil')
        except:
            if self.verbose:
                print "ERROR (NEOData): Fatal error!  Missing out.neo.equil."
            return

        n_spec = self.grid['n_species']
        if len(equil.shape)==1:
            equil = equil[None,:equil.shape[0]]

        self.equil['r_over_a']      = equil[:,0]
        self.equil['dphidr']        = equil[:,1]
        self.equil['q']             = equil[:,2]
        self.equil['rho_star']      = equil[:,3]
        self.equil['R0_over_a']     = equil[:,4]
        self.equil['omega0']        = equil[:,5]
        self.equil['domega0dr']     = equil[:,6]
        self.equil['n']             = equil[:,7+0*n_spec:7+1*n_spec]
        self.equil['T']             = equil[:,7+1*n_spec:7+2*n_spec]
        self.equil['a_over_ln']     = equil[:,7+2*n_spec:7+3*n_spec]
        self.equil['a_over_lt']     = equil[:,7+3*n_spec:7+4*n_spec]
        self.equil['tauinv_norm']   = equil[:,7+4*n_spec:7+5*n_spec]

    #-------------------------------------------------------------------------#

    def read_theory(self):
        """Reads out.neo.theory."""

        import sys
        import numpy as np

        try:
            data = np.loadtxt(self.dirname+'/out.neo.theory')
        except:
            if self.verbose:
                print "ERROR (NEOData): Fatal error!  Missing out.neo.theory."
            return

        n_spec = self.grid['n_species']

        if len(data.shape)==1:
            data = data[None,:data.shape[0]]

        self.theory['r_over_a'] = data[:,0]
        self.theory['HHGamma']  = data[:,1]
        self.theory['HHQi']     = data[:,2]
        self.theory['HHQe']     = data[:,3]
        self.theory['HHjparB']  = data[:,4]
        self.theory['HHk']      = data[:,5]
        self.theory['HHuparB']  = data[:,6]
        self.theory['HHvtheta'] = data[:,7]
        self.theory['CHQi']     = data[:,8]
        self.theory['TGQi']     = data[:,9]
        self.theory['SjparB']   = data[:,10]
        self.theory['Sk']       = data[:,11]
        self.theory['SuparB']   = data[:,12]
        self.theory['Svtheta']  = data[:,13]
        self.theory['HRphisq']  = data[:,14]
        self.theory['HSGamma']  = data[:,15+0*n_spec:15+1*n_spec]
        self.theory['HSQ']      = data[:,15+1*n_spec:15+2*n_spec]
        self.theory['KjparB']   = data[:,15+2*n_spec]

    #-------------------------------------------------------------------------#

    def read_transport(self):
        """Reads out.neo.transport."""

        import sys
        import numpy as np

        try:
            data = np.loadtxt(self.dirname+'/out.neo.transport')
        except:
            if self.verbose:
                print "ERROR (NEOData): Fatal error!  Missing out.neo.transport."
            return

        n_spec = self.grid['n_species']
        if len(data.shape)==1:
            data = data[None,:data.shape[0]]

        self.transport['r_over_a'] = data[:,0]
        self.transport['phisq']    = data[:,1]
        self.transport['jparB']    = data[:,2]
        self.transport['vtheta0']  = data[:,3]
        self.transport['uparB0']   = data[:,4]
        self.transport['Gamma']    = data[:,5+0*n_spec:5+1*n_spec]
        self.transport['Q']        = data[:,5+1*n_spec:5+2*n_spec]
        self.transport['Pi']       = data[:,5+2*n_spec:5+3*n_spec]
        self.transport['uparB']    = data[:,5+3*n_spec:5+4*n_spec]
        self.transport['k']        = data[:,5+4*n_spec:5+5*n_spec]
        self.transport['K']        = data[:,5+5*n_spec:5+6*n_spec]
        self.transport['vtheta']   = data[:,5+6*n_spec:5+7*n_spec]
        self.transport['vphi']     = data[:,5+7*n_spec:5+8*n_spec]

    #-------------------------------------------------------------------------#

    def read_transport_gv(self):
        """Reads out.neo.transport_gv."""

        import sys
        import numpy as np

        try:
            data = np.loadtxt(self.dirname+'/out.neo.transport_gv')
        except:
            if self.verbose:
                print "ERROR (NEOData): Fatal error!  Missing out.neo.transport_gv."
            return

        n_spec = self.grid['n_species']
        if len(data.shape)==1:
            data = data[None,:data.shape[0]]

        self.transport_gv['r_over_a'] = data[:,0]
        self.transport_gv['Gamma_gv'] = data[:,1+0*n_spec:1+1*n_spec]
        self.transport_gv['Q_gv']     = data[:,1+1*n_spec:1+2*n_spec]
        self.transport_gv['Pi_gv']    = data[:,1+2*n_spec:1+3*n_spec]

    #-------------------------------------------------------------------------#

    def read_transport_exp(self):
        """Reads out.neo.transport_exp."""

        import sys
        import numpy as np

        try:
            data = np.atleast_2d(np.loadtxt(self.dirname+'/out.neo.transport_exp'))
        except:
            if self.verbose:
                print "ERROR (NEOData): Fatal error!  Missing out.neo.transport_exp."
            return

        n_spec = self.grid['n_species']

        self.transport_exp['r']       = data[:,0]
        self.transport_exp['phisq']   = data[:,1]
        self.transport_exp['jparB']   = data[:,2]
        self.transport_exp['vtheta0'] = data[:,3]
        self.transport_exp['uparB0']  = data[:,4]
        self.transport_exp['Gamma']   = data[:,5+0*n_spec:5+1*n_spec]
        self.transport_exp['Q']       = data[:,5+1*n_spec:5+2*n_spec]
        self.transport_exp['Pi']      = data[:,5+2*n_spec:5+3*n_spec]
        self.transport_exp['uparB']   = data[:,5+3*n_spec:5+4*n_spec]
        self.transport_exp['k']       = data[:,5+4*n_spec:5+5*n_spec]
        self.transport_exp['K']       = data[:,5+5*n_spec:5+6*n_spec]
        self.transport_exp['vtheta']  = data[:,5+6*n_spec:5+7*n_spec]
        self.transport_exp['vphi']    = data[:,5+7*n_spec:5+8*n_spec]

    #-------------------------------------------------------------------------#

    def read_rotation(self):
        """Reads out.neo.rotation."""

        import sys
        import numpy as np

        n_spec   = self.grid['n_species']
        n_theta  = self.grid['n_theta']
        n_radial = self.grid['n_radial']

        try:
            data = np.loadtxt(self.dirname+'/out.neo.rotation')
            self.rotation['r_over_a'] = data[:,0]
            self.rotation['dphi_ave'] = data[:,1]
            self.rotation['n_ratio']  = data[:,2+0*n_spec:2+1*n_spec]
            self.rotation['dphi']     = data[:,2+1*n_spec:2+1*n_spec+1*n_theta]
        except:
            if self.verbose:
                print "Warning (NEOData): Missing out.neo.rotation."
            self.rotation['r_over_a'] = self.grid['r_over_a']
            self.rotation['dphi_ave'] = np.zeros([n_radial])
            self.rotation['n_ratio']  = np.zeros([n_radial,n_spec])
            self.rotation['dphi']     = np.zeros([n_radial,n_theta])

    #-------------------------------------------------------------------------#

    def read_vel(self):
        """Reads out.neo.vel."""

        import sys
        import numpy as np

        try:
            data = np.loadtxt(self.dirname+'/out.neo.vel')
        except:
            if self.verbose:
                print "ERROR (NEOData): Missing out.neo.vel."
            return

        self.vel = data.reshape((self.grid['n_radial'],
                                 self.grid['n_species'],
                                 self.grid['n_theta']))

        self.veltag = ['n_radial','n_species','n_theta']
