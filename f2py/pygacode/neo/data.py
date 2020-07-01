import numpy as np

class NEOData:
    """A class for management of NEO output data.

    Data:

    Example Usage:
        >>> from gacode.python.neo.data import NEOData
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
        self.read_theory()
        self.read_transport()
        self.read_transport_gv()
        self.read_equil()
        self.read_vel()
        self.read_rotation()
        self.read_neoinputgeo()
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
        self.phi           = {}
        self.vel           = []
        self.veltag        = []
        self.rotation      = {}
        self.geo           = {}
        self.flow          = {}

        
    #-------------------------------------------------------------------------#
    
    def join(self,obj):
        from copy import deepcopy

        copy = deepcopy(obj)
        
        copy.grid['n_radial'] += 1
        copy.grid['r_over_a'] =  np.hstack((self.grid['r_over_a'],obj.grid['r_over_a']))
        copy.vel = np.vstack(( self.vel , obj.vel ))
        
        dicts = 'equil','theory','theory_nclass','transport','transport_gv','transport_exp','rotation','flow','geo'
        
        for d in dicts:
            keys = getattr(copy,d)
            for k in keys:
               getattr(copy,d)[k] = np.vstack((getattr(self,d)[k],getattr(obj,d)[k]))
        
        return copy
        
        
    #-------------------------------------------------------------------------#

    def set_directory(self, path):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        self.dirname = expanduser(expandvars(path))

    #-------------------------------------------------------------------------#
    
    def read_neoinputgeo(self):
        """Load shape of magnetics flux surfaces"""                
        
        def mom2rz(rcos,rsin,zcos,zsin,theta):
            #composition of the flux surfaces from the Fourier moments
            #rcos - jt,jrho,jmom  (or just jrho,jmom, or jmom)

            nmom = np.size(rcos,0)
            angle = np.outer(np.arange(nmom),theta )
            cos = np.cos(angle)
            sin = np.sin(angle)
            cos[0]/= 2
            
            r_plot = np.tensordot(rcos,cos,axes=([0,0]))
            r_plot+= np.tensordot(rsin,sin,axes=([0,0]))
            z_plot = np.tensordot(zcos,cos,axes=([0,0]))
            z_plot+= np.tensordot(zsin,sin,axes=([0,0]))
            return r_plot,z_plot


        

        try:
            data = np.loadtxt(self.dirname+'/input.geo')
            n_fourier = int(data[0])
            fourier = data[1:] #rcos,rsin,zcos,zsin are dimensionalles!!!
            
            rcos,rsin,zcos,zsin,drcos,drsin,dzcos,dzsin=fourier.reshape((8,n_fourier+1), order='F')
 
        except:
            try:
                data = np.loadtxt(self.dirname+'/input.profiles.geo')
            
            
                n_fourier = int(data[0])
                fourier = data[1:] 
                
                rcos,rsin,zcos,zsin = fourier.reshape((4,n_fourier+1 ,-1), order='F')
  
            except:
                print("ERROR (NEOData): Fatal error!  Missing input.geo.")
                return
 

        self.geo = {'rcos':rcos,'rsin':rsin,'zcos':zcos,'zsin':zsin}
        

        theta = np.linspace(-np.pi,np.pi,1000)
        r,z = mom2rz(rcos,rsin,zcos,zsin, self.grid['theta'])

        self.geo['R'] = r
        self.geo['z'] = z

    #-------------------------------------------------------------------------#

    def read_grid(self):
        """Reads out.neo.grid"""


        try:
            data = np.loadtxt(self.dirname+'/out.neo.grid')
        except:
            if self.verbose:
                print("ERROR (NEOData): Fatal error!  Missing out.neo.grid.")
            return

        self.grid['n_species'] = int(data[0])
        self.grid['n_energy']  = int(data[1])
        self.grid['n_xi']      = int(data[2])
        self.grid['n_theta']   = int(data[3])  #BUG it is wrong!!!
        self.grid['theta']     = data[4:4+int(data[3])]
        self.grid['n_radial']  = int(data[4+int(data[3])])
        self.grid['r_over_a']  = data[5+int(data[3]):]

    #-------------------------------------------------------------------------#

    def read_equil(self):
        """Read out.neo.equil."""

        try:
            equil = np.loadtxt(self.dirname+'/out.neo.equil')
        except:
            if self.verbose:
                print("ERROR (NEOData): Fatal error!  Missing out.neo.equil.")
            return

        if len(equil.shape)==1:
            equil = equil[None,:equil.shape[0]]

        self.equil['r_over_a']      = equil[:,0]
        self.equil['dphidr']        = equil[:,1]
        self.equil['q']             = equil[:,2]
        self.equil['rho_star']      = equil[:,3]
        self.equil['R0_over_a']     = equil[:,4]
        self.equil['omega0']        = equil[:,5]
        self.equil['domega0dr']     = equil[:,6]
        self.equil['n']             = equil[:,7+0::5]
        self.equil['T']             = equil[:,7+1::5]
        self.equil['a_over_ln']     = equil[:,7+2::5]
        self.equil['a_over_lt']     = equil[:,7+3::5]
        self.equil['tauinv_norm']   = equil[:,7+4::5]

    #-------------------------------------------------------------------------#

    def read_theory(self):
        """Reads out.neo.theory."""

        try:
            data = np.loadtxt(self.dirname+'/out.neo.theory')
        except:
            if self.verbose:
                print("ERROR (NEOData): Fatal error!  Missing out.neo.theory.")
            return


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
        self.theory['HSGamma']  = data[:,15+0::3]
        self.theory['HSQ']      = data[:,15+1::3]
        self.theory['KjparB']   = data[:,15+2::3]

    #-------------------------------------------------------------------------#

    def read_transport(self):
        """Reads out.neo.transport."""

        try:
            data = np.loadtxt(self.dirname+'/out.neo.transport')
        except:
            if self.verbose:
                print("ERROR (NEOData): Fatal error!  Missing out.neo.transport.")
            return

        if len(data.shape)==1:
            data = data[None,:data.shape[0]]

        self.transport['r_over_a'] = data[:,0]
        self.transport['phisq']    = data[:,1]
        self.transport['jparB']    = data[:,2]
        self.transport['vtheta0']  = data[:,3]
        self.transport['uparB0']   = data[:,4]
        self.transport['Gamma']    = data[:,5::8]
        self.transport['Q']        = data[:,6::8]
        self.transport['Pi']       = data[:,7::8]
        self.transport['uparB']    = data[:,8::8]
        self.transport['k']        = data[:,9::8]
        self.transport['K']        = data[:,10::8]
        self.transport['vtheta']   = data[:,11::8]
        self.transport['vphi']     = data[:,12::8]


    #-------------------------------------------------------------------------#

    def read_transport_gv(self):
        """Reads out.neo.transport_gv."""

        try:
            data = np.loadtxt(self.dirname+'/out.neo.transport_gv')
        except:
            if self.verbose:
                print("ERROR (NEOData): Fatal error!  Missing out.neo.transport_gv.")
            return

        if len(data.shape)==1:
            data = data[None,:data.shape[0]]
        
        self.transport_gv['r_over_a'] = data[:,0]
        self.transport_gv['Gamma_gv'] = data[:,1+0::3]
        self.transport_gv['Q_gv']     = data[:,1+1::3]
        self.transport_gv['Pi_gv']    = data[:,1+2::3]

    #-------------------------------------------------------------------------#

    def read_transport_exp(self):
        """Reads out.neo.transport_exp."""


        try:
            data = np.atleast_2d(np.loadtxt(self.dirname+'/out.neo.transport_exp'))
        except:
            if self.verbose:
                print("ERROR (NEOData): Fatal error!  Missing out.neo.transport_exp.")
            return

        self.transport_exp['r']       = data[:,0]
        self.transport_exp['phisq']   = data[:,1]
        self.transport_exp['jparB']   = data[:,2]
        self.transport_exp['vtheta0'] = data[:,3]
        self.transport_exp['uparB0']  = data[:,4]
        self.transport_exp['Gamma']   = data[:,5::8]
        self.transport_exp['Q']       = data[:,6::8]
        self.transport_exp['Pi']      = data[:,7::8]
        self.transport_exp['uparB']   = data[:,8::8]
        self.transport_exp['k']       = data[:,9::8]
        self.transport_exp['K']       = data[:,10::8]
        self.transport_exp['vtheta']  = data[:,11::8]
        self.transport_exp['vphi']    = data[:,12::8]

    #-------------------------------------------------------------------------#

    def read_rotation(self):
        """Reads out.neo.rotation."""

        n_spec   = self.grid['n_species']
        n_theta  = self.grid['n_theta']
        n_radial = self.grid['n_radial']
 
        try:
            data = np.loadtxt(self.dirname+'/out.neo.rotation',ndmin=2)
            self.rotation['r_over_a'] = data[:,0]
            self.rotation['dphi_ave'] = data[:,1]
            self.rotation['n_ratio'] = data[:,2:n_spec*2+2:2]
            self.rotation['V_conv']  = data[:,3:n_spec*2+3:2]
            
            N = n_spec*2+2
            self.rotation['phi_theta'] = data[:,N:N+n_theta]
            N += n_theta
            self.rotation['n_ov_n0'] = data[:,N:].reshape(n_radial,n_spec,n_theta)

            
            
        except:
            print("Warning (NEOData): Missing out.neo.rotation.")
            self.rotation['r_over_a'] = self.grid['r_over_a']
            self.rotation['dphi_ave'] = np.zeros([n_radial])
            self.rotation['n_ratio']  = np.ones([n_radial,n_spec])
            self.rotation['dphi']     = np.zeros([n_radial,n_theta])
            self.rotation['V_conv']   = np.zeros([n_radial,n_spec])
            self.rotation['n_ov_n0']  = np.ones([n_radial,n_spec,n_theta])
            self.rotation['phi_theta']= np.zeros([n_radial,n_theta])


            
   
    #-------------------------------------------------------------------------#

    def read_vel(self):
        """Reads out.neo.vel."""

        try:
            data = np.loadtxt(self.dirname+'/out.neo.vel')
        except:
            if self.verbose:
                print("ERROR (NEOData): Missing out.neo.vel.")
            return

        self.vel = data.reshape((self.grid['n_radial'],
                                 self.grid['n_species'],
                                 self.grid['n_theta']))

        self.veltag = ['n_radial','n_species','n_theta']
        
        
        try:
            data = np.loadtxt(self.dirname+'/out.neo.vel_fourier',ndmin=2)
        except:
            print("ERROR (NEOData): Missing out.neo.vel_fourier.")
            raise Exception('no data') 



        n_theta = self.grid['n_theta']
        nbspecies = self.grid['n_species']
        n_radial = self.grid['n_radial']
        theta = self.grid['theta']
        
        
        m_theta = (n_theta-1)//2-1

        data = data.reshape(n_radial, nbspecies, 6, m_theta+1)

        cos_par_flow = data[:,:,0]
        sin_par_flow = data[:,:,1]

        cos_pol_flow = data[:,:,2]
        sin_pol_flow = data[:,:,3]

        cos_tor_flow = data[:,:,4]
        sin_tor_flow = data[:,:,5]


        s = np.sin(theta*np.arange(m_theta+1)[:,None])
        c = np.cos(theta*np.arange(m_theta+1)[:,None])
        
        
        self.flow['par_flow'] = np.einsum('ijk,kl',cos_par_flow, c)+np.einsum('ijk,kl',sin_par_flow, s)
        self.flow['pol_flow'] = np.einsum('ijk,kl',cos_pol_flow, c)+np.einsum('ijk,kl',sin_pol_flow, s)
        self.flow['tor_flow'] = np.einsum('ijk,kl',cos_tor_flow, c)+np.einsum('ijk,kl',sin_tor_flow, s)

     
     
     
     
