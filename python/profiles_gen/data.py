vars_input_profiles = [
    ['rho', 'rmin', 'polflux', 'q', 'omega0'],
    ['rmaj', 'zmag', 'kappa', 'delta', 'zeta'],
    ['ne', 'Te', 'ptot', 'z_eff', 'NULL'],
    ['ni_1', 'ni_2', 'ni_3', 'ni_4', 'ni_5'],
    ['ni_6', 'ni_7', 'ni_8', 'ni_9', 'ni_10'],
    ['Ti_1', 'Ti_2', 'Ti_3', 'Ti_4', 'Ti_5'],
    ['Ti_6', 'Ti_7', 'Ti_8', 'Ti_9', 'Ti_10'],
    ['vtor_1', 'vtor_2', 'vtor_3', 'vtor_4', 'vtor_5'],
    ['vtor_6', 'vtor_7', 'vtor_8', 'vtor_9', 'vtor_10'],
    ['vpol_1', 'vpol_2', 'vpol_3', 'vpol_4', 'vpol_5'],
    ['vpol_6', 'vpol_7', 'vpol_8', 'vpol_9', 'vpol_10'],
    ['flow_beam', 'flow_wall', 'flow_mom', 'NULL', 'NULL'],
    ['pow_e', 'pow_i', 'pow_ei', 'pow_e_aux', 'pow_i_aux'],
    ['pow_e_fus', 'pow_i_fus', 'pow_e_sync', 'pow_e_brem', 'pow_e_line'],
    ['sbeame','sbcx','sscxl','NULL','NULL'],
]

vars_input_profiles_extra = [
    'bunit', 's', 'drmaj', 'dzmag', 'sdelta', 'skappa', 'szeta', 
    'dlnnedr', 'dlntedr',
    'dlnnidr_1', 'dlnnidr_2', 'dlnnidr_3', 'dlnnidr_4', 'dlnnidr_5',
    'dlnnidr_6', 'dlnnidr_7', 'dlnnidr_8', 'dlnnidr_9', 'dlnnidr_10',
    'dlntidr_1', 'dlntidr_2', 'dlntidr_3', 'dlntidr_4', 'dlntidr_5',
    'dlntidr_6', 'dlntidr_7', 'dlntidr_8', 'dlntidr_9', 'dlntidr_10',
    'dlnptotdr', 'drdrho', 'w0p', 'vol', 'volp', 'cs', 'rhos', 'ni_new', 'dlnnidr_new',
    'grad_r0', 'ave_grad_r', 'bp0', 'bt0', 'gamma_e', 'gamma_p', 'mach', 'ip',
    'sdlnnedr', 'sdlntedr',
    'sdlnnidr_1', 'sdlnnidr_2', 'sdlnnidr_3', 'sdlnnidr_4', 'sdlnnidr_5',
    'sdlnnidr_6', 'sdlnnidr_7', 'sdlnnidr_8', 'sdlnnidr_9', 'sdlnnidr_10',
    'sdlntidr_1', 'sdlntidr_2', 'sdlntidr_3', 'sdlntidr_4', 'sdlntidr_5',
    'sdlntidr_6', 'sdlntidr_7', 'sdlntidr_8', 'sdlntidr_9', 'sdlntidr_10',
    'nuee'
]

vars_input_profiles_jbs = ['expro_rho', 'jbs_err', 'jbs_neo', 'jbs_sauter', 'jbs_nclass', 'jbs_koh']

#tuple: 1) latex 2) units 3) profiles column name
fancyNames = \
             {'NULL'       : (''                          ,''             ,'[null]'),
              'rho'        : ('{\\hat \\rho}'             ,''             ,'rho(-)'),
              'rmin'       : ('a'                         ,'m'            ,'rmin(m)'),
              'rmaj'       : ('R_0'                       ,'m'            ,'rmaj(m)'),
              'q'          : ('q'                         ,''             ,'q(-)'),
              'kappa'      : ('\\kappa'                   ,''             ,'kappa(-)'),
              'delta'      : ('\\delta'                   ,''             ,'delta(-)'),
              'Te'         : ('T_e'                       ,'keV'          ,'Te(keV)'),
              'ne'         : ('n_e'                       ,'10^{19}/m^3'  ,'ne(10^19/m^3)'),
              'z_eff'      : ('Z_\mathrm{eff}'            ,''             ,'zeff(-)'),
              'omega0'     : ('\\omega_0'                 ,'1/s'          ,'omega0(1/s)'),
              'flow_mom'   : ('S_\mathrm{\\omega}'        ,'Nm'           ,'flow_mom(Nm)'),
              'sbcx'       : ('sbcx'                      ,'1/m^3/s'      ,'sbcx(/m^3/s)'),
              'sbeame'     : ('sbeame'                    ,'1/m^3/s'      ,'sbeame(/m^3/s)'),
              'sscxl'      : ('sscxl'                     ,'1/m^3/s'      ,'sscxl(/m^3/s)'),
              'pow_e'      : ('P_e'                       ,'MW'           ,'pow_e(MW)'),
              'pow_i'      : ('P_i'                       ,'MW'           ,'pow_i(MW)'),
              'pow_ei'     : ('P_{ei}'                    ,'MW'           ,'pow_ei(MW)'),
              'zeta'       : ('\\zeta'                    ,''             ,'zeta(-)'),
              'flow_beam'  : ('S_\mathrm{n,beam}'         ,'kW/eV'        ,'flow_beam(kW/eV)'),
              'flow_wall'  : ('S_\mathrm{n,wall}'         ,'kW/eV'        ,'flow_wall(kW/eV)'),
              'zmag'       : ('Z_0'                       ,'m'            ,'zmag(m)'),
              'ptot'       : ('p_\mathrm{total}'          ,'Pa'           ,'ptot(Pa)'),
              'polflux'    : ('\\psi'                     ,'Wb/rad'       ,'polflux(Wb/rad)'),
              'pow_e_aux'  : ('P_{e,\\rm aux}'            ,'MW'           ,'pow_e_aux(MW)'),
              'pow_i_aux'  : ('P_{i,\\rm aux}'            ,'MW'           ,'pow_i_aux(MW)'),
              'pow_e_fus'  : ('P_{e,\\rm fus}'            ,'MW'           ,'pow_e_fus(MW)'),
              'pow_i_fus'  : ('P_{i,\\rm fus}'            ,'MW'           ,'pow_i_fus(MW)'),
              'pow_e_sync' : ('P_{e,\\rm sync}'           ,'MW'           ,'pow_e_sync(MW)'),
              'pow_e_brem' : ('P_{e,\\rm brem}'           ,'MW'           ,'pow_e_brem(MW)'),
              'pow_e_line' : ('P_{e,\\rm line}'           ,'MW'           ,'pow_e_line(MW)'),
              # extra (needs units)
              'bunit'      : ('B_\mathrm{unit}'           ,'T'            ,''),
              's'          : ('s'                         ,''             ,''),
              'drmaj'      : ('dR_0/dr'                   ,''             ,''),
              'dzmag'      : ('dZ_0/dr'                   ,''             ,''),
              'sdelta'     : ('s_\\delta'                  ,''             ,''),
              'skappa'     : ('s_\\kappa'                  ,''             ,''),
              'szeta'      : ('s_\\zeta'                   ,''             ,''),
              'dlnnedr'    : ('-dln(n_e)/dr'              ,'1/m'          ,''),
              'dlntedr'    : ('-dln(T_e)/dr'              ,'1/m'          ,''),
              'dlnptotdr'  : ('-dln(p_\mathrm{tot})/dr'   ,'1/m'          ,''),
              'drdrho'     : ('dr/d\\rho'                 ,''             ,''),
              'w0p'        : ('d(\\omega_0)/dr'           ,'1/s/m'        ,''),
              'vol'        : ('V'                         ,'m^3'          ,''),
              'volp'       : ('dV/dr'                     ,'m^2'          ,''),
              'cs'         : ('c_\mathrm{s}'              ,'m/s'          ,''),
              'rhos'       : ('\\rho_\mathrm{s,unit}'     ,'m'            ,''),
              'ni_new'     : ('n_i'                       ,'10^{19}/m^3'  ,''), #[Corrected for quasin.]
              'dlnnidr_new': ('-dln(n_i)/dr'              ,'1/m'          ,''), #[Corrected for quasin.]
              'grad_r0'    : ('|\\nabla_r|_{\\theta=0}'   ,''             ,''),
              'ave_grad_r' : ('<|\\nabla_r|>'             ,''             ,''),
              'bp0'        : ('B_p|_{\\theta=0}'          ,'T'            ,''),
              'bt0'        : ('B_t|_{\\theta=0}'          ,'T'            ,''),
              'gamma_e'    : ('r/q d(\\omega_0)/dr'       ,'1/s'          ,''),
              'gamma_p'    : ('R_0 d(\\omega_0)/dr'       ,'1/s'          ,''),
              'mach'       : ('R_0 \\omega_0/c_s'         ,''             ,''),
              'ip'         : ('I_p'                       ,'A'            ,''),
              'sdlnnedr'   : ("-n_e''/n_e"                ,'1/m^2'        ,''),
              'sdlntedr'   : ("-T_e''/T_e"                ,'1/m^2'        ,''),
              #jbs
              'expro_rho'  : ('\\rho'                     ,''             ,''),
              'jbs_err'    : ('j_\mathrm{bs,err}'         ,'MA/m^2'       ,''),
              'jbs_neo'    : ('j_\mathrm{bs,neo}'         ,'MA/m^2'       ,''),
              'jbs_sauter' : ('j_\mathrm{bs,sauter}'      ,'MA/m^2'       ,''),
              'jbs_nclass' : ('j_\mathrm{bs,nclass}'      ,'MA/m^2'       ,''),
              'jbs_koh'    : ('j_\mathrm{bs,koh}'         ,'MA/m^2'       ,''),
              }

for _k in range(1, 11):
    fancyNames['ni_%d'%_k]      =('n_{i,%d}'%_k            ,'10^{19}/m^3'  ,'ni_%d(10^19/m^3)'%_k)
    fancyNames['Ti_%d'%_k]      =('T_{i,%d}'%_k            ,'keV'          ,'Ti_%d(keV)'%_k)
    fancyNames['vtor_%d'%_k]    =('v_{\\varphi,%d}'%_k     ,'m/s'          ,'vtor_%d(m/s)'%_k)
    fancyNames['vpol_%d'%_k]    =('v_{\\theta,%d}'%_k      ,'m/s'          ,'vpol_%d(m/s)'%_k)
    fancyNames['dlntidr_%d'%_k] =('-dln(T_{i,%d})/dr'%_k   ,'1/m'          ,'')
    fancyNames['dlnnidr_%d'%_k] =('-dln(n_{i,%d})/dr'%_k   ,'1/m'          ,'')
    fancyNames['sdlntidr_%d'%_k]=("-T_{i,%d}''/T_{i,%d}"%(_k,_k)  ,'1/m^2'        ,'')
    fancyNames['sdlnnidr_%d'%_k]=("-n_{i,%d}''/n_{i,%d}"%(_k,_k)  ,'1/m^2'        ,'')

definition = \
             {'NULL' : 'Empty Column',
              'rho'  : 'Square root of the normalized toroidal flux\n'+
              '- Should span the region [0,1]\n'+
              '- Grid spacing can be nonuniform',
              #
              'rmin' : 'The generalized minor radius, r, in units of metres.\n'+
              '- see https://fusion.gat.com/theory/Gyrogeometry',
              #
              'polflux' : 'The polflux!',
              #
              'q' : 'The safety factor q (dimensionless)',
              #
              'omega0' : 'The omega0!',
              #
              'rmaj' : 'The generalized major radius, R_0(r), in units of metres.\n'+
              '- see https://fusion.gat.com/theory/Gyrogeometry',
              #
              'zmag' : 'The omega0!',
              #
              'kappa'     : 'The omega0!',
              #
              'delta'     : 'The omega0!',
              #
              'zeta'     : 'The omega0!',
              #
              'ne'     : 'The omega0!',
              #
              'Te'     : 'The omega0!',
              #
              'ptot'     : 'The omega0!',
              #
              'z_eff'     : 'The omega0!',
              }

for _k in range(1, 11):
    definition['ni_%d'%_k] = 'Ion %d density (10^19/m^3)'%_k
    definition['Ti_%d'%_k] = 'Ion %d temperature in keV'%_k
    definition['vtor_%d'%_k] = 'Ion %d toroidal velocity (m/s)'%_k
    definition['vpol_%d'%_k] = 'Ion %d poloidal velocity (m/s)'%_k



class profiles_genData:
    """
        profiles_gen output data class.
        """

    # ---------------------------------------------------------------------------#
    # Methods

    def __init__(self, infile):
        """
        Constructor reads in data from directory and creates new object.
        """

        import numpy as np

        # Initialize data
        self.infile = infile
        self.data = {}
        self.n_exp = 0
        self.geo = {}
        self.fancy = fancyNames
        
        # Read input.profiles
        try:
            tmp = profiles_gen(infile)
            self.data.update(tmp.data)
            self.n_exp = tmp.n_exp
            print('(INFO): (profiles_genData) ' + infile + ' found.')
        except Exception as E:
            raise (IOError('(ERROR): (profiles_genData) ' + infile + ' not found: ' + str(E)))

        # OPTIONAL: Read input.profiles.extra if it exists
        try:
            self.data.update(profiles_gen_extra(infile + '.extra').data)
            print('(INFO): (profiles_genData) ' + infile + '.extra found.')
        except Exception as E:
            print('(INFO): (profiles_genData) ' + infile + '.extra NOT loaded: ' + str(E))

        # OPTIONAL: Read input.profiles.geo if it exists
        try:
            tmp = profiles_gen_geo(infile + '.geo')
            self.geo.update(tmp.geo)
            self.nfourier = tmp.nfourier
            print('(INFO): (profiles_genData) ' + infile + '.geo found.')
        except Exception as E:
            print('(INFO): (profiles_genData) ' + infile + '.geo NOT loaded: ' + str(E))

        # OPTIONAL: Read input.profiles.jbs if it exists
        try:
            self.data.update(profiles_gen_jbs(infile + '.jbs').data)
            print('(INFO): (profiles_genData) ' + infile + '.jbs found.')
        except Exception as E:
            print('(INFO): (profiles_genData) ' + infile + '.jbs NOT loaded: ' + str(E))

class profiles_gen:
    def __init__(self, infile):
        import numpy as np

        self.data = {}
        self.shot = 0
        self.n_ion = 0
        self.n_exp = 0
        self.bt_exp = 0.0
        self.arho_exp = 0.0

        row = 0
        for line in open(infile, 'r').readlines():
            row = row + 1
            if line[0:5] == 'SHOT':
                self.shot = int(line.split('=')[1])
            if line[0:5] == 'N_ION':
                self.n_ion = int(line.split('=')[1])
            if line[0:5] == 'N_EXP':
                self.n_exp = int(line.split('=')[1])
            if line[0:6] == 'BT_EXP':
                self.bt_exp = float(line.split('=')[1])
            if line[0:8] == 'ARHO_EXP':
                self.arho_exp = float(line.split('=')[1])
            if line[0:7] == '#rho(-)':
                break

        data = np.loadtxt(infile, skiprows=row)

        for k1, line in enumerate(vars_input_profiles):
            if self.n_exp * (k1 + 1) > data.shape[0]:
                break
            for k2, var in enumerate(line):
                self.data[var] = data[self.n_exp * k1:self.n_exp * (k1 + 1), k2]
        if 'NULL' in self.data:
            del self.data['NULL']

class profiles_gen_extra:
    def __init__(self, infile):
        import numpy as np

        data = np.loadtxt(infile, comments='#')
        data = data.reshape((-1, len(vars_input_profiles_extra)), order='F')

        self.data = {}
        for k, var in enumerate(vars_input_profiles_extra):
            self.data[var] = data[:, k]

class profiles_gen_geo:
    def __init__(self, infile):
        import numpy as np

        # First, get number of Fourier modes
        fp = open(infile)
        for i, line in enumerate(fp):
            if i == 11:
                self.nfourier = int(line)
                break
        fp.close()
        
        data = np.loadtxt(infile, skiprows=12)

        # Dimension 9 assumes nfourier=8
        self.geo = {}
        data = data.reshape((4, self.nfourier + 1, -1), order='F')
        self.geo['ar'] = data[0, :, :]
        self.geo['br'] = data[1, :, :]
        self.geo['az'] = data[2, :, :]
        self.geo['bz'] = data[3, :, :]

class profiles_gen_jbs:
    def __init__(self, infile):
        import numpy as np

        data = np.loadtxt(infile, comments='#')
        data = data.reshape((-1, len(vars_input_profiles_jbs)), order='F')

        self.data = {}
        for k, var in enumerate(vars_input_profiles_jbs):
            self.data[var] = data[:, k]

class profiles_gen_format:
    def __init__(self):

        self.tags = vars_input_profiles
        self.desc = fancyNames
        self.defs = definition
        


if __name__ == '__main__':
    # check that all fancy names have been assigned
    import numpy

    tmp = set(numpy.array(vars_input_profiles).flatten())
    tmp = tmp.union(set(vars_input_profiles_jbs))
    tmp = tmp.union(set(vars_input_profiles_extra))
    print tmp.difference(set(fancyNames.keys()))

