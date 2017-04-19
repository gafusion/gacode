import os
import sys
import numpy as np

class tgyrodata:
    """TGYRO output data class"""

    #---------------------------------------------------------------------------#
    # Methods

    def __init__(self, sim_directory, verbose=True):
        """Constructor reads data from sim_directory"""

        self.verbose      = verbose
        self.loc_n_ion    = 0
        self.n_iterations = 0
        self.n_fields     = 0
        self.n_r          = 0
        self.data         = {}
        self.pedflag      = 0
    
        self.dir = sim_directory
        self.n_ion = int(self.get_tag_value("LOC_N_ION"))
        self.getcontrol()

        # Manage new (or old) data output formats
        if os.path.isfile(self.dir+'/out.tgyro.profile_i1'):
            if self.verbose: print "INFO: (data.py) Detected new multi-ion version of TGYRO [GOOD]."
            self.getdata()
        else:
            if self.verbose: print "INFO: (data.py) Detected old static ion version of TGYRO."
            if self.verbose: print "      --------> Using getolddata function ..."
            self.getolddata() 
        
    #---------------------------------------------------------------------------#

    def getolddata(self):
        """Reconstruct new datafiles from old dataset"""

        self.fileparser('out.tgyro.alpha')
        self.fileparser('out.tgyro.geometry.1')
        self.fileparser('out.tgyro.geometry.2')
        self.fileparser('out.tgyro.gyrobohm')
        self.fileparser('out.tgyro.nu_rho',onerow=True)
        self.fileparser('out.tgyro.power_e')
        self.fileparser('out.tgyro.power_i')

        self.fileparser('out.tgyro.flux_e')
        self.fileparser('out.tgyro.flux_i')
        self.fileparser('out.tgyro.flux_target')
        self.fileparser('out.tgyro.gradient')
        self.fileparser('out.tgyro.mflux_e')
        self.fileparser('out.tgyro.mflux_i')
        self.fileparser('out.tgyro.profile')
        self.data['ni1'] = self.data['ni']
        self.data['a/Lni1'] = self.data['a/Lni']
        self.data['ti1'] = self.data['ti']
        self.data['a/Lti1'] = self.data['a/LTi']
        self.data['pflux_i1_neo'] = self.data['pflux_i_neo']
        self.data['pflux_i1_tur'] = self.data['pflux_i_tur']
        self.data['pflux_i1_tot'] = self.data['pflux_i_neo']+self.data['pflux_i_tur']
        self.data['pflux_i1_target'] = 0.0*self.data['pflux_i_neo']
        self.data['eflux_i1_neo'] = self.data['eflux_i_neo']
        self.data['eflux_i1_tur'] = self.data['eflux_i_tur']
        for i in range(2,self.n_ion+1):
            self.fileparser('out.tgyro.profile'+str(i))
            self.data['ni'+str(i)] = self.data['ni']
            self.data['a/Lni'+str(i)] = self.data['a/Lni']
            self.fileparser('out.tgyro.flux_i'+str(i))
            self.fileparser('out.tgyro.mflux_i'+str(i))
            self.data['pflux_i'+str(i)+'_neo'] = self.data['pflux_i_neo']
            self.data['pflux_i'+str(i)+'_tur'] = self.data['pflux_i_tur']
            self.data['eflux_i'+str(i)+'_neo'] = self.data['eflux_i_neo']
            self.data['eflux_i'+str(i)+'_tur'] = self.data['eflux_i_tur']
            self.data['pflux_i'+str(i)+'_tot'] = 0.0*self.data['pflux_i_neo']
            self.data['pflux_i'+str(i)+'_target'] = 0.0*self.data['pflux_i_neo']
        self.data['a/Lte']  = self.data['a/LTe']

    def getdata(self):
        """Read all tgyro-format datafiles"""
        
        self.fileparser('out.tgyro.alpha')
        self.fileparser('out.tgyro.geometry.1')
        self.fileparser('out.tgyro.geometry.2')
        self.fileparser('out.tgyro.gyrobohm')
        self.fileparser('out.tgyro.nu_rho')
        self.fileparser('out.tgyro.power_e')
        self.fileparser('out.tgyro.power_i')

        self.fileparser('out.tgyro.residual')
        self.fileparser('out.tgyro.evo_er')
        self.fileparser('out.tgyro.evo_te')
        self.fileparser('out.tgyro.evo_ti')
        self.fileparser('out.tgyro.profile')

        # Species series 
        self.fileparser('out.tgyro.evo_ne')
        self.fileparser('out.tgyro.flux_e')
        self.fileparser('out.tgyro.profile_e')
        self.pedflag = self.fileparser('out.tgyro.ped',ped=True)
        for i in range(self.n_ion):
            for root in ['evo_n','flux_i','profile_i']:
                self.fileparser('out.tgyro.'+root+str(i+1))

        if self.verbose: print self.data.keys()

    #---------------------------------------------------------------------------#

    def get_tag_value(self, tag):
        """
        Return the specified variable from input.tgyro.gen.
        tag = input.tgyro tag
        """

        datafile = file(self.dir+'/input.tgyro.gen','r')

        for line in datafile:
            try:
                if line.split()[1] == tag:
                    return float(line.split()[0])
            except IndexError:
                print "WARNING: Cannot find specified input parameter: ", tag
                return 0

    #---------------------------------------------------------------------------#

    def getcontrol(self):
        """
        Read control.out to set resolutions.
        """

        import numpy as np
    
        data = np.loadtxt(self.dir+'/out.tgyro.control')

        self.n_r = int(data[0])
        self.n_evolve     = int(data[1])
        self.n_iterations = int(data[2])
   
    #---------------------------------------------------------------------------#

    def fileparser(self,file,onerow=False,ped=False):
        """Generic parser for standard TGYRO iteration-based file format"""

        import os
        import string
        import numpy as np

        if not os.path.exists(self.dir+'/'+file) :
            print "WARNING: (data.py) "+file+" does not exist."
            return 0

        data = open(self.dir+'/'+file,'r').readlines()

        if not len(''.join(data).strip()):
            print "WARNING: (data.py) "+file+" is empty."
            return 0

        # Data dimensions
        if onerow == False:
            ix = 2
        else:
            ix = 1

        if ped == False:
            nrad = self.n_r
        else:
            nrad = 1

        nr = nrad+ix
        nc = len(string.split(data[0]))
        nb = self.n_iterations+1

        numdata = np.zeros((nc,nb,nrad))
        for ib in range(nb):
            try:
                tags = string.split(data[ib*nr])
            except:
                print "WARNING: (data.py) "+file+" shorter than expected."
                return 0

            for ir in range(nrad):
                row = string.split(data[ib*nr+ir+ix])
                numdata[:,ib,ir] = row[:]

        # Populate data list
        for ic in range(nc):
            self.data[tags[ic]] = numdata[ic,:,:]
        
        if self.verbose: print 'INFO: (data.py) Read data in '+file
        return 1
