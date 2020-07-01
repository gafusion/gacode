import os
import sys
import numpy as np
import re

# class for managing tgyro output data (out.tgyro.*)

class tgyrodata:

    def __init__(self, sim_directory, verbose=True):

        self.verbose      = verbose
        self.loc_n_ion    = 0
        self.n_iterations = 0
        self.n_fields     = 0
        self.n_r          = 0
        self.data         = {}
        self.pedflag      = 0

        self.dir = sim_directory
        self.n_ion = int(self.get_tag_value("LOC_N_ION"))

        data = np.genfromtxt(self.dir+'/out.tgyro.control')
        self.n_r         = int(data[0])
        self.n_evolve     = int(data[1])
        self.n_iterations = int(data[2])

        self.getdata()

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

        if self.verbose: print(list(self.data.keys()))


    def get_tag_value(self, tag):

        datafile = open(self.dir+'/input.tgyro.gen','r').readlines()

        for line in datafile:
            try:
                if line.split()[1] == tag:
                    return float(line.split()[0])
            except IndexError:
                print("WARNING: Cannot find specified input parameter: ", tag)
                return 0


    def fileparser(self,file,onerow=False,ped=False):
        """Generic parser for standard TGYRO iteration-based file format"""

        import os
        import numpy as np

        if not os.path.exists(self.dir+'/'+file) :
            print("WARNING: (data.py) "+file+" does not exist.")
            return 0

        data = open(self.dir+'/'+file,'r').readlines()

        if not len(''.join(data).strip()):
            print("WARNING: (data.py) "+file+" is empty.")
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
        nc = len(data[0].split())
        nb = self.n_iterations+1

        numdata = np.zeros((nc,nb,nrad))
        for ib in range(nb):
            try:
                tags = data[ib*nr].split()
            except:
                print("WARNING: (data.py) "+file+" shorter than expected.")
                return 0

            for ir in range(nrad):
                row = data[ib * nr + ir + ix].split()
                for ic in range(nc):
                    try:
                        numdata[ic, ib, ir] = float(row[ic])
                    except ValueError:
                        # parse badly formatted numbers like this 3.952525-323
                        tmp = re.match('([\-\+]*[\.0-9]+)([\-\+])([0-9]+)', row[ic])
                        if tmp.groups():
                            numdata[ic, ib, ir] = float('%sE%s%s' % (tmp.groups()[0], tmp.groups()[1], tmp.groups()[2]))
                        else:
                            raise

        # Populate data list
        for ic in range(nc):
            self.data[tags[ic]] = numdata[ic,:,:]

        if self.verbose: print('INFO: (data.py) Read data in '+file)
        return 1
