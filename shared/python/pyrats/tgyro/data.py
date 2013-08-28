import numpy as np
class TGYROData:
    """TGYRO output data class.

    Data:

    loc_n_ion
    tgyro_mode
    n_iterations
    n_fields
    n_radial
    dirname = ""
    data

    Example Usage:
    >>> from pyrats.tgyro.data import TGYROData
    >>> sim = TGYROData('examples')
    """

    #---------------------------------------------------------------------------#
    # Methods

    def __init__(self, sim_directory):
        """Constructor reads in data from sim_directory and creates new object.
        """
        
        self.init_data()
        self.set_directory(sim_directory)
        self.read_data()

    #---------------------------------------------------------------------------#

    def init_data(self):
        """Initialize object data."""

        self.loc_n_ion    = 0
        self.tgyro_mode   = 0
        self.n_iterations = 0
        self.n_fields     = 0
        self.n_radial     = 0
        self.data         = {}

    #---------------------------------------------------------------------------#

    def set_directory(self, path):
        """Set the simulation directory."""

        from os.path import expanduser, expandvars
        self.dirname = expanduser(expandvars(path))

    #---------------------------------------------------------------------------#

    def read_data(self):
        """Read in object data."""
        
        self.tgyro_mode = self.get_tag_value("TGYRO_MODE")
        if self.tgyro_mode == 2:
            self.read_stabilities()
        else:
            self.loc_n_ion = int(self.get_tag_value("LOC_N_ION"))
            self.read_control()
            self.fileparser('gyrobohm.out')
            self.fileparser('flux_e.out')
            self.fileparser('flux_i.out')
            self.fileparser('flux_target.out')
            self.fileparser('mflux_target.out')
            self.fileparser('chi_e.out')
            self.fileparser('chi_i.out')
            self.fileparser('profile.out')
            self.fileparser('gradient.out')
            self.fileparser('out.tgyro.geometry.1')
            self.fileparser('out.tgyro.geometry.2')
            self.fileparser('out.tgyro.power')
            for i in range(2,self.loc_n_ion+1):
                si = '%d'%i
                for fn in ['profile','chi_i','mflux_i','flux_i']:
                    self.fileparser(fn+si+'.out',spec_num=i)
            self.get_residual()
    #---------------------------------------------------------------------------#

    def get_tag_value(self, tag):
        """
        Return the specified variable from input.tgyro.gen.
        tag = input.tgyro tag
        """

        datafile = file(self.dirname+'/input.tgyro.gen','r')

        for line in datafile:
            try:
                if line.split()[1] == tag:
                    return float(line.split()[0])
            except IndexError:
                print "Cannot find specified input parameter: ", tag
                return 0

    #---------------------------------------------------------------------------#

    def read_control(self):
        """
        Read control.out to set resolutions.
        """

        import numpy as np
    
        data = np.loadtxt(self.dirname+'/control.out')

        self.n_radial     = int(data[0])
        self.n_fields     = int(data[1])
        self.n_iterations = int(data[2])

    #---------------------------------------------------------------------------#

    def fileparser(self,file,spec_num=''):
        """
        Generic parser for standard TGYRO file format.
        """

        import string
        import numpy as np

        data = open(self.dirname+'/'+file,'r').readlines()

        # Data dimensions
        nr = self.n_radial+2
        nc = len(string.split(data[0]))
        nb = self.n_iterations+1

        numdata = np.zeros((nc,nb,nr-2))
        for ib in range(nb):
            tags=string.split(data[ib*nr])
            null=string.split(data[ib*nr+1])
            for ir in range(nr-2):
                row=string.split(data[ib*nr+ir+2])
                for ic in range(nc):
                    numdata[ic,ib,ir] = row[ic]

        for ic in range(nc):
            if tags[ic]=='r/a' and spec_num!='':
                continue
            self.data[tags[ic]+str(spec_num)] = numdata[ic,:,:]

    #---------------------------------------------------------------------------#

    def get_residual(self):
        """Read residual.out
        """
        import string
        import numpy as np

        fn = 'residual.out'
        data = open(self.dirname+'/'+fn,'r').readlines()
        
        # Data dimensions
        nr = self.n_radial
        nb = self.n_iterations+1
        nc = int(1+2*self.get_tag_value("LOC_NE_FEEDBACK_FLAG") \
              +2*self.get_tag_value("LOC_TE_FEEDBACK_FLAG") \
              +2*self.get_tag_value("LOC_TI_FEEDBACK_FLAG"))
        
        numdata = np.zeros((nc,nb,nr-1),dtype=float)
        
        for ib in range(nb):
            tags=string.split(data[ib*nr]) # Contains overall residual
            for ir in range(nr-1):
                row=string.split(data[ib*nr+ir+1])
                for ic in range(nc):
                    numdata[ic,ib,ir] = row[ic]

        self.data['residual'] = numdata
        
    #---------------------------------------------------------------------------#

    def read_stab_file(self, file_name):
        """Read files generated with stability analysis mode.

        file_name - stability file, eg "wi_elec.out" (Default)

        Call with read_stab_file("wr_elec.out").

        Returns [r, ks, freq]

        """

        from numpy import loadtxt

        full_path = self.dirname + '/' + file_name
        raw_data = loadtxt(file(full_path), skiprows=2, unpack=True)
        r = raw_data[0,:]
        freqs = raw_data[1:,:]
        num_ks = freqs.shape[1]
        ks = loadtxt(file(full_path), skiprows=1, usecols=range(1,num_ks+1))[0]

        return [r, ks, freqs]

    def read_stabilities(self):
        """Read output files from TGYRO_METHOD=2, store into variables.

        Files read:
            wr_ion.out   ->    self.wr_ion
            wi_ion.out   ->    self.wi_ion
            wr_elec.out  ->    self.wr_elec
            wi_elec.out  ->    self.wi_elec

        """

        self.data['wr_ion'] = self.read_stab_file("wr_ion.out")
        self.data['wi_ion'] = self.read_stab_file("wi_ion.out")
        self.data['wr_elec'] = self.read_stab_file("wr_elec.out")
        self.data['wi_elec'] = self.read_stab_file("wi_elec.out")

        self.data['r/a'] = self.data['wr_ion'][0]
        self.n_radial = len(self.data['r/a'])
    
    # ----------------------------------------- #
    # Get data back

    def get_local_res(self, field=1, iteration=-1):
        try:
            return self.data['local_res'][iteration][:,2*field-1]
        except IndexError:
            print "Invalid field residual."
            print self.data['local_res']
            return None

    def get_stability_at_radius(self, radius=0, frequency='r', direction='ion'):
        """Get frequency vs. ky at specified radius from stability analysis.

        Parameters:

            radius        -       index of requested radius (Default: 0)
                                  integer ranging from 0 : n_r-1

            frequency     -       real or imaginary spectrum
                                  'r'     : real (Default)
                                  'i'     : imaginary

            direction     -       direction of spectrum
                                  'ion'   : ion direction (Default)
                                  'elec'  : electron direction

        Usage:

            [ky, wi_ion] = get_stability_at_radius(0, 'i', 'ion')

        """

        # Check for valid inputs


        if radius < 0 or radius > self.n_radial-1:
            print "Error in get_stability_at_radius. Invalid radial point: ", radius
            print "Valid range: 0 - ", self.n_radial-1
            return

        if frequency is not 'r' and frequency is not 'i':
            print "Error in get_stability_at_radius. Invalid frequency: ", frequency
            print "Valid options are 'r' and 'i'"
            return

        if direction is not 'ion' and direction is not 'elec':
            print "Error in get_stability_at_radius. Invalid frequency: ", frequency
            print "Valid options are 'ion' and 'elec'"
            return

        # Return requested data
        if direction is 'ion':
            if frequency is 'r':
                k = self.data['wr_ion'][1]
                f = self.data['wr_ion'][2][:, radius]
            else:
                k = self.data['wi_ion'][1]
                f = self.data['wi_ion'][2][:, radius]
        else:
            if frequency is 'r':
                k = self.data['wr_elec'][1]
                f = self.data['wr_elec'][2][:, radius]
            else:
                k = self.data['wi_elec'][1]
                f = self.data['wi_elec'][2][:, radius]

        return [k, f]

    def get_most_unstable_at_radius(self, radius=0, direction='ion'):
        """Return the most unstable mode at specified radius in the specified
        direction.

        Parameters:

           radius     -   requested radial index (Default: 0)
                          integer between 0 and n_r-1

           direction  -   'ion'  most unstable ion mode
                          'elec' most unstable electron mode

        Usage:

           [k_e, omega_e, gamma_e] = get_most_unstable_at_radius(0, 'elec')

        """

        # Check input sanity

        if radius < 0 or radius > self.n_radial-1:
            print "Error in get_most_unstable_at_radius. Invalid radial point: ", radius
            print "Valid range: 0 - ", self.n_radial-1
            return

        if direction is not 'ion' and direction is not 'elec':
            print "Error in get_most_unstable_at_radius. Invalid frequency: ", frequency
            print "Valid options are 'ion' and 'elec'"
            return

        # Get ion and electron data at requested radius

        ks, gammas = self.get_stability_at_radius(radius, 'i', direction)
        k2, omegas = self.get_stability_at_radius(radius, 'r', direction)

        gamma_max = max(gammas)
        max_index = gammas.argmax()
        k_max = ks[max_index]
        omega_max = omegas[max_index]
        
        return k_max, omega_max, gamma_max

# Analysis Tools

    def make_gradient_vs_field_space(self, r=1, evolve_field=1, grad='a/LTi', profile='ti'):
        """Return matrix of Space[iteration][gradient,profile,residual] at given radial point.

        eg: Space[[4.0, 1.7, 2.5], [3.6, 1.5, 0.1]] if grad 4.0->3.6 while Ti->1.5
            and local residual went from 2.5 -> 0.1

        """
        import numpy
        space=[]
        for iteration in range(self.n_iterations + 1):
            #print "iteration", iteration, "r", r
            loc_grad = self.data[grad][iteration][r]
            loc_prof = self.data[profile][iteration][r]
            loc_resi = self.get_local_res(evolve_field, iteration)[r]
            #print loc_grad, loc_prof, loc_resi
            space.append([loc_grad, loc_prof, loc_resi])

        return numpy.array(space)             

    def make_gradient_vs_field_plot(self, r=1, evolve_field=1, grad='a/LTi', profile='ti', arrows=False):
        """Return a pyplot figure of gradient_vs_field_space."""
        import matplotlib.pyplot as plt
        from pylab import Arrow
        space = self.make_gradient_vs_field_space(r, evolve_field, grad, profile)
        fig = plt.figure()
        walk = fig.add_subplot(111)
        x=space[:,0]
        y=space[:,1]
        z=space[:,2]
        cax = walk.scatter(x,y,c=z,s=40,vmin=min(z),vmax=max(z))
        walk.set_xlabel(grad)
        walk.set_ylabel(profile)
        cbar = plt.colorbar(cax)
        cbar.set_label("Local Residual")
        walk.set_title("Radial Point " + str(r))
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        deltax = xmax - xmin
        deltay = ymax - ymin

        if arrows:
            arrow_width = 0.05 * deltax
            for point in range(len(x)-1):
                dx = x[point+1] - x[point]
                dy = y[point+1] - y[point]
                arrow = Arrow(x[point], y[point], dx, dy, alpha=0.75, width=arrow_width, \
                              linewidth=0, edgecolor='none', facecolor='black')
                walk.add_patch(arrow)

        for i in range(len(x)-1):
            if i % 5 == 0:
                walk.annotate(str(i), (x[i], y[i]), xycoords='data', size=15)           
        walk.annotate("Start", (x[0], y[0]), xycoords='data', xytext=((xmax - 0.3 * deltax), (ymin + 0.05 * deltay)), \
                      textcoords='data', size=15, \
                      arrowprops=dict(width=3, headwidth=5, shrink=0.05),)
        walk.annotate("End, Iteration "+ str(len(x)-1), (x[-1], y[-1]), xycoords='data', xytext=((xmin + 0.05 * deltax), \
                      (ymax - 0.1 * deltay)), \
                      textcoords='data', size=15, \
                      arrowprops=dict(width=3, headwidth=5, shrink=0.05, facecolor="red"),)
        minpos = 0
        for j in range(len(z-1)):
            if z[j] == min(z):
                minpos = j
        walk.annotate("Min Residual = " + str(z[minpos]), (x[minpos], y[minpos]), xycoords='data', \
                      xytext=((xmin + 0.05 * deltax), (ymin + 0.05 * deltay)), \
                      textcoords='data', size=15, arrowprops=dict(width=3, headwidth=5, shrink=0.05, facecolor="black"))

        return fig
    
    def sprofile(self, what, nf0=201, x='r/a'):
        """
        This function returns smooth profiles on a uniform r/a grid
        @param what: what profile to return ['r/a', 'rho', 'rmin', 'rmaj/a', 'te', 'a/LTe', 'ne', 'a/Lne' , 'ti', 'a/LTi', 'ni', 'a/Lni', 'M=wR/cs', 'M=wR/cs']
        @param nf0: number of points
        @param x: return profiles equally spaced in 'r/a', 'rho', 'rmin', 'rmaj/a'
        @param verbose: plotting of the `what` quantity
        @return: `what` quantity (niterations x nf0) or `x` at the locations
        """
        bc_offset = self.get_tag_value('LOC_BC_OFFSET')
        nr = len(self.data[x][0,:])
        r_ref = self.data[x][0,nr-bc_offset-1]
        n   = self.n_iterations+1
        rf00 = np.linspace(self.data[x][0,0],self.data[x][0,-1],int(nf0))
        rf0 = np.array(sorted(set(list(rf00)+list(self.data[x][0]))))
        rf00_ind = np.zeros(len(rf0),dtype=bool)
        for i in range(len(rf0)):
            if rf0[i] in rf00:
                rf00_ind[i] = True
        r_ref_ind = np.nonzero(rf0==r_ref)[0]
        if what == x:
            return rf0[rf00_ind]
        elif what in ['r/a', 'rho','rmin','rmaj/a']:
            return np.interp(rf0[rf00_ind],self.data[x][0],self.data[what][0])

        quantity={'ne':'a/Lne', 'ni':'a/Lni', 'te':'a/LTe', 'ti':'a/LTi'}#, 'M=wR/cs':'M=wR/cs'}
        quantityScaleLen={'a/Lne':'ne', 'a/Lni':'ni', 'a/LTe':'te', 'a/LTi':'ti'}#, 'M=wR/cs':'M=wR/cs'}
        for spec in range(2,self.loc_n_ion+1):
            for quant in ['T','n']:
                quantity['%si%d'%(quant,spec)] = 'a/L%si%d'%(quant,spec)
                quantityScaleLen['a/L%si%d'%(quant,spec)] = '%si%d'%(quant,spec)
        
        if what in quantity:
            whatScale=quantity[what]
        elif what in quantityScaleLen:
            whatScale=what
            what=None
        else:
            raise(Exception("Quantity '"+what+"' is not in %s"%repr(quantity.keys()+quantity.values())))

        # Linearly interpolate gradient scale lenghts
        zf0 = np.zeros((n,len(rf0)))
        for l in range(n):
            zf0[l,:] = np.interp(rf0,self.data[x][0],self.data[whatScale][l,:])

        if what is not None:
            pf0 = np.zeros((n,len(rf0)))
            for l in range(n):
                # Set boundary condition at pivot point
                pf0[l,r_ref_ind] = self.data[what][l,nr-bc_offset-1]
            # Exponential integration to obtain smooth profiles
            for i in range(r_ref_ind,0,-1):
                pf0[:,i-1] = pf0[:,i]*np.exp( 0.5*(rf0[i]-rf0[i-1])*(zf0[:,i]+zf0[:,i-1]))
            for i in range(r_ref_ind+1,len(rf0)):
                pf0[:,i] = pf0[:,i-1]*np.exp(-0.5*(rf0[i]-rf0[i-1])*(zf0[:,i]+zf0[:,i-1]))
            
        else:
            what=whatScale
            pf0=zf0

        return pf0[:,rf00_ind]
