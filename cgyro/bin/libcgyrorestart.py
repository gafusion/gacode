# restart_resize module

import sys,os
import struct

header_size = 1024*16 # 2xDouble precisio
restart_fname="bin.cgyro.restart"

class CGyroGrid:
    def __init__(self):
        # nc = n_radial*n_theta
        self.n_theta = 0
        self.n_radial = 0

        # nv = n_energy*n_xi*n_species
        self.n_xi = 0
        self.n_energy = 0
        self.n_species = 0

        # toroidal
        self.n_toroidal = 0

    def __str__(self):
        nc_str = "=== nc ===\nn_theta=%i\nn_radial=%i"%(self.n_theta,self.n_radial)
        nv_str = "=== nv ===\nn_xi=%i\nn_energy=%i\nn_species=%i"%(self.n_xi,self.n_energy,self.n_species)
        tor_str= "=== tor ===\nn_toroidal=%i"%self.n_toroidal
        return "%s\n%s\n%s"%(nc_str,nv_str,tor_str)

    def load_from_dict(self,adict):
        self.n_theta =    int(adict["N_THETA"])
        self.n_radial =   int(adict["N_RADIAL"])
        self.n_xi =       int(adict["N_XI"])
        self.n_energy =   int(adict["N_ENERGY"])
        self.n_species =  int(adict["N_SPECIES"])
        self.n_toroidal = int(adict["N_TOROIDAL"])

    def get_nc(self):
        return self.n_radial*self.n_theta

    def get_nv(self):
        return self.n_energy*self.n_xi*self.n_species

    def isSame(self,other,exclude_species=False):
        isSame=( (self.n_theta==other.n_theta) and
                 (self.n_xi==other.n_xi) and
                 (self.n_energy==other.n_energy) )
        if (not exclude_species):
           isSsame = isSame and (self.n_species==other.n_species)
        isSame= isSame and (
                 (self.n_radial==other.n_radial) and
                 (self.n_toroidal==other.n_toroidal) )
        return isSame

class CGyroRestartFormat:
    def __init__(self):
        self.el_size = 16 # 2xDouble precision (complex)
        # nv split used in restart formatting
        self.nv_loc = 0
        # toroidals_per_proc
        self.nt_loc = 0
        # velocity_order
        self.velocity_order = 0

    def __str__(self):
        return "nv_loc=%i\nnt_loc=%i\nvelocity_order=%i"%(self.nv_loc,self.nt_loc,self.velocity_order)

    def load_from_dict(self,adict):
	# self.nv_loc not in dict, has to be computed
        if 'VELOCITY_ORDER' in adict:
           # some old versions did not have the option of changing the velocity order
           self.velocity_order = int(adict["VELOCITY_ORDER"])
        if 'TOROIDALS_PER_PROC' in adict:
           # some old versions did not have the option of changing toroidals_per_proc
           self.nt_loc = int(adict["TOROIDALS_PER_PROC"])

class CGyroRestartHeader:
    def __init__(self):
        self.grid = CGyroGrid()
        self.fmt = CGyroRestartFormat()
        self.info_only_keys = [ "nc_loc", "n_toroidal_procs", "n_proc" , "mpi_rank_order",
                   "delta_t_method", "nonlinear_flag", "n_jtheta", "nsplit", "nsplitA", "nsplitB",
                   "nup_theta", "nv", "nc" ]
        self.info_only = {}
        self.reset_info()

    def __str__(self):
        out="%s\n=== formatting ==="%self.grid
        out += "\n%s\n=== info ==="%self.fmt
        for k in self.info_only_keys:
            out += "\n%s=%i"%(k,self.info_only[k])
        return out

    def reset_info(self):
        for k in self.info_only_keys:
           self.info_only[k] = 0

    def load(self,fdir, allow_old=False):
        fname = os.path.join(fdir,restart_fname)
        with open(fname,"rb") as fd:
            magic_b = fd.read(4)
            [magic] = struct.unpack('i',magic_b)
            if (magic==140906808):
                self.fmt.velocity_order = 1
                read_velocity_order = False
            elif (magic==140916753):
                self.fmt.velocity_order = 2
                read_velocity_order = False
            elif (magic==140974129):
                read_velocity_order = True
            else:
                raise IOError("Wrong CGyroRestartHeader magic number %i"%magic)
            magic1=magic
            version_b = fd.read(4)
            [version] = struct.unpack('i',version_b)
            if (not(( (magic==140974129) and (version==3) ) or
                    ( (magic!=140974129) and (version==2) )
                   ) ):
                raise IOError("Unsupported CGyroRestartHeader version %i for magic %i"%(version,magic))
            if ((version!=3) and (not allow_old)):
                raise IOError("Not latest CGyroRestartHeader version")

            grid_b = fd.read(6*4)
            [self.grid.n_theta,self.grid.n_radial,
             self.grid.n_species,
             self.grid.n_xi,self.grid.n_energy,
             self.grid.n_toroidal] = struct.unpack('6i',grid_b)
           
            self.reset_info() # will keep to 0 anything I do not know
            if (version==3):
               # new format
               grid_c = fd.read(3*4)
               [self.fmt.velocity_order,
                self.fmt.nt_loc, self.fmt.nv_loc] = struct.unpack('3i',grid_c)
               read_velocity_order = False
               opt_b = fd.read(len(self.info_only_keys)*4)
               opt_up = struct.unpack('%ii'%len(self.info_only_keys),opt_b)
               i=0
               for k in self.info_only_keys:
                  self.info_only[k] = opt_up[i]
                  i+=1
            else:
               # assume old format 
               self.fmt.nt_loc = 0 # unknown
               self.fmt.nv_loc = 0 # unknown
               mpi_b = fd.read(2*4)
               [self.info_only["mpi_rank_order"],self.info_only["n_proc"]] = struct.unpack('2i',mpi_b)

            magic_b = fd.read(4)
            [magic] = struct.unpack('i',magic_b)
            if (magic!=magic1):
                raise IOError("Wrong CGyroRestartHeader magic(2) number %i!=%i"%(magic,magic1))

    def savev3(self,fd):
            fd.seek(0)
            magic = 140974129
            magic_b= struct.pack('i',magic)
            fd.write(magic_b)
            version_b= struct.pack('i',3)
            fd.write(version_b)
            gridarr=[self.grid.n_theta,self.grid.n_radial,
                     self.grid.n_species,
                     self.grid.n_xi,self.grid.n_energy,
                     self.grid.n_toroidal,
                     self.fmt.velocity_order,
                     self.fmt.nt_loc, self.fmt.nv_loc]
            for g in gridarr:
                g_b = struct.pack('i',g)
                fd.write(g_b)
            for k in self.info_only_keys:
                v_b = struct.pack('i',self.info_only[k])
                fd.write(v_b)
            fd.write(magic_b)

    def get_ncbytes(self):
        return self.fmt.el_size*self.grid.get_nc()

    def get_total_bytes(self):
        return ( header_size +
                self.get_ncbytes()*self.grid.get_nv()*self.grid.n_toroidal)

    def get_iv(self, i_s, i_x, i_e):
        iv = 0
        if self.fmt.velocity_order==1:
            iv = (i_e*self.grid.n_xi+i_x)*self.grid.n_species+ i_s;
        else:
            iv = (i_s*self.grid.n_energy+i_e)*self.grid.n_xi+ i_x;
        return iv

    # offset to start of nc block
    # not including the header
    def nc_offset(self, i_s, i_x, i_e, i_t):
        nv = self.grid.get_nv()
        nv_loc = self.fmt.nv_loc
        nt_loc = self.fmt.nt_loc
        iv = self.get_iv(i_s,i_x,i_e)
        off_block = (iv%nv_loc) + (i_t%nt_loc)*nv_loc
        off_row = (iv//nv_loc)*(nv_loc*nt_loc)
        off_col = (i_t//nt_loc)*(nv*nt_loc)
        offset = self.get_ncbytes()*(off_block+off_row+off_col)
        return offset

