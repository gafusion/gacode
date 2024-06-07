# restart_resize module

import sys,os
import struct

header_size = 1024*16 # 2xDouble precisio
restart_fname="bin.cgyro.restart"

class CGyroGrid:
    def __init__(self):
        # nc
        self.n_theta = 0
        self.n_radial = 0

        # nv
        self.n_xi = 0
        self.n_energy = 0

        # toroidal
        self.n_toroidal = 0

        # velocity_order
        self.velocity_order = 1

    def load_from_dict(self,adict):
        self.n_theta =    int(adict["N_THETA"])
        self.n_radial =   int(adict["N_RADIAL"])
        self.n_xi =       int(adict["N_XI"])
        self.n_energy =   int(adict["N_ENERGY"])
        self.n_toroidal = int(adict["N_TOROIDAL"])
        if 'VELOCITY_ORDER' in adict:
           # some old versions did not have the option of changing the velocity order
           self.velocity_order = int(adict["VELOCITY_ORDER"])

    def get_nc(self):
        return self.n_radial*self.n_theta

    def get_nvg(self):
        return self.n_energy*self.n_xi

    def isSame(self,other):
        return ( (self.n_theta==other.n_theta) and
                 (self.n_radial==other.n_radial) and
                 (self.n_xi==other.n_xi) and
                 (self.n_energy==other.n_energy) and
                 (self.n_toroidal==other.n_toroidal) and
                 (self.velocity_order==other.velocity_order))

class CGyroRestartHeader:
    def __init__(self):
        self.grid = CGyroGrid()
        self.n_species = 0
        self.mpi_rank_order = 0
        self.n_proc = 0

    def load(self,fdir):
        fname = os.path.join(fdir,restart_fname)
        with open(fname,"rb") as fd:
            magic_b = fd.read(4)
            [magic] = struct.unpack('i',magic_b)
            if (magic==140906808):
                self.grid.velocity_order = 1
            elif (magic==140916753):
                self.grid.velocity_order = 2
            else:
                raise IOError("Wrong CGyroRestartHeader magic number %i"%magic)
            magic1=magic
            version_b = fd.read(4)
            [version] = struct.unpack('i',version_b)
            if (version!=2):
                raise IOError("Unsupported CGyroRestartHeader version %i"%version)

            grid_b = fd.read(6*4)
            [self.grid.n_theta,self.grid.n_radial,
             self.n_species,
             self.grid.n_xi,self.grid.n_energy,
             self.grid.n_toroidal] = struct.unpack('6i',grid_b)
            
            mpi_b = fd.read(2*4)
            [self.mpi_rank_order,self.n_proc] = struct.unpack('2i',mpi_b)

            magic_b = fd.read(4)
            [magic] = struct.unpack('i',magic_b)
            if (magic!=magic1):
                raise IOError("Wrong CGyroRestartHeader magic(2) number %i!=%i"%(magic,magic1))

    def savev2(self,fdir):
        fname = os.path.join(fdir,restart_fname)
        with open(fname,"rb+") as fd:
            magic = 140906808 if (self.grid.velocity_order == 1) else 140916753
            magic_b= struct.pack('i',magic)
            fd.write(magic_b)
            version_b= struct.pack('i',2)
            fd.write(version_b)
            gridarr=[self.grid.n_theta,self.grid.n_radial,
                     self.n_species,
                     self.grid.n_xi,self.grid.n_energy,
                     self.grid.n_toroidal,
                     self.mpi_rank_order,self.n_proc]
            for g in gridarr:
                g_b = struct.pack('i',g)
                fd.write(g_b)
            fd.write(magic_b)

def add_species(org_dir, new_dir, grid, org_pre_species, org_post_species, new_species):
    org_fname = os.path.join(org_dir,restart_fname)
    new_fname = os.path.join(new_dir,restart_fname)

    tag_fname="out.cgyro.tag"
    new_tag_fname = os.path.join(new_dir,tag_fname)


    # Data structure (in fortran terms)
    # if velocity_order==1
    #   doublex2 h_x[nc,n_species,nvg,n_toroidal]
    # else
    #   doublex2 h_x[nc,nvg,n_species,n_toroidal]

    el_size = 16 # 2xDouble precision
    nc =      grid.get_nc()
    ncbytes = nc*el_size

    nvg =     grid.get_nvg()

    zerobuf = bytearray(ncbytes*new_species)
    for i in range(ncbytes*new_species):
        zerobuf[i] = 0  # 0.0 is also a binary 0

    org_header = CGyroRestartHeader()
    org_header.load(org_dir)
    # TODO: verify consistency
    if (not grid.isSame(org_header.grid)):
        raise IOError("Wrong CGyroRestartHeader grid content")
    if (org_header.n_species!=(org_pre_species+org_post_species)):
        raise IOError("Wrong CGyroRestartHeader number of species")

    with open(org_fname,"rb") as org_fd:
        with  open(new_fname,"wb") as new_fd:
            tmp=org_fd.read(header_size)
            new_fd.write(tmp)

            for t in range(grid.n_toroidal):
              if org_header.grid.velocity_order==1:
                if (org_pre_species>0):
                    tmp=org_fd.read(ncbytes*org_pre_species)
                    new_fd.write(tmp)

                new_fd.write(zerobuf)
                for n in  range(nvg-1):
                    # write the pre and post together
                    tmp=org_fd.read(ncbytes*(org_pre_species+org_post_species))
                    new_fd.write(tmp)

                    new_fd.write(zerobuf)

                if (org_post_species>0):
                    tmp=org_fd.read(ncbytes*org_post_species)
                    new_fd.write(tmp)
              else: # velocity_order==2
                if (org_pre_species>0):
                    # keep memory footprint small
                    for n in range(nvg):
                        tmp=org_fd.read(ncbytes*org_pre_species)
                        new_fd.write(tmp)

                # keep memory footprint small
                for n in range(nvg):
                    new_fd.write(zerobuf)

                if (org_post_species>0):
                    # keep memory footprint small
                    for n in range(nvg):
                        tmp=org_fd.read(ncbytes*org_post_species)
                        new_fd.write(tmp)


    new_header = org_header
    new_header.n_species = org_header.n_species + new_species
    new_header.savev2(new_dir)

    return
