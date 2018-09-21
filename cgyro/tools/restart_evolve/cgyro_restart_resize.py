#
# cgyro_restart_resize module
#
import sys,os
import struct

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

    def load_from_dict(self,adict):
        self.n_theta =    int(adict["N_THETA"])
        self.n_radial =   int(adict["N_RADIAL"])
        self.n_xi =       int(adict["N_XI"])
        self.n_energy =   int(adict["N_ENERGY"])
        self.n_toroidal = int(adict["N_TOROIDAL"])

    def get_nc(self):
        return self.n_radial*self.n_theta

    def get_nvg(self):
        return self.n_energy*self.n_xi

    def isSame(self,other):
        return ( (self.n_theta==other.n_theta) and
                 (self.n_radial==other.n_radial) and
                 (self.n_xi==other.n_xi) and
                 (self.n_energy==other.n_energy) and
                 (self.n_toroidal==other.n_toroidal))

class CGyroRestartHeader:
    def __init__(self):
        self.grid = CGyroGrid()
        self.n_species = 0

    def load(self,fdir):
        fname = os.path.join(fdir,restart_fname)
        with open(fname,"rb") as fd:
            magic_b = fd.read(4+12)
            [magic] = struct.unpack('i12x',magic_b)
            if (magic!=140906808):
                raise IOError("Wrong CGyroRestartHeader magic number %i"%magic)
            version_b = fd.read(4+12)
            [version] = struct.unpack('i12x',version_b)
            if (version!=2):
                raise IOError("Unsupported CGyroRestartHeader version %i"%version)

            grid_b = fd.read(6*(4+12))
            [self.grid.n_theta,self.grid.n_radial,
             self.n_species,
             self.grid.n_xi,self.grid.n_energy,
             self.grid.n_toroidal] = struct.unpack('i12xi12xi12xi12xi12xi12x',grid_b)
            


def add_species(org_dir, new_dir, grid, org_pre_species, org_post_species, new_species):
    header_size = 1024

    org_fname = os.path.join(org_dir,restart_fname)
    new_fname = os.path.join(new_dir,restart_fname)

    tag_fname="out.cgyro.tag"
    new_tag_fname = os.path.join(new_dir,tag_fname)


    # Data structure (in fortran terms)
    # double h_x[nc,n_species,nvg,2*n_toroidal]

    el_size = 8 # Double precision
    nc =      grid.get_nc()
    ncbytes = nc*el_size

    nvg =     grid.get_nvg()

    tor2 =    grid.n_toroidal*2

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

            for t in range(tor2):
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

    with open(new_tag_fname,"w") as tag_fd:
        tag_fd.write("           0\n 0.0000E+00\n")


    return
