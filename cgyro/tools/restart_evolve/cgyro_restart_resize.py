#
# cgyro_restart_resize module
#
import sys,os

class CGyroGrid:
    def __init__(self):
        # nc
        self.n_radial = 0
        self.n_theta = 0

        # nv
        self.n_energy = 0
        self.n_xi = 0

        # toroidal
        self.n_toroidal = 0

    def load_from_dict(self,adict):
        self.n_radial =   int(adict["N_RADIAL"])
        self.n_theta =    int(adict["N_THETA"])
        self.n_energy =   int(adict["N_ENERGY"])
        self.n_xi =       int(adict["N_XI"])
        self.n_toroidal = int(adict["N_TOROIDAL"])

    def get_nc(self):
        return self.n_radial*self.n_theta

    def get_nvg(self):
        return self.n_energy*self.n_xi


def add_species(org_dir, new_dir, grid, org_pre_species, org_post_species, new_species):

    restart_fname="out.cgyro.restart"
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

    with open(org_fname,"rb") as org_fd:
        with  open(new_fname,"wb") as new_fd:
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
