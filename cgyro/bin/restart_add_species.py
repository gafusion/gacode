#!/usr/bin/env python3
#
# This tool evolves the cgyro restart file
# Given an original (input.gen,restart) pair
# and a new input.gen,it generates a new restart file
# that keeps the old info
#
# Currently limited to adding an additional species
#

import sys,os
import argparse
import libcgyrorestart

class CgyroInput:
    """Input parser for input.cgyro.gen file"""
    def __init__(self):
        self.user_dict = {}
        self.def_fname = "input.cgyro.gen"

    def readInput(self,fdir):
        inputfile=os.path.join(fdir,self.def_fname)
        with open(inputfile,'r') as fin:
            in_lines = fin.readlines()
        for line in in_lines:
            # Remove leading and trailing whitespace from line
            line = line.strip()

            # Skip blank lines
            if len(line) > 0 and line[0] != '#':
                (val,arg) = line.split()
                self.user_dict[arg] = val

    def isSameGrid(self, other):
         cmpArgs=["N_ENERGY","N_XI","N_THETA","N_RADIAL","N_TOROIDAL"]
         isSame=True
         for arg in cmpArgs: isSame = isSame and (self.user_dict[arg]==other.user_dict[arg])
         return isSame

    # i is 0-based
    def isSameSpecEl(self, other, myi, otheri):
        isSame=True

        myarg="Z_%i"%(myi+1)
        otherarg="Z_%i"%(otheri+1)
        if (myarg in self.user_dict and otherarg in other.user_dict):
            isSame = isSame and (self.user_dict[myarg]==other.user_dict[otherarg])
        else:
            isSame=False

        myarg="MASS_%i"%(myi+1)
        otherarg="MASS_%i"%(otheri+1)
        if (myarg in self.user_dict and otherarg in other.user_dict):
            isSame = isSame and (self.user_dict[myarg]==other.user_dict[otherarg])
        else:
            isSame=False

        return isSame

    def isSameSpecies(self, other):
        my_n_species = int(self.user_dict["N_SPECIES"])
        other_n_species = int(other.user_dict["N_SPECIES"])
        if (my_n_species!=other_n_species):
            return False

        # Number not enough, check also type
        isSame=True
        for i in range(my_n_species):
            isSame = isSame and self.isSameSpecEl(other,i,i)

        return isSame

    # if true, also return pre, post and new species
    def isSpeciesSuperset(self, other, diff_species):
        my_n_species = int(self.user_dict["N_SPECIES"])
        other_n_species = int(other.user_dict["N_SPECIES"])
        if (my_n_species<=other_n_species):
            return False

        # find out how many of the species we had before
        org_pre_species=0
        for i in range(other_n_species):
            if self.isSameSpecEl(other,i,i):
                org_pre_species+=1
            else:
                break

        # find out how many new species we have
        new_species=0
        for i in range(org_pre_species,my_n_species):
            if ((i>=other_n_species) or (not self.isSameSpecEl(other,i,i))):
                new_species+=1
            else:
                break
        
        if (new_species!=(my_n_species-other_n_species)):
            return False # all new species must be together... looks like they are not

        # now check that all remaining species are the same
        org_post_species=0
        for i in range(org_pre_species+new_species,my_n_species):
            if (((i-new_species)<other_n_species) and self.isSameSpecEl(other,i,i-new_species)):
                org_post_species+=1
            else:
                return False # all post should have been the same

        diff_species["org_pre"]     = org_pre_species
        diff_species["org_post"]    = org_post_species
        diff_species["new_species"] = new_species

        return True

def get_arguments():

   parser=argparse.ArgumentParser(description="Utility to add species to CGYRO restart file")

   parser.add_argument('-o',
                       metavar='ORIG',
                       help="Original directory containing input.cgyro.gen and bin.cgyro.restart",
                       type=str,
                       required=True)
    
   parser.add_argument('-n',
                       metavar='NEW',
                       help="New directory containing input.cgyro.gen with added species",
                       type=str,
                       required=True)
    
   args=parser.parse_args()

   return args.o,args.n


def add_species(org_dir, new_dir, org_grid, new_grid, org_pre_species, org_post_species, new_species):
    org_fname = os.path.join(org_dir,libcgyrorestart.restart_fname)
    new_fname = os.path.join(new_dir,libcgyrorestart.restart_fname)

    header_size = libcgyrorestart.header_size

    org_header = libcgyrorestart.CGyroRestartHeader()
    org_header.load(org_dir)
    if (not org_grid.isSame(org_header.grid)):
        raise IOError("Wrong CGyroRestartHeader grid content")
    if (not new_grid.isSame(org_header.grid, exclude_species=True)):
        raise IOError("Incompatible CGyroRestartHeader grid content")
    if (org_header.grid.n_species!=(org_pre_species+org_post_species)):
        raise IOError("Incompatible CGyroRestartHeader number of species")
    ncbytes = org_header.get_ncbytes()
    org_fsize = org_header.get_total_bytes()
    if os.stat(org_fname).st_size!=org_fsize:
        raise IOError("Wrong restart file size")

    zerobuf = bytearray(ncbytes)
    for i in range(ncbytes):
        zerobuf[i] = 0  # 0.0 is also a binary 0

    new_header = libcgyrorestart.CGyroRestartHeader()
    new_header.load(org_dir)  # keep the same format as the old one
    new_header.grid.n_species += new_species # just add the species
    new_fsize = new_header.get_total_bytes()
    with open(org_fname,"rb") as org_fd:
        with  open(new_fname,"wb") as new_fd:
          new_fd.truncate(new_fsize)
          new_header.savev3(new_fd)

          for i_t in range(org_grid.n_toroidal):
           for i_e in range(org_grid.n_energy):
            for i_x in range(org_grid.n_xi):
              # copy the initial species
              for i_s in range(org_pre_species):
                org_off = org_header.nc_offset(i_s, i_x, i_e, i_t)
                new_off = new_header.nc_offset(i_s, i_x, i_e, i_t)
                org_fd.seek(header_size+org_off)
                new_fd.seek(header_size+new_off)
                tmp=org_fd.read(ncbytes)
                new_fd.write(tmp)

              # fill additional species with 0
              for i in range(org_post_species):
                i_s = org_pre_species+new_species + i
                new_off = new_header.nc_offset(i_s, i_x, i_e, i_t)
                new_fd.seek(header_size+new_off)
                new_fd.write(zerobuf)

              # copy the final species
              for i in range(org_post_species):
                i_s = org_pre_species+new_species + i
                org_off = org_header.nc_offset(i_s, i_x, i_e, i_t)
                new_off = new_header.nc_offset(i_s, i_x, i_e, i_t)
                org_fd.seek(header_size+org_off)
                new_fd.seek(header_size+new_off)
                tmp=org_fd.read(ncbytes)
                new_fd.write(tmp)

    return

old_dir,new_dir=get_arguments()

try:
    old_cfg = CgyroInput()
    new_cfg = CgyroInput()

    old_cfg.readInput(old_dir)
    new_cfg.readInput(new_dir)
except IOError as err:
    print("IO error: {0}".format(err))
    sys.exit(21)
except:
    print("Unexpected error:", sys.exc_info()[0])
    sys.exit(22)

if (not new_cfg.isSameGrid(old_cfg)):
    print("ERROR: Grids not the same")
    sys.exit(11)

if (new_cfg.isSameSpecies(old_cfg)):
    print("INFO: No changes detected. No restart file generated")
    sys.exit(1)

diff_species={"org_pre":0,"org_post":0,"new_species":0}
if (not new_cfg.isSpeciesSuperset(old_cfg,diff_species)):
    print("ERROR: New species are not a superset of the old ones")
    sys.exit(11)

if diff_species["org_pre"]>0:
  print("INFO: Keeping first %i species"%diff_species["org_pre"])
print("INFO: Adding %i species"%diff_species["new_species"])
if diff_species["org_post"]>0:
  print("INFO: Keeping last %i species"%diff_species["org_post"])


old_grid_obj =  libcgyrorestart.CGyroGrid()
new_grid_obj =  libcgyrorestart.CGyroGrid()
old_grid_obj.load_from_dict(old_cfg.user_dict)
new_grid_obj.load_from_dict(new_cfg.user_dict)

try:
  add_species(old_dir, new_dir, old_grid_obj, new_grid_obj,
              diff_species["org_pre"],diff_species["org_post"],diff_species["new_species"])
except IOError as err:
    print("IO error: {0}".format(err))
    sys.exit(21)
