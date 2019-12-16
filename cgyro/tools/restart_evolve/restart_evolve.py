#!/usr/bin/env python

#
# This tool evolves the cgyro restart file
# Given an original (input.gen,restart) pair
#   and a new input.gen,
#   it generates a new restart file that keeps the old info
#
# Currently limited to adding an additional species
#

import sys,os

import cgyro_restart_resize

class CgyroInput:
    """Input parser for input.cgyro.gen file"""
    def __init__(self):
        self.user_dict = {}
        self.def_fname = "input.cgyro.gen"

    def readInput(self,fdir):
        inputfile=os.path.join(fdir,self.def_fname)
        for line in open(inputfile,'r').readlines():
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

old_dir="."
new_dir="t1/"

if len(sys.argv)!=3:
    print("ERROR: Wrong argument count")
    print("Usage:")
    print("  restart_evolve.py org_dir new_dir")
    sys.exit(10)

old_dir=sys.argv[1]
new_dir=sys.argv[2]
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

grid_obj =  cgyro_restart_resize.CGyroGrid()
grid_obj.load_from_dict(new_cfg.user_dict)

if (new_cfg.isSameSpecies(old_cfg)):
    print("INFO: No changes detected. No restart file generated")
    sys.exit(1)

diff_species={"org_pre":0,"org_post":0,"new_species":0}
if (not new_cfg.isSpeciesSuperset(old_cfg,diff_species)):
    print("ERROR: New species are not a superset of the old ones")
    sys.exit(11)

print("INFO: Adding %i species"%diff_species["new_species"])
try:
  cgyro_restart_resize.add_species(old_dir, new_dir, grid_obj,
                                   diff_species["org_pre"],diff_species["org_post"],diff_species["new_species"])
except IOError as err:
    print("IO error: {0}".format(err))
    sys.exit(21)
