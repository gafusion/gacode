#!/usr/bin/env python

#
# This tool upgrades the cgyro restart file version
# Given an original (input.gen,restart) pair
#   and the n_proc used to generate it
#   it generates a new restart file that keeps the old info
#
# Currently limited to v1 -> v2 upgrade
#

import string
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
            line = string.strip(line)

            # Skip blank lines
            if len(line) > 0 and line[0] != '#':
                (val,arg) = line.split()
                self.user_dict[arg] = val

    def getNSpecies(self):
        return int(self.user_dict["N_SPECIES"])
    
    def getMpiRankOrder(self):
        return int(self.user_dict["MPI_RANK_ORDER"])

old_dir="."
new_dir="t1/"

if len(sys.argv)!=4:
    print("ERROR: Wrong argument count")
    print("Usage:")
    print("  restart_upgrade.py org_dir new_dir n_proc")
    sys.exit(10)

old_dir=sys.argv[1]
new_dir=sys.argv[2]
n_proc =int(sys.argv[3])

try:
    old_cfg = CgyroInput()
    old_cfg.readInput(old_dir)
except IOError as err:
    print("IO error: {0}".format(err))
    sys.exit(21)
except:
    print("Unexpected error:", sys.exc_info()[0])
    sys.exit(22)

grid_obj =  cgyro_restart_resize.CGyroGrid()
grid_obj.load_from_dict(old_cfg.user_dict)

n_species = old_cfg.getNSpecies()
mpi_rank_order =  old_cfg.getMpiRankOrder()

if mpi_rank_order==1:
    try:
        cgyro_restart_resize.upgrade_v1v2_ro1(old_dir, new_dir, grid_obj, n_species, n_proc)
    except IOError as err:
        print("IO error: {0}".format(err))
        sys.exit(21)
    except:
        print("Unexpected error:", sys.exc_info()[0])
        sys.exit(22)
else:
    print("ERROR: MPI_RANK_ORDER %i not supported"%mpi_rank_order)
    sys.exit(21)
   


try:
  cgyro_restart_resize.add_species(old_dir, new_dir, grid_obj,
                                   diff_species["org_pre"],diff_species["org_post"],diff_species["new_species"])
except IOError as err:
    print("IO error: {0}".format(err))
    sys.exit(21)
