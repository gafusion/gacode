#!/usr/bin/env python3
#
# This tool evolves the cgyro restart file
# Given an original (input.gen,restart) pair
# and a new input.gen,it generates a new restart file
# that keeps the old info
#
# Currently limited to adding radial and toroidal dimensions
#

import sys,os
import argparse
import libcgyrorestart
import array

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
         cmpArgs=["N_ENERGY","N_XI","N_THETA","N_SPECIES"]
         isSame=True
         for arg in cmpArgs: isSame = isSame and (self.user_dict[arg]==other.user_dict[arg])
         return isSame

    def isSameDims(self, other):
         cmpArgs=["N_RADIAL","N_TOROIDAL"]
         isSame=True
         for arg in cmpArgs: isSame = isSame and (self.user_dict[arg]==other.user_dict[arg])
         return isSame

    # if true, also return pre, post and new species
    def isDimsSuperset(self, other):
        my_n_radial = int(self.user_dict["N_RADIAL"])
        my_n_toroidal = int(self.user_dict["N_TOROIDAL"])
        other_n_radial = int(other.user_dict["N_RADIAL"])
        other_n_toroidal = int(other.user_dict["N_TOROIDAL"])
        if (my_n_radial<other_n_radial):
            return False
        if (my_n_toroidal<other_n_toroidal):
            return False
        return True

def get_arguments():

   parser=argparse.ArgumentParser(description="Utility to add dimensions to CGYRO restart file")

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
    
   parser.add_argument('-s',
                       metavar='SCALE',
                       help="Scaling factor for exended elements (Default: 1.0)",
                       type=float,
                       default=1.0,
                       required=False)
    
   args=parser.parse_args()

   return args.o,args.n,args.s


def add_dims(org_dir, new_dir, org_grid, new_grid, dim_offset, scale_val):
    org_fname = os.path.join(org_dir,libcgyrorestart.restart_fname)
    new_fname = os.path.join(new_dir,libcgyrorestart.restart_fname)

    header_size = libcgyrorestart.header_size

    org_header = libcgyrorestart.CGyroRestartHeader()
    org_header.load(org_dir)
    if (not org_grid.isSame(org_header.grid)):
        raise IOError("Wrong CGyroRestartHeader grid content")
    thetabytes = org_header.get_thetabytes()
    org_fsize = org_header.get_total_bytes()
    if os.stat(org_fname).st_size!=org_fsize:
        raise IOError("Wrong restart file size")

    new_header = libcgyrorestart.CGyroRestartHeader()
    new_header.load(org_dir)  # keep the same format as the old one
    # Update the relevant dims
    new_header.grid = new_grid
    # nt_loc may not be compatible with new n_toroidal, just set to 1 to be safe
    new_header.fmt.nt_loc = 1
    # invalidate optional info, to maintain consistency
    new_header.reset_info()
    new_fsize = new_header.get_total_bytes()
    with open(org_fname,"rb") as org_fd:
        with  open(new_fname,"wb") as new_fd:
          new_fd.truncate(new_fsize)
          new_header.savev3(new_fd)

          for i_t in range(new_grid.n_toroidal):
           j_t = min(i_t, org_grid.n_toroidal-1) # will extend last element past the limit
           for i_e in range(new_grid.n_energy):
            for i_x in range(new_grid.n_xi):
             for i_s in range(new_grid.n_species):
              for i_r in range(new_grid.n_radial):
                new_off = new_header.theta_offset(i_r,
                                                  i_s, i_x, i_e,
                                                  i_t)
                new_fd.seek(header_size+new_off)
                if i_r<dim_offset:
                  j_r = 0 # use the lowest element for the left pad
                else:
                  j_r = min(i_r-dim_offset, org_grid.n_radial-1) # will extend the last eelement past the limit
                org_off = org_header.theta_offset(j_r,
                                                  i_s, i_x, i_e,
                                                  j_t)
                org_fd.seek(header_size+org_off)
                tmp=org_fd.read(thetabytes)
                if ( (j_t!=i_t) and ((j_r+dim_offset)!=i_r) and # only scale new elements
                     (scale_val!=1.0) ):  # and only if the scale is not 1
                  tarr=array.array('d', tmp)
                  tlist=tarr.tolist()
                  for i in range(len(tlist)):
                    tlist[i] *= scale_val
                  tarr=array.array('d')
                  tarr.fromlist(tlist)
                  tmp=tarr.tobytes()
                new_fd.write(tmp)

    return

old_dir,new_dir,scale_val=get_arguments()

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

if (new_cfg.isSameDims(old_cfg)):
    print("INFO: No changes detected. No restart file generated")
    sys.exit(1)

if (not new_cfg.isDimsSuperset(old_cfg)):
    print("ERROR: New dimensions are not a superset of the old ones")

if (scale_val!=1.0):
  print("INFO: Scaling additional elements from border by %f"%scale_val)

old_grid_obj =  libcgyrorestart.CGyroGrid()
new_grid_obj =  libcgyrorestart.CGyroGrid()
old_grid_obj.load_from_dict(old_cfg.user_dict)
new_grid_obj.load_from_dict(new_cfg.user_dict)

dim_offset = (new_grid_obj.n_radial-old_grid_obj.n_radial)//2
if (dim_offset>0):
    print("INFO: Shifting radial by %i"%dim_offset)


try:
  add_dims(old_dir, new_dir, old_grid_obj, new_grid_obj, dim_offset, scale_val)
except IOError as err:
    print("IO error: {0}".format(err))
    sys.exit(21)
