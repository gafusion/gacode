from gacodeinput import *
import sys,os

# Usage:
#  xgyro_parse.py n_mpi

total_MPI = int(sys.argv[1])

# Read the fixed parameters first
base = SimpleInput()
base.add('MPI_RANK_ORDER','1')
base.add('N_DIRS','0')
# Perform the parsing
base.read_input('input.xgyro')

n_dirs = int(base.data_dict['N_DIRS'])
if (n_dirs<1):
  print("ERROR: At least one CGYRO dir needed in input.xgyro.", file=sys.stderr)
  sys.exit(1)

del base

# now do the full parsing
x = SimpleInput()

# XGYRO input parameters
x.add('MPI_RANK_ORDER','1')
x.add('N_DIRS',"%s"%n_dirs)
for i in range(n_dirs):
  fi = i+1
  x.add('MIN_MPI_%i'%fi, '1') 
  x.add('DIR_%i'%fi, 'cgyro_dir_%i'%fi)

# Perform the parsing
x.read_input('input.xgyro',write_ext=False)
if (x.error!=0):
  x.printmsg()
  sys.exit(x.error)

# convert MIN_MPI into N_MPI, to get the same sum as total_MPI
min_mpi=0
for i in range(n_dirs):
  fi = i+1
  min_mpi += int(x.data_dict['MIN_MPI_%i'%fi])

if min_mpi>total_MPI:
  print("ERROR: Need at least %i MPI processes."%min_mpi, file=sys.stderr)
  sys.exit(1)

if (total_MPI%min_mpi)!=0:
  print("ERROR: MPI rank must be multiple of %i."%min_mpi, file=sys.stderr)
  sys.exit(1)

mpi_mult = total_MPI//min_mpi

for i in range(n_dirs):
  fi = i+1
  old_k = 'MIN_MPI_%i'%fi
  new_k = 'N_MPI_%i'%fi
  x.rename(old_k,new_k)
  x.data_dict[new_k] = '%i'%(int(x.data_dict[new_k])*mpi_mult)

# now deal with cgyro setup
for i in range(n_dirs):
  fi = i+1
  cgyro_dir = x.data_dict['DIR_%i'%fi]
  if not os.path.isfile(cgyro_dir+"/input.cgyro"):
    print("ERROR: Cound not find %s/input.cgyro"%cgyro_dir, file=sys.stderr)
    sys.exit(1)

# now we can write out the gen file
x.write_parsed('input.xgyro.gen')

# now that everything works, print out the CGYRO subdirs for further processing
for i in range(n_dirs):
  fi = i+1
  cgyro_dir = x.data_dict['DIR_%i'%fi]
  print(cgyro_dir)

