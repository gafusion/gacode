#!/usr/bin/env python
# Usage:
#
#   python ufile2.py <datafile> <time>
#
#   NOTES:
#   - DOC: http://tokamak-profiledb.ccfe.ac.uk/DOCS/PR08MAN/pdbman.html

import sys
import numpy as np

def extract2d(infile):

    vt = np.zeros(1)
    data_region = 0
    for line in open(infile,'r').readlines():
        if line.count("END-OF-DATA") == 1:
            data_region = 0
            # Compute average if time-point in range
            fxt  = vy.reshape((nx,nt),order='F')
            print('Converting '+var)
            # Write the averaged data for current profile (var)
            # Output filename: "out.TAG.ave"
            np.savetxt(var+'.txt',np.transpose((vx,fxt[:,0])),fmt='%1.6e')

        if data_region == 0:

            if line.count("-DEP") == 1:
                # Extract current variable name
                var=line.split('  ')[0].strip()

            if line.count("X PTS") == 1:
                # Get length of radial grid
                nx=int(line.split(";-# OF X")[0].strip())
                vx = np.zeros(nx)

            if line.count("Y PTS") == 1:
                # Get size of time grid
                nt=int(line.split(";-# OF Y")[0].strip())
                vy = np.zeros(nx*nt)
                vt = np.zeros(nt)
                # Signal header end
                data_region = 1
                ix = 0
                it = 0
                iy = 0
                # Write dimensions
                #if var=="NE":
                f=open('out.dim','w')
                f.write(str(nx)+'\n')
                f.close()

        else:

            # Here we are reading numbers

            # Number of columns
            n = np.max([line.count("E"),line.count("e")])

            if ix < nx:
                # Read the radial grid [nx points]
                for i in range(n):
                    vx[ix]=float(line[1+13*i:1+13*i+13])
                    ix = ix+1
            elif it < nt:
                # Read the time grid [nt points]
                vt[0]=float(line)
                it = 1
            else:
                # Read the radial-time grid [(nx+1)*ny points]
                for i in range(n):
                    vy[iy]=float(line[1+13*i:1+13*i+13])
                    iy = iy+1

#----------------------------------------------------------------------------
# Manage input parameters

# <datafile>
try:
    infile=sys.argv[1]
except:
    print('Usage: python ufile2.py <datafile> <time>')
    sys.exit()

varlist = []
t0 = -1.0
# Generate list of included profiles (tags)
for line in open(infile,'r').readlines():
    if line.count("-DEP") == 1:
        # Extract current variable name
        var=line.split('  ')[0].strip()
        varlist.append(var)
# Print list of included profiles (tags) in blocks of 10
print('INFO: (ufile2) tags ->')
for i in np.arange(start=0,stop=len(varlist),step=10):
    print('       '+' '.join(varlist[i:i+10]))

extract2d(infile)
