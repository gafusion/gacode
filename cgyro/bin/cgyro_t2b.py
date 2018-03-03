import os
import sys
import numpy as np

if len(sys.argv) == 1:
    print 'python cgyro_t2b.py <text-source> <binary-dest>'
    sys.exit()
    
afile = sys.argv[1]
bfile = sys.argv[2]

if not os.path.isfile(afile):
    print "INFO: (cgyro_t2b) "+afile+' not found'
    sys.exit()
    
arr = np.fromfile(afile,dtype='float32',sep=' ')

print "INFO: (cgyro_t2b) Read ASCII data in "+afile

arr.tofile(bfile)

print "INFO: (cgyro_t2b) Wrote binary data to "+bfile
