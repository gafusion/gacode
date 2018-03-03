import os
import sys
import numpy as np

if len(sys.argv) == 1:
    print 'python cgyro_b2b.py <float64> <float32>'
    sys.exit()
    
afile = sys.argv[1]
bfile = sys.argv[2]

if not os.path.isfile(afile):
    print "INFO: (cgyro_b2b) "+afile+' not found'
    sys.exit()
    
arr = np.fromfile(afile,dtype='float64')
brr = arr.astype('float32')

print "INFO: (cgyro_b2b) Read float64 data in "+afile

brr.tofile(bfile)

print "INFO: (cgyro_t2b) Wrote float32 data to "+bfile
