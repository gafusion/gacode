# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import sys
import numpy as np
from gacodeplotdefs import *

infile = 'out.profiles_3d.geoall'
ftype  = 'screen'

data=np.loadtxt(infile)
n_t = int(data[0])+1
n_p = int(data[1])+1
n = len(data)-2
vec = data[2:len(data)]

print('Sanity check: ',n-n_t*n_p*4*4)

vec = np.reshape(vec,(4,4,n_p,n_t),order='F')

d = np.zeros([n_p,n_t])
vp = np.arange(n_p)
vt = np.arange(n_t)

for i in range(2):
    for j in range(4):
        d = d+vec[i,j,:,:]**2

d = np.sqrt(d)

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.set_xlabel(r'$n$')
ax.set_ylabel(r'$\delta$')
ax.set_yscale('log')

for i in vt:
    ax.plot(vp[0:],d[0:,i],'o',markersize=8,label='m='+str(i))
    ax.plot(vp[0:],d[0:,i],'k',linewidth=1,alpha=0.4)


ax.legend()
if ftype == 'screen':
    plt.show()
else:
    outfile = ftype
    plt.savefig(outfile)
