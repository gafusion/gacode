# file processed by 2to3
import sys
import numpy as np
from gacodeplotdefs import *

ftype = sys.argv[1]
ie    = int(sys.argv[2])
n0    = int(sys.argv[3])
dir   = sys.argv[4]

#-------------------------------------------------------
# Read grid dimension and axes
#
data = np.loadtxt('out.neo.grid_3d')

n_species = int(data[0])
n_energy  = int(data[1])+1
nxi       = 100
msize     = int(data[3])
#-------------------------------------------------------

#-------------------------------------------------------
# (m,n) indices
#
indx = np.loadtxt('out.le3.geoindx',dtype=int)
#-------------------------------------------------------

#-------------------------------------------------------
# Read H
#
data = np.loadtxt('out.neo.gxi_3d')
fx = np.reshape(data,(nxi,msize,n_energy,n_species),'F')
fx = fx/np.max(fx)
#-------------------------------------------------------

#-------------------------------------------------------
# Read xi grid
#
xi = np.loadtxt('out.neo.gxi_3d_x')
#-------------------------------------------------------

fig = plt.figure(figsize=(14,14))
fig.subplots_adjust(left=0.07,right=0.95,top=0.94,bottom=0.06,wspace=0.25,hspace=0.32)
fig.suptitle(r'$\mathbf{dir}$')

nm = int(max(indx[:,1]))+1
nn = int(max(indx[:,2]))+1

if nm <= 4:
    k=2
elif nm <= 9:
    k=3
elif nm <= 16:
    k=4

colors = ['black','red','magenta','blue']

init = np.zeros([nm],dtype=int)

for p in range(msize):
    i = indx[p,0]
    m = indx[p,1]
    n = indx[p,2]
    if n == n0:
        if init[m] == 0:
            ax = fig.add_subplot(k,k,m+1)
            ax.grid(which="both",ls=":")
            ax.grid(which="major",ls=":")
            ax.set_title(r'$(m,n)=('+str(m)+','+str(n)+')$')
            ax.set_xlabel(r'$\xi = v_\parallel/v$')
            init[m] = 1

        fp = np.array(fx[:,p,ie,0])
        ax.plot(xi,fp,'-',color=colors[i-1],markersize=2)

if ftype == 'screen':
    plt.show()
else:
    outfile = 'h.'+ftype
    plt.savefig(outfile)
