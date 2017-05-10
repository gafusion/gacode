import sys
import numpy as np
from gacodeplotdefs import *
from mpl_toolkits.mplot3d import Axes3D
from cgyro.data import cgyrodata

ftype = sys.argv[1]
w     = float(sys.argv[2])
field = sys.argv[3]

sim = cgyrodata('./')
sim.getbigfield()

#-----------------------------------------------------------------
# Note array structure
# self.phi = np.reshape(data,(2,self.n_radial,self.n_n,nt),'F')

t = sim.t
nx=sim.n_radial
ny=sim.n_n

f = np.zeros([nx-1,ny])
n = sim.n_time

itheta=0

imin = int((1.0-w)*n)
for i in np.arange(imin,n):
    f = f+sim.phisq[1:,itheta,:,i]

# Fix (0,0)
i0 = nx/2-1
f[i0,0] = 1e-6

# Reverse y order for image plotting
f = f[:,::-1]

# Scale data
f = np.log(f)

x0 = max(abs(sim.kx))*0.25
y0 = max(abs(sim.ky))

asp=y0/(2*x0)

fig = plt.figure(figsize=(15,17*asp))
fig.subplots_adjust(left=0.05,right=0.96,top=0.91,bottom=0.14)
ax = fig.add_subplot(111)

windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+']$'

ax.set_xlabel(r'$k_x \rho_s/4$',fontsize=GFONTSIZE)
ax.set_ylabel(r'$k_y \rho_s$',fontsize=GFONTSIZE)
ax.set_title(r'$\mathrm{Average~fluctuation~intensity} \quad $'+windowtxt)

ax.imshow(np.transpose(f),extent=[-x0,x0,0,y0],interpolation='none')

if ftype == 'screen':
    plt.show()
else:
    outfile = 'kxky_phisq.'+ftype
    plt.savefig(outfile)
