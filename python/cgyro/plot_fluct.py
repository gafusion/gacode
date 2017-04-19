import sys
import numpy as np
import os.path

from gacodeplotdefs import *
from gacodefuncs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
moment = sys.argv[2]
species = int(sys.argv[3])
ymin = sys.argv[4]
ymax = sys.argv[5]
nx = int(sys.argv[6])
ny = int(sys.argv[7])
istr = sys.argv[8]
fmin = sys.argv[9]
fmax = sys.argv[10]
colormap = sys.argv[11]

sim = cgyrodata('./')
nt = sim.n_time
nr = sim.n_radial
nn = sim.n_n
ns = sim.n_species

epx = np.zeros([nx,nr],dtype=np.complex)
eny = np.zeros([ny,nn],dtype=np.complex)
x = np.zeros([nx])
y = np.zeros([ny])

#------------------------------------------------------------------------
# Fourier arrays
#
for i in range(nx):
    x[i] = i*2*np.pi/(nx-1)
    for p in range(nr):    
        epx[i,p]=np.exp(1j*(p-nr/2)*x[i])

for j in range(ny):
    y[j] = j*2*np.pi/(ny-1)
    for n in range(nn):    
        eny[j,n]=np.exp(-1j*n*y[j])

# factor of 1/2 for n=0
eny[:,0] = 0.5*eny[:,0]
#------------------------------------------------------------------------

#------------------------------------------------------------------------
# Real-space field resonstruction
#
def maptoreal(nr,nn,nx,ny,c):

    import numpy as np
    import time

    start = time.time()

    # This needs to be fast, so we use numpy.outer
    f = np.zeros([nx,ny])
    for p in range(nr):
        for n in range(nn):
            f[:,:] = f[:,:]+np.real(c[p,n]*np.outer(epx[:,p],eny[:,n]))
    
    end = time.time()
  
    return f,str(end-start)
#-----------------------------------------------------------------------

if istr == '-1':
    ivec = range(nt)
else:
    ivec = str2list(istr)

# Set filename and title
if (moment == 'n'):
    fdata = 'out.cgyro.kxky_n'
    title = r'$\delta \mathrm{n}$'
elif (moment == 'e'):
    fdata = 'out.cgyro.kxky_e'
    title = r'$\delta \mathrm{E}$'
elif (moment == 'phi'):
    fdata = 'out.cgyro.kxky_phi'
    title = r'$\delta\phi$'

# Check to see if it exists
if not os.path.isfile(fdata):
    print fdata+' does not exist.  Try -moment phi'
    sys.exit()

#(2,n_radial,n_species,n_n,nt)
if (moment == 'phi'):
    n_chunk = 2*nr*nn
else:
    n_chunk = 2*nr*ns*nn
    
m = 0
i = 0

aa = np.zeros([n_chunk])
for line in open(fdata):
    aa[m] = float(line)
    m = m+1
    if m == n_chunk:
        i = i+1
        m = 0
        print 'INFO: (plot_fluct) Time index '+str(i) 
        if i in ivec:
            if (moment == 'phi'):
                a = np.reshape(aa,(2,nr,nn),'F')
                c = a[0,:,:]+1j*a[1,:,:]
            else:
                a = np.reshape(aa,(2,nr,ns,nn),'F')
                c = a[0,:,species,:]+1j*a[1,:,species,:]
                
            f,t = maptoreal(nr,nn,nx,ny,c)
            if fmin == 'auto':
                f0=np.min(f)
                f1=np.max(f)
            else:
                f0=float(fmin)
                f1=float(fmax)

            fig = plt.figure(figsize=(10,10))
            fig.subplots_adjust(left=0.08,right=0.96,top=0.94,bottom=0.08)
            ax = fig.add_subplot(111)
            ax.set_title(title,fontsize=16)
            ax.set_xlabel(r'$x/\rho_s$',fontsize=16)
            ax.set_ylabel(r'$y/\rho_s$',fontsize=16)
        
            levels = np.arange(f0,f1,(f1-f0)/256)
            ax.contourf(x/2/np.pi*sim.length,y/sim.ky[1],np.transpose(f),levels,cmap=plt.get_cmap(colormap))
            print 'INFO: (plot_fluct) min,max = ',f0,f1

            if ftype == 'screen':
                plt.show()
            else:
                fname = fdata+str(i)
                # Filename uses number padded with zeros
                plt.savefig(str(i)+'.png')
                # Close each time to prevent memory accumulation
                plt.close()

            if i == max(ivec):
                sys.exit()





