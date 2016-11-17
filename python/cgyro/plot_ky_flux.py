import sys
import numpy as np
from gacodeplotdefs import *
from gacodefuncs import *
from cgyro.data import cgyrodata

ftype  = sys.argv[1]
w = float(sys.argv[2])
moment = sys.argv[3]
ymin   = sys.argv[4]
ymax   = sys.argv[5]

sim = cgyrodata('./')
sim.getbigflux()

ns = sim.n_species
t  = sim.t

ky  = sim.ky
ave = np.zeros((sim.n_n,ns))

if moment == 'n':
    ntag = 'Density~flux'
    mtag = '\Gamma'
    y = np.sum(sim.kxky_flux_n,axis=0)
    ftag = 'ky_flux_n'
elif moment == 'e':
    ntag = 'Energy~flux'
    mtag = 'Q'
    y = np.sum(sim.kxky_flux_e,axis=0)
    ftag = 'ky_flux_e'
else:
    print 'ERROR (plot_ky_flux.py) Invalid moment.'
    sys.exit()

# Determine tmin
imin=iwindow(t,w)

fname = 'out.cgyro.'+ftag

#============================================================
# Dump and exit if desired:
if ftype == 'dump':
    arr = np.zeros([len(ky),ns+1])
    arr[:,0] = ky
    stag = '# (k_y rho_s'
    for ispec in range(ns):
        for j in range(sim.n_n):
            ave = average(y[ispec,j,:],sim.t,w)
            arr[j,ispec+1] = ave
        stag = stag+' , s'+str(ispec)
            
    fid = open(fname,'w')
    fid.write('# Moment  : '+mtag+'\n')
    fid.write('# Time    : '+str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])+'\n')
    fid.write(stag+')\n')
    np.savetxt(fid,arr,fmt='%.5e')
    fid.close()

    print 'INFO: (plot_ky_flux) Created '+fname
    sys.exit()

#============================================================
# Otherwise plot
fig = plt.figure(figsize=(6*ns,6))
fig.subplots_adjust(left=0.06,right=0.96,top=0.92,bottom=0.12)

color = ['m','k','b','c']

windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+']$'

if ky[1] < 0.0:
    ky = -ky
    xlabel=r'$-k_\theta \rho_s$'
else:
    xlabel=r'$k_\theta \rho_s$'
    
dk = ky[1]-ky[0]
    
for ispec in range(ns):
    for j in range(sim.n_n):
        ave[j,ispec] = average(y[ispec,j,:],sim.t,w)

for ispec in range(ns):
    stag = str(ispec)
    ax = fig.add_subplot(1,ns,ispec+1)
    ax.set_xlabel(xlabel)
    u = specmap(sim.mass[ispec],sim.z[ispec])
    ax.set_ylabel(r'$'+mtag+'_'+u+'$',color='k')
    ax.set_title(windowtxt)
    ax.bar(ky-dk/2.0,ave[:,ispec],width=dk/1.1,color=color[ispec],alpha=0.5,edgecolor='black')
    ax.set_xlim([0,ky[-1]+dk])
    if ymax != 'auto':
        ax.set_ylim([0,float(ymax)])
    # Dissipation curve             
    ax.plot(ky,sim.alphadiss*ax.get_ylim()[1]*0.5,linewidth=2,color='k',alpha=0.2)

if ftype == 'screen':
   plt.show()
else:
    fname=fname+'.'+ftype
    print 'INFO: (plot_ky_flux) Created '+fname
    plt.savefig(fname)
