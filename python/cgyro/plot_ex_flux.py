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
sim.getexflux()

ns = sim.n_species
nl = sim.n_global+1
t  = sim.t

ky  = sim.ky
ave = np.zeros((sim.n_n,ns))

#self.lky_flux_n =-> (2,ng,self.n_species,self.n_n,nt),'F')

if moment == 'n':
    ntag = 'Density~flux'
    mtag = '\Gamma'
    z = np.sum(sim.lky_flux_n,axis=3)
    ftag = 'ex_flux_n'
elif moment == 'e':
    ntag = 'Energy~flux'
    mtag = 'Q'
    z = np.sum(sim.lky_flux_e,axis=3)
    ftag = 'ex_flux_e'
elif moment == 'v':
    ntag = 'Momentum~flux'
    mtag = '\Pi'
    z = np.sum(sim.lky_flux_v,axis=3)
    ftag = 'ex_flux_v'
else:
    print 'ERROR (plot_ex_flux.py) Invalid moment.'
    sys.exit()

xr = np.zeros((ns,nl))
xi = np.zeros((ns,nl))
for ispec in range(ns):
    for l in range(nl):
        xr[ispec,l] = average(z[0,l,ispec,:],sim.t,w)
        xi[ispec,l] = average(z[1,l,ispec,:],sim.t,w)


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

    print 'INFO: (plot_ex_flux) Created '+fname
    sys.exit()

#============================================================
# Otherwise plot
fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.05,right=0.96,top=0.92,bottom=0.12)
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(r'$r/L_x$')

color = ['k','m','b','c','g','r']

windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+']$'

ax.set_title(r'$\mathrm{'+ntag+'} \quad $'+windowtxt)
    
t = 2*np.pi*np.arange(0.0,1.0,0.001)

ax.set_title(r'$\mathrm{'+ntag+'} \quad $'+windowtxt)

e=0.0
a1 = 0.0
a2 = 2*np.pi

for ispec in range(ns):
    stag = str(ispec)
    u = specmap(sim.mass[ispec],sim.z[ispec])
    # Flux curve
    g = np.zeros(len(t))
    g = xr[ispec,0] 
    for l in range(1,nl):
        g = g+2*(np.cos(l*t)*xr[ispec,l]-np.sin(l*t)*xi[ispec,l])
    ax.plot(t/(2*np.pi),g)
    # Flux integral
    gi = xr[ispec,0]
    for l in range(1,nl):
        gi = gi+2*((np.sin(l*a2)-np.sin(l*a1))*xr[ispec,l]-(np.cos(l*a2)-np.cos(l*a1))*xi[ispec,l])/l/(a2-a1)
    label = r'$'+mtag+'_'+u+'/'+mtag+'_\mathrm{GB}: '+str(round(gi,3))+'$'
    ax.plot([e,1-e],[gi,gi],label=label)
    ax.axvspan(0,e,facecolor='g',alpha=0.1)
    ax.axvspan(1-e,1,facecolor='g',alpha=0.1)

if ymax != 'auto':
    ax.set_ylim([0,float(ymax)])

ax.legend(loc=2)

if ftype == 'screen':
   plt.show()
else:
    fname=fname+'.'+ftype
    print 'INFO: (plot_ex_flux) Created '+fname
    plt.savefig(fname)
