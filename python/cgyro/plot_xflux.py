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


# NOTE: lky_flux_* -> [ 2, nl , ns , n_n , nt ]
#                       0  1    2     3    4 

if moment == 'n':
    ntag = 'Density~flux'
    mtag = '\Gamma'
    z = np.sum(sim.lky_flux_n,axis=3)
    ftag = 'xflux_n'
elif moment == 'e':
    ntag = 'Energy~flux'
    mtag = 'Q'
    z = np.sum(sim.lky_flux_e,axis=3)
    ftag = 'xflux_e'
elif moment == 'v':
    ntag = 'Momentum~flux'
    mtag = '\Pi'
    z = np.sum(sim.lky_flux_v,axis=3)
    ftag = 'xflux_v'
else:
    print 'ERROR (plot_xflux.py) Invalid moment.'
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
    fid = open(fname,'w')
    fid.write('# Moment  : '+mtag+'\n')
    fid.write('# Time    : '+str(sim.t[imin])+' < (c_s/a) t < '+str(sim.t[-1])+'\n')
    #fid.write(stag+')\n')
    np.savetxt(fid,xr[0,:],fmt='%.5e')
    np.savetxt(fid,xi[0,:],fmt='%.5e')
    fid.close()
    print 'WARNING: (plot_xflux) Data dump not fully implemented in '+fname
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

for ispec in range(ns):
    # Flux curve
    g = np.zeros(len(t))
    g = xr[ispec,0] 
    for l in range(1,nl):
        g = g+2*(np.cos(l*t)*xr[ispec,l]-np.sin(l*t)*xi[ispec,l])
    # Flux average
    g0 = xr[ispec,0]
    u     = specmap(sim.mass[ispec],sim.z[ispec])
    label = r'$'+mtag+'_'+u+'/'+mtag+'_\mathrm{GB}: '+str(round(g0,3))+'$'
    # Plot
    ax.plot(t/(2*np.pi),g)
    ax.plot([0,1],[g0,g0],label=label)

if ymax != 'auto':
    ax.set_ylim([float(ymin),float(ymax)])

ax.axvspan(0.25,0.75,facecolor='g',alpha=0.1)

ax.legend(loc=2)

if ftype == 'screen':
   plt.show()
else:
    fname=fname+'.'+ftype
    print 'INFO: (plot_xflux) Created '+fname
    plt.savefig(fname)
