import sys
import numpy as np
from gacodeplotdefs import *
from gacodefuncs import *
from cgyro.data import cgyrodata

ftype  = sys.argv[1]
w      = float(sys.argv[2])
moment = sys.argv[3]
ymin   = sys.argv[4]
ymax   = sys.argv[5]
field  = int(sys.argv[6])
fc     = int(sys.argv[7])

sim = cgyrodata('./')
sim.getflux()

ns = sim.n_species
t  = sim.t

field_tag = '\mathrm{Total}'

if hasattr(sim,'ky_flux'):
    # New flux format
    if fc == 0:
        ys = np.sum(sim.ky_flux,axis=(2,3))
    else:
        ys = np.sum(sim.ky_flux[:,:,field,:,:],axis=2)
        if field == 0:
            field_tag = '\phi'
        elif field == 1:
            field_tag = 'A_\parallel'
        else:
            field_tag = 'B_\parallel'
else:
    # Old/large flux format
    sim.getbigflux()          
    ys = np.zeros([ns,3,sim.n_time])
    if moment == 'n':
        ys[:,0,:] = np.sum(sim.kxky_flux_n,axis=(0,2))
    if moment == 'e':
        ys[:,1,:] = np.sum(sim.kxky_flux_e,axis=(0,2))

# Now, ys -> {n_species,3,nt}

if moment == 'n':
    ntag = 'Density~flux'
    mtag = '\Gamma'
    ttag = 'G'
    ftag = 'flux_n'
    y = ys[:,0,:]
elif moment == 'e':
    ntag = 'Energy~flux'
    mtag = 'Q'
    ttag = 'Q'
    ftag = 'flux_e'
    y = ys[:,1,:]
elif moment == 'v':
    ntag = 'Momentum~flux'
    mtag = '\Pi'
    ttag = 'Pi'
    ftag = 'flux_v'
    y = ys[:,2,:]
else:
    print 'ERROR (plot_flux.py) Invalid moment.'
    sys.exit()

# Get index for average window
imin=iwindow(t,w)

fname = 'out.cgyro.'+ftag

#============================================================
# Dump and exit if desired:
if ftype == 'dump':
    data = np.column_stack((sim.t,y[0,:]))
    head = '(cs/a) t     '+ttag+'_1/'+ttag+'_GB'
    for ispec in range(1,ns,1):
        head = head+'         '+ttag+'_'+str(ispec+1)+'/'+ttag+'_GB'
        data = np.column_stack((data,y[ispec,:]))
    np.savetxt(fname,data,fmt='%.8e',header=head)
    print 'INFO: (plot_flux) Created '+fname
    sys.exit()
#============================================================

#============================================================
# Otherwise plot
fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.05,right=0.96,top=0.92,bottom=0.12)
ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(TIME)

color = ['k','m','b','c','g','r']

windowtxt = '['+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+']'

ax.set_title(r'$\mathrm{'+ntag+'} \quad '+windowtxt+'\quad ['+field_tag+']$')

for ispec in range(ns):
    ave = average(y[ispec,:],t,w)
    y_ave = ave*np.ones(len(t))
    u = specmap(sim.mass[ispec],sim.z[ispec])
    label = r'$'+mtag+'_'+u+'/'+mtag+'_\mathrm{GB}: '+str(round(ave,3))+'$'
    # Average
    ax.plot(t[imin:],y_ave[imin:],'--',color=color[ispec])
    # Time trace
    ax.plot(sim.t,y[ispec,:],label=label,color=color[ispec])

ax.legend(loc=2)

if ymax != 'auto':
    ax.set_ylim([0,float(ymax)])

if ftype == 'screen':
   plt.show()
else:
   fname = fname+'.'+ftype
   print 'INFO: (plot_flux) Created '+fname
   plt.savefig(fname)
