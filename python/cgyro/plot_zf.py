import sys
import numpy as np
from gacodeplotdefs import *
from gacodefuncs import *
from cgyro.data import cgyrodata

ftype  = sys.argv[1]
w      = float(sys.argv[2])
ifield = int(sys.argv[3])

sim = cgyrodata('./')

if sim.n_n > 1:
    print "ERROR: (plot_zf.py) This plot option valid for ZF test only."
    sys.exit()
elif ifield > 0:
    print "ERROR: (plot_zf.py) Only Phi can be plotted."
else:
    sim.getbigfield()

t    = sim.t
phic = sim.kxky_phi[0,0,0,:]
y    = phic[:]*1e6*(1-np.i0(sim.kx[0]**2)*np.exp(-sim.kx[0]**2))/(np.i0(sim.kx[0]**2)*np.exp(-sim.kx[0]**2))
#initialization in code is with 1e-6*besselj0 # phic[0]

#----------------------------------------------------
# Average calculations
imin=iwindow(t,w)

ave = average(y[:],t,w)
print 'INFO: (plot_zf) Integral time-average = %.6f' % ave

ave_vec = ave*np.ones(len(t))
#----------------------------------------------------

#=====================================================================
fig = plt.figure(figsize=(10,6))
fig.subplots_adjust(left=0.1,right=0.95,top=0.95,bottom=0.12)

ax = fig.add_subplot(111)
ax.grid(which="majorminor",ls=":")
ax.grid(which="major",ls=":")
ax.set_xlabel(TIME)
ax.set_ylabel(r'$\mathrm{Re}\Phi/\Phi_0$')

for i in range(sim.n_radial):
    ax.plot(t,sim.kxky_phi[0,i,0,:]*1e6*(1-np.i0(sim.kx[i]**2)*np.exp(-sim.kx[i]**2))/(np.i0(sim.kx[i]**2)*np.exp(-sim.kx[i]**2)),label=r'$k_x=%.3g$' % sim.kx[i])


    
ax.plot(t[imin:],ave_vec[imin:],color='b',
        label=r'$\mathrm{Average}$',linewidth=1)

theory = 1./(1.0+1.6*sim.q**2/np.sqrt(sim.rmin/sim.rmaj))
ax.plot([0,max(t)],[theory,theory],color='grey',
        label=r'$\mathrm{RH \; theory}$',alpha=0.3,linewidth=4)

theory2 = 1./(1.0+2.*sim.q**2)
ax.plot([0,max(t)],[theory2,theory2],color='m',
        label=r'$\mathrm{fluid \; theory}$',alpha=0.3,linewidth=4)

ax.legend(loc=1,prop={'size':10})


if ftype == 'screen':
    plt.show()
else:
    outfile = 'out.cgyro.zf.'+str(ifield)+'.'+ftype
    plt.savefig(outfile)

