import sys
import numpy as np
from gacodeplotdefs import *
from gacodefuncs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
w = float(sys.argv[2])
ymin  = sys.argv[3]
ymax  = sys.argv[4]
nstr  = sys.argv[5]

sim = cgyrodata('./')
sim.getbigfield()

t   = sim.t
kx  = sim.kx
ave = np.zeros(sim.n_radial)

imin=iwindow(sim.t,w)
    
dk = kx[1]-kx[0]
x0 = kx[-1]+dk

#============================================================
# Set axes and plot
fig = plt.figure(figsize=(12,6))
fig.subplots_adjust(left=0.08,right=0.96,top=0.92,bottom=0.12)
ax = fig.add_subplot(1,1,1)

color = ['m','k','b','c']
xlabel=r'$k_x \rho_s$'
windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+']$'

ax.set_title(r'$\mathrm{Average~fluctuation~intensity} \quad $'+windowtxt)
ax.set_xlabel(xlabel)

if nstr == 'null':
    y = np.sum(sim.phisq,axis=1)
    for j in range(sim.n_radial):
        ave[j] = average(y[j,:],sim.t,w)
    ax.set_ylabel(r'$\overline{\delta \phi_\mathrm{total}}$',color='k')
    ax.plot(kx,np.sqrt(ave[:]),color=color[0],ls='steps')
else:
    y = np.zeros([sim.n_radial,sim.n_time])
    nvec = str2list(nstr)
    print 'INFO: (plot_kx_phi) n = '+str(nvec)
    ax.set_ylabel(r'$\overline{\Phi_n}$',color='k')
    for n in nvec:
        num = r'$n='+str(n)+'$'
        y[:] = sim.phisq[:,n,:]
        for j in range(sim.n_radial):
            ave[j] = average(sim.phisq[j,n,:],sim.t,w)
        ax.plot(kx+dk/2,np.sqrt(ave[:]),ls='steps',label=num)
            
ax.set_xlim([-x0,x0])
ax.set_yscale('log')

ymin,ymax=setlimits(ax.get_ylim(),ymin,ymax)
ax.set_ylim([ymin,ymax])

# Dissipation curve             
ax.plot(kx,sim.radialdiss*ax.get_ylim()[1]*0.5,linewidth=2,color='k',alpha=0.2)

if sim.n_n > 16:
    ax.legend(loc=4, ncol=5, prop={'size':12})
else:
    ax.legend(loc=4, ncol=6, prop={'size':12})

fname = 'out.cgyro.kx_phi'

if ftype == 'screen':
   plt.show()
else:
    fname=fname+'.'+ftype
    print 'INFO: (plot_kx_phi) Created '+fname
    plt.savefig(fname)
