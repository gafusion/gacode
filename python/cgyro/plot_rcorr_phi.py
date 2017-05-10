import sys
import numpy as np
import scipy as scipy
import scipy.signal as signal
from scipy.optimize import curve_fit

def absexp(x, tau):
    return np.exp(- np.abs(x)/tau)

from gacodeplotdefs import *
from gacodefuncs import *
from cgyro.data import cgyrodata

ftype = sys.argv[1]
w = float(sys.argv[2])
ymin  = sys.argv[3]
ymax  = sys.argv[4]
#nstr  = sys.argv[5]

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
xlabel=r'$r / \rho_s$'
windowtxt = r'$['+str(t[imin])+' < (c_s/a) t < '+str(t[-1])+']$'

ax.set_title(r'$\mathrm{Average~radial~correlation} \quad $'+windowtxt)
ax.set_xlabel(xlabel)

itheta=0

#y = np.sum(sim.phisq,axis=1)
y = np.sum(sim.phisq[:,itheta,1:,:],axis=1)
for j in range(sim.n_radial):
  ave[j] = average(y[j,:],sim.t,w)
#corr = np.zeros(sim.n_radial)
ave = np.roll(ave,-sim.n_radial/2)
ave[0] = 0.
corr = np.fft.fft(ave,sim.n_radial)
corr = np.fft.fftshift(corr)
corr /= np.max(np.abs(corr))
corr = corr.real
delta_r = np.fft.fftfreq(sim.n_radial)
delta_r = np.fft.fftshift(delta_r)
Lx = 2*np.pi/dk
delta_r *= Lx

#calculate envleope
corr_hilbert = signal.hilbert(corr)
corr_env = np.abs(corr_hilbert)
ax.set_ylabel(r'$C_{\delta \phi}(\Delta r)$',color='k')
ax.plot(delta_r, 0*delta_r, color='k', ls='--')
ax.plot(delta_r, corr, color=color[0])

#ax.plot(delta_r, corr_env, color=color[0],ls=':')
l_corr, pcov = curve_fit(absexp, delta_r, corr_env, p0=10.)
ax.plot(delta_r, absexp(delta_r,l_corr), color=color[1],ls='-.')
print 'l_corr = ', l_corr

ax.set_xlim([np.min(delta_r),np.max(delta_r)])
ax.set_ylim(-1,1)

#if sim.n_n > 16:
#    ax.legend(loc=4, ncol=5, prop={'size':12})
#else:
#    ax.legend(loc=4, ncol=6, prop={'size':12})

fname = 'out.cgyro.rcorr_phi'

if ftype == 'screen':
   plt.show()
else:
    fname=fname+'.'+ftype
    print 'INFO: (plot_rcorr_phi) Created '+fname
    plt.savefig(fname)
