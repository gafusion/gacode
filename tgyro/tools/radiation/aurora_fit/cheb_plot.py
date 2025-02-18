import numpy as np
import numpy.polynomial.chebyshev as pc
import matplotlib.pyplot as plt
from matplotlib import rc

rc('font',size=24)
rc('text',usetex=True)

i = 'Kr'
nc = 12

# OLD
x = np.loadtxt('aurora/x.txt')
t = np.loadtxt('aurora/te.txt')
y = np.loadtxt('aurora/'+i+'.txt')
s = np.log(y*10**13)
# 1MW/cm^3 = 10^13 erg/cm^3/s
c = pc.chebfit(x,s,nc)

# New
xn = np.loadtxt('aurora_2025/x.txt')
tn = np.loadtxt('aurora_2025/te.txt')
yn = np.loadtxt('aurora_2025/'+i+'.txt')
sn = np.log(yn*10**13)
# 1MW/cm^3 = 10^13 erg/cm^3/s
cn = pc.chebfit(xn,sn,nc)

fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(1,1,1)
ax.grid(which="both",ls=":")
ax.set_xscale('log')

ax.set_xlabel(r'$T [\mathrm{keV}]$')
ax.set_ylabel(r'$\mathrm{cooling~rate~[erg/cm^3/s]}$')

ax.plot(t/1e3,np.exp(s),color='k',label=i)
ax.plot(t/1e3,np.exp(pc.chebval(x,c)),color='m',label=i+' fit')

ax.plot(tn/1e3,np.exp(sn),color='k',linestyle='--',label=i+' (2025)')
ax.plot(tn/1e3,np.exp(pc.chebval(xn,cn)),color='m',linestyle='--',label=i+' fit (2025)')
ax.legend()
plt.tight_layout(pad=0.3)

#if ofile1 == 'screen':
plt.savefig(i+'.png')
#else:
#   plt.savefig(ofile1)
