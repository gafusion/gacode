import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from profiles_gen.data import profiles_genData

rc('text',usetex=True)
rc('font',size=18)

# Number of theta-points for plotting
narc = 256

infiles = sys.argv[1]
surf    = sys.argv[2]
i0      = int(sys.argv[3])
ftype   = sys.argv[4]

filevec = string.splitfields(infiles,',')

# Support only one input.profiles file
prof = profiles_genData(filevec[0])

fig = plt.figure(figsize=(8,12))
ax = fig.add_subplot(111,aspect='equal')
ax.set_xlabel(r'$R$')
ax.set_ylabel(r'$Z$')

t = 2*np.pi*np.arange(0,narc)/float(narc-1)

if i0 < 0:
    i1=0
    i2=prof.n_exp
else:
    i1=i0
    i2=i1+1

if surf == 'miller' or surf == 'both':
    # Miller geometry flux-surfaces
    for i in np.arange(i1,i2):
        bigr = prof.data['rmaj'][i]  
        bigz = prof.data['zmag'][i]  
        r    = prof.data['rmin'][i]
        d    = prof.data['delta'][i]
        k    = prof.data['kappa'][i]
        z    = prof.data['zeta'][i]

        x = bigr+r*np.cos(t+np.arcsin(d)*np.sin(t))
        y = bigz+k*r*np.sin(t+z*np.sin(2*t))

        if i == i1:
            ax.plot(x,y,'-k',linewidth=1,label=r'$\mathrm{Miller}$')
        else:
            ax.plot(x,y,'-k',linewidth=1)


if surf == 'fourier' or surf == 'both':
    # Fourier geometry flux-surfaces 
    for i in np.arange(i1,i2):
        a0r = prof.geo['ar'][0,i]
        a0z = prof.geo['az'][0,i]

        x = a0r/2
        y = a0z/2

        for j in range(1,1+prof.nfourier,1):
            ar = prof.geo['ar'][j,i]
            br = prof.geo['br'][j,i]
            az = prof.geo['az'][j,i]
            bz = prof.geo['bz'][j,i]
            x = x+ar*np.cos(j*t)+br*np.sin(j*t)
            y = y+az*np.cos(j*t)+bz*np.sin(j*t)

        if i == i1:
            ax.plot(x,y,'-b',linewidth=1,label=r'$\mathrm{Fourier}~'+str(prof.nfourier)+'$')
        else:
            ax.plot(x,y,'-b',linewidth=1)


# LCFS
i = prof.n_exp-1
a0r = prof.geo['ar'][0,i]
a0z = prof.geo['az'][0,i]

x = a0r/2
y = a0z/2

for j in range(1,1+prof.nfourier,1):
    ar = prof.geo['ar'][j,i]
    br = prof.geo['br'][j,i]
    az = prof.geo['az'][j,i]
    bz = prof.geo['bz'][j,i]
    x = x+ar*np.cos(j*t)+br*np.sin(j*t)
    y = y+az*np.cos(j*t)+bz*np.sin(j*t)

ax.plot(x,y,'-m',linewidth=2)

ax.legend()
if ftype == 'screen':
    plt.show()
else:
    outfile = ftype
    plt.savefig(outfile)
