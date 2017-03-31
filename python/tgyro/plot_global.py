import sys
import string
import numpy as np
from gacodeplotdefs import *
from tgyro.data import tgyrodata
from profiles_gen.data import profiles_genData

simdir = sys.argv[1]
f    = sys.argv[2]
ext  = sys.argv[3]
rvar = sys.argv[4]
loc  = int(sys.argv[5])

fig = plt.figure(figsize=(14,8))
fig.subplots_adjust(left=0.07,right=0.96,top=0.94,bottom=0.09)
ax  = fig.add_subplot(111)

def plotcurve(sim,prof):
    # Look for f in keys
    keys  = sorted(prof.data.keys())
    for i in range(len(keys)):
        if f == keys[i].split()[0]:
            success = 1
            fulltag = keys[i]
    if success == 0:
        print "ERROR: (plot_global.py) Bad profile : "+f
        sys.exit()

    ax.plot(x,prof.data[fulltag],color=color,alpha=alpha,linewidth=width,
            label=r'$\mathbf{'+label+'}$')

sim  = tgyrodata(simdir)
prof = profiles_genData('input.profiles')

if rvar == "r":
    ax.set_xlabel(r"$r~\mathrm{(m)}$",fontsize=GFONTSIZE)
    x = prof.data['rmin']
    r_star = sim.data['r_*/a'][0]*max(x)
    r_top  = sim.data['r_top/a'][0]*max(x)
    ax.axvspan(r_star,r_top,facecolor='purple',alpha=0.1)

if rvar == "r/a":
    ax.set_xlabel(r"$r/a$",fontsize=GFONTSIZE)
    x = prof.data['rmin']
    x = x/max(x)
    r_star = sim.data['r_*/a'][0]
    r_top  = sim.data['r_top/a'][0]
    ax.axvspan(r_star,r_top,facecolor='purple',alpha=0.1)

if rvar == "rho":
    ax.set_xlabel(r"$\rho$",fontsize=GFONTSIZE)
    x = prof.data['rho']

ftag = r'$'+prof.fancy[f][0]+'\;[\mathrm{'+prof.fancy[f][1]+'}]$'
ax.set_ylabel(ftag)

color='black' ; width=2 ; alpha = 0.5 ; label='experiment'
plotcurve(sim,prof)

i=0
prof = profiles_genData('input.profiles.'+str(i))
color='blue' ; width=1 ; alpha = 1 ; label='initial'
plotcurve(sim,prof)

i=sim.n_iterations
color='m' ; width=1 ; alpha = 1 ; label = 'final'
prof = profiles_genData('input.profiles.'+str(i))
plotcurve(sim,prof)

ax.legend(loc=loc)

if ext == 'screen':
    plt.show()
else:
    outfile = f+'.'+ext
    print "Saving plot to "+outfile
    plt.savefig(outfile)

