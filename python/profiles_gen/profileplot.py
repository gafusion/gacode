"""Generate a plot of a profile function"""

import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from profiles_gen.data import profiles_genData

rc('text',usetex=True)
rc('font',size=18)

rvar    = sys.argv[1]
infiles = sys.argv[2]
plots   = sys.argv[3]
ftype   = sys.argv[4]
loc     = int(sys.argv[5])
t       = sys.argv[6]
title   = sys.argv[7]

plotvec = string.splitfields(plots,',')
filevec = string.splitfields(infiles,',')

n = len(plotvec)

for j in range(n):

    if len(t.split(',')) == n:
        tlabel = r'$\;\mathrm{'+t.split(',')[j]+'}$'
    else:
        tlabel = ''

    if filevec[j] == '.':
        filevec[j] = filevec[j-1]

    if plotvec[j] == '.':
        plotvec[j] = plotvec[j-1]

    prof = profiles_genData(filevec[j])
    tag  = plotvec[j]
    keys  = sorted(prof.data.keys())

    success = 0
    for i in range(len(keys)):
        if tag == keys[i].split()[0]:
            success = 1
            fulltag = keys[i]

    if success == 0:
        print "ERROR: (profiles_gen_plot) Bad profile = "+tag
        sys.exit()

    if j==0:
        fig = plt.figure(figsize=(12,7))
        ax  = fig.add_subplot(111)

        if title != 'null':
            ax.set_title(r'$'+title+'$',size=18)

        if rvar == "r":
            ax.set_xlabel(r"$r \, \mathrm{(m)}$")
        if rvar == "r/a":
            ax.set_xlabel(r"$r/a$")
        if rvar == "rho":
            ax.set_xlabel(r"$\rho$")

        ax.grid(which="majorminor",ls=":")
        ax.grid(which="major",ls=":")

    if rvar == "r":
        x = prof.data['rmin']

    if rvar == "r/a":
        x = prof.data['rmin']
        x = x/max(x)

    if rvar == "rho":
        x = prof.data['rho']

    ftag = r'$'+prof.fancy[tag][0]+'\;[\mathrm{'+prof.fancy[tag][1]+'}]$'
    ax.plot(x,prof.data[fulltag],'o-',label=ftag+tlabel)

ax.legend(loc=loc)

if ftype == 'screen':
    plt.show()
else:
    print "Saving plot to "+ftype
    outfile = ftype
    plt.savefig(outfile)
