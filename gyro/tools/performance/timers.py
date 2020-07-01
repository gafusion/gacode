# file processed by 2to3
from __future__ import print_function, absolute_import
from builtins import map, filter, range
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.text as text
#from matplotlib.lines import Line2D
import matplotlib.font_manager
import matplotlib.font_manager
prop = matplotlib.font_manager.FontProperties(size=14)

fig = plt.figure(figsize=(5*3.5,3*3.5))

fig.subplots_adjust(left=0.05,right=0.95)
fig.subplots_adjust(bottom=0.06,top=0.95)
fig.subplots_adjust(wspace=0.35,hspace=0.3)

titles=['Field interp.a',
        'Field_interp.b',
        'Velocity sum',
        'Field explicit',
        'Field implicit',
        'Gyroave h',
        'Implicit he',
        'RHS total',
        'Coll. step',
        'Coll comm',
        'Nonlinear step',
        'Nonlinear comm',
        'Diagnos. allstep',
        'Diagnos. datastep']

counts=[64,128,256,512]

m = len(titles)
n = len(counts)

tick_locs = counts
tick_lbls = []
for i in range(n):
    tick_lbls.append(str(counts[i]))

t = np.zeros(n)

#------------------------------------------------------------
for p in range(m):
    plot = fig.add_subplot(3,5,p+1)
    plot.set_title(titles[p])
    plot.set_xlabel('Number of cores')
    plot.set_ylabel('Elapsed time (s)')
    plot.set_xscale('log',basex=2)
    plot.set_yscale('log')
    plot.grid(True)

    for i in range(n):
        data = np.loadtxt(str(counts[i])+"/out.gyro.timing",skiprows=3)
        t[i] = np.sum(data[2:,p])

    plot.plot(counts,t,'-',color='black')
    plot.set_xticks(tick_locs)
    plot.set_xticklabels(tick_lbls)
#------------------------------------------------------------

#plt.savefig('fig1_flops.pdf')
plt.show()
#------------------------------------------------------------





