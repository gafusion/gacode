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
fig.subplots_adjust(wspace=0.2,hspace=0.2)

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
    tick_lbls.append(str(titles[i]))

t = np.zeros(m)

#------------------------------------------------------------
for p in range(n):
    plot = fig.add_subplot(2,2,p+1)
    plot.set_title('n='+str(counts[p]))
    plot.set_xlabel('Kernel')
    plot.set_ylabel('Elapsed time (s)')
    plot.grid(True)

    for i in range(m):
        data = np.loadtxt(str(counts[p])+"/out.gyro.timing",skiprows=3)
        t[i] = np.sum(data[2:,i])

    plot.bar(np.arange(m)-0.45+1.0,t,width=0.9,alpha=0.4,edgecolor='black')
    plot.set_xlim([0.5,14.5])
    plot.set_xticks(np.arange(m)+1)
#------------------------------------------------------------

#plt.savefig('fig1_flops.pdf')
plt.show()
#------------------------------------------------------------





