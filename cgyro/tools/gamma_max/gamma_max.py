# USAGE:
# $ python gamma_max.py <kmin> <kmax> nk
# EXAMPLE:
# $ python gamma_max.py 0.1 0.7 7

import sys
import os
import numpy as np
from gacodefuncs import *

kmin = float(sys.argv[1])
kmax = float(sys.argv[2])
nk   = int(sys.argv[3])

nmpi = 4

g = [] ; x = []

list=kmin+np.linspace(0,1,nk)*(kmax-kmin)

os.system('cgyro -clean')

for i,k in enumerate(list):
    os.system('echo "KY='+str(k)+'" >> input.cgyro')
    os.system('rm -f out.cgyro.tag ; cgyro -e . -n '+str(nmpi)+' > out')
    data = np.genfromtxt('out.cgyro.freq')
    gamma = data[-1][-1]
    g.append(gamma)
    x.append(k)
    print('{:6.4f} {:6.4f}'.format(k,gamma))
    if i > 1:
        d = g[-1]-g[-2]
        if d < 0:
            break

xs,fs = quadratic_max(x[-3:],g[-3:])

print('root: {:6.4f} {:6.4f}'.format(xs,fs))
