# USAGE:
# $ python peak.py <kmin> <kmax> nk
# EXAMPLE:
# $ python peak.py 0.1 0.7 7

import sys
import os
import numpy as np

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

# f(x) at 3 points
f1 = g[-3] ; f2 = g[-2] ; f3 = g[-1]
x1 = x[-3] ; x2 = x[-2] ; x3 = x[-1]

# Extrema fs=f(xs) based on 3-point fit to parabola 
xs = (f1-f3)/2.0/(f3-2*f2+f1)
xs = x2+xs*(x3-x2)

fs = f2+(f3-f1)**2/8.0/(2*f2-f3-f1)

print('root: {:6.4f} {:6.4f}'.format(xs,fs))
