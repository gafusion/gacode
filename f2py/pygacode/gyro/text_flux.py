import sys
import numpy as np
from ..gacodefuncs import *
from .data import GYROData

sim       = GYROData('./')
w         = float(sys.argv[1])
wmax      = 0.0
field     = sys.argv[2]
i_moment  = int(sys.argv[3])

n_field   = int(sim.profile['n_field'])
n_kinetic = int(sim.profile['n_kinetic'])

t = sim.t['(c_s/a)t']

# Read data in gbflux_i and make gbflux
sim.read_gbflux_i()

flux = sim.gbflux

# b is collection of all arrays to be plotted
b = np.zeros((len(t),n_kinetic+1))

b[:,0] = t

# Manage field
if field == 's':
    flux0 = np.sum(flux,axis=1)
    ftag = 'TOT   '
else:
    i_field = int(field)
    flux0 = flux[:,i_field,:,:]
    if i_field == 0: 
        ftag = 'ES    '
    if i_field == 1: 
        ftag = 'EM    '
    if i_field == 2: 
        ftag = 'COM   '

# Manage moment
if i_moment == 0: 
    mtag = 'GAMMA [GB]  '
if i_moment == 1: 
    mtag = 'Q [GB]      '
if i_moment == 2: 
    mtag = 'PI [GB]     '
if i_moment == 3: 
    mtag = 'S [GB]      '

ul = '----------  '

line1 = '                '
line2 = '    (cs/a)t     '
line3 = '    '+ul

tag = []

# Manage species
for i in range(n_kinetic):

    b[:,i+1] = flux0[i,i_moment,:]

    if sim.profile['electron_method'] == 2 or  sim.profile['electron_method'] == 4:
        if i == n_kinetic-1:
            stag = 'elec  '
        else:
            stag = 'ion-'+str(i+1)+' '
    if sim.profile['electron_method'] == 1:
        stag = 'ion-'+str(i+1)+' '
    if sim.profile['electron_method'] == 3:
        stag = 'elec  '

    line1 = line1+mtag
    line2 = line2+stag+ftag
    line3 = line3+ul
    tag.append(mtag.strip()+' '+stag+ftag.strip()+': ')

np.set_printoptions(precision=3,suppress=False,threshold=100000)

print(line1)
print(line2)
print(line3)
print(b)

# Determine tmin
imin,imax=iwindow(t,w,wmax)

print

if imin == len(t)-1:
    print("Averaging Window too small.")
else:
    print('Average Window:',str(t[imin])+' < (c_s/a) t < '+str(t[-1]))
    print('')
    for i in range(n_kinetic):
        print(tag[i],average(b[:,i+1],t,w,wmax))

if sim.profile['boundary_method'] == 2:
    print
    print('Buffers have been properly omitted from average.')
