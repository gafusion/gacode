import sys
import string
import numpy as np
from gyro.data import GYROData

#---------------------------------------------------------------
def average(f,t,window):
 
    n_time = len(t)

    # Manage eigenvalue case (2 time points)
    if len(t) == 2:
        tmin = t[n_time-1]
        tmax = tmin
        ave  = f[n_time-1]
        return ave
    
    tmin = (1.0-window)*t[n_time-1]
    tmax = t[n_time-1]

    t_window = 0.0
    ave      = 0.0
    for i in range(n_time-1):
        if t[i] > tmin: 
            ave = ave+0.5*(f[i]+f[i+1])*(t[i+1]-t[i])
            t_window = t_window+t[i+1]-t[i]

    ave = ave/t_window

    return ave
#---------------------------------------------------------------
 
sim       = GYROData(sys.argv[1])
field     = sys.argv[2]
i_moment  = int(sys.argv[3])
window    = float(sys.argv[4])

n_field   = int(sim.profile['n_field'])
n_kinetic = int(sim.profile['n_kinetic'])

t    = sim.t['(c_s/a)t']

# Read data in gbflux_i and make gbflux
sim.read_gbflux_i()
sim.make_gbflux()

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
    tag.append(string.strip(mtag)+' '+stag+string.strip(ftag)+': ')

np.set_printoptions(precision=3,suppress=False,threshold=100000)

print line1
print line2
print line3
print b

# Determine tmin
imin=0
for i in range(len(t)):
    if t[i] < (1.0-window)*t[len(t)-1]:
        imin = i+1

print

if imin == len(t)-1:
    print "Averaging Window too small." 
else:
    print 'Average Window:',str(t[imin])+' < (c_s/a) t < '+str(t[-1])
    print
    for i in range(n_kinetic):
        print tag[i],average(b[:,i+1],t,window)

if sim.profile['boundary_method'] == 2:
    print
    print 'Buffers have been properly omitted from average.'
