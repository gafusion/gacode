import numpy as np
from scipy import integrate,optimize,interpolate

try:
   import gapy
   hasgapy = True
except:
   print('WARNING: (banana_width.py) gapy.so not found.  Using circular geometry.')
   hasgapy = False

data=np.loadtxt('out.locpargen')

# Define theta grid
n = 32
theta = np.linspace(0,2*np.pi,n)-np.pi

if hasgapy:
   gapy.geo.geo_q_in = data[0]
   gapy.geo.geo_rmin_in = data[1]
   gapy.geo.geo_rmaj_in = data[2]
   gapy.geo.geo_drmaj_in = data[3]
   gapy.geo.geo_kappa_in = data[4]
   gapy.geo.geo_interp(theta,True)
   b_gapy = gapy.geo.geo_b
   s_gapy = gapy.geo.geo_gsin*gapy.geo.geo_g_theta*gapy.geo.geo_grad_r
   s_gapy[0] = 0.0 ; s_gapy[-1] = 0.0
else:
   eps = data[1]/data[2]
   b_gapy = 1.0/(1.0+eps*np.cos(theta))
   s_gapy = np.sin(theta)
    
b_s = interpolate.splrep(theta,b_gapy,per=1)
s_s = interpolate.splrep(theta,s_gapy,per=1)

def func(t):

   b = interpolate.splev(t,b_s)
   s = interpolate.splev(t,s_s)
   f = (1+l*b/2)/np.sqrt(1.0-l*b)/b*s
   return f

def bounce(t):
   b = interpolate.splev(t,b_s)
   return 1-l*b

lb   = 1.0/interpolate.splev(np.pi,b_s)
lmax = 1.0/interpolate.splev(0.0,b_s)
l    = lb*0.9999

if l < lb:
   t0 = np.pi
elif l < lmax:
   t0 = optimize.fsolve(bounce,0.01)[0]
else:
   raise ValueError('l too large')

y,err = integrate.quad(func,0,t0)
db = data[5]
print('# orbit width/a = {:.3f}'.format(np.abs(y*db)))

