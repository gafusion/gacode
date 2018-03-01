from mpmath import mp
import numpy as np
import sys

b0   = mp.sqrt(mp.mpf(sys.argv[1]))
m    = int(sys.argv[2])
type = sys.argv[3]
imom = int(sys.argv[4])

mp.dps = 50

mu = mp.matrix(2*m,1)
a = mp.matrix(m,1)
b = mp.matrix(m,1)
q2  = mp.matrix(m,1)
q2x = mp.matrix(m,1)
p = mp.matrix(m+1,1)

c = mp.matrix(m+1,m+1)

cx = mp.mpf(1)
cb = mp.mpf(1)

if type == 'leg':
   for i in range(2*m):
      mu[i] = mp.mpf(1.0)/(i+1)
elif type == 'lag':
   for i in range(2*m):
      mu[i] = mp.gamma((i+1)/mp.mpf(2))/2 
elif type == 's0':
   for i in range(2*m):
      arg = (i+1)/mp.mpf(2)
      mu[i] = (mp.gamma(arg)-mp.gammainc(arg,b0**2))/2
elif type == 's1':
   for i in range(2*m):
      arg = (i+2)/mp.mpf(2)
      mu[i] = (mp.gamma(arg)-mp.gammainc(arg,b0**2))/2
elif type == 's2':
   for i in range(2*m):
      arg = (i+3)/mp.mpf(2)
      mu[i] = (mp.gamma(arg)-mp.gammainc(arg,b0**2))/2
      
c[0,0] = mp.mpf(1.0)

for n in range(m):
   q2[n]  = mp.mpf(0)
   q2x[n] = mp.mpf(0)
   for i in range(n+1):
      q2[n] = q2[n] + c[i,n]**2*mu[2*i]
      q2x[n] = q2x[n] + c[i,n]**2*mu[2*i+1]
      if i > 0:
         for j in range(i):
           q2[n] = q2[n] + 2*c[i,n]*c[j,n]*mu[i+j] 
           q2x[n] = q2x[n] + 2*c[i,n]*c[j,n]*mu[i+j+1] 
           
   a[n] = q2x[n]/q2[n]
   if n > 0:
      b[n] = q2[n]/q2[n-1]
   else:
      b[0] = mp.mpf(0)
   
   for i in range(n+1):
      if i > 0:
         cx = c[i-1,n]
      else:
         cx = mp.mpf(0.0)
      if n > 0:
         cb = c[i,n-1]
      else:
         cb = mp.mpf(0.0)

      c[i,n+1] = cx-a[n]*c[i,n]-b[n]*cb

   c[n+1,n+1] = mp.mpf(1.0)

# Reverse
for n in range(m+1):
   p[n] = c[m-n,m]
   
x = mp.matrix(mp.polyroots(p,maxsteps=1000,extraprec=200))

# Nodes
if type == 'leg':
   x = b0*x
   
mat = mp.matrix(m,m)
d1 = mp.matrix(m,m)
d2 = mp.matrix(m,m)
w = mp.matrix(m,1)

for i in range(m):
   for j in range(m):
      mat[i,j] = x[i]**j
      d1[i,j] = j*x[i]**(j-1)
      d2[i,j] = j*(j-1)*x[i]**(j-2)

# Inverse
mi = mat**(-1)
d1 = d1*mi
d2 = d2*mi

# True weight function (v^2) for any previous orthogonal poly nodes
mom = mp.matrix(m,1)
for i in range(m):
   mom[i] = mp.gamma((i+1+imom)/mp.mpf(2))/mp.mpf(2)

for i in range(m):
   w[i] = 0.0
   for j in range(m):
      w[i] = w[i]+mom[j]*mi[j,i]

# WEIGHTS
for i in range(m):
   w[i] = w[i]*4/mp.sqrt(mp.pi)*x[i]**(2-imom)

# Write to datafile
fout = open('out.cgyro.egrid','w')
for k in range(m):
    fout.write(mp.nstr(x[k]**2,17)+' '+mp.nstr(w[k],17)+'\n')
for k in range(m):
    for kp in range(m):
        fout.write(mp.nstr(b0*d1[k,kp],17)+' ')
    fout.write('\n')
for k in range(m):
    for kp in range(m):
        fout.write(mp.nstr(b0**2*d2[k,kp],17)+' ')
    fout.write('\n')
    
fout.close()
               
