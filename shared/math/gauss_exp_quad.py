import sys
from mpmath import mp

emax = sys.argv[1]
n    = int(sys.argv[2])

# Precision set to huge value (>= 50)
prec = 50

b = mp.sqrt(emax)
prec = 50
n    = n+1

#-----------------------
nplot=n*8
mp.dps=prec

pa = mp.matrix(n,1)
pb = mp.matrix(n,1)
qa = mp.matrix(n,1)
ra = mp.matrix(n,1)
alpha = mp.matrix(n,1)
beta = mp.matrix(n,1)
gamma = mp.matrix(n,1)
xl = mp.matrix(n,1)
xr = mp.matrix(n,1)
xk = mp.matrix(n-1,1)
xf = mp.matrix(n-1,n-1)
wk = mp.matrix(n-1,1)

z  = mp.matrix(n-1,n-1)
zp = mp.matrix(n-1,n-1)
zpp = mp.matrix(n-1,n-1)

alpha[0] = (mp.exp(-b*b)-1)/(mp.sqrt(mp.pi)*mp.erf(b))
beta[0] = 0
gamma[0] = mp.sqrt(mp.pi)/2*mp.erf(b)

pa[0] = 1
pb[0] = 1

# Compute recursion coefficients

for k in range(n-1):
    if k > 0:
        pa[k+1] = alpha[k]*pa[k]+beta[k]*pa[k-1]
        pb[k+1] = (b+alpha[k])*pb[k]+beta[k]*pb[k-1]
    else:
        pa[k+1] = alpha[k]*pa[k]
        pb[k+1] = (b+alpha[k])*pa[k]
    gamma[k+1] = ((k+1)*gamma[k]-(mp.exp(-b*b)*pb[k+1]*pb[k]-pa[k+1]*pa[k]))/2
    alpha[k+1] = 1/(2*gamma[k+1])*(mp.exp(-b*b)*pb[k+1]*pb[k+1]-pa[k+1]*pa[k+1])
    beta[k+1] = -gamma[k+1]/gamma[k]

#for k in range(n):
#    print k,alpha[k],beta[k],g[k]

# Print polynomial data

def poly(x0):
    global alpha
    global beta
    global n
    global pa
    global qa
    global ra

    pa[0] = 1
    qa[0] = 0
    ra[0] = 0
    for k in range(n-1):
        if k > 0:
            pa[k+1] = (x0+alpha[k])*pa[k]+beta[k]*pa[k-1]
            qa[k+1] = (x0+alpha[k])*qa[k]+beta[k]*qa[k-1]+pa[k]
            ra[k+1] = (x0+alpha[k])*ra[k]+beta[k]*ra[k-1]+2*qa[k]
        else:
            pa[k+1] = (x0+alpha[k])*pa[k]
            qa[k+1] = (x0+alpha[k])*qa[k]+pa[k]
            ra[k+1] = (x0+alpha[k])*ra[k]+2*qa[k]

    f = pa[n-1]
    g = qa[n-1]
    return f,g
    
x0 = mp.mpf(0)
f,g = poly(x0)
k = 0
for i in range(nplot):
    x0 = i*b/(nplot-1)
    f0 = f
    f,g = poly(x0)

    if f*f0 < 0:
        xl[k] = x0-b/(nplot-1)
        xr[k] = x0
        k = k+1 
        
if k < n-1:
    print 'Not all roots bracketed'
    sys.exit()
    
# Solve for roots using Newton's method
for k in range(n-1):
    x0 = (xl[k]+xr[k])/2
    f,g = poly(x0)
    while abs(f/g) > 10.0**(2-prec):
        x0  = x0-f/g
        f,g = poly(x0)
    xk[k] = x0
    wk[k] = gamma[n-2]/(g*pa[n-2])
    # Capture data for derivative matrix
    for kp in range(n-1):
        z[k,kp]   = pa[kp]
        zp[k,kp]  = qa[kp]
        zpp[k,kp] = ra[kp]
      
zp  = zp*mp.powm(z,-1)*b
zpp = zpp*mp.powm(z,-1)*b*b

# Redefine weight and renormalize last element

wg = mp.matrix(n-1,1)
s = mp.mpf(0)
for k in range(n-1):
   wg[k] = wk[k]*4/mp.sqrt(mp.pi)*xk[k]**2
   s = s+wg[k]

wg[n-2] = wg[n-2]+mp.mpf(1)-s

# Write to datafile
fout = open('out.cgyro.egrid','w')
for k in range(n-1):
    fout.write(mp.nstr(xk[k]**2,17)+' '+mp.nstr(wg[k],17)+'\n')
for k in range(n-1):
    for kp in range(n-1):
        fout.write(mp.nstr(zp[k,kp],17)+' ')
    fout.write('\n')
for k in range(n-1):
    for kp in range(n-1):
        fout.write(mp.nstr(zpp[k,kp],17)+' ')
    fout.write('\n')
    
fout.close()
               
