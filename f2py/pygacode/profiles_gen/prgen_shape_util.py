import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# Function to compute Fourier integrals
# f,w are periodic
def moment(n,f,w,d):

   s0 = np.sum((f[:-1]*w[:-1]+f[1:]*w[1:])*d[:-1])
   s1 = np.sum((w[:-1]*w[:-1]+w[1:]*w[1:])*d[:-1])
   
   return s0/s1
    
def extrap(x,u):
    p = 2
    m = (u[p]-u[p-1])/(x[p]-x[p-1])
    b = u[p]-m*x[p]
    u[0] = b
    return u

def zero(x,u):
    u[0] = 0.0
    return u

def plot_ang(r,z,x,vr,xr,cr,sr,outfile):

    nf = len(cr)-1
    
    rmin = xr[0] 
    rmaj = xr[1]
    zmin = xr[2]*rmin
    zmaj = xr[3] 

    # Angles approximated by Fourier series
    n_arc = len(x)
    pr = np.ones(n_arc)*cr[0]
    for p in range(1,nf+1):
        pr = pr+cr[p]*np.cos(p*x)+sr[p]*np.sin(p*x)

    # Latex fonts
    rc('text',usetex=True)
    rc('font',size=18)

    fig = plt.figure(figsize=(18,9))

    # PLOT contour
    ax = fig.add_subplot(121,aspect='equal')
    ax.set_xlabel(r'$R$')
    ax.set_ylabel(r'$Z$')
    ax.grid(which="both",ls=":")

    # Data
    ax.plot(r,z,'-k',linewidth=2,alpha=0.3)
    # Parameterized contour
    rp = rmaj+rmin*np.cos(x+pr)
    zp = zmaj+zmin*np.sin(x)
    ax.plot(rp,zp,'-m',linewidth=1)

    # PLOT angle 
    ax = fig.add_subplot(122)
    ax.set_xlabel(r'$x = \theta_Z/(2\pi)$')
    ax.set_ylabel(r'$\theta_R-\theta_Z$')
    ax.grid(which="both",ls=":")

    x=x/(2*np.pi)

    ax.set_xlim([0,1])

    ax.plot(x,pr,'-m',linewidth=1)

    # Shift elements so that x<1.0
    s=-2
    for i in range(n_arc):
        if x[i] > 1.0:
            s = i-1
            break

    if s > -2:
        x[0:-1]  = np.roll(x[0:-1],-s)  ; x[-1] = x[0] 
        vr[0:-1] = np.roll(vr[0:-1],-s) ; vr[-1] = vr[0]
        pr[0:-1] = np.roll(pr[0:-1],-s) ; pr[-1] = pr[0]

    for i in range(n_arc-1):
        if x[n_arc-2-i] > x[n_arc-1-i]:
            x[n_arc-2-i] = x[n_arc-2-i]-1.0

    # Angles from data
    ax.plot(x,vr,'-k',linewidth=2,alpha=0.3)

    plt.tight_layout()
    print('INFO: (plot_ang) Writing '+outfile)
    plt.savefig(outfile)
    plt.close()

    return

def plot_coef(pnorm,ci,si,xi):
   
   # Latex fonts
   rc('text',usetex=True)
   rc('font',size=18)

   fig = plt.figure(figsize=(14,12))

   label=['r','R','\kappa','Z']
   for i in range(4):
      ax = fig.add_subplot(3,4,i+1)
      ax.set_xlabel(r'$\psi$')
      ax.set_title(r'$'+label[i]+'$')
      ax.grid(which="both",ls=":")
      ax.set_xlim([0,1])
      u = xi[i,:] 
      ax.plot(pnorm,u,'-r',linewidth=1,alpha=1)

   label=['c_0','c_1','c_2','c_3']
   for i in range(4):
      ax = fig.add_subplot(3,4,i+5)
      ax.set_xlabel(r'$\psi$')
      ax.set_title(r'$'+label[i]+'$')
      ax.grid(which="both",ls=":")
      ax.set_xlim([0,1])
      u = ci[i,:]
      ax.plot(pnorm,u,'-r',linewidth=1,alpha=1)

   label=['-','\delta','-\zeta','s_3']
   for i in range(4):
      ax = fig.add_subplot(3,4,i+9)
      ax.set_xlabel(r'$\psi$')
      ax.set_title(r'$'+label[i]+'$')
      ax.grid(which="both",ls=":")
      ax.set_xlim([0,1])
      if i > 0:
         u = si[i,:]
         ax.plot(pnorm,u,'-r',linewidth=1,alpha=1)
            
   plt.tight_layout()
   plt.savefig('plot_coef.png')
   plt.close()

   return
