import sys
import string
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

from ..gacodefuncs import *

def prgen_shape(r,z,narc,nf,xplot):

   # Number of theta-points for plotting
   ur = np.zeros(narc) ; uz = np.zeros(narc)
   vr = np.zeros(narc) ; vz = np.zeros(narc)

   # Extrema ; definitions of rmin, rmaj, etc.
   n1 = np.argmax(z) ; m1 = np.argmin(z)
   n2 = np.argmax(r) ; m2 = np.argmin(r)
   za = z[n1] ; zb = z[m1]
   ra = r[n2] ; rb = r[m2]
   zmaj = 0.5*(za+zb) ; zmin = 0.5*(za-zb)
   rmaj = 0.5*(ra+rb) ; rmin = 0.5*(ra-rb)

   if z[1] < z[0]:
      # Reverse order (may be needed)
      r = np.flip(r,0) ; z = np.flip(z,0)

   # Compute generalized angles (new method July 2023)
   eps = 1.0-1e-11
   # (ur,uz): principle angles (discontinuous)
   uz[:] = np.arcsin(eps*(z[:]-zmaj)/zmin)
   ur[:] = np.arccos(eps*(r[:]-rmaj)/rmin)
   
   # Determine correct branches (vr,vz) via extrema
   i0 = np.zeros(4,dtype=int)

   i0[0] = np.argmin(ur)
   i0[1] = np.argmax(uz)
   i0[2] = np.argmax(ur)
   i0[3] = np.argmin(uz)
  
   j1 = i0[1] ; j2 = i0[2] ; j3 = i0[3]

   # Array notation for speed
   vz[:j1]   = uz[:j1]         ; vr[:j1]   = ur[:j1]
   vz[j1:j2] = np.pi-uz[j1:j2] ; vr[j1:j2] = ur[j1:j2]
   vz[j2:j3] = np.pi-uz[j2:j3] ; vr[j2:j3] = 2*np.pi-ur[j2:j3]
   vz[j3:]   = 2*np.pi+uz[j3:] ; vr[j3:]   = 2*np.pi-ur[j3:]

   # Define vr as deviation from vz
   vr[:] = vr[:]-vz[:] ; vr[-1] = vr[0]
   
   x = vz
   # Now:
   #
   # R = rmaj + rmin*cos(x+vr)
   # Z = zmaj + zmin*sin(x)
   #
   # NOTE: kappa=zmin/rmin

   xr = np.zeros(4)
   xr[0] = rmin
   xr[1] = rmaj
   xr[2] = zmin/rmin
   xr[3] = zmaj

   dx = np.diff(x)

   # Compute expansion coefficients (cr):
   #  vr = sum cr*cos(nx)+sr*sin(nx)
   cr = np.zeros(nf+1)
   sr = np.zeros(nf+1)

   cr[0] = moment(narc,vr,np.ones(narc),dx)
   for p in range(1,nf+1):
      cr[p] = moment(narc,vr,np.cos(p*x),dx)
      sr[p] = moment(narc,vr,np.sin(p*x),dx)

   if xplot > 0.0:
      outfile = '{:.3f}'.format(xplot)
      plot_ang(r,z,x,vr,xr,cr,sr,i0,outfile)

   return cr,sr,xr

def prgen_fshape(rd,zd,nf):

    nd = len(rd)

    # Construct equally-spaced poloidal angle
    theta  = np.linspace(0,1,nd)*2*np.pi
    dtheta = theta[1]-theta[0]

    ar = np.zeros(nf+1) ; br = np.zeros(nf+1)
    az = np.zeros(nf+1) ; bz = np.zeros(nf+1)

    ds = dtheta/np.pi

    # Trapezoidal integration spectrally-accurate
    for n in range(nf+1):
       y = np.cos(n*theta[:])*rd[:] ; ar[n] = ds*np.sum(y[:-1])
       y = np.sin(n*theta[:])*rd[:] ; br[n] = ds*np.sum(y[:-1])
       y = np.cos(n*theta[:])*zd[:] ; az[n] = ds*np.sum(y[:-1])
       y = np.sin(n*theta[:])*zd[:] ; bz[n] = ds*np.sum(y[:-1])

    return ar,br,az,bz

#-----------------------------------------------------------------------------
# Old Fourier expansion
#-----------------------------------------------------------------------------
def oldfourier(ri,zi,nf,rnorm):

   print("INFO: (oldfourier) Generating legacy Fourier coefficients")
   npsi = len(rnorm)

   ari = np.zeros([nf+1,npsi]) ; bri = np.zeros([nf+1,npsi])
   azi = np.zeros([nf+1,npsi]) ; bzi = np.zeros([nf+1,npsi])

   for i in range(npsi-1):
      r=ri[:,i+1] ; z=zi[:,i+1]
      ar,br,az,bz = prgen_fshape(r,z,nf)
      ari[:,i+1] = ar[:] ; bri[:,i+1] = br[:]
      azi[:,i+1] = az[:] ; bzi[:,i+1] = bz[:]

   # Repair origin
   ari[0,:] = extrap(rnorm,ari[0,:])
   azi[0,:] = extrap(rnorm,azi[0,:])

   ari[1:,0] = 0.0 ; bri[1:,0] = 0.0
   azi[1:,0] = 0.0 ; bzi[1:,0] = 0.0
      
   u = ari
   u = np.append(u,bri)
   u = np.append(u,azi)
   u = np.append(u,bzi)
   u.tofile('fluxfit.geo')

   return

# Function to compute Fourier integrals
# f,w are periodic
def moment(n,f,w,d):

   s0 = np.sum((f[:-1]*w[:-1]+f[1:]*w[1:])*d[:])
   s1 = np.sum((w[:-1]*w[:-1]+w[1:]*w[1:])*d[:])
   
   return s0/s1
    
def extrap(x,u):
    p = 2
    m = (u[p]-u[p-1])/(x[p]-x[p-1])
    b = u[p]-m*x[p]
    u[0] = b
    return u

def plot_ang(r,z,x,vr,xr,cr,sr,i0,outfile):

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
    #rc('text',usetex=True)
    rc('font',size=18)

    fig = plt.figure(figsize=(18,10))

    # PLOT contour
    ax = fig.add_subplot(121,aspect='equal')
    ax.set_title(r'$\psi='+outfile+'$')
    ax.set_xlabel(r'$R$')
    ax.set_ylabel(r'$Z$')
    ax.grid(which="both",ls=":")

    # Data
    ax.plot(r,z,'--k',linewidth=1,label='data')
    ax.plot(r[i0],z[i0],'o')
    
    # Parameterized contour
    rp = rmaj+rmin*np.cos(x+pr)
    zp = zmaj+zmin*np.sin(x)
    ax.plot(rp,zp,'-m',linewidth=1,label='mapped')
    ax.legend()
    #-------------------------------------------------------
    
    # PLOT angle 
    ax = fig.add_subplot(122)
    ax.set_title(r'$n_\mathrm{harm} = '+str(nf)+'$')
    ax.set_xlabel(r'$x = \theta_Z/(2\pi)$')
    ax.set_ylabel(r'$\theta_R-\theta_Z$')
    ax.grid(which="both",ls=":")

    x=x/(2*np.pi)

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

    # Angle from data
    ax.plot(x,vr,'-k',linewidth=2,alpha=0.3)

    # Angle from parameterization
    ax.plot(x,pr,'-m',linewidth=1)

    ax.set_xlim([0,1])

    plt.tight_layout()
    ofile = 'pnorm_'+outfile+'.pdf'
    print('INFO: (prgen_shape) Writing '+ofile)
    plt.savefig(ofile)
    plt.close()

    err = 100*np.average(abs(vr-pr))

    print('INFO: (prgen_shape) Mapper error (%) = {:.4f}'.format(err)) 
    return

def plot_coef(pnorm,ci,si,xi):
   
   # Latex fonts
   #rc('text',usetex=True)
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
