#----------------------------------------------------------------------
# gacodefuncs.py
#
# PURPOSE:
#  Functions used for computing averages, manipulating strings, etc.
#----------------------------------------------------------------------

import numpy as np
import time

# Useful labels
TIME=r'$(c_s/a)\,t$'

TEXPHI  = r'\delta\phi'
TEXAPAR = r'{A_\parallel}'
TEXBPAR = r'{B_\parallel}'
TEXDN   = r'\delta n'
TEXDE   = r'\delta E'
TEXDV   = r'\delta v'

#---------------------------------------------------------------
# Legacy average function (should go away eventually)
def average(f,t,w,wmax):
    imin,imax = time_index(t,str(w))
    return time_average(f,t,imin,imax)
#---------------------------------------------------------------

#---------------------------------------------------------------
# Determine index imin,imax for time-averaging window
def time_index(t,w):

    tvec = np.array(w.split(',')).astype(float)

    if len(tvec) == 1:
        imax = len(t)-1
        imin = np.argmin(abs(t-t[-1]*(1-tvec[0])))
    else:
        imax = np.argmin(abs(t-tvec[1]))
        imin = np.argmin(abs(t-tvec[0]))
    
    return imin,imax
#---------------------------------------------------------------

#---------------------------------------------------------------
# Compute time average (of 1D, 2D or 3D array)
def time_average(f,t,imin,imax):
    
    T = t[imax]-t[imin]
    dt = np.diff(t)[imin:imax]/T

    if f.ndim == 1:
        # 1D array
        sf = 0.5*(f[imin:imax]+f[imin+1:imax+1])
        ave = np.sum(sf*dt)
    elif f.ndim ==2:
        # 2D array (average over last dimension)
        sf = 0.5*(f[:,imin:imax]+f[:,imin+1:imax+1])
        ave = np.sum(sf*dt,axis=-1)
    else:
        # 3D array (average over last dimension)
        sf = 0.5*(f[:,:,imin:imax]+f[:,:,imin+1:imax+1])
        ave = np.sum(sf*dt,axis=-1)
        
    return ave
#---------------------------------------------------------------

#------------------------------------------------------
# Construct an explicit integer list based on string
def str2list(str):

    nvec = []
    for i in str.split(','):
        if '-' in i:
            v = i.split('-')
            for j in range(int(v[0]),int(v[1])+1):
                nvec.append(j)
        else:
            nvec.append(int(i))

    return nvec
#------------------------------------------------------

#------------------------------------------------------
# Set axis limits
def setlimits(a,fmin,fmax):

    fmin0=a[0]
    fmax0=a[1]

    if fmin != 'auto':
        fmin0=float(fmin)
    if fmax != 'auto':
        fmax0=float(fmax)

    return fmin0,fmax0
#------------------------------------------------------

#---------------------------------------------------------------
# Determine species name (returnval) from mass and charge
def specmap(m_in,z_in):

  # Assume Deuterium normalization
  m = int(m_in*2)
  z = int(z_in)

  if z < 0:
    name = 'e'
  elif m == 1:
     name = 'H'
  elif m == 2:
      name = 'D'
  elif m == 3:
     if z == 1:
        name = 'T'
     elif z == 2:
        name = 'He3'
     else:
        name = '?'
  elif m == 4:
     name = 'He'
  elif m == 7:
     name = 'Li'
  elif m == 9:
     name = 'Be'
  elif m == 12:
     name = 'C'
  elif m == 14:
     name = 'N'
  elif m == 16:
     name = 'O'
  elif m == 20:
     name = 'Ne'
  elif m == 40:
     name = 'Ar'
  elif m > 180:
     name = 'W'
  else:
     name = '?'

  return name

#---------------------------------------------------------------
# Generate time window text string
def wintxt(imin,imax,t,usec=0,fc=0,field=0):

    if usec:
        cstr = 'domain/half'
    else:
        cstr = 'domain/full'

    if fc == 0:
        fstr = 'field/all'
    elif field == 0:
        fstr = '\phi'
    elif field == 1:
        fstr = 'A_\parallel'
    else:
        fstr = 'B_\parallel'

    pre = '['+cstr+','+fstr+']'
    win = '[{:.1f} < (c_s/a) t < {:.1f}]'.format(t[imin],t[imax])
   
    mpre = r'$\mathrm{'+pre+'}$'
    mwin = r'$'+win+'$'

    print('INFO: (wintxt.py) Average Window: '+win)
      
    return mpre,mwin

#---------------------------------------------------------------
# Generate smooth profile using exponential integration
def smooth_pro(x,z,p,n,type='log'):

    nx = len(x)
    xf = np.zeros((nx-1)*n+1)
    zf = np.zeros((nx-1)*n+1)
    pf = np.zeros((nx-1)*n+1)
    j = 0
    for i in range(nx-1):
        for m in range(n):
            u = m/(1.0*n)
            xf[j] = x[i]*(1-u)+x[i+1]*u
            zf[j] = z[i]*(1-u)+z[i+1]*u
            j = j+1

    xf[j] = x[nx-1]
    zf[j] = z[nx-1]
    pf[j] = p[nx-1]

    # Exponential integration to obtain smooth profiles
    if type == 'log':
        for i in np.arange(j,0,-1):
            pf[i-1] = pf[i]*np.exp(0.5*(xf[i]-xf[i-1])*(zf[i]+zf[i-1]))
    else:
        for i in np.arange(j,0,-1):
            pf[i-1] = pf[i]-0.5*(xf[i]-xf[i-1])*(zf[i]+zf[i-1])
            
    return xf,pf
#---------------------------------------------------------------

#---------------------------------------------------------------
def extract(d,sd,key,w,spec,moment,norm=False,wmax=0.0,cflux='auto',dovar=False):

   import os
   import re
   from cgyro.data import cgyrodata

   # d        = directory
   # sd       = prefix of subdirectory ('a' for a1,a2,a3)
   # key      = key to scan (for example, 'GAMMA_P')
   # w        = time-averaging width
   # spec     = (0 ...)
   # moment   = (0 ...)
   # norm     = True (density normalization)
   # wmax     = time-averaging minimum
   # cflux    = 'on'/'off'/'auto'
   # dovar    = True/False (variance calculation)

   x = []
   f = []
   for i in range(64):
      sub = d+'/'+sd+str(i)+'/'
      if os.path.isdir(sub) == True:
         # If this is a directory, get the key value
         with open(sub+'/input.cgyro') as in_cgyro:
            for line in in_cgyro.readlines():
               if re.match(key,line):
                  found = float(line.split('=')[1])
         x.append(found)
         # Get the corresponding flux
         sim = cgyrodata(sub+'/')
         sim.getflux(cflux)
         y = np.sum(sim.ky_flux,axis=(2,3))
         # Flux for input (spec,moment) window w
         ave,var = variance(y[spec,moment,:],sim.t,w,wmax)
         if norm == True:
             var = var/sim.dens[spec]
             ave = ave/sim.dens[spec]
         if dovar:
            f.append(var)
         else:
            f.append(ave)
         print('INFO: (extract) Processed data in '+sub)

   # return (scan parameter, flux, variance)
   return np.array(x),np.array(f)
#---------------------------------------------------------------

#---------------------------------------------------------------
# Determine species name (returnval) from mass and charge
def tag_helper(mass,z,moment):

  u=specmap(mass,z)

  # Set filename root and title
  isfield = True
  if (moment == 'n'):
      fdata = '.cgyro.kxky_n'
      title = r'${\delta \mathrm{n}}_'+u+'$'
      isfield = False
  elif (moment == 'e'):
      fdata = '.cgyro.kxky_e'
      title = r'${\delta \mathrm{E}}_'+u+'$'
      isfield = False
  elif (moment == 'v'):
      fdata = '.cgyro.kxky_v'
      title = r'${\delta \mathrm{v}}_'+u+'$'
      isfield = False
  elif (moment == 'phi'):
      fdata = '.cgyro.kxky_phi'
      title = r'$\delta\phi$'
  elif (moment == 'apar'):
      fdata = '.cgyro.kxky_apar'
      title = r'$\delta A_\parallel$'
  elif (moment == 'bpar'):
      fdata = '.cgyro.kxky_bpar'
      title = r'$\delta B_\parallel$'

  return fdata,title,isfield
#---------------------------------------------------------------

#---------------------------------------------------------------
# Get time vector from commmand line option
def time_vector(istr,nt):

   if istr == '-1':
      ivec = [nt]
   elif istr == 'all':
      ivec = list(range(nt))
   else:
      ivec = str2list(istr)

   return ivec
#---------------------------------------------------------------

#---------------------------------------------------------------
def mkfile(ext):

    s=ext.split('.')
    if len(s) == 2:
        pre   = s[0]
        ftype = s[1]
    else:
        pre = ''
        ftype = s[0]

    return pre,ftype
#---------------------------------------------------------------

#---------------------------------------------------------------
def parameter_val(infile,par):

   val = 'null'
   with open(infile,'r') as f:
      for line in f:
         x = line.split()
         if x[1] == par:
            val = x[0]

   return val
#---------------------------------------------------------------

#---------------------------------------------------------------
def quadratic_max(x,g):
    
    # f(x) at 3 points
    f1 = g[-3] ; f2 = g[-2] ; f3 = g[-1]
    x1 = x[-3] ; x2 = x[-2] ; x3 = x[-1]

    # Extrema fs=f(xs) based on 3-point fit to parabola 
    xs = (f1-f3)/2.0/(f3-2*f2+f1)
    xs = x2+xs*(x3-x2)

    fs = f2+(f3-f1)**2/8.0/(2*f2-f3-f1)

    return xs,fs
#---------------------------------------------------------------

#---------------------------------------------------------------
def indx_theta(i,n):

   # Select theta index
   if n == 1:
      itheta = 0
      thetapi = 0.0
   elif i == -1:
      itheta = n//2
      thetapi = 0.0
   else:
      itheta = i
      thetapi = -1+2.0*itheta/n

   print('INFO: (indx_theta) Selected theta index {:d} of {:d}-{:d} : theta={:.2f}pi'.
         format(itheta,0,n-1,thetapi))

   return itheta,thetapi
#---------------------------------------------------------------

def shift_fourier(f,imin,imax):

    nx = f.shape[0]+1
    nn = f.shape[1]
    nt = f.shape[2]
    
    y1 = np.zeros([nn])
    y2 = np.zeros([nn])

    ephi  = np.zeros([2*nx,nt],dtype=complex)
    ephip = np.zeros([2*nx,nt],dtype=complex)

    wpos = np.zeros([2*nx])
    wneg = np.zeros([2*nx])
    wneg[nx//2+1:3*nx//2] = 1.0
    wneg[nx//2]  =  0.5
    wneg[3*nx//2] = 0.5
    wpos[:] = 1-wneg[:]
        
    phi  = np.zeros([nx,nt],dtype=complex)
    phip = np.zeros([nx,nt],dtype=complex)

    for n in range(nn):
        for p in range(1,nx):
            phi[p,:] = f[p-1,n,:]
            phip[p,:] = -(p-nx//2)*f[p-1,n,:]
            
        ephi[nx//2:3*nx//2,:]  = phi[:,:]
        ephip[nx//2:3*nx//2,:] = phip[:,:]            

        # NOTE: We use *inverse* FFT (ifft) for correct +sign convention of
        #       the exponent. Also note order convention:
        #       - a[0] = p=0
        #       - a[1:nx/2] = p > 0
        #       - a[nx/2:n] = p < 0
         
        phi_T = np.fft.ifft(np.fft.ifftshift(ephi,axes=0),axis=0)
        phip_T = np.fft.ifft(np.fft.ifftshift(ephip,axes=0),axis=0)

        # Not quite correct time-average (correct only for fixed dt)
        pn_t = np.zeros([2*nx])
        pd_t = np.zeros([2*nx])
        for jt in np.arange(imin,imax+1):
            pn_t[:] = pn_t[:] + np.real(np.conj(phi_T[:,jt])*phip_T[:,jt])
            pd_t[:] = pd_t[:] + np.real(np.conj(phi_T[:,jt])*phi_T[:,jt])
    
        # Shift in -gamma domain (standard order: p=0 is 0th index)
        pn = np.sum(pn_t[:]*wneg[:])
        pd = np.sum(pd_t[:]*wneg[:])
            
        y2[n] = pn/pd

        # Shift in central domain 
        pn = np.sum(pn_t[:]*wpos[:])
        pd = np.sum(pd_t[:]*wpos[:])
         
        y1[n] = pn/pd

    return y1,y2


def shift_legendre(f,imin,imax):

    import scipy.special as sp

    nx = f.shape[0]+1
    nn = f.shape[1]
    nt = f.shape[2]
    n0 = nx//2
    nk = 2*n0

    y1 = np.zeros([nn])
    y2 = np.zeros([nn])

    c1 = np.zeros([nk],dtype=complex)
    c2 = np.zeros([nk],dtype=complex)

    mat1 = np.zeros([nk,n0-1])
    mat2 = np.zeros([nk,n0-1])
    kvec = np.arange(nk)
    pvec = np.arange(1,n0)
    z = pvec*np.pi/2
    for k in kvec:
        mat1[k,:] = sp.spherical_jn(k,z)
        mat2[k,:] = mat1[k,:]*(-1)**pvec[:]

    ai = 1j**kvec   
    ak = 2*kvec+1   

    c1m = np.zeros(nk,dtype=complex)
    c2m = np.zeros(nk,dtype=complex)

    for n in range(nn):
        n_all1 = n_all2 = 0.0
        d_all1 = d_all2 = 0.0
        for jt in np.arange(imin,imax+1):

            y = f[:,n,jt]

            phim = np.flip(y[0:n0-1])
            phi0 = y[n0-1]
            phip = y[n0:]

            #print(len(phim),len(phip))

            mp1 = np.matmul(mat1,phip)
            mm1 = np.matmul(mat1,phim)
            mp2 = np.matmul(mat2,phip)
            mm2 = np.matmul(mat2,phim)

            c1[:] = ak[:]*(mp1[:]*ai[:]+mm1[:]*np.conj(ai[:]))
            c2[:] = ak[:]*(mp2[:]*ai[:]+mm2[:]*np.conj(ai[:]))
            
            c1[0] = c1[0]+phi0
            c2[0] = c2[0]+phi0

            c1m[0:2] = c1[0:2]
            c2m[0:2] = c2[0:2]
            for m in range(2,nk):
                c1m[m] = c1m[m-2]+c1[m]
                c2m[m] = c2m[m-2]+c2[m]
               
            n1 = 2*np.sum(c1m[0:nk-1]*np.conj(c1[1:nk]))
            n2 = 2*np.sum(c2m[0:nk-1]*np.conj(c2[1:nk]))
            
            d1 = 2*np.sum((np.abs(c1[:]))**2/ak[:])
            d2 = 2*np.sum((np.abs(c2[:]))**2/ak[:])

            n_all1 = n_all1 + np.imag(n1)
            n_all2 = n_all2 + np.imag(n2)
            d_all1 = d_all1 + d1
            d_all2 = d_all2 + d2

            # Derivative
            # demoninator is <phi | phi> = sum_m 2/(2m+1) |cm|^2 

        y1[n] = n_all1/d_all1*(2.0/np.pi)
        y2[n] = n_all2/d_all2*(2.0/np.pi)

    return y1,y2

