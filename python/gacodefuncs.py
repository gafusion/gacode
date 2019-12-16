#----------------------------------------------------------------------
# gacodefuncs.py
#
# PURPOSE:
#  Functions used for computing averages, manipulating strings, etc.
#----------------------------------------------------------------------

# Useful labels
TIME=r'$(c_s/a)\,t$'

#---------------------------------------------------------------
# Generalization of average routine to include variance
def variance(f,t,wmin,wmax):

    import numpy as np
    
    n_time = len(t)

    # Manage case with 2 time points (eigenvalue)
    if len(t) == 2:
        tmin = t[-1]
        tmax = tmin
        ave  = f[-1]
        return ave

    tmin = (1.0-wmin)*t[-1]
    tmax = (1.0-wmax)*t[-1]

    t_window = 0.0
    ave      = 0.0 ; av2 = 0.0
    for i in range(n_time-1):
        if t[i] >= tmin and t[i] <= tmax: 
            ave = ave+0.5*(f[i]+f[i+1])*(t[i+1]-t[i])
            av2 = av2+0.5*(f[i]**2+f[i+1]**2)*(t[i+1]-t[i])
            t_window = t_window+t[i+1]-t[i]

    ave = ave/t_window
    var = np.sqrt((av2/t_window-ave**2))
   
    return ave,var
#---------------------------------------------------------------
#---------------------------------------------------------------
def average(f,t,wmin,wmax):
 
    n_time = len(t)

    # Manage case with 2 time points (eigenvalue)
    if len(t) == 2:
        tmin = t[-1]
        tmax = tmin
        ave  = f[-1]
        return ave

    tmin = (1.0-wmin)*t[-1]
    tmax = (1.0-wmax)*t[-1]

    t_window = 0.0
    ave      = 0.0
    for i in range(n_time-1):
        if t[i] >= tmin and t[i] <= tmax: 
            ave = ave+0.5*(f[i]+f[i+1])*(t[i+1]-t[i])
            t_window = t_window+t[i+1]-t[i]

    ave = ave/t_window

    return ave
#---------------------------------------------------------------
#---------------------------------------------------------------
def average_n(f,t,wmin,wmax,n):
 
    import numpy as np

    ave = np.zeros(n)

    n_time = len(t)
    tmin = (1.0-wmin)*t[-1]
    tmax = (1.0-wmax)*t[-1]

    t_window = 0.0
    for i in range(n_time-1):
        if t[i] >= tmin and t[i] <= tmax: 
            ave[:] = ave[:]+0.5*(f[:,i]+f[:,i+1])*(t[i+1]-t[i])
            t_window = t_window+t[i+1]-t[i]

    ave = ave/t_window

    return ave
#---------------------------------------------------------------
#---------------------------------------------------------------
# Determine index imin,imax for time-averaging window
def iwindow(t,wmin,wmax):
 
    for i in range(len(t)):
        if t[i] < (1.0-wmin)*t[-1]:
            imin = i+1
        if t[i] <= (1.0-wmax)*t[-1]:
            imax = i

    return imin,imax
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
     name = 'He4'
  elif m == 7:
     name = 'Li'
  elif m == 9:
     name = 'Be'
  elif m == 12:
     name = 'C'
  elif m > 180:
     name = 'W'
  else:
     name = '?'

  return name
#---------------------------------------------------------------

#---------------------------------------------------------------
def smooth_pro(x,z,p,n):

    import numpy as np

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
    for i in np.arange(j,0,-1):
        pf[i-1] = pf[i]*np.exp(0.5*(xf[i]-xf[i-1])*(zf[i]+zf[i-1]))


    return xf,pf
#---------------------------------------------------------------

#---------------------------------------------------------------
def extract(d,sd,key,w,spec,moment,norm=False,wmax=0.0,cflux='auto',dovar=False):

   import os
   import re
   import numpy as np
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
         for line in open(sub+'/input.cgyro').readlines():
            if re.match(key,line):
               found = float(line.split('=')[1])
         x.append(found)
         # Get the corresponding flux
         sim = cgyrodata(sub+'/')
         sim.getflux(cflux)
         y = np.sum(sim.ky_flux,axis=(2,3))
         # Flux for input (spec,moment) window w
         ave,var = variance(y[spec,moment,:],sim.t,w,wmax)
         if dovar:
            f.append(var)
         else:
            f.append(ave)
         print('INFO: (extract) Processed data in '+sub)

   # return (scan parameter, flux, variance)
   if norm == True:
      return np.array(x),np.array(f)/sim.dens[spec]
   else:
      return np.array(x),np.array(f)
#---------------------------------------------------------------

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
     name = 'He4'
  elif m == 7:
     name = 'Li'
  elif m == 9:
     name = 'Be'
  elif m == 12:
     name = 'C'
  elif m > 180:
     name = 'W'
  else:
     name = '?'

  return name
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
def theta_indx(theta,theta_plot):

   # Select theta index
   if theta_plot == 1:
      itheta = 0
   else:
       # theta=0 check just to be safe
       if theta == 0.0:
           itheta = theta_plot/2
       else:
           itheta = int((theta+1.0)/2.0*theta_plot)

   print('INFO: (theta_indx) Selected index',itheta+1,'of',theta_plot)
   return itheta
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
   f = open(infile,'r')
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
