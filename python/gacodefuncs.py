#----------------------------------------------------------------------
# gacodeplotfuncs.py
#
# PURPOSE:
#  Functions used for computing averages, manipulating strings, etc.
#----------------------------------------------------------------------

#---------------------------------------------------------------
def average(f,t,window):
 
    n_time = len(t)

    # Manage case with 2 time points (eigenvalue)
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
#---------------------------------------------------------------
def average_n(f,t,window,n):
 
    import numpy as np

    ave = np.zeros(n)

    n_time = len(t)
    tmin = (1.0-window)*t[n_time-1]
    tmax = t[n_time-1]

    t_window = 0.0
    for i in range(n_time-1):
        if t[i] > tmin:
            ave[:] = ave[:]+0.5*(f[:,i]+f[:,i+1])*(t[i+1]-t[i])
            t_window = t_window+t[i+1]-t[i]

    ave = ave/t_window

    return ave
#---------------------------------------------------------------
#---------------------------------------------------------------
# Determine index imin for time-averaging window
def iwindow(t,window):
 
    imin=0
    for i in range(len(t)):
        if t[i] < (1.0-window)*t[-1]:
            imin = i+1

    return imin
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
# Determine index imin for time-averaging window
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
