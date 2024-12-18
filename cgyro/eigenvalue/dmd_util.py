import numpy as np
from scipy.signal import butter,filtfilt

# Map from (kx,theta) to theta_* 
def map1d(f2d,q):
    nr    = f2d.shape[0]
    nt    = f2d.shape[1]
    ntime = f2d.shape[2]
    px = np.arange(nr)-nr//2
    f1d = np.zeros([nr,nt,ntime],dtype=complex)
    anorm = f1d[nr//2,nt//2,:]

    for ir in range(nr):
        f1d[ir,:,:] = f2d[ir,:,:]*np.exp(-2*np.pi*1j*px[ir]*q)
    f1d = f1d.reshape(nr*nt,ntime)

    return f1d,anorm

# downsample with averaging
def downsample(f,n):
    npar = f.shape[0]
    nave = f.shape[1]

    # array end trim
    m = (nave//n)*n

    # downsample
    fdown = np.mean(f[:,:m].reshape(npar,-1,n),axis=2)

    return fdown
