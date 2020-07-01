import numpy as np
from scipy import interpolate
import matplotlib.path as mp
import matplotlib._contour as _contour

def contourPaths(x, y, Z, levels):
    '''
    :param x: 1D x coordinate
    :param y: 1D y coordinate
    :param Z: 2D data
    :param levels: levels to trace
    :return: list of segments
    '''

    [X,Y]=np.meshgrid(x,y)
    contour_generator = _contour.QuadContourGenerator(X, Y, Z, None, True, 0)

    allsegs = []
    for level in levels:
        segs_=contour_generator.create_contour(level)
        segs=list(map(mp.Path,segs_))
        allsegs.append(segs)

    return allsegs

def paraboloid(x, y, z):
    '''
    z = ax*x^2 + bx*x + ay*y^2 + by*y + c

    NOTE: This function uses only the first 5 points of the x, y, z arrays
    to evaluate the paraboloid coefficients

    :return: ax,bx,ay,by,c
    '''
    if np.any(np.isnan(x.flatten())) or np.any(np.isnan(y.flatten())) or np.any(np.isnan(z.flatten())):
        raise(np.linalg.LinAlgError('paraboloid could not be fit with x=%s y=%s z=%s'%(x,y,z)))
    A=[]
    for k in range(5):
        A.append([x[k]**2,x[k],y[k]**2,y[k],1])
    A=np.array(A)
    ax,bx,ay,by,c=np.dot(np.linalg.inv(A),np.array(z[:5]))
    return ax,bx,ay,by,c

class RectBivariateSplineNaN:
    def __init__(self, Z, R, Q, *args, **kw):
        tmp=Q.copy()
        nan=np.isnan(tmp)
        self.thereAreNaN=False
        if nan.any():
            self.thereAreNaN=True
            tmp[nan]=0
            self.mask=interpolate.RectBivariateSpline(Z, R, nan ,kx=1, ky=1)
        self.spline=interpolate.RectBivariateSpline(Z, R, tmp, *args, **kw)

    def __call__(self,*args):
        tmp=self.spline(*args)
        if self.thereAreNaN:
            mask=self.mask(*args)
            tmp[mask>0.01]=np.nan
        return tmp

    def ev(self,*args):
        tmp=self.spline.ev(*args)
        if self.thereAreNaN:
            mask=self.mask.ev(*args)
            tmp[mask>0.01]=np.nan
        return tmp

def prgen_contour(geqdsk,nrz,levels,psinorm,narc,quiet):

    Rin      = np.linspace(0,geqdsk['RDIM'],geqdsk['NW'])+geqdsk['RLEFT']
    Zin      = np.linspace(0,geqdsk['ZDIM'],geqdsk['NH'])-geqdsk['ZDIM']/2.0+geqdsk['ZMID']
    psirz    = geqdsk['PSIRZ']
    efitpsi0 = geqdsk['SIMAG']
    efitpsi1 = geqdsk['SIBRY']
    efitp    = geqdsk['PRES']
    efitq    = geqdsk['QPSI']
    efitf    = geqdsk['FPOL']
    raxis    = geqdsk['RMAXIS']
    zaxis    = geqdsk['ZMAXIS']
    rlim     = geqdsk['RLIM']
    zlim     = geqdsk['ZLIM']

    if any(np.isnan(rlim)) or any(np.isnan(zlim)):
        print('ERROR: (prgen_contour) rlim/zlim arrays contain NaNs')
        return

    #-----------------------------------------------------------------
    # Change resolution
    #-----------------------------------------------------------------
    if not quiet:
        print('INFO: (prgen_contour) Levels based on psi')

    nrzz = int(nrz*(min(Zin)-max(Zin))/(min(Rin)-max(Rin)))
    
    r2d   = np.linspace(min(Rin),max(Rin),nrz)
    z2d   = np.linspace(min(Zin),max(Zin),nrzz)
    psi2d = RectBivariateSplineNaN(Zin,Rin,psirz)(z2d,r2d)

    if not quiet:
       dres=np.sqrt((r2d[1]-r2d[0])**2+(z2d[1]-z2d[0])**2)
       print('INFO: (prgen_contour) Grid diagonal resolution [m] = {:.5f}'.format(dres))

    #-----------------------------------------------------------------
    # Crop 
    #-----------------------------------------------------------------
    if len(rlim) and len(zlim):
        bbox=[min(rlim),max(rlim),min(zlim),max(zlim)]
        limits=[max([np.argmin(abs(r2d-bbox[0]))-1,0]),
                min([np.argmin(abs(r2d-bbox[1]))+1,len(r2d)-1]),
                max([np.argmin(abs(z2d-bbox[2]))-1,0]),
                min([np.argmin(abs(z2d-bbox[3]))+1,len(z2d)-1])]
        psi2d=psi2d[limits[2]:limits[3],limits[0]:limits[1]]
        r2d=r2d[limits[0]:limits[1]]
        z2d=z2d[limits[2]:limits[3]]

    #-----------------------------------------------------------------
    # Find axis
    #-----------------------------------------------------------------
    if not quiet:
        print('INFO: (prgen_contour) Finding magnetic axis and separatrix')

    dmax = (r2d[-1]-r2d[0])/2.0/1.5

    RR,ZZ = np.meshgrid(r2d-raxis,z2d-zaxis)
    DD = np.sqrt(RR**2+ZZ**2) < dmax
    tmp = psi2d.copy()
    tmp[np.where(DD == 0)] = np.nan

    # figure out sign
    ax = 0
    for k in range(1,int(min([len(r2d)//2,len(z2d)//2])))[::-1]:
        ri = (len(r2d)//2 + np.array([-k, 0, +k, 0, 0])).astype(int)
        zi = (len(z2d)//2 + np.array([0, 0, 0, +k, -k])).astype(int)
        try:
            ax, bx, ay, by, c = paraboloid(r2d[ri], z2d[zi], tmp[zi, ri])
            break
        except np.linalg.LinAlgError:
            pass
    if ax > 0:
        # look for the minimum
        m = np.nanargmin(tmp)
    else:
        # look for the maximum
        m = np.nanargmax(tmp)
        
    Zmi = int(m / psi2d.shape[1])
    Rmi = int(m - Zmi * psi2d.shape[1])

    # fit paraboloid in the vicinity of the grid-based center
    ri = (Rmi + np.array([-1, 0, +1, 0, 0])).astype(int)
    zi = (Zmi + np.array([0, 0, 0, +1, -1])).astype(int)
    ax, bx, ay, by, c = paraboloid(r2d[ri], z2d[zi], psi2d[zi,ri])
    raxis_new = -bx/(2*ax)
    zaxis_new = -by/(2*ay)

    # set as the center value of PSI (based on R0 and Z0)
    psi0 = RectBivariateSplineNaN(z2d, r2d, psi2d).ev(zaxis_new,raxis_new)
   
    #-----------------------------------------------------------------
    # Find separatrix
    #-----------------------------------------------------------------
    accuracy = 9
        
    # separatrix is found by looking for the largest closed path enclosing the magnetic axis

    flxm=np.nanmin(psi2d)
    flxM=np.nanmax(psi2d)
        
    kdbgmax=100
    psi1 = None
    for kdbg in range(kdbgmax):

        flx=0.5*(flxM+flxm)
    
        line=[]
        paths=contourPaths(r2d,z2d,psi2d,[flx])[0]
        for path in paths:
            if not np.isnan(path.vertices[:]).any() and np.allclose(path.vertices[0,0],path.vertices[-1,0]) and np.allclose(path.vertices[0,1],path.vertices[-1,1]):
                path.vertices[0,0] = path.vertices[-1,0] = (path.vertices[0,0]+path.vertices[-1,0])*0.5
                path.vertices[0,1] = path.vertices[-1,1] = (path.vertices[0,1]+path.vertices[-1,1])*0.5
                simplePath=mp.Path(path.vertices[::len(path.vertices[:,0])//10+1,:])
                if (np.max(simplePath.vertices[:,0]) > raxis_new and
                    np.min(simplePath.vertices[:,0]) < raxis_new and
                    np.max(simplePath.vertices[:,1]) > zaxis_new and
                    np.min(simplePath.vertices[:,1]) < zaxis_new):
                    if simplePath.contains_point((raxis_new,zaxis_new)):
                        dR = path.vertices[1,0]-path.vertices[0,0]
                        dZ = path.vertices[1,1]-path.vertices[0,1]
                        orientation = int(np.sign((path.vertices[0,1]-zaxis_new)*dR-(path.vertices[0,0]-raxis_new)*dZ))
                        line=path.vertices[::orientation,:]
                        break
            
        if len(line):
            try:
                # stop condition
                np.testing.assert_array_almost_equal(sep/raxis_new,line/raxis_new,accuracy)
                break
            except Exception:
                pass
            finally:
                sep = line
                psi1 = flx
                flxm = flx
        else:
            open_sep=paths
            flxM=flx

    if kdbg == kdbgmax-1:
        print('WARNING: (prgen_contour) Finding LCFS aborted after %d iterations!'%kdbgmax)

    print('INFO  (prgen_contour) dpsi = {:.9f} [EFIT] {:.9f} [new]'.format(efitpsi1-efitpsi0,psi1-psi0))

    # Separatrix hits edges of computation domain
    if ((np.abs(np.min(sep[:,0])-np.min(r2d)) < 1e-3) or
        (np.abs(np.max(sep[:,0])-np.max(r2d)) < 1e-3) or
        (np.abs(np.min(sep[:,1])-np.min(z2d)) < 1e-3) or
        (np.abs(np.max(sep[:,1])-np.max(z2d)) < 1e-3)):
        print('WARNING: (prgen_contour) New separatrix hits computation boundary: using EFIT separatrix')
        psi0 = efitpsi0 ; psi1 = efitpsi1
            
    #-----------------------------------------------------------
    # Find surfaces
    #-----------------------------------------------------------
        
    if not quiet:
        print('INFO: (prgen_contour) Contour levels = {:d} | n_arc = {:d} | nrz = {:d}'.format(levels,narc,nrz))

    # absolute psi levels to interpolate (psinorm is normalized)
    # packsep < 1.0 increases separatrix packing
    packsep = 1.0
    
    out_psi = (np.linspace(0,psinorm,levels))**packsep*(psi1-psi0)+psi0

    contours = contourPaths(r2d,z2d,psi2d,out_psi)

    RI = np.zeros([narc,levels]) ; ZI = np.zeros([narc,levels])
    r1 = np.zeros([narc])        ; z1 = np.zeros([narc])

    for k,item1 in enumerate(contours):
        if k==0:
            # axis
            r1[:] = raxis_new*np.ones([narc]) ; z1[:] = zaxis_new*np.ones([narc])
        else:
            # all others
            path=item1[-1]
            # Reverse order (compared to original OMFIT order)
            r=path.vertices[::-1,0] ; r[-1] = r[0]
            z=path.vertices[::-1,1] ; z[-1] = z[0]
            # check for contour above separatrix
            if np.average(z) > np.max(sep[:,1]) or np.average(z) < np.min(sep[:,1]):
                path=item1[-2]
                r=path.vertices[::-1,0] ; r[-1] = r[0]
                z=path.vertices[::-1,1] ; z[-1] = z[0]
                
            if any(np.isnan(r*z)):
                print('ERROR: (prgen_contour) NaN encountered')

            # Arc length
            dl = np.sqrt(np.diff(r)**2+np.diff(z)**2)
            larc = np.zeros([len(r)]) ; larc[1:] = np.cumsum(dl)
               
            # Cubic interpolation from fine contour mesh to coarse t-mesh
            t = np.linspace(0,1,narc)*larc[-1]
            cs = interpolate.splrep(larc,r,per=True) ; r1=interpolate.splev(t,cs) 
            cs = interpolate.splrep(larc,z,per=True) ; z1=interpolate.splev(t,cs)

            # Shift elements so that first index at max(R).
            s = np.argmax(r1)
            r1[0:-1] = np.roll(r1[0:-1],-s) ; r1[-1] = r1[0] 
            z1[0:-1] = np.roll(z1[0:-1],-s) ; z1[-1] = z1[0]
            
        RI[:,k] = r1[:] ; ZI[:,k] = z1[:]

        
    efitpsi = np.linspace(out_psi[0],out_psi[-1],len(efitp))
    cs = interpolate.interp1d(efitpsi,efitp,kind='quadratic') ; out_p = cs(out_psi)
    cs = interpolate.interp1d(efitpsi,efitf,kind='quadratic') ; out_f = cs(out_psi)
    #cs = interpolate.interp1d(efitpsi,efitq,kind='quadratic') ; out_q = cs(out_psi)

    # Recalculate q based on definition (and some identities)
    loopint = np.zeros([levels])
    for i in range(narc-1):
       loopint[:] = loopint[:]+(RI[i+1,:]-RI[i,:])*(ZI[i+1,:]+ZI[i,:])/(RI[i+1,:]+RI[i,:])
      
    cs = interpolate.splrep(out_psi,loopint) ; out_q = interpolate.splev(out_psi,cs,der=1)
    out_q = out_f*out_q/(2*np.pi)
    
    return RI,ZI,out_psi,out_q,out_p,out_f

