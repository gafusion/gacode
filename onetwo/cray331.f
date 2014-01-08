      subroutine absorb (ngrid, peout, piout)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
c
      complex*16      skxx, skxy, skxz, skyy, skyz, skzz
      complex*16      efld, dielt, xkie, ey, eytot,
     .                pdsx, pdsy, pdsz, plex
      real*8          peout(kich), piout(kich), pds(kich,ks),
     .                pdsm(kich,ks,4)
c
      common /diel  / skxx(kich,ks,4),skxy(kich,ks,4),skxz(kich,ks,4),
     .                skyy(kich,ks,4),skyz(kich,ks,4),skzz(kich,ks,4)
      common /rfp   / efld(kich,3,4),dielt(6,2),exy2(kich,2),
     .                ezy2(kich,2),elr2(kich,2)
      common /ichflx/ xkie(kich,4),ey(kich,4),eytot(kich)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      do i=1,ngrid
      do j=1,neq
      pds(i,j) = 0.0
      do k=1,4
        pdsx = CMPLX (0.0, 0.0)
        pdsy = CMPLX (0.0, 0.0)
        pdsz = CMPLX (0.0, 0.0)
        do kk=1,4
          pdsx = pdsx
     .         + ey(i,kk) * (PIMAG (skxx(i,j,kk))*efld(i,1,kk)
     .         - CMPLX (0.0, REAL (skxy(i,j,kk)))
     .         + PIMAG (skxz(i,j,kk))*efld(i,3,kk))
          pdsy = pdsy
     .         + ey(i,kk)*(CMPLX (0.0, REAL (skxy(i,j,kk)))*efld(i,1,kk)
     .         + PIMAG (skyy(i,j,kk))
     .         + CMPLX (0.0, -REAL (skyz(i,j,kk)))*efld(i,3,kk))
          pdsz = pdsz
     .         + ey(i,kk)*(PIMAG (skxz(i,j,kk))*efld(i,1,kk)
     .         - CMPLX (0.0, -REAL (skyz(i,j,kk)))
     .         + PIMAG (skzz(i,j,kk))*efld(i,3,kk))
        end do
        plex = conjg(ey(i,k)*efld(i,1,k))*pdsx+conjg(ey(i,k))*pdsy
     .       + conjg(ey(i,k)*efld(i,3,k))*pdsz
        pdsm(i,j,k) = REAL (plex)
        pds(i,j) = pds(i,j)+pdsm(i,j,k)
      end do
      end do
      peout(i) = pds(i,1)
      piout(i) = 0.0
      do m=2,neq
        piout(i) = piout(i) + pds(i,m)
      end do
      end do
      return
c
      end

      subroutine akim (fpi, rhoi, nsp, nslb, vi, li, zims, frq, wkpl,
     .                 gam, kxr2, li0)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c gam(nsp,i,nslb),i = 1,2,3--1:(1,1);2:(2,2),3:(1,2)--are the imaginary
c parts of the ion polarization tensor for each species
c
      include 'ichp2.i'
c
      parameter (nx = 2, nz = 3, lgam = 3)
      real*8     kxr2
      dimension  fpi(nx,ny),rhoi(nx,ny),vi(nx,ny),li(nx),
     .           zims(nx,ny,nz),gam(nx,lgam,ny),work(10),kxr2(ny)
c
      data pi/3.14159265/
c
      rpi = SQRT (pi)
      vph = frq / wkpl
      do i=1,nslb
        do n=1,nsp
          fac1 = rpi*vph*(fpi(n,i)/frq)**2
          fac1 = fac1/vi(n,i)
          fac2 = 0.5*kxr2(i)*rhoi(n,i)**2
          call modbes (fac2, li(n)+2, work)
          do j=1,3
            gam(n,j,i) = 0.0
          end do
          do l=1,li(n)
            l1   = l+1
            l2   = l+2
            fac4 = zims(n,i,l)**2
            if (fac4 .le. 100.0) then
              fac3 = EXP (-fac4)
            else
              fac3 = 0.0
            end if
            gam(n,1,i) = gam(n,1,i)+fac1*fac3* l**2*work(l1)/fac2
            gam(n,2,i) = gam(n,2,i)+fac1*fac3*(l**2*work(l1)/fac2
     .                 + 2.0*fac2*(work(l1)-0.5*work(l)-0.5*work(l2)))
            gam(n,3,i) = gam(n,3,i)+fac1*fac3*l*(work(l1)
     .                 - work(l)*0.5-work(l2)*0.5)
          end do
        end do
      end do
      return
c
      end

      subroutine asyr (x, zr, dzr, n, rem, tol)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      rem1 = 1.0
      s    = 1.0
      dzr  = 0.0
      b    = 1.0 / x
c
      do 10 n=1,40
        x1  = FLOAT (n)
        p   = b**2*(2.0*x1-1.0)*0.5
        if (n .eq. 1)  go to 11
        p1  = p1*p
        go to 12
   11   p1  = p
   12   rem = ABS (p1) / ABS (s)
        if (rem  .lt. tol)  go to 13
        if (rem1 .lt. rem)  go to 13
        rem1 = rem
        dzr  = dzr + p1
   10 s = s + p1
c
   13 zr  =  -b * s
      dzr = 2.0 * (dzr+p1)
      return
c
      end

      subroutine beslci (x, y, nb, ize, br, bi, ncalc)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     Journal of Research   National Bureau of Standards
c     b-mathematical sciences   vol.77-b  nos.3-4  july-dec. 1973  p118
c     routine to calculate Bessel functions j and i
c     of complex argument and integer order
c                  input variables
c     x real      part of complex argument
c     y imaginary part of complex argument
c     nb a positive integer designating highest order to be calculated
c     ize     zero for j*s    one for i*s
c     br  for normal exit, br contains real      part of solution vector
c     bi  for normal exit, bi contains imaginary part of solution vector
c     normal exit if  ncalc = nb
c     br-bi-ncalc need not be initialized
c                  machine dependent constants
c     nsig  decimal significance desired. set to INT (LOG10 (2)*nbit+1)
c     nbit  number of bits in the mantissa
c     the relative truncation error is limited to t = 0.5*10**-nsig
c     for order greater than ABS (z).
c     for order less    than ABS (z) (general test), the relative error is
c     limited to t for function values of magnitude at least 1.
c     the absolute error is limited to t for smaller values.
c     nten  largest integer k such that 10**k is machine representable.
c     largez  upper limit on the magnitude of z. if ABS (z) = n, at least
c     n iterations of the backward recursion will be executed.
c     exparg  largest argument that the library exp routine can handle.
c                  error returns
c     let g denote either i or j.
c     in case of an error, ncalc .ne. nb and not all g*s are calculated
c     to the desired accuracy.
c     if ncalc .lt. 0, an argument is out of range. nb .le. 0 or ize is
c     neither 0 nor 1 or ize = 0 and ABS (y) > exparg, or ize = 1 and
c     ABS (x) > exparg.  in this case,the vectors br and bi are not
c     calculated, and ncalc is set to MIN0 (nb,0)-1 so ncalc .ne. nb.
c     nb .gt. ncalc .gt. 0 will occur if nb .gt. magz and ABS (g-sub-nb-of-z/g
c     -sub-magx-of-z) .lt. 10.0**(nten/2), i.e. nb is much greater than
c     magz. in this case, br(n) and bi(n) are calculated to the desired
c     accuracy for n .le. ncalc, but for ncalc .lt. n .le. nb, precision  is
c     lost. if n .gt. ncalc and ABS (g(ncalc-1)/g(n-1)) .eq. 10**-k, then
c     the last k significant figures of g(n-1) ( = br(n)+i*bi(n)) are
c     erroneous. if the user wishes to calculate g(n-1) to higher
c     accuracy, he should use an asymptotic formula for large order.
c
      real*4     singloid ! single precision variable
      dimension  br(*), bi(*)
      data       nsig, nten, largez, exparg
     .          /  15,  307,  10000,  700.0/
c
      tempar   = SQRT (x*x + y*y)
      singloid =       tempar
      magz     = INT (singloid)
      if (nb .gt. 0 .and.   magz  .le. largez  .and.
     . ((ize .eq. 0 .and. ABS (y) .le. exparg)  .or.
     .  (ize .eq. 1 .and. ABS (x) .le. exparg)))  go to 1
c
c     error return    z,nb,or ize is out of range
c
      ncalc = MIN0 (nb, 0) - 1
      return
c
    1 sign  = 1 - 2 * ize
      ncalc = nb
c
c     use 2-term ascending series for small z
c
      if (tempar**4 .lt. 0.1**nsig)  go to 50
c
c     initialize the calculation of the p*s
c
      nbmz   = nb-magz
      n      = magz+1
      if (ABS (x) .lt. ABS (y))  go to 2
      zinvr  = 1.0/(x+y*y/x)
      zinvi  = -y*zinvr/x
      go to 3
    2 zinvi  = -1.0/(y+x*x/y)
      zinvr  = -x*zinvi/y
    3 plastr = 1.0
      plasti = 0.0
      pr     = sign*(n+n)*zinvr
      pi     = sign*(n+n)*zinvi
      test   = 2.0 * (10.0)**nsig
      m      = 0
      if (nbmz .lt. 3)  go to 6
c
c     calculate p*s until n = nb-1.  check for possible overflow.
c
      tover  = 10.0**(nten-nsig)
      nstart = magz+2
      nend   = nb - 1
      do 5 n=nstart,nend
        poldr  = plastr
        poldi  = plasti
        plastr = pr
        plasti = pi
        pr     = sign*((n+n)*(plastr*zinvr-plasti*zinvi)-poldr)
        pi     = sign*((n+n)*(plasti*zinvr+plastr*zinvi)-poldi)
        if ((pr/tover)**2+(pi/tover)**2-1.0)  5, 5, 7
    5 continue
      n = nend
c
c     calculate special significance test for nbmz .gt. 2.
c
      tempbi = MAX (ABS (pr), ABS (pi))
      tempbi = tempbi * SQRT (2.0 * (10.0)**nsig * SQRT (((pr/tempbi)**2
     .       + (pi / tempbi)**2)*((plastr/tempbi)**2
     .                          + (plasti/tempbi)**2)))
      test   = MAX (test,tempbi)
c
c     calculate p*s until significance test is passed
c
    6 n      = n + 1
      poldr  = plastr
      poldi  = plasti
      plastr = pr
      plasti = pi
      pr     = sign*((n+n)*(plastr*zinvr-plasti*zinvi)-poldr)
      pi     = sign*((n+n)*(plasti*zinvr+plastr*zinvi)-poldi)
      if ((pr/test)**2+(pi/test)**2 .lt. 1.0)  go to 6
      if (m .eq. 1)  go to 12
c
c     calculate strict variant of significance test,
c     and calculate p*s until this test is passed
c
      m      = 1
      tempbi = MAX (ABS (pr), ABS (pi))
      tempbr = SQRT (((pr/tempbi)**2+(pi/tempbi)**2)/
     .              ((plastr/tempbi)**2+(plasti/tempbi)**2))
      tempbi = (n+1) / tempar
      if (tempbr+1.0/tempbr .gt. 2.0*tempbi)
     .    tempbr = tempbi + SQRT (tempbi**2-1.0)
      test       = test / SQRT (tempbr-1.0/tempbr)
      if ((pr/test)**2+(pi/test)**2-1.0) 6,12,12
    7 nstart     = n + 1
c
c     to avoid overflow, normalize p*s by dividing by tover.
c     calculate p*s until unnormalized p would overflow.
c
      pr     = pr/tover
      pi     = pi/tover
      plastr = plastr/tover
      plasti = plasti/tover
      psaver = pr
      psavei = pi
      tempcr = plastr
      tempci = plasti
      test   = 10.0**(2*nsig)
    8 n      = n+1
      poldr  = plastr
      poldi  = plasti
      plastr = pr
      plasti = pi
      pr     = sign*((n+n)*(plastr*zinvr-plasti*zinvi)-poldr)
      pi     = sign*((n+n)*(plasti*zinvr+plastr*zinvi)-poldi)
      if (pr**2+pi**2 .le. test)  go to 8
c
c     calculate backward test,and find ncalc,the highest n
c     such that the test is passed.
c
      tempbr = SQRT ((plastr**2+plasti**2)/(poldr**2+poldi**2))
      tempbi = n / tempar
      if (tempbr+1.0/tempbr .gt. 2.0*tempbi)
     .    tempbr = tempbi + SQRT (tempbi**2-1.0)
      test = 0.5*(1.0-1.0/tempbr**2)/10.0**nsig
      test = ((plastr**2+plasti**2)*test)*((poldr**2+poldi**2)*test)
      pr   = plastr*tover
      pi   = plasti*tover
      n    = n-1
      nend = MIN0 (nb,n)
      do 9 ncalc=nstart,nend
      poldr  = tempcr
      poldi  = tempci
      tempcr = psaver
      tempci = psavei
      psaver = sign*((n+n)*(tempcr*zinvr-tempci*zinvi)-poldr)
      psavei = sign*((n+n)*(tempci*zinvr+tempcr*zinvi)-poldi)
      if ((psaver**2+psavei**2)*(tempcr**2+tempci**2)-test)  9, 9, 10
    9 continue
      ncalc  = nend+1
   10 ncalc  = ncalc-1
c
c     the coefficient of b(n) in the normalized sum is
c     m * SQRT (-1)**imag,where m = -2,0,or2,and imag is 0 or 1.
c     calculate recursion rules for m and imag,and initialize them.
c
   12 n      = n+1
      tempbr = ize*x+(1-ize)*y
      ipos   = 0
      if (tempbr)  13, 14, 13
   13 singloid = 1.1 * tempbr / ABS (tempbr)
      ipos   = INT (singloid)
   14 mrecur = 4*((2+ize+ipos)/2)-3-2*(ize+ipos)
      k      = 2+ipos+2*ize*ipos**2-ize
      l      = n-4*(n/4)
      mlast  = 2+8*((k*l)/4)-4*((k*l)/2)
      if (ipos .eq. 0 .and. (l .eq. 1 .or. l .eq. 3)) mlast = 0
      l      = l+3-4*((l+3)/4)
      m      = 2+8*((k*l)/4)-4*((k*l)/2)
      if (ipos .eq. 0 .and. (l .eq. 1 .or. l .eq. 3)) m = 0
      imrecr = (1-ize)*ipos**2
      imag   = imrecr*(l-2*(l/2))
c
c     initialize the backward recursion and the normalization sum
c
      tempbr = 0.0
      tempbi = 0.0
      if (ABS (pi) .gt. ABS (pr))  go to 15
      tempar = 1.0 / (pr+pi*(pi/pr))
      tempai = -(pi*tempar)/pr
      go to 16
   15 tempai = -1.0/(pi+pr*(pr/pi))
      tempar = -(pr*tempai)/pi
   16 if (imag .ne. 0)  go to 17
      sumr   = m*tempar
      sumi   = m*tempai
      go to 18
   17 sumr   = -m * tempai
      sumi   =  m * tempar
   18 nend   =  n - nb
      if (nend) 26, 22, 19
c
c     recur backward via difference equation calculating (but not
c     storing) br(n) and bi(n) until n = nb
c
   19 do 21 l=1,nend
        n      = n-1
        tempcr = tempbr
        tempci = tempbi
        tempbr = tempar
        tempbi = tempai
        pr     = (n+n)*zinvr
        pi     = (n+n)*zinvi
        tempar = pr*tempbr-pi*tempbi-sign*tempcr
        tempai = pr*tempbi+pi*tempbr-sign*tempci
        imag   = (1-imag)*imrecr
        k      = mlast
        mlast  = m
        m      = k*mrecur
        if (imag .ne. 0)  go to 20
        sumr   = sumr+m*tempar
        sumi   = sumi+m*tempai
        go to 21
   20   sumr   = sumr-m*tempai
        sumi   = sumi+m*tempar
   21 continue
c
c     store  br(nb),bi(nb)
c
   22 br(n) = tempar
      bi(n) = tempai
      if (n .gt. 1)  go to 23
c
c     nb = 1.  since 2*tempar and 2*tempai were added to sumr and sumi
c     respectively,we must subaxis(i)ract tempar and tempai
c
      sumr = sumr-tempar
      sumi = sumi-tempai
      go to 35
c
c     calculate and store br(nb-1),bi(nb-1)
c
   23 n     = n-1
      pr    = (n+n)*zinvr
      pi    = (n+n)*zinvi
      br(n) = pr*tempar-pi*tempai-sign*tempbr
      bi(n) = pr*tempai+pi*tempar-sign*tempbi
      if (n .eq. 1)  go to 34
      imag  = (1-imag)*imrecr
      k     = mlast
      mlast = m
      m     = k*mrecur
      if (imag .ne. 0)  go to 24
      sumr  = sumr+m*br(n)
      sumi  = sumi+m*bi(n)
      go to 30
   24 sumr  = sumr-m*bi(n)
      sumi  = sumi+m*br(n)
      go to 30
c
c     n < nb, so store br(n),bi(n) and set higher orders zero
c
   26 br(n)   = tempar
      bi(n)   = tempai
      nend    = -nend
      do 27 l=1,nend
      br(n+l) = 0.0
   27 bi(n+l) = 0.0
   30 nend    = n-2
      if (nend .eq. 0)  go to 33
c
c     calculate via difference equation and store br(n),bi(n) until n = 2
c
      do 32 l=1,nend
      n     = n-1
      pr    = (n+n)*zinvr
      pi    = (n+n)*zinvi
      br(n) = pr*br(n+1)-pi*bi(n+1)-sign*br(n+2)
      bi(n) = pr*bi(n+1)+pi*br(n+1)-sign*bi(n+2)
      imag  = (1-imag)*imrecr
      k     = mlast
      mlast = m
      m     = k*mrecur
      if (imag .ne. 0)  go to 31
      sumr  = sumr+m*br(n)
      sumi  = sumi+m*bi(n)
      go to 32
   31 sumr  = sumr-m*bi(n)
      sumi  = sumi+m*br(n)
   32 continue
c
c     calculate and store br(1),bi(1)
c
   33 br(1) = 2.0 * (br(2)*zinvr-bi(2)*zinvi)-sign*br(3)
      bi(1) = 2.0 * (br(2)*zinvi+bi(2)*zinvr)-sign*bi(3)
   34 sumr  = sumr+br(1)
      sumi  = sumi+bi(1)
c
c     calculate normalization factor. tempar+i*tempai
c
   35 if (ize .eq. 1)  go to 36
      tempcr = ipos*y
      tempci = -ipos*x
      go to 37
   36 tempcr = ipos*x
      tempci = ipos*y
   37 tempcr = EXP (tempcr)
      tempbr = COS (tempci)
      tempbi = SIN (tempci)
      if (ABS (sumr) .lt. ABS (sumi))  go to 38
      tempci = sumi/sumr
      tempcr = (tempcr/sumr)/(1.0+tempci*tempci)
      tempar = tempcr*(tempbr+tempbi*tempci)
      tempai = tempcr*(tempbi-tempbr*tempci)
      go to 39
   38 tempci = sumr/sumi
      tempcr = (tempcr/sumi)/(1.0+tempci*tempci)
      tempar = tempcr*(tempbr*tempci+tempbi)
      tempai = tempcr*(tempbi*tempci-tempbr)
c
c     normalize
c
   39 do 40 n=1,nb
      tempbr = br(n)*tempar-bi(n)*tempai
      bi(n)  = br(n)*tempai+bi(n)*tempar
   40 br(n)  = tempbr
      return
c
c     two-term ascending series for small z
c
   50 tempar = 1.0
      tempai = 0.0
      tempcr = 0.25*(x*x-y*y)
      tempci = 0.5*x*y
      br(1)  = 1.0-sign*tempcr
      bi(1)  = -sign*tempci
      if (nb .eq. 1)  go to 52
      do 51 n=2,nb
      tempbr = (tempar*x-tempai*y)/(n+n-2)
      tempai = (tempar*y+tempai*x)/(n+n-2)
      tempar = tempbr
      tempbr = n
      br(n)  = tempar*(1.0-sign*tempcr/tempbr)+tempai*tempci/tempbr
   51 bi(n)  = tempai*(1.0-sign*tempcr/tempbr)-tempar*tempci/tempbr
   52 return
c
      end


      subroutine call_genray ( model, qrfe, qrfi, currf,
     .                         totgrpe, totgrpi, totgrc, totgrp)
c
c******   
c******   Text file interface to GENRAY  (Bob Harvey, Jan 24, 2005)
c******   Added netCDF interface to GENRAY (BH, Sept 30, 2005)
c
c
c ----------------------------------------------------------------------
c
c  This is an interface routine between ONETWO and the all frequencies
c  ray tracing code GENRAY.
c  [This interface routine is modified from subroutine raytrace, above.]
c
c  The broad sequence of operations is:
c  (1) One file, genray_profs_in.txt will be written as a text
c      and netCDF file (genray_profs_in.nc)
c      file by ONETWO to pass plasma profile data to GENRAY.
c      [Later, the text file will be dropped.]
c      A separate namelist file, genray.dat gives launch
c      parameters for the rays, and must be specified by the user.
c  (2) The ray tracing code will be run, producing a text or netCDF
c      output file, genray_profs_out/genray_profs_out.nc, 
c      containing power deposition and current drive data,
c      [later, the text file will be dropped]  to be
c      read in by this subroutine
c
c  Output variables are:
c     qrfe (j), j=1,nj    power density to electrons (W/cm**3)
c     qrfi (j), j=1,nj    power to (all) ion species (W/cm**3)
c     currf(j), j=1,nj    RF current density <j.B/B_0> (A/cm**2)
c     totgrpe, totgrpi, totgrc,  totgrp, are respectively total power to
c       electrons (W), to ions (W), and total toroidal RF current (Amps),
c       input power to genray (W) for this mode.
c ----------------------------------------------------------------------
c
c     Variables in rfmod.f90, MODULE rf
      USE rf,only : rgenray,pgre,pgrit,pgrc,pgri,totgrps,
     .              save_genray_io,genray_fi,nspecgr,
                    !nspecgr is the number of ion species passed to genray,
                    !including possible hot beam species and 'dt' species
                    !(bulk and beam) which have been split into two.
     .              charge_nc,dmass_nc,en_nc,temp_nc,
     .              genraydat,gfilename,nfwsmth,nicsmth,necsmth,
     .              char1,char2,char_name1,char_name2,char_name
!BH: Folowing statement is in numbrs.f90      USE param,only : kj
c      USE param  !Commented, since parameters are in modules
      USE ext_prog_info, only : get_genray
      USE ename ,only : eqdskfilename,eqfile,eqdsk_tdem
      USE io, only : eqdskin,ncrt,nout  !in io.f90
      USE machin, only : rmajor
      USE mesh, only : r
      USE geom, only : codeid,hcap
      USE soln, only : ene,en,te,ti
      USE solcon, only : time,irfcnt
      USE ions, only : namep,namei,nameb,atw,atomno,zeff
      USE nub2, only : enbeam,wbeam,tbeam
      USE numbrs, only : nprim,nimp,nion,nj !in numbrs.f90
      USE fusion, only : fd,fdbeam
      implicit  integer (i-n), real*8 (a-h, o-z)

      character(len=256)      :: eqdsk_name
      character(len=256)      :: genraydatt
      character(len=256),save :: genray_path_out  !From ext_prog_info
      character(len=256)      :: command
      dimension qrfe(nj),qrfi(nj),currf(nj)
      integer istat

c --- include file for netCDF declarations 
c --- (obtained from NetCDF distribution)
      include 'netcdf.inc'

c --- some stuff for netCDF file ---
      integer ncid,vid,istatus
      character filenc*128,ltitle*128,name*128
      integer dims_id(2),start(2),count(2),countm(2),char256dim
      integer nbulk,nbulkm,nprofs
      integer nj_id,nspecgr_id
      integer nprofs_id,nbulk_id,nbulkm_id



      logical,save ::  first_time 
      data first_time / .true./
      save :: pemassrat
      data pemassrat /1836.2/    !proton-to-electron mass ratio
c
      data start/1,1/

c     
c     Check whether there is a beam ion component
c     [Beam is time-dependent, and enbeam 0'ed at each entry to source.]
      ibeamcmpt=0
      do j=1,nj
         if (enbeam(j).gt.0.) ibeamcmpt=1
      enddo

      nspecgr=1+nprim+nimp+ibeamcmpt  !electrons, ions, & beam component

c     If one of the primary ion species is 'dt', then split it up into
c     two species
      ipsplit=0
      do i=1,nprim
         if (namep(i).eq.'dt') ipsplit=ipsplit+1
      enddo
      if (ipsplit.gt.1) STOP 'CALL_GENRAY not configrd for two dt ions'
      nspecgr=nspecgr+ipsplit

      nion1=nion+ipsplit   !total ion species, except beam components.

c     If the beam species is 'dt', then split it up into two species
      ibsplit=0
c      if (nameb.eq.'dt') ibsplit=1 HSj changed  to following 8/26/11
      if (nameb.eq.'dt' .AND. ibeamcmpt == 1 ) ibsplit=1
      nspecgr=nspecgr+ibsplit   !total species, including any split ion
                                !'dt' species.
c      print *,'ndt beam  nspecgr,ibsplit =', nspecgr,ibsplit ! 888889999
c     The beam component must be the first or second primary species,
c     and is restricted to 'h', 'd', 't', 'dt'.
c     If ibeamcmpt.ne.0, give it the primary species component number.
      if (ibeamcmpt.ne.0) then
         do i=1,nprim
            if (nameb.eq.namep(i)) ibeamcmpt=i
         enddo
      endif
c     

c     Write out number of plasma species, dimension of the
c     equispaced radial mesh, radial mesh, and zeff. 
c     Then for each species,
c     write out charge number (pos for electron and ion), 
c     atomic weight in units of electron mass, 
c     then density (/cm**3) and temperature (keV). Start with electrons.
c
c     NOTE:  Charge of ions is taken to be fully stripped value,
c            that is, the atomic number.  A more careful treatment
c            would account for the ionization state of Z>1 atoms
c            as a function of plasma radius.  A volume average
c            charge state might be sufficient for most purposes.
c            (BH: 060202)

c     Also, allocate some storage to facilitate netCDF write:
      if(.not. allocated(charge_nc))then
         allocate (charge_nc(nspecgr),STAT = istat)
         if(istat .ne. 0)
     .        call allocate_error("charge_nc,sub call_genray",0,istat)
      endif
      if(.not. allocated(dmass_nc))then
         allocate (dmass_nc(nspecgr),STAT = istat)
         if(istat .ne. 0)
     .        call allocate_error("dmass_nc,sub call_genray",0,istat)
      endif
      if(.not. allocated(en_nc))then
         allocate (en_nc(nj,nspecgr),STAT = istat)
         if(istat .ne. 0)
     .        call allocate_error("en_nc,sub call_genray",0,istat)
      endif
      if(.not. allocated(temp_nc))then
         allocate (temp_nc(nj,nspecgr),STAT = istat)
         if(istat .ne. 0)
     .        call allocate_error("temp_nc,sub call_genray",0,istat)
      endif




c
c     ---    Create genray_profs_in.txt
c
      call getioun(igenray,43)  ! local unit no. for files
c     
      open (unit=igenray,file='genray_profs_in.txt',status='UNKNOWN')
c     Write itendification header
      write(igenray,*)'Written by ONETWO: In subroutine call_genray'

      write (igenray, 1003)  nspecgr,nj
      write (igenray, 1001)  (r(j), j=1,nj)     !cms
      write (igenray, 1001)  (zeff(j), j=1,nj)
      one=1.
      write (igenray, 1001) one,one
      write (igenray, 1001)  (ene(j), j=1,nj)   !/cm**3
      write (igenray, 1001)  (te(j), j=1,nj)   ! keV
        ispecgr=1
        charge_nc(ispecgr)=one
        dmass_nc(ispecgr)=one
        do j=1,nj
           en_nc(j,ispecgr)=ene(j)
           temp_nc(j,ispecgr)=te(j)
        enddo

      do i=1,nprim
         ispecgr=ispecgr+1
         if(namep(i).eq.'dt') then

            write (igenray, 1001) one,2.*pemassrat
            write (igenray, 1001)  (fd*en(j,i), j=1,nj) !/cm**3
            write (igenray, 1001)  (ti(j), j=1,nj) ! keV

            charge_nc(ispecgr)=one
            dmass_nc(ispecgr)=2.*pemassrat
            do j=1,nj
               en_nc(j,ispecgr)=fd*en(j,i)
               temp_nc(j,ispecgr)=ti(j)
            enddo

            write (igenray, 1001) one,3.*pemassrat
            write (igenray, 1001)  ((1.-fd)*en(j,i), j=1,nj) !/cm**3
            write (igenray, 1001)  (ti(j), j=1,nj) ! keV

            ispecgr=ispecgr+1
            charge_nc(ispecgr)=one
            dmass_nc(ispecgr)=3.*pemassrat
            do j=1,nj
               en_nc(j,ispecgr)=(1.-fd)*en(j,i)
               temp_nc(j,ispecgr)=ti(j)
            enddo

         else

            write (igenray, 1001) atomno(i),atw(i)*pemassrat
            write (igenray, 1001)  (en(j,i), j=1,nj) !/cm**3
            write (igenray, 1001)  (ti(j), j=1,nj) ! keV

            charge_nc(ispecgr)=atomno(i)
            dmass_nc(ispecgr)=atw(i)*pemassrat
            do j=1,nj
               en_nc(j,ispecgr)=en(j,i)
               temp_nc(j,ispecgr)=ti(j)
            enddo
         endif
      enddo
      
      do i=nprim+1,nprim+nimp
         ispecgr=ispecgr+1

         write (igenray, 1001) atomno(i),atw(i)*pemassrat
         write (igenray, 1001)  (en(j,i), j=1,nj) !/cm**3
         write (igenray, 1001)  (ti(j), j=1,nj) ! keV

         charge_nc(ispecgr)=atomno(i)
         dmass_nc(ispecgr)=atw(i)*pemassrat
         do j=1,nj
            en_nc(j,ispecgr)=en(j,i)
            temp_nc(j,ispecgr)=ti(j)
         enddo
      enddo

      if (ibeamcmpt.ne.0) then
c        Beam temperature taken to be (2./3)[avg beam energy (keV/cm**3)]/
c        [beam density (/cm**3)]:
         do j=1,nj
            tbeam(j)=(2./3.)*wbeam(j)/enbeam(j)
         enddo

         if(nameb.eq.'dt') then
            write (igenray, 1001) one,2.*pemassrat
            write (igenray, 1001)  (fd*enbeam(j), j=1,nj) !/cm**3
            write (igenray, 1001)  (tbeam(j), j=1,nj) ! keV

            ispecgr=ispecgr+1
            charge_nc(ispecgr)=one
            dmass_nc(ispecgr)=2.*pemassrat
            do j=1,nj
               en_nc(j,ispecgr)=fd*enbeam(j)
               temp_nc(j,ispecgr)=tbeam(j)
            enddo
              
            write (igenray, 1001) one,3.*pemassrat
            write (igenray, 1001)  ((1.-fd)*enbeam(j), j=1,nj) !/cm**3
            write (igenray, 1001)  (tbeam(j), j=1,nj) ! keV

            ispecgr=ispecgr+1
            charge_nc(ispecgr)=one
            dmass_nc(ispecgr)=3.*pemassrat
            do j=1,nj
               en_nc(j,ispecgr)=(1.-fd)*enbeam(j)
               temp_nc(j,ispecgr)=tbeam(j)
            enddo
        else
            write (igenray, 1001) atomno(ibeamcmpt),
     .                            atw(ibeamcmpt)*pemassrat
            write (igenray, 1001)  (enbeam(j), j=1,nj) !/cm**3
            write (igenray, 1001)  (tbeam(j), j=1,nj) ! keV
            
            ispecgr=ispecgr+1
            charge_nc(ispecgr)=atomno(ibeamcmpt)
            dmass_nc(ispecgr)=atw(ibeamcmpt)*pemassrat
            do j=1,nj
               en_nc(j,ispecgr)=enbeam(j)
               temp_nc(j,ispecgr)=tbeam(j)
            enddo
         endif
      endif
         
 1000 format(10a8)
 1001 format (5(1pe16.9))
 1003 format (10i5)
c     
c     Close genray_profs_in.txt
      call giveupus(igenray)
      close (unit = igenray)

c     Get eqdsk_name and create eqdsk, as needed.
c     eqdsk_name will be passed to genray.

      if (codeid .eq. 'onedee') then
         call giveupus(igenray)
         close (unit =  igenray)
         
         write(*,*)'call_genray: GENRAY not set up for 1D'
         STOP
c     
      else
c
c     2D equilibrium case [following coding in subroutine wrt_curray_in]
c     Add blank at end of file name, for detection when read by NetCDF.

         eqfile = eqdskfilename
         if(eqdsk_tdem .ne. 'tdem' ) then
            lenge = LEN_TRIM(eqdskfilename)
            eqdsk_name  = eqdskfilename(1:lenge)//' '
         else
            call wrt_tdem_eqdsk(time,eqdsk_name)
            eqdsk_name=eqdsk_name//' '
            eqfile = eqdsk_name
            print *,'eqdsk file created in sub wrt_tdem_eqdsk'
         endif
         print *,'call_genray: eqdsk_name =', eqdsk_name
         print *,'call_genray: eqfile =',eqfile
         print *,'call_genray: eqdsk_tdem =',eqdsk_tdem
         if (eqfile .eq. 'none')  eqdsk_name =  'eqdskin '
 
      endif


c-----------------------------------------------------------------------
c     Set up and write data to netcdf file
c-----------------------------------------------------------------------

c
c     Open netCDF file for output
      filenc='genray_profs_in.nc'
      ncid=nccre(filenc,NCCLOB,istatus)
      call check_err(istatus)
      write(*,*)'In call_genray after nccre, istatus=',istatus
c     Brief description to be added to file:
      ltitle='Profile data passed from ONETWO to GENRAY'

c
      call ncaptc(ncid,NCGLOBAL,'title',NCCHAR,length_char(ltitle),
     +     ltitle,istatus)


      count(1)=nj
      count(2)=nspecgr

cl    Define dimensions
      nj_id=ncddef(ncid,'nj_dim',nj,istatus)
      nspecgr_id=ncddef(ncid,'nspecgr_dim',nspecgr,istatus)
      char256dim=ncddef(ncid,'char256dim',256,istatus)

c     Define vector of dimensions
      dims_id(1)=nj_id
      dims_id(2)=nspecgr_id

c     Define variable names for netcdf file

      vid=ncvdef(ncid,'eqdsk_name',NCCHAR,1,char256dim,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,39,
     +          'Name of input eqdsk; may vary with time',istatus)

      vid=ncvdef(ncid,'nj',NCLONG,0,0,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,32,
     +           'Transport radial mesh dimension',istatus)
      call check_err(istatus)

      vid=ncvdef(ncid,'nspecgr',NCLONG,0,0,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,40,
     +           'Number of plasma species, incl electrons',istatus)
      call check_err(istatus)

      vid=ncvdef(ncid,'r',NCDOUBLE,1,nj_id,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,21,
     +           'transport radial mesh',istatus)
      call ncaptc(ncid,vid,'units',NCCHAR,3,
     +           'cms',istatus)
      call check_err(istatus)

      vid=ncvdef(ncid,'zeff',NCDOUBLE,1,nj_id,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,29,
     +           'zeff on transport radial mesh',istatus)
      call check_err(istatus)

      vid=ncvdef(ncid,'charge',NCDOUBLE,1,nspecgr_id,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,40,
     +           'Charge number of e (pos) and ion species',istatus)
      call check_err(istatus)

      vid=ncvdef(ncid,'dmass',NCDOUBLE,1,nspecgr_id,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,42,
     +           'Atomic wt of each spec, normd to elec mass',istatus)
      call check_err(istatus)

      vid=ncvdef(ncid,'en',NCDOUBLE,2,dims_id,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,36,
     +           'Electron (first) and ion densities',istatus)
      call ncaptc(ncid,vid,'units',NCCHAR,6,
     +           '/cm**3',istatus)
      call check_err(istatus)
 
      vid=ncvdef(ncid,'temp',NCDOUBLE,2,dims_id,istatus)
      call ncaptc(ncid,vid,'long_name',NCCHAR,43,
     +           'Electron (first) and ion temperatures',istatus)
      call ncaptc(ncid,vid,'units',NCCHAR,3,
     +           'keV',istatus)
      call check_err(istatus)


c     End the define-mode and start the data-mode

      call ncendf(ncid,istatus)
      call check_err(istatus)

c     Write data into netcdf file
c$$$      write (igenray, 1003)  nspecgr,nj
c$$$      write (igenray, 1001)  (r(j), j=1,nj)     !cms
c$$$      write (igenray, 1001)  (zeff(j), j=1,nj)


      vid=ncvid(ncid,'eqdsk_name',istatus)
c     Add 1 to get blank into variable [See NetCDF manual].
      ll=length_char(eqdsk_name)+1
      if (ll.gt.256) then
         write(*,*)'call_genray: eqdsk_name length too great'
         stop
      endif
      call ncvptc(ncid,vid,1,ll,eqdsk_name,ll,istatus)
      call check_err(istatus)

      vid=ncvid(ncid,'nj',istatus)
      call ncvpt(ncid,vid,1,1,nj,istatus)
      call check_err(istatus)

      vid=ncvid(ncid,'nspecgr',istatus)
      call ncvpt(ncid,vid,1,1,nspecgr,istatus)

      vid=ncvid(ncid,'r',istatus)
      call ncvpt(ncid,vid,1,count(1),r(1),istatus)

      vid=ncvid(ncid,'zeff',istatus)
      call ncvpt(ncid,vid,1,count(1),zeff(1),istatus)

      vid=ncvid(ncid,'charge',istatus)
      call ncvpt(ncid,vid,1,count(2),charge_nc,istatus)

      vid=ncvid(ncid,'dmass',istatus)
      call ncvpt(ncid,vid,1,count(2),dmass_nc,istatus)

      vid=ncvid(ncid,'en',istatus)
      call ncvpt(ncid,vid,start,count,en_nc,istatus)

      vid=ncvid(ncid,'temp',istatus)
      call ncvpt(ncid,vid,start,count,temp_nc,istatus)
      call check_err(istatus)
    
      call ncclos(ncid,istatus)
      call check_err(istatus)

c
c     Copy full-path-specified genray input file to genray.dat.
c     (default is genray.dat in present directory.)
      genraydatt=genraydat(model)
      command = 'cp '//genraydatt(1:LEN_TRIM(genraydatt))//
     .               ' '//'genray.dat'
      write(*,*)'cray309: command= ',command
      if ( genraydatt(1:LEN_TRIM(genraydatt)) .ne. 'genray.dat') then
           if (ISHELL (command) .lt. 0)
     .     call STOP ('sub source: failure of spawned cp  command', 67)
      endif	


c     SAVE COPY OF INPUT FILES (save_genray_io .eq. 1)

      if (save_genray_io .eq. 1) then

c        Construct ascii 'model.irfcnt' designator for files
c        (ASCII character 48 is a '0')
         if (model.gt.99) then
            write(*,*) 'source:  STOP in model char_name creation'
            STOP
         endif
c        Set up char1 for 1st digit in 10:99, char2 for second digit
         char1=achar(model/10 + 48)
         char2=achar(mod(model,10) + 48)
         write(*,*)'source:  model,char1,char2 = ',model,char1,char2
         char_name1=char1//char2

c        Similary for irfcnt(model):
c        Since irfcnt() not incremented in sub source yet, add 1:
         irfcntt=irfcnt(model)+1
         if (irfcntt.gt.99) then
            write(*,*) 'source:  STOP in irfcnt char_name creation'
            write(*,*) 'source:  Add a little more coding...'
            STOP
         endif
         char1=achar(irfcntt/10 + 48)
         char2=achar(mod(irfcntt,10) + 48)
         char_name2=char1//char2
         char_name=char_name1//'.'//char_name2
         write(*,*)'source: irfcntt,char1,char2,char_name = ',
     .       irfcntt,char1,char2,char_name 

c        Save genray Input files
         gfilename='genray.'//char_name//'.dat'
         command = 'cp '//'genray.dat '//gfilename
         write(*,*)'source: command= ',command
         if (ISHELL (command) .lt. 0)   
     .   call STOP ('sub call_genray: failure of spawned cp'//
     .     ' genray.dat command',67)

         gfilename='genray_profs_in.'//char_name//'.nc'
         command = 'cp '//'genray_profs_in.nc '//gfilename
         write(*,*)'source: command= ',command
         if (ISHELL (command) .lt. 0)   
     .   call STOP ('sub call_genray: failure of spawned cp'//
     .     'genray_profs_in.nc command',67)

         gfilename='genray_eqdsk'//'.'//char_name
         write(*,*)'source: eqdskfilename = ',eqdskfilename
         command = 'cp '//eqdskfilename(1:LEN_TRIM(eqdskfilename))//
     .        ' '//gfilename
         write(*,*)'source: command= ',command
         if (ISHELL (command) .lt. 0)   
     .   call STOP ('sub call_genray: failure of spawned cp'//
     .     ' eqdskfilename command',67)

      endif



      
c     GET GENRAY EXECUTABLE

      if(first_time)then
         !get fully qualified name of genray to run:
          call get_genray(ncrt,nout,genray_path_out,len_str)
          first_time = .false.
      endif
      command = ADJUSTL(genray_path_out(1:LEN_TRIM( genray_path_out)))



c     EXECUTE RAY-TRACING CODE (GENRAY) 
c     and then connect to file genray_profs_out
c
      write  (6, '(/ '' ---- GENRAY starting'')')
      write(6, FMT ='(" running : ",a)')command
c
      if (ISHELL (command) .lt. 0)
     .  call STOP ('subroutine CALL_GENRAY: failure of spawned GENRAY',
     .  67)
c
      write  (6, '(/ '' ---- GENRAY finished'')')
c
c

c     SAVE COPY OF OUTPUT FILES (save_genray_io .eq. 1)

      if (save_genray_io .eq. 1) then

         gfilename='genray_profs_out.'//char_name//'.nc'
         command = 'cp '//'genray_profs_out.nc '//gfilename
         write(*,*)'source: command= ',command
         if (ISHELL (command) .lt. 0)   
     .   call STOP ('sub call_genray: failure of spawned cp'//
     .     'genray_profs_out.nc command',67)

      endif



c$$$c
c$$$c  -----------------------------------------------------------------------
c$$$c     Read data produced by ray-tracing code,
c$$$c     and put into form for output variables in this subroutine
c$$$c  -----------------------------------------------------------------------
c$$$c
c$$$      call getioun(igenray,42)
c$$$      open (unit = igenray, file = 'genray_profs_out', status = 'OLD')
c$$$c
c$$$      read   (igenray, '(3i5)') nbulk,nprofs,indexrho
c$$$                               !#species,grid size,radial coord type
c$$$      write(*,*)'call_genr:nbulk,nprofs,indexrho',nbulk,nprofs,indexrho
c$$$c     Check # species agrees with expected #, # of radial mesh points:
c$$$      if (nbulk.ne.1 .and. nbulk.ne.nspecgr) then
c$$$         write(*,*)'call_genray:nbulk,nspecgr= ',nbulk,nspecgr
c$$$         write(*,*)'call_genray: problem with nbulk'
c$$$         STOP
c$$$      endif
c$$$      if (nprofs.ne.nj) then
c$$$         write(*,*)'call_genray: problem with nj'
c$$$         STOP
c$$$      endif
c$$$c     Check radial coord of genray data is sqrt(tor flux)
c$$$      if (indexrho.ne.2) then
c$$$         write(*,*)'call_genray: problem with genray radial coord type'
c$$$         STOP
c$$$      endif
c
c
c  -----------------------------------------------------------------------
c     Read netCDF data produced by ray-tracing code,
c     and put into form for output variables in this subroutine
c  -----------------------------------------------------------------------
c
      ncid = ncopn('genray_profs_out.nc',NCNOWRIT,istatus)
      write(*,*)'after ncopn ncid=',ncid,'istatus',istatus
c.......................................................................
c     read in dimension IDs
      nprofs_id = ncdid(ncid,'nprofs_dim',istatus)
      write(*,*)'after ncdid nprofs_id=',nprofs_id,'istatus',istatus
      nbulk_id = ncdid(ncid,'nbulk_dim',istatus)
      write(*,*)'after ncdid nbulk_id=',nbulk_id,'istatus',istatus

c --- inquire about dimension sizes:#species,grid size---
      call ncdinq(ncid,nprofs_id,name,nprofs,istatus)
      call ncdinq(ncid,nbulk_id,name,nbulk,istatus)

      if (nbulk.gt.1) then
         nbulkm_id = ncdid(ncid,'nbulkm_dim',istatus)
         write(*,*)'after ncdid nbulkm_id=',nbulkm_id,'istatus',istatus
         call ncdinq(ncid,nbulkm_id,name,nbulkm,istatus)
         countm(1)=nprofs
         countm(2)=nbulkm
      endif

      vid = ncvid(ncid,'indexrho',istatus)   !radial coord type
      call ncvgt(ncid,vid,1,1,indexrho,istatus)

      write(*,*)'call_genray:nbulk,nprofs,indexrho ',
     .     nbulk,nprofs,indexrho

c     Check # species agrees with expected #, # of radial mesh points:
      write(*,*)'call_genray:nbulk,nspecgr=  ',nbulk,nspecgr
      if (nbulk.ne.1 .and. nbulk.ne.nspecgr) then
         write(*,*)'call_genray: problem with nbulk'
         STOP
      endif
      if (nprofs.ne.nj) then
         write(*,*)'call_genray: problem with nj'
         STOP
      endif
c     Check radial coord of genray data is sqrt(tor flux)
      if (indexrho.ne.2) then
         write(*,*)'call_genray: problem with genray radial coord type'
         STOP
      endif

      !allocate storage and initialize to zero:
c  -----------------------------------------------------------------------

      if(.not. allocated(rgenray))then
         allocate (rgenray(nprofs),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("rgenray,sub call_genray",0,istat)
      endif
      rgenray=(/(0d0,i=1,nprofs)/)

      if(.not. allocated(pgre))then
         allocate (pgre(nprofs),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("pgre,sub call_genray",0,istat)
      endif
      pgre=(/(0d0,i=1,nprofs)/)

      if(.not. allocated(pgrc))then
         allocate (pgrc(nprofs),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("pgrc,sub call_genray",0,istat)
      endif
      pgrc=(/(0d0,i=1,nprofs)/)

      if(.not. allocated(pgrit))then
         allocate (pgrit(nprofs),STAT = istat)
         if(istat .ne. 0)
     .        call allocate_error("pgrit,sub call_genray",0,istat)
      endif
      pgrit=(/(0d0,i=1,nprofs)/)
      
      if(.not. allocated(pgri))then
         allocate (pgri(nprofs,nspecgr-1),STAT = istat)
         if(istat .ne. 0)
     .     call allocate_error("pgri,sub call_genray",0,istat)
      endif
      do j=1,nspecgr-1
         do i=1,nprofs	 
            pgri(i,j)=0d0	
         enddo			
      enddo			
      
      if(.not. allocated(totgrps))then
         allocate (totgrps((nspecgr-1)),STAT = istat)
         if(istat .ne. 0)
     .        call allocate_error("totgrps,sub call_genray",0,istat)
      endif
      totgrps=(/(0d0,i=1,nspecgr-1)/)

      totgrp=0d0
      totgrpe=0d0
      totgrpi=0d0
      totgrc=0d0
      

c  ----------------------------------------------------------------------

c     continue reading data:
c$$$      read   (igenray,   3010 )  (rgenray(j), j=1,nprofs) !normalized rho ?
c$$$      read   (igenray,   3010 )  (pgre(j), j=1,nprofs) !power to electrons
c$$$      read (igenray, 3020)  totgrp, totgrpe, totgrc
c$$$         read (igenray, 3010)  (pgrc(j),j=1,nprofs)  !rf current density

c     For a pure electron heating EC case (nbulk=1), no power calculation
c     results are passed from genray.
c$$$      if (nbulk.gt.1) then
c$$$         read   (igenray,   3010 )  (pgrit(j), j=1,nprofs)
c$$$                                    !total power to ions
c$$$         do n=1,nbulk-1
c$$$           read (igenray, 3010 )  (pgri(j,n),j=1,nprofs)
c$$$                                  !power to ion specie
c$$$         end do
c$$$         read (igenray, 3020)  totgrpi, totgrps(1:nbulk-1)
c$$$      endif
c$$$      
c$$$ 3010 format (5(1pe19.6))
c$$$ 3020 format (4(1pe19.6))

      vid = ncvid(ncid,'rgenray',istatus)   
      call ncvgt(ncid,vid,1,nprofs,rgenray,istatus) !normalized rho

      vid = ncvid(ncid,'pgre',istatus)   
      call ncvgt(ncid,vid,1,nprofs,pgre,istatus)    !power den to electrons

      vid = ncvid(ncid,'totgrpe',istatus)           !tot pwr to electrons
      call ncvgt(ncid,vid,1,1,totgrpe,istatus)

      if (nbulk.gt.1) then
         vid = ncvid(ncid,'pgri',istatus)       !pwr den to each ion specie 
         call ncvgt(ncid,vid,start,countm,pgri,istatus)
         
         vid = ncvid(ncid,'pgrit',istatus)   
         call ncvgt(ncid,vid,1,nprofs,pgrit,istatus) !power density to ions
         
         vid = ncvid(ncid,'totgrps',istatus)   
         call ncvgt(ncid,vid,1,nbulk-1,totgrps,istatus)!total pwr to each ion
         
         vid = ncvid(ncid,'totgrpi',istatus)           !total pwr to ions 
         call ncvgt(ncid,vid,1,1,totgrpi,istatus)
      endif

      vid = ncvid(ncid,'totgrp',istatus)            !total input power to  
      call ncvgt(ncid,vid,1,1,totgrp,istatus)       !genray for this mode

      vid = ncvid(ncid,'pgrc',istatus)              !rf current density
      call ncvgt(ncid,vid,1,nprofs,pgrc,istatus)    !<j.B/Beqd>

      vid = ncvid(ncid,'totgrc',istatus)            !total tor current
      call ncvgt(ncid,vid,1,1,totgrc,istatus)

c
c  ----------------------------------------------------------------------

c     renormalize power fraction to ions and electrons if called for:
      if(genray_fi .gt.  0.0 .and. genray_fi .lt. 1.0)then
        do j=1,nprofs
           pgrtoti =0d0
           do n=1,nspecgr-1
              pgrtoti    = pgrtoti + pgri(j,n)
           enddo
           pgrtot   = pgrtoti  + pgre(j)
           pgre(j) = (1.-genray_fi)*pgrtot
           do n=1,nspecgr-1
             pgri(j,n) = genray_fi*pgrtot*pgri(j,n)/pgrtoti
           enddo
           pgrit(j) = genray_fi*pgrtot
        enddo
      endif

c
      do j=1,nprofs
        rgenray(j) = rgenray(j) * r(nj)                 !convert rgenray to cm
      enddo

      if (nbulk .gt. 1) then
      write(*,*)' '
      write(*,*)'******************************************************'
      write(*,*)'call_genray:  For the time being, power to a hot'
      write(*,*)'call_genray:  Maxwellian ion distn representing beam'
      write(*,*)'call_genray:  ions is calculated in genray, returned'
      write(*,*)'call_genray:  to ONETWO, but not further processed.'
      write(*,*)'call_genray:  One approx possibility is to divvy it'
      write(*,*)'call_genray:  as additional source power to  e and i'
      write(*,*)'call_genray:  based on Coulomb rates.'
      write(*,*)'******************************************************'
      write(*,*)' '
      endif


c
c     Cast powers and current onto ONETWO grid
c
cBH  The rgenray should be the same as the r(), for genray      
cBH      call intrp (0, 1, rgenray, pwe , nprofs, r, qrfe , nj)
cBH      call intrp (0, 1, rgenray, pwit, nprofs, r, qrfi , nj)
cBH      call intrp (0, 1, genray, pwc , nprofs, r, currf, nj)

      do j=1,nj
         qrfe(j)=pgre(j)
         qrfi(j)=pgrit(j)
         currf(j)=pgrc(j)
      enddo

c
c     Smooth out profiles qrfe, qrfi, currf. Note that
c     no smoothing is the default(nfwsmth,nicsmth,necsmth = 0) so this
c     is in effect only if user requested through inone.
c
      nsmth = MAX(nfwsmth,nicsmth) ! may separate these in the future
      call smthsj (r, qrfe , nj, nsmth)
      call smthsj (r, qrfi , nj, nsmth)
      call smthsj (r, currf, nj, nsmth)
      if ( necsmth.ne.0 .and. nsmth.eq.0 ) then
         call smthsj (r, qrfe , nj, nsmth)
         call smthsj (r, currf, nj, nsmth)
      endif
c
c
c     Total up powers and currents
c
      call trapv (r, qrfe , hcap, nj, totgrpe)
      call trapv (r, qrfi , hcap, nj, totgrpi)
      call trapv (r, currf, hcap, nj, totgrc)
c
      pi      = ATAN2 (0d0, -1d0)
      totgrpe = 4.0 * pi**2 * 100.0 * rmajor * totgrpe
      totgrpi = 4.0 * pi**2 * 100.0 * rmajor * totgrpi
      totgrc  = 2.0 * pi * totgrc
      print *,'call_genray:tot pe,pi,jc =',totgrpe,totgrpi,totgrc
c
	    
      call giveupus(igenray)
      close (unit = igenray)
      if (nbulk.ne.1)
     +     deallocate (pgri)    ! nspecgr may change from call to call


      return
c
      end   !end call_genray
c     
c     
      integer function length_char(string)
c     Returns length of string, ignoring trailing blanks,
c     using the fortran intrinsic len().
      character*(*) string
      do i=len(string),1,-1
         if(string(i:i) .ne. ' ') goto 20
      enddo
 20   length_char=i
      return
      end
c
c
      subroutine cgbco (abd, lda, n, ml, mu, ipvt, rcond, z)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      integer    lda, n, ml, mu, ipvt(*)
      real*8     rcond
      complex*16 abd(lda,*), z(*), cdotc, ek, t, wk, wkm
      real*8     anorm, s, scasum, sm, ynorm
      integer    is, info, j, ju, k, kb, kp1, l, la, lm, lz, m, mm
      complex*16 zdum, zdum1, zdum2, csign1
      real*8     cabs1
c
      cabs1 (zdum       ) = ABS (REAL (zdum)) + ABS (PIMAG (zdum))
      csign1(zdum1,zdum2) = cabs1 (zdum1) * (zdum2 / cabs1 (zdum2))
c
c     compute 1-norm of a
c
      zero  = 0.0
      anorm = 0.0e0
      L     = ML + 1
      IS    = L + MU
      do 10 j=1,n
         anorm = MAX (anorm, scasum(l, abd(is,j), 1))
         if (IS .gt. ML + 1)  IS = IS - 1
         if ( J .le. MU    )  L  = L + 1
         if ( J .ge. N - ML)  L  = L - 1
   10 continue
c
c     FACTOR
c
      call cgbfa (abd,lda,n,ml,mu,ipvt,info)
c
c     rcond = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
c     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E .
c     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A .
c     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
c     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E .
c     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
c
c     SOLVE CTRANS(U)*W = E
c
      EK = (1.0e0,0.0e0)
      do j=1,n
        Z(J) = (0.0e0,0.0e0)
      end do
      M  = ML + MU + 1
      JU = 0
      do 100 k=1,n
         if (cabs1(Z(K)) .ne. 0.0e0)  ek = csign1 (ek,-z(k))
         if (cabs1(EK-Z(K)) .le. cabs1(abd(M,K)))  go to 30
            s  = cabs1 (abd(M,K)) / cabs1 (ek-z(k))
            call csscal(n,s,z,1)
            ek = CMPLX (s, zero) * ek
   30    continue
         WK  = EK - Z(K)
         WKM = -EK - Z(K)
         S   = cabs1(WK)
         SM  = cabs1(WKM)
         if (cabs1(abd(M,K)) .eq. 0.0e0)  go to 40
            WK  = WK/CONJG(abd(M,K))
            WKM = WKM/CONJG(abd(M,K))
         go to 50
   40    continue
            WK  = (1.0e0,0.0e0)
            WKM = (1.0e0,0.0e0)
   50    continue
         KP1 = K + 1
         JU  = MIN0 (MAX0 (JU,MU+IPVT(K)),N)
         MM  = M
         if (KP1 .gt. JU)  go to 90
            do 60 j=kp1,ju
               MM   = MM - 1
               SM   = SM + cabs1(Z(J)+WKM*CONJG(abd(MM,J)))
               Z(J) = Z(J) + WK*CONJG(abd(MM,J))
               S    = S + cabs1(Z(J))
   60       continue
            if (S .ge. SM)  go to 80
               T  = WKM - WK
               WK = WKM
               MM = M
               do 70 j=kp1,ju
                  MM   = MM - 1
                  Z(J) = Z(J) + T*CONJG(abd(MM,J))
   70          continue
   80       continue
   90    continue
         Z(K) = WK
  100 continue
      S = 1.0e0/scasum(N,Z,1)
      call csscal(n,s,z,1)
c
c     SOLVE CTRANS(L)*Y = W
c
      do kb=1,n
         K  = N + 1 - KB
         LM = MIN0 (ML,N-K)
         if (K .lt. N) Z(K) = Z(K) + cdotc(LM,abd(M+1,K),1,Z(K+1),1)
         if (cabs1(Z(K)) .le. 1.0e0)  go to 110
            S = 1.0e0/cabs1(Z(K))
            call csscal(n,s,z,1)
  110    continue
         L    = IPVT(K)
         T    = Z(L)
         Z(L) = Z(K)
         Z(K) = T
      end do
      S = 1.0e0/scasum(N,Z,1)
      call csscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     SOLVE L*V = Y
c
      do k=1,n
         L    = IPVT(K)
         T    = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         LM   = MIN0 (ML,N-K)
         if (K .lt. N) call caxpy(LM,T,abd(M+1,K),1,Z(K+1),1)
         if (cabs1(Z(K)) .le. 1.0e0)  go to 130
           S = 1.0e0/cabs1(Z(K))
           call csscal(n,s,z,1)
           ynorm = s * ynorm
  130    continue
      end do
      S     = 1.0e0/scasum(N,Z,1)
      call csscal(n,s,z,1)
      ynorm = s * ynorm
c
c     SOLVE  U*Z = W
c
      do kb=1,N
         K = N + 1 - KB
         if (cabs1(Z(K)) .le. cabs1(abd(M,K)))  go to 150
            S     = cabs1(abd(M,K))/cabs1(Z(K))
            call csscal(n,s,z,1)
            ynorm = S*ynorm
  150    continue
         if (cabs1(abd(M,K)) .ne. 0.0e0) Z(K) = Z(K)/abd(M,K)
         if (cabs1(abd(M,K)) .eq. 0.0e0) Z(K) = (1.0e0,0.0e0)
         LM = MIN0 (K,M) - 1
         LA = M - LM
         LZ = K - LM
         T  = -Z(K)
         call caxpy(LM,T,abd(LA,K),1,Z(LZ),1)
      end do
c
c     MAKE ZNORM = 1.0
c
      S     = 1.0e0/scasum(N,Z,1)
      call csscal(n,s,z,1)
      ynorm = S*ynorm
c
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
c
      end

      subroutine cgbfa (abd, lda, n, ml, mu, ipvt, info)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      integer    lda,n,ml,mu,ipvt(*),info
      complex*16 abd(lda,*), t
      integer    i,icamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
      complex*16 zdum
      real*8     cabs1
c
      cabs1 (zdum) = ABS (REAL (zdum)) + ABS (PIMAG (zdum))
c
      m    = ml + mu + 1
      info = 0
c
c     ZERO INITIAL FILL-IN COLUMNS
c
      j0 = mu + 2
      j1 = MIN0 (n, m) - 1
      if (j1 .lt. j0)  go to 30
      do 20 jz=j0,j1
         i0 = m + 1 - jz
         do 10 i=i0,ml
            abd(i,jz) = (0.0e0, 0.0e0)
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
c
      nm1 = n - 1
      if (nm1 .lt. 1)  go to 130
      do 120 k=1,nm1
         kp1 = k + 1
c
c        ZERO NEXT FILL-IN COLUMN
c
         jz = jz + 1
         if (JZ .gt. N)  go to 50
         if (ML .lt. 1)  go to 50
            do 40 i=1,ml
               abd(I,JZ) = (0.0e0,0.0e0)
   40       continue
   50    continue
c
c        FIND L = PIVOT INDEX
c
         LM      = MIN0   (ML,N-K)
         L       = icamax (LM+1, abd(M,K), 1) + M - 1
         IPVT(K) = L + K - M
c
c        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
c
         if (cabs1(abd(L,K)) .eq. 0.0e0)  go to 100
c
c           INTERCHANGE IF NECESSARY
c
            if (L .eq. M)  go to 60
               T        = abd(L,K)
               abd(L,K) = abd(M,K)
               abd(M,K) = T
   60       continue
c
c           COMPUTE MULTIPLIERS
c
            T = -(1.0e0,0.0e0)/abd(M,K)
            call cscal(LM,T,abd(M+1,K),1)
c
c           ROW ELIMINATION WITH COLUMN INDEXING
c
            JU = MIN0 (MAX0 (JU,MU+IPVT(K)),N)
            MM = M
            if (JU .lt. KP1)  go to 90
            do 80 j=kp1,ju
               L  = L - 1
               MM = MM - 1
               T  = abd(L,J)
               if (L .eq. MM)  go to 70
                  abd( L,J) = abd(MM,J)
                  abd(MM,J) = T
   70          continue
               call caxpy(LM,T,abd(M+1,K),1,abd(MM+1,J),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = K
  110    continue
  120 continue
  130 continue
      IPVT(N) = N
      if (cabs1(abd(M,N)) .eq. 0.0e0) info = N
      return
c
      end

      subroutine cgbsl (abd, lda, n, ml, mu, ipvt, b, job)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      integer lda,n,ml,mu,ipvt(*),job
      complex*16 abd(lda,*), b(*)
c
      complex*16 cdotc,T
      integer k,kb,l,la,lb,lm,m,nm1
c
      M   = MU + ML + 1
      NM1 = N - 1
      if (JOB .ne. 0)  go to 50
c
c        JOB = 0 , SOLVE  A * X = B
c        FIRST SOLVE L*Y = B
c
         if (ML .eq. 0)  go to 30
         if (NM1 .lt. 1)  go to 30
            do 20 k=1,nm1
               LM = MIN0 (ML,N-K)
               L  = IPVT(K)
               T  = B(L)
               if (L .eq. K)  go to 10
                  B(L) = B(K)
                  B(K) = T
   10          continue
               call caxpy(LM,T,abd(M+1,K),1,B(K+1),1)
   20       continue
   30    continue
c
c        NOW SOLVE  U*X = Y
c
         do 40 kb=1,n
            K    = N + 1 - KB
            B(K) = B(K)/abd(M,K)
            LM   = MIN0 (K,M) - 1
            LA   = M - LM
            LB   = K - LM
            T    = -B(K)
            call caxpy(LM,T,abd(LA,K),1,B(LB),1)
   40    continue
      go to 100
   50 continue
c
c        JOB = NONZERO, SOLVE  CTRANS(A) * X = B
c        FIRST SOLVE  CTRANS(U)*Y = B
c
         do k=1,n
            LM   = MIN0 (K,M) - 1
            LA   = M - LM
            LB   = K - LM
            T    = cdotc(LM,abd(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/CONJG(abd(M,K))
         end do
c
c        NOW SOLVE CTRANS(L)*X = Y
c
         if ( ML .eq. 0)  go to 90
         if (NM1 .lt. 1)  go to 90
            do 80 kb=1,nm1
               K    = N - KB
               LM   = MIN0 (ML,N-K)
               B(K) = B(K) + cdotc(LM,abd(M+1,K),1,B(K+1),1)
               L    = IPVT(K)
               if (L .eq. K)  go to 70
                  T    = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          continue
   80       continue
   90    continue
  100 continue
      return
c
      end
c       
c
      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
         write(*,*)  'check_err:  netCDF error'
         stop
      endif
      return
      end



           subroutine cspline (x, y, n, yp1, ypn, y2)
c------------------------------------------------------------
c           from gafit_lib
c-----------------------------------------------------------HSJ

c
      implicit none
c
      character rcs_id*63
      save rcs_id
      data rcs_id/
     &"$Id: cray331.f,v 1.101 2013/09/03 18:53:21 stjohn Exp $"/
c
      integer i, n, nmax, k
      doubleprecision p, sig, qn, un
      parameter (nmax = 129)
      doubleprecision x(n), y(n), y2(n), u(nmax), yp1, ypn
c
      if (n .gt. nmax) call STOP (
     &'subroutine CSPLINE: dimensional error', 201)
c
      if (yp1 .gt. 0.99d30) then
      y2(1) = 0.d0
      u(1) = 0.d0
      else
      y2(1) = -0.5d0
      u(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
      sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
      p = sig*y2(i-1)+2.d0
      y2(i) = (sig-1.d0)/p
      u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) /(x(i)-x(i
     &-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn .gt. 0.99d30) then
      qn = 0.d0
      un = 0.d0
      else
      qn = 0.5d0
      un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k=n-1,1,-1
      y2(k) = y2(k)*y2(k+1)+u(k)
      enddo
      return
c
      end





      complex*16 function det (xx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      complex*16 xx,kxx,kxy,kxz,kyy,kyz,kzz
      complex*16 zlam,besi,besd,cxq
c
      common /rfzl/ zlam(5),besi(7,5),besd(7,5)
      common /rf4 / xnparz,xnpar,dbg(99)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      common /rffq/ wc(5),wp2(5),wp(5),vtx(5),vtz(5)
c
      dimension  besr(7), besa(7)
c
      zero = 0.0
      one  = 1.0
      det  = CMPLX (0.0, 0.0)
      cxq  = SQRT (xx)
      if (REAL (cxq) .lt. 0.0)  cxq = -cxq
      do ii=1,7
        besr(ii) = 0.0
        besa(ii) = 0.0
      end do
      dbg(11) =  REAL (xx)
  100 dbg(12) = PIMAG (xx)
      do 1 k=1,neq
      zlam(k) = 0.5*xx*(vtx(k)/wc(k))**2
      dbg(50) =  REAL (xx)
      go to 122
  122 dbg(51) = PIMAG (xx)
      zlr     =  REAL (zlam(k))
      go to 101
  101 zli = PIMAG (zlam(k))
      if (zlam(k) .eq. CMPLX (0.0, 0.0))  go to 223
      call beslci (zlr, zli, nhp, 1, besr, besa, ncalc)
      besi(nhp,k) = CMPLX (besr(1), zero) + CMPLX (zero, one)*besa(1)
      besd(nhp,k) = CMPLX (besr(2), zero) + CMPLX (zero, one)*besa(2)
      dbg(1) = besr(1)
      dbg(2) = besa(1)
      dbg(3) = besr(2)
      go to 102
  102 dbg(4) = besa(2)
c
  223 do 2 mm=2,nhigh+1
      m = mm - 1
      if (zlam(k) .eq. CMPLX (0.0, 0.0))  go to 202
      besi(m,k) = CMPLX (besr(mm), zero) + CMPLX (zero, one)*besa(mm)
      dbg(5)    = besr(mm)
      go to 103
  103 dbg(6) = besa(mm)
      if (m .ne. 1)  besd(m,k) = besi(m-1,k)-m*besi(m,k)/zlam(k)
      if (m .eq. 1)  besd(1,k) = besi(nhp,k)-besi(1,k)/zlam(k)
      go to 201
  202 besd(m,k)   = CMPLX (0.0, 0.0)
      besd(1,k)   = CMPLX (0.5, 0.0)
      besd(nhp,k) = CMPLX (0.0, 0.0)
      besi(m,k)   = CMPLX (0.0, 0.0)
      besi(nhp,k) = CMPLX (1.0, 0.0)
c
  201 dbg(7) =  REAL (besd(m,k))
      go to 104
  104 dbg(8) = PIMAG (besd(m,k))
    2 continue
    1 continue
      xn = (clight/(2.0 * pi * freq))**2
      xnz = xkpar*xkpar*xn
      xnpar = xn
      xnparz = xnz
      go to 107
  107 det = ((kxx(xx)-xnz)*((kyy(xx)-xnz-xx*xn)
     .    * (kzz(xx)-xx*xn)+kyz(xx)*kyz(xx))
     .      +kxy(xx)*(kxy(xx)*(kzz(xx)-xx*xn)
     .      +kyz(xx)*(kxz(xx)+cxq*xnz/xkpar))
     .     +(kxz(xx)+cxq*xnz/xkpar)*(kxy(xx)*kyz(xx))
     .     -(kyy(xx)-xnz-xx*xn)*(kxz(xx)+cxq*xnz/xkpar))
      dbg( 9) =  REAL (det)
      go to 106
  106 dbg(10) = PIMAG (det)
      return
c
      end
      subroutine ech (freq, rfpow, wrfo, idamp, gafsep, necsmth,
     .                codeid, xsrc, zsrc, nray, thet, phai, hlwec,
     .                ratwec, nj, r, psir, ene, te, q, hcap, dr,
     .                qrfe, currf,jtor_rf, totecpe, totecc, ifixshap,
     .                time,iplotec, r0, bt0, zeff, rf_output)
c
      USE param
      USE ext_prog_info, ONLY : echin_save,toray_version,
     .                          toray_version_switch,get_toray
      USE mhdpar
      USE io 
      USE ename
      USE constnts,       ONLY : twopi
      use rf,             ONLY : ech_input
      USE mhdcom
      USE file_proc,      ONLY : append_time_to_filename,
     .                           rename_disk_file,
     .                           append_gyro_id_to_filename,
     .                           delete_file
     
      USE geom, ONLY : rcapi


      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      CHARACTER ech_filename*256
      CHARACTER torayinpt*256
      save      rcs_id
      data      rcs_id /
     ."$Id: cray331.f,v 1.101 2013/09/03 18:53:21 stjohn Exp $"/
c
c       parameter (nx_param = 129, ny_param = 129) ! HSJ 10/24/96
c

c     do not include mesh.i as there are conflicts below, see for example ra, - HSJ
c
      external     LENGTH, FLUSH
      integer      LENGTH, rf_output
      logical      first_time,ex
      character*8  codeid, ntitle(5)
      character*12 cgmfile, glogfile, tlogfile, gafit_log, toray_log,
     .              operator*4
      character (len = 256) program, toray_to_run,gafiterr
      integer len_toray_to_run,read_toray_in
      dimension    r(nj), ene(nj), te(nj), q(nj), qrfe(nj),
     .             currf(nj), psir(nj), hcap(nj), dr(nj), zeff(nj),
     .             pdata(nw,nh), dumaray5(5)
      REAL*8 jtor_rf(nj)
c
      dimension      rp_eqd(2001),   zp_eqd(2001),
     .             rlim_eqd(2001), zlim_eqd(2001)


      dimension  xbouni(ktoray),vb(ktorm1),qrfeb(ktoray),tpowde(ktorm1),
     .         wiecrt(ktoray),tidept(ktorm1),psinrm(kj),dxbouni(ktorm1)
      REAL *8,dimension(:),allocatable :: currf_toray,dumy
c
      character (len = 24) task
      save         first_time        , igafit   , itoray,
     .             toray_to_run,len_toray_to_run
c
      data         first_time/.true./, igafit/0/, itoray/0/,
     .             iecall/0/, icount/0/,read_toray_in /0/
c
      if(ktoray .ne. kj)then 
         print *,'ktoray =',ktoray
         print *,'kj =',kj
         call STOP('Subroutine ech, ktoray must equal kj',0)
      endif
      nx_param = nw
      ny_param = nh
      nbfld  = 1
      jtor_rf(:) = 0.0d0
      task ='write psiin'

c     get fully qualified name of toray to run:
c                          toray_to_run(1:len_toray_to_run)

      if(first_time)
     . call get_toray(nx_param,ny_param,ktoray,ncrt,nout,
     .          toray_to_run,len_toray_to_run)

          !required for echin:
          psidif = psir(1)-psir(nj)
          do   j=1,nj
             psinrm(j) = (psir(1)-psir(j))/psidif
          enddo
          xmaxis = rma/100.
          rmajor = r0/100.                           !added 03/10/03 HSJ
c-----------------------------------------------------------------------  
c     decide if gafit should  be  called, label 19 bypasses gafit section:
c------------------------------------------------------------------------
c

      if(toray_version .lt. toray_version_switch)then
           !old version of toray does not accept igafit input in
           !file toray.in. Hence if toray.in exists we must check to see
           !that toray.in is appropriate.
           call rwrt_toray_in
c           print * ,'gafsep =',gafsep
           if (gafsep .eq. 0.0 .or.    codeid .eq. 'onedee')  go to 19
           nbfld  = 3
           !if this is a fixed shape plasma (ifixshap = 1) and
           !we have called gafit before then skip it this time
           !because file psiin exists already:
           if (iecall .ne. 0   .and. ifixshap .eq.  1      )  go to 19
           iecall = iecall + 1
           ldata  = 2
           call rf_mhddat(task)   !writes psiin file
      else     !new version of toray, do not call gafit directly.
               !if gafit is used at all it will be called from
               !toray directly
         nbfld  = 3
         iecall = iecall + 1
         ldata  = 2
         !if(read_toray_in .eq. 0) call rwrt_toray_in

         call rwrt_toray_in

         read_toray_in = 1
          call rf_mhddat(task)   !writes psiin file

c         call rf_mhd_interface        !write file  mhddat (read  by toray)directly
          print*,'skipping gafit spawn from Onetwo'
          go to 19
      endif


c
c ----------------------------------------------------------------------
c --- read an existing eqdsk. this will be the eqdsk given in inone
c --- if the equilibrium is not changed (ifixshap = 1). otherwise
c --- if new equilibria are being calculated (ifixshap = 0) then
c --- the eqdsk to be read is the latest one written by the code itself.
c --- see subroutine WREQDSK .
c ----------------------------------------------------------------------
c

       print *,'done rf_mhddat, eqdsk_tdem = ',
     .                           eqdsk_tdem(1:LEN_TRIM(eqdsk_tdem))
c
       irftest = 1  !disable the folowing, remove the dead code later HSJ 3/06/03
      if (eqdsk_tdem .ne. 'tdem' ) then
c
        eqfile = eqdskfilename

        if (eqfile .eq. 'none')  eqfile = eqdskin
c
        call getioun(neq,neq)
        open (unit = neq, file = eqfile, status = 'OLD', err = 500)
c
        print *,'eqdsk for toray =',eqfile(1:LEN_TRIM(eqfile))
        go to 11
c
  500   write  (6, 12) eqfile
   12   format (' subroutine ECH reports:'         /
     .          '   unable to open eqdsk file ', a /
     .          '   ONETWO must stop'               )
        call giveupus(neq)
        call STOP ('subroutine ECH: cannot open EQDSK file', 50)
c
   11   read (neq, '(6a8, 3i4)') (ntitle(i), i=1,5), vid, ipestg, nx, nz
        read (neq, 1006)  xdim, zdim, rmajor, redgem, zmid
        read (neq, 1006)  xmaxis, zmaxis, psimax, psilim, btorus
        print *,'idamp =',idamp
        if (idamp .lt. 0)  go to 19
        read (neq, 1006)  dumaray5
        read (neq, 1006)  dumaray5
c
c        call getioun(nscr,nscr)
****    open (unit = nscr, file = 'psiin', status = 'UNKNOWN')
****    pssiin_old goes with rf_mhddat interface
****    (whereas psiin goes with rf_mhd_interface)
c        open (unit = nscr, file = 'psiin_old', status = 'UNKNOWN')
c
c        write (nscr, 1005)  ldata, nj, nx, nz
        ymid = 0.0
c        write (nscr, 1006)  xdim,zdim,rmajor,redgem,ymid
c        write (nscr, 1006)  xmaxis,zmaxis,psimax,psilim,btorus
c        write (nscr, 1006)  gafsep
         read  (neq , 1006)  (pdata(i,1), i=1,nx)
c        write (nscr, 1006)  (pdata(i,1), i=1,nx)
        read  (neq , 1006)  (pdata(i,1), i=1,nx)
        read  (neq , 1006)  (pdata(i,1), i=1,nx)
        read  (neq , 1006)  (pdata(i,1), i=1,nx)
        read  (neq , 1006) ((pdata(i,j), i=1,nx), j=1,nz)
c        write (nscr, 1006) ((pdata(i,j), i=1,nx), j=1,nz)
c        write (nscr, 1006)  (psir(i), i=1,nj)
****    close (unit = nscr)
        read  (neq , 1006)  (pdata(i,1), i=1,nx)
c        write (nscr, 1006)  (pdata(i,1), i=1,nx)
        read  (neq , '(2i5)')  np_eqd, nlim_eqd
c        write (nscr, '(2i5)')  np_eqd, nlim_eqd
        read  (neq , 1006)  (rp_eqd(i), zp_eqd(i), i=1,np_eqd)
c        write (nscr, 1006)  (rp_eqd(i), zp_eqd(i), i=1,np_eqd)
        read  (neq , 1006)  (rlim_eqd(i), zlim_eqd(i), i=1,nlim_eqd)
c        write (nscr, 1006)  (rlim_eqd(i), zlim_eqd(i), i=1,nlim_eqd)
c        call giveupus(nscr)
c        close (unit = nscr)
cyr99
        psidif = psir(1)-psir(nj)
        do 221 j=1,nj
 221       psinrm(j) = (psir(1)-psir(j))/psidif
cyr99
c
c       if q is used takes absolute value here
c
        call giveupus(neq)
        close (unit = neq)        ! close file each time
c
      else  ! get equivalent info from netCDF file
c
        nx = nx_param
        nz = ny_param
c        print *,'netcdf section,nx,ny =',nx,nz
        call ech_netcdf_interface (pdata,nx,nz,rlim_eqd,zlim_eqd,
     .                             nlim_eqd,rp_eqd,zp_eqd,np_eqd,
     .                             xdim,zdim,rmajor,redgem,zmid,xmaxis,
     .                             zmaxis,psimax,psilim,btorus,time)
        write (6, '("ECH data read from netCDF eqdsk file")')

        if(irftest .eq. 0)then
        call getioun(nscr,nscr)
****    open (unit = nscr, file = 'psiin'    , status = 'UNKNOWN')
****    pssiin_old goes with rf_mhddat interface
****    (whereas psiin goes with rf_mhd_interface)
        open (unit = nscr, file = 'psiin_old', status = 'UNKNOWN')
c
        write (nscr, 1005)  ldata,nj,nx,nz
        write (nscr, 1006)  xdim,zdim,rmajor,redgem
c        write   (nunit,   200  )  rdim,zdim,rcenter,rinside
        write (nscr, 1006)  xmaxis,zmaxis,psimax,psilim,btorus
c       write   (nunit,   200  )  rmaxis,zmaxis,psimax,psilim,b0
        write (nscr, 1006)  gafsep                !gasep
        write (nscr, 1006)  (pdata(i,1), i=1,nx)  !sf(1..nx)
        write (nscr, 1006) ((pdata(i,j), i=1,nx), j=1,nz) !psi(i,j)
        write (nscr, 1006)  (psir(i), i=1,nj)      
        write (nscr, 1006)  (pdata(i,1), i=1,nx)  !qpsi,i=1,nx
        write (nscr, '(2i5)')  np_eqd, nlim_eqd
        write (nscr, 1006)  (  rp_eqd(i),   zp_eqd(i), i=1,np_eqd  )
        write (nscr, 1006)  (rlim_eqd(i), zlim_eqd(i), i=1,nlim_eqd)
        call giveupus(nscr)
        close (unit = nscr)
cyr99
        endif
        psidif = psir(1)-psir(nj)
        do 22 j=1,nj
   22   psinrm(j) = (psir(1)-psir(j))/psidif
cyr99
c

      end if  ! netCDF file, or eqdsk read branch

 



c
c ----------------------------------------------------------------------
c comments below (from INIT) show the role of RF_OUTPUT input variable
c ----------------------------------------------------------------------
c  rf_output 0:  discard standard output of TORAY and GAFIT
c            1:    allow standard output of TORAY and GAFIT to flow <-- DEFAULT
c            2: redirect standard output of TORAY and GAFIT to
c                 "toray.log"    and "gafit.log"    respectively (  overwriting)
c            3: redirect standard output of TORAY and GAFIT to
c                 "toray.log"    and "gafit.log"    respectively (concatenating)
c            4: redirect standard output of TORAY and GAFIT to
c                 "toray_NN.log" and "gafit_NN.log" respectively, NN = 01,02,...
c ----------------------------------------------------------------------
c
c -------- begin GAFIT section -----------------------------------------
c
      print *,'Entering gafit section of Onetwo'
c      program = ADJUSTL('gafit')
c      proper gafit is in same location as toray:
       program =toray_to_run(1:len_toray_to_run-5)//'gafit'
c      is gafit available ??
       inquire(FILE =program(1:len_toray_to_run) ,EXIST  = ex)
       if(.not. ex)then
          gafiterr='ERROR '//program(1:len_toray_to_run)//
     .                                     ' not available'
          write(ncrt,'(a)')gafiterr
          write(nout,'(a)')gafiterr
          call STOP('SUB ECH Requested gafit not found ',0)
       endif
       print *,'Onetwo has determined that gafit to run is'
       print *,program(1:len_toray_to_run)
      igafit  = igafit + 1
      write  (gafit_log, 1002)  'gafit_', igafit, '.log'
 1002 format (a6, i2.2, a4)
c
      if (first_time) then
        if      (rf_output .eq. 0) then
          operator = ' >  '
          glogfile = '/dev/null'
        else if (rf_output .eq. 1) then
          operator = '    '
          glogfile = ' '
        else if (rf_output .eq. 2) then
          operator = ' >  '
          glogfile = 'gafit.log'
        else if (rf_output .eq. 3) then
          operator = ' >> '
          glogfile = 'gafit.log'
          call DESTROY (glogfile)
        else if (rf_output .eq. 4) then
          operator = ' >  '
          glogfile =  gafit_log
        else
          operator = ' >  ' 
          glogfile = gafit_log
          write  (ncrt, 10)  rf_output
   10     format (/
     .  ' WARNING: user input of  RF_OUTPUT =', i3, '  is invalid.'   /
     .       10x, 'RF_OUTPUT may only be 0, 1, 2, 3, or 4,'           /
     .       10x, 'so ONETWO will continue but treat RF_OUTPUT as 4.' /)
          rf_output = 4
        end if
      else if (rf_output .eq. 4) then
        glogfile = gafit_log
      end if
c
      if (rf_output .eq. 1) then
        write  (ncrt , 1001) program(1:len_toray_to_run)
        write  (nitre, 1001) program(1:len_toray_to_run)
 1001   format (/ ' ---- ', a,' started, output not redirected' /)
      else
        write  (ncrt , 1000)  program(1:len_toray_to_run), 
     .                             glogfile(1:LENGTH(glogfile))
        write  (nitre, 1000)  program(1:len_toray_to_run),
     .                             glogfile(1:LENGTH(glogfile))
 1000   format (/ ' ---- ', a,' started, output directed to "', a, '"')
      end if
c
      if (ISHELL (program(1:len_toray_to_run) // operator //
     .                                            glogfile) .ne. 0)
     .  call STOP ('subroutine ECH: failure of spawned GAFIT', 64)
c
      write  (ncrt , 1003)  program(1:5)
      write  (nitre, 1003)  program(1:5)
 1003 format (  ' ---- ', a, ' finished' /)
c
c -------- end GAFIT section -------------------------------------------
c

c----------------------------Echin  Section-----------------------------
c     enter here if gafit is not run directly. In Newer versions
c     of toray gafit may be run by calling toray (rather than gafit)
c     This depends on switch igafit in file toray.in.
c---------------------------------------------------------------------HSJ
   19 call getioun(nscr,nscr)

      open (unit = nscr, file = 'echin', status = 'UNKNOWN')
      if (idamp .lt. 0)  idamp = -idamp

      print *,'idamp,nbfld =',idamp,nbfld
      ra = r(nj)
      write  (nscr, 1006)  time
 1006 format (5e16.9)
      write  (nscr, 1005)  idamp, nj, nray, nbfld
 1005 format (20i4)
      write  (nscr, 1006)  freq,wrfo,xsrc,zsrc,thet,phai,
     .                     hlwec,ratwec,r0,bt0,ra
c
cyr99      if (nbfld .eq. 3) then
cyr99        psidif = psir(1)-psir(nj)
cyr99        do 22 j=1,nj
cyr99   22   psinrm(j) = (psir(1)-psir(j))/psidif
cyr99      else
c
      if (nbfld .eq. 1) then
        do j=1,nj
          psinrm(j) = r(j) / ra   !note ra is defined above as r(nj)
        end do
      end if
c
      write (nscr, 1006)  (psinrm(j), j=1,nj)
      write (nscr, 1006)  (zeff  (j), j=1,nj)
      write (nscr, 1006)  (ene   (j), j=1,nj)
      write (nscr, 1006)  (te    (j), j=1,nj)
****  write (nscr, 1006)  (q     (j), j=1,nj)
      call giveupus(nscr)
      close (unit = nscr)
c
c -------- begin TORAY section -----------------------------------------
c
c       program = 'toray'

       program =toray_to_run(1:len_toray_to_run)
      itoray  = itoray + 1
      write (toray_log, 1002)  'toray_', itoray, '.log'
c
      if (first_time) then
        first_time = .false.
        if      (rf_output .eq. 0) then
          tlogfile = '/dev/null'
        else if (rf_output .eq. 1) then
          tlogfile = ' '
        else if (rf_output .eq. 2) then
          tlogfile = 'toray.log'
        else if (rf_output .eq. 3) then
          tlogfile = 'toray.log'
          call DESTROY (tlogfile)
        else if (rf_output .eq. 4) then
          tlogfile =  toray_log
        end if
      else if (rf_output .eq. 4) then
        tlogfile  = toray_log
      end if
c
      if (rf_output .eq. 1) then
        write  (ncrt , 1001)  program(1:len_toray_to_run)
        write  (nitre, 1001)  program(1:len_toray_to_run)
      else
        write  (ncrt , 1000)  program(1:len_toray_to_run), 
     .                         tlogfile(1:LENGTH(tlogfile))
        write  (nitre, 1000)  program(1:len_toray_to_run),
     .                         tlogfile(1:LENGTH(tlogfile))
      end if
c
      if (ISHELL (program(1:len_toray_to_run)
     .                            // operator // tlogfile) .ne. 0)
     .  call STOP ('subroutine ECH: failure of spawned TORAY', 65)
c
      write  (ncrt , 1003)  program(1:len_toray_to_run)
      write  (nitre, 1003)  program(1:len_toray_to_run)
c
c -------- end TORAY section -------------------------------------------
c

c ---------Rename toray.nc to toray_time.nc---------------------------------------------
      ech_filename = 'toray.nc'
 

!      CALL append_time_to_filename(time,ech_filename) ! creates new name stored in module
       CALL append_gyro_id_to_filename(time,ech_filename) ! creates new name stored in module

      CALL rename_disk_file(ech_filename) !new_name,toray_XXXXX.nc  was loaded in 
c                                          append_time_to_filename
c                                          file_proc.f90 contains more info
c                                          Assumes current working directory is proper one

c ---------Rename toray.in to toray.in_xx_time---------------------------------------------
c   if toray.in was created we have to rename it so that the next ech channel gets 
c   the correct set of values defined in the echin.nc file:
c   the model id is taken from model_globl in echdat_module.f90
c   dont save toray.in, just delete it instead:
      IF(ech_input .NE. 'none')THEN 
        CALL delete_file('toray.in')
!       INQUIRE(FILE ='toray.in',EXIST  = ex)
!       IF(ex)THEN
!          torayinpt ='toray.in'
!          CALL append_gyro_id_to_filename(torayinpt)
!          CALL rename_disk_file(torayinpt)
!       ENDIF
      ENDIF

c ----------------------------------------------------------------------------------------

c ------- read file echout, created by toray ---------------------------
c
      call getioun(nscr,nscr)
      open (unit = nscr, file = 'echout', status = 'OLD')
c
c   Read data from TORAY
c     ledge              = number of radial bin boundaries
c     xbouni(1:ledge)    = normalized poloidal flux of bin boundaries
c                          from 0(center) to edge(1.0)
c     voltoray           = plasma volume (cm**3)
c     vb(1:ledge-1)      = bin volumes, as fraction of plasma volume
c     qrfeb(1:ledge-1)   = RF pwr dens in each bin, per unit incident power.
c     tpowde(1,ledge-1)  = integrated power starting at plasma center,
c                          for unit incident power.
c     wiecrt(1:ledge-1)  = amps per cm**3 in each bin *voltoray,
c                          per incident watt.
c     tidept(1:ledge-1)  = integrated current (amps),starting at plasma
c                          center, per incident watt.
c     (Probably correct, but should check with Gary Smith (Bob H. 7/26/88))
c
c     (Some data is passed for cross-checking.)
c
 790  format (5e16.9)
 800  format (i4)
      read   (nscr, 800)  ledge
      if (ledge .gt. ktoray) then
        call giveupus(nscr)
        call STOP ('subroutine ECH: ledge > ktoray', 52)
      end if
      read   (nscr, 790)  (xbouni(j), j=1,ledge)
      read   (nscr, 790)   voltoray
      read   (nscr, 790)  (vb(j)    , j=1,ledge-1)
      read   (nscr, 790)  (qrfeb(j) , j=1,ledge-1)
      read   (nscr, 790)  (tpowde(j), j=1,ledge-1)
      read   (nscr, 790)  (wiecrt(j), j=1,ledge-1)
      read   (nscr, 790)  (tidept(j), j=1,ledge-1)


c     new values in echout, read by Onetwo as of version 3.94, HSJ 02/14/06
c     h_factr < sqrt(1-B/Bmax) >   ( Not used in Onetwo )
c     bsq_avg <B^2/B0^2>           ( Not used in Onetwo )
c     b_avg  <B/B0>                ( Not used in Onetwo )
c     r0rinv <R0/R>                ( Not used in Onetwo )
c     currf_toray <J dot B/Bt0> accounts for  total power input 
c                                  eq not per incident power
c     rjpdrho = < j_parallel > per incident watt ( Not used in Onetwo )
      IF( .NOT. ALLOCATED(currf_toray))ALLOCATE(currf_toray(ledge))
      IF( .NOT. ALLOCATED(dumy))ALLOCATE(dumy(ledge))
      read   (nscr, 790)  (dumy(j), j=1,ledge)    !h_factr
      read   (nscr, 790)  (dumy(j), j=1,ledge)    !bsq_avg
      read   (nscr, 790)  (dumy(j), j=1,ledge)    !b_avg
      read   (nscr, 790)  (dumy(j), j=1,ledge)    !r0rinv
      read   (nscr, 790)  (currf_toray(j), j=1,ledge)
      read   (nscr, 790)  (jtor_rf(j), j=1,ledge)    !rjpdrho
!          print *,'ledge,voltoray =',ledge,voltoray
!          print *,'xbouni =',xbouni
!          print *,'currf_toray =',currf_toray
!          print *,'qrfeb =', qrfeb
!          print *,'tpowde =',tpowde
!          print *,'tidept =',tidept
!           print *,'read rjpdrhop =',jtor_rf
!          do j=1,nj
!             print *,'j,r(j),psinrm(j) =',j,r(j),psinrm(j)
!          enddo
!        call STOP('temp,line 1828,ech',1)
      call giveupus(nscr)
      close  (unit = nscr)
      write  (6, '(a)')  ' closed file echout'
c
c     Convert wiecrt to amp/cm**2
c       darea = vb(j)*<rmajor/R>/(2 * pi * rmajor)
c       We approximate the flux surface average <rmajor/R> = 1.0
c
      if (tpowde(ledge-1) .le. 1.0e-4)  go to 90  ! if no accumulated deposited power skip out
      pi = ATAN2 (0.0, -1.0) 
      if (nbfld .eq. 1)  xmaxis = r0 / 100.0      !note r0 = rmajor
      write (6, '(a, f14.2)')  ' voltoray = ', voltoray
      write (6, '(a, f14.2)')  ' mag axis toray = ', 100.*xmaxis 
      do j=1,ledge-1
        wiecrt(j) = wiecrt(j)*(2.0 * pi * 100.0*xmaxis) / voltoray
      end do
c
c     Interpolate qrfeb and wiecrt onto ONETWO mesh.
c     qrfeb and wiecrt are given at mesh centers.
c     Use MESCON to convert to mesh boundaries.
c
      do j=1,ledge-1 
        dxbouni(j) = xbouni(j+1)-xbouni(j)
      end do
      call mescon(qrfeb,dxbouni,ledge)
      do j=1,ledge
        if (qrfeb(j) .lt. 0.0)  qrfeb(j) = 0.0
      end do
      call mescon (wiecrt, dxbouni, ledge)
c
c     Set boundary values
c
      qrfe(1) = qrfeb(1)
      qrfe(nj) = qrfeb(ledge)
!      currf(1) = wiecrt(1)
!      currf(nj) = wiecrt(ledge)
!      if (wiecrt(ledge-1)*currf(nj) .lt. 0.0)  currf(nj) = 0.0

      !NEW for toray version 1.8 convert currf_toray to cell edges:
      !note: ledge may not equal nj:
      call mescon(currf_toray,dxbouni,ledge)
      currf(1)  = currf_toray(1)
      currf(nj) = currf_toray(ledge)

      dumy(:) = jtor_rf(:)
      Call mescon(dumy,dxbouni,ledge)
      jtor_rf(1) = dumy(1) 
      jtor_rf(nj) = dumy(ledge)
      


c
c     Fill in rest of mesh (the normalized, nonuniform psi mesh in xbouni
c     does  should  match the onetwo normalized psi mesh,psinrm. But 
c     its possible that toray has recalcualted it so to be sure
c     we do the interpolation:
      do j=2,nj-1
        call find(j1,j2,psinrm(j),xbouni,ledge)
        if (j1 .ne. j2) then
          del = (psinrm(j)-xbouni(j1))/(xbouni(j2)-xbouni(j1))
        else
          del =  0.0
        end if
        qrfe (j) = (qrfeb (j1) + (qrfeb (j2) - qrfeb (j1)) * del)
c        currf(j) = (wiecrt(j1) + (wiecrt(j2) - wiecrt(j1)) * del)
         currf(j) = currf_toray(j1) +(currf_toray(j2)-
     .                                        currf_toray(j1))*del
         jtor_rf(j) = dumy(j1) +(dumy(j2)-dumy(j1))*del
      end do
c

      call FLUSH (ncrt, iflush)
c
c     doing the profile smoothing using SMTHSJ (HSJ and YRLL, 9 Jun 94)
c     conserve the total power and current (necsmth is default 0
c     (meaning no smoothing) otherwise set by user input.  
c
      call smthsj (r, qrfe , nj, necsmth)
      call smthsj (r, currf, nj, necsmth)
c
      ptot = tpowde(ledge-1)
      ctot = tidept(ledge-1)
c
c     integrate profiles and compare with TORAY results
c
      call trapv (r, qrfe , hcap, nj, ptotn)
      call trapv (r, currf, hcap, nj, ctotn)
****  ptotn = 4.0 * pi**2 * 100.0 * xmaxis * ptotn
      ptotn = 4.0 * pi**2 * 100.0 * rmajor * ptotn  ! 31 Oct 94 HSJ/YRLL
      ctotn = 2.0 * pi    * ctotn
      write  (ncrt, '("ledge = ",i5," ptot = ",e16.9,
     .        " ptotn = ",e16.9,/," ctot = ",e16.9,
     .        "ctotn = ",e16.9 )')  ledge, ptot, ptotn, ctot, ctotn
c
c     renormalize power and current densities
c


c       changed the following /03/10/03 HSJ:
c      if (ptotn .eq. 0.0)  ptotn = 1.0e-10
c      if (ctotn .eq. 0.0)  ctotn = 1.0e-12
c       changed the following /11/07/03 HSJ:
      if (ptotn .eq. 0.0)  ptotn = 1.0
      if (ctotn .eq. 0.0)  ctotn = 1.0

       if(ptotn .eq. 0.0 .or. ctotn .eq. 0)then
          write(ncrt,FMT = '("rmajor =",1pe12.4)')rmajor
          call STOP('ech, ptotn  or ctotn =0',0)
       endif
c
      do j=1,nj
        qrfe (j)  = ptot/ptotn*qrfe (j)*rfpow
!        currf(j) = ctot/ctotn*currf(j)*rfpow
        currf(j)  = currf(j)*rfpow !changed 04/14/06 per Prater recom., HSJ
      end do
      print *,'ptot,ptotn,rfpow =',ptot,ptotn,rfpow
c      call trapv (r, qrfe , hcap, nj, ptotnn)
c      print *,'ptotnn =',ptotnn
c      ptotnnn = 4.0 * pi**2 * 100.0 * rmajor * ptotnn 
c      print *,'rmajor ,ptotnnn =',rmajor,ptotnnn
c
      jtor_rf(:) = jtor_rf(:)*rfpow

      totecpe = ptot * rfpow
      totecc  = ctot * rfpow
c
      if (iplotec .eq. 0)  go to 2
****  if (te(1) .ne. told .and. ABS (te(1)-told)/te(1) .lt. 0.2) go to 2
c
c --- generate sequential names for XPLOT's CGM output files
c --- to avoid file clobberage
c
      icount = icount + 1
      write  (cgmfile, '(a6, i2.2, a4)')  'xplot_', icount, '.cgm'
      write  (ncrt, 70)  cgmfile
   70 format (/ ' ---- XPLOT started, plots sent to CGM file "', a, '"')
c
c --- run XPLOT via shell script "run_xplot" to allow remote execution
c --- this change was implemented 17 Nov 98 for Tim Luce, by Joe Freeman
c
c      if (ISHELL ('run_xplot ' // cgmfile) .ne. 0)
c     .  call STOP ('subroutine ECH: failure of spawned XPLOT', 66)


c     xplot is not available everywhere, let the command fail if necessary HSJ
       j = ISHELL('run_xplot ' // cgmfile)
c
      write  (ncrt, 80)
   80 format (  ' ---- XPLOT finished' /)
c
    2 told = te(1)

      return
c
   90 do j=1,nj
        qrfe (j) = 0.0
        currf(j) = 0.0
      end do
      totecpe    = 0.0
      totecc     = 0.0
      told       = te(1)
      return
c
      end

      subroutine eqdskrt
c
      USE param
      USE io 
      USE mhdpar
      USE ename
      USE extra
      USE numbrs
      USE machin
      USE rhog
      USE gpsi
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- This replaces an old version of eqdskrt  (BH  10/7/91)
c ----------------------------------------------------------------------
c
c    Make an extended eqdsk.  Used with ray-tracing code.
c    Use the latest eqdsk available.
c
c    The following INCLUDE files are for purposes of obtaining 2D eqdsk data:
c
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'numbrs.i'
c      include 'io.i'
c      include 'machin.i'
c      include 'small.i'
c      include 'rhog.i'
c      include 'extra.i'
c      include 'ename.i'
c
      character  ntitle(5)*8, eqdsk*64
      dimension  fdum(5)
c
c     no MHD calculations done, only eqdsk available is eqdskin
c
      eqdsk = eqdskin
c
c     get most recent eqdsk. this is not eqdskin because we have
c     recalculated equilibria and generated new eqdsks
c
      if (eqdskfilename .ne. 'none')  eqdsk = eqdskfilename
c
      write  (6, 36)  eqdskfilename, eqdskin, eqdsk
   36 format (/ ' eqdskfilename: ', a /
     .          '       eqdskin: ', a /
     .          '         eqdsk: ', a )
c
c --- open eqdsk file and create eqdskex (extended eqdsk)
c    (name used to be eqdsk.ex but was changed to eqdskex by JF for portability)
c
      iguess = 19 ! myopen calls getioun, open, and giveupus if open failed
      call myopen (iguess, eqdsk, 2, length, iread)
c
      if (iread .ne. 1) then
        write  (ncrt, 81)
        write  (nout, 81)
   81   format (' FATAL ERROR: could not open guess file')
        call STOP ('subroutine EQDSKRT: could not open guess file', 54)
      end if
c
      call getioun(nscr,nscr)
      open (unit = nscr, file = 'eqdskex', status = 'UNKNOWN')
c
c     Transcribe eqdsk to eqdskex
c
      read   (iguess, 89)  (ntitle(i), i=1,5), dat, ipestg, nxeqd, nyeqd
      write  (nscr  , 89)  (ntitle(i), i=1,5), dat, ipestg, nxeqd, nyeqd
   89 format (6a8, 3i4)
c
      do 10 ii=1,4
        read   (iguess, 82)  (fdum(i), i=1,5)
   10   write  (nscr  , 82)  (fdum(i), i=1,5)
   82   format (5e16.9)
c
      nlinesx  = (nxeqd         - 1) / 5 + 1
      nlinesxy = (nxeqd * nyeqd - 1) / 5 + 1
c
      do 20 jj=1,4
      do 20 ii=1,nlinesx
        read  (iguess, 82)  (fdum(i),i=1,5)
   20   write (nscr  , 82)  (fdum(i),i=1,5)
c
      do ii=1,nlinesxy
        read  (iguess, 82)  (fdum(i), i=1,5)
        write (nscr  , 82)  (fdum(i), i=1,5)
      end do
c
      do ii=1,nlinesx
        read  (iguess, 82)  (fdum(i), i=1,5)
        write (nscr  , 82)  (fdum(i), i=1,5)
      end do
c
      call giveupus(iguess)
      close (unit = iguess)
c
c add on eqdsk extension
c
      call eqdskext (nscr)
c
      call giveupus(nscr)
      close (unit = nscr)
      return
c
      end

      subroutine fcreal (ar, zr, zi, dzr, dzi, n, rem, tol)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      in    = 0
      sqrpi = 1.772454
      x     = ABS (ar)
      if (x .le. 1.5)  go to 11
      if (x .le. 4.0)  go to 12
   14 call asyr(x,zr,dzr,n,rem,tol)
      in  = 1
   31 asq = -x**2
      if (asq .lt. -8.0e1)  go to 81
      zi  = sqrpi * EXP (asq)
      if (in .ne. 1)  go to 82
      dzi = -2.0*x*zi
      go to 20
   82 dzr = -2.0 * (1.0 + x*zr)
      dzi = -2.0 * x * zi
      go to 20
   81 zi = 0.0
      dzi = 0.0
   20 if (ar .gt. 0.0)  return
      zr = -zr
      dzi = -dzi
      return
   11 call powr(x,zr,n,rem,tol)
      go to 31
   12 call intr(x,zr,n,rem,tol)
      go to 31
c
      end

      subroutine fischb (ifb, lmode, rfpow, rfcur, r, nj, ene, te,
     .                   zeff, rnp, rmajor, hcap, iprof, alpha, x,
     .                   qrfe, currf, itrapech)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension  alpha(3), currf(nj), ene(nj), hcap(nj), qrfe(nj),
     .           r(nj), te(nj), x(nj), zeff(nj)
c
c     rfpow   watts
c     rfcur  amps
c     qrfe   w/cm**3
c     qrfi   w/cm**3
c     currf  amps/cm**2
c     ene    cm**-3
c     te     keV,  (but not explicitly needed below).
c     ifb    1, specify total power and calc profile and associated RF cur.
c            2, specity total current and calc prof and associated power.
c     iprof  1,2, or 3.  a profile disignator.
c     lmode  0, LH current drive efficency.
c            1,2, or 3.  ECH lmode-harmonic efficency.
c
c     fjdp0 is in amp*meters/watt, as in cordey, edlington and start,
c             plasma physics, 1982 (i.e. 4 times fisch-boozer normalization).
c
      fjdp0(j) = 4.0 * 1.28e12 * te(j) / ene(j)
      c        = 3.0e8
c
      if (iprof .gt. 3  .or.  lmode .gt. 3)
     .  call STOP ('subroutine FISCHB: unspecified problem', 57)
c
c     profile calculation
c     make RF profile zero at outer nb mesh points (nb .ge. 1)
c
      nb = 4
c
c     introduce smoothing of profiles.
c     ns = the number of 3-point smoothing operations performed
c
      ns  = 5
      njb = nj - nb + 1
      go to (25, 50, 75), iprof
c
c     parabolic profile
c
   25 do 30 j=1,njb
   30 x(j) = (1.0 - (r(j)*r(j)/r(njb)/r(njb))**alpha(2))**alpha(1)
      if (njb .eq. nj)  go to 32
      do 31 j=njb+1,nj
      x(j) = 0.0
   31 continue
   32 continue
      go to 100
c
c     gaussian profile
c
   50 r0   = r(njb)*alpha(1)
      x01  = r(njb)*r(njb)/r0/r0
      do 55 j=1,njb
      x02  = r(j)*r(j)/r0/r0
   55 x(j) = EXP (-x02) - EXP (-x01)
      if (njb .eq. nj)  go to 57
      do 56 j=njb+1,nj
      x(j) = 0.0
   56 continue
   57 continue
      go to 100
c
c     modified maxwellian - hollow profile
c
   75 r0  = alpha(1)*r(njb)
      r1  = alpha(2)*r(njb)
      x02 = (r(njb) - r1)*(r(njb) - r1)/r0/r0
      x01 = r(njb)*r(njb)/r0/r0
c
c --- the following corrections were copied from cray33.f  19 Jan 93 HSJ
c
****  xedge = x01 * EXP (-x02)                             **bh** 930118
      xedge =       EXP (-x02)
      n = INT (alpha(3)*njb) + 1
      do 90 j=n,njb
****    x01 =  r(j)*r(j)/r0/r0                             **bh** 930118
        x02 = (r(j) - r1) * (r(j) - r1) / r0 / r0
***90   x(j) = x01 * EXP (-x02) - xedge                    **bh** 930118
   90   x(j) =       EXP (-x02) - xedge
      n1 = n
      if (n .gt. 1) n1 = n - 1
      do 95 j=1,n1
   95 x(j) = x(n)
      if (njb .eq. nj)  go to 100
      do j=njb+1,nj
        x(j) = 0.0
      end do
c
c     smoothing.  use currf for temporary storage
c
  100 if (ns .eq. 0)  go to 101
c
      do is=1,ns
        do j=2,nj-1
          currf(j) = 0.25*(x(j-1)+2.0*x(j)+x(j+1))
        end do
        do j=2,nj-1
          x(j) = currf(j)
        end do
      end do
c
c     integrate profile to determine normalization
c     xtot = integral(dr) from 0 to r(nj) of r*hcap*x
c
  101 call trapv (r, x, hcap, nj, xtot)
      x02 = 0.5/3.1416/xtot
      do 105 j=1,nj
  105 x(j) = x02*x(j)
c
      if (ifb .eq. 2)  go to 125
c
c     have been given total dissipated power,rfpow; want to calculate
c     heating profile, qrfe, current profile, currf, and total current, rfcur
c
c     current density in amps/cm**2
c
      do j=1,nj
        qrfe(j) = rfpow*x(j)/6.2832/rmajor
        if (rnp .le. 0.0)  u0 = ABS (rnp)
        if (rnp .gt. 0.0)  u0 = c/rnp/(1.325e7 * SQRT (te(j)))
        fjdp1    = fjdp0(j)*fjdp(u0,lmode,zeff(j))
        aspctinv = r(j)/rmajor
        trapfrc  = 1.0
        if (itrapech .eq. 1)  trapfrc = (1.0-aspctinv**0.7)**0.73
        currf(j) = trapfrc*qrfe(j)*fjdp1*10**2
      end do
c
      call trapv (r, currf, hcap, nj, rfcur)
      rfcur    = rfcur*6.2832
      go to 150
c
c     have been given total current,rfcur; want to calculate current profile,
c     currf, heating profile, qrfe, and total dissipated power, rfpow.
c
  125 do 135 j=1,nj
      aspctinv = r(j)/rmajor
      trapfrc = 1.0
      if (itrapech .eq. 1)  trapfrc = (1.0 - aspctinv**0.7)**0.73
      currf(j) = trapfrc*rfcur*x(j)
      if (rnp .le. 0.0)  u0 = ABS (rnp)
      if (rnp .gt. 0.0)  u0 = c/rnp/(1.325e7 * SQRT (te(j)))
      fjdp1 = fjdp0(j)*fjdp(u0,lmode,zeff(j))
  135 qrfe(j) = currf(j)/(fjdp1*10**2)
      call trapv (r,qrfe,hcap,nj,rfpow)
      rfpow = rfpow*rmajor*39.4784
c
  150 return
c
      end

      real*8 function fjdp (u0, lmode, zeff)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c calculate normalized j/p from Cordey et al.
c
      save
      dimension alphal(5,3), betal(5,3), z(5)
      data      alphal/1.5, 1.4, 0.65, 0.37, 0.2, 1.05, 1.0, 0.66,
     .                 0.4, 0.21, 1.4, 1.0, 0.7, 0.45, 0.27/,
     .           betal/0.06, 0.13, 0.14, 0.13, 0.09, 0.14, 0.17, 0.17,
     .                 0.13, 0.1, 0.16, 0.19, 0.17, 0.13, 0.08/,
     .               z/1.0, 2.0, 4.0, 8.0, 16.0/
c
c lmode = 0, LH
c
c 1-3, harmonic of ECH
c
      if (lmode .eq. 0  )  go to 30
      if (zeff  .lt. 1.0)
     .  call STOP ('function FJDP: unspecified problem', 58)
c
      do j=2,5
        jz = j-1
        if (zeff .lt. z(j))  go to 20
      end do
c
   20 fjdp = alphal(jz,lmode)*u0+betal(jz,lmode)*u0**2
     .       + (zeff-z(jz))/(z(jz+1)-z(jz))*
     .         ((alphal(jz+1,lmode)-alphal(jz,lmode))*u0
     .       +(betal(jz+1,lmode)-betal(jz,lmode))*u0**2)
      return
c
   30 fjdp = 1.2 * (6.0/(zeff+5.0))*(1.4/u0+0.7*u0**2)
      return
c
      end





      subroutine gafit_write(ipsi,ntoray,runid)
c------------------------------------------------------------------------
c
c     write a simplified version of gafit.in suitable
c     for use in Onetwo/toray:
c
c-------------------------------------------------------------HSJ-/03/06/03
       implicit none
       integer ipsi,ntoray
       character *(*) runid
         namelist /fitdat/ ipsi
         open (unit = ntoray, file = 'gafit.in', status = 'NEW')
         write (unit = ntoray, fmt = '(3x, a72)') runid
         write (unit = ntoray, nml = fitdat)
         close (unit = ntoray)
       return
       end
      subroutine get_newpsir (psimax, psilim, psiout, psir, j0, nj,
     .                        gasep, psirnew, xbouni, ltab, njnew)
c
      implicit none
c
      real*8    psimax, psilim, psiout, psir(*), gasep
      real*8    psirnew(*), xbouni(*)
      integer j0, nj, njnew, ltab(*)
c
      integer n, j
c
      psirnew(j0-1)=psimax
      xbouni(j0)=0.
      ltab(j0)=1
c
      n=j0
      do j=j0,nj
         xbouni(n+1)=(psir(j)-psimax)/(psilim-psimax)
         if (xbouni(n+1)-xbouni(n) .ge. gasep) then
            psirnew(n)=psir(j)
            ltab(n+1)=j
            n=n+1
         end if
      end do
      njnew=n-1
c
      if (psirnew(njnew) .ne. psiout) then
         psirnew(njnew)=psiout
         ltab(njnew+1)=nj
         xbouni(njnew+1)=(psiout-psimax)/(psilim-psimax)
      end if
      return
c
      end

      subroutine icherr (ierr)
c
      USE param
      USE io 
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- ICHERR handles error conditions for subroutine ICH
c
c      include 'param.i'
c      include 'io.i'
c
      if (ierr .ne. 1)  go to 102
      write  (nout, 8001)
      write  (nqik, 8001)
      write  (ncrt, 8001)
 8001 format (/' ICH error: RF matrix is singular')
      call STOP ('subroutine ICHERR: singular RF matrix', 55)
  102 if (ierr .ne. 2)  go to 103
      write  (nout, 8002)
      write  (nqik, 8002)
      write  (ncrt, 8002)
 8002 format (/ ' ICH error: array dimensions are too small')
      call STOP ('subroutine ICHERR: array dimensions too small', 56)
c
  103 return
c
      end








      subroutine ich (freq0,xkpar0,nhigh0,ykprp0,iside,navg,ichmod,
     .                betalm,nrfrad,rfrad,rfrow,jrfmin,jrfmax,kj,nj,r,
     .                ene,en,te,ti,nion,atw,zion,r0,bt0,hcap,rfpow,qrfe,
     .                qrfi)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          ich                      t.k.mau     9/30/82
c         (ONETWO interface:        b.stockdale)
c
c          Ion cyclotron heating model.
c
c          inputs:
c          freq0   - frequency of applied wave (1/sec)
c          xkpar0  - k-parallel of applied wave (1/cm)
c          nhigh0  - number of terms retained in K
c          ykprp0  - ykperp
c          iside   - +1 for high-field launch (inside)
c                    +2  or -1 for low-field launch (outside)
c          navg    - number of adjacent points to be used in moving
c                    average of computed power absorbtion profiles
c          ichmod  - 1 to use t.k.mau code
c                  - 2 to use s.c.chiu code
c                  - 3 to use t.k.mau code when betai0 < betalm, and
c                    use s.c.chiu code when betai0 >= betalm
c          betalm  - limit on betai0 when ichmod = 3
c          nrfrad  - number of points in rfrad
c          rfrad(i)- major radius grid points (cm) to be used
c          rfrow(i)- rho value at rfrad(i) (cm)
c          jrfmin  - rfrad(jrfmin) is innermost point to use
c          jrfmax  - rfrad(jrfmax) is outermost point to use
c          kj      - first dimension of transport code arrays
c          nj      - number of points in transport code mesh
c          r(i)    - radial grid point i in transport mesh (cm)
c          ene(i)  - electron density at r(i)  (1/cm**3)
c          en(i,j) - ion species j density at r(i)  (1/cm**3)
c          te(i)   - electron temperature at r(i)  (keV)
c          ti(i)   - ion temperature at r(i)  (keV)
c          nion    - number of ion species (up to 3)
c          atw(j)  - atomic weight of ion species j (h = 1.0)
c          zion(i,j)-charge of ion species j at r(i)
c          r0      - major radius of transport grid axis (cm)
c          bt0     - toroidal field at r0  (Gauss)
c          hcap(i) - geometric factor H from transport code
c          rfpow   - total RF power launched (W)
c
c          outputs:
c          qrfe(i) - ich power source density  (keV/cm**3-sec)
c          qrfi(i) - ich power source density  (keV/cm**3-sec)
c
      dimension rfrad(nrfrad),rfrow(nrfrad)
      dimension r(nj),ene(nj),en(kj,nion),te(nj),ti(nj)
      dimension atw(nion),zion(kj,nion),hcap(nj)
      dimension qrfe(nj),qrfi(nj)
c
c          nx,ny must be same as in trpic2. ny <= kich.
c
      include 'ichp2.i'
      parameter (ikj = 101, nx = 2)
      include 'ichpar.i'
      include 'ichcom.i'
c
      dimension prfe(kich),prfi(kich)
c
      common /rf1/ dene(kich),den1(kich),den2(kich),den3(kich),
     .             terf(kich),tirf(kich),brf(kich)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      common /nn / nstart,i3,j1
c
      complex*16 abd(16,kich4),work(kich4),bcol(kich4)
      integer ipvt(kich4)
      dimension dvol(ikj),qrfrad(kich),work1(ikj)
      dimension rrow(ny),dni(nx,ny),pabs(nx,ny),pabsie(ny)
c
c check that dimensions are adequate
c
      if (nrfrad .gt. kich)  call icherr(2)
      if (  nion .gt. ks-1)  call icherr(2)
      if (    nj .gt. ikj )  call icherr(2)
c
c initialize
c
      freq     = freq0
      xkpar    = xkpar0
      nhigh    = nhigh0
      ykperp   = ykprp0
      kinc     = iside
      if (kinc .eq. -1)
     .kinc     = 2
      ngrid    = jrfmax - jrfmin + 1
      neq      = nion  + 1
      nhp      = nhigh + 1
      nh2p     = nhigh*2 + 1
      pi       = 3.14159
      wfreq    = 2.0 * pi * freq
      clight   = 2.998e10
      ergtkev  = 1.6022e-9
      charge   = 4.8032e-10
      pmass    = 1.6726e-24
      anumb(1) = 1.0
      emass(1) = 9.1095e-28
      do k=2,neq
        anumb(k) = zion(1,k-1)
        emass(k) = atw(k-1)*pmass
      end do
c
c decide which ich code to call
c
      icode = 1
      if (ichmod .eq. 1)  go to 90
      if (ichmod .eq. 2)  icode = 2
      if (ichmod .eq. 2)  go to 90
      betai0 = 0.0
      do 80 k=1,nion
   80 betai0 = betai0 + en(1,k)*ti(1)
      betai0 = betai0*8.0 * pi * 1.6e-9 / (bt0**2)
      if (betai0 .ge. betalm)  icode = 2
   90 if (icode .eq. 2 .and. ngrid .gt. 2*ny)  call icherr(2)
      if (icode .eq. 2 .and.  nion .gt.   nx)  call icherr(2)
c
c left-justify rfrad grid in rmrf.
c obtain plasma profiles at rmrf grid points by linear interpolation.
c
      do 105 i=1,ngrid
        if (kinc .eq. 1) jj = jrfmin + i - 1
        if (kinc .eq. 2) jj = jrfmax - i + 1
        if (icode .eq. 1) ii = i
        if (icode .eq. 1)  go to 101
        if (ngrid .le. ny) ii = i
        if (ngrid .le. ny)  go to 101
        if (i .ne. ((i-1)/2)*2+1)  go to 105
        ii = (i-1)/2 + 1
  101   rmrf(ii) = rfrad(jj)
        if (icode .eq. 2) rrow(ii) = rfrow(jj)
        call interp(rfrow(jj),r,nj,ene,dene(ii))
        call interp(rfrow(jj),r,nj,en,den1(ii))
        if (nion .ge. 2)
     .     call interp(rfrow(jj),r,nj,en(1,2),den2(ii))
        if (nion .ge. 3)
     .     call interp(rfrow(jj),r,nj,en(1,3),den3(ii))
        call interp(rfrow(jj),r,nj,te,terf(ii))
        call interp(rfrow(jj),r,nj,ti,tirf(ii))
        brf(ii) = bt0*r0/rmrf(ii)
  105 continue
c
      if (icode .eq. 2)  go to 180
c
c          solve for ich power deposition profile
c
      nstart = ngrid/2
      call rfheat(nstart)
c
      do j=nstart-1,1,-1
        call rfheat(j)
      end do
c
      do jj=nstart+1,ngrid
        call rfheat(jj)
      end do
c
      do i=1,ngrid
        call param0(i)
      end do
c
      n41 = 4*(ngrid-1)
      i3 = n41-3
      j1 = n41-1
      call rfwave(n41,abd,bcol,ipvt,work,icheck)
      if (icheck .ne. 0)  call icherr(1)
      call absorb(ngrid,prfe,prfi)
      call rflux(ngrid,fabs)
      go to 200
c
c call s.c.chiu model when absorption is nearly 100%
c the other code has trouble at high beta
c
  180 nslb = ii
      do 185 i=1,nslb
      dni(1,i) = den1(i)
      dni(2,i) = den2(i)
      tirf(i) = tirf(i)*1.0e3
  185 terf(i) = terf(i)*1.0e3
      x0 = rmrf(1)-r0
      xf = rmrf(nslb)-r0
      freqr = 2.0 * pi * freq
      call trpic2(dni,anumb(2),atw,tirf,terf,brf,
     .            x0,xf,nion,nslb,xkpar,freqr,pabs,prfi,prfe,
     .            pabsie,fabs,pleft,pitot,petot)
c
c          compute volume of flux shells and
c          convert to power source density for transport code
c
  200 fact = 4.0 * pi**2 * r0
      do 205 i=2,nj-1
  205 dvol( i) = fact*hcap( i)*r(i)*0.5*(r(i+1)-r(i-1))
      dvol( 1) = fact*hcap( 1)*0.125*r(2)**2
      dvol(nj) = fact*hcap(nj)*0.125*(3.0*r(nj)+r(nj-1))*(r(nj)-r(nj-1))
      if (icode .eq. 2)  go to 250
      do 210 i=1,ngrid
      if (kinc .eq. 1) jj = jrfmin + i - 1
      if (kinc .eq. 2) jj = jrfmax - i + 1
  210 qrfrad(jj) = prfe(i)
      icentr = 1
      call qcon(icentr,rfrad,rfrow,jrfmin,jrfmax,qrfrad,
     .          r0,nj,r,dvol,qrfe)
      do 220 i=1,ngrid
      if (kinc .eq. 1) jj = jrfmin + i - 1
      if (kinc .eq. 2) jj = jrfmax - i + 1
  220 qrfrad(jj) = prfi(i)
      icentr = 1
      call qcon(icentr,rfrad,rfrow,jrfmin,jrfmax,qrfrad,
     .          r0,nj,r,dvol,qrfi)
      go to 280
c
c          do conversion for trpic2 case
c
  250 if (kinc .eq. 2)  call switch(rmrf,nslb,work)
      if (kinc .eq. 2)  call switch(rrow,nslb,work)
      do i=1,nslb
        if (kinc .eq. 1) jj = i
        if (kinc .eq. 2) jj = nslb - i + 1
        qrfrad(jj) = prfe(i)
      end do
c
      icentr = 1
      call qcon (icentr,rmrf,rrow,1,nslb,qrfrad,r0,nj,r,dvol,qrfe)
c
      do i=1,nslb
        if (kinc .eq. 1) jj = i
        if (kinc .eq. 2) jj = nslb - i + 1
        qrfrad(jj) = prfi(i)
      end do
c
      icentr = 1
      call qcon(icentr,rmrf,rrow,1,nslb,qrfrad,r0,nj,r,dvol,qrfi)
c
c "post-processing": apply moving average if navg > 0
c
  280 do 290 j=1,nj
      if (qrfe(j) .lt. 0.0) qrfe(j) = 0.0
  290 if (qrfi(j) .lt. 0.0) qrfi(j) = 0.0
      if (navg .eq. 0)  go to 300
      call qnavg(navg,nj,qrfe,work1)
      call qnavg(navg,nj,qrfi,work1)
c
c          normalize to total power
c
  300 ptot = 0.0
      do 310 j=1,nj
  310 ptot = ptot + (qrfe(j)+qrfi(j))*dvol(j)
      ptot = fabs*rfpow*0.625e16/ptot
      do 320 j=1,nj
      qrfe(j) = qrfe(j)*ptot
  320 qrfi(j) = qrfi(j)*ptot
      return
c
      end

      subroutine intqdr (rm1, rm2, rfrad, j1, j2,
     .                   qrf1, qrf2, qrfrad, pwr)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c Intqdr integrates the linear power density profile
c from rm1 to rm2, using the trapezoidal rule.
c rm1 and rm2 must lie between rfrad(j1) and rfrad(j2).
c qrfrad(i) is the linear power density at rfrad(i).
c qrf1 and qrf2 are the (given) interpolated linear power
c densities at rm1 and rm2, respectively.
c The output is pwr, the total power between rm1 and rm2.
c
      dimension  rfrad(*), qrfrad(*)
c
      jj = j1
  100 if (rfrad(jj) .ge. rm1)  go to 110
      jj = jj+1
      go to 100
c
c Now rfrad(jj) is the next grid point after rm1.
c Check for rm2 < rfrad(jj).
c
  110 if (rfrad(jj) .lt. rm2)  go to 120
      pwr = 0.5 * (qrf1+qrf2)*(rm2-rm1)
      return
c
c Begin with integral from rm1 to rfrad(jj).
c Add rfrad intervals until rfrad(jj) >= rm2.
c
  120 pwr = 0.5*(qrf1+qrfrad(jj))*(rfrad(jj)-rm1)
  130 jj = jj+1
      if (rfrad(jj) .ge. rm2)  go to 140
      pwr = pwr + 0.5*(qrfrad(jj-1)+qrfrad(jj))*(rfrad(jj)-rfrad(jj-1))
      if (jj .eq. j2)  return
      go to 130
c
c Add the integral from rfrad(jj-1) to rm2.
c
  140 pwr = pwr + 0.5*(qrfrad(jj-1)+qrf2)*(rm2-rfrad(jj-1))
      return
c
      end

      subroutine intr (x, zr, n, rem, tol)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension f(10)
c
      i1 = 0
      b  = 1.5
      call powr(b,f(1),n,rem,tol)
   35 n = 0
   13 f(2) = -2.0*(1.0+b*f(1))
      do 10 i=1,7
      x1 = FLOAT (i)
   10 f(i+2) = -2.0*x1*f(i)-2.0*b*f(i+1)
      h = (4.032e3*tol / ABS (f(9)))**0.125
      r = b + h
      if (r .lt. x)  go to 19
      h = x - b
      i1 = 1
   19 b = b + h
      p = 1.0
      do 12 i=2,8
      x1 = FLOAT (i-1)
      p = p*h/x1
   12 f(1) = f(1) + p*f(i)
      n = n + 1
      if (i1 .eq. 0)  go to 13
   14 zr = f(1)
      rem = 0.1 * tol * FLOAT (n) + rem
      return
c
      end

      complex*16 function kxx (xx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      complex*16    sxx   ,sxy   ,sxz   ,syy   ,syz   ,szz
      common /rf2 / sxx(5),sxy(5),sxz(5),syy(5),syz(5),szz(5)
      complex*16    zlam   ,besi     ,besd     ,fthet,fphi,xx,temp(5)
      common /rfzl/ zlam(5),besi(7,5),besd(7,5)
      common /rfz / zeta(11,5),ck(5),fthet(11,5),fphi(11,5)
      common /rf4 / xnparz,xnpar,dbg(99)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      kxx = CMPLX (0.0, 0.0)
      do 200 ii=1,neq
        sxx (ii) = CMPLX (0.0, 0.0)
        temp(ii) = CMPLX (0.0, 0.0)
  200 continue
      do 1 k=1,neq
      if (zlam(k) .eq. CMPLX (0.0, 0.0))  go to 222
      do 2 m=1,nhigh
        mn = m + nhigh
        temp(k) = ck(k) * EXP (-zlam(k)) * m * m * BESI (m,k) *
     .            (fthet(m,k)+fthet(mn,k))/zlam(k)
        kxx = kxx+temp(k)
        sxx(k) = sxx(k)+temp(k)
        dbg(29) =  REAL (kxx)
        dbg(30) = PIMAG (kxx)
    2 continue
c
      go to 103
  222 temp(k) = 0.5*ck(k)*(fthet(1,k)+fthet(nhp,k))
      kxx = kxx+temp(k)
      sxx(k) = sxx(k)+temp(k)
  103 kxx = kxx + CMPLX (1.0/neq,0.0)
      sxx(k) = sxx(k) + CMPLX (1.0/neq,0.0)
      dbg(31) =  REAL (kxx)
      dbg(32) = PIMAG (kxx)
    1 continue
c
      return
c
      end

      complex*16 function kxy (xx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      complex*16    xx,temp(5),zlam,besi,besd,fthet,fphi
      complex*16    sxx,sxy,sxz,syy,syz,szz
      common /rf2 / sxx(5),sxy(5),sxz(5),syy(5),syz(5),szz(5)
      common /rfzl/ zlam(5),besi(7,5),besd(7,5)
      common /rfz / zeta(11,5),ck(5),fthet(11,5),fphi(11,5)
      common /rf4 / xnparz,xnpar,dbg(99)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      kxy = CMPLX (0.0, 0.0)
      do 200 ii=1,neq
      sxy(ii) = CMPLX (0.0, 0.0)
  200 continue
      do 1 k=1,neq
      do 2 m=1,nhigh
      mn = m+nhigh
      go to 101
  101 temp(k) = CMPLX (0.0, 1.0)*ck(k) * EXP (-zlam(k))*m*((besd(m,k)
     .        - besi(m,k))*(fthet(m,k)-fthet(mn,k)))
      kxy = kxy+temp(k)
      sxy(k) = sxy(k)+temp(k)
    2 continue
      dbg(27) =  REAL (kxy)
      go to 102
  102 dbg(28) = PIMAG (kxy)
    1 continue
      return
c
      end

      complex*16 function kxz (xx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      complex*16 xx,temp(5),cxq,zlam,besi,besd,fthet,fphi
      complex*16 sxx,sxy,sxz,syy,syz,szz
      common /rf2 / sxx(5),sxy(5),sxz(5),syy(5),syz(5),szz(5)
      common /rfzl/ zlam(5),besi(7,5),besd(7,5)
      common /rfz / zeta(11,5),ck(5),fthet(11,5),fphi(11,5)
      common /rf4 / xnparz,xnpar,dbg(99)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      common /rffq/ wc(5),wp2(5),wp(5),vtx(5),vtz(5)
c
      kxz = CMPLX (0.0, 0.0)
      cxq = SQRT (xx)
      if (REAL (cxq) .lt. 0.0)  cxq = -cxq
      do 200 ii=1,neq
        sxz(ii) = CMPLX (0.0, 0.0)
  200 continue
      do 1 k=1,neq
      do 2 m=1,nhigh
      mn = m+nhigh
      go to 101
  101 if (zlam(k) .ne. CMPLX (0.0, 0.0))
     .  temp(k) = -ck(k) * EXP (-zlam(k))*
     .             cxq*m*besi(m,k)*(fphi(m,k)-fphi(mn,k))/zlam(k)/wc(k)
      kxz     = kxz    + temp(k)
      sxz(k)  = sxz(k) + temp(k)
      dbg(25) =  REAL (kxz)
      go to 102
  102 dbg(26) = PIMAG (kxz)
    2 continue
    1 continue
      return
c
      end

      complex*16 function kyy (xx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      complex*16    sxx,sxy,sxz,syy,syz,szz
      common /rf2 / sxx(5),sxy(5),sxz(5),syy(5),syz(5),szz(5)
      common /rfzl/ zlam(5),besi(7,5),besd(7,5)
      complex*16    fthet,fphi,xx,temp(5),zlam,besi,besd
      common /rfz / zeta(11,5),ck(5),fthet(11,5),fphi(11,5)
      common /rf4 / xnparz,xnpar,dbg(99)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      kyy = CMPLX (0.0, 0.0)
      do 200 ii=1,neq
      syy(ii) = CMPLX (0.0, 0.0)
  200 continue
      do 1 k=1,neq
      if (zlam(k) .eq. CMPLX (0.0, 0.0))  go to 231
      go to 101
  101 temp(k) = CMPLX (1.0/neq,0.0)-2.0*ck(k)*zlam(k) * EXP (-zlam(k))*
     .(besd(nhp,k)-besi(nhp,k))*fthet(nh2p,k)
      kyy = kyy+temp(k)
      syy(k) = syy(k)+temp(k)
      dbg(21) =  REAL (kyy)
      go to 102
  102 dbg(22) = PIMAG (kyy)
      do 2 m=1,nhigh
      mn = m+nhigh
      go to 103
  103 temp(k) = ck(k) * EXP (-zlam(k))*(m*m*besi
     .      (m,k)/zlam(k)-2.0*zlam(k)*(besd(m,k)-besi(m,k)))*(fthet(m,k)
     .     + fthet(mn,k))
      kyy = kyy+temp(k)
      syy(k) = syy(k)+temp(k)
    2 continue
      go to 230
  231 kyy = kyy + CMPLX (1.0/neq,0.0)+.5*ck(k)*(fthet(1,k)+fthet(nhp,k))
  230 dbg(23) =  REAL (kyy)
      go to 104
  104 dbg(24) = PIMAG (kyy)
    1 continue
      return
c
      end

      complex*16 function kyz (xx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      complex*16 xx,temp(5),cxq,fthet,fphi,zlam,besi,besd
      complex*16 sxx,sxy,sxz,syy,syz,szz
      common /rf2 / sxx(5),sxy(5),sxz(5),syy(5),syz(5),szz(5)
      common /rfzl/ zlam(5),besi(7,5),besd(7,5)
      common /rfz / zeta(11,5),ck(5),fthet(11,5),fphi(11,5)
      common /rf4 / xnparz,xnpar,dbg(99)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      common /rffq/ wc(5),wp2(5),wp(5),vtx(5),vtz(5)
c
      kyz = CMPLX (0.0, 0.0)
      cxq = SQRT (xx)
      if (REAL (cxq) .lt. 0.0)  cxq = -cxq
      do 200 ii=1,neq
      syz(ii) = CMPLX (0.0, 0.0)
  200 continue
      do 1 k=1,neq
      go to 101
  101 temp(k) = ck(k)*(EXP (-zlam(k))*cxq*(besd(nhp,k)-besi(nhp,k))
     .             *fphi(nh2p,k)/wc(k)) * CMPLX (0.0, 1.0)
      kyz = kyz+temp(k)
      syz(k) = syz(k)+temp(k)
      dbg(17) =  REAL (kyz)
      go to 102
  102 dbg(18) = PIMAG (kyz)
      do 2 m=1,nhigh
      mn = m+nhigh
      go to 103
  103 temp(k) = ck(k)*(EXP (-zlam(k)) * SQRT (xx)*(besd(m,k)-besi(m,k))
     .    * (fphi(m,k)+fphi(mn,k))/wc(k))*CMPLX (0.0, 1.0)
      kyz = kyz+temp(k)
      syz(k) = syz(k)+temp(k)
      dbg(19) =  REAL (kyz)
      go to 104
  104 dbg(20) = PIMAG (kyz)
    2 continue
    1 continue
      return
c
      end

      complex*16 function kzz (xx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      complex*16 xx,temp(5),zlam,besi,besd,fthet,fphi
      complex*16 sxx,sxy,sxz,syy,syz,szz
      common /rf2 / sxx(5),sxy(5),sxz(5),syy(5),syz(5),szz(5)
      common /rfzl/ zlam(5),besi(7,5),besd(7,5)
      common /rfz / zeta(11,5),ck(5),fthet(11,5),fphi(11,5)
      common /rf4 / xnparz,xnpar,dbg(99)
      common /rffq/ wc(5),wp2(5),wp(5),vtx(5),vtz(5)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      kzz = CMPLX (0.0, 0.0)
      do ii=1,neq
        szz(ii) = CMPLX (0.0, 0.0)
      end do
      do 1 k=1,neq
      go to 101
  101 temp(k) = CMPLX (1.0/neq, 0.0)
     .        - wp2(k) * EXP (-zlam(k))*besi(nhp,k)*fphi(nh2p,k)
     .         /xkpar/xkpar/vtz(k)/vtz(k)
      kzz = kzz+temp(k)
      szz(k)  = szz(k)+temp(k)
      dbg(13) =  REAL (kzz)
      dbg(14) = PIMAG (kzz)
      do 2 m=1,nhigh
      mn = m+nhigh
      temp(k) = -wp2(k) * EXP (-zlam(k))*besi(m,k)*(zeta(m,k)*fphi(m,k)
     .          +zeta(mn,k)*fphi(mn,k))/(2.0 * pi * freq*xkpar*vtz(k))
      kzz     = kzz    + temp(k)
      szz(k)  = szz(k) + temp(k)
      dbg(15) =  REAL (kzz)
      dbg(16) = PIMAG (kzz)
    2 continue
    1 continue
      return
c
      end

      subroutine lintrp (rt,rfrad,jrfmin,jrfmax,qrfrad,qt)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          LINTRP performs a linear interpolation on qrfrad to
c          estimate the linear power density (qt) at rt.
c          (rfrad(jrfmin) < rt < rfrad(jrfmax),
c          rfrad(i) < rfrad(i+1)  and
c          qrfrad(i) = power density at rfrad(i))
c
      dimension  rfrad(*), qrfrad(*)
c
      qt = 0.0
      if (rt .gt. rfrad(jrfmin))  go to 110
      qt = qrfrad(jrfmin)
      return
  110 if (rt .lt. rfrad(jrfmax))  go to 120
      qt = qrfrad(jrfmax)
      return
  120 do 130 i=jrfmin+1,jrfmax
      if (rt .le. rfrad(i))  go to 140
  130 continue
  140 qt = qrfrad(i-1) + (qrfrad(i)-qrfrad(i-1))*
     .             (rt-rfrad(i-1))/(rfrad(i)-rfrad(i-1))
      return
c
      end

      subroutine modbes (x, l, bi)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     modified Bessel function i from order 0 to order l-1
c     ebi(j) is i of order j-1
c
      dimension  bi(l), bim(5)
      real*8     zero
c
      zero = 0.0
      call beslci (x, zero, l, 1, bi, bim, ncalc)
      do i=1,l
        bi(i) = bi(i) * EXP (-x)
      end do
      if (l .ne. ncalc)  ier = 999
      return
c
      end

      subroutine param0 (i)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c subroutine to calculate wave and plasma parameters to set up matrix equation
c
      include 'ichpar.i'
      include 'ichcom.i'
c
c *** kinc = 1    high field fast wave incidence
c          = 2    low  field fast wave incidence
c
      complex*16 skxx,skxy,skxz,skyy,skyz,skzz
      complex*16 ey,eytot,xkperp,efld,dielt,xkp2
      complex*16 rfa,rfb,rfc,kie(4),xkie,sie(4)
      common /matx/rfa(kich,4),rfb(kich,4),rfc(kich,4)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      common /rf3/xrcp(kich),xkperp(kich,2),xkp2(kich,2)
      common /rfp/efld(kich,3,4),dielt(6,2),exy2(kich,2),
     .           ezy2(kich,2),elr2(kich,2)
      common /ichflx/xkie(kich,4),ey(kich,4),eytot(kich)
      common /diel/skxx(kich,ks,4),skxy(kich,ks,4),skxz(kich,ks,4),
     .            skyy(kich,ks,4),skyz(kich,ks,4),skzz(kich,ks,4)
      dbg67 = REAL (xkperp(i,1))
      dbg48 = REAL (xkperp(i,2))
c
      if (i .eq. 1)  go to 102
      if (i .eq. ngrid)  go to 104
      rdel = 0.5*(rmrf(i+1)-rmrf(i-1))
      go to 103
  102 rdel = 0.5*(rmrf(i+1)-rmrf(i))
      go to 103
  104 rdel = 0.5*(rmrf(i)-rmrf(i-1))
  103 go to (61,62),kinc
   61 continue
      if (dbg67 .gt. 0.0)  go to 100
      go to 63
   62 if (dbg67 .le. 0  )  go to 100
   63 kie(3) =  xkperp(i,1)
      kie(1) = -kie(3)
      go to 101
  100 kie(1) =  xkperp(i,1)
      kie(3) = -kie(1)
  101 go to (65, 66), kinc
   65 if (dbg48 .gt. 0.0)  go to 200
      go to 67
   66 if (dbg48 .le. 0.0)  go to 200
   67 kie(4) =  xkperp(i,2)
      kie(2) = -kie(4)
      go to 201
  200 kie(2) =  xkperp(i,2)
      kie(4) = -kie(2)
  201 go to (71, 72), kinc
   72 do k=1,4
        efld(i,3,k) = -efld(i,3,k)
        do j=1,neq
          skxz(i,j,k) = -skxz(i,j,k)
          skyz(i,j,k) = -skyz(i,j,k)
        end do
      end do
c
   71 do k=1,4
        sie (  k) = CMPLX (0.0, 1.0)*kie(k)*rdel
        xkie(i,k) = kie(k)
        rfa (i,k) = xkpar*efld(i,1,k)-kie(k)*efld(i,3,k)
        rfb (i,k) = kie(k)-ykperp*efld(i,1,k)
        rfc (i,k) = EXP (sie(k))
      end do
      return
c
      end

      subroutine param1 (frq,vi,ve,li,le,nsp,nslb,wkpl,fci,fce,
     .                   zips,zi0,zims,zeps,zems,ze0,li0)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c PARAM1 calculates parameters for z-functions up to li(i) harmonics
c for ions and le harmonics for electrons
c
      include 'ichp2.i'
c
      parameter (nx = 2, nz = 3, nne = 1)
      dimension vi(nx,ny),li(nx),fci(nx,ny),ve(ny),
     .          fce(ny),zips(nx,ny,nz),zi0(nx,ny),zims(nx,ny,nz),
     .          zeps(ny,nne),zems(ny,nne),ze0(ny),work(10)
c
      do 100 j=1,nslb
      do 110 l=1,le
      zeps(j,l) = frq+l*fce(j)
      zems(j,l) = frq-l*fce(j)
      zeps(j,l) = zeps(j,l)/(wkpl*ve(j))
      zems(j,l) = zems(j,l)/(wkpl*ve(j))
  110 continue
      ze0(j) = frq/(wkpl*ve(j))
      do 120 i=1,nsp
      work(i) = wkpl*vi(i,j)
      do 130 l=1,li(i)
      zips(i,j,l) = (frq+l*fci(i,j))/work(i)
      zims(i,j,l) = (frq-l*fci(i,j))/work(i)
  130 continue
      zi0(i,j) = frq/work(i)
  120 continue
  100 continue
      return
c
      end

      subroutine pfunc (x, y, u, v)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     PFUNC calculates  z / (SQRT (pie) * i) = (u,v)  at point (x,y)
c
      dimension w287(4),w283(4),et2(63)
      data      pie/3.14159265358979/, osqpi/0.56418958354776/
c
c     8 point hermite gaussian quadrature  (+,-) x(i)
c
      data w283/0.381186990207322,1.157193712446780,1.981656756695843,
     .          2.930637420257244/
c
c     8 point hermite gaussian quadrature w(i)
c
      data w287/0.6611470125582,0.2078023258149,1.707798300741e-2,
     .          1.996040722114e-4/
      data itest/0/
c
      eighty = 80.0
c
      if (itest .gt. 0)  go to 5
c
c     set up array of e(-t*t)
c     values of EXP (-t*t) for values of 0.08 le x le 5.04 by steps of 0.08
c
      do 3 i=1,63
        t      = FLOAT (i)*0.08
        et2(i) = EXP (-t*t)
    3 continue
c
      itest = 1
c
c     find which quadrant point u = (x,y) is in
c
    5 ii = 1
      assign 244 to j
      c5 = x
      c6 = y
      if (c5 .lt. 0.0)  go to 8
      if (c6 .lt. 0.0)  go to 287
      go to 11
    8 if (c6 .ge. 0.0)  go to 14
      assign 245 to i
      go to 20
   11 assign 257 to i
      go to 46
   14 assign 255 to i
      go to 46
c
c     point has y < 0, calculate analytic continuation part of
c     expression = 2 * EXP (-u*u) / SQRT (pie)
c
   20 z  = c6*c6-c5*c5
      z  = MIN (z,  eighty)
      z  = MAX (z, -eighty)
      co = EXP (z)
      c7 = co + co
      co = c5 * c6
      c9 = co + co
      c8 = -c7 * SIN (c9)
      c7 =  c7 * COS (c9)
   46 c5 = ABS (c5)
      c6 = ABS (c6)
c
c     test for magnitude of y to pick step size to do
c     integral 1.0/(pie*i) * EXP (-t*t)/(t-u)
c
      if (c5 .ge. 6.0)  go to 219
   50 if (c6 .le. 0.5)  go to 65
      if (c6 .le. 3.0)  go to 61
      if (c6 .gt. 6.0)  go to 219
      ns = 6
      go to 73
   61 if (c6 .le. 1.5)  go to 71
      ns = 3
      go to 73
c
c     will do integral for y=0.5 even if y is actually < 0.5
c     will use taylor series to get final answer.
c
   65 c10 = c6
      c6  = 0.5
      assign 128 to j
   71 ns = 1
c
c     evaluate integrand at zero
c
   73 c62 = c6*c6
      c19 = c5*c5
      c20 = c62+c19
      c19 = 1.0/c20
      c17 = c19*c6
      c18 = c19*c5
      c9  = 0.08
c
c     do trapezoidal integration from -5 < x < 5
c
      do 100 m=ns,63,ns
      ws = FLOAT (m)*c9
      c1m = c5-ws
      c1p = c5+ws
      c20p = et2(m)/(c62+c1p*c1p)
      c20m = et2(m)/(c62+c1m*c1m)
      c17 = c6*(c20p+c20m)+c17
      c18 = c1m*c20m+c1p*c20p+c18
  100 continue
      ws = c9 / pie * FLOAT (ns)
      c17 = c17*ws
      c18 = c18*ws
      go to j,(128,244)
c
c     do Taylor series expansion for case y < 0.5 to get correct answer
c
  128 c11 = c17
      c12 = c18
      c9 = 2.0
      c6 = c10-0.5
      c6 = c6+c6
      c10 = c11/2.0
      c13 = (c5*c12+c10-osqpi)*c6
      c10 = c12/2.0
      c14 = (-c5*c11+c10)*c6
      c17 = c11+c13
      c18 = c12+c14
  165 c10 = c6/c9
      c19 = c13/2.0
      c19 = c5*c14+c19
      c15 = (c6/2.0*c11+c19)*c10
      c17 = c15+c17
      t1 = c5*c13
      c19 = (c6*c12+c14)/2.0
      c16 = (-t1+c19)*c10
      c18 = c16+c18
      t1 = c17+c15
      if (ABS (t1-c17) .gt. 1.0e-7 * ABS (t1+c17))  go to 207
      t1 = c18+c16
      if (ABS (t1-c18) .le. 1.0e-7 * ABS (t1+c18))  go to 244
  207 c11 = c13
      c12 = c14
      c13 = c15
      c14 = c16
      c9 = c9+1.0
      go to 165
c
c     do 8 point hermite gaussian quadrature for values of x or y gt 6
c
  219 c17 = 0.0
      c18 = 0.0
      c62 = c6*c6
      do 230 m=1,4
      c1m = c5-w283(m)
      c1p = c5+w283(m)
      c2m = w287(m)/(c62+c1m*c1m)
      c2p = w287(m)/(c62+c1p*c1p)
      c17 = c6*(c2m+c2p)+c17
      c18 = c1m*c2m+c1p*c2p+c18
  230 continue
      c17 = c17/pie
      c18 = c18/pie
c
c     reset value of z for correct quadrant
c
  244 go to i,(245,249,255,257)
  245 c8 = -c8
      c18 = -c18
  249 c17 = c7-c17
      c18 = c8-c18
  255 c18 = -c18
  257 u = c17
      v = c18
      return
  287 c5 = -c5
      assign 249 to i
      go to 20
c
      end




      subroutine plsmap (dni,dne,zi,gmi,ti,te,frq,b0,nsp,nslb,
     .                   fci,fce,vi,ve,rhoi,fpi,fpe)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  param calculates the plasma parameters necessary for evaluation
c  of the dispersion relation
c
c  dni(i,j) = ion density of ith species at jth slab
c  dne(j) = electron density
c  zi(i) = charge state of ith species
c  gmi(i) = mass number of ith species
c  nsp = number of species
c  ti(i,j) = temperature of ith species at jth slab
c  te(j) = electron temperature
c  frq = frequency, b0(j) = toroidal b-field in jth slab
c  fci(i,j) = ion cyclotron frequency of ith species
c  fce(j) = electron cyclotron frequency
c  vi(i,j) = ith species thermal velocity
c  ve(j) = electron thermal velocity
c  rhoi(i,j) = gyroradius of ith species
c  fpi(i,j) = plasma frequency of ith ion species
c  fpe(j) = electron plasma frequency
c  nslb = number of slabs
c
      include 'ichp2.i'
c
      parameter (nx = 2, nz = 3)
      dimension dni(nx,ny),dne(ny),zi(nx),gmi(nx),
     .          ti(nx,ny),te(ny),b0(ny),fci(nx,ny),fce(ny),
     .          vi(nx,ny),ve(ny),rhoi(nx,ny),fpi(nx,ny),fpe(ny)
c
      data fec/1.76e7/, fic/9.58e3/,
     .     fip/1.32e3/, fep/5.64e4/, fte/5.9255e7/, fti/1.3845e6/
c
      do j=1,nslb
        dne(j) = 0.0
        do i=1,nsp
          fpi (i,j) = fip*zi(i) * SQRT (dni(i,j)/gmi(i))
          dne (  j) = dne(j)+zi(i)*dni(i,j)
          fci (i,j) = fic*zi(i)*b0(j)/gmi(i)
          vi  (i,j) = fti * SQRT (ti(i,j)/gmi(i))
          rhoi(i,j) = vi(i,j)/fci(i,j)
        end do
        fce(j) = fec*b0(j)
        ve (j) = fte * SQRT (te(j))
        fpe(j) = fep * SQRT (dne(j))
      end do
c
      return
c
      end

      subroutine powr (x, zr, n, rem, tol)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      s = 1.0
      do 10 n=1,40
      x1 = FLOAT (n)
      p = -x**2*2.0/(2.0*x1+1.0)
      if (n .eq. 1)  go to 11
      p2 = p2*p
      go to 12
   11 p2 = p
   12 rem = ABS (p2) / ABS (s)
      if (rem .lt. tol)  go to 13
   10 s = s+p2
   13 zr = -2.0*x*s
      return
c
      end

      subroutine qcon (icentr,rfrad,rfrow,jrfmin,jrfmax,qrfrad,
     .                 r0,nj,r,dvol,qrf)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c QCON converts a linear power density profile produced
c by an RF code into a radial volume power density profile.
c qrfrad(i) contans the linear power density at rfrad(i)
c (which corresponds to rho = rfrow(i)).  The profile is
c given at rfrad(jrfmin) <= r-major <= rfrad(jrfmax).
c The linear power density is integrated along the line
c segments falling within a flux shell, and the resulting
c power is divided by the volume of that shell.
c Set icentr = 1 for 'onedee' cases with ylaunch = 0.0, and
c icentr = 0 for other cases.
c
      dimension rfrad(*),rfrow(*),qrfrad(*),r(nj),dvol(nj),qrf(nj)
c
      do 110 j=1,nj
  110   qrf(j) = 0.0
      if (icentr .eq. 0)  go to 200
c
c icentr = 1 for 'onedee' case with ylaunch = 0.0
c Begin at rfrad(jrfmin) or r0-r(nj), whichever is greater.
c Integrate linear power density to compute volume power density in each shell.
c
      r2     = r(nj)
      rm2    = rfrad(jrfmin)
      qrf2   = qrfrad(jrfmin)
      ifirst = 0
c
      do 150 j=nj,2,-1
        r1   = r2
        rm1  = rm2
        qrf1 = qrf2
        r2   = 0.5*(r(j)+r(j-1))
        if (rfrad(jrfmin) .ge. r0-r2)  go to 150
        if (rfrad(jrfmin) .ge. r0-r1)  ifirst = 1
        if (ifirst .eq. 0) then
          rm1 = r0 - r1
          call lintrp (rm1, rfrad, jrfmin, jrfmax, qrfrad, qrf1)
        end if
        ifirst = 1
        rm2    = r0 - r2
        call lintrp (rm2, rfrad, jrfmin, jrfmax, qrfrad, qrf2)
        call intqdr (rm1, rm2, rfrad, jrfmin, jrfmax, qrf1, qrf2,
     .               qrfrad, pwr)
        qrf(j) = pwr / dvol(j)
  150 continue
c
c          Integrate within center shell
c
      rm1  = rm2
      qrf1 = qrf2
      rm2  = r0 + r2
      call lintrp(rm2,rfrad,jrfmin,jrfmax,qrfrad,qrf2)
      call intqdr(rm1,rm2,rfrad,jrfmin,jrfmax,
     .            qrf1,qrf2,qrfrad,pwr)
      qrf(1) = pwr / dvol(1)
c
c Integrate outwards to r0+r(nj) or rfrad(jrfmax), whichever is lesser.
c
      do j=2,nj
        r1   = r2
        rm1  = rm2
        qrf1 = qrf2
        if (j .lt. nj)  r2 = 0.5*(r(j)+r(j+1))
        if (j .eq. nj)  r2 = r(nj)
        if (rfrad(jrfmax) .le. r0+r2)  go to 170
        rm2 = r0 + r2
        call lintrp(rm2,rfrad,jrfmin,jrfmax,qrfrad,qrf2)
        call intqdr(rm1,rm2,rfrad,jrfmin,jrfmax,
     .              qrf1,qrf2,qrfrad,pwr)
        qrf(j) = qrf(j) + pwr/dvol(j)
      end do
c
      return
c
  170 rm2 = rfrad(jrfmax)
      qrf2 = qrfrad(jrfmax)
      call intqdr(rm1,rm2,rfrad,jrfmin,jrfmax,
     .            qrf1,qrf2,qrfrad,pwr)
      qrf(j) = qrf(j) + pwr/dvol(j)
      return
c
c icentr = 0 cases
c
  200 return
c
      end

      subroutine qnavg (navg, nj, qrf, qt)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c Qnavg smoothes the power source profile with
c a moving average over navg adjacent points
c
      dimension qrf(nj),qt(nj)
c
      do 100 j=1,nj
      qt(j) = qrf(j)
 100  qrf(j) = 0.0
c
c do points within navg of either end
c
      do 120 j=1,navg
      npav = j+navg
      do 110 jj=1,npav
      qrf(j) = qrf(j) + qt(jj)
 110  qrf(nj-j+1) = qrf(nj-j+1) + qt(nj-jj+1)
      qrf(j) = qrf(j)/npav
 120  qrf(nj-j+1) = qrf(nj-j+1)/npav
c
c do points in between
c
      npav = 2*navg+1
      do 140 j=navg+1,nj-navg
      do 130 jj=j-navg,j+navg
  130 qrf(j) = qrf(j) + qt(jj)
  140 qrf(j) = qrf(j)/npav
      return
c
      end

      subroutine raytrace (freqcy,ptot0,nalfa,ebkev,
     .                     atmf,azf,nrt,nthin,idrive,thgril, 
     .                     height,psi0,maxref, islofa,
     .                     nscr, codeid, rmajor, rminor, btor, curtot,
     .                     kappa, qsaf0, qsafa, r, hcap, nfwsmth, qrfe,
     .                     qrfi, currf, totrtpe, totrtpi, totrtc)
c
c******   TKM   9/99   'thgril' added in argument list
c******   New interface for  Onetwo v3.0, HSJ 05/20/03
c
c
c ----------------------------------------------------------------------
c
c  This is an interface routine between ONETWO and the LH/FW ray
c  tracing routine curray.
c  The broad sequence of operations is:
c  (1) Two files will be written by ONETWO to pass data from ONETWO
c      to the raytracing code.  One file (curray_in) will contain
c      data for setting up the ray configuration.
c       The other file,trxpl.out,will provide profile data .
c  (2) The ray tracing code will be run, producing output file
c      (raytrout) containing power deposition and current drive data.
c  (3) The data from raytrout will be used to give the output
c      variables from the subroutine.
c
c  Output variables are:
c     qrfe (j), j=1,nj    power density to electrons (W/cm**3)
c     qrfi (j), j=1,nj    power to (all) ion species (W/cm**3)
c     currf(j), j=1,nj    RF current density (A/cm**2)
c     totrtpe, totrfpi, totrtc, are respectively total power
c        to electron (W), to ions (W), and total RF current (Amps).
c ----------------------------------------------------------------------
c
c      parameter   (kjrt = 101, kspec = 10)  removed HSJ 5/20/03
       USE rf,only : nrayptrt,powers,nminor,wrt_trxpl_out,incrt,
     .                  nnkparrt, rcurray,pwe,irayiort,bmaxrt,
     .                  pwc,pwit,pwi,nnkpolrt,anpinfrt,anpsuprt,
     .                wrt_curray_in,runid_curray,anzinfrt,anzsuprt,
     .                curray_in_spawn,trxplout_spawn,save_curray_input,
     .                curray_fi

       USE param,only : kj
       USE ext_prog_info , only : get_curray
       USE io, only : ncrt,nout
       USE solcon,only : time 
       USE numbrs
       implicit  integer (i-n), real*8 (a-h, o-z)
c       include 'numbrs.i'
c       include 'io.i'
c       parameter ( kspec =10)

c
      real*8       kappa,pwtotc
      character*8  codeid
      character(len = 256),save :: curray_path_out
      character(len =256)       :: command
      dimension    r(nj),hcap(nj),qrfe(nj),qrfi(nj),currf(nj)
      integer istat
      logical,save ::  first_time 
      data first_time / .true./
c





c --- create file curray_in for passing input to curray,
c --- save file curray_in if requested by user:
      itrans = 1
      nprofs = nj





c       set nspect to number of bins in powers that have non zero power
c       allow fraction curdrive with irfcur:

      runid_curray ='Curray input file created by Onetwo'
      call wrt_curray_in(nprim,nion,nimp,ptot0,ebkev,nalfa,atmf,azf,
     .                   freqcy,nthin,islofa,psi0,maxref,
     .                   irayiort,incrt,bmaxrt,idrive,nprofs,thgril,
     ,                   nminor,height,time)



      interface = 1

      if(interface .eq. 0)then  !this block is obsolete  HSJ 06/19/03
c
c ---    Create raytrin
c ---    next two lines added 21 feb 92 for running CURRAY (raytrace / s.c.chiu)
c
         call getioun(irtrace,42) ! local unit no. for files raytrin, raytrout,trxpl.out
c
         open (unit = irtrace, file = 'raytrin', status = 'UNKNOWN')
c
         jrt = 1
         write (irtrace, 1000)  freqcy, nrayptrt, nrt, codeid
         write (irtrace, 1001) (powers(j),j=1,nrt)
         write (irtrace, 1003) (nnkparrt(j,jrt),j=1,nrt)
         write (irtrace, 1001) (anzinfrt(j,jrt),j=1,nrt)
         write (irtrace, 1001) (anzsuprt(j,jrt),j=1,nrt)
c        following were added t match TK mau version of raytrace
c        these variables are also in curray_in
         write (irtrace, 1003) (nnkpolrt(j,jrt),j=1,nrt)
         write (irtrace, 1001) (anpinfrt(j,jrt),j=1,nrt)
         write (irtrace, 1001) (anpsuprt(j,jrt),j=1,nrt)

c******    TKM   9/99   'thgril' added in write command
         write (irtrace, 1002)  nthin,maxref,islofa,height,thgril,
     .                       curtot,btor

c        write (irtrace, 1002)  nthin,maxref,islofa,height,curtot,btor
         write (irtrace, 1001)  rmajor,rminor,kappa,qsaf0,qsafa
c
c --- next two lines added 21 feb 92 for running CURRAY (raytrace / s.c.chiu)
c
         write (irtrace, 1003)  nprofs, itrans
         write (irtrace, 1001)  (r(j), j=1,nj) 
 1000    format (1pe16.9,2i5,2x,a8)
 1001    format (5(1pe16.9))
 1002    format (3i5,4(1pe16.9))
 1003    format (10i5)
c
         if (codeid .eq. 'onedee') then
c
c          write eqdsk extension into file raytrin, and close file
c
           call eqdskext (irtrace)
           call giveupus(irtrace)
           close (unit =  irtrace)
c
         else
c
c          2D equilibrium case.  Close raytrin, write extended eqdsk
c          (file named eqdskex) and close the file
c
           call giveupus(irtrace)
           close (unit = irtrace)
           call eqdskrt
         end if
      else   !interface .eq. 1, new 05/029/03 interface to curray using trxpl.out:
          print *,'nj,nion ,r(nj) =',nj,nion,r(nj)
          call wrt_trxpl_out(nj,nion,r)

      endif


       



      if(first_time)then
         !get fully qualified name of curray to run:
          call get_curray(ncrt,nout,curray_path_out,len_str)
          first_time = .false.
      endif
      command = ADJUSTL(curray_path_out(1:LEN_TRIM( curray_path_out)))
      command = command(1:LEN_TRIM(command))//' '//
     .      curray_in_spawn(1:LEN_TRIM(curray_in_spawn))
c     execute ray-tracing code (CURRAY) and then connect to file raytrout
c
      write  (6, '(/ '' ---- CURRAY started'')')
      write(6, FMT ='(" running : ",a)')command
c
      if (ISHELL (command) .lt. 0)
     .  call STOP ('subroutine RAYTRACE: failure of spawned CURRAY', 67)
c
      write  (6, '(/ '' ---- CURRAY finished'')')
      if(save_curray_input .eq. 1)then
      command = 'mv '//'trxpl.out '//
     .            trxplout_spawn(1:LEN_TRIM(trxplout_spawn))
          if (ISHELL (command) .lt. 0)   
     .   call STOP ('sub wrt_curray: failure of spawned mv command', 67)
      endif
c
c     Read data produced by ray-tracing code,
c     and put into form for output variables in this subroutine
c
       
      call getioun(irtrace,42)
      open (unit = irtrace, file = 'raytrout', status = 'OLD')
c
      read   (irtrace, '(2i5)')   nprofs, nspec           !grid size,#species

      !allocate storage accordingly:
c  -----------------------------------------------------------------------
      if(.not. allocated(rcurray))then
         allocate (rcurray(nprofs),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("rcurray,sub raytrace",0,istat)
      endif

      if(.not. allocated(pwe))then
         allocate (pwe(nprofs),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("pwe,sub raytrace",0,istat)
      endif

      if(.not. allocated(pwc))then
         allocate (pwc(nprofs),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("pwc,sub raytrace",0,istat)
      endif

      if(.not. allocated(pwit))then
         allocate (pwit(nprofs),STAT = istat)
        if(istat .ne. 0)
     .        call allocate_error("pwit,sub raytrace",0,istat)
      endif

      allocate (pwi(nprofs,nspec),STAT = istat)
      if(istat .ne. 0)
     .        call allocate_error("pwi,sub raytrace",0,istat)

c  ----------------------------------------------------------------------

c     continue reading data:
      read   (irtrace,   3010 )  (rcurray(j), j=1,nprofs) !normalized rho ?
      read   (irtrace,   3010 )  (    pwe(j), j=1,nprofs) !power to electrons


      do n=1,nspec
        read (irtrace,   3010 )  (pwi(j,n),j=1,nprofs)     !power to ion specie
      end do
      read (irtrace, 3010)  (pwc(j),j=1,nprofs)
 3010 format (5(1pe19.6))


c

c     renormalize fraction to ions and electron if called for:
      if(curray_fi .gt.  0.0 .and. curray_fi .lt. 1.0)then
        do j=1,nprofs
           pwtoti =0.0
           do n=1,nspec
              pwtoti    = pwtoti + pwi(j,n)
           enddo
           pwtot   = pwtoti  + pwe(j)
           pwe(j) = (1.-curray_fi)*pwtot
           do n=1,nspec
             pwi(j,n) =  curray_fi*pwtot*pwi(j,n)/pwtoti
           enddo
        enddo
      endif

      do j=1,nprofs
        rcurray(j) = rcurray(j) * r(nj)                 !convert rcurray to cm
          pwit(j) = 0.0
        do n=1,nspec
          pwit(j) = pwit(j) + pwi(j,n)
        end do
      end do


c
c     Cast powers and current onto ONETWO grid, here nprofs may not
c     be eqaul to nj:
c
      
      call intrp (0, 1, rcurray, pwe , nprofs, r, qrfe , nj)
      call intrp (0, 1, rcurray, pwit, nprofs, r, qrfi , nj)
      call intrp (0, 1, rcurray, pwc , nprofs, r, currf, nj)

c
c     Smooth out profiles qrfe, qrfi, currf. Note that
c     no smoothing is the default(nfwsmth = 0) so this
c     is in effect only if user requested it through inone.
c
      call smthsj (r, qrfe , nj, nfwsmth)
      call smthsj (r, qrfi , nj, nfwsmth)
      call smthsj (r, currf, nj, nfwsmth)
c
c
c     Total up powers and currents
c
      call trapv (r, qrfe , hcap, nj, totrtpe)
      call trapv (r, qrfi , hcap, nj, totrtpi)
      call trapv (r, currf, hcap, nj, totrtc)
c
      pi      = ATAN2 (0.0, -1.0)
      totrtpe = 4.0 * pi**2  * rmajor * totrtpe
      totrtpi = 4.0 * pi**2  * rmajor * totrtpi
      totrtc  = 2.0 * pi * totrtc
      print *,'tot pe,pi,jc =',totrtpe,totrtpi,totrtc
c
      call giveupus(irtrace)
      close (unit = irtrace)
      deallocate (pwi) ! nspec may change from call to call


      return
c
      end





      subroutine relkx2 (fpi, nsp, nslb, fpe, fci, fce, kxr2, frq,
     .                   s, d, wkpl)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  relkx2 calculates the square of the real part of kx
c  using the cold plasma dispersion.
c  output is kxr2(l),s(l),d(l) for each slab
c
      include 'ichp2.i'
c
      parameter (nx = 2)
      real*8     kxr2
      dimension  fpi(nx,ny), fpe(ny), fci(nx,ny), fce(ny), kxr2(ny),
     .           s(ny), d(ny), wk(10)
      data       c/2.9979e10/
c
      wk0  = frq / c
      wk02 = wk0**2
      wnz  = wkpl/wk0
      wnz2 = wnz**2
      do l=1,nslb
        s(l) = 1.0
        d(l) = 0.0
        do i=1,nsp
          wk(1) = fpi(i,l)**2
          wk(2) = frq**2-fci(i,l)**2
          wk(3) = wk(1)/wk(2)
          s(l)  = s(l)-wk(3)
          d(l)  = d(l)+wk(3)*frq/fci(i,l)
        end do
        wk(1)   = fpe(l)**2/fce(l)**2
        s(l)    = s(l)+wk(1)
        kxr2(l) = (s(l)-wnz2)**2-d(l)**2
        kxr2(l) = kxr2(l)/(s(l)-wnz2)
        kxr2(l) = kxr2(l)*wk02
      end do
      return
c
      end

      subroutine rffld (i)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
c
      complex*16   efld,           dielt, denom
      common /rfp/ efld(kich,3,4), dielt(6,2), exy2(kich,2),
     .             ezy2(kich,2), elr2(kich,2)
c
      do l=1,2
        denom = dielt(1,l)-dielt(3,l)*dielt(3,l)/dielt(6,l)
        efld(i,3,l) =  (dielt(1,l)*dielt(5,l)+dielt(3,l)*dielt(2,l)) /
     .                  dielt(6,l)/denom
        efld(i,2,l) = CMPLX (1.0, 0.0)
        efld(i,1,l) = -(dielt(2,l)+dielt(3,l)*dielt(5,l)/dielt(6,l)) /
     .                  denom
        efld(i,3,l+2) = -efld(i,3,l)
        efld(i,2,l+2) = efld(i,2,l)
        efld(i,1,l+2) = efld(i,1,l)
        exy2(i,l)     = ABS (efld(i,1,l))**2
        ezy2(i,l)     = ABS (efld(i,3,l))**2
        exey          = exy2(i,l)+1.
        elr2(i,l)     = (exey+2.0 * PIMAG (efld(i,1,l)))/(exey-2.0
     .                            * PIMAG (efld(i,1,l)))
      end do
c
      return
c
      end

      subroutine rfheat (i)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      external det
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      real*8     zero, ten_to_the_minus_7th
      complex*16 det,rr,s,el,aq,bq,cq,detnt,root1,root2
      complex*16 z(2),z1(11,5),dz1(11,5),fthet,fphi,xkperp,xkp2
      dimension eperp(5),infer(2),reden(5),zx(5),zxm(5)
c
      common /rfz/zeta(11,5),ck(5),fthet(11,5),fphi(11,5)
      common /rf1/dene(kich),den1(kich),den2(kich),den3(kich),
     .           terf(kich),tirf(kich),brf(kich)
      common /rf4/xnparz,xnpar,dbg(99)
      common /rf3/xrcp(kich),xkperp(kich,2),xkp2(kich,2)
      common /rffq/wc(5),wp2(5),wp(5),vtx(5),vtz(5)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      common /nn/ nstart,i3,j1
c
      ten_to_the_minus_7th = 1.0e-7
      zero     = 0.0
      w        = wfreq
      rr       = 0.0
      s        = 0.0
      el       = 0.0
      reden(1) = dene(i)
      reden(2) = den1(i)
      reden(3) = den2(i)
      reden(4) = den3(i)
      eperp(1) = terf(i)
      eperp(2) = tirf(i)
      eperp(3) = tirf(i)
      eperp(4) = tirf(i)
c
      do 1 k=1,neq
      wc(k) = charge*anumb(k)*brf(i)/emass(k)/clight
      if (k .eq. 1) wc(1) = -wc(1)
      vtx(k) = SQRT (2.0*eperp(k)*ergtkev/emass(k))
      vtz(k) = SQRT (2.0*eperp(k)*ergtkev/emass(k))
      zx(k)  = 1.0
      zxm(k) = 0.0
      wp2(k) = 4.0 * pi*reden(k)*(charge*anumb(k))**2/emass(k)
      wp(k)  = SQRT (wp2(k))
      do 4 nn=1,nhigh
      mn = nn+nhigh
      zeta(nn,k)  = (w-nn*wc(k))/xkpar/vtz(k)
      zeta(mn,k)  = (w+nn*wc(k))/xkpar/vtz(k)
      dbg(33)     = zeta(nn,k)
      dbg(34)     = zeta(mn,k)
      call fcreal (zeta(nn,k), zr, zi, dzr, dzi, ntrm, rem,
     .             ten_to_the_minus_7th)
      z1(nn,k)    = CMPLX (zr , zero) + CMPLX (0.0, 1.0)*zi
      dz1(nn,k)   = CMPLX (dzr, zero) + CMPLX (0.0, 1.0)*dzi
      call fcreal (zeta(mn,k), zr1, zi1, dzr1, dzi1, ntm, rem,
     .             ten_to_the_minus_7th)
      z1(mn,k)    = CMPLX (zr1 , zero) + CMPLX (0.0, 1.0)*zi1
      dz1(mn,k)   = CMPLX (dzr1, zero) + CMPLX (0.0, 1.0)*dzi1
      fthet(nn,k) = 2.0 * zx(k)*z1(nn,k)/vtz(k)-xkpar*zxm(k)*dz1(nn,k)/w
      fthet(mn,k) = 2.0 * zx(k)*z1(mn,k)/vtz(k)-xkpar*zxm(k)*dz1(mn,k)/w
      fphi(nn,k)  = (1.0 - nn*wc(k)*zxm(k)/w)*dz1(nn,k)
      fphi(mn,k)  = (1.0 + nn*wc(k)*zxm(k)/w)*dz1(mn,k)
      dbg(35)     =  REAL (fthet(nn,k))
      dbg(36)     = PIMAG (fthet(nn,k))
      dbg(37)     =  REAL (fthet(mn,k))
      dbg(38)     = PIMAG (fthet(mn,k))
      dbg(39)     =  REAL (fphi(nn,k))
      dbg(40)     = PIMAG (fphi(nn,k))
      dbg(41)     =  REAL (fphi(mn,k))
      dbg(42)     = PIMAG (fphi(mn,k))
    4 continue
      zeta(nh2p,k)  = w/xkpar/vtz(k)
      dbg(43)       = zeta(nh2p,k)
      call fcreal (zeta(nh2p,k), zr0, zi0, dzr0, dzi0, itm, rem,
     .             ten_to_the_minus_7th)
       z1(nh2p,k)   = CMPLX ( zr0, zero) + CMPLX (0.0, 1.0)* zi0
      dz1(nh2p,k)   = CMPLX (dzr0, zero) + CMPLX (0.0, 1.0)*dzi0
      fthet(nh2p,k) =
     .  2.0*zx(k)*z1(nh2p,k)/vtz(k)-xkpar*zxm(k)*dz1(nh2p,k) / w
      fphi(nh2p,k)  = dz1(nh2p,k)
      dbg(44)       =  REAL (fthet(nh2p,k))
      dbg(45)       = PIMAG (fthet(nh2p,k))
      dbg(46)       =  REAL (fphi(nh2p,k))
      dbg(47)       = PIMAG (fphi(nh2p,k))
      ck(k)         = 0.5*wp2(k)/zx(k)/w/xkpar
c
c     calculation cold plasma dispersion relation for (kperp)**2
c
      if (i .ne. nstart)  go to 1
      rr = rr+0.5*ck(k)*(fthet(1,k)+fthet(nhigh+1,k))
      s = s + CMPLX (0.0, 0.5)*ck(k)*(fthet(1,k)-fthet(nhigh+1,k))
      el = el-2.0*ck(k)*zeta(nh2p,k)*fphi(nh2p,k)/vtz(k)
    1 continue
      if (i .eq.  nstart   )  go to 220
      if (i .ne. (nstart+1))  go to 120
      z(1)  = root1
      z(2)  = root2
      go to 120
  220 isgn  = 0
      ak    = (w/clight)**2
      akp   = xkpar**2/ak
      aq    = rr
      bq    = (akp-rr)*(rr+el)-s*s
      cq    = el*((akp-rr)**2+s*s)
      detnt = bq*bq-4.0 * aq*cq
      root2 = 0.5*(-bq+SQRT (detnt))/aq
      root1 = 0.5*(-bq-SQRT (detnt))/aq
      z(1)  = root1*ak
      z(2)  = root2*ak
      abz1  = ABS (z(1))
      abz2  = ABS (z(2))
      if (abz1 .lt. abz2)  go to 261
      fact1 =  0.001
      fact2 = -1.0
      go to 262
  261 fact1 =  1.0
      fact2 = -0.001
  262 z(1)  = fact1 * z(1)
      z(2)  = fact2*z(2)
  120 xrcp(i) = REAL (z(1))
      call zanlyt1 (det, ten_to_the_minus_7th,
     .              6, 0, 2, 2, z, 200, infer, ier, i)
      if (isgn .ne. 0)  go to 224
      isgn  = 1
      root1 = z(1)
      root2 = z(2)
c
  224 do nn=1,2
        xkp2(i,nn)   =        z(nn)
        xkperp(i,nn) =  SQRT (z(nn))
        dbg(48)      =  REAL (xkperp(i,nn))
        dbg(49)      = PIMAG (xkperp(i,nn))
      end do
c
      dbg(65) =  REAL (z(1))
      dbg(66) = PIMAG (z(1))
      dbg(67) =  REAL (xkperp(i,1))
      dbg(68) = PIMAG (xkperp(i,1))
      call rffld(i)
      return
c
      end

      subroutine rflux (np1, fabs)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      complex*16 ey,efld,dielt,xkie
      complex*16 hz,chz,ehyz,hy,chy,ez,ehzy,eyt,eytot,chzt,ezt,chyt,eh
      real*8     pxx(kich,4),pxt(kich)
      common /rfp/efld(kich,3,4),dielt(6,2),exy2(kich,2),ezy2(kich,2),
     .            elr2(kich,2)
      common /ichflx/xkie(kich,4),ey(kich,4),eytot(kich)
      common /plx/ px1(kich),px2(kich),px3(kich),px4(kich),
     .             pxtot(kich)
      equivalence (pxx,px1),(pxt,pxtot)
c
      do l=1,np1
        eyt  = CMPLX (0.0, 0.0)
        chzt = CMPLX (0.0, 0.0)
        ezt  = CMPLX (0.0, 0.0)
        chyt = CMPLX (0.0, 0.0)
        do i=1,4
          hz = ey(l,i)*(xkie(l,i)-ykperp*efld(l,1,i))
          hy = ey(l,i)*(xkpar*efld(l,1,i)-xkie(l,i)*efld(l,3,i))
          ez = efld(l,3,i)*ey(l,i)
          chz = conjg(hz)
          ehyz = ey(l,i)*chz
          chy = conjg(hy)
          ehzy = ez*chy
          eh = ehyz-ehzy
          pxx(l,i) = 0.5 * REAL (eh)
          eyt = eyt+ey(l,i)
          chzt = chzt+chz
          ezt = ezt+ez
          chyt = chyt+chy
        end do
        pxt(l) = REAL (eyt*chzt-ezt*chyt) * 0.5
      end do
c
      px0 = ABS (pxx(1,1))
c
      do k=1,np1
        pxtot(k) = pxt(k)/px0
        px1(k) = pxx(k,1)/px0
        px2(k) = pxx(k,2)/px0
        px3(k) = pxx(k,3)/px0
        px4(k) = pxx(k,4)/px0
      end do
c
      fabs = (pxtot(1)-pxtot(np1))/pxtot(1)
      return
c
      end

      complex*16 function rfmat (i, j, kflag)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      integer       kflag(2)
      complex*16    rfa,rfb,rfc,efld,dielt
      common /matx/ rfa(kich,4),rfb(kich,4),rfc(kich,4)
      common /rfp / efld(kich,3,4),dielt(6,2),exy2(kich,2),ezy2(kich,2),
     .              elr2(kich,2)
      common /nn  / nstart,i3,j1
c
      if (j .gt. 2 .and. j .lt. j1)  go to 10
      if (j .le. 2)  go to 20
      if (i .lt. i3)  go to 50
      if1 = j-j1+1
      if3 = i-i3+1
      nn  = 4
      if (kflag(2) .gt. 0)  nn = 2
      go to (21, 22, 23, 24), if3
   21 go to (101, 102), if1
  101 rfmat = -CMPLX (1.0, 0.0)
      return
  102 rfmat = -CMPLX (1.0, 0.0)
      return
   22 go to (201, 202), if1
  201 rfmat = -efld(ngrid,3,1)
      return
  202 rfmat = -efld(ngrid,3,nn)
      return
   23 go to (301, 302), if1
  301 rfmat = -rfa(ngrid,1)
      return
  302 rfmat = -rfa(ngrid,nn)
      return
   24 go to (401, 402), if1
  401 rfmat = -rfb(ngrid,1)
      return
  402 rfmat = -rfb(ngrid,nn)
      return
   20 if (i .gt. 4)  go to 50
      nn = 1
      if (kflag(1) .gt. 0)  nn = 2
      go to (1, 2, 3, 4), i
    1 rfmat = rfc(1,j+nn)
      return
    2 rfmat = efld(1,3,j+nn)*rfc(1,j+nn)
      return
    3 rfmat = rfa(1,j+nn)*rfc(1,j+nn)
      return
    4 rfmat = rfb(1,j+nn)*rfc(1,j+nn)
      return
   10 k  = (j-3)/4+2
      ii = 4*(k-2)
      lr = i-ii
      if (lr .le. 0 .or. lr .gt. 8)  go to 50
      lc = j-ii-2
      go to (11, 12, 13, 14, 15, 16, 17, 18), lr
   11 rfmat = -CMPLX (1.0, 0.0)
      return
   12 rfmat = -efld(k,3,lc)
      return
   13 rfmat = -rfa(k,lc)
      return
   14 rfmat = -rfb(k,lc)
      return
   15 rfmat = rfc(k,lc)
      return
   16 rfmat = rfc(k,lc)*efld(k,3,lc)
      return
   17 rfmat = rfc(k,lc)*rfa(k,lc)
      return
   18 rfmat = rfc(k,lc)*rfb(k,lc)
      return
   50 rfmat = 0.0
      return
c
      end

      subroutine rf_mhddat(task)
c
c------------------------------------------------------------------------
c-----This subroutine generates the mhd information required by the
c-----rf codes,using either the input eqdsk,or information form the
c-----time dependent eqdsk mode,or information from the fixed boundary
c-----calculations,or 1-D runs.
c--------------------------------------------------------HSJ-2-14-00-----
c
      USE param
      USE contour
      USE io 
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE rf, only:  gafsep
      USE extra
      USE numbrs
      USE mesh


      USE machin
      USE geom ! rcapi
      USE constnts
      USE soln2d
      USE rhog
      USE mhdcom
      USE bicube
      USE flxav
      USE gpsi
      USE  replace_imsl,              ONLY : my_ibcccu
      implicit  integer (i-n), real*8 (a-h, o-z)

c      include 'small.i'    !p(nw,nh)(= psi(R,Z))
      include 'storage.i'  !work,wdum,zdum,xdum,vdum
c
      real *8 ,dimension(:),allocatable :: psi_rf   !the psi grid for rf calcs

      real *8 zcontour_rf(kstore),rcontour_rf(kstore),
     .        bpol_rf(kstore),rmajor_rf(kstore) ! all these LOCAL storage only
      real *8 rmaxis,zmaxis,psimax_rf,psilim_rf,
     .        sgnpsi,rgmax,rgmin,zgmax,zgmin,psi2d_rf(nw,nh),
     .        rpout,rpin,zptop,zpbot,vol_psi_rf(kj),volume_rf,
     .        xbouni(kj),q_rf(kj),pprim_rf(kj),
     .        fpsi_rf(kj),bmini_rf(kj),bmaxi_rf(kj),dx2i(kj),
     .        cmtom,rgrid_rf(nw),zgrid_rf(nh),psir_rf(kj),
     .        kg_to_tesla
      REAL *8 b_avg_rf(nw),bsq_avg_rf(nw),h_factr_rf(nw),r0rinv_rf(nw)
      equivalence (zcontour_rf(1),zdum(1))
      equivalence (rcontour_rf(1),wdum(1))
      equivalence (bpol_rf(1),xdum(1))
      equivalence (rmajor_rf(1),vdum(1))
      equivalence (psi2d_rf(1,1) ,zdum(1))

      real*8 rpbdry_rf(nrplas), zpbdry_rf(nrplas)
      real*8 sf_rf(kpsi), qpsi_rf(kpsi)
      character (len = *) task



      ngrid_rf = nj      !can modify the grid passed to rf codes here
                         !but logic below not yet consistent
      allocate (psi_rf(1:ngrid_rf),STAT = istat)
      if(istat .ne. 0)
     .          call allocate_error("psi_rf",0,istat) 
      ngrid_rfm1 = ngrid_rf-1
c
      if (codeid .eq. 'onedee') then
c
c ----------------------------------------------------------------------
c calculate psi_rf, an array of flux zone boundaries used in rf codes.
c for a 1-D run the zones are of equal width in minor radius.
c ----------------------------------------------------------------------
c
          do 10 j=1,nj
   10     psir(j) = r(j)
          dpsi = r(nj)/ngrid_rfm1
          do 20 i=1,ngrid_rf
   20     psi_rf(i) = (i-1)*dpsi
c
c ----------------------------------------------------------------------
c calculate volumes of flux zones for elliptical cross sections
c ----------------------------------------------------------------------
c


          factor = 2.0 * pi**2 * rmajor*kappa*(rminor/ngrid_rfm1)**2
          do 30 i=1,ngrid_rfm1
   30     vol_psi_rf(i) = factor*(i**2-(i-1)**2)
c
c         oned coding not done
c
          call STOP ('subroutine RF_MHDDAT: 1D not implemented', 287)
c
      else        ! calculations for 1-1/2 d cases
c
c ----------------------------------------------------------------------
c calculate psi_rf, an array of flux zone boundaries used in rf codes.
c ----------------------------------------------------------------------
c
c          dpsi = (psir(nj)-psir(1))/ngrid_rfm1
c          psi_rf(1)=psir(1)
c          do 100 i=2,ngrid_rfm1
c  100     psi_rf(i) = psi_rf(i-1) + dpsi
c          psi_rf(ngrid_rf)=psir(nj)               !avoid roundoff
          call copya(psir,psi_rf,ngrid_rf)         !assumes same size
          call copya(psir,psir_rf,kj)
          call multpl1(psir_rf,kj,psimks)
c          print *,'psi-rf in rf_mhddat ='
c          print *,psi_rf
c
c         get fpsi on rf grid from Fcap (= f(psilim)/f(psi))
          call copya(fcap,xdum,nj)        !xdum is temp storage
          cconst = rmajor*btor*1.e-3      !btor is in gauss
          do j=1,nj
             fpsi_rf(j)= cconst/xdum(j) !kg-cm
          enddo
 

c
c ----------------------------------------------------------------------
c flux values are assumed to be in  Kgauss-cm**2
c (This will be the case if this subroutine is called from sub source)
c get bicubic representation of psi ( = p):
c (reverse the sign of psi so that CNTOUR will work)
c ----------------------------------------------------------------------
c
          cconst = -1.0
          call multpl1 (p,nwh,cconst)
          call my_ibcccu  (p,rmhdgrid,nw,zmhdgrid,nh,cspln,
     .                     nw,wnoperm,ier)
          call multpl1 (p,nwh,cconst)
c
c ----------------------------------------------------------------------
c generate the psi contours corresponding to the rho grid.
c ----------------------------------------------------------------------
c
          bperr     = 0.05
          iconvg    = 0
          arcl      = 2.0    ! 2 cm arclength increment
          taxis     = 5.0
          tlim      = 30.0
          a         = (tlim-taxis)/((psi_rf(1)-psi_rf(ngrid_rf)))
          bincp     = taxis
          delta_psi = (-psi_rf(ngrid_rf-1)+psi_rf(ngrid_rf))
          do j=1,ngrid_rfm1    ! loop from plasma edge to axis
             i = ngrid_rf-j+1
             psi_psi_rf = psi_rf(i)
             if (ifixshap .eq. 1) then
c
c               the actual contour to be traced has the original (at
c               time t=time0) values of psi associated with the (R,Z)
c               grid IF the equilibrium is not evolved. Under this
c               condition the values of psi in psi_rf will not map
c               back to the (R,Z) contours which are fixed in time.
c               Here we assume that the poloidal flux fraction inside each
c               (R,Z) contour remains invariant as a fucntion of time.
c               get normalized poloidal flux fraction at current time,psi_rfn:
c
                psi_rfn = (psi_psi_rf -psi_rf(1))/
     .                                  ( psi_rf(ngrid_rf)-psi_rf(1))
c
c               get psi_psi_rf, the psi value that represents the same
c               normalized poloidal flux fraction as psi_rfn but at the time the
c               eqdsk was read. Then get the (R,Z) contour points that
c               correspondt to psi_psi_rf on the eqdsk  and ASSUME that
c               this contour corresponds to
c               the current value of psi, psi_rf(i)
                psi_psi_rf= psi_rfn*(psir_at_t0(nj)-psir_at_t0(1)) +
     .                                          psir_at_t0(1)
             endif
c
             ptrace = -psi_psi_rf

c
c --- prevent searching for plasma boundary by moving psilim in slightly:
c
             if (j .eq. 1)  ptrace = ptrace +0.01*delta_psi
c
             iauto  = 0
             iconvg = 0
             dang   = a*(ptrace+psi_rf(1))+bincp
             drx    = 0.0
             dry    = 0.0
             if (j .eq. 1 .and. mhdmethd .ne. 'tdem')  go to 501
             if (j .eq. 1 .and. mhdmethd .eq. 'tdem') then
               iauto  = 1
               iconvg = 0
             end if
c
 
             call cntour (xmagn1,ymagn1,ptrace,rcmin,rcmax,zcmin,
     .                    zcmax,zrcmin,zrcmax,rzcmin,rzcmax,dang,arcl,
     .                    bperr,drx,dry,100.0*xlimiter(nlimiter+1),
     .                    100.0*xlimiter(nlimiter+2),
     .                    100.0*ylimiter(nlimiter+1),
     .                    100.0*ylimiter(nlimiter+2),
     .                    iauto,iautoc,rcontour_rf,zcontour_rf,
     .                    ncontour_rf,rmhdgrid,nw,
     .                    zmhdgrid,nh,cspln,n2cspln,nh2,ncrt,kstore,
     .                    ierr,bpol_rf,iconvg,delta_psi)
            if (ierr .ne. 0)
     .          call STOP ('subroutine RF_MHDDAT: CNTOUR error', 289)
            go to 502
c
  501       call fixedcntour (rplasbdry,zplasbdry,nplasbdry,
     .                        rcontour_rf,zcontour_rf,ncontour_rf,
     .                        rcmin,rcmax,zcmin,zcmax,
     .                        rzcmin,rzcmax,zrcmin,zrcmax,
     .                        rmhdgrid,zmhdgrid,nw,nh,
     .                        bpol_rf,cspln,n2cspln,nh2,pds)
        call limiter_check(rcmin/100.,rcmax/100.,zcmin/100.,zcmax/100.,
     .                    xlimiter,ylimiter,nlimiter)
c
c           get the flux zone volumes:
c
 502        call volcalc (rcontour_rf,zcontour_rf,ncontour_rf,
     .                                 xmagn1,ymagn1,udum(j),area)
            if (j .eq. 1) then
                npbdry_rf=ncontour_rf
                do i=1,npbdry_rf
                   rpbdry_rf(i)=rcontour_rf(i)*0.01
                   zpbdry_rf(i)=zcontour_rf(i)*0.01
                enddo
                rpout = rcmax*.01  ! max outside radius of plasma,m
                rpin  = rcmin *.01 ! max inside radius of plasma,m
                zptop = zcmax*.01  ! max height of plasma, m
                zpbot = zcmin*.01  ! min height of plasma, m
            endif
c     get min/max mod b on contour:
            btor_sq = (fpsi_rf(j)/rcontour_rf(1))**2
            btotal = bpol_rf(1)**2 + btor_sq
            bmod_max = btotal
            bmod_min = btotal
            do k=2,ncontour_rf  !loop over points on contour
               btor_sq = (fpsi_rf(j)/rcontour_rf(k))**2
c              bpol is bpoloidal at each contour point
               btotal = bpol_rf(k)**2 + btor_sq
               bmod_max = MAX (bmod_max, btotal)
               bmod_min = MIN (bmod_min, btotal)
            end do                 ! end loop over contour points
            bmaxi_rf(ngrid_rf-j+1) = SQRT (bmod_max)     ! kgauss
            bmini_rf(ngrid_rf-j+1) = SQRT (bmod_min)     ! kgauss
          end do                  ! end loop over plasma contours
          bmaxi_rf(1)= ABS(fpsi_rf(1))/xmagn1
          bmini_rf(1)= bmaxi_rf(1)
          vol_psi_rf(ngrid_rf) = 0.0
          volume_rf    = 0.0
          udum(ngrid_rf) = 0.0    !volume at magnetic axis
c         vol_psi_rf(ngrid_rfm1)= plasma volume of shell at psilim
c         vol_psi_rf(1)= plasma volume of tube centered at magnetic axis
          dx2i(1)=0.0
          do j=1,ngrid_rfm1 !
            vol_psi_rf(ngrid_rfm1-j+1) = udum(j)-udum(j+1)
            dx2i(ngrid_rf-j+1)=vol_psi_rf(ngrid_rfm1-j+1)*1.e-6 !m**3
            volume_rf = volume_rf+vol_psi_rf(ngrid_rfm1-j+1)
          end do
          volume_rf = volume_rf*1.e-6 !m**3
      endif
c
c     the q profile on the psi_rf grid is obtained by interpolation
c     from the q profile on the psir grid,qpsir (which corresponds
c     to the rho grid as well). Note that q will be negative if btor is negative
c      do j=2,ngrid_rfm1
c         call interp(psi_rf(j),psir,nj,qpsir,dumy)
c         q_rf(j)=dumy
c      enddo
c      q_rf(1)=SNGL(qpsir(1))
c      q_rf(ngrid_rf)= qpsir(nj)
      do j=1,kj
         q_rf(j) = ABS(q(j))
      enddo
c
c     get p-prime on rf grid from pprim:
c      do j=1,ngrid_rf
c         call interp(psi_rf(j),psir,nj,pprim,dumy)
c         pprim_rf(j)= dumy
c      enddo
      call copya(pprim,pprim_rf,kj)
c
c     convert to meters
c
      rmaxis = xmagn1*.01
      zmaxis = ymagn1*.01
      psimax_rf = psir(1)*psimks
      psilim_rf = psir(nj) *psimks
c     get equispaced (in psi) normalized psi:
      do j=1,ngrid_rf
         xbouni(j) = (psi_rf(j)-psi_rf(1))/
     .                     (psi_rf(ngrid_rf)-psi_rf(1))
      enddo
c
      nr_rf = nw    ! nr_rf is nr in subroutine outputb (nr is used in onetwo)
      nz_rf = nh    ! nz_rf     nz
      sgnpsi =-1.
      b0 = btor *1.e-4  !tesla
      cmtom = 0.01
      rgmax = rmhdgrid(nw)*cmtom
      rgmin = rmhdgrid(1)*cmtom
      zgmax = zmhdgrid(nh)*cmtom
      zgmin = zmhdgrid(1)*cmtom
      call copya(p,psi2d_rf,nwh)
      cconst = psimks
      call multpl1 (psi2d_rf,nwh,cconst)
      call copya(rmhdgrid,rgrid_rf,nw)
      call copya(zmhdgrid,zgrid_rf,nh)
      call multpl1(rgrid_rf,nw,cmtom)
      call multpl1(zgrid_rf,nh,cmtom)
      kg_to_tesla = 0.1
      call multpl1(bmini_rf,kj,kg_to_tesla)
      call multpl1(bmaxi_rf,kj,kg_to_tesla)
      cconst=1.e-3
      call multpl1(fpsi_rf,kj,cconst)


      if(task .eq. 'write psiin')then
c           write the binary file for the rf codes:
c           call STOP('wrt_rf_mhd  no longer works called from here',0)
c           call wrt_rf_mhd(n_rf, nw, nh, nj,
c        .     rgmax, rgmin, zgmax, zgmin, rpout, rpin, zptop, zpbot,
c        .     b0, rmaxis, zmaxis, psimax_rf, psilim_rf,
c        .     rgrid_rf, nr_rf, zgrid_rf, nz_rf, psi2d_rf,
c        .     sgnpsi, psir_rf, fpsi_rf, xbouni, dx2i, q_rf,
c        .     bmini_rf, bmaxi_rf)
c
c          write psiin
           lcentr=2
           rdim_rf=rgmax-rgmin
           zdim_rf=zgmax-zgmin
           rzero_rf=rmajor*cmtom
           gasep_rf=gafsep
c          convert some profiles form nj grid to nr_rf grid
c          (which is what psiiin demands:
           sf_rf(1)=fpsi_rf(1)
           sf_rf(nr_rf)=fpsi_rf(nj)
           qpsi_rf(1)=q_rf(1)
           qpsi_rf(nr_rf)=q_rf(nj)
           h_factr_rf(1) = h_factr_cap(1)
           h_factr_rf(nr_rf) = h_factr_cap(nj)
           b_avg_rf(1) = b_avg_cap(1)
           b_avg_rf(nr_rf) = b_avg_cap(nj)
           bsq_avg_rf(1) = bsq_avg_cap(1)
           bsq_avg_rf(nr_rf) = bsq_avg_cap(nj)
           r0rinv_rf(1) = rcapi(1)        ! rcapi = <1/R> , 1/cm
           r0rinv_rf(nr_rf) = rcapi(nj)
           dpsi=(psir(nj)-psir(1))/(nr_rf-1)
           do i=2,nr_rf-1
              psival=psir(1)+(i-1)*dpsi
              call interp(psival, psir, nj, fpsi_rf, sf_rf(i))
              call interp(psival, psir, nj, q_rf,    qpsi_rf(i))
              call interp(psival, psir, nj, b_avg_cap,b_avg_rf(i))
              call interp(psival, psir, nj, bsq_avg_cap,bsq_avg_rf(i))
              call interp(psival, psir, nj, h_factr_cap,h_factr_rf(i))
              call interp(psival, psir, nj, rcapi,r0rinv_rf(i))
           enddo
           r0rinv_rf(1:nr_rf) = rma*r0rinv_rf(1:nr_rf) !<Rmag_axis/R>
   
 
c
           call wrt_psiin(n_rf, lcentr, nj, nr_rf, nz_rf,
     $          rdim_rf, zdim_rf, rzero_rf, rgmin,
     $          rmaxis, zmaxis, psimax_rf, psilim_rf, b0,
     $          gasep_rf, sf_rf, psi2d_rf, nw,
     $          psir, qpsi_rf, npbdry_rf, nlimiter,
     $          rpbdry_rf, zpbdry_rf, xlimiter, ylimiter,
     $          h_factr_rf,b_avg_rf,bsq_avg_rf,r0rinv_rf)
           print *,'file psiin created in sub rf_mhddat'
c

c
      else
         call STOP('task  error in sub rf_mhddat',1)
      endif !task switch 

      deallocate (psi_rf, STAT = istat)
      if(istat .ne. 0)
     .          call deallocate_error("psi_rf",0,istat) 
c 

      return
      end

      subroutine rf_mhd_interface
c  a modified version of this routine, rf_mhddat is being
c  used. rf_mhddat relies on the GAFIT .
c------------------------------------------------------------------------
c-----This subroutine generates the mhd information required by toray
c-----using either the input eqdsk,or information form the
c-----time dependent eqdsk mode,or information from the fixed boundary
c-----calculations,or 1-D runs. File mhddat is written directly
c-----so that file psiin and prgoram Gafit are not needed.
c--------------------------------------------------------HSJ-1/07/03-----
c
      USE param
      USE contour
      USE io 
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE extra
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE constnts
      USE soln2d
      USE rhog
      USE mhdcom
      USE bicube
      USE flxav
      USE gpsi
      USE replace_imsl,           ONLY :  my_ibcccu,my_dbcevl1
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'constnts.i'
c      include 'io.i'       !ncrt,n_rf
c      include 'mhdpar.i'   !nw,nh,nwh
c      include 'bicube.i'   !cspln,wnoperm,pds,nh2,n2cspln
c      include 'contour.i'  !rplasbdry,zplasbdry,nplasbdry
c      include 'extra.i'    !q profile
c      include 'flxav.i'    !xmagn1,ymagn1
c      include 'geom.i'     !codeid,fcap(nj)
c      include 'limiter.i'  ! xlimiter,ylimiter,nlimiter
c      include 'numbrs.i'   !nj
c      include 'machin.i'   ! rmajor,rminor,kappa,btor
c      include 'mhdgrid.i'  ! rmhdgrid,zmhdgrid
c      include 'mhdcom.i'   !mhdmethd
c      include 'mesh.i'     !r(1..nj)
c      include 'rhog.i'     !psir(1..nj),psir_at_t0,qpsir
c      include 'soln2d.i'   !ifixshap
c      include 'small.i'    !p(nw,nh)(= psi(R,Z))
      include 'storage.i'  !work,wdum,zdum,xdum,vdum
c
      real *8 ,dimension(:),allocatable :: psi_rf   !the psi grid for rf calcs
      real *8 zcontour_rf(kstore),rcontour_rf(kstore),
     .        bpol_rf(kstore)     !all these are LOCAL storage only
      real *8 rmaxis,zmaxis,psimax_rf,psilim_rf,
     .        sgnpsi,rgmax,rgmin,zgmax,zgmin,psi2d_rf(nw,nh),
     .        rpout,rpin,zptop,zpbot,vol_psi_rf(kj),volume_rf,
     .        xbouni(kj),q_rf(kj),pprim_rf(kj),
     .        fpsi_rf(kj),bmini_rf(kj),bmaxi_rf(kj),dx2i(kj),
     .        cmtom,rgrid_rf(nw),zgrid_rf(nh),
     .        kg_to_tesla
      equivalence (zcontour_rf(1),zdum(1))
      equivalence (rcontour_rf(1),wdum(1))
      equivalence (bpol_rf(1),xdum(1))
c
      ngrid_rf = nj      !can modify the grid passed to rf codes here
      allocate (psi_rf(1:ngrid_rf),STAT = istat)
      if(istat .ne. 0)
     .          call allocate_error("psi_rf",0,istat) 
c                        
      ngrid_rfm1 = ngrid_rf-1
      
c
      if (codeid .eq. 'onedee') then
c
c ----------------------------------------------------------------------
c calculate psi_rf, an array of flux zone boundaries used in rf codes.
c for a 1-D run the zones are of equal width in minor radius.
c ----------------------------------------------------------------------
c
          do 10 j=1,nj
   10     psir(j) = r(j)
          dpsi = r(nj)/ngrid_rfm1
          do 20 i=1,ngrid_rf
   20     psi_rf(i) = (i-1)*dpsi
c
c ----------------------------------------------------------------------
c calculate volumes of flux zones for elliptical cross sections
c ----------------------------------------------------------------------
c
          factor = 2.0 * pi**2 * rmajor*kappa*(rminor/ngrid_rfm1)**2
          do 30 i=1,ngrid_rfm1
   30     vol_psi_rf(i) = factor*(i**2-(i-1)**2)
c
c         oned coding not done
c
      call STOP ('subroutine RF_MHD_INTERFACE: 1D not implemented', 286)
c
      else        ! calculations for 1-1/2 d cases
c
c ----------------------------------------------------------------------
c calculate psi_rf, a uniform  array of flux zone boundaries
c ----------------------------------------------------------------------
c
c          dpsi = (psir(nj)-psir(1))/ngrid_rfm1
c          psi_rf(1)=psir(1)
c          do 100 i=2,ngrid_rfm1
c  100     psi_rf(i) = psi_rf(i-1) + dpsi
c          psi_rf(ngrid_rf)=psir(nj)             !avoid roundoff
           call copya(psir,psi_rf,ngrid_rf)

c
c         get fpsi on rf grid from Fcap (= f(psilim)/f(psi))
c          call copya(fcap,xdum,kj)            !xdum is temp storage
          cconst = rmajor*btor*1.e-3           !btor is in gauss here
          do j=1,nj
             fpsi_rf(j)= cconst/fcap(j)        !f(psi) in kgauss-cm
          enddo
c
c ----------------------------------------------------------------------
c flux values are assumed to be in  Kgauss-cm**2
c (This will be the case if this subroutine is called from sub source)
c get bicubic representation of psi ( = p):
c (reverse the sign of psi so that CNTOUR will work)
c ----------------------------------------------------------------------
c
          cconst = -1.0
          call multpl1 (p,nwh,cconst)
c          print *,'psi(middle) =',p(nw/2,nh/2) !postive with psi_rf negative
c          print *,'psi_rf(1)(ngrid_rf) =',psi_rf(1),psi_rf(ngrid_rf)
c          print *,'rmhdgrid(nw) =',rmhdgrid(nw)
          call my_ibcccu  (p,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                     wnoperm,ier)
          call multpl1 (p,nwh,cconst)
c
c ----------------------------------------------------------------------
c generate the psi contours corresponding to the rho grid.
c ----------------------------------------------------------------------
c
          bperr     = 0.05 
          iconvg    = 0
          arcl      = 2.0    ! 2 cm arclength increment
          taxis     = 5.0
          tlim      = 30.0
          a         = (tlim-taxis)/((psi_rf(1)-psi_rf(ngrid_rf)))
          bincp     = taxis
          delta_psi = (-psi_rf(ngrid_rf-1)+psi_rf(ngrid_rf))
          do j=1,ngrid_rfm1    ! loop from plasma edge to axis
             i = ngrid_rf-j+1
             psi_psi_rf = psi_rf(i)
             if (ifixshap .eq. 1) then
c
c               the actual contour to be traced has the original (at
c               time t=time0) values of psi associated with the (R,Z)
c               grid IF the equilibrium is not evolved. Under this
c               condition the values of psi in psi_rf will not map
c               back to the (R,Z) contours which are fixed in time.
c               Here we assume that the poloidal flux fraction inside each
c               (R,Z) contour remains invariant as a fucntion of time.
c               get normalized poloidal flux fraction at current time,psi_rfn:
c
                psi_rfn = (psi_psi_rf -psi_rf(1))/
     .                                  ( psi_rf(ngrid_rf)-psi_rf(1))
c
c               get psi_psi_rf, the psi value that represents the same
c               normalized poloidal flux fraction as psi_rfn but at the time the
c               eqdsk was read. Then get the (R,Z) contour points that
c               correspondt to psi_psi_rf on the eqdsk  and ASSUME that
c               this contour corresponds to
c               the current value of psi, psi_rf(i)
                psi_psi_rf= psi_rfn*(psir_at_t0(nj)-psir_at_t0(1)) +
     .                                          psir_at_t0(1)
             endif
c
             ptrace = -psi_psi_rf
c
c --- prevent searching for plasma boundary by moving psilim in slightly:
c
             if (j .eq. 1)  ptrace = ptrace +0.01*delta_psi
c
             iauto  = 0
             iconvg = 0
             dang   = a*(ptrace+psi_rf(1))+bincp
             drx    = 0.0
             dry    = 0.0
             if (j .eq. 1 .and. mhdmethd .ne. 'tdem')  go to 501
             if (j .eq. 1 .and. mhdmethd .eq. 'tdem') then
               iauto  = 1
               iconvg = 0
             end if
c
             call cntour (xmagn1,ymagn1,ptrace,rcmin,rcmax,zcmin,
     .                    zcmax,zrcmin,zrcmax,rzcmin,rzcmax,dang,arcl,
     .                    bperr,drx,dry,100.0*xlimiter(nlimiter+1),
     .                    100.0*xlimiter(nlimiter+2),
     .                    100.0*ylimiter(nlimiter+1),
     .                    100.0*ylimiter(nlimiter+2),
     .                    iauto,iautoc,rcontour_rf,zcontour_rf,
     .                    ncontour_rf,rmhdgrid,nw,
     .                    zmhdgrid,nh,cspln,n2cspln,nh2,ncrt,kstore,
     .                    ierr,bpol_rf,iconvg,delta_psi)
c
            if (ierr .ne. 0)
     .      call STOP ('subroutine RF_MHD_INTERFACE: CNTOUR error', 288)
            go to 502
c
  501       call fixedcntour (rplasbdry,zplasbdry,nplasbdry,
     .                        rcontour_rf,zcontour_rf,ncontour_rf,
     .                        rcmin,rcmax,zcmin,zcmax,
     .                        rzcmin,rzcmax,zrcmin,zrcmax,
     .                        rmhdgrid,zmhdgrid,nw,nh,
     .                        bpol_rf,cspln,n2cspln,nh2,pds)
        call limiter_check(rcmin/100.,rcmax/100.,zcmin/100.,zcmax/100.,
     .                    xlimiter,ylimiter,nlimiter)
c
c             get the flux zone volumes:
c
 502        call volcalc (rcontour_rf,zcontour_rf,ncontour_rf,
     .                                 xmagn1,ymagn1,udum(j),area)
c
            if (j .eq. 1) then
                rpout = rcmax*.01  ! max outside radius of plasma,m
                rpin  = rcmin *.01 ! max inside radius of plasma,m
                zptop = zcmax*.01  ! max height of plasma, m
                zpbot = zcmin*.01  ! min height of plasma, m
            endif
c     get min/max mod b on contour:
            btor_sq = (fpsi_rf(j)/rcontour_rf(1))**2
            btotal = bpol_rf(1)**2 + btor_sq
            bmod_max = btotal
            bmod_min = btotal
            do k=2,ncontour_rf  !loop over points on contour
               btor_sq = (fpsi_rf(j)/rcontour_rf(k))**2
c              bpol is bpoloidal at each contour point
               btotal = bpol_rf(k)**2 + btor_sq
               bmod_max   = MAX(bmod_max,btotal)
               bmod_min   = MIN(bmod_min,btotal)
            end do ! end loop over contour points
            bmaxi_rf(ngrid_rf-j+1) = SQRT (bmod_max)     ! kgauss
            bmini_rf(ngrid_rf-j+1) = SQRT (bmod_min)     ! kgauss
          end do   ! end loop over plasma contours
c
          bmaxi_rf(1)= ABS(fpsi_rf(1))/xmagn1
          bmini_rf(1)= bmaxi_rf(1)
          vol_psi_rf(ngrid_rf) = 0.0
          volume_rf            = 0.0
          udum(ngrid_rf)       = 0.0       ! volume at magnetic axis
c
c         vol_psi_rf(ngrid_rfm1)= plasma volume of shell at psilim
c         vol_psi_rf(1)= plasma volume of tube centered at magnetic axis
c
          dx2i(1)=0.0
          do j=1,ngrid_rfm1
            vol_psi_rf(ngrid_rfm1-j+1) = udum(j)-udum(j+1)
            dx2i(ngrid_rf-j+1)=vol_psi_rf(ngrid_rfm1-j+1)*1.e-6 ! m**3
            volume_rf = volume_rf+vol_psi_rf(ngrid_rfm1-j+1)
          end do
          volume_rf = volume_rf*1.e-6 !m**3
      endif
c
c     the q profile on the psi_rf grid is obtained by interpolation
c     from the q profile on the psir grid,qpsir (which corresponds
c     to the rho grid as well). Note that q will be negative if btor is negative
c      do j=2,ngrid_rfm1
c         call interp(psi_rf(j),psir,nj,qpsir,dumy)
c         q_rf(j)=dumy
c      enddo
c      q_rf(1)= SNGL(qpsir(1))
c      q_rf(ngrid_rf)= qpsir(nj)
      do j=1,kj
         q_rf(j) = ABS(q(j))
      enddo
c
c     get p-prime on rf grid from pprim:
c      do j=1,ngrid_rf
c         call interp(psi_rf(j),psir,nj,pprim,dumy)
c         pprim_rf(j)= dumy
c      enddo
      call copya(pprim,pprim_rf,kj)
c
c     convert to meters
c
      rmaxis = xmagn1*.01
      zmaxis = ymagn1*.01
      psimax_rf = psir(1)*psimks
      psilim_rf = psir(nj) *psimks
      btor_rf   = 1.e-3*fpsi_rf(1)/rmaxis
c     get equispaced (in psi) normalized psi:
      do j=1,ngrid_rf
         xbouni(j) = (psi_rf(j)-psi_rf(1))/
     .                     (psi_rf(ngrid_rf)-psi_rf(1))
      enddo
c
      nr_rf = nw    ! nr_rf is nr in subroutine outputb (nr is used in onetwo)
      nz_rf = nh    !nz_rf     nz
      sgnpsi =-1.
      b0 = btor *1.e-4  !tesla
      cmtom = 0.01
      rgmax = rmhdgrid(nw)*cmtom
      rgmin = rmhdgrid(1)*cmtom
      zgmax = zmhdgrid(nh)*cmtom
      zgmin = zmhdgrid(1)*cmtom
      cconst = psimks
      call copya(rmhdgrid,rgrid_rf,nw)
      call copya(zmhdgrid,zgrid_rf,nh)
      call multpl1(rgrid_rf,nw,cmtom)
      call multpl1(zgrid_rf,nh,cmtom)
      kg_to_tesla = 0.1
      call multpl1(bmini_rf,kj,kg_to_tesla)
      call multpl1(bmaxi_rf,kj,kg_to_tesla)
      cconst=1.e-3
      call multpl1(fpsi_rf,kj,cconst)
      call multpl1(psi_rf,kj,psimks)       !convert psi_rf from
                                                 !kg/cm**2 to  volt-sec/rad
c     write the binary file for the rf codes:
           call wrt_rf_mhd(n_rf,nj,
     .     rgmax, rgmin, zgmax, zgmin, rpout, rpin, zptop, zpbot,
     .     btor_rf, rmaxis, zmaxis, psimax_rf, psilim_rf,
     .     rgrid_rf, nr_rf, zgrid_rf, nz_rf,
     .     sgnpsi, psi_rf, fpsi_rf, xbouni, dx2i, q_rf,
     .     bmini_rf, bmaxi_rf)
c
      deallocate (psi_rf, STAT = istat)
      if(istat .ne. 0)
     .          call deallocate_error("psi_rf",0,istat) 
c                        

      return
c
      end








      subroutine rfwave (n41, abd, bcol, ipvt, work, icheck)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     subroutine to calculate the field distribution along the major radius
c     by inversion of band matrix
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      integer ipvt(n41),kflag(2)
      complex*16 bcol(n41),work(n41),abd(16,n41),rfmat
      complex*16 rfa,rfb,rfc,ey,eytot,xkie,efld,dielt
      common /matx/rfa(kich,4),rfb(kich,4),rfc(kich,4)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      common /rfp/efld(kich,3,4),dielt(6,2),exy2(kich,2),
     .           ezy2(kich,2),elr2(kich,2)
      common /ichflx/xkie(kich,4),ey(kich,4),eytot(kich)
c
c assuming incident fast wave only
c assuming reflected wave energy from the n0n-incident side to ne zero
c setting up matrix equation to solve for the wave frofile
c rewriting the matrix in band storage mode
c
      icheck   = 0
      kflag(1) = 0
      kflag(2) = 0
      rik1     = PIMAG (xkie(1,2))
      rik2     = PIMAG (xkie(ngrid,2))
      go to (301, 302), kinc
c
  302 if (rik1 .lt. 0.0)  kflag(1) = 1
      if (rik2 .lt. 0.0)  kflag(2) = 1
      go to 303
c
  301 if (rik1 .gt. 0.0)  kflag(1) = 1
      if (rik2 .gt. 0.0)  kflag(2) = 1
c
  303 ml = 5
      mu = 5
      m  = ml + mu + 1
      do j=1,n41
        i1 = MAX0 (  1, j-mu)
        i2 = MIN0 (n41, j+ml)
        do i=i1,i2
          k        = i - j + m
          abd(k,j) = rfmat(i,j,kflag)
        end do
      end do
      lda = 2*ml + mu + 1
      call cgbco (abd, lda, n41, ml, mu, ipvt, rcond, work)
      rcond1 = 1.0+rcond
****  if (rcond1 .eq. 1.0)  go to 50
c
c setting up the column vector on the right hand side
c
      bcol(1) = -rfc(1,1)
      bcol(2) = -rfc(1,1)*efld(1,3,1)
      bcol(3) = -rfa(1,1)*rfc(1,1)
      bcol(4) = -rfb(1,1)*rfc(1,1)
      do i=5,n41
        bcol(i) = CMPLX (0.0, 0.0)
      end do
c
c calling subroutine to solve band matrix equation
c
      call cgbsl (abd, lda, n41, ml, mu, ipvt, bcol, 0)
c
c setting field profile from solution
c
      ey(1,1) = CMPLX (1.0, 0.0)
      if (kflag(1) .gt. 0)  go to 101
      ey(1,2) = bcol(1)
      ey(1,3) = bcol(2)
      ey(1,4) = CMPLX (0.0, 0.0)
      go to 100
  101 ey(1,2) = CMPLX (0.0, 0.0)
      ey(1,3) = bcol(1)
      ey(1,4) = bcol(2)
c
  100 do i=2,ngrid-1
        ii = 4*(i-2)+2
        do k=1,4
          ey(i,k) = bcol(ii+k)
        end do
      end do
      ey(ngrid,1) = bcol(n41-1)
      ey(ngrid,3) = CMPLX (0.0, 0.0)
      if (kflag(2) .gt. 0)  go to 201
      ey(ngrid,2) = CMPLX (0.0, 0.0)
      ey(ngrid,4) = bcol(n41)
      go to 200
  201 ey(ngrid,2) = bcol(n41)
      ey(ngrid,4) = CMPLX (0.0, 0.0)
c
  200 do i=2,ngrid
        do k=1,4
          ey(i,k) = ey(i,k) * EXP (CMPLX (0.0, 1.0) * xkie(i,k) *
     .                            (rmrf(i) - rmrf(i-1)) * 0.5)
        end do
      end do
      return
c
   50 icheck = 1
      return
c
      end


      subroutine rwrt_gafit_in(ntoray,runid,nout,ncrt)
c-------------------------------------------------------------------------
c     this rotuine is called only if toray_version > toray_version_switch
c --------------------------------------------------------------HSJ-03/06/03
        USE ext_prog_info, ONLY : toray_version,toray_version_switch
        implicit none
        integer ntoray,nout,ncrt,rwrt,ipsi,npts
        character *(*) runid
        real *8 gafsep,dsrat,tolmap,percenflux
        logical newbdry,ex
        namelist /fitdat/ ipsi, npts, gafsep, dsrat, tolmap, percenflux,
     &                   newbdry

      ipsi = -100 
      ntoray = 101
      call getioun(ntoray,ntoray)
      inquire(FILE ='gafit.in',EXIST  = ex)
      if(ex)then
            print *,'gafit.in  exists'
            !read gafit.in
            open (unit = ntoray, file = 'gafit.in', status = 'OLD')
            read (unit = ntoray, fmt = '(3x, a72)') runid
            read (unit = ntoray, nml = fitdat)
            close (unit = ntoray)
            print *,'read gafit.in in Onetwo'
            if(toray_version .lt. toray_version_switch)then
               if(ipsi .ne. -100)then
                 write(nout,2)
                 write(ncrt,2)
 2               format(2x,'Input in gafit.in not consistent with ',
     .               ' required Onewo setup.',/,
     .               2x,'File gafit.in has  ipsi set ',/,
     .               2x,'But the selected version of gafit',/,
     .               2x,'does not accept this parameter')
                     call STOP(
     .                 'ERROR, gafit.in  namelist not compatible',0)
               endif
            else     !toray_version toray_version_switch
               If (ipsi .ne. 0)then
                 write(nout,1)
                 write(ncrt,1)
 1               format(2x,'Input in gafit.in not consistent with ',
     .               ' required Onewo setup.',/,
     .               2x,'File gafit.in must specify ipsi = 0 ',/,
     .               2x,'(Because we have to read file psiin)',/,
     .               2x,' We do not overide the missing  ',/,
     .               2x,'value since that might lead to user '
     .               'confusion' )
                     call STOP(
     .                 'ERROR, gafit.in  namelist not compatible',0)

                endif
            endif
      else    

           print *,'gafit.in doesnt exist'
           if(toray_version .gt. toray_version_switch)then
              !create a stripped down version of gafit.in
              !that just has ipsi =0 set:
              ipsi = 0
              call gafit_write(ipsi,ntoray,runid)
              print *,'file gafit.in created'
           endif
      endif
      call giveupus(ntoray)
      return
      end


      subroutine rwrt_toray_in1(igafit,ntoray,runid,
     .                          nout,ncrt,rwrt)
       USE precision_mod    ! this file and (comnam.i below) are
                            ! copied from the toray cvs tree
c-----------------------------------------------------------
c---     this subroutine isolates  the toray include files 
c-----------------------------------------------------------
         integer ntoray,istat,ISHELL,rwrt,nout,ncrt
         character *(*)runid
         !these *.i file are  copied from toray distribution
         include 'paradm.i'
         include 'comio.i'
         include 'comnam.i'    !edata namelist from here ,from toray cvs tree     
         raypatt=' '
         if(rwrt .eq. 1)then
            open (unit = ntoray, file = 'toray.in', status = 'OLD')
            read (unit = ntoray, fmt = '(3x, a72)') runid
            read (unit = ntoray, nml = edata)
            close (unit = ntoray)
            print *,'read toray.in in Onetwo'
            !return the value if igafit (if it was set in namelist) 
         else
            write(ncrt,1)igafit
            write(nout,1)igafit
 1           format(2x,'Error, file toray.in must have igafit = ',i3,
     .              /,'for this case to work correctly.')
             call STOP('Onetwo run terminated by bad Toray input',0)
            open  (unit = ntoray, file = 'toray.in', status = 'NEW')
            write (unit = ntoray, fmt = '(3x, a72)') runid
            write (unit = ntoray, nml = edata)
            close (unit = ntoray)
         endif

         return
      end




      subroutine  rwrt_toray_in
c-------------------------------------------------------------------
c     This subroutine checks to see if toray.in exists in the
c     local path. if it does, the namelist is read.

c      For toray versions that expect gafit to be external to toray
c      the variable igafit is not part of the namelist.
c      If it is set in the namelist inform the user and stop.
c      If it is not set in the namelist or toray.in doesnt 
c      exist  then we continue.


c      For toray versions that expect gafit to be called from
c      within toray igafit must be set to 1. (The default
c      in the new versions of toray is 0). So if toray .in exists
c      read it, and check the value of igafit. If igafit is 0,
c      inform the user and quit. If igafit =1 then continue.
c      If toray.in doesnt exist, create a stripped down version that
c      contains only the statment "igafit =1" and continue.
c----------------------------------------------------------------HSJ

      USE param
      USE ext_prog_info, ONLY : toray_version,toray_version_switch
      USE io 
      implicit  integer (i-n), real*8 (a-h, o-z)
      integer iostat,ntoray,istat,ISHELL,rwrt
      logical ex
      character(len = 72)runid


 100  igafit = -100           !see if this value changes by namelist read below

      ntoray = 101
      call getioun(ntoray,ntoray)
      inquire(FILE ='toray.in',EXIST  = ex)
      if(ex)then
         print *,'toray.in exists,toray_version =',toray_version
         rwrt = 1
         call rwrt_toray_in1(igafit,ntoray,runid,nout,ncrt,rwrt)
         if(toray_version .lt. toray_version_switch )then
              if(igafit .gt. -100  )then
                 !toray.in exists and igafit was set . This means that
                 !toray.in was created for toray versions > 1.4 
                 write(nout,1)
                 write(ncrt,1)
 1               format(2x,'Input in toray.in not consistent with ',
     .             ' version of toray being run:',/,
     .             2x,'File toray.in has specified igafit',/,
     .             2x,'But  only toray versions greater than 1.4 ',/,
     .             2x,'accept this namelist variable' )
                 call STOP(
     .           'ERROR,toray.in  namelist not compatible',0) 
              endif
               !igafit not set so toray.in is ok
               !skip down to check gafit.in
         else  !toray_version >= toray_version_switch
              if(igafit .gt. -100  )then         !igafit was set in namelist read
                                                 !if value is not 1 then error
                  if(igafit .ne. 1)then
                     write(nout,2)igafit
                     write(ncrt,2)igafit
 2                   format(2x,'Input in toray.in not consistent with ',
     .               ' required Onewo setup.',/,
     .               2x,'File toray.in has specified igafit =',i5,/,
     .               2x,'But we have to run gafit and hence ',
     .               'igafit =1 '/,
     .               2x,'is required. We do not overide the '
     .               'incorrect ',/,
     .               2x,'value since that might lead to user '
     .               'confusion' )
                     call STOP(
     .               'ERROR, toray.in  namelist not compatible',0)
                  endif
               else    !igafit .eq.-100
                       !namelist toray.in exists but igafit is not set
                       !this is an error condition:
                     write(nout,3)
                     write(ncrt,3)
 3                   format(2x,'Input in toray.in not consistent with ',
     .               ' required Onewo setup.',/,
     .               2x,'File toray.in must specify the required '
     .               'input igafit = 1 ',/,
     .               2x,'Because we have to run gafit',/,
     .               2x,' We do not overide the missing  ',/,
     .               2x,'value since that might lead to user '
     .               'confusion' )
                       call STOP(
     .                 'ERROR, toray.in  namelist not compatible',0)
               endif
         endif
       else       !file toray.in doesnt exist. 
          print *,'toray.in doesnt exist,toray_version=',toray_version
c          if(toray_version .gt. toray_version_switch)then  
          if(toray_version .GE. toray_version_switch)then  ! changed 7/24/2011 HSJ
                  !create a stripped down toray.in that will have igafit = 1
            igafit = 1
            runid ='toray.in created by Onetwo'      !  if changed also change above
            call wrt_toray_in2(igafit,ntoray,runid)  !  want to redefine  the namelist 
            print *,'file toray.in created'  
          endif
       endif

       call giveupus(ntoray)
C      now check gafit.in:
       
       rwrt = 1
       call rwrt_gafit_in(ntoray,runid,nout,ncrt)

      return
      end

      subroutine switch (r, n, work)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- subroutine SWITCH reverses the order of the contents of array r
c
      dimension  r(n), work(n)
c
      do 100 i=1,n
  100 work(i) = r(i)
      do 110 i=1,n
  110 r(i) = work(n-i+1)
      return
c
      end


      subroutine TorGA_get_fspt2 (psi1d, farray, npsi1d, fspt2)
cProlog

c
      implicit none
c
      integer npsi1d
      real *8  psi1d(*), farray(*), fspt2(*)
c
c     uses subroutine CSPLINE
c
c     local variables:
c
      doubleprecision yp1, yp2
c use natural b.c.
      data yp1, yp2 /1.0d30, 1.0d30/ 
c
      call cspline (psi1d, farray, npsi1d, yp1, yp2, fspt2)
      return
c
      end

       subroutine TorGA_psigrid (psimax, psilim, npsi1d, psi1d)
cProlog
 
c
       implicit none
c
       integer npsi1d
       real*8  psi1d(*), psimax, psilim
c
       integer j
       real*8 dpsi
c
       dpsi=(psilim-psimax)/(npsi1d-1.d0)
       psi1d(1)=psimax
       psi1d(npsi1d)=psilim
       do j=2,npsi1d-1
       psi1d(j)=psimax+(j-1)*dpsi
       enddo
       return
c
       end

      subroutine trpic2 (dni, zi, gmi, temi, te, b0, x0, xf, nsp, nslb,
     .                   wkpl, frq, pabs, pabsi, pabse, pabsie, pabst,
     .                   pleft, pitot, petot)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  subroutine trpic2 is an alternative heating model to subplement
c  the wave operator heating code at high beta.
c  it uses cold plasma dispersion to calculate the real part of
c  the perpendicular wave number, then uses this kx to calculate
c  the imaginary part of kx
c  isw = 2 is low field incidence; then xslb=thickness of slab is negative;
c  also kxr is negative;
c  isw = 1 is high field incidence, the signs are positive.
c
c   input parameters:
c      dni(nsp,nslb),zi(nsp),gmi(nsp),ti(nsp,nslb),te(nslb),
c      b0(nslb) = b-field at each slab;
c      x0 = starting position with respect to plasma axis,(low field side is
c               positive),
c      xf = final position with respect to plasma axis,
c      nsp = number of ion-species
c      nslb = number of slabs
c
c   output quantities:
c      pabs(nsp,nslb) = fractional power absorbed by each ion species
c                       at each slab
c      pabsi(nslb) = fractional power absorbed by all ion species at each slab
c      pabse(nslb) = fractional power absorbed by electrons at each slab
c      pabsie(nslb) = fractional power consumed at each slab
c      pitot = total power absorbed by the ions
c      petot = total power absorbed by the electrons
c      pabst = total power absorbed
c      pleft = fractional power left
c
      include 'ichp2.i'
c
      parameter (nx = 2, nz = 3, nne = 1, nkxi = 3, ngam = 3)
      real*8     kxr2(ny),kxr(ny),kxi(nkxi,ny),kxim(ny)
      dimension  temi(ny),pabsi(ny)
      dimension  dni(nx,ny),dne(ny),zi(nx),gmi(nx),ti(nx,ny),
     .           te(ny),b0(ny),fci(nx,ny),fce(ny),vi(nx,ny),
     .           ve(ny),rhoi(nx,ny),fpi(nx,ny),fpe(ny),li(nx),
     .           zips(nx,ny,nz),zi0(nx,ny),zims(nx,ny,nz),zeps(ny,nne),
     .           zems(ny,nne),ze0(ny),s(ny),d(ny),
     .           gam(nx,ngam,ny),game(ny)
      dimension  sumk(ny),rm(ny),pabs(nx,ny),pabse(ny),pabsie(ny)
      data       pi/3.14159265/, c/2.9979e10/
c
      anslb = nslb
      xslb  = (xf - x0) / anslb
      if (xslb .lt. 0.0)  isw = 2
      if (xslb .gt. 0.0)  isw = 1
      le = 1
c
      do n=1,nsp
        li(n) = 3
      end do
c
      do l=1,nslb
        do n=1,nsp
          ti(n,l) = temi(l)
        end do
      end do
c
      wk0  = frq/c
      wk02 = wk0**2
      wnz  = wkpl/wk0
      wnz2 = wnz**2
      rpi  = SQRT (pi)
      call plsmap (dni,dne,zi,gmi,ti,te,frq,b0,nsp,nslb,
     .             fci,fce,vi,ve,rhoi,fpi,fpe)
      call param1 (frq,vi,ve,li,le,nsp,nslb,wkpl,fci,fce,
     .             zips,zi0,zims,zeps,zems,ze0,3)
      call relkx2(fpi,nsp,nslb,fpe,fci,fce,kxr2,frq,s,d,wkpl)
      call akim(fpi,rhoi,nsp,nslb,vi,li,zims,frq,wkpl,gam,kxr2,3)
c
      do i=1,nslb
        kxr(i)  = SQRT (kxr2(i))
        if (isw .eq. 2)  kxr(i) = -kxr(i)
        game(i) = -0.5*rpi*(ve(i)*fpe(i)/(c*fce(i)))**2*ze0(i)*
     .             EXP (-ze0(i)**2)
      end do
c
c the nsp+1-th species in kxi(n,nslb) is the electron contribution
c
      do l=1,nslb
        snz = s(l)-wnz2
        do n=1,nsp
          kxi(n,l) = wk0**2*((gam(n,1,l)+gam(n,2,l))*snz-2.0*gam(n,3,l)*
     .               d(l))-kxr2(l)*gam(n,1,l)
          kxi(n,l) = kxi(n,l)*0.5/(kxr(l)*snz)
        end do
        kxi(nsp+1,l) = -kxr(l)*game(l)
      end do
c
c kxim(i) is the sum of kxi of all species
c
      do l=1,nslb
        kxim(l) = 0.0
        do n=1,nsp+1
          kxim(l) = kxim(l) + kxi(n,l)
        end do
      end do
c
c sumk(nslb) is sum of damping exponent
c
      sumk(1) = 0.0
      rm(1)   = 1.0
      do l=2,nslb
        sumk(l) = sumk(l-1)+2.0*kxim(l-1)*xslb
        rm(l)   = EXP (-sumk(l))
      end do
      pitot = 0.0
      petot = 0.0
      pabst = 0.0
c
      do n=1,nsp
        pabs(n,1) = 0.0
      end do
c
      pabse(1) = 0.0
c
      do l=2,nslb
        pabsie(l) = 0.0
        do n=1,nsp
          pabs(n,l) = kxi(n,l-1)*(rm(l-1)-rm(l))/kxim(l-1)
          pabsie(l) = pabsie(l)+pabs(n,l)
          pabst = pabst+pabs(n,l)
        end do
        pabsi(l) = pabsie(l)
        pabse(l) = kxi(nsp+1,l-1)*(rm(l-1)-rm(l))/kxim(l-1)
        pitot = pitot+pabsi(l)
        petot = petot+pabse(l)
        pabsie(l) = pabsie(l)+pabse(l)
        pabst = pabst+pabse(l)
      end do
c
      if (pabst .gt. 1.0) then
        do l=1,nslb
          do n=1,nsp
            pabs(n,l) = pabs(n,l)/pabst
          end do
          pabse (l) = pabse (l) / pabst
          pabsi (l) = pabsi (l) / pabst
          pabsie(l) = pabsie(l) / pabst
        end do
        pabst = 1.0
      end if
c
      pleft = 1.0 - pabst
      return
c
      end

      subroutine wrt_psiin (nunit,lcentr, ledge, nr, nz, rdim, zdim,
     .                       rcenter, rinside, rmaxis, zmaxis,
     .                       psimax, psilim, b0, gasep, sf, psi,
     .                       nw, psir, qpsi, npbdry, nlimtr, rpbdry,
     .                       zpbdry, rlim, zlim,h_factr_rf,b_avg_rf,
     .                       bsq_avg_rf,r0rinv_rf)
c
      implicit none
c
      integer nunit
      integer lcentr, ledge, nr, nz, npbdry, nlimtr
      real*8    rdim, zdim, rcenter, rinside, rmaxis, zmaxis
      real*8    psimax, psilim, b0
      integer nw
      real*8    gasep
      real*8    sf(*), psi(nw,*), psir(*)
      real*8    qpsi(*)
      real*8    rpbdry(*), zpbdry(*), rlim(*), zlim(*)
      REAL *8   h_factr_rf(*),b_avg_rf(*),bsq_avg_rf(*)
      REAL *8   r0rinv_rf(*)
****  real*8    rlimmin, rlimmax, zlimmin, zlimmax
c
c --- write 'psiin' (an input file for GAFIT)
c
      integer i, j
c
      open (unit = nunit, file = 'psiin', status = 'unknown')
c
      write   (nunit, '(4i4)')  lcentr,ledge,nr,nz
      write   (nunit,   200  )  rdim,zdim,rcenter,rinside
      write   (nunit,   200  )  rmaxis,zmaxis,psimax,psilim,b0
      write   (nunit,   200  )  gasep
      write   (nunit,   200  )  (sf(i),i=1,nr)
      write   (nunit,   200  )  ((psi(i,j),i=1,nr),j=1,nz)
      write   (nunit,   200  )  (psir(i),i=1,ledge)
      write   (nunit,   200  )  (qpsi(i),i=1,nr)
      write   (nunit, '(2i5)')   npbdry,nlimtr
      write   (nunit,   200  )  (rpbdry(i),zpbdry(i),i=1,npbdry)
      write   (nunit,   200  )  (rlim(i),zlim(i),i=1,nlimtr)
      write   (nunit,   200  )  (h_factr_rf(i),i =1,nr) !<SQRT(1-B/Bmax)>
      write   (nunit,   200  )  (bsq_avg_rf(i),i =1,nr) !<B882/B0**2>
      write   (nunit,   200  )  (b_avg_rf(i),i =1,nr)   !<B/B0>, B0 == Btor
      write   (nunit,   200  )  (r0rinv_rf(i),i =1,nr)   !<R/R0> 
      close  (unit = nunit)

  200 format (5e16.9)
      return
      end


      subroutine wrt_rf_mhd (n_rf,nj,
     .     rgmax, rgmin, zgmax, zgmin, rpout, rpin, zptop, zpbot,
     .     bzero, rmaxis, zmaxis, psimax, psilim,
     .     rgrid, nr, zgrid, nz, 
     .     sgnpsi, psir, fpsi, xbouni, dx2i, qsafetyi, bmini, bmaxi)
c----------------------------------------------------------------HSJ----
      USE ext_prog_info, only : nw_rf,nh_rf,kj_rf,nj_rf
      USE mhdpar
      USE bicube
      implicit none
c
      integer n_rf,npsimax,nj,npsi1d,ledge,j
c      include 'mhdpar.i'  !dims for bicube.i
c      include 'bicube.i'  !get cspln from here
      real*8    rgmax,rgmin,zgmax,zgmin,rpout,rpin,zptop,zpbot
      real*8    bzero,rmaxis,zmaxis,psimax,psilim,psiout
      integer   nr,nz
      integer,  allocatable,dimension(:) :: ltab
      real*8    rgrid(nw),zgrid(nh),psi2d_rf(nw,nh),sgnpsi
      real*8    psir(nj), fpsi(nj)
      real*8    xbouni(nj),dx2i(nj),qsafetyi(nj),
     .          bmini(nj),bmaxi(nj)
      real*8,   dimension(:), allocatable :: psi1d,fspt2
      character *8  outfile
      logical   exists

      allocate(psi1d(nw))
      allocate(fspt2(nw))
      allocate(ltab(kj_rf))


      outfile = 'mhddat'
      psiout =  psilim
      ledge = nw-1 ! # of psi zones  ??
c      call get_newpsir (psimax, psilim, psiout, psir, lcentr, ledge,
c     .                  gasep, psirnew, xbouni, ltab, ledgenew)

      nj_rf = nj
      npsimax = 65                     !npsimax follows toray convention
      if(nh_rf .ge. 129) npsimax =129
      npsi1d = nw
      call TorGA_psigrid (psimax, psilim, npsi1d, psi1d)
c     fspt2 are the spline coefficients for fpsi on the psi1d grid
c     which is of size npsi1d = nw.
      call TorGA_get_fspt2 (psi1d, fpsi, npsi1d, fspt2)
c---open a binary file for output that will be read by the rf codes
c
      inquire (file = outfile, exist = exists)
      if (exists)   call DESTROY(outfile)
      call getioun(n_rf,n_rf)
      open  (unit = n_rf, file = outfile, status = 'NEW',
     .        form='UNFORMATTED')
       write (unit = n_rf) nw_rf,nh_rf,2*nh_rf,npsimax,nj_rf
       write (unit = n_rf) rgmax,rgmin,zgmax,zgmin,rpout,rpin,
     .                                                  zptop,zpbot
       write (unit = n_rf) rmaxis,zmaxis,bzero,psimax,psilim,psiout
       write (unit = n_rf) rgrid,nr,zgrid,nz,cspln,sgnpsi
       write (unit = n_rf) psi1d,fpsi,fspt2,npsi1d
       write (unit = n_rf) ledge,ltab
       write (unit = n_rf) xbouni,dx2i,qsafetyi,bmini,bmaxi
      call giveupus(n_rf)
      close (unit = n_rf)
c      print *,'psi1d =',psi1d
      deallocate(psi1d)
      deallocate(fspt2)
      deallocate(ltab)

      call STOP('temp',0)
      return
c
      end


      subroutine wrt_toray_in2(igafit,ntoray,runid)
c--------------------------------------------------------------
c -- write toray.in file ( because none was found
c -- in the current path.
c----------------------------------------------------------HSJ-
       USE rf,                           ONLY : ech_input,rfpow    
   
       USE echdat_module,         ONLY : gauszone_lcl  => gauszone,
     .                                   nharm_lcl     => nharm,    
     .                                   netcdfdat_lcl => netcdfdat,
     .                                   modelc_lcl    => modelc,
     .                                   model_globl

       USE io,                    ONLY : io_toray_hist
       IMPLICIT NONE

       REAL*8   ergspj,powinc
       INTEGER  igafit,ntoray
       LOGICAL  nlout(2)
       INTEGER  gauszone,nharm,netcdfdat,modelc
       character*(*) runid
 


       namelist /edata/ igafit,gauszone,nharm,netcdfdat,
     .                  powinc,modelc
!     .                  nlout causes problems with toray.nc write
!      set namelist variables. if ech_input = none then the **_lcl values
!      are not known and are defaulted to likely values!!!!!!!!!!!
       ergspj    = 10000000.
       igafit    = 1
       nlout(:)  = .FALSE.
       IF(ech_input .NE. 'none')THEN
          gauszone  = gauszone_lcl(model_globl)
          nharm     = nharm_lcl(model_globl)  
          netcdfdat = netcdfdat_lcl(model_globl)
          powinc    = rfpow(model_globl) * ergspj ! powinc in ergs/sec
          modelc    = modelc_lcl(model_globl)
       ELSE
          gauszone  = 4
          nharm     = 2
          powinc    = rfpow(model_globl)*ergspj
          netcdfdat = 1
          modelc    = 5
       ENDIF


       ! load edata namelist variables from echin.nc if available.
       ! otherwise just put igafit =1 in namelist
       open (unit = ntoray, file = 'toray.in', status = 'NEW')
       write (unit = ntoray, fmt = '(3x, a72)') runid
       IF(ech_input .EQ. 'none')write(UNIT = ntoray, 
     .      fmt = '(3x,"file ech_input not present",/,
     .      " Therefore the following namelist is a ",/
     .      "  generic default,except for rfpow" )')


       write (unit = ntoray, nml = edata)
       close (unit = ntoray)

       write(io_toray_hist,FMT='("toray.in loaded with:")')
       write(io_toray_hist,FMT='("model,rfpow =",i5,x,1pe14.6)')
     .                            model_globl,powinc
       write(io_toray_hist,FMT='("---------------------------",/)')

       return
       end







      subroutine zanlyt1 (f, eps, nsig, kn, nguess, n, x, itmax, infer,
     .                    ier, ia)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c -------------- ZANLYT1 ------- s ------- library 3 -------------------
c
c   function            - determination of zeros of an analytic complex
c                           function using muller's method with
c                           deflation
c   usage               - call zanlyt1 (f,eps,nsig,kn,nguess,n,x,itmax,infer,
c                                       ier, i)
c   parameters   f      - a complex function subprogram, f(z), written
c                           by the user specifying the equation whose
c                           roots are to be found.  f must appear in
c                           an external statement in the calling pro-
c                           gram.
c                eps    - first stopping criterion.  let fp(z) = f(z)/p
c                           where p = (z-z(1))*(z-z(2))*,,,*(z-z(k-1))
c                           and z(1),...,z(k-1) are previously found
c                           roots.  if ((ABS (f(z)) .le. eps)  .and.
c                           (ABS (fp(z)) .le. eps)), then z is accepted
c                           as a root.
c                nsig   - 2nd stopping criterion.  a root is accepted
c                           if two successive approximations to a given
c                           root agree in the first nsig digits. (input)
c                             note. if either or both of the stopping
c                             criteria are fulfilled, the root is
c                             accepted.
c                kn     - the number of known roots which must be stored
c                           in x(1),...,x(kn), prior to entry to zanlyt1
c                nguess - the number of initial guesses provided. these
c                           guesses must be stored in x(kn+1),...,
c                           x(kn+nguess).  nguess must be set equal
c                           to zero if no guesses are provided. (input)
c                n      - the number of new roots to be found by
c                           zanlyt1 (input)
c                x      - a complex vector of length kn+n.  x(1),...,
c                           x(kn) on input must contain any known
c                           roots.  x(kn+1),..., x(kn+n) on input may,
c                           at the user]s option, contain initial
c                           guesses for the n new roots which are to be
c                           computed.  if the user does not provide
c                           an initial guess, zero is used.
c                           on output, x(kn+1),...,x(kn+n) contain the
c                           approximate roots found by zanlyt1.
c                itmax  - the maximum allowable number of iterations
c                           per root (input)
c                infer  - an integer vector of length kn+n.  on
c                           output infer(j) contains the number of
c                           iterations used in finding the j-th root
c                           when convergence was achieved.  if
c                           convergence was not obaxis(i)ained in itmax
c                           iterations, infer(j) will be greater than
c                           itmax (output).
c                ier    - error parameter (output)
c                         warning error
c                           ier = 33 failure to converge within itmax
c                           iterations for at least one of the (n) new
c                           roots.
c               i      - grid number
c   precision           - single
c   reqd. IMSL routines - uertst
c   language            - fortran
c ----------------------------------------------------------------------
c   latest revision     - february 5, 1974
c
      external F
c
      include 'ichpar.i'
      include 'ichcom.i'
c
      logical       false                         ! added 7 Aug 96 by JF
      dimension     infer(*), x(*)
      complex*16    x,d,dd,den,fprt,frt,h,rt,t1,t2,t3,
     .              tem,z0,z1,z2,bi,f,xx,xl,y0,y1,y2,x0,
     .              zero,p1,one,four,p5,rtt
      complex*16    efld,dielt
      complex*16    skxx,skxy,skxz,skyy,skyz,skzz
      complex*16    sxx,sxy,sxz,syy,syz,szz
      common /diel/ skxx(kich,ks,4),skxy(kich,ks,4),skxz(kich,ks,4),
     .              skyy(kich,ks,4),skyz(kich,ks,4),skzz(kich,ks,4)
      common /rf2 / sxx(5),sxy(5),sxz(5),syy(5),syz(5),szz(5)
      common /rfp / efld(kich,3,4),dielt(6,2),exy2(kich,2),
     .              ezy2(kich,2),elr2(kich,2)
c
      common /rfl   / pi, wfreq, clight, ergtkev, charge, emass(5),
     .                anumb(5),
     .                neq, nhp, nh2p
c
      common /rf4 / xnparz, xnpar, dbg(99)
      data          zero/(0.0,0.0)/, p1/(0.1,0.0)/,
     .              one/(1.0,0.0)/, four/(4.0,0.0)/,
     .              p5/(0.5,0.0)/
      data          rzero/0.0/, rten/10.0/, rhun / 100.0/, rp01/0.01/
      data          ax/0.1/,ickmax/3/
      data          false /.false./               ! added 7 Aug 96 by JF
c
      real_zero = 0.0
      ier   = 0
      if (n .lt. 1)  go to 9005
  101 eps1  = rten**(-nsig)
      eps1  = MIN (eps1, rp01)
c                                  set number of iterations
      knp1  = kn + 1
      knpn  = kn + n
      knpng = kn + nguess
      do i=1,knpn
        infer(i) = 0
        if (i .gt. knpng)  x(i) = zero
      end do
      l   = knp1
   10 jk  = 0
      ick = 0
      xl  = x(l)
   15 ic  = 0
      h   = ax
      h   = p1*h
      if (ABS (xl) .gt. ax)  h = p1*xl
c
c     The IF block below was added 7 Aug 96 by JF as a cheap workaround:
c
      if (false) then       ! NOTE that variable FALSE is always .false.
        assign 50 to nn
      end if
c
c     first three points are xl+h, xl-h, xl
c
      rt = xl+h
      assign 20 to nn
      go to 50
   20 z0 = fprt
      y0 = frt
      x0 = rt
      rt = xl-h
      assign 25 to nn
      go to 50
   25 z1 = fprt
      y1 = frt
      h = xl-rt
      d = h/(rt-x0)
      rt = xl
      assign 30 to nn
      go to 50
   30 z2 = fprt
      y2 = frt
c                                  begin main algorithm
   35 dd = one + d
      t1 = z0*d*d
      t2 = z1*dd*dd
      xx = z2*dd
      t3 = z2*d
      bi = t1-t2+xx+t3
      den = bi*bi-four*(xx*t1-t3*(t2-xx))
c
c     use denominator of maximum amplitude
c
      t1 = SQRT (den)
      qz = rhun * MAX (ABS (bi),ABS (t1))
      t2 = bi + t1
      if (ABS (t2)+qz .eq. qz) t2 = zero
      t3 = bi - t1
      if (ABS (t3)+qz .eq. qz) t3 = zero
      den = t2
      qz = ABS (t3)-ABS (t2)
      if (qz .gt. rzero)  den = t3
c
c     test for zero denominator
c
      assign 30 to nn
      if (ABS (den) .eq. rzero)  go to 65
      d = -xx/den
      d = d+d
      h = d*h
      rt = rt + h
c
c     check convergence of the first kind
c
      if (ABS (h) .le. eps1 * MAX (ABS (rt),ax))  go to 70
      if (ic .ne. 0)  go to 15
      assign 40 to nn
      go to 50
   40 qz = ABS (fprt)-ABS (z2)*rten
      if (qz .ge. rzero)  go to 45
      z0 = z1
      z1 = z2
      z2 = fprt
      y0 = y1
      y1 = y2
      y2 = frt
      go to 35
c
c     take remedial action to induce convergence
c
   45 d    = d*p5
      h    = h*p5
      rt   = rt-h
   50 jk   = jk+1
      if (jk .gt. itmax)  go to 75
      frt  = f(rt)
      fprt = frt
c
c     test to see if first root is being determined
c
      if (l .eq. 1)  go to 60
c
c     compute deflated function
c
      lm1 = l-1
      do 55 i=1,lm1
         tem = rt - x(i)
         if (ABS (tem) .eq. rzero)  go to 65
   55 fprt = fprt / tem
c
c     check convergence of the second kind
c
   60 if (ABS (fprt) .le. eps .and. ABS (frt) .le. eps)  go to 80
c
      go to nn, (20, 25, 50, 40)
c
   65 if (ic .ne. 0)  go to 15
      tem = rten*eps1
      if (ABS (rt) .gt. ax) tem = tem*rt
      rt  = rt+tem
      d   = (h+tem)*d/h
      h   = h+tem
      go to 50
c                                  check solution
   70 if (ic .ne. 0)  go to 80
      ic  = 1
      z0  = y1
      z1  = y2
      z2  = f(rt)
      xl  = rt
      ick = ick+1
      if (ick .le. ickmax)  go to 35
c                                         warning error, itmax = maximum
      jk  = itmax + jk
   75 ier = 33
c                                         a root has been found
   80 x(l)       = rt
      infer(l)   = jk
      err1       =  REAL (rt)
      go to 108
  108 err2       = PIMAG (rt)
      rtt        =  SQRT (rt)
      if (REAL (rtt) .le. 0.0)  rtt = -rtt
      dielt(1,l) = CMPLX (dbg(31), real_zero) + CMPLX (0.0, 1.0)*dbg(32)
     .           - CMPLX (xnparz , real_zero)
      dbg(71)    =  REAL (dielt(1,l))
      dbg(52)    = PIMAG (dielt(1,l))
      dielt(2,l) = CMPLX (dbg(27), real_zero) + CMPLX (0.0, 1.0)*dbg(28)
      dielt(3,l) = CMPLX (dbg(25), real_zero) + CMPLX (0.0, 1.0)*dbg(26)
     .           + rtt*xnparz/xkpar
      dbg(55)    =  REAL (dielt(3,l))
      dbg(56)    = PIMAG (dielt(3,l))
      dielt(4,l) = CMPLX (dbg(23), real_zero) + CMPLX (0.0, 1.0)*dbg(24)
     .           - rt*xnpar
     .           - CMPLX (xnparz , real_zero)
      dbg(57)    =  REAL (dielt(4,l))
      dbg(58)    = PIMAG (dielt(4,l))
      dielt(5,l) = CMPLX (dbg(19), real_zero) + CMPLX (0.0, 1.0)*dbg(20)
      dielt(6,l) = CMPLX (dbg(15), real_zero) + CMPLX (0.0, 1.0)*dbg(16)
     .           - rt*xnpar
      dbg(53)    =  REAL (dielt(6,l))
      dbg(54)    = PIMAG (dielt(6,l))
c
      do k=1,neq
        skxx(ia,k,l)   = sxx(k)
        skxy(ia,k,l)   = sxy(k)
        skxz(ia,k,l)   = sxz(k)
        skyy(ia,k,l)   = syy(k)
        skyz(ia,k,l)   = syz(k)
        skzz(ia,k,l)   = szz(k)
        skxx(ia,k,l+n) = skxx(ia,k,l)
        skxy(ia,k,l+n) = skxy(ia,k,l)
        skxz(ia,k,l+n) = -skxz(ia,k,l)
        skyy(ia,k,l+n) = skyy(ia,k,l)
        skyz(ia,k,l+n) = -skyz(ia,k,l)
        skzz(ia,k,l+n) = skzz(ia,k,l)
      end do
c
      l = l + 1
      if (l   .le. knpn)  go to 10
      if (ier .eq.    0)  go to 9005
****  call uertst (ier, 'zanlyt1')
 9005 return
c
      end

      subroutine zfunc (w, z, zp)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     calculates plasma dispersion function z = (zr,zi) and its
c     derivative zprime = (zpr,zpi) for the point u = (x,y)
c
      complex*16  z, zp, w
      data        sqrtpi /1.7724538509055/
c
      x   =  REAL (w)
      y   = PIMAG (w)
      call pfunc (x, y, zr, zi)
      xx  =  zr
      zr  = -sqrtpi * zi
      zi  =  sqrtpi * xx
      zpr = -2 * ( 1.0 + x*zr - y*zi)
      zpi = -2 * (y*zr + x*zi)
      z   = CMPLX (zr , zi )
      zp  = CMPLX (zpr, zpi)
      return
c
      end
