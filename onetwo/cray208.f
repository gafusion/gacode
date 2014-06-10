 

      subroutine allocate_error(var,myid,istat)
c------------------------------------------------------------------
      character *(*) var
      integer istat,myid
      write(*,1)var,istat,myid
 1    format(2x,"Memory Allocation error encountered",/,
     .       2x,"Unable to allocate ",a,/,
     .       2x,"status =",i5)
      istat =0 !reset for next case
      return
      end



      subroutine deallocate_error(var,myid,istat)
c---------------------------------------------------------------------
      character *(*) var
      integer istat, myid
      write(*,1)var,istat,myid
 1    format(2x,"Memory DE-Allocation error encountered",/,
     .       2x,"Unable to deallocate ",a,/,
     .       2x,"status =",i5, " process rank =",i5)
      istat =0 !reset for next case
      return
      end
      subroutine extrapbd (f1, f2, f3, f4, x1, y1, x2, y2, xt, yt,
     .                     xt1, yt1, xt2, yt2, psivl, area, dx, dy)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray208.f,v 1.69 2014/04/29 18:15:38 stjohn Exp $"/
c
c***********************************************************************
c**                                                                   **
c**     MAIN PROGRAM:  MHD FITTING CODE                               **
c**                                                                   **
c**     SUBPROGRAM DESCRIPTION:                                       **
c**          extrapbd extrapolates across a (x,y) cell.               **
c**                                                                   **
c**     RECORD OF MODIFICATION (Beware, this is very incomplete!):    **
c**          26/04/83..........first created                          **
c**          07/09/83..........replaced                               **
c**                                                                   **
c******************************************************************* HSJ
c
      dimension  xp(2), yp(2)
c
      ip    = 0
      fmin1 = MIN (f1,f2)
      fmax1 = MAX (f1,f2)
      fmin2 = MIN (f2,f3)
      fmax2 = MAX (f2,f3)
      fmin3 = MIN (f3,f4)
      fmax3 = MAX (f3,f4)
      fmin4 = MIN (f1,f4)
      fmax4 = MAX (f1,f4)
      a     = y2*(f2-f1)+y1*(f4-f3)
      b     = x1*(f2-f3)+x2*(f4-f1)
      c     = f1-f2+f3-f4
c
c ----------------------------------------------------------------------
c vertical asymptote is x = -b/c
c horizontal asymptote is y = -a/c
c there can only be vertical and horizontal asymptotes
c ----------------------------------------------------------------------
c
      if (c .eq. 0.0)  go to 950
      xasym = -b/c
      yasym = -a/c
      if ((   x1 .le. xasym) .and. (xasym .le. x2))  go to 840
      if ((yasym .lt. y1   ) .or.  (yasym .gt. y2))  go to 950
c
c ----------------------------------------------------------------------
c there is a horizontal asymptote
c find psi value on asymptote at x1 and x2
c ----------------------------------------------------------------------
c
      psias4 = ((f4-f1)*(yasym-y1)/dy)+f1
      psias2 = ((f3-f2)*(yasym-y1)/dy)+f2
      if (yt .lt. yasym)  go to 810
c
c ----------------------------------------------------------------------
c upper section of cell
c ----------------------------------------------------------------------
c
      fpmin4 = MIN (psias4,f4)
      fpmax4 = MAX (f4,psias4)
      if ((psivl .lt. fpmin4) .or. (psivl .gt. fpmax4))  go to 790
      ip = ip+1
      xp(ip) = x1
      yp(ip) = ((psivl-f1)/(f4-f1))*dy+y1
      if (psivl .eq. f4)  yp(ip) = y2
  790 fpmin2 = MIN (f3,psias2)
      fpmax2 = MAX (f3,psias2)
      if ((psivl .lt. fpmin2) .or. (psivl .gt. fpmax2))  go to 800
      ip = ip+1
      xp(ip) = x2
      yp(ip) = ((psivl-f2)/(f3-f2))*dy+y1
      if (psivl .eq. f3)  yp(ip) = y2
      if (ip .eq. 2)  go to 990
  800 ip = ip+1
      yp(ip) = y2
      xp(ip) = ((psivl-f4)/(f3-f4))*dx+x1
      go to 990
c
c ----------------------------------------------------------------------
c lower section of cell
c ----------------------------------------------------------------------
c
  810 fpmin4 = MIN (f1,psias4)
      fpmax4 = MAX (f1,psias4)
      if ((psivl .lt. fpmin4) .or. (psivl .gt. fpmax4))  go to 820
      ip = ip+1
      xp(ip) = x1
      yp(ip) = ((psivl-f1)/(f4-f1))*dy+y1
  820 fpmin2 = MIN (f2,psias2)
      fpmax2 = MAX (f2,psias2)
      if ((psivl .lt. fpmin2) .or. (psivl .gt. fpmax2))  go to 830
      ip = ip+1
      xp(ip) = x2
      yp(ip) = ((psivl-f2)/(f3-f2))*dy+y1
      if (ip .eq. 2)  go to 990
  830 ip = ip+1
      yp(ip) = y1
      xp(ip) = ((psivl-f1)/(f2-f1))*dx+x1
      go to 990
  840 if ((y1 .lt. yasym) .and. (yasym .lt. y2))  go to 900
c
c ----------------------------------------------------------------------
c vertical asymptote
c     find psi value on asymptote at y1 and y2
c ----------------------------------------------------------------------
c
      psias1 = ((f2-f1)*(xasym-x1)/dx)+f1
      psias3 = ((f3-f4)*(xasym-x1)/dx)+f4
      if (xt .lt. xasym)  go to 870
c
c ----------------------------------------------------------------------
c right side of cell
c ----------------------------------------------------------------------
c
      fpmin1 = MIN (psias1,f2)
      fpmax1 = MAX (psias1,f2)
      if ((psivl .lt. fpmin1) .or. (psivl .gt. fpmax1))  go to 850
      ip = ip+1
      yp(ip) = y1
      xp(ip) = ((psivl-f1)/(f2-f1))*dx+x1
      if (psivl .eq. f2)  xp(ip) = x2
  850 fpmin3 = MIN (f3,psias3)
      fpmax3 = MAX (f3,psias3)
      if ((psivl .lt. fpmin3) .or. (psivl .gt. fpmax3))  go to 860
      ip = ip+1
      yp(ip) = y2
      xp(ip) = ((psivl-f4)/(f3-f4))*dx+x1
      if (psivl .eq. f3)  xp(ip) = x2
      if (ip .eq. 2)  go to 990
  860 ip = ip+1
      xp(ip) = x2
      yp(ip) = ((psivl-f2)/(f3-f2))*dy+y1
      go to 990
c
c ----------------------------------------------------------------------
c left side of cell
c ----------------------------------------------------------------------
c
  870 fpmin1 = MIN (f1,psias1)
      fpmax1 = MAX (f1,psias1)
      if ((psivl .lt. fpmin1) .or. (psivl .gt. fpmax1))  go to 880
      ip = ip+1
      yp(ip) = y1
      xp(ip) = ((psivl-f1)/(f2-f1))*dx+x1
  880 fpmin3 = MIN (f4,psias3)
      fpmax3 = MAX (f4,psias3)
      if ((psivl .lt. fpmin3) .or. (psivl .gt. fpmax3))  go to 890
      ip = ip+1
      yp(ip) = y2
      xp(ip) = ((psivl-f4)/(f3-f4))*dx+x1
      if (ip .eq. 2)  go to 990
  890 ip = ip+1
      xp(ip) = x1
      yp(ip) = ((psivl-f1)/(f4-f1))*dy+y1
      go to 990
c
c ----------------------------------------------------------------------
c both horizontal and vertical asymptotes are present
c find psi on asymptotes
c ----------------------------------------------------------------------
c
  900 psias1 = ((f2-f1)*(xasym-x1)/dx)+f1
      psias2 = ((f3-f2)*(yasym-y1)/dy)+f2
      psias4 = ((f4-f1)*(yasym-y1)/dy)+f1
      psias3 = ((f3-f4)*(xasym-x1)/dx)+f4
      if (xt .gt. xasym)  go to 920
      if (yt .gt. yasym)  go to 910
c
c ----------------------------------------------------------------------
c xt .lt. xasym and yt .lt. yasym
c ----------------------------------------------------------------------
c
      fpmin1 = MIN (f1,psias1)
      fpmax1 = MAX (f1,psias1)
      yp(1)  = y1
      xp(1)  = ((psivl-f1)/(f2-f1))*dx+x1
      fpmin4 = MIN (f1,psias4)
      fpmax4 = MAX (f1,psias4)
      xp(2)  = x1
      yp(2)  = ((psivl-f1)/(f4-f1))*dy+y1
      go to 990
c
c ----------------------------------------------------------------------
c xt .lt. xasym and yt .gt. yasym
c ----------------------------------------------------------------------
c
  910 fpmin4 = MIN (f4,psias4)
      fpmax4 = MAX (f4,psias4)
      xp(1)  = x1
      yp(1)  = ((psivl-f1)/(f4-f1))*dy+y1
      fpmin3 = MIN (psias3,f4)
      fpmax3 = MAX (psias3,f4)
      yp(2) = y2
      xp(2) = ((psivl-f4)/(f3-f4))*dx+x1
      go to 990
  920 if (yt .gt. yasym)  go to 930
c
c ----------------------------------------------------------------------
c xt .gt. xasym and yt .lt. yasym
c ----------------------------------------------------------------------
c
      fpmin1 = MIN (f2,psias1)
      fpmax1 = MAX (f2,psias1)
      yp(1) = y1
      xp(1) = ((psivl-f1)/(f2-f1))*dx+x1
      fpmin2 = MIN (f2,psias2)
      fpmax2 = MAX (f2,psias2)
      xp(2) = x2
      yp(2) = ((psivl-f2)/(f3-f2))*dy+y1
      go to 990
c
c ----------------------------------------------------------------------
c xt > xasym and yt > yasym
c ----------------------------------------------------------------------
c
  930 fpmin2 = MIN (f3,psias2)
      fpmax2 = MAX (f3,psias2)
      xp(1) = x2
      yp(1) = ((psivl-f2)/(f3-f2))*dy+y1
      fpmin3 = MIN (f3,psias3)
      fpmax3 = MAX (f3,psias3)
      yp(2) = y2
      xp(2) = ((psivl-f4)/(f3-f4))*dx+x1
      go to 990
c
c ----------------------------------------------------------------------
c no asymptotes
c ----------------------------------------------------------------------
c
  950 if ((psivl .lt. fmin4) .or. (psivl .gt. fmax4))  go to 960
      ip = ip+1
      xp(ip) = x1
      yp(ip) = ((psivl-f1)/(f4-f1))*dy+y1
      if (psivl .eq. f4)  yp(ip) = y2
  960 if ((psivl .lt. fmin2) .or. (psivl .gt. fmax2))  go to 970
      ip = ip+1
      xp(ip) = x2
      yp(ip) = ((psivl-f2)/(f3-f2))*dy+y1
      if (psivl .eq. f3)  yp(ip) = y2
      if (ip .eq. 2)  go to 990
  970 if ((psivl .le. fmin1) .or. (psivl .gt. fmax1))  go to 980
      ip = ip+1
      yp(ip) = y1
      xp(ip) = ((psivl-f1)/(f2-f1))*dx+x1
      if (ip .eq. 2)  go to 990
  980 ip = ip+1
      yp(ip) = y2
      xp(ip) = ((psivl-f4)/(f3-f4))*dx+x1
  990 xt1 = xp(1)
      xt2 = xp(2)
      yt1 = yp(1)
      yt2 = yp(2)
      return
c
      end

      subroutine fillp (x, y, nold, nnew, dist, xwork, ywork, maxnew)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subroutine expands a set of x,y points defining a closed
c curve using linear interpolation so that the maximum distance
c between points is < = dist.
c ----------------------------------------------------------------------
c
      dimension  x(*), y(*), xwork(*), ywork(*)
c
      ncrt = 6
c
      do i=1,nold
        xwork(i) = x(i)
        ywork(i) = y(i)
      end do
c
      nold = nold+1
      xwork(nold) = xwork(1)
      ywork(nold) = ywork(1)
      nnew = 0
c
      do 4000 i=1,nold-1
        nnew = nnew+1
        if (nnew .gt. maxnew)  go to 9000
        x(nnew) = xwork(i)
        y(nnew) = ywork(i)
        del = (xwork(i+1)-xwork(i))**2+(ywork(i+1)-ywork(i))**2
        del = SQRT (del)
        if (del .le. dist)  go to 4000
        nadd = del/dist
        do k=1,nadd
          c = (k*dist)/del
          nnew = nnew+1
          if (nnew .gt. maxnew)  go to 9000
          x(nnew) = (1.0 - c)*xwork(i)+c*xwork(i+1)
          y(nnew) = (1.0 - c)*ywork(i)+c*ywork(i+1)
        end do
 4000 continue
c
      return
c
 9000 write  (ncrt, 8000)
 8000 format (/ ' sorry, ran out of space while expanding points',
     .          ' in subroutine FILLP')
      call STOP ('subroutine FILLP: unspecified problem', 13)
c
      end

      subroutine find2 (table, valin, n, iguess, ilow, ihigh, ncrt,nout)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subroutine uses a binary search algorithm to locate
c valin in a table
c ----------------------------------------------------------------------
c
      dimension table(*)
c
      valuse = valin
      if (valuse .gt. table(n))  valuse = table(n)
      if (valuse .lt. table(1))  valuse = table(1)
c
      ilow  = 0
      ihigh = n + 1
      if (iguess .le. 0)  iguess = 0.5 * (ilow+ihigh)
      if (iguess .gt. n)  iguess = 0.5 * (ilow+ihigh)
c
 2100 if (ilow .ge. ihigh-1)  return
      i = iguess
      if (table(i) .gt. valuse)  go to 2300
      if (table(i) .lt. valuse)  go to 2200
      ilow = iguess
      if (ilow .eq. n)  ilow = n - 1
      ihigh = ilow + 1
      return
c
 2200 ilow   = i
      iguess = 0.5 * (ilow+ihigh)
      go to 2100
c
 2300 ihigh  = i
      iguess = 0.5 * (ilow+ihigh)
      go to 2100
c
      end

      subroutine findeqlin (xaxd, yaxd, thet, a, bincp, iflg, isgn)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- find eq(uation)lin(e) finds the equation of the straight line
c --- passing through (xaxd,yaxd) and with slope whose angle is thet radians.
c --- input
c
c   xaxd
c   yaxd
c   thet
c
c --- output
c
c   a
c   bincp
c   iflg  =0 or 1 depending on how the line is defined, see below
c   isgn  +1 or -1,gives direction in which to move along the line to go
c           from (xaxd,yaxd) toward the boundary in direction thet
c
c ----------------------------------------------------------------------
c
      data pi,piov2,piov4,fpiov4,spiov4
     .    /3.141592654,1.570796327,0.7853981634,3.926990818,5.497787145/
      data twopi, tpiov4, tpiov2
     .    /6.283185307,2.356194491,4.712388904/
c
c --- get equation of ray emanating from (xaxd,yaxd)
c
      if (( piov4 .le. thet) .and. (thet .le. tpiov4))  go to 20
      if ((fpiov4 .le. thet) .and. (thet .le. spiov4))  go to 20
c
c --- y as a function of x            y = a*x+bincp
c
      isgn = -1
      if ((thet .lt. piov4) .or. (thet .gt. spiov4))  isgn = 1
      a = TAN (thet)
      iflg = 0
      bincp = yaxd-a*xaxd
      return
c
c --- x as a function of y            x = a*y+bincp
c
   20 isgn = 1
      if (thet .gt. pi)  isgn = -1
      if (isgn .eq. -1)  go to 22
      thet1 = piov2 - thet
      if (thet .gt. piov2)  thet1 = twopi - ABS (thet1)
      go to 25
   22 thet1 = tpiov2-thet
      if (thet .gt. tpiov2)  thet1 = pi - ABS (thet1)
   25 a     = TAN (thet1)
      iflg  = 1
      bincp = xaxd-a*yaxd
      return
c
c --- now have y = a*x+bincp    (iflg=0)   or
c --- x = a*y+bincp             (iflg=1)
c
      end

      subroutine fittocur (tocur, iounit)
c



c
c --- this subroutine samples the xchisq surface as a function of tocur
c --- near the given total current. The object is to modify tocur so
c --- that chisq is as small as possible. Note that only tocur is
c --- changed but chisq is calculated on the basis of the magnetic signals,
c --- with a small contribution from the current deviation. Note also
c --- that we do not try to converge to the minimum chisq very precisely.
c --- hence the simplistic approach given below is adequate.
c --- first we bracket a (local) minimum in chisq. Then we
c --- use parabolic extrapolation to zero in on the minimum more precisely.
c


      USE mhdpar      
      USE mhdbcdtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'mhdpar.i'
c      include 'mhdbcdtn.i'
c
      parameter (gold = 1.618304,glimit = 2.0,tiny = 1.0e-20,
     .           itmaxfit = 4,cgold=0.3819660,zeps = 1.0e-05)
c
c --- bracket the (local) minimum:
c
      tol     =  0.03
      indchi  =  0
      maxiteq = 20              ! max iterations allowed in FREEBDRY
      tocursv = tocur
      tocur   = 0.97*tocursv    ! assume 3% deviation increases xchisq
      ax      = tocur
      call freebdry (maxiteq, iounit, indchi, relerr)
      fa      = xchisq+((ax-tocursv)/(0.05*tocursv))**2
      tocur   = tocursv
      bx      = tocur
      call freebdry (maxiteq, iounit, indchi, relerr)
      fb      = xchisq+((bx-tocursv)/(0.05*tocursv))**2
      if (fb .gt. fa) then
        dum = ax
        ax  = bx
        bx  = dum
        dum = fb
        fb  = fa
        fa  = dum
      end if
      cx    = 2.0 * bx - ax
      tocur = cx
      call freebdry (maxiteq, iounit, indchi, relerr)
      fc    = xchisq+((cx-tocursv)/(0.05*tocursv))**2
    1 if (fb .ge. fc) then
        r = (bx-ax)*(fb-fc)
        q = (bx-cx)*(fb-fa)
        u = bx -
     .   ((bx-cx)*q-(bx-ax)*r) / (2.0 * SIGN (MAX (ABS (q-r),tiny),q-r))
        ulim = bx+glimit*(cx-bx)
        if ((bx-u)*(u-cx) .gt. 0.0) then
          tocur = u
          call freebdry(maxiteq,iounit,indchi,relerr)
          fu = xchisq+((u-tocursv)/(0.05*tocursv))**2
          if (fu .lt. fc) then
            ax = bx
            fa = fb
            bx = u
            fb = fu
            go to 1
          else if (fu .gt. fb) then
            cx = u
            fc = fu
            go to 1
          end if
          u = cx+gold*(cx-bx)
          tocur = u
          call freebdry(maxiteq,iounit,indchi,relerr)
          fu = xchisq+((u-tocursv)/(0.05*tocursv))**2
        else if ((cx-u)*(u-ulim) .gt. 0.0) then
          tocur = u
          call freebdry(maxiteq,iounit,indchi,relerr)
          fu = xchisq+((u-tocursv)/(0.05*tocursv))**2
          if (fu .lt. fc) then
            bx = cx
            cx = u
            u = cx+gold*(cx-bx)
            fb = fc
            fc = fu
            tocur = u
            call freebdry(maxiteq,iounit,indchi,relerr)
            fu = xchisq+((u-tocursv)/(0.05*tocursv))**2
          end if
        else if ((u-ulim)*(ulim-cx) .ge. 0.0) then
          u = ulim
          tocur = u
          call freebdry(maxiteq,iounit,indchi,relerr)
          fu = xchisq+((u-tocursv)/(0.05*tocursv))**2
        else
          u = cx+gold*(cx-bx)
          tocur = u
          call freebdry(maxiteq,iounit,indchi,relerr)
          fu = xchisq+((u-tocursv)/(0.05*tocursv))**2
        end if
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
        go to 1
      end if
c
c --- min chisq now equals fb with corresponding current bx
c --- next use parabolic extrapolation to refine the estimate:
c
      tol = 0.01    ! try for about 1% precision
      a = MIN (ax, cx)
      b = MAX (ax, cx)
      v = bx
      w = v
      x = v
      e = 0.0
      fx = fb
      fv = fx
      fw = fx
      do 11 iter=1,itmaxfit
        xm = 0.5 * (a+b)
        tol1 = tol * ABS (x)+zeps
        tol2 = 2.0 * tol1
        if (ABS (x-xm) .le. (tol2 - 0.5 * (b-a)))  go to 3
        if (ABS (e) .gt. tol1) then
          r = (x-w)*(fx-fv)
          q = (x-v)*(fx-fw)
          p = (x-v)*q-(x-w)*r
          q = 2.0 * (q-r)
          if (q .gt. 0.0)  p = -p
          q = ABS (q)
          etemp = e
          e     = d
          if (ABS (p) .ge. ABS (0.5 * q*etemp) .or. p .le. q*(a-x) .or.
     .    p .ge. q*(b-x))  go to 30
          d = p/q
          u = x+d
          if (u-a .lt. tol2 .or. b-u .lt. tol2)  d = SIGN (tol1,xm-x)
          go to 2
        end if
   30   if (x .ge. xm) then
          e = a-x
        else
          e = b-x
        end if
        d = cgold*e
    2   if (ABS (d) .ge. tol1) then
          u = x + d
        else
          u = x + SIGN (tol1, d)
        end if
        tocur = u
        call freebdry(maxiteq,iounit,indchi,relerr)
        fu = xchisq+((u-tocursv)/(0.05*tocursv))**2
        if (fu .le. fx) then
          if (u .ge. x) then
            a = x
          else
           b = x
        end if
        v = w
        fv = fw
        w = x
        fw = fx
        x = u
        fx = fu
      else
        if (u .lt. x) then
          a = u
        else
          b = u
        end if
        if (fu .le. fw .or. w .eq. x) then
          v = w
          fv = fw
          w = u
          fw = fu
        else if (fu .le. fv .or. v .eq. x .or. v .eq. w) then
          v = u
          fv = fu
        end if
      end if
   11 continue
    3 xmin = x
      tocur = xmin    ! best estimate of the current for min chisq
      return
c
      end

      subroutine fluxav12 (p, x, nx, y, ny, psival, npsi, fpsi,
     .                   ffppsival, qpsival,
     .                   bax, curaxis, xmagax, ymagax, ixcal, btor,
     .                   rmajor, ndisk, iscr, iounit, rmin, rmax,
     .                   zmin, zmax, elong, eps, vprime, rm2, psivolp,
     .                   xhm2, xi11,xi33, xips, rm2i, ravg,ravgi,
     .                   bsqinvavg,bsq_avg,b_avg,h_factr,ierr,
     .                   widep, hitep, kappa, rminor, rgeom, btgeom,
     .                   circum, rvloop, zvloop, psivloop, tocur,
     .                   ali, mhdmode, mhdmethd, psifctr, cxareanpsi,
     .                   grho1npsi, grho2npsi, rmajavnpsi, rminavnpsi,
     .                   triangnpsi,triangnpsi_l,pindentnpsi,sfareanpsi,
     .                   torfluxnpsi, ifixshap)
c
c
c ----------------------------------------------------------------------
c This subroutine calculates flux surface averages of quantities
c listed below.  The contours used are defined by the vector psival.
c note that this routine assumes psival(1) = value at plasma edge
c and psival(npsi) = value at magnetic axis.
c all input and output is strictly MKS!
c --- input through argument list:
c p(nx,ny)       psi (poloidal flux function)
c x(nx)
c y(ny)            psi is defined on the (x,y) mhdgrid
c psival           list of psi values for contours to be found
c npsi             number of elements in psival
c fpsi(npsi)       f(psi) defined over psival grid
c qpsi(npsi)       safety factor defined on psival grid
c bax              magnetic field (tesla) on magnetic axis
c curaxis          current density,amps/m**2 on magnetic axis
c xmagax
c ymagax           location of magnetic axis
c ndisk            Fortran unit number for file "scratch" or "scratch1"
c set ndisk = 0 if these files are not required
c iscr             = 1 implies create file "scratch1",
c                      write contour points in it and leave file open
c                      on return from FLUXAV12
c                  # 1 means same as above but file "scratch" is closed
c                      on return from FLUXAV12
c iounit           Fortran unit number for error messages
c                  (set iounit = 0 to suppress all error messages)
c btor             b-toroidal at R0
c rmajor           R0
c ixcal            switch for special calculations (see below)
c
c rvloop
c zvloop                  (r,z) of voltage loop,meters
c tocur                     total toroidal current,amps
c
c mhdmode              character variable, either mhdmode = 'coils'
c or mhdmode = 'no coils'
c mhdmethd
c
c ifixshap           used to set psifctr
c
c --- input through INCLUDE files:
c --- INCLUDE file mhdpar
c nw,nh           (required to define parameters in INCLUDE file bicube.i only)
c --- INCLUDE file constnts
c pi,pisq,u0              (u0 = 4 * pi * 1.0e-07)
c --- INCLUDE file spare
c aspare                   no idea on this one
c --- INCLUDE file storage
c xdum,ydum,zdum,wdum      work vectors (min size npsi)
c sdum,tdum,udum
c --- INCLUDE file bicube:
c cspln(n2cspln,nw,nh2) work array for bicubic spline coeff.
c wnoperm(nwork)        work array,min length =2*nw*nh+2*nh
c pds                   work vector (min length 6)
c
c --- INCLUDE file contour:
c rplasbdry
c zplasbdry
c nplasbdry              plasma boundary used, if mhdmode = 'no coils'
c rcontr
c zcontr
c ncontr                 plasma flux surface contours
c INCLUDE file zerocom:
c zero(nwnh)             used for diagnostic printout in dump_psi_values
c wzero(nwnh)            used in toroidal flux calculations
c
c INCLUDE file io.i:
c  use_efit_cntr         selects use of plasma boundary from eqdsk if=1
c --- output
c
c --- through argument list:
c rmin
c rmax
c zmin
c zmax               min/max extent of plasma
c widep                plasma width
c hitep                       height
c circum                      circumference
c rminor                      minor radius
c kappa                       vertical elongation (real number)
c rgeom                       geometric center
c btgeom                       btor at rgeom
c psivloop          twopi*(psi(rvloop,zvloop)-psi(plasma boundary))
c (used for inductive voltage calculations)
c ali                plasma inductance,defined as
c                    ali = integral(bp**2 dv)/(v*<bp**2>)
c                    this can be written as
c                    ali = vprime*integral(dpsi(integral bp*dl))/
c                    (volume*u0*tocur)
c                    However vprime at the plasma surface can be
c                    difficult to calculate accurately if an x point
c                    is present(due to the divergence of dl/bp).
c                    hence we use <bp**2> = (u0*I/circumfrence)**2 instead.
c cxareanpsi(j)  cross sectional area of plasma surfaces
c
c --- the following output vectors are all of length
c --- npsi and defined so that the first value corresponds to the
c --- plasma edge and the npsi'th value to the magnetic axis:
c
c elong(i)     elongation of flux surface i, = (zmax-zmin)/(rmax-rmin)
c eps(i)       horizontal inverse aspect ratio = (rmax-rmin)/(rmax+rmin)
c vprime(i)    dvolume/dpsi
c psivolp(i)   volume, psivolp(1) = total plasma volume
c rm2(i)              = <1/R**2>
c rm2i(i)             = <R**2>
c ravg(i)             = <R>
c ravgi(i)             = <1/R>
c bsqinvavg(i)        = < (Bt0/B)**2>
c bsq_avg(i)          = < B**2/Bt0**2>
c b_avg(i)            = <B/Bt0>
c h_factr             = < sqrt(1- B/Bmax).
c rmajavnpsi(i)       = avg rmajor (defined at magnetic axis elevation)
c rminavnpsi(i)       = avg rminor (              "                   )
c triangnpsi(i)       = (upper)triangularity
c triangnpsi_l(i)     = (lower) triangularity
c pindentnpsi(i)      = indentation
c sfareanpsi(i)       = surface area
c torfluxnpsi(i)
c grho1npsi(npsi)     = < ABS (grad psi)>
c grho2npsi(npsi)     = <(grad psi)**2>
c                     NOTE THAT GRHO1 AND GRHO2 ARE CHANGED INTO
c                     <abs (grad rho) > and < (grad rho)**2 > after
c                     we return from this suroutine
c --- the following flux surface average vectors are calculated only
c --- if ixcal = 1 . They are required for neoclassical transport models.
c --- Note that these are all dimensionless quantitites.
c xhm2(i)         <btotl**2/bax**2>
c xi11(i)           eq. 4.2-53   of ga-a16178
c xi33(i)           eq. 4.2-55
c xips(i)           eq. 4.2-59
c ierr     error flag. if no error, ierr = 0. if an error occurred
c ierr = 1 is returned and the above output can not be assumed
c available (or correct).
c --- there is no output to INCLUDE files !
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE contour
      USE io
      USE limiter
      USE mhdpar 
      USE tfact
      USE constnts
      USE ifs
      USE bicube
      USE metrics
      USE neo2dp
      USE cer
      USE numbrs,             ONLY : nj
      USE mhdcom,             ONLY : zero,wzero
      USE replace_imsl,       ONLY : my_ibcccu,my_dbcevl1 


      implicit  integer (i-n), real*8 (a-h, o-z)
D     include 'mpif.h'

      include 'shape.i'
      include 'spare.i'
      include 'storage.i'

      include 'imsl.i'

c      include 'zerocom.i'

c
      parameter (nlam = 150, nxxp = 5 )
      dimension  xxp(nxxp),yyp(nxxp)
      dimension  p(nx,ny),x(nx),y(ny),fpsi(npsi),psival(npsi)
      dimension  ffppsival(npsi), qpsival(npsi)
      dimension  elong(npsi),vprime(npsi),psivolp(npsi),rm2(npsi)
      dimension  rm2i(npsi),ravg(npsi),cxareanpsi(npsi)
      dimension  ravgi(npsi)
      dimension  bsqinvavg(npsi),bsq_avg(npsi),b_avg(npsi)
      dimension  h_factr(npsi)
      dimension  eps(npsi),xp(nconmax),yp(nconmax),arclen(nconmax)
      dimension  xhm2(*),xi11(*),xi33(*),xips(*)
      dimension  h(kstore),suml(nlam),bpinv(kstore)
      dimension  gp(nwork),psivlcpy(kstore),bpl(kstore)
      dimension  rmajavnpsi(npsi),rminavnpsi(npsi),triangnpsi(npsi),
     .           pindentnpsi(npsi),grho1npsi(npsi),grho2npsi(npsi),
     .           sfareanpsi(npsi),torfluxnpsi(npsi),triangnpsi_l(npsi)
      REAL *8    mag_bsq(nconmax),mag_b(nconmax)      
      REAL *8    max_b,factr
      dimension     ftlin(npsi), hlin(nconmax)

D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_elong, t_eps
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_psivolp
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_cxareanpsi
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_sfareanpsi, t_bpl,t_vprime
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_ravg,t_rm2i,t_rm2,t_ravgi
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_bsq_avg,t_bsqinvavg
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_grho1npsi,t_grho2npsi
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_grth,t_bsq,t_bmsq
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_grbmsq
D      real *8 ,DIMENSION(:,:),ALLOCATABLE:: t_gfm
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_cer_btdr,t_cer_bp
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_xhm2,t_xips,t_ftlin
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_ftnclp, t_xi11, t_xi33
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_rhod_psi_ifs
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_rminavnpsi
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_rmajavnpsi
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_triangnpsi
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_torfluxnpsi
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_pindentnpsi
      logical       change_cntour
      integer       mp_low, mp_high
      real*8        kappa, yp_dif, xp_mid
      character*(*) mhdmode, mhdmethd
      equivalence (psivlcpy(1),xdum(1))
      equivalence (bpinv(1),zdum(1))
      equivalence (h(1),wdum(1))
      equivalence (xp(1),sdum(1))
      equivalence (yp(1),tdum(1))
      equivalence (arclen(1),udum(1))
      equivalence (bpl(1),vdum(1))
c

c-------------------------------------------MPI related---------------
c --- for running with just one processor:
      myid=0        !must be set to zero if only 1 processor is used 
      master = 0    !must be 0 for do loop 50 to work properly 
      numprocs =1
c --- for running with multiple processors. The D in col 1
c     activates these lines if compiled with -Mdlines
c     to ignore these lines compile with -Mnodilines
D      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr) !get processor id
D      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr) !get num processors 
       !print *,"myid =",myid
       !print *,"numprocs = ",numprocs
      if(numprocs .gt. 1)then
D            !use f90 to zero the following arrays. These arrays are filled
D            !on each processer with a stride of numprocs. Hence we get
D            !filled arrays( with a stride of one)  by simply adding 
D            !together the non zero elemnts that each processor provides. 
D            allocate(t_elong(size(elong)))
D            allocate(t_eps(npsi))
D            allocate(t_psivolp(npsi)) 
D            allocate(t_cxareanpsi(npsi))
D            allocate(t_sfareanpsi(npsi))
D            allocate(t_bpl(npsi))
D            allocate(t_vprime(npsi))
D            allocate(t_rm2(npsi),t_rm2i(npsi),t_ravg(npsi))
D            allocate(t_ravgi(npsi))
D            allocate(t_bsq_avg(npsi),t_bsqinvavg(npsi))
D            allocate(t_b_avg(npsi))
D            allocate(t_h_fctr(npsi))
D            allocate(t_grho1npsi(npsi),t_grho2npsi(npsi))
D            allocate( t_grth(npsi), t_bsq(npsi))
D            allocate( t_bmsq(npsi),t_grbmsq(npsi))
D            allocate( t_gfm(3,npsi), t_cer_btdr(npsi))
D            allocate( t_cer_bp(npsi),t_ftnclp(npsi),t_xi11(npsi))
D            allocate( t_xi33(npsi), STAT = istat)
D            if(istat .ne. 0)
D    .          call allocate_error("t_xi33(npsi)",myid,istat)
D            allocate( t_torfluxnpsi(npsi), STAT = istat)
D            if(istat .ne. 0)
D    .          call allocate_error("t_torfluxnpsi(npsi)",myid,istat)
D            allocate( t_pindentnpsi(npsi), STAT = istat)
D            if(istat .ne. 0)
D    .          call allocate_error("t_pindentnpsi(npsi)",myid,istat)
D            allocate( t_rhod_psi_ifs(npsi))
D            allocate( t_rmajavnpsi(npsi),STAT = istat)
D            if(istat .ne. 0)
D    .          call allocate_error("t_rmajavnpsi(npsi)",myid,istat)
D            allocate (t_rminavnpsi(npsi))
D            if(istat .ne. 0)
D    .          call allocate_error("t_rminavnpsi(npsi)",myid,istat)
D            allocate( t_triangnpsi(npsi) )
D            t_elong(:) = 0.0 ; t_eps(:) = 0.0 
D            t_psivolp(:) = 0.0 ;  t_cxareanpsi(:) = 0.0 
D            t_sfareanpsi(:) =0.0; t_rhod_psi_ifs(:) = 0.0
D            t_bpl(:) = 0.0 ;  t_vprime(:) =0.0
D            t_rm2(:) =0.0 ; t_rm2i(:) =0.0
D            t_ravg(:) =0.0 ; t_bsq_avg(:) =0.0 ;t_ravgi(:) =0.0
D            t_b_avg(:) =0.0
D            t_h_factr(:) =0.0
D            t_bsqinvavg(:) = 0.0 ; t_grho1npsi(:) =0.0
D            t_grho2npsi(:) = 0.0
D            t_grth(:) =0.0
D            t_bsq(:) = 0.0 ; t_bmsq(:) =0.0 ; t_grbmsq(:) =0.0
D            t_gfm(:,:) =0.0 ;  t_cer_btdr(:) = 0.0 ; t_cer_bp(:)=0.0
D            t_triangnpsi(:) = 0.0
D            t_pindentnpsi(:) = 0.0
D            t_rmajavnpsi(:) = 0.0 
D            t_rminavnpsi(:) = 0.0
D            if(ixcal .ne. 0)then
D               allocate(t_xhm2(npsi),t_xips(npsi),t_ftlin(npsi))
D               t_xhm2(:) =0.0 ; t_xips(:) =0.0 ; t_ftlin(:) = 0.0
D               t_ftnclp(:) = 0.0 ; t_xi11(:) =0.0 ; t_xi33(:) =0.0
D            endif
      endif


c----------------------------------------end MPI related-------------------







      change_cntour = .false.
      one    = 1.0
      isetw  =  0          ! calc wzero on first call to torflux2
      mpmin  = 50   ! # points on plasma surface is < this => error exit
      imslmd = 'fluxav'
      ierr   = -1
      if (nx .ne. nw .or. ny .ne. nh)  go to 1000    ! error exit
      ierr   =  0
****  nlam   = 15   ! param sets grid size for xi11 and xi33 calculation
      isignn =  1
      if (bax .lt. 0.0)  isignn = -1
c
c if  iscr .ne. 1  create file 'scratch'  and write flux surface points into it,
c then close 'scratch'
c if  iscr .eq. 1  create file 'scratch1' and write flux surface points into it,
c then leave file open for subsequent addition by RHOSET
c
      ndisko = ndisk
      if (ndisk .gt. 0 .and. myid .eq. master ) then
        call getioun(ndisko,ndisk)
        if( iscr .ne. 1) then
          open (unit = ndisko, file = 'scratch' , status = 'UNKNOWN')
        else
          open (unit = ndisko, file = 'scratch1', status = 'UNKNOWN')
        endif
      endif
      ndisk = ndisko
c
c ----------------------------------------------------------------------
c initialize for contour subroutines
c need max psi at magnetic axis for CNTOUR
c ----------------------------------------------------------------------
c
      call copya (psival, psivlcpy, npsi)    ! copy psival into psivlcpy
c
c --- psifctr is used to pull in the outer boundary for stability
c
      psivlcpy(1) = psifctr * (psivlcpy(2) - psivlcpy(1)) + psivlcpy(1)
      cconst = -1.0
      call multpl1 (p       , nx*ny, cconst)
      call multpl1 (psivlcpy,  npsi, cconst)

c
c --- set up bicubic spline coefficient array cspln
c
      imslmd='997 c208'
      call my_ibcccu (p, x, nx, y, ny, cspln, nx, wnoperm, ierr)



c
****  call multpl1 (p,nx*ny,cconst) ! p is changed back after we're..
c                                   ..done with torflux2
c
c --- all contours to be found must be inside the box defined by
c --- xmin, xmax, ymin, ymax
c
   10 xmin  = xlimiter(nlimiter+1)
      xmax  = xlimiter(nlimiter+2)
      ymin  = ylimiter(nlimiter+1)
      ymax  = ylimiter(nlimiter+2)
      dx    =  0.0
      dy    =  0.0
      arcl  =  0.02
      arcl  =  0.001 ! 999,
      taxis =  5.0
      tlim  = 30.0
      a     = (tlim - taxis) / (psivlcpy(1) - psivlcpy(npsi))
      bincp =         taxis
c
c ----------------------------------------------------------------------
c for npsi-1 values of psi in psivlcpy do in sequence
c from the plasma edge to the contour next to the magnetic axis:
c a) find the contour (xp(i),yp(i)),i = 1...mp
c b) form the required averages
c for the magnetic axis,j = npsi,results are obtained after were done here
c ----------------------------------------------------------------------
c
      delta_psi=psivlcpy(2)-psivlcpy(1)
c      do 50 j=1,npsi-1               ! npsi is magnetic axis
      do 50 j= myid+1,npsi-1,numprocs !parallel version also works for non mpi
         isetarcl = 0
 55      iauto    = 0
         if (j .eq. 1)  iauto  = 1   ! j = 1 is plasma edge
         iconvg   = 0
         if (j .eq. 1)  iconvg = 1
         dang     =    a * (psivlcpy(j) - psivlcpy(npsi)) + bincp
         dang     = dang * (isetarcl + 1)
         dang     = 0.1  ! 999,
         bperr    = 0.01
         bperr    = 1.e-3 ! 999,
c
         if (j .eq. 1) then
            if (use_efit_cntr .eq. 1) then
              write (iounit, '("using eqdsk contour for outermost ",
     .                         "flux surface ")')
              go to 15
            else
              write  (iounit, 14)  'psivlcpy(   1) =', psivlcpy(1)
              write  (iounit, 14)  'psivlcpy(npsi) =', psivlcpy(npsi)
              write  (iounit, 14)  '        xmagax =', xmagax
              write  (iounit, 14)  '        ymagax =', ymagax
   14         format ('  in FLUXAV: ', a, f18.12)
            end if
         end if
c
c         if (j .eq. 1 .and. mhdmode  .eq. 'no coils'
c     .                .and. mhdmethd .eq. 'sorpicrd')  go to 15
c


         call cntour (xmagax,ymagax,psivlcpy(j),rcmin,rcmax,zcmin,zcmax,
     .                zrcmin,zrcmax,rzcmin,rzcmax,dang,arcl,bperr,
     .                dx,dy,xmin,xmax,ymin,ymax,iauto,iautoc,xp,yp,
     .                mp,x,nx,y,ny,cspln,n2cspln,nh2,iounit,nconmax,
     .                ierr,gp,iconvg,delta_psi)


         if (ierr .eq. 0  .and.  j .eq. 1)
     .   write (iounit, *) ' number of points on plasma boundary = ', mp
c
c --- too many points were found; reduce the accuracy requirement and try again:
c
         if (ierr .gt. 1) then
           isetarcl = isetarcl + 1
               arcl =     arcl * 1.3
           if (isetarcl .lt. 5)  go to 55
         end if
c
         if (rcmin .eq. rcmax)  ierr = 1
         if ( ierr .ne. 0    )  go to 1000    ! error exit
         if (j .eq. 1 .and. mhdmode .ne. 'no coils') then
           call copya (xp, rplasbdry, mp)
           call copya (yp, zplasbdry, mp)
           nplasbdry = mp
         end if
         go to 16
c
c --- mhdmode = 'no coils', fixed boundary case. don't trace the contour
c --- instead use the known contour:
c
   15    call fixedcntour (rplasbdry,zplasbdry,nplasbdry,
     .                     xp,yp,mp,rcmin,rcmax,zcmin,zcmax,
     .                     rzcmin,rzcmax,zrcmin,zrcmax,x,y,nx,ny,
     .                     gp,cspln,n2cspln,nh2,pds)
            call limiter_check(rcmin,rcmax,zcmin,zcmax,
     .                    xlimiter,ylimiter,nlimiter) 
   16    arclen(1) = 0.0
         if (j .eq. 1 .and. mp .lt. mpmin)  go to 1000
         do i=1,mp-1
           arclen(i+1) = arclen(i) + SQRT ((xp(i+1)-xp(i))**2
     .                             + (yp(i+1)-yp(i))**2)
         end do
c
c --- elongation and (horizontal) inverse aspect ratio
c
         elong(j) = (zcmax-zcmin)/(rcmax-rcmin)
         eps  (j) = (rcmax-rcmin)/(rcmax+rcmin)
         rhod_psi_ifs(j) = 0.5*(rcmax-rcmin)   ! rhod_ifs on npsi grid
         if(numprocs .gt. 1)then
D             t_rhod_psi_ifs(j) =  rhod_psi_ifs(j)
D             t_elong(j)=elong(j)
D             t_eps(j)=eps(j)
         endif
         if (j .gt. 1)  go to 110
c
c --- save the following parameters (on the plasma edge,j = 1) for output
c --- for more than one processor this info is stored only on the master
c --- and must be broadcast to the other processors.
         rmin      = rcmin
         rmax      = rcmax
         zmin      = zcmin
         zmax      = zcmax
         widep     = rmax-rmin
         hitep     = zmax-zmin
         kappa     = hitep/widep
         rminor    = 0.5 * widep
         rgeom     = 0.5 * (rmax+rmin)
         btgeom    = btor*rmajor/rgeom
         circum    = arclen(mp)
         bpminimum = gp(1)
         do 105 nc=1,mp
           bpminimum = MIN (bpminimum,gp(nc))
           if (bpminimum .eq. gp(nc)) then
             rbpminimum = xp(nc)
             zbpminimum = yp(nc)
           end if
           rcontr(nc) = xp(nc)
  105      zcontr(nc) = yp(nc)
c
         ncontr = mp
c
c --- calculate the plasma volume ( = twopi*integral (R*Z*dR)).
c
         call volcalc (rcontr, zcontr, ncontr,
     .                 xmagax, ymagax, psivolp(j), cxareanpsi(j))
         if(numprocs .gt. 1)then
D             t_psivolp(j)=psivolp(j)
D             t_cxareanpsi(j) =cxareanpsi(j)
         endif
c
c ----------------------------------------------------------------------
c calculate volume of flux surface j
c ----------------------------------------------------------------------
c
  110    call volcalc (xp, yp, mp, xmagax, ymagax,
     .                 psivolp(j), cxareanpsi(j))

c
c ----------------------------------------------------------------------
c problems in contour tracing can be detected here
c ----------------------------------------------------------------------
c
      if (j .gt. 1 .and. numprocs .eq. 1 ) then
         if (psivolp(j)-psivolp(j-1) .ge. 0.0) then
            if (.not. change_cntour) then
               change_cntour = .not. change_cntour
               use_cnt1      = .not. use_cnt1
               use_cnt2      = .not. use_cnt2
c
****           if file scratch or scratch1 is open, rewind to start over
****
****           inquire (unit=ndisk, iostat=ios, err=11, opened=lopen)
****           if (lopen)  rewind ndisk
c
               if (ndisk .gt. 0)  rewind ndisk
               call copya (psival, psivlcpy, npsi)
               psivlcpy(1) = psifctr * (psivlcpy(2) - psivlcpy(1))
     .                                               + psivlcpy(1)
               cconst = -1.0
               call multpl1 (psivlcpy,  npsi, cconst)
               go to 10
   11          call STOP ('subroutine FLUXAV: scratch file prob', 204)
            else
               write (*,'(" subroutine FLUXAV reports:"                /
     .                    " error in determination of plasma contours" /
     .                    " Using input variable PSIFCTR may help"     /
     .                    " See file cray102.f for description")')
               if (nw .lt. 65)
     .         write (*,'(" using a finer mesh eqdsk"                  /
     .                    " with the APPROPRIATE VERSION of ONETWO"    /
     .                    " may also eliminate the problem")')
               call STOP ('subroutine FLUXAV: contour problem', 79)
            end if
         end if
       end if
c
c ----------------------------------------------------------------------
c calculate surface area of flux surface j
c ----------------------------------------------------------------------
c
         call surfarea (xp, yp, mp, sfareanpsi(j))

c
c ----------------------------------------------------------------------
c calculate the toroidal flux inside the surface
c as used here wzero has to be redefined for each contour
c so isetw is never set to 1. These calculations are based
c on the cell size of the MHD grid and become progressively more
c inaccurate as the cell size increases or (the contour shrinks
c around the magnetic axis). we use it only for the outermost surface
c even then, the 33 by 65 grid may be a bit too coarse.
c ----------------------------------------------------------------------
c
         if (j .eq. 1) then
           call torflux2 (p, nx, ny, x, y, wzero, xp, yp, mp,
     .                    torfluxnpsi(j), fpsi, psivlcpy, npsi, isetw)
         else
           torfluxnpsi(j) = 0.0  ! remaining values calculated in RHOSET
         end if

c
c ----------------------------------------------------------------------
c calculate line integral of mag(bpoloidal)
c ----------------------------------------------------------------------
c
         call integ (arclen, gp, mp, bpl(j))

c
c ----------------------------------------------------------------------
c calculate vprime and line integral of 1.0 / bpoloidal
c ----------------------------------------------------------------------
c
            do i=1,mp
               if (gp(i) .eq. 0.0) then
                 write  (6, 3210)  i, mp
 3210            format (' ERROR: gp(i) = 0, i, mp =', 2(2x,i5))
                 call STOP ('subroutine FLUXAV: problem #1', 150)
               end if
               bpinv(i) = 1.0 / gp(i)
            end do
c
            call integ (arclen,bpinv,mp,bpinteg)
            vprime(j) = 2.0 * pi * bpinteg
c            if (j .eq. 1)  bpintsv = bpinteg ! save for use below
            if (bpinteg .eq. 0.0) then
              write  (6, 3212) mp
 3212         format (' ERROR: bpinteg = 0 ,mp =',i5)
              call STOP ('subroutine FLUXAV: problem #2', 151)
            end if

c
c ----------------------------------------------------------------------
c calculate rm2 = <1/R**2>,rm2i=<R**2>,ravg=<R>,ravgi=<1/R>
c bsqinvavg=<Bt0**2/B**2>
c ----------------------------------------------------------------------
c
            do i=1,mp
              if (xp(i) .eq. 0.0) then
                write  (6, 3211)  i, mp
 3211           format (' error,xp(i) = 0, i, mp =', 2(2x,i5))
                call STOP ('subroutine FLUXAV: problem #3', 152)
              end if
              ydum(i) = bpinv(i)/(xp(i)**2)
            end do
c
            call integ (arclen, ydum, mp, sum2)
            rm2(j) = sum2/bpinteg
            do 140 i=1,mp
  140       ydum(i) = bpinv(i)*xp(i)**2
            call integ (arclen,ydum,mp,sum2)
            rm2i(j) = sum2/bpinteg

            do 150 i=1,mp
  150       ydum(i) = bpinv(i)*xp(i)
            call integ (arclen,ydum,mp,sum2)
            ravg(j) = sum2/bpinteg

            do  i=1,mp
               ydum(i) = bpinv(i)/xp(i)
            enddo
            call integ (arclen,ydum,mp,sum2)
            ravgi(j) = sum2/bpinteg
c
c
c           bsq_avg,mag_bs,mag_b,max_b:
c
            max_b = 0.0D0
            do i=1,mp
              mag_bsq(i) =    (fpsi(j)/xp(i))**2 + gp(i)**2
              mag_b(i)   =    SQRT(mag_bsq(i))
              max_b      =    MAX(max_b,mag_b(i))
              ydum(i)    =    mag_bsq(i)*bpinv(i)
            end do
            call integ (arclen,ydum,mp,sum2)
            bsq_avg(j) = sum2/bpinteg/Btor**2
c
c
c           b_avg :
c
            do i=1,mp
              ydum(i) = mag_b(i)*bpinv(i)
            end do
            call integ (arclen,ydum,mp,sum2)
            b_avg(j) = sum2/bpinteg/ABS(Btor)

c
c           h_factr = <sqrt(1-B/Bmax)>:
c
            do i =1,mp
               factr  = 1.d0-mag_b(i)/max_b
               ydum(i) =0.0D0
               if(factr .gt. 0.0D0)
     .           ydum(i) = SQRT(factr)*bpinv(i)
            enddo
            call integ (arclen,ydum,mp,sum2)
            h_factr(j)  = sum2/bpinteg/ABS(Btor)

c
c           bsqinvavg is the same as btor*btor*bmsq
c           it is recalculated here because bmsq is special
c           to the Houlberg model (and may not be avaialble)  HSJ 8/22/98
c
            do i=1,mp
              ydum(i) = bpinv(i) / ((fpsi(j)/xp(i))**2 + gp(i)**2)
            end do
c
            call integ (arclen,ydum,mp,sum2)
            bsqinvavg(j) = btor*btor*sum2/bpinteg

c
            call integ(arclen,xp,mp,sum2)
            grho1npsi(j) = sum2/bpinteg             ! <ABS (grad psi)>

            do 156 i=1,mp
  156       ydum(i) = xp(i)*xp(i)*gp(i)
            call integ(arclen,ydum,mp,sum2)
            grho2npsi(j) = sum2/bpinteg             ! <(grad psi)**2>

c
c find the index of the inboard and outboard points closest
c to the magnetic axis elevation, ir and il:
c
            almin = y(nh)
            armin = y(nh)
            do i=1,mp-1
              advert = ABS (ymagax-yp(i))
              if (xp(i) .le. xmagax) then
                almin = MIN (almin, advert)
                if (almin .eq. advert)  il = i
              else
                armin = MIN (armin, advert)
                if (armin .eq. advert)  ir = i
              end if
            end do
            if (ir*il .eq. 0  .or.  ir .eq. il)  go to 1000
            rmajavnpsi (j) =  0.5 * (xp(ir) + xp(il))
            rminavnpsi (j) =  0.5 * (xp(ir) - xp(il))
            triangnpsi (j) = (0.5 * (rcmax+rcmin)-rzcmax)/rminavnpsi(j)
          triangnpsi_l (j) = (0.5 * (rcmax+rcmin)-rzcmin)/rminavnpsi(j)
            pindentnpsi(j) = (xp(il) - rcmin) / (2.0 * rminavnpsi(j))
            if(numprocs .gt. 1)then
D             t_vprime(j) =vprime(j)
D             t_bpl(j)=bpl(j)
D             t_torfluxnpsi(j) = torfluxnpsi(j)
D             t_psivolp(j) = psivolp(j)
D             t_cxareanpsi(j) =cxareanpsi(j)
D             t_sfareanpsi(j)= sfareanpsi(j)
D             t_pindentnpsi(j) = pindentnpsi(j)
D             t_triangnpsi(j) = triangnpsi(j)
D             t_rmajavnpsi(j) = rmajavnpsi(j)
D             t_rminavnpsi(j) = rminavnpsi(j)
D             t_grho1npsi(j) = grho1npsi(j)
D             t_grho2npsi(j) =  grho2npsi(j)
D             t_bsqinvavg(j) = bsqinvavg(j)
D             t_rm2(j) = rm2(j)
D             t_rm2i(j) = rm2i(j)
D             t_ravg(j) = ravg(j)
D             t_ravgi(j) = ravgi(j)
D             t_bsq_avg(j) =  bsq_avg(j)
D             t_b_avg(j)   = b_avg(j)
D             t_h_factr(j) = h_factr(j)
            endif
c
c ----------------------------------------------------------------------
c Metrics for Houlberg/NCLASS bootstrap models
c ----------------------------------------------------------------------
c 
            call getbsmodel (jhirsh0)
            if (jhirsh0 .eq. 95 .or. jhirsh0 .eq. 96
     .         .or. jhirsh0 .eq. 99 .or. jhirsh0 .eq. 100) then
              if (j .eq. 1 .and. jhirsh0 .lt. 99)  write (iounit, *)
     .       ' calculating bootstrap current using Houlberg Model'
            if (j .eq. 1 .and. (jhirsh0 .eq. 99 .or. jhirsh0 .eq. 100))
     .         write (iounit, *)
     .       ' calculating bootstrap current using NCLASS Model'
              call nclmetrics (xmagax, ymagax, rmajor, mp, xp, yp, gp,
     .                         eps(j),fpsi(j),qpsival(j), ffppsival(j),
     .                         x, nx, y, ny, cspln, grth(j), bsq(j),
     .                         bmsq(j), grbmsq(j), gfm(1,j))
              if(numprocs .gt. 1)then
D                t_grth(j)   = grth(j)
D                t_bsq(j)    = bsq(j)
D                t_bmsq(j)   = bmsq(j)
D                t_grbmsq(j) = grbmsq(j)
D                t_gfm(:,j)  = gfm(:,j)
              endif
              if (j_cer .ne. 0) then
c
c                 compute Bp and Bt/R along CER cord at Z=0
c
                  mp_low = 0
                  xp_mid = xmax
                  cer_btdr(j) = 0.0
                  cer_bp(j)   = 0.0
                  if(numprocs .gt. 1)then
D                    t_cer_bp(j) = cer_bp(j)
D                    t_cer_btdr(j) = cer_btdr(j)
                  endif
c
c                 find the Z=0 crossings
c
                  do i=1,mp-1
                     if (yp(i)*yp(i+1) .le. 0) then
                        if (mp_low .eq. 0) then
                          mp_low = i
                          xp_mid = xp(i)
                        end if
                        if (mp_low .ne. 0 .and. xp(i) .gt. xp_mid)
     .                      mp_low = i
                     end if
                  end do
                  if (yp(mp_low+1) .eq. 0.0)  mp_low = mp_low+1
                  if (mp_low .ne. 0) then
c
c                   a crossing was found proceed
c
                    mp_high = mp_low
                    if (yp(mp_low)*yp(mp_low+1).lt.0) mp_high = mp_low+1
                    if (yp(mp_low)*yp(mp_low-1).lt.0) mp_high = mp_low-1
                    yp_dif = yp(mp_high) - yp(mp_low)
                    if (yp_dif .ne. 0.0)  yp_dif = -yp(mp_low)/yp_dif
                    xp_mid = xp(mp_low) +
     .                      (xp(mp_high) - xp(mp_low))*yp_dif
                    cer_btdr(j) = fpsi(j)/xp_mid**2
                    cer_bp(j) = gp(mp_low) +
     .                         (gp(mp_high) - gp(mp_low))*yp_dif
                    if(numprocs .gt. 1)then
D                      t_cer_btdr(j) = cer_btdr(j)
D                      t_cer_bp(j) = cer_bp(j)
                    endif
                  end if
              end if
            end if
c
c ----------------------------------------------------------------------
c special calculations
c ----------------------------------------------------------------------
c
             if (ixcal .eq. 0)  go to 160

c
c --- get max total b field on contour
c
             bmax = -1.0e30
c
             do i=1,mp
               gp(i) = (xp(i)*gp(i))**2 ! convert to grad psi squared
               bb = SQRT (fpsi(j)**2+gp(i))/xp(i)
               if (bb .gt. bmax) bmax = bb
             end do
c
             hmin =     bax / (isignn * bmax)
             dl2  = hmin**2 / (nlam - 1)
c
             do i=1,mp
                h(i)    = (isignn * SQRT (fpsi(j)**2+gp(i))/xp(i))
                if (h(i) .eq. 0.0) then
                   write  (6, 3215)  i, mp
 3215              format (' ERROR: h(i)=0, i, mp =', 2(2x,i5))
                   call STOP ('subroutine FLUXAV: problem #4', 153)
                 end if
               h(i)    = bax / (isignn * SQRT (fpsi(j)**2+gp(i))/xp(i))
               ydum(i) = h(i)**2 * bpinv(i)
               hlin(i) = ABS (bax/(bmax*h(i)))
               hlin(i) = MIN (hlin(i), one)   ! required due to roundoff
             end do
c
             call integ (arclen, ydum, mp, sum2)
             xhp2 = sum2 / bpinteg
             if (j .eq. 1)  xhp2sv = xhp2 ! save for use below
c
             do i=1,mp
               ydum(i) = bpinv(i) / h(i)**2
             end do
c
             call integ (arclen, ydum, mp, sum2)
             xhm2(j) = sum2 / bpinteg

             if (xhm2(j) .eq. 0.0) then
               write  (6, 3217)  j, mp
 3217          format (' ERROR: xhm2(j) = 0   j, mp =', 2(2x,i5))
               call STOP ('subroutine FLUXAV: problem #5', 154)
             end if
c
             xips(j) = xhp2 - 1.0 / xhm2(j)

c

             do l=1,nlam
               al = SQRT ((l-1)*dl2)
               do i=1,mp
                 val     =  ABS (1.0-al/h(i))
                 ydum(i) = SQRT (val) * bpinv(i)
               end do
               call integ (arclen, ydum, mp, suml(l))
               suml(l) = suml(l) / bpinteg
             end do
c
            xi0 = 0.0
            do l=1,nlam-1
              if (suml(l) .eq. 0.0 .or. suml(l+1) .eq. 0.0) then
                write  (6, 3218)  l, mp
 3218           format (' ERROR: sumll(l),sumll(l+1) = 0.0,l,mp =',
     .                    2(2x,i5))
                call STOP ('subroutine FLUXAV: problem #6', 155)
              end if
              xi0 = xi0 + 0.5 * (1.0/suml(l)+1.0/suml(l+1))*dl2 * 0.5
            end do
c
c           use approximation of Liu, GAA-21820 for ft
c
            do i=1,mp
              ydum(i) = (hlin(i)**2) * bpinv(i)
            end do
            call integ (arclen, ydum, mp, hlin2avg)
            hlin2avg = hlin2avg/bpinteg           ! <h**2>
            do i=1,mp
              ydum(i) = hlin(i) * bpinv(i)
            end do
            call integ (arclen, ydum, mp, hlinavg)
            hlinavg = hlinavg/bpinteg             ! <h>
            do i=1,mp
              ydum(i) = (1.0 - (SQRT (1.0-hlin(i))*(1.0+0.5*hlin(i))))
     .                                 / (hlin(i)**2)
              ydum(i) = ydum(i)*bpinv(i)
            end do
            call integ (arclen, ydum, mp, hcomavg)
            hcomavg   = hcomavg / bpinteg
            ftl       = 1.0 - hlin2avg * hcomavg
            ftu       = 1.0 + 0.5 * hlinavg
            ftu       = 1.0 - ftu * SQRT (1.0 - hlinavg)
            ftu       = ftu * hlin2avg / (hlinavg * hlinavg)
            ftu       = 1.0 - ftu
            ftlin (j) = 0.75 * ftu + 0.25 * ftl

            if (ftcalc .eq. 'analytic')  ftnclp(j) = ftlin(j)

            if (j .eq. 1)  xi0sv = xi0    ! save for use below
            xi11(j)   =  1.33333*xhp2-xi0
            xi33(j)   = (1.33333-xhm2(j)*xi0)*xhm2(j)
            if(numprocs .gt. 1)then
D               t_xi11(j) =  xi11(j)
D               t_xi33(j) = xi33(j)
D               t_ftlin(j) = ftlin(j)
D               t_ftnclp(j) = ftlin(j)
D               t_xhm2(j) = xhm2(j) 
D               t_xips(j) = xips(j)
            endif


  160       if (ndisk .le. 0 .or. myid .ne. master)  go to 50
c
c ----------------------------------------------------------------------
c write points to file ndisk in CGS (Gaussian) units
c ----------------------------------------------------------------------
c

               write (ndisk, 8100)  mp
               write (ndisk, 8120)  (100.0*xp(i), 100.0*yp(i), i=1,mp)
               write (ndisk, 8120)  (1.0e-04*bpinv(i),i=1,mp)
               write (ndisk, 8120)  (100.0*arclen(i),i=1,mp)
 8100          format (i6,2x,i6)
 8120          format (6e12.5)
c
   50        continue         ! end loop on psi values to be contoured


D         call MPI_Barrier( MPI_COMM_WORLD,ierr)



        if(numprocs .gt. 1)then !merge results
C            note that the receiving array,elong, is zeroed automatically before the reduction
D            call MPI_Allreduce(t_elong,elong,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_elong)


D            call MPI_Allreduce(t_eps,eps,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_eps)

D            call MPI_Allreduce(t_psivolp,psivolp,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_psivolp)

D            call MPI_Allreduce(t_cxareanpsi,cxareanpsi,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_cxareanpsi)

D            call MPI_Allreduce(t_bpl,bpl,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_bpl)

D            call MPI_Allreduce(t_vprime,vprime,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_vprime)

D            call MPI_Allreduce(t_rm2,rm2,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_rm2)

D            call MPI_Allreduce(t_rm2i,rm2i,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_rm2i)

D            call MPI_Allreduce(t_ravg,ravg,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_ravg)

D            call MPI_Allreduce(t_ravgi,ravgi,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_ravgi)

D            call MPI_Allreduce(t_bsq_avg,bsq_avg,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_bsq_avg)

D            call MPI_Allreduce(t_b_avg,b_avg,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_b_avg)

D            call MPI_Allreduce(t_h_factr,b_avg,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_h_factr)

D            call MPI_Allreduce(t_bsqinvavg,bsqinvavg,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_bsqinvavg)

D            call MPI_Allreduce(t_grho1npsi,grho1npsi,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_grho1npsi)

D            call MPI_Allreduce(t_bsq,bsq,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_bsq)

D            call MPI_Allreduce(t_bmsq,bmsq,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_bmsq)

D            call MPI_Allreduce(t_grbmsq,grbmsq,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_grbmsq)

D            call MPI_Allreduce(t_gfm,gfm,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_gfm)

D            call MPI_Allreduce(t_cer_btdr,cer_btdr,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_cer_btdr)

D            call MPI_Allreduce(t_cer_bp,cer_bp,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_cer_bp)


D            call MPI_Allreduce(t_sfareanpsi,sfareanpsi,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D            deallocate(t_sfareanpsi)
D               if(istat .ne. 0) call deallocate_error (
D    &                               "t_sfareanpsi",myid,myid,istat)

D            if(ixcal .ne. 0)then
D               call MPI_Allreduce(t_xhm2,xhm2,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D               deallocate(t_xhm2)

D               call MPI_Allreduce(t_xips,xips,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D               deallocate(t_xips)

D               call MPI_Allreduce(t_ftlin,ftlin,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D               deallocate(t_ftlin)

D               if(ftcalc .eq. 'analytic')
D    .               call MPI_Allreduce(t_ftnclp,ftnclp,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D               deallocate(t_ftnclp)

D               call MPI_Allreduce(t_xi11,xi11,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D               deallocate(t_xi11)

D               call MPI_Allreduce(t_xi33,xi33,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D               deallocate(t_xi33,STAT = idealc)
D               if(idealc .ne. 0)print *,'idealcxi33  =',idealc

D               call MPI_Allreduce(t_rmajavnpsi,rmajavnpsi,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D               deallocate(t_rmajavnpsi,STAT = istat)
D               if(istat .ne. 0) call deallocate_error (
D    &                               "t_rmajavnpsi",myid,istat)



D               call MPI_Allreduce(t_rminavnpsi,rminavnpsi,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D               deallocate(t_rminavnpsi,STAT = istat)
D                if(istat .ne. 0) call deallocate_error (
D    &                               "t_rmajavnpsi",myid,istat)

D               call MPI_Allreduce(t_triangnpsi,triangnpsi,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D               deallocate(t_triangnpsi)
D                if(istat .ne. 0) call deallocate_error (
D    &                               "t_rmajavnpsi",myid,istat)


D               call MPI_Allreduce(t_pindentnpsi,pindentnpsi,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D                deallocate(t_pindentnpsi, STAT = istat)
D                if(istat .ne. 0) call deallocate_error (
D    &                               "t_pindentnpsi",myid,istat)

D               call MPI_Allreduce(t_torfluxnpsi,torfluxnpsi,npsi,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
D               deallocate(t_torfluxnpsi, STAT = istat)
D                if(istat .ne. 0) call deallocate_error (
D    &                               "t-_torfluxpsi",myid,istat)
D            endif

c         broadcast quantities that were only calculated for j =1
c         (which is done by process 0)


D           if(nplasbdry .ne. 0 )then
D             call MPI_BCAST(nplasbdry,1,MPI_Integer,
D    .                      master,MPI_COMM_WORLD,ierr)
c             nplasbdry must be BCAST first so following will work
c             The unwritten rule seems to be that the number of 
c             of elements broadcast must exist on the listener
c             apriori.
D             call MPI_BCAST(rplasbdry,nplasbdry,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D             call MPI_BCAST(zplasbdry,nplasbdry,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)


D           endif
D           call MPI_Barrier( MPI_COMM_WORLD,ierr)
D            if(myid .eq. 0)
D    .       print *,"nplasbdry, id 0 =",nplasbdry
D           if(myid .eq. 1)
D    .       print *,"nplasbdry, id 1 =",nplasbdry
D          call MPI_BCAST(rmin,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)

D          call MPI_BCAST(rmax,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D          call MPI_BCAST(zmin,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)



D           call MPI_Barrier( MPI_COMM_WORLD,ierr)
D            if(myid .eq. 0)
D    .       print *,"zmin id 0 =",zmin
D           if(myid .eq. 1)
D    .       print *,"zmin, id 1 =",zmin







D          call MPI_BCAST(zmax,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D          call MPI_BCAST(widep,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D          call MPI_BCAST(hitep,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D          call MPI_BCAST(kappa,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D          call MPI_BCAST(rminor,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D          call MPI_BCAST(rgeom,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)



D           call MPI_Barrier( MPI_COMM_WORLD,ierr)
D            if(myid .eq. 0)
D    .       print *,"rgeom id 0 =",rgeom
D           if(myid .eq. 1)
D    .       print *,"rgeom, id 1 =",rgeom



D          call MPI_BCAST(btgeom,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D          call MPI_BCAST(circum,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D          call MPI_BCAST(bpminimum,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D          call MPI_BCAST(rpminimum,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D          call MPI_BCAST(zbpminimum,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)




D           call MPI_Barrier( MPI_COMM_WORLD,ierr)



D             call MPI_BCAST(ncontr,1,MPI_integer,
D    .                      master,MPI_COMM_WORLD,ierr)
D             call MPI_BCAST(rcontr,ncontr,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D             call MPI_BCAST(zcontr,ncontr,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)


	 endif   !numprocs .gt. 1


c
        if (use_efit_cntr .eq. 1) then
          vprime(1) = (psivolp(1)-psivolp(2))/delta_psi
          rm2   (1) = ABS (twopi*twopi*qpsival(1) / (fpsi(1)*vprime(1)))
        end if
c
c --- calculate elongation on axis based on general quadratic form of
c --- ellipse   a*x**2+2b*x*y+c*y**2+d*x+e*y+f   the coefficients a,b,
c --- c,d,e,f are obtained from the Taylor Series expansion of psi
c --- about the magnetic axis
c
       call my_dbcevl1 (x,nx,y,ny,cspln,nx,xmagax,ymagax,pds,ierr,6)
       a = 0.5 * pds(5)
       b = 0.5 * pds(4)
       c = 0.5 * pds(6) 
c
c --- rotate coordinates
c
       thet    = 0.5 * ATAN2 (2*b, a-c)
       ap      = a * COS (thet)**2  +  2.0 * b * SIN (thet) * COS (thet)
     .         + c * SIN (thet)**2
       cp      = a * SIN (thet)**2  -  2.0 * b * SIN (thet) * COS (thet)
     .         + c * COS (thet)**2
       discrim = b**2-a*c ! < 0 ellipse, = 0 parabola, > 0 hyperbola
c
       if (ap/cp .gt. 0.0 .and. discrim .lt. 0.0) then
         elong(npsi) = SQRT (cp/ap)
       else
c
c      elongation on axis is not used for other calculations,
c      so the following crude approximation is OK
c
         elong(npsi) = elong(npsi-1)
       end if
c

       rm2        (npsi) = 1.0 / xmagax**2
       rm2i       (npsi) = 1.0 / rm2(npsi)
       ravg       (npsi) = xmagax
       ravgi      (npsi) = 1./xmagax
       psivolp    (npsi) = 0.0
       bpl        (npsi) = 0.0
       eps        (npsi) = 0.0
       rmajavnpsi (npsi) = xmagax
       rminavnpsi (npsi) = 0.0
       triangnpsi (npsi) = 0.0
       triangnpsi_l (npsi) = 0.0
       pindentnpsi(npsi) = 0.0
       grho1npsi  (npsi) = 0.0
       grho2npsi  (npsi) = 0.0
       sfareanpsi (npsi) = 0.0
c
c ----------------------------------------------------------------------
c Houlberg's Metrics on axis defined here. Since bootstrap is zero on
c axis, it is ok to use metric(npsi)=metric(npsi-1) for axis values
c ----------------------------------------------------------------------
c
      bsq    (npsi) = bax * bax
      bmsq   (npsi) = 1.0 / bsq(npsi)
      grbmsq (npsi) = 0.0
      grth   (npsi) = 0.0
      do m=1,3
        gfm(m,npsi) = 0.0
      end do
c
      cer_bp  (npsi) = 0.0
      cer_btdr(npsi) = bax / xmagax
c
****  grth (  npsi) = grth(npsi-1)
****  do m=1,3
****    gfm(m,npsi) = gfm(m,npsi-1)
****  end do
c
c ----------------------------------------------------------------------
c to get vprime at the magnetic axis we do the following
c (q is not known at this point so 4 * pi**2 * q/(f*<1.0/R**2> is
c not used for vprime at this stage)
c the geometric form is:
c   vprime(npsi) = pisq * (1.0 + elong(npsi)**2)
c .              / (u0 * ABS (curaxis) * 2.0 * elong(npsi))
c the following form for vprime(npsi) is obtained by
c integrating vprime from the magnetic axis out to the first
c flux surface wrt psi and setting the result equal to the
c volume at that point (which is known):
c ----------------------------------------------------------------------
c

      vprime(npsi) = 2.0 * psivolp(npsi-1) /
     .              (psivlcpy(npsi)-psivlcpy(npsi-1)) - vprime(npsi-1)
      if (ixcal .eq. 0)  go to 350
      xhm2(npsi)   = 1.0
      xi11(npsi)   = 0.0
      xi33(npsi)   = 0.0
      xips(npsi)   = 0.0
c
c    get bsqinvavg on axis by polynomial extrapolation:
c     use a cubic polynomial
c
      xxp(1)=psivlcpy(npsi-1)
      xxp(2)=psivlcpy(npsi-2)
      xxp(3)=psivlcpy(npsi-3)
      xxp(4)=psivlcpy(npsi-4)
      yyp(1)=bsqinvavg(npsi-1)
      yyp(2)=bsqinvavg(npsi-2)
      yyp(3)=bsqinvavg(npsi-3)
      yyp(4)=bsqinvavg(npsi-4)
      xxp(5)=psivlcpy(npsi)          !evaluation point
      call polint(xxp,yyp,nxxp-1,xxp(5),bsqinvavg(npsi),dyyp)
c
c     get bsq_avg on axis by polynomial extrapolation:
c     use a cubic polynomial
c
      xxp(1)=psivlcpy(npsi-1)
      xxp(2)=psivlcpy(npsi-2)
      xxp(3)=psivlcpy(npsi-3)
      xxp(4)=psivlcpy(npsi-4)
      yyp(1)=bsq_avg(npsi-1)
      yyp(2)=bsq_avg(npsi-2)
      yyp(3)=bsq_avg(npsi-3)
      yyp(4)=bsq_avg(npsi-4)
      xxp(5)=psivlcpy(npsi)          !evaluation point

      call polint(xxp,yyp,nxxp-1,xxp(5),bsq_avg(npsi),dyyp)

      yyp(1)=b_avg(npsi-1)
      yyp(2)=b_avg(npsi-2)
      yyp(3)=b_avg(npsi-3)
      yyp(4)=b_avg(npsi-4)
      xxp(1)=psivlcpy(npsi-1)
      xxp(2)=psivlcpy(npsi-2)
      xxp(3)=psivlcpy(npsi-3)
      xxp(4)=psivlcpy(npsi-4)
      xxp(5)=psivlcpy(npsi)         !evaluation point
      call polint(xxp,yyp,nxxp-1,xxp(5),b_avg(npsi),dyyp)




      yyp(1)=h_factr(npsi-1)
      yyp(2)=h_factr(npsi-2)
      yyp(3)=h_factr(npsi-3)
      yyp(4)=h_factr(npsi-4)
      xxp(1)=psivlcpy(npsi-1)
      xxp(2)=psivlcpy(npsi-2)
      xxp(3)=psivlcpy(npsi-3)
      xxp(4)=psivlcpy(npsi-4)
      xxp(5)=psivlcpy(npsi)         !evaluation point
      call polint(xxp,yyp,nxxp-1,xxp(5),h_factr(npsi),dyyp)

c
c ----------------------------------------------------------------------
c integrate bpl over psi, form the inductance li
c ----------------------------------------------------------------------
c
 350  call integ (psivlcpy,bpl,npsi,bpsq)
      if (u0 .eq. 0.0 .or. tocur .eq. 0.0 .or. psivolp(1) .eq. 0.0) then
        write  (6, 3233)  u0,tocur,psivolp(1)
 3233   format (' u0, tocur, psivolp(1) =', 3(2x, 1pe14.6))
        call STOP ('subroutine FLUXAV: problem #7', 156)
      end if
      ali1 = vprime(1)*bpsq / (u0*tocur*psivolp(1)) ! not used
      ali  = twopi*bpsq*circum*circum/(psivolp(1)*u0*tocur*u0*tocur)
c
      mdum = 0
      if (ndisk .gt. 0  .and. myid .eq. master )then
         write (ndisk, 8100) mdum
         write (ndisk,8100)npsi,nj
         write (ndisk, 8120)  (psivlcpy(i),i=1,npsi)
         write (ndisk, 8120)  100.*xmagax,100.*ymagax
         if ( iscr .ne. 1)  close (unit = ndisk)
      endif
c
c --- pick up psival (the edge value could have changed slightly)
c
      psival(1) = -psivlcpy(1)
c
c --- get psi at (rvloop, zvloop):
c

      if ((x(1) .le. rvloop .and. rvloop .le. x(nx)) .and.
     .    (y(1) .le. zvloop .and. zvloop .le. y(ny))) then
        call my_dbcevl1(x,nx,y,ny,cspln,nx,rvloop,zvloop,pds,ierr,1)
        psivloop = -twopi * (pds(1) + psival(1))
      else
        psivloop = 0.0
      end if
c
c     reset wzero (it was changed above in the torflux2 calculations)
c
      iflag = 2
     
      call zlim (wzero, nx, ny, ncontr, rcontr, zcontr, x, y, iflag)
      cconst = -1.0
      call multpl1 (p, nx*ny, cconst) ! change p back to original sign
c
c     apply psifctr only once if MHD calcs are not done:
c
      if (ifixshap .eq. 1)  psifctr = 0.0

      return
c
 1000 if (iounit .ne. 0) then
        write  (iounit, 1010) ierr,myid
 1010   format (' CNTOUR, called by FLUXAV, returned with ierr =', i5 /
     .          ' program must stop',/,' process id = ',i5,/)
        write  (iounit, 1011)  rcmin, rcmax, zcmin, zcmax
 1011   format ( '  rcmin,rcmax,zcmin,zcmax =',4(2x,1pe14.6))
        write  (iounit, 1012)  isetarcl, j, mp
 1012   format ( '  isetarcl,j,mp =',3(2x,i5))
        write  (iounit, 1013)  (psivlcpy(i), i=1,npsi)
 1013   format ( '  psivlcpy =' / (5(2x,1pe14.6)))
      end if
c
c if an error exit by way of label 1000 is taken then we still have to
c change the sign of psi back to what it used to be
c
      cconst = -1.0
      call multpl1 (p, nx*ny, cconst)
      call dump_psi_values (p, nx, ny, x, y, n77, -1, zero,
     .                      0, dumy, .false., dumy, 0, ptrace)
      call STOP ('subroutine FLUXAV: problem #8', 14)
      return
c
      end

      subroutine fofpsi (flimarg)
c
c
c ----------------------------------------------------------------------
c  given ffprim on psir grid,calculate f(psi) on psival grid
c  flimarg is the known boundary value of f on the plasma surface:
c --- input assumes:
c  psir(1) = mag axis,psir(nj) = plasma edge
c  psival(1) = edge,pisval(npsi) = axis
c --- on output we have:
c   fpsi(1)  =   f(psival(1)) = flimarg=edge value
c   fpsi(npsi) = f(psival(npsi)) = axis value
c ----------------------------------------------------------------------
c
      USE param
      USE mhdpar 
      USE numbrs
      USE machin
      USE psig
      USE rhog
      USE flxav
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'flxav.i'
c      include 'machin.i'
c      include 'numbrs.i'
c      include 'psig.i'
c      include 'rhog.i'
      include 'storage.i'
c
      isignn = 1
      if (btor .lt. 0.0)  isignn = -1
c
c --- interpolate ffprim(psir) to get ffprim(psival) ( = ffppsival)
c
      call intrp (1, 1, psir, ffprim, nj, psival, ffppsival, npsi)
       
c
c --- integrate from plasma edge to plasma axis
c
      call trap2 (psival, ffppsival, fpsi, npsi)
      do j=1,npsi
        fpsi(j) = isignn * SQRT (flimarg * flimarg + 2.0 * fpsi(j))
      end do


c
c --- Now interpolate qpsir(psir) to get qpsival(psival)
c
      call intrp (1, 1, psir, qpsir, nj, psival, qpsival, npsi)
      return
c
      end

      subroutine fqlin (x1, y1, x2, y2, f1, f2, f3, f4, x, y,
     .                  area, psivl)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      a1    = (x2-x) * (y2-y)
      a2    = (x-x1) * (y2-y)
      a3    = (x-x1) * (y-y1)
      a4    = (x2-x) * (y-y1)
      psivl = (a1*f1 + a2*f2 + a3*f3 + a4*f4) / area
      return
c
      end

      subroutine freebdry (maxit, iounit, ind, relerr)
c

c
c ----------------------------------------------------------------------
c   routine does free boundary equilibrium calculation
c   The solution is iterated until converged. As determined in
c   subroutine CONCEK, the convergence condition is
c      (sum over i,j (psi(i,j,iter-1)-psi(i,j,iter))/
c                  (psiaxis(iter)-psibdry(iter)))**2 < toleq
c   if convergence is not achieved in maxiter iterations the solution is
c   terminated and ind = 2 is returned. If the solution is converged before
c   maxiter iterations are exhausted,ind = 1 is returned.
c   The iterations are required because j toroidal depends on psi.
c   This introduces three sources of non-linearity:
c               a) the support set of j is unknown a priori
c               b) the boundary conditions on the border of the
c                  MHD grid are unknown a priori
c               c) ffprim and pprim are functions of psi
c   The fixed point iteration scheme used here to handle these problems
c   is dynamically relaxed using the krylov sequence of residuals to
c   discover the spectral radius of the iteration matrix and then define
c   the relaxation parameter in such a way as to annihilate the largest
c   eigenvalue (of course this will change the spectral radius of the
c   iteration scheme each time the relaxation parameter is redefined,
c   independent of the fact that the iterations are not stationary,
c   consequnetly and updating scheme,with a trust region,is employed).
c   both p prime and ff prime are renormalized to get proper total
c   current
c
c --- input
c --- through argument list:
c  maxiter      max number of iterations allowed (if ii = 1)
c  iounit       > 0 write pertinent output as a function of iteration
c                   number to unit number iounit
c              =0 do not write this information
c --- input through INCLUDE files,only those variables in the INCLUDE files that
c --- are used in this subroutine are listed below:
c  INCLUDE file soln2d:
c  toleq      tolerance parameter for setting convergence definition (see above)
c  ieq          =0 implies initial equilibrium, .ne. 0 implies subsequent one
c  INCLUDE file etc:
c  tocur             total toroidal current at this time,in amps
c  INCLUDE file mhdcom:
c  rma,zma              (r,z) coords of magnetic axis (input/output)
c  psibdry              psi value on plasma boundary  (input/output)
c  psiaxis              psi value on plasma magnetic axis (input/output)
c                       stored at the end of these vectors (see subroutine INIT)
c  psi(i,j)             i = 1..nw,j=1...nh,(input/output)
c  psi1d(k = 1,..nwh)     1d version of psi,indexed as k = (i-1)*nh+j
c
c INCLUDE file limiter.i:
c  xlimiter(i)            i = 1...nlimiter,the limiter points vectors
c  ylimiter(i)            the extremes of the limiter are assumed
c                          (which is not the fortran default)
c
c  INCLUDE file mhdgrid:
c  rmhdgrid(1..nw)
c  zmhdgrid(1..nh)         r,z vectors which define mhdgrid
c
c  INCLUDE file mhdpar:
c  nw,nh,nwh               MHD grid specifiers
c  n2cspln,nh2,nwork       parameters for work and bicubic spline
c
c  INCLUDE file contour:
c  nconmax                 max allowed points in rcontr,zcontr
c
c  INCLUDE file zerocom:
c  zero(i)        i = 1,...nwh,the interior/exterior indicator for limiter
c
c  INCLUDE file io:
c  ncrt
c  nout
c  nqik           fortran unit numbers for output writes
c
c  INCLUDE file storage:
c  xdum(i)
c  wdum(i)
c  zdum(i)          i = 1,...nj, temporary storage vectors
c
c  INCLUDE file numbrs:
c  nj                size of vectors pprim,ffprim,psir
c
c  INCLUDE file bicube:
c  cspln(n2cspln,nw,nh2)   bicubic spline quantitites,need not be set
c  wnoperm(nwork)          on input and do not contain anything useful
c  pds(6)                  on output.
c
c  INCLUDE file rhog:
c  pprim(i)
c  ffprim(i)            i = 1,...nj   dpress/dpsi,f*df/dpsi
c  psir(i)    psi vector over which ppprim,ffprim are defined,input/output
c
c  INCLUDE file mhdbcdtn:
c  fixfcoil              switch for fcoil fitting,may change dynamically
c
c INCLUDE file imsl
c imslmd                   error control for local IMSL routines
c
c --- output:
c --- through argument list
c  ind     =0 implies convergence to toleq was not achieved and
c             iterations were terminated due to vanishing of plasma
c             current (see also ieqfail below)
c          =1 implies convergence to toleq in less than maxiter iterations
c             was achived
c          =2 implies convergence to toleq was not achived in maxiter iterations
c --- through INCLUDE files:
c  INCLUDE file mhdcom:
c  psi(i,j)          the new computed psi
c  psibdry
c  psiaxis           boundary and axis values of psi
c  pcurrent(i)       the toroidal plasma current density,amps/m**2
c
c  INCLUDE file soln2d:
c  ieqfail  = 0 implies free boundary calculation was successful
c  ieqfail  = 1 implies support set for plasma current vanished (i.e., during
c               the iteration process the plasma shrunk to nothing. We
c               have not seen this happen using subroutine bound as the
c               current support set generator. But it is known to occur
c               in other scrape-off algorithms).
c  ieqfail  = 2 implies magnetic axis could not be found
c
c  INCLUDE file zerocom:
c  wzero(k)          k = 1,..nwh, wzero(k) = 1 if plasma current exists at
c                    rmhdgrid(i),zmhdgrid(j). wzero(k) = 0 otherwise
c                    (k = (i-1)*nh+j)
c
c  INCLUDE file contour:
c  rcontr(i)        i = 1,...ncontr    the plasma boundary points
c  zcontr(i)
c  ncontr
c  rplasmin
c  rplasmax
c  zplasmin
c  zplasmax             min/max extension of plasma
c  rsep
c  zsep                    (r,z) coords of xpoint (if none xsep = ysep = 0.0)
c
c  INCLUDE file small:
c  xax(1),yax(1)        the magnetic axis (same as rma,zma)
c  elongax             elongation (of psi contour) at magnetic axis
c
c  INCLUDE file rhog:
c  psir(j)              j = 1,..nj   new psi vector for pprim,ffprim
c  pprim(j)
c  ffprim(j)            newly normalized pprim,ffprim
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE io
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE numbrs
      USE mhdgrid    
      USE mesh
      USE machin
      USE geom
      USE constnts
      USE soln2d
      USE psig
      USE rhog
      USE mhdcom
      USE bicube
      USE flxav
      USE neo2dp
      USE etc
      USE gpsi
      USE mhdbcdtn
      USE replace_imsl,       ONLY : my_ibcccu,my_dbcevl1 
      implicit  integer (i-n), real*8 (a-h, o-z)

c      include 'small.i'
      include 'storage.i'
c      include 'zerocom.i'
c      include 'mhdbcdtn.i'
      include 'imsl.i'
c
      parameter   (maxdrift = 5)   ! used in vertical stabilization only
      dimension    ivertidx(maxdrift)
      character*6  intfl, icurflg*5
      dimension    resid1(nwh), resid2(nwh), wzerosv(nwh), pcurold(nwh)
c
c ----------------------------------------------------------------------
c initialization
c if ieq = 0 we know the intial psi from initmhd or from reading an
c eqdsk. if ieq .ne. 0 then the inital psi is the converged value
c of psi calculated at the previous time point. in either case
c we have an initial psi to work with.
c ----------------------------------------------------------------------
c
      itmin    = iteq / 3
      maxiter  = maxit
      itmin    = MAX0 (itmin, 5)
      xchisqmn = 1.0d100
      ifxwzero = 0
      iteqsv   = iteq
      noxpoint =   0 ! =0 may search for xpoint,=1 don't
      isignpsi = -1 ! sign of psi on axis
      isgngren = -1 ! multiplier for some Green's functions, see getfcur
      mfcinv   =  0 ! means calc inverse of design matrix in getfcur
      kstart   =  1
      kend     =  nwh  ! use j=kstart,kend NOT j=1,nwh, adj. dynamically
      iavfcoil =  0 ! switch for averaging f coil currents,see GETFCUR
      idifsum  =  0 ! sensor for vertical instability
      ifixbdry =  0 ! switch for fixed boundary, see subroutine SOLVEGS
      itervert = 10 ! start vert. instab. sensing after itervert iters
      itencoil =  0 ! iteration ctr for no coils option
      fixfcopy = fixfcoil   ! save for restoration at exit from freebdry
      kvertsbl =  0         ! index for ivertidx
      if (mhdmode .eq. 'no coils' .or. mhdmethd .eq. 'green')
     .ivertsbl = 0
      drmhdgrd = (rmhdgrid(2)-rmhdgrid(1))
      dzmhdgrd = (zmhdgrid(2)-zmhdgrid(1))
      darea    =  drmhdgrd*dzmhdgrd
      drhalf   = 0.5 * drmhdgrd
      dzhalf   = 0.5 * dzmhdgrd
      drdzhalf = drhalf*dzhalf
      alpha    = 1.0
      alphao   = 1.0     ! pprim and ffprim normalization factors
      alph2    = 1.0     ! cumulative multiplier, for printout only
      if (omeq .eq. 0.0) then
        omeqmin = 0.3
        omeqmax = 1.5
        omeq    = 1.0
        iomeq   = 1
        iomeqit = 3
        iomeqct = 0
      else
        iomeq   = 0
      end if
      ind     = 0
      iter    = 0
      relerr  =  1.0e6        ! initial guess of relative error
      ichklim = -1            ! check limiter flag for subroutine XPOINT
      modbdry =  0    ! flag sets plasma boundary refinement calculation
      if (rma .eq. 0.0)  rma = 0.5 * (rmhdgrid(nw)+rmhdgrid(1))
      if (zma .eq. 0.0)  zma = 0.5 * (zmhdgrid(nh)+zmhdgrid(1))
      xsep    = 0.0
      ysep    = 0.0           ! initial guess of xpoint
c
c --- following calcs are to set up wzero, which is not known if ieq = 0
c
      psiaxlocl=psiaxis
      if (ieq .eq. 0) then
        if (isignpsi .lt. 0) then
          cconst = -1.0
          call multpl1 (psi, nwh, cconst)
          psibdry = -psibdry
          psiaxlocl=-psiaxis
        end if
        delta_psi=(psiaxlocl-psibdry)/(nw-1)
c
c --- get bicubic spline representation of psi:
c
         imslmd = 'ibcccu'
         imslmd='2419c208'
         call my_ibcccu (psi,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                   wnoperm,ierr)
c
c --- trace the contour corresponding to the plasma boundary
c --- psi(r,z) = psibdry:
c
         dang    = 5.0
         dang    = 0.1 ! 999,
         arcl    = 0.02
         arcl    = 0.001 ! 999,
         bperr   = 0.05
         bperr    = 1.e-3 ! 999,
         iauto   = 1
         xsearch = 0.0
         ysearch = 0.0
c
         call cntour (rma,zma,psibdry,rplmin,rplmax,zplmin,zplmax,
     .         zrplmin,zrplmax,rzplmin,rzplmax,dang,arcl,bperr,xsearch,
     .         ysearch,
     .         xlimiter(nlimiter+1),xlimiter(nlimiter+2),
     .         ylimiter(nlimiter+1),ylimiter(nlimiter+2),
     .         iauto,iautoc,
     .         rcontr,zcontr,ncontr,rmhdgrid,nw,zmhdgrid,nh,cspln,
     .         n2cspln,nh2,ncrt,nconmax,ierr,wzero,0,delta_psi)
c
         kstart = 0
         do 20 i=1,nw
           do 20 j=1,nh    ! require this order in i,j
c                            for kstart,kend logic
             k0        = (i-1)*nh+j
             wzero(k0) = 0.0
             if (psi(i,j)    .lt. psibdry)  go to 20
             if (zmhdgrid(j) .lt.  zplmin)  go to 20
             if (zmhdgrid(j) .gt.  zplmax)  go to 20
             if (rmhdgrid(i) .lt.  rplmin)  go to 20
             if (rmhdgrid(i) .gt.  rplmax)  go to 20
             kend = k0
             if (kstart .eq. 0)  kstart = k0
             wzero(k0) = 1.0
   20    continue
        if (isignpsi .lt. 0) then
          cconst = -1.0
          call multpl1 (psi,nwh,cconst)
          psibdry = -psibdry
        end if
      end if    ! end ieq = 0 branch
c
      write (ncrt,10)
   10 format (//'   free boundary equilibrium calculation' //
     .   '  it',4x,'reler',3x,'current',2x,
     .  '# pts',3x,'psibdry',3x,'psiaxis',3x,'rplmin',
     .  5x,'zplmin',3x,'alpha',5x,'curaxis',7x,'  xax',10x,'yax')
      write (nout  , 10)
      write (nqik  , 10)
      if (iounit .ne. 0)
     .write (iounit, 10)
c
c ----------------------------------------------------------------------
c begin iteration
c these are the "inner" iterations. Even though the current
c profile is known as a function of rho,the support set of
c the current as a function of (r,z) is not known and must be
c determined by evaluating psi(r,z),in such a manner that psi
c is consistent with J(rho). In other words,the mapping
c rho(psi(r,z))is the real object of this calculation.
c ----------------------------------------------------------------------
c
      tsolveg = 0.0
      tsolvav = 0.0
      tbavg   = 0.0
      txptavg = 0.0
      ixptcl  = 0
   30 iter    = iter + 1
****  write (6, 99) iter-1,relerr,tsolveg,imax,jmax,omeq,rma,zma,xchisq,
**** .              tocur
** 99 FORMAT ( '  ITER = ',I5,'  RELERR=',1PE10.2                /
**** .         '  tsolveg =',1pe12.6,' IMAX JMAX =',2(2X,I3)   /
**** .         '  OMEQ= ',1PE12.4, '  RMA ZMA =',2(1PE12.5,1X) /
**** .         '  XCHISQ =',1PE12.5,' TOTAL CURRENT,AMPS ',1PE12.5)
      if (fixfcoil .ne. fixfcopy)  itencoil = itencoil + 1
      do 45 j=1,nh
      do 45 i=1,nw
      k0 = (i-1)*nh+j
      psi1d(k0) = psi(i,j)
   45 p(i,j)    = psi(i,j)    ! save previous psi
c
c ----------------------------------------------------------------------
c calculate the plasma current in AMPS
c NOTE THAT MULTIPLACTION BY DAREA IS DONE HERE TO AVOID HAVING
c TO DO IT IN SEVERAL OTHER PLACES. THUS PCURRENT IS NOT A CURRENT
c DENSITY. IT IS THE TOTAL CURRENT ASSOCIATED WITH A FILAMENT OF CROSS
c SECTIONAL AREA DAREA.
c pprim and ffprim are defined on the psir(1...nj) psi grid.
c to get j(r,z) we use a natural cubic spline representation
c of pprim and ffprim to interpolate onto (r,z) by going through psi(r,z):
c ----------------------------------------------------------------------
c
      do 50 j=1,nj
        pprim(j) = alpha*pprim(j)
   50   ffprim(j) = alpha*ffprim(j)
c
      if (ifxwzero .eq. 1)  call copya (pcurrent, pcurold, nwh)
      call zeroa (pcurrent, nwh)
      call zeroa (wnoperm , nwh)
      k0 = kend - kstart + 1
c
****  imslmd = 'freebd1'
****  call intrp (1,1,psir,pprim,nj,psi1d(kstart),
**** .            pcurrent(kstart),k0)
****  imslmd = 'freebd2'
****  call intrp (1,1,psir,ffprim,nj,psi1d(kstart),
**** .            wnoperm(kstart),k0)
c
c --- try linear interpolation
c
      imslmd = 'freebdry'
      do 51 k=kstart,kend
        if (wzero(k) .eq. 0.0)  go to 51
        psiwnt = psi1d(k)
        call find(kjl,kju,psiwnt,psir,nj)
****    if (kjl .eq. 0)  go to 51
        if (kjl .eq. 0) then
          kjl = nj
          kju = nj
        end if
        if (kjl .eq. kju) then
          pcurrent(k) = pprim (kjl)
          wnoperm (k) = ffprim(kjl)
        else
          delpsi      = (psiwnt-psir(kjl))/(psir(kju)-psir(kjl))
          delp        = pprim(kju)-pprim(kjl)
          delf        = ffprim(kju)-ffprim(kjl)
          pcurrent(k) = pprim(kjl)+delp*delpsi
          wnoperm (k) = ffprim(kjl)+delf*delpsi
        end if
   51 continue
c
      do k0=kstart,kend
        i = k0 / nh + 1
        j = k0 - (i - 1) * nh
        if (j .eq. 0)  i = i - 1
        pcurrent(k0) = -(pcurrent(k0)*rmhdgrid(i)
     .                 + wnoperm(k0)/(u0*rmhdgrid(i)))*wzero(k0)*darea
      end do
c
c ----------------------------------------------------------------------
c --- average the plasma current if ifxwzero = 1
c ----------------------------------------------------------------------
c
      if (ifxwzero .eq. 1) then
        call adda(pcurold,pcurrent,nwh)
        cconst = 0.5
        call multpl1 (pcurrent,nwh,cconst)
      end if
c
c ----------------------------------------------------------------------
c --- subroutine SOLVEGS sets up and solves Grad-Shafranov equation
c ----------------------------------------------------------------------
c
****  if (iter .gt. 20)  iavfcoil = 1   ! avg old and new fcoil currents
      call SECOND (tsolvegs)
      call solvegs (iounit, ierr, kstart, kend, ifixbdry, iavfcoil,
     .              nout, ncrt)
      call SECOND (tsolveg1)
      tsolveg = tsolveg1-tsolvegs
      tsolvav = tsolvav+tsolveg
      if (ierr .ne. 0)  go to 1100
      do 100 j=1,nh
         do 100 i=1,nw
            k = (i-1)*nh+j
            psi(i,j) = psi1d(k)
  100 continue
c
c ----------------------------------------------------------------------
c underrelax solution, estimate eigenvalue of residuals to predict omeq:
c ----------------------------------------------------------------------
c
      if (iter .gt. 5) then
        do 110 j=1,nh
        do 110 i=1,nw
        psi(i,j) = omeq*psi(i,j) + (1.0-omeq)*p(i,j)
        k0 = (i-1)*nh+j
  110   psi1d(k0) = psi(i,j)
      end if
c
      if (iomeq .eq. 1) then
      do 111 j=1,nh
      do 111 i=1,nw
        k0 = (i-1)*nh+j
        resid2(k0) = resid1(k0)        ! residual for previous iteration
  111 resid1(k0) = psi(i,j) - p(i,j)   ! current residual
      if (iter .gt. 5) then
        sum1 = 0.0
        sum2 = 0.0
        do 113 j=1,nwh
          sum1 = sum1+resid1(j)*resid2(j)
  113     sum2 = sum2+resid2(j)*resid2(j)
        eigstar  = sum1/sum2            ! estimate of largest eigenvalue
c                                         of iteration matrix
        omeqpred = omeq/(1.0-eigstar) ! to annihilate largest eigenvalue
c
c --- the trust reqion for omeq is given by (omeqmin,omeqmax):
c
        omeqpred = MIN (omeqpred,omeqmax)
        omeqpred = MAX (omeqpred,omeqmin)
        iomeqct = iomeqct+1
        if (iomeqct .ge. iomeqit) then
          iomeqct = 0
          omeq = omeqpred
        end if
      end if
      end if
c
c ----------------------------------------------------------------------
c find magnetic axis (there is only one for non-doublets)
c note: MAGAX destroys wnoperm
c iknowax = 1 is returned if magax found an appropriate magnetic axis
c iknowax = 0 is returned if magax failed to find the magnetic axis
c in this case we abort the calculations:
c ----------------------------------------------------------------------
c
  130 ispln   = 0  ! set up cspln in subroutine MAGAX..
c                  ..cspln is used below also
      iknowax = 0
      if (iter .gt. 5)  iknowax = 1
      zmaold = zma
      call magax (psi,rmhdgrid,zmhdgrid,nw,nh,isignpsi,iknowax,
     .            ncrt,ispln,zero,rma,zma,cspln,n2cspln,nh2,
     .            psiaxis,elongax)
      if (iknowax .eq. 0) then
        call dump_psi_values (psi,nw,nh,rmhdgrid,zmhdgrid,n77,
     .                isignpsi,zero,0,map,.false.,dumy,0,ptrace)
        call STOP ('subroutine Freebdry: problem #8', 14)
        ieqfail = 2   ! magnetic axis not found..
c                     ..so set no equilibrium failure flag
        return        ! abort the calculation
      end if
c
c ----------------------------------------------------------------------
c get the free boundary,rcontr(i),zcontr(i),i = 1..ncontr
c ----------------------------------------------------------------------
c
      dpsi     = 1.0e10
      rminsrch = xlimiter(nlimiter+1) - (rmhdgrid(2)-rmhdgrid(1))
      itrace   = 1
      limfag   = 2
      isignn   = isignpsi
      if (ieq .eq. 0 .and. iter .eq. 1)  radold = 0.5 * (rma+rminsrch)
      call SECOND (tbound)
      call bound (psi1d,nw,nh,nwh,zero,rmhdgrid,zmhdgrid,rma,zma,
     .            rminsrch,itrace,nlimiter,xlimiter,ylimiter,isignn,
     .            limfag,radold,nconmax,rcontr,zcontr,ncontr,dpsi,
     .            rplmin,rplmax,zplmin,zplmax,
     .            rzplmin,rzplmax,zrplmin,zrplmax,psibdry)
      call SECOND (tbound1)
      tbound = tbound1-tbound
      tbavg  = tbavg+tbound
      if (isignn .ne. 0)  go to 1000    ! error exit
      if ((isignpsi .gt. 0 .and. psibdry .gt. psiaxis)   .or.
     .    (isignpsi .lt. 0 .and. psibdry .lt. psiaxis)) then
        write (nout, 1101)  isignpsi, psibdry, psiaxis
        write (ncrt, 1101)  isignpsi, psibdry, psiaxis
        if (iounit .ne. 0) then
          write  (iounit, 1101)  isignpsi, psibdry, psiaxis
 1101     format (' subroutine FREEBDRY detects an ERROR'             /
     .            ' isignpsi, psibdry, psiaxis =', i5, 2(2x, 1pe14.6) /
     .            ' program must stop')
        end if
        call STOP ('subroutine FREEBDRY: problem #3', 94)
      end if
      if (noxpoint .ne. 0)  dpsi = 0.0
      if (dpsi .lt. 0.01 * ABS (psiaxis-psibdry))  go to 122  ! x point?
c
c ----------------------------------------------------------------------
c subroutine bound claims a diverted plasma was found
c use refined two-d search to locate x point exactly
c ichklim = 1 means retrace the contour found  by bound using the
c bi-cubic spline representation of psi and subroutine CNTOUR.
c once ichklim is set to 1 we leave it that way to avoid oscillations
c ----------------------------------------------------------------------
c
      ispln    = 1                              ! cspln was set in MAGAX
      iknowxpt = 0
      if (iter .gt. 5)  iknowxpt = 1   ! XPOINT will correct if not true
      call SECOND (txpoint)
      ixptcl = ixptcl + 1
      call xpoint (psi,rmhdgrid,zmhdgrid,nw,nh,xsep,ysep,iknowxpt,
     .             cspln,n2cspln,nh2,wnoperm,nwork,
     .             ispln,xlimiter,ylimiter,nlimiter,0,
     .             nconmax,ichklim,rcontr,zcontr,ncontr,rma,zma,
     .             xdum,ydum,psisep,rplmin,rplmax,zplmin,zplmax,
     .             rzplmin,rzplmax,zrplmin,zrplmax,wdum)
      call SECOND (txpoint1)
      txpoint = txpoint1-txpoint
      txptavg = txptavg+txpoint
      if (iknowxpt .eq. 1)  psibdry = psisep
****  if (isignpsi .gt. 0 .and. psibdry .gt. psiaxis) go to 1100 ! abort
****  if (isignpsi .lt. 0 .and. psibdry .lt. psiaxis) go to 1100 ! abort
      if ((isignpsi .gt. 0 .and .psibdry .gt. psiaxis) .or.
     .    (isignpsi .lt. 0 .and. psibdry .lt. psiaxis)) then
            write (nout, 1102)  isignpsi, psibdry, psiaxis
            write (ncrt, 1102)  isignpsi, psibdry, psiaxis
            if (iounit .ne. 0)
     .      write  (iounit, 1102)  isignpsi, psibdry, psiaxis
 1102       format ('  after call to XPOINT an error exists'           /
     .              '  isignpsi,psibdry,psiaxis =', i5, 2(2x, 1pe14.6) /
     .              '  program must stop')
            call STOP ('subroutine FREEBDRY: problem #4', 96)
      end if
c
c ----------------------------------------------------------------------
c wzero is now set like array zero except the support set of wzero
c is reduced to the box defined by rplmin,rplmax and zplmin,zplmax.
c ----------------------------------------------------------------------
c
  122 kstart = 0
      do 200 i=1,nw
      do 200 j=1,nh       ! require this i,j order for kstart,kend logic
      k0        = (i-1)*nh+j
      wzero(k0) = 0.0
      if ((isignpsi .gt. 0) .and. (psi(i,j) .lt. psibdry))  go to 200
      if ((isignpsi .lt. 0) .and. (psi(i,j) .gt. psibdry))  go to 200
      if (zmhdgrid(j) .lt. zplmin)  go to 200
      if (zmhdgrid(j) .gt. zplmax)  go to 200
      if (rmhdgrid(i) .lt. rplmin)  go to 200
      if (rmhdgrid(i) .gt. rplmax)  go to 200
      if (kstart .eq. 0) kstart = k0
      kend      = k0
      wzero(k0) = 1.0
  200 continue
c
c ----------------------------------------------------------------------
c  we need to rescale the psir grid given the new psi values.
c  the scaling is as follows. let pnew be the new psgrid and
c  pold be the old psi grid.
c  then  pnew = a*pold+b  is the linear transformation that is used.
c  note that,for example,dpress/dpnew = (1.0/a)*dpress/dpold
c  with this scaling. hence the effect is to multiply the
c  current density by (1.0/a). since the total current must
c  be kept constant,the support set of the current density
c  must change to compensate for this (1/a) factor.
c  note that
c       a = (psinew(edge)-psinew(axis))/(psiold(edge)-psiold(axis))
c  a is a measure of how much the current must be changed due to
c  changes in magnitude of psi. below we renormalize again,using
c  alpha. alpha takes care of the changes in the support set
c  of the current density. the ratio of a/alpha is a combined
c  measure of how much the values of psi versus the geometry is having.
c ----------------------------------------------------------------------
c
      if (isignpsi .lt. 0) then
        pmin = psiaxis
        pmax = psibdry
      else
        pmin = psibdry
        pmax = psiaxis
      end if
        psi1 = psir(1)
        do 124 j=1,nj
        percnf = (psir(j)-psi1)/(psir(nj)-psi1)
  124   psir(j) = psiaxis+percnf*isignpsi*(pmin-pmax)
      pcurext = 0.0
      if (modbdry .ne. 0 .and. iter .gt. 5) then
c
c ----------------------------------------------------------------------
c --- adjust wzero for the boundary points to get a smoother transition
c --- pcurext is the plasma current (in amps) which is due to cells whose
c --- center (at rmhdgrid(i),zmhdgrid(j) ) is outside the plasma but part
c --- of the cell intersects the plasma. this current is not assigned to
c --- any grid point but is used in the current normalization calculation.
c ----------------------------------------------------------------------
c
      ksep  = 0
      ksep1 = 0
      if (xsep .ne. 0) then
        isep  = (xsep-rmhdgrid(1))/drmhdgrd
        jsep  = (ysep-zmhdgrid(1))/dzmhdgrd
        ksep  = (isep-1)*nh+jsep
        ksep1 = ksep+nh
      end if
      do 500 i=2,nw-1
      do 500 j=2,nh-1
        k = (i-1)*nh+j
        if (k .eq. ksep )  go to 500   ! don't do for x point
        if (k .eq. ksep1)  go to 500
        wsum = wzero(k-nh-1)+wzero(k-nh)+wzero(k-nh+1)+
     .         wzero(k+1)+wzero(k+nh+1)+wzero(k+nh)+
     .         wzero(k+nh-1)+wzero(k-1)
        if (wsum .gt. 7.95)  go to 500    ! not a boundary point
        if (wsum .lt. 0.05)  go to 500    ! not a boundary point
c
c --- other points are possible boundary points so run them through areaclc:
c
        psinw = 0.25*(psi(i-1,j)+psi(i-1,j+1)+psi(i,j+1)+psi(i,j))
        psine = 0.25*(psi(i,j)+psi(i,j+1)+psi(i+1,j+1)+psi(i+1,j))
        psise = 0.25*(psi(i,j-1)+psi(i,j)+psi(i+1,j)+psi(i+1,j-1))
        psisw = 0.25*(psi(i-1,j-1)+psi(i-1,j)+psi(i,j)+psi(i,j-1))
c
c --- using the above definition of psi values you can identify the cell
c --- which is under investigation.
c
        psik = psi(i,j)
        call areaclc (psisw,psinw,psine,psise,psibdry,psik,
     .                drmhdgrd,dzmhdgrd,drhalf,dzhalf,drdzhalf,ncrt,
     .                dareain,dareaout,psictr,rctr)
        rctr = rctr + rmhdgrid(i)   ! change rctr to MHD grid coordinate
        if (dareain .eq. 0.0) then      ! point k is outside the plasma
          if (dareaout .ne. 0.0) then   ! part of cell intersects plasma
          call find (kjl, kju, psictr, psir, nj)           ! get kjl,kju
          if (kjl .eq. 0)  go to 500
          if (kjl .eq. kju) then
            ppextra  = pprim(kjl)
            ffpextra = ffprim(kjl)
          else
            delpsi   = (psictr-psir(kjl))/(psir(kju)-psir(kjl))
            delp     =  pprim(kju)-pprim(kjl)
            delf     =  ffprim(kju)-ffprim(kjl)
            ppextra  =  pprim(kjl)+delp*delpsi
            ffpextra =  ffprim(kjl)+delf*delpsi
          end if
            pcurext  = -(ppextra*rctr+ffpextra/(u0*rctr))*dareaout
     .                 + pcurext
          end if
        else   ! point k inside the plasma; assign current to grid point
          if (psictr .ne. psik) then        ! cell is only partly inside
          call find (kjl, kju, psictr, psir, nj)          ! get kjl, kju
          if (kjl .eq. 0)  go to 500
            if (kjl .eq. kju) then
              ppextra = pprim(kjl)
              ffpextra = ffprim(kjl)
            else
              delpsi = (psictr-psir(kjl))/(psir(kju)-psir(kjl))
              delp = pprim(kju)-pprim(kjl)
              delf = ffprim(kju)-ffprim(kjl)
              ppextra = pprim(kjl)+delp*delpsi
              ffpextra = ffprim(kjl)+delf*delpsi
            end if
            pcurr1 = -(ppextra*rctr+ffpextra/(u0*rctr))*dareain
          call find (kjl, kju, psik, psir, nj)            ! get kjl, kju
          if (kjl .eq. 0)  go to 500
            if (kjl .eq. kju) then
              ppextra = pprim(kjl)
              ffpextra = ffprim(kjl)
            else
              delpsi = (psik-psir(kjl))/(psir(kju)-psir(kjl))
              delp = pprim(kju)-pprim(kjl)
              delf = ffprim(kju)-ffprim(kjl)
              ppextra = pprim(kjl)+delp*delpsi
              ffpextra = ffprim(kjl)+delf*delpsi
            end if
            pcurr2 = -(ppextra*rmhdgrid(i)+ffpextra/(u0*rmhdgrid(i)))
     .                *darea
            if (pcurr2 .ne. 0.0)
     .        wzero(k) = pcurr1/pcurr2
          end if
        end if
  500 continue
      end if    ! end modbdry boundary calculation branch
c
c ----------------------------------------------------------------------
c calculate total current (in amps)
c ----------------------------------------------------------------------
c
      pcuriter = 0.0
      icurpts  = 0
      area     = 0.0
      do 140 k0=kstart,kend
****  pcurrent(k0) = pcurrent(k0)*wzero(k0)
      if (wzero(k0) .gt. 0.0)  icurpts = icurpts+1
      area     = area+wzero(k0)
  140 pcuriter = pcuriter+pcurrent(k0)
      area     = area*darea
      pcuriter = pcuriter+pcurext
c
c --- try a different approach to pcuriter
c
      pcuriter = 0.0
      sumdlbar = 0.0
      do j=1,ncontr-1
        rbar = 0.5 * (rcontr(j)+rcontr(j+1))
        zbar = 0.5 * (zcontr(j)+zcontr(j+1))
        call my_dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,nw,rbar,
     .                                          zbar,pds,ier,3)
        dlbar    = SQRT ((rcontr(j+1)-rcontr(j))**2
     .                 + (zcontr(j+1)-zcontr(j))**2)
        bpbar    = SQRT (pds(2)**2+pds(3)**2)/rbar
        pcuriter = pcuriter+bpbar*dlbar
        sumdlbar = sumdlbar+dlbar
      end do
      pcuriter = pcuriter/u0
      if (pcuriter .gt. 0.0)  go to 144
      write  (ncrt , 4447)  ncontr, pcuriter
      write  (nout , 4447)  ncontr, pcuriter
      write  (nitre, 4447)  ncontr, pcuriter
 4447 format ('  ncontr =',i5, '  pcuriter= ',1pe12.4)
      ieqfail = 1
      return
  144 continue
      alphao = alpha
      alpha  = tocur/pcuriter
      if ((1.0-alpha)*(1.0-alphao) .lt. 0. .and. iter .gt. 5)
     .        alpha = 0.25*alphao+0.75*alpha
      alph2 = alph2*alpha
c
c --- get the relative error, relerr, and check for convergence
c
  240 relerro = relerr
c
c --- current density at the magnetic axis
c
      curaxis = -rma*pprim(1)-ffprim(1)/(u0*rma)
      call concek (psi,p,psiaxis,psibdry,nw,nh,toleq,ind,relerr,
     .             imax,jmax)
      if (iter .lt. 5)  ind=0  ! do at least min(5,maxiter) iterations
      if (ind .ne. 1 .and. iter .gt. maxiter) ind = 2
      if (ind .eq. 2)  itencoil = 0   ! if maxiters exceeded enable exit
      if (ind .gt. 0 .and. itencoil .eq. 0)  go to 160   ! no exit until
c                                                          itencoil = 0
      if (fixfcoil .eq. 0) then
        write  (intfl, 1200) icurpts
 1200   format ('+',i4)
      else
        write  (intfl, 1202) icurpts
 1202   format ('-',i4)
      end if
      read   (intfl, 1201) icurflg
 1201 format (a)
      if (iounit .ne. 0 .and. ieqprt .eq. 1) then
       write (iounit,150) iter,relerr,pcuriter,icurflg,psibdry,
     .psiaxis,rplmin,zplmin,alpha,curaxis,rma,zma,imax,jmax,area,pcurext
     .  ,xchisq,tocur,xsep,ysep,sumdlbar,pcurptsv
       write (nitre,150) iter,relerr,pcuriter,icurflg,psibdry,
     .psiaxis,rplmin,zplmin,alpha,curaxis,rma,zma,imax,jmax,area,pcurext
     .  ,xchisq,tocur,xsep,ysep,sumdlbar,pcurptsv
  150 format (' ',i3,1pe8.1,1pe14.7,a,2(1pe15.7),3(1pe15.7) /
     . 1x,3(1pe15.7),i5,i5,2x,4(1pe15.7) / 2x,6(1pe15.7))
      end if
      if (ivertsbl .eq. 1) then    ! vertical stability logic
        idifsign = 1
        if (zmaold .gt. zma)  idifsign = -1
        kvertsbl = kvertsbl+1
        if (kvertsbl .gt. maxdrift)  kvertsbl = 1
        ivertidx(kvertsbl) = idifsign
        idifsum = 0
        do 350 j=1,kvertsbl
  350   idifsum = idifsum+ivertidx(j)
c
c --- determine if fcoil currents should be frozen at present value.
c --- Error must be increasing (rather than decreasing),last maxdrift
c --- iterations must yield a drift in zma in the same direction,
c --- and at least itervert iterations must have been done (need
c --- itervert .ge. maxdrift):
c
        if (iter .gt. itervert .and. IABS (idifsum) .eq. maxdrift
     .                         .and. relerro       .lt.   relerr) then
            fixfcoil = 1
            ifixbdry = 1
          end if
        if (relerr .lt. toleq) then   ! switch back after convergence
           itencoil = 0
           fixfcoil = fixfcopy
           ifixbdry = 0
        end if
      end if
c
c ----------------------------------------------------------------------
c --- save the minimum chi-sq support set for the current
c ----------------------------------------------------------------------
c
      xchisq = 0.0
      if (mhdmode .eq. 'coils')  call exptlclc(kstart,kend)
      if (iter .gt. itmin) then
      xchisqmn = MIN (xchisq, xchisqmn)
      if (xchisqmn .eq. xchisq) then
         kstartsv = kstart
         kendsv   = kend
         pcurptsv = icurpts
         call copya (wzero,wzerosv,nwh)
      end if
      end if
      if ((iter .gt. 5+itmin .and. relerr .gt. 1.0e-04)  .or.
     .                               (ifxwzero .eq. 1)) then
         ifxwzero = 1
         kstart   = kstartsv
         kend     = kendsv
         icurpts  = pcurptsv
         call copya (wzerosv,wzero,nwh)
         if (iteq .eq. iteqsv)  iteq = iteq + itmin   ! allow additional
c                                                       itmin iterations
         maxiter = iteq
      end if
      go to 30
c
c ----------------------------------------------------------------------
c iterations over - do cleanup
c ----------------------------------------------------------------------
c
  160 iteq    = iteqsv
      tsolvav = tsolvav/iter
      tbavg   = tbavg/iter
      if (ixptcl .ne. 0)  txptavg = txptavg/ixptcl
      write  (ncrt, 151)  tsolvav, tbavg, txptavg
  151 format (
     . ' average time per iteration in SOLVEGS: ', 1pe12.6, ' seconds' /
     . ' average time per iteration in   BOUND: ', 1pe12.6, ' seconds' /
     . ' average time per iteration in  XPOINT: ', 1pe12.6, ' seconds' )
      if (iounit .ne. 0)
     .  write (iounit, 151)  tsolvav, tbavg, txptavg
      fixfcoil = fixfcopy    ! set back to original value
      write (ncrt,150) iter,relerr,pcuriter,icurflg,psibdry,psiaxis,
     .       rplmin,zplmin,alpha,curaxis,rma,zma,imax,jmax,area,pcurext,
     .       xchisq,tocur,xsep,ysep,sumdlbar,pcurptsv
      write (nout,150) iter,relerr,pcuriter,icurflg,psibdry,psiaxis,
     .       rplmin,zplmin,alpha,curaxis,rma,zma,imax,jmax,area,pcurext,
     .       xchisq,tocur,xsep,ysep,sumdlbar,pcurptsv
      write (nqik,150) iter,relerr,pcuriter,icurflg,psibdry,psiaxis,
     .       rplmin,zplmin,alpha,curaxis,rma,zma,imax,jmax,area,pcurext,
     .       xchisq,tocur,xsep,ysep,sumdlbar,pcurptsv
      if (iounit .ne. 0)
     .write (iounit,150) iter,relerr,pcuriter,icurflg,psibdry,psiaxis,
     .       rplmin,zplmin,alpha,curaxis,rma,zma,imax,jmax,area,pcurext,
     .       xchisq,tocur,xsep,ysep,sumdlbar,pcurptsv
  300 xax(1)   = rma
      yax(1)   = zma
      rplasmin = rplmin
      rplasmax = rplmax
      zplasmin = zplmin
      zplasmax = zplmax
      rsep     = xsep
      zsep     = ysep
      call copya (psi,p,nwh)
      if (iomeq .eq. 1)  omeq = 0.0
c
c --- calculate psiloop and probe values based on calculated
c --- currents,compare with exptl. data:
c
      if (mhdmode .eq. 'coils')  call exptlclc(kstart,kend)
c
c --- calculate integral of bpsquare dv:
c
      sum = 0.0
      do k=kstart,kend
        sum = sum+pcurrent(k)*wzero(k)*(psibdry-psi1d(k))
      end do
      bpsqtot = twopi*u0*sum
c
c --- change pcurrent into current density:
c
      do 400 k0=kstart,kend
  400   pcurrent(k0) = pcurrent(k0)/darea
      return
c
 1000 write (ncrt, 1010)  isignn
      write (nout, 1010)  isignn
      write (nqik, 1010)  isignn
 1010 format (/ ' BOUND unable to determine the free boundary'        /
     .          ' reason unknown, ONETWO must stop, error flag =', i5  )
      call STOP ('subroutine FREEBDRY: problem #1', 15)
 1100 write  (ncrt, 1110)
      write  (nout, 1110)
      write  (nqik, 1110)
 1110 format (' subroutine SOLVEGS detected an ERROR' /
     .        ' ONETWO must stop'                      )
      call STOP ('subroutine FREEBDRY: problem #2', 16)
c
      end

      subroutine getadjv (voladj, volume, volwant, volnudge, klim,
     .                    klimold, limpos, vadjmin, vadjmax)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character limpos*8
      data      xfact /1.0/
c
      nout    = 7
      holdadj = voladj
      if (limpos .eq. 'left' )  go to 3000
      if (limpos .eq. 'right')  go to 4000
      return
c
c ======================================================================
c we want plasma to limit on the left
c ======================================================================
c
 3000 if (klim .eq. 1 .and. klimold .eq. 1)  go to 3200
      if (klim .eq. 1 .and. klimold .eq. 0)  go to 3300
      if (klim .eq. 1 .and. klimold .eq. 2)  go to 3400
      if (klim .eq. 1 .and. klimold .eq. 3)  go to 3400
      if (klim .ne. 1 .and. klimold .eq. 1)  go to 3500
      if (klim .ne. 1 .and. klimold .eq. 0)  go to 3600
      if (klim .ne. 1 .and. klimold .eq. 2)  go to 3700
      if (klim .ne. 1 .and. klimold .eq. 3)  go to 3700
c
c ----------------------------------------------------------------------
c klim = 1    klimold = 1
c ----------------------------------------------------------------------
c
 3200 xfact = 1.0
      if (volume .lt. volwant .and. voladj .gt. vadjmin)
     .vadjmin = voladj
      if (volume .gt. volwant .and. voladj .lt. vadjmax)
     .vadjmax = voladj
      deladj  = ovoladj-voladj
      if (deladj .eq. 0.0)  go to 9000
      deriv   = (oldvol-volume)/deladj
      if (deriv .eq. 0.0)  go to 9000
      voladj  = voladj - (volume-volwant)/deriv
      if (vadjmin .lt. -1.0e20)  go to 3220
      if ( voladj .lt. vadjmin)  voladj  = 0.5 * (voladj+vadjmin)
 3220 if (vadjmax .gt. 1.0e30 )  go to 6000
      if ( voladj .gt. vadjmax)  vadjmax = 0.5 * (voladj+vadjmax)
      go to 6000
c
c ----------------------------------------------------------------------
c klim = 1    klimold = 0
c ----------------------------------------------------------------------
c
 3300 xfact = 1.0
      if (volume .lt. volwant .and. voladj .gt. vadjmin)
     .vadjmin = voladj
      if (volume .gt. volwant .and. voladj .lt. vadjmax)
     .vadjmax = voladj
      voladj  = voladj - SIGN (volnudge, (volume-volwant))
      go to 6000
c
c ----------------------------------------------------------------------
c klim = 1  and  (klimold = 2  or  klimold = 3)
c ----------------------------------------------------------------------
c
 3400 xfact = 1.0
      if (volume .lt. volwant .and. voladj .gt. vadjmin)
     .vadjmin = voladj
      if (volume .gt. volwant .and. voladj .lt. vadjmax)
     .vadjmax = voladj
      if (vadjmin .lt. -1.0e20)  go to 3410
      voladj  = voladj-volnudge
      go to 6000
 3410 voladj  = 0.5 * (vadjmin+vadjmax)
      go to 6000
c
c ----------------------------------------------------------------------
c klim # 1  and  klimold = 1
c ----------------------------------------------------------------------
c
 3500 if (voladj .lt. vadjmax) vadjmax = voladj
      voladj = 0.5 * (vadjmin+vadjmax)
      go to 6000
c
c ----------------------------------------------------------------------
c klim # 1  and  klimold = 0
c ----------------------------------------------------------------------
c
 3600 if (voladj .lt. vadjmax) vadjmax = voladj
      xfact  = xfact*2.0
      voladj = voladj-xfact*volnudge
      go to 6000
c
c ----------------------------------------------------------------------
c klim # 1  and  (klimold = 2  or  klimold = 3)
c ----------------------------------------------------------------------
c
 3700 if ( voladj .lt. vadjmax)  vadjmax = voladj
      if (vadjmin .lt. -1.0e20)  go to 3720
      voladj = 0.5 * (vadjmin+vadjmax)
      go to 6000
 3720 xfact  = xfact*2.0
      voladj = voladj-xfact*volnudge
      go to 6000
c
c ==============================================================================
c we want plasma to limit on the right
c ==============================================================================
c
 4000 if (klim .eq. 3 .and. klimold .eq. 3)  go to 4200
      if (klim .eq. 3 .and. klimold .eq. 0)  go to 4300
      if (klim .eq. 3 .and. klimold .eq. 1)  go to 4400
      if (klim .eq. 3 .and. klimold .eq. 2)  go to 4400
      if (klim .ne. 3 .and. klimold .eq. 3)  go to 4500
      if (klim .ne. 3 .and. klimold .eq. 0)  go to 4600
      if (klim .ne. 3 .and. klimold .eq. 1)  go to 4700
      if (klim .ne. 3 .and. klimold .eq. 2)  go to 4700
      if (klim .ne. 3 .and. klimold .eq. 3)  go to 4700
c
c ----------------------------------------------------------------------
c klim = 3    klimold = 3
c ----------------------------------------------------------------------
c
 4200 xfact = 1.0
      if (volume .gt. volwant .and. voladj .gt. vadjmin)
     .vadjmin = voladj
      if (volume .lt. volwant .and. voladj .lt. vadjmax)
     .vadjmax = voladj
      deladj  = ovoladj-voladj
      if (deladj .eq. 0.0)  go to 9000
      deriv   = (oldvol-volume)/deladj
      if (deriv .eq. 0.0)  go to 9000
      voladj  = voladj - (volume-volwant)/deriv
      if (vadjmin .lt. -1.0e20)  go to 4220
      if ( voladj .lt. vadjmin)  voladj  = 0.5 * (voladj+vadjmin)
 4220 if (vadjmax .gt. +1.0e30)  go to 6000
      if ( voladj .gt. vadjmax)  vadjmax = 0.5 * (voladj+vadjmax)
      go to 6000
c
c ----------------------------------------------------------------------
c klim = 3    klimold = 0
c ----------------------------------------------------------------------
c
 4300 xfact = 1.0
      if (volume .gt. volwant .and. voladj .gt. vadjmin)
     .vadjmin = voladj
      if (volume .lt. volwant .and. voladj .lt. vadjmax)
     .vadjmax = voladj
      voladj  = voladj + SIGN (volnudge, (volume-volwant))
      go to 6000
c
c ----------------------------------------------------------------------
c klim = 3  and  (klimold = 1  or  klimold = 2)
c ----------------------------------------------------------------------
c
 4400 xfact = 1.0
      if ( volume .gt. volwant .and. voladj .gt. vadjmin)
     .vadjmin = voladj
      if ( volume .lt. volwant .and. voladj .lt. vadjmax)
     .vadjmax = voladj
      if (vadjmax .lt. 1.0e20)  go to 4410
      voladj  = voladj + volnudge
      go to 6000
c
 4410 voladj = 0.5 * (vadjmin + vadjmax)
      go to 6000
c
c ----------------------------------------------------------------------
c klim # 3  and  klimold = 3
c ----------------------------------------------------------------------
c
 4500 if (voladj .gt. vadjmin) vadjmin = voladj
      voladj = 0.5 * (vadjmin+vadjmax)
      go to 6000
c
c ----------------------------------------------------------------------
c klim # 3  and  klimold = 0
c ----------------------------------------------------------------------
c
 4600 if (voladj .gt. vadjmin)  vadjmin = voladj
      xfact  = xfact*2.0
      voladj = voladj-xfact*volnudge
      go to 6000
c
c ----------------------------------------------------------------------
c klim # 3  and  (klimold = 1  or  klimold = 2)
c ----------------------------------------------------------------------
c
 4700 if (voladj .gt. vadjmin) vadjmin = voladj
      if (vadjmax .gt. 1.0e20)  go to 4720
      voladj = 0.5 * (vadjmin+vadjmax)
      go to 6000
c
 4720 xfact  = xfact*2.0
      voladj = voladj+xfact*volnudge
      go to 6000
c
c ==============================================================================
c
 6000 ovoladj = holdadj
      oldvol = volume
      klimold = klim
      if (vadjmin .lt. -1.0e+20)  return
      if (vadjmax .gt.  1.0e+20)  return
      xfact = 1.0
      return
c
 9000 write  (nout, 8000)
 8000 format (' FATAL ERROR in width control subroutine')
      write  (nout, 8010)  ovoladj, oldvol, holdadj, volume
 8010 format (' ovoladj,oldvol,holdadj,volume',4e12.3)
      call STOP ('subroutine GETADJV: unspecified problem', 17)
c
      end

      subroutine getfcur (iounit, ierr, kstart, kend, iavfcoil)
c
c
c ----------------------------------------------------------------------
c --- subroutine calculates the f coil currents,curfcoil,in a least
c --- squares sense. Psi loops, magnetic probes
c --- and experimental fcoil currents are used to determine
c --- the calculated fcoil currents. The effect of vessel and
c --- ecoil currents,both assumed known from the input file inone,
c --- can be accounted for on option.
c --- input(argument list):
c   iounit            fortran unit number for diagnostic output,
c                     set to 0 to suppress all diagnostic output.
c    kstart
c    kend                to avoid going over whole grid
c    iavfcoil           if = 1 average the new fcoil cur. with the old
c                          = 0 do not average,
c                              return the new fcoil currents
c
c --- input (through INCLUDE files):
c --- INCLUDE file param
c   kstore            size of temporary work vectors xdum,ydum,zdum,wdum
c                     (no other parameters from param.i are used)
c
c --- INCLUDE file mhdpar
c   nfcoil           these are the parameters set in mhdpar.i
c   nsilop           which are used in this subroutine
c   magpr2
c   nwh
c   nesum
c   nvessel
c   nwork
c
c --- common block mhdcom (INCLUDE file mhdpar):
c   pcurrent(nwh)  toroidal plasma filament current (amps)
c   ifitpsi       =1 use psi loop values in fcoil cur determination
c                 =0 do not use psi loop values in fcoil cur determination
c   ifitprob      =1 use mag probe values in fcoil cur determination
c                 =0 do not use mag probe values in fcoil cur determination
c   iecurr        =1 include ecoil current effects
c                 =0 exclude ecoil current effects
c   ivessel       =1 include vessel current effects
c                 =0 exclude vessel current effects
c   fitfcoil      =1 include experimental fcoil currents
c                 =0 exclude experimental fcoil currents from fit
c    mfcinv       =0 set up and calc inverse of design matrix
c                 =1 inverse of design matrix is allready set up
c
c --- from common blocks associated with Green's table (INCLUDE file mhdcom):
c --- these are the various inductive coupling factors used in this subroutine
c   rsilfc(nsilop,nfcoil)         fcoil to psi loop
c   rmp2fc(magpr2,nfcoil)         fcoil to magnetic probe
c   rsilec(nsilop,nesum)          ecoil to psi loop
c   rmp2ec(magpr2,nesum)          ecoil to magnetic probe
c   rsilvs(nsilop,nvessel)        vessel to psi loop
c   rmp2vs(magpr2,nvessel)        vessel to magnetic probe
c   rsilpc(nsilop,nwh)            plasma to psi loop
c   rmp2pc(magpr2,nwh)            plasma to magnetic probe
c
c --- common block mhdbctdn(INCLUDE file mhdbctdn)
c   psiloop(nsilop)        psi loop values at the current time
c                          (not used if ifitpsi = 0)
c   probeval(magpr2)       probe values at the current time
c                          (not used if ifitprob = 0)
c   ecoilcur(nesum)        e coil current at the current time
c                          (not used if iecurr = 0)
c   vescur(nvessel)        vessel current at the current time
c                          (not used if ivessel = 0)
c   expfcoil(nfcoil)       experimental fcoil currents,amps
c                          a value of exactly 0.0 for the j'th value of
c                          psiloop(j),probeval(j),ecoilcur(j),expfcoil(j) or
c                          vescur(j) means the corresponding quantity
c                          will not be used.
c  fwtpsilp(nsilop)
c  fwtmprbe(1..magpri)
c  fwtfcoil(1..nfcoil)
c  fwtecoil(1...nesum)
c  fwtvescr(1..nvessel)       the weighting vectors ( = 1/sigma)
c                             fwtecoil and fwtvescr are currently not used
c  isgngren                   set to -1 if  del-star psi =  u0*r*jphi
c                             set to +1                  = -u0*r*jphi
c
c --- common block bicube (INCLUDE file bicube)
c   wnoperm           temporary storage area
c --- common block storage (INCLUDE file storage)
c   xdum              temporary storage
c   ydum
c   zdum
c   wdum
c
c common block IMSL (INCLUDE file imsl)
c imslmd                an unecessary complication because error
c                       control (useing uertst) for IMSL errors
c                       was made local to ONETWO.
c --- output
c --- through argument list
c   ierr            =0 if no error occurred
c                    .gt. 0 if a fatal error ocurred
c --- through common block mhdbcdtn,INCLUDE file mhdbcdtn:
c   curfcoil(1..nfcoil)    the required fcoil currents (amps)
c
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE mhdpar                                
      USE mhdcom
      USE bicube
      USE mhdbcdtn
      USE replace_imsl,            ONLY : my_lginf
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'mhdbcdtn.i'
c      include 'mhdcom.i'
c      include 'bicube.i'
      include 'storage.i'
      include 'imsl.i'
c
      dimension    desgnm(maxobsrv,nfcoil),rhs(maxobsrv),idum(maxobsrv),
     .             singval(kstore)
c      equivalence (desgnm(1,1),wnoperm(1))
      equivalence (rhs(1),xdum(1))
      equivalence (singval(1),ydum(1))
      equivalence (idum(1),zdum(1))
c
      if (mfcinv .eq. 1)  go to 650
c
c --- do a couple of cross checks
c
      ierr    = 0
      idesize = maxobsrv * nfcoil
      if (idesize .gt. nwork) then
        ierr = 1
        if (iounit .ne. 0)  write (iounit, 10)  nwork, idesize
   10   format ('  subroutine GETFCUR reports nwork is too small' /
     .          '  nwork =', i5, ' required = ', i5               /
     .          '  program must be recompiled')
      end if
      if (2*maxobsrv .gt. kstore) then
        ierr = 1
        if (iounit .ne. 0)  write (iounit, 20) kstore, maxobsrv
   20   format (' subroutine GETFCUR reports kstore is too small' /
     .          '   kstore =', i5, ' required = ', i5             /
     .          '   program must be recompiled')
      end if
      if (ierr .eq. 1)  return
c
c --- set up the design matrix for the fcoil current problem:
c             desgnm(i,j)*curfcoil(j) = dpsiloop(i) (sum on j)
c             desgnm(l,j)*curfcoil(j) = dprobe(l)   (sum on j)
c --- where dpsiloop is the deficit in psiloop values which has to be
c --- made up by fcoil currents and dprobe is the deficit in mag probe
c --- values, also to be curtailed by f coil currents.
c
      tol0 = 1.0e-10      ! CRAY only
      tol  = 0.0
      itry = 0
   30 do nk=1,nfcoil
        nj = 0
        if (ifitpsi .eq. 1) then
c
c --- fcoil current to psi loop
c
          do m=1,nsilop
            if (psiloop(m) .ne. 0.0) then
              nj = nj+1
              desgnm(nj,nk) = isgngren*rsilfc(m,nk)*fwtpsilp(m)
            end if
          end do
        end if
        if (ifitprob .eq. 1) then
c
c --- fcoil current to probe contribution:
c
          do m=1,magpr2
            if (probeval(m) .ne. 0.0) then
              nj = nj+1
              desgnm(nj,nk) = rmp2fc(m,nk)*fwtmprbe(m)
            end if
          end do
        end if
        if (fitfcur .eq. 1) then
c
c --- experimental fcoil currents are included in the fit
c
          do m=1,nfcoil
            if (expfcoil(m) .ne. 0.0) then
              nj            = nj + 1
              desgnm(nj,nk) = 0.0
              if (nk .eq. m)
     .        desgnm(nj,nk) = fwtfcoil(m)
            end if
          end do
        end if
      end do
c
      neqn = nj    ! neqn equations, nfcoil unknowns
c
c --- check to see if we have enough data
c
      if (neqn .lt. nfcoil) then
        ierr = 1
        if (iounit .ne. 0)
     .  write  (iounit, 550)  neqn, nfcoil
  550   format (' FATAL ERROR from subroutine GETFCUR'                 /
     .          ' number of equations =',i5, '  number of unknowns =',
     .            i5 / ' underdetermined system is not allowed'        /
     .          ' fcoil currents must be unique')
        return
      end if
c
c --- get the generalized inverse of the (neqn,nfcoil) matrix using
c --- singular value decomposition:
c --- (recall that the least squares solution of y = A*x is given by
c --- x = (At*A)**-1*At*y,where At is A transpose. The generalized
c --- inverse of matrix A is, by definition, (At*A)**-1*At, hence that's
c --- all we need.)
c
      ier    = 0
      imslmd = 'lginf'
      call my_lginf (desgnm,maxobsrv,neqn,nfcoil,tol,fcoilinv,
     .            nfcoil,singval,wdum,ier)
c
c --- make up for ill conditioning of desgnm if necessary
c
      if (ier  .ne. 0) then
        if (itry .gt. 2)  go to 5000
        itry = itry + 1
        tol  = tol0 * 100.0   ! may make sense only on CRAY
        go to 30              ! necessary because lginf destroyed desgnm
      end if
      mfcinv = 1              ! set flag to indicate fcoilinv exists
c
c --- fcoilinv(nfcoil,neqn) is now the generalized inverse of desgnm
c --- we use this inverse below to get the least squares solution.
c --- note that fcoilinv is recalculated only when mfcinv = 0
c --- next get the rhs:
c --- enter here if inverse of design matrix is known.
c
  650 nj = 0
      if (ifitpsi .eq. 1) then
        do m=1,nsilop
          sum = 0.0
          if (psiloop(m) .ne. 0.0) then
            nj = nj+1
c
c --- get contribution to psi loop m from all sources (other than fcoil)
c --- a)plasma current:
c
            do 600 k=kstart,kend
  600       sum = sum+isgngren*rsilpc(m,k)*pcurrent(k)
            if (iecurr .ne. 0.0) then
c
c --- b)from ecoils:
c
              do 700 k=1,nesum
  700         sum = sum+isgngren*rsilec(m,k)*ecoilcur(k)
            end if
            if (ivessel .ne. 0) then
c
c --- c)from vessel currents:
c
              do k=1,nvessel
                sum = sum + isgngren * rsilvs(m,k) * vescur(k)
              end do
            end if
            rhs(nj) = fwtpsilp(m)*(psiloop(m)-sum) ! here is the deficit
          end if
        end do
      end if
      if (ifitprob .eq. 1) then
        do m=1,magpr2
          sum = 0.0
          if (probeval(m) .ne. 0.0) then
            nj = nj + 1
c
c --- a) plasma current contribution
c
            do k=kstart,kend
              sum = sum + rmp2pc(m,k) * pcurrent(k)
            end do
            if (iecurr .ne. 0) then
c
c --- b) ecoil contribution
c
              do k=1,nesum
                sum = sum+rmp2ec(m,k)*ecoilcur(k)
              end do
            end if
c
            if (ivessel .ne. 0) then
c
c --- c) vessel contribution
c
              do 1300 k=1,nvessel
 1300         sum = sum+rmp2vs(m,k)*vescur(k)
            end if
            rhs(nj) = fwtmprbe(m)*(probeval(m)-sum)      ! probe deficit
          end if
        end do
      end if
      if (fitfcur .eq. 1) then
c
c --- experimental fcoil current contributions
c
      do m=1,nfcoil
         if (expfcoil(m) .ne. 0.0) then
           nj = nj+1
           rhs(nj) = expfcoil(m)*fwtfcoil(m)
         end if
      end do
c
      end if
c
c --- now solve the overdetermined system (neqn>nfcoil) for the fcoil currents
c
      do 1500 j=1,nfcoil
         wdum(j) = curfcoil(j)     ! save old one for average below
         curfcoil(j) = 0.0
         do 1500 i=1,neqn
 1500       curfcoil(j) = curfcoil(j)+fcoilinv(j,i)*rhs(i)
         if (iavfcoil .eq. 1) then
            do 1600 j=1,nfcoil
 1600       curfcoil(j) = 0.5 * (curfcoil(j)+wdum(j))
         end if
c
      return
c
 5000 ierr = ier
      if (iounit .ne. 0)
     .  write  (iounit, 5200)  ier
 5200   format (' subroutine GETFCUR reports ERROR in least squares' /
     .          ' solution for fcoil currents, ier from lginf =', i5)
      return
c
      end

      subroutine getg (ii, elongax)
c
c
c ----------------------------------------------------------------------
c  the flux surface average of the toroidal component of Ampere's law,
c  which in MKS units reads (eq. 2.4-1):
c             curden = (1/u0)*d(gcap*hcap*rho*bp0)/d(rho)/(hcap*rho),
c  is integrated over rho to get
c             gcap = flux surface average of (grad(rho)*rmajor/r)**2.
c  this subroutine assumes hcap*rho*gcap*bp0 evaluated at rho => 0 = 0
c  subroutine also calculates the rbp(j) profile. (the initial
c  condition for Faraday's Law, each time a new transport cycle is started).
c --- input
c  through argument list:
c  ii                          not used right now
c
c --- through INCLUDE files
c  INCLUDE file param
c  kj                        size of vectors indexed by nj below
c
c  INCLUDE file etc:
c  bp(1..nj)                 bp(rho) is defined by bp = (1/rmajor)*dpsi/drho
c
c  INCLUDE file geom:
c  fcap(1...nj)              f(psi)/f(psilim)
c  hcap(1..nj)               fcap(psi)/<R0**2/R**2>
c
c  INCLUDE file constnts
c  u0,twopi
c
c  INCLUDE file mesh:
c  r(1..nj)                  rho grid,meters
c
c  INCLUDE file soln:
c  rbp(nj)                   fcap*gcap*hcap*rho*bp at rho = rhomax,
c                            only the nj'th value is used (for
c                            normalization . rbp(nj) is set by
c                            the total current boundary condition
c  curden(1...n)             <jtoroidal*R0/R>
c
c  INCLUDE file storage
c  xdum
c  ydum(1...nj)              work area,min length nj
c
c --- output
c  INCLUDE file geom:
c  gcap(1..nj)                 <(grad rho *R0/R)**2>
c
c  INCLUDE file soln:
c  rbp(1...nj)               fcap*gcap*hcap*rho*bp
c
c ----------------------------------------------------------------------
c
      USE param
      USE soln
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE constnts
      USE etc
      implicit  integer (i-n), real*8 (a-h, o-z)

c      include 'etc.i'

      include 'storage.i'
c
      dimension qconsist(kj)
c
c ----------------------------------------------------------------------
c calculate integral of hcap*rho*curden and store in ydum
c (Thus ydum is the toroidal current divided by two pi, since the
c element of cross sectional area da is given by da = 2 * pi * hcap*rho,
c note that curden =<jtoroidal*R0/R> and rbp(nj) = fcap(nj)*hcap(nj)*
c gcap(nj)*rho(nj)*bp0(nj) and fcap(nj) = 1.0 )
c ----------------------------------------------------------------------
c



      call curpro(nj)   ! here nj is used to indicate argument > 0

      u0otpi  = u0/twopi
      do 10 j=1,nj
   10 xdum(j) = hcap(j)*r(j)*curden(j)
      call trap2(r,xdum,ydum,nj)
      xnorm   = rbp(nj)/(u0otpi*ydum(nj))

c
c ----------------------------------------------------------------------
c now normalize ydum so it is the toroidal current
c ----------------------------------------------------------------------
c
      do 15 j=1,nj
   15 ydum(j) = xnorm*ydum(j)
c
c ----------------------------------------------------------------------
c calculate gcap; note that bp here is in tesla
c ----------------------------------------------------------------------
c
      do  j=2,nj
         gcap(j) = ydum(j)*u0otpi/(r(j)*bp(j)*hcap(j))
         if(rbp_save(nj) .ne. 0.0)
     .         gcap (j) = rbp_save(j)*1.e-6/(fcap(j)*hcap(j)*r(j)*bp(j))
      enddo
      call cubicextrp (gcap(2),gcap(3),gcap(4),r(2),r(3),r(4),gcap(1),2)
   
c
c ----------------------------------------------------------------------
c define the correct rbp as used in transport (in MKS units however)
c ----------------------------------------------------------------------
c

         DO  j=1,nj
            rbp(j) = ydum(j)*fcap(j)*u0otpi
            if(rbp_save(nj) .ne. 0.0)
     .            rbp(j) = rbp_save(j)*1.e-6
         ENDDO
 
c
c --------------------------------------------------------------------------
c     as a consistency check we calculate q again, now using the formula
c       q = twopi * (btor / rmajor) * (gcap * hcap * rho**2) / (I(rho))
c     qconsist(1) = q(1)  not calculated here
c
      do j=2,nj
        qconsist(j) = btor * gcap(j) * hcap(j) * r(j)**2
     .                     / ydum(j) / rmajor / u0otpi
      end do



c
      return
c
      end

      subroutine getmhdbc (ltest_code)
c

c
c ----------------------------------------------------------------------
c --- subroutine gets the  boundary conditions for psi
c --- at time point t = time.These are the values of psi
c --- at the flux loops and/or probes  for the coils option,
c --- and the values of psi
c --- on the boundary of the MHD grid for the no coils option.
c --- note that the array flxeqbcd serves double duty. for the no coils
c --- option flxeqbcd contains the input boundary values of psi.
c --- for the coils option flxeqbcd contains the psi loop values if
c --- ifitpsi = 1,otherwise flxeqbcd is not used.
c --- Sets up f coil current vector,curfcoil(1,..nfcoil),if
c --- the fcoil currents are assumed known at the boundary condition
c --- times. If  the fcoil currents are to be calculated consistently
c --- with the psi loop and probe data then the experimental fcoil current
c --- values are loaded into expfcoil(j).(curfcoil will be defined in
c --- subroutine GETFCUR for this latter case).
c --- the ecoil and vessel currents are also interpolated in time
c --- if respective options are set.
c --- get the toroidal current,pcurmhdt,the toroidal b field,btormhdt
c --- and the loop voltage,vlopmhdt,all at time point t = time.
c --- (these quantities may or may not be used,but are loaded at any rate).
c --- finally the magnetic probe and partial rogowski values are determined
c --- appropriate weighting vectors used in least squares fitting
c --- of fcoil currents (see subroutine GETFCUR) are constructed here.
c
c --- if this is a test case,(ltest_code = 1),then the first time we call
c --- this routine we take the boundary conditions found and write them
c --- into the other time slots so that the MHD boundary conditions
c --- become constant in time.
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE io
      USE solcon
      USE mhdpar 
      USE extra   
      USE mhdcom
      USE shapctr
      USE etc
      USE mhdbcdtn
      implicit  integer (i-n), real*8 (a-h, o-z)

c      include 'mhdbcdtn.i'
c      include 'etc.i'

c      include 'shapctr.i'
c      include 'mhdcom.i'
c
      data ltest_pass /0/
c
c --- find first time pt in bc cond. vector which is gt current time:
c --- should allow .ge. in lin 4050 for final time
c
      do i=1,mxtbcmhd
        i2 = i
        if (timeqbcd(i) .gt. time)  go to 2200
      end do
c
c ----------------------------------------------------------------------
c at the final time,the parameter time may be slightly greater than
c timmax (because subroutine CHEKDT adds a tolerance,time_tol to time.)
c Hence we do the following check:
c (timeqbcd(itbcmhd) gives the largest time at which boundary
c conditions are available.)NOTE the parameter time_tol used below should
c have the same value as in subroutine CHEKDT.
c We do not pass it by INCLUDE file however.
c ----------------------------------------------------------------------
c
c      time_tol = 1.0e-09   defined in solcon
      if (time-timeqbcd(itbcmhd) .lt. 10.0*time_tol) then
        i2 = itbcmhd
        go to 2200
      end if
      write  (ncrt, 8000)
      write  (nout, 8000)
 8000 format (' FATAL ERROR in subroutine GETMHDBC:  ',
     .        ' time > timeqbcd(i)  for  i = 1,mxtbcmhd ' )
      call STOP ('subroutine GETMHDBC: problem #1', 18)
c
 2200 i1 = i2-1
      if (mhdmode .eq. 'no coils') then
        call load_psincbcd
      else if (mhdmode .eq. 'coils') then
c
        do i=1,nsilop
          if (i1 .gt. 0) then
            psiloop(i) = ((flxeqbcd(i,i2)-flxeqbcd(i,i1))/(timeqbcd(i2)
     .               -timeqbcd(i1)))*(time-timeqbcd(i1))+flxeqbcd(i,i1)
          end if
          psiloop(i)  = -psiloop(i)    ! psi convention in ONETWO
****      fwtpsilp(i) = ABS (errpsilp*psiloop(i))
****      if (fwtpsilp(i) .ne. 0.0)  fwtpsilp(i) = 1.0 / fwtpsilp(i)
        end do
c
      fwtmax = 0.0
      do 2510 i=1,nsilop
 2510   fwtmax = MAX (fwtmax, ABS (psiloop(i)))
      if (fwtmax .ne. 0.0)  fwtmax = 1.0/(errpsilp*fwtmax)
      do i=1,nsilop
        fwtpsilp(i) = 0.0
        if (psiloop(i) .ne. 0.0)  fwtpsilp(i) = fwtmax
      end do
      if (fixfcoil .eq. 1) then   ! fcoil currents fixed at input values
        do i=1,nfcoil
          if (i1 .gt. 0) then
          curfcoil(i) = ((fcoilcur(i,i2)-fcoilcur(i,i1))/(timeqbcd(i2)
     .                -timeqbcd(i1)))*(time-timeqbcd(i1))+fcoilcur(i,i1)
          else
             curfcoil(i)=fcoilcur(i,i2)
          end if
          curfcoil(i) = mhdmultp*curfcoil(i)
        end do
      end if
c
c --- load exptl values of fcoil currents
c
          do i=1,nfcoil
            if (i1 .gt. 0) then
            expfcoil(i) = ((fcoilcur(i,i2)-fcoilcur(i,i1))/(timeqbcd(i2)
     .                -timeqbcd(i1)))*(time-timeqbcd(i1))+fcoilcur(i,i1)
            else
               expfcoil(i)=fcoilcur(i,i2)
            end if
            expfcoil(i) = mhdmultp*expfcoil(i)
****        fwtfcoil(i) = ABS (errfcoil*expfcoil(i))
****        if (fwtfcoil(i) .ne. 0.0)  fwtfcoil(i) = 1.0 / fwtfcoil(i)
          end do
      fwtmax = 0.0
      do 2360 i=1,nfcoil
 2360   fwtmax = MAX (fwtmax, ABS (expfcoil(i)))
      if (fwtmax .ne. 0.0)  fwtmax = 1.0 / (errfcoil*fwtmax)
      do i=1,nfcoil
        fwtfcoil(i) = 0.0
        if (expfcoil(i) .ne. 0.0)  fwtfcoil(i) = fwtmax
      end do
        if (iecurr .eq. 1) then
          do i=1,nesum
            if (i1 .gt. 0) then
              ecoilcur(i) = ((ecurrt(i,i2)-ecurrt(i,i1))/(timeqbcd(i2)
     .                -timeqbcd(i1)))*(time-timeqbcd(i1))+ecurrt(i,i1)
            else
              ecoilcur(i)=ecurrt(i,i2)
            end if
            ecoilcur(i) = mhdmultp*ecoilcur(i)
****        fwtecoil(i) = ABS (errecoil*ecoilcur(i))
****        if (fwtecoil(i) .ne. 0.0)  fwtecoil(i) = 1.0 / fwtecoil(i)
          end do
      fwtmax = 0.0
      do 2660 i=1,nesum
 2660   fwtmax = MAX (fwtmax, ABS (ecoilcur(i)))
      if (fwtmax .ne. 0.0)  fwtmax = 1.0/(errecoil*fwtmax)
      do i=1,nesum
        fwtecoil(i) = 0.0
        if (ecoilcur(i) .ne. 0.0)  fwtecoil(i) = fwtmax
      end do
        end if
        if (ivessel .eq. 1) then
          do i=1,nvessel
            if (i1 .gt. 0) then
            vescur(i) = ((vescurrt(i,i2)-vescurrt(i,i1))/(timeqbcd(i2)
     .              -timeqbcd(i1)))*(time-timeqbcd(i1))+vescurrt(i,i1)
            else
               vescur(i)=vescurrt(i,i2)
            end if
            vescur(i) = vescur(i)*mhdmultp
****        fwtvescr(i) = ABS (errvescr*vescur(i))
****        if (fwtvescr(i) .ne. 0.0)  fwtvescr(i) = 1.0 / fwtvescr(i)
          end do
      fwtmax = 0.0
      do i=1,nvessel
        fwtmax = MAX (fwtmax, ABS (vescur(i)))
      end do
      if (fwtmax .ne. 0.0)  fwtmax = 1.0 / (errvescr*fwtmax)
      do i=1,nvessel
        fwtvescr(i) = 0.0
        if (vescur(i) .ne. 0.0)  fwtvescr(i) = fwtmax
      end do
      end if
c
c --- load experimental probe values
c
          do 2800 i=1,magpr2
            if (i1 .gt. 0) then
          probeval(i) = ((expmp2(i,i2)-expmp2(i,i1))/(timeqbcd(i2)
     .           -timeqbcd(i1)))*(time-timeqbcd(i1))+expmp2(i,i1)
          else
             probeval(i)=expmp2(i,i2)
          end if
          fwtmprbe(i) = ABS (errmprbe*probeval(i))
          probeval(i) = probeval(i)*mhdmultp
          if (fwtmprbe(i) .ne. 0.0)  fwtmprbe(i) = 1.0 / fwtmprbe(i)
 2800     continue
      fwtmax = 0.0
      do 2860 i=1,magpr2
 2860   fwtmax = MAX (fwtmax, ABS (probeval(i)))
      if (fwtmax .ne. 0.0)  fwtmax = 1.0 / (errmprbe*fwtmax)
      do 2870 i=1,magpr2
        fwtmprbe(i) = 0.0
        if (probeval(i) .ne. 0.0)  fwtmprbe(i) = fwtmax
 2870 continue
c
c --- other MHD data
c
          if (i1 .gt. 0) then
          pcurmhdt = ((pcurmhd(i2)-pcurmhd(i1))/(timeqbcd(i2)
     .           -timeqbcd(i1)))*(time-timeqbcd(i1))+pcurmhd(i1)
          pcurmhdt = pcurmhdt*mhdmultp
          btormhdt = ((btormhd(i2)-btormhd(i1))/(timeqbcd(i2)
     .           -timeqbcd(i1)))*(time-timeqbcd(i1))+btormhd(i1)
          btormhdt = btormhdt*mhdmultp
          vlopmhdt = ((vloopmhd(i2)-vloopmhd(i1))/(timeqbcd(i2)
     .           -timeqbcd(i1)))*(time-timeqbcd(i1))+vloopmhd(i1)
          vlopmhdt = vlopmhdt*mhdmultp
          else
             pcurmhdt =pcurmhd(i2)
             pcurmhdt = pcurmhdt*mhdmultp
             btormhdt = btormhd(i2)
             btormhdt = btormhdt*mhdmultp
             vlopmhdt = vloopmhd(i2)
             vlopmhdt = vlopmhdt*mhdmultp
          end if
      else
        write  (nout, 2550) mhdmode
        write  (ncrt, 2550) mhdmode
 2550   format (' subroutine GETMHDBC detects incorrect setting' /
     .          ' of mhdmode = ', a                              /
     .          ' program must stop'                              )
        call STOP ('subroutine GETMHDBC: problem #2', 19)
      end if
      if (i1 .gt. 0) then
        volwant =
     .    ((volaray(i2)-volaray(i1))/(timeqbcd(i2)-timeqbcd(i1)))*
     .    (time-timeqbcd(i1))+volaray(i1)
      else
        volwant = volaray(i2)
      end if
c
      if (ltest_code .eq. 0)  return
      if (ltest_pass .ne. 0)  return
      if (mhdmode .ne. 'coils') then
        write  (nout, 3000)
        write  (ncrt, 3000)
 3000   format (' ERROR: mhdmode = coils must be set' /
     .              8x, 'if ltest_code = 1.  ONETWO must stop.')
        call STOP ('subroutine GETMHDBC: problem #3', 20)
      end if
c
c --- do the following only if ltest_code=1 and ltest_pass =0
c --- load the MHD boundary condition vectors with the current
c --- value. further interpolations in time,done above,will
c --- thus always return the same values for testing purposes.
c
      ltest_pass = 1
      do j=1,mxtbcmhd
            do i=1,nsilop
                 flxeqbcd(i,j) = -psiloop(i)
            end do
            do i=1,nfcoil
                 fcoilcur(i,j) = expfcoil(i)/mhdmultp
            end do
            if (iecurr .eq. 1) then
                  do i=1,nesum
                       ecurrt(i,j) = ecoilcur(i)/mhdmultp
                  end do
            end if
            if (ivessel .eq. 1) then
                  do i=1,nvessel
                       vescurrt(i,j) = vescur(i)/mhdmultp
                  end do
            end if
            do i=1,magpr2
                 expmp2(i,j) = probeval(i)/mhdmultp
            end do
            pcurmhd(j) = pcurmhdt/mhdmultp
            btormhd(j) = btormhdt/mhdmultp
            vloopmhd(j) = vlopmhdt/mhdmultp
      end do
c
      return
c
      end

      subroutine getrmaj1 (iounit, nj)
c
c
c ----------------------------------------------------------------------
c --- subroutine gets values of rmajor at the magnetic axis elevation,zma,
c --- that correspond to the rho values given in r,on the outboard side
c --- of the plasma. this routine is not set up to handle cases where
c --- zma is not the magnetic axis!!!!
c --- the rho values are not used directly. instead,since the psir
c --- vector gives values of psi that are in one to one correspondence
c --- with the rho values,we use the psir grid.
c
c --- INPUT(argument list):
c  iounit           unit number for error message (if any)
c  nj               size of r,psir
c
c --- INPUT INCLUDE files:
c
c --- INCLUDE file param.i:
c  kj               size of various vectors (used in rhog.i)
c
c --- INCLUDE file mhdpar.i:
c  nw               size of vector rmhdgrid
c  nh               size of vector zmhdgrid
c
c --- INCLUDE file bicube.i:
c  cspln            stores bicubic spline coefficients
c  wnoperm          work storage
c  pds(6)           work storage
c
c --- INCLUDE file mhdcom.i:
c  psi(nw,nh)       psi array (on rmhdgrid,zmhdgrid)
c  rma              major radius to magnetic axis,m
c  zma              elevation of magnetic axis,m
c
c --- INCLUDE file mhdgrid.i:
c  rmhdgrid(i)      i = 1,2,..nw horizontal grid,m
c  zmhdgrid(j)      j = 1,2,..nh vertical grid,m
c
c --- INCLUDE file rhog.i:
c  psir(j)          psi values corresponding to r(j) (volt-sec)
c
c --- INCLUDE file contour.i:
c  rplasmax         max radial extent of palsma,m
c
c --- OUTPUT:
c
c --- INCLUDE file rhog.i:
c  rmajorvec(j)     j = 1,2,..nj major radius vector,m, corresponding to
c                   rho values in r and psi values in psir
c  bprmaj(j)        true (i.e., not flux surf. avg.) magnitude of
c                   poloidal b field at rmajorvec(j) (which corresponds
c                   to r(j) and psir(j))
c  btotrmaj(j)      j=1,2,..nj same for total field
c
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE contour
      USE mhdpar  
      USE mhdgrid  
      USE rhog
      USE mhdcom
      USE geom,   ONLY : fcap
      USE machin,  ONLY : rmajor,btor
      USE bicube
      USE replace_imsl,       ONLY : my_ibcccu,my_dbcevl1
      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'imsl.i'
c
      logical psiminonax
c
      rsearchup = MIN (rplasmax+0.01,rmhdgrid(nw))
      imslmd    = 'getrmaj1'
c
c --- define the bicubic spline coefficients,cspln
c
      ier = 0
      imslmd='4140c208'
      call my_ibcccu (psi,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                wnoperm,ier)
      if (ier .ne. 0)  go to 20
c
c --- set up to get rmajorvec
c
      psiminonax = .false.
      if (psir(1) .lt. psir(nj))  psiminonax = .true.
      rmajorvec(1) = rma
      bprmaj   (1) = 0.0
      btotrmaj(1)  = rmajor*btor/(fcap(1)*rma)
      do j=2,nj
        itry = 0
        rlow = rmajorvec(j-1)
        rup = rsearchup
   10   rtry = 0.5 * (rlow+rup)
        itry = itry+1
        ier1 = 0
        call my_dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,nw,rtry,zma,
     .                   pds,ier1,6)
        if (ier1 .ne. 0)  go to 20
        if (psiminonax) then
          if (pds(1) .lt. psir(j)) then
              rlow = MAX (rtry,rma)
          else
              rup = MIN (rtry,rsearchup)
          end if
        else
          if (pds(1) .ge. psir(j)) then
              rlow = MAX (rtry,rma)
          else
              rup = MIN (rtry,rsearchup)
          end if
        end if
        if (itry .gt. 150)  go to 20
        if (rlow .lt. rma .or. rlow .gt. rsearchup)  go to 20
        if (rup .lt. rma .or. rup .gt. rsearchup)  go to 20
        if (ABS (rlow-rup) .gt. 1.0e-6)  go to 10
        rmajorvec(j) = rtry
        bprmaj   (j) =
     .  SQRT (pds(2)**2 + pds(3)**2) / rtry ! b-poloidal at rmajorvec(j)
        btotrmaj(j)  = SQRT((rmajor*btor/rmajorvec(j)*fcap(j))**2 +
     .                 bprmaj(j)**2)
      end do

 
c
      return
   20 write  (iounit, 1)
    1 format (' subroutine GETRMAJ1 reports:'     /
     .        ' could not determine major radius' /
     .        ' diagnostics are:')
      if (ier  .ne. 0)  write (iounit, 2) ier
      if (ier1 .ne. 0)  write (iounit, 3) ier1
    2 format (' IBCCCU error, ier =', i5)
    3 format (' DBCEVL error, ier =', i5)
      write  (iounit, 4)  itry, j
    4 format ('  itry,j =',2(2x,i5))
      write  (iounit, 5)  rma,zma,rsearchup,rtry
    5 format ('  rma,zma,rsearchup,rtry =' / 4(2x,1pe12.6))
      write  (iounit, 6)  pds(1),psir(j)
    6 format ('  pds(1),psir(j) =', 2(2x,1pe12.6))
      write  (iounit, 7)
    7 format (' ONETWO is terminated')
      call STOP ('subroutine GETRMAJ1: unspecified problem', 21)
c
      end

      subroutine getrmaj (zelev, nf, ierr, rmajor,
     .                    psirmaj, isigncur, psiedge)
c
c
c --- at elevation zelev,get a uniformly-spaced major radius array, rmajor,
c --- starting at the inside edge of the plasma and working to the
c --- outside edge. for each value of rmajor(j),get a corresponding
c --- value of psi,stored in psirmaj(j). if zelev is above the
c --- top of the plasma or below the bottom of the plasma return
c --- with error code ierr = 1.0
c --- ROUTINE ASSUMES BICUBIC SPLINE ARRAY CSPLN IS SET UP AND READY FOR USE
c --- BY DBCEVL1 !
c ------------------------------------------------------------------ HSJ
c
      USE contour
      USE mhdpar
      USE mhdgrid    
      USE mhdcom
      USE bicube
      USE replace_imsl,                      ONLY : my_dbcevl1 
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'mhdpar.i'
c      include 'mhdgrid.i'
c      include 'mhdcom.i'
c      include 'bicube.i'
c      include 'contour.i'
      include 'imsl.i'
c
      dimension  rmajor(nf), psirmaj(nf)
c
      imslmd = 'getrmaj'
      ierr   = 0
      if ((zelev .ge. zplasmax) .or. (zelev .le. zplasmin))  ierr = 1
      if (ierr .eq. 1)  return
      toler  = 1.0e-06
c
c --- first get inside radius
c
      ra = rmhdgrid(1)
      rb = rma
c
  100 rg = 0.5 * (ra+rb)
      call chkinout (rg, zelev, psirzp, ind, isigncur, psiedge)
      if (ABS (ra-rb) .lt. toler)  go to 120
      if (ind .eq. 0) then
        ra = rg
      else
        rb = rg
      end if
      go to 100
c
  120 rmajor(1) = rg
c
c --- next get outside radius
c
      ra = rma
      rb = rmhdgrid(nw)
c
  140 rg = 0.5 * (ra+rb)
      call chkinout (rg, zelev, psirzp, ind, isigncur, psiedge)
      if (ABS (ra-rg) .lt. toler)  go to 130
      if (ind .eq. 0) then
        rb = rg
      else
        ra = rg
      end if
      go to 140
c
  130 rmajor(nf) = rg
      dr = (rmajor(nf)-rmajor(1))/(nf-1)
      do 150 j=2,nf-1
        rmajor(j) = rmajor(j-1)+dr
        call my_dbcevl1 (rmhdgrid, nw, zmhdgrid, nh, cspln, nw, 
     .                    rmajor(j),zelev, pds, ier, 1)
  150 psirmaj(j)  = pds(1)
      psirmaj(1)  = psibdry
      psirmaj(nf) = psibdry
c
      return
c
      end

      subroutine headeq (iunit)
c
      USE solcon
      USE param
      USE io
      USE soln2d
      USE etc
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subroutine is to generate the top line of each page of MHD output
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'etc.i'
c      include 'io.i'
c      include 'solcon.i'
c      include 'soln2d.i'
c
      if (iunit .ne. nout)  return
      write (iunit, 10)  time, n, ieq, itre
   10 format ('1time = ',3pf14.2,' msec,  time point=', i4,
     . '   equilibrium point = ', i3,
     . '   equilibrium iteration number = ', i2,
     .  '  equilibrium')
      return
c
      end

      subroutine initmhd (rmhdgrid,nw,zmhdgrid,nh,zero,rmagax,
     .                    zmagax,rscale,zscale,rmajor,u0,totcur,cspln,
     .                    icsplnd1,icsplnd2,icsplnd3,pds,psi1d,
     .          wnoperm,nlimiter,xlimiter,ylimiter,rcontr,zcontr,maxpts,
     .          torcurnt,nj,iounit,psi,pprim,ffprim,rm2inv,psir,press)

c
c ----------------------------------------------------------------------
c --- Note this routine assumes input is in MKS units!
c --- routine is used for startup of equilibrium calculations.
c --- It assumes that psi has the form
c      psi(i,j) = a * EXP (-(rmhdgrid(i)-rmagax)**2/rscale
c                     -(zmhdgrid(j)-zmagax)**2/zscale)
c --- where a >0.0
c --- and uses this psi together with the input current and
c --- kinetic profiles to get pprim and ffprim starting guesses.
c --- input
c  rmhdgrid(1,..nw)
c  zmhdgrid(1...nh)        gives the MHD grid (in meters)
c  zero(1..nw,1...nh)      the limiter scrapeoff array
c  xlimiter(1...nlimiter)   the limiter points vectors
c  ylimiter(1...nlimiter)
c  rmagax
c  zmagax
c  rscale
c  zscale                  shape control factors,see above definition
c                          of psi. by setting these factors appropriately
c                          limited control over the initial plasma
c                          shape is possible. if any of these factors
c                          is not set on input appropriate values will
c                          be determined by this subroutine,see coding
c                          below.
c  rmajor                  factor used to calculate <rmajor**2/R**2>
c  u0                      4.0 * pi*1.0e-07
c  totcur                  total toroidal current,amps
c  torcurnt(1...nj)        current vector <jtoroidal*rmajor/R>,amps/m**2
c  iounit                  fortran unit no. for error diagnostics,set to
c                          zero to supress diagnostic output.
c --- temporary storage required:
c  cspln(icsplnd1,icsplnd2,icsplnd3)
c                          cspln is storage array for cubic spline coefs.
c                          icsplnd1 must be exactly =2
c                          icsplnd2 must be exactly equal to nw
c                          icsplnd3 must be .ge.  2*nh
c  pds(6)                  storage vector
c  wnoperm(1)              storage vector of min length 2*nw*nh+2*max(nw,nh)
c  rcontr(1...maxpts)
c  zcontr(1...maxpts)     storage for contour points found by
c                        subroutines bound and cntour.
c                        maxpts is the dimension of the vectors rcontr,
c                        zcontr as defined in the calling program. If
c                        the number of points found exceeds maxpts,an
c                        error exit is taken.
c --- output
c  psi(1...nw,1...nh)      psi array defined over rmhdgrid,zmhdgrid
c                          in units of volt sec per radian.
c  psi1d(kk)            kk = 1,..nw*nh  is a 1 d bersion of psi
c                       psi1d is arranged so that psi(i,j) = psi1d(kk)
c                       where kk = (i-1)*nh+j (note that this is NOT
c                       the standard fortran convention,since j indexes
c                       fastest here. This is required for compatability
c                       with the Green's table and the scrapeoff package!
c  pprim(1..nj)
c  ffprim(1..nj)             the dp/dpsi and f*df/dpsi vectors defined
c                            on the rho grid (rho grid is implicit in
c                            this routine)
c  rm2inv(1....nj)           flux surface average,<rmajor**2/R**2>
c  psir(1...nj)              psi on the rho grid,volt-sec/radian
c  press(1...nj)             thermal pressure  nt/m**2
c ------------------------------------------------------------------ HSJ
c
      USE replace_imsl,       ONLY : my_ibcccu
c
      
      implicit  integer (i-n), real*8 (a-h, o-z)
      include 'imsl.i'

      dimension rmhdgrid(nw),zmhdgrid(nh),zero(nw*nh),pds(*),
     .  torcurnt(nj),psi(nw,nh),pprim(nj),ffprim(nj),rm2inv(nj),
     .  cspln(icsplnd1,icsplnd2,icsplnd3),psir(nj),
     .  xlimiter(nlimiter),ylimiter(nlimiter),psi1d(*),
     .  rcontr(*),zcontr(*),wnoperm(*),press(nj)
c
c --- set and check scale factors:
c
      if (rmagax .eq. 0.0)  rmagax = 0.5 * (rmhdgrid(nw)+rmhdgrid(1))
      if (zmagax .eq. 0.0)  zmagax = 0.5 * (zmhdgrid(nh)+zmhdgrid(1))
      if (rscale .eq. 0.0)  rscale = 0.25*(rmhdgrid(nw)-rmhdgrid(1))**2
      if (zscale .eq. 0.0)  zscale = 0.25*(zmhdgrid(nh)-zmhdgrid(1))**2
      if (rmagax .lt. rmhdgrid(3) .or. rmagax .gt. rmhdgrid(nw-2))
     .  go to 2000
      if (zmagax .lt. zmhdgrid(3) .or. zmagax .gt. zmhdgrid(nh-2))
     .  go to 2000
c
c --- get non normalized psi first:
c
      do 10 j=1,nh
        dzscsq = (zmhdgrid(j)-zmagax)**2/zscale
        do 10 i=1,nw
          drscsq = (rmhdgrid(i)-rmagax)**2/rscale
           psi(i,j) = EXP (-drscsq-dzscsq)
        kk = (i-1)*nh+j
        psi1d(kk) = psi(i,j)
   10 continue
c
c --- next find the plasma boundary with above assumed psi:
c
        xlmin  = rmhdgrid(2)
        ix     = 0
        nerr   = 0
        isignn = 1
        limfag = 2
        radold = xlimiter(nlimiter+1)
        if (totcur .lt. 0.0)  isignn = -1
        call bound (psi1d,nw,nh,nw*nh,zero,rmhdgrid,zmhdgrid,rmagax,
     .        zmagax,xlmin,ix,nlimiter,xlimiter,ylimiter,nerr,limfag,
     .        radold,maxpts,rcontr,zcontr,ncontr,dpsi,rplasmin,rplasmax,
     .        zplasmin,zplasmax,rzmin,rzmax,zrmin,zrmax,psimin)
c
c --- now get estimate for magnitude of psi based on support set for
c --- j toroidal. we do this by evaluating del-star and
c --- setting result equal to j-toroidal:
c
      sum  = 0.0
      drdz = (rmhdgrid(2)-rmhdgrid(1))*(zmhdgrid(2)-zmhdgrid(1))
      do 20 j=1,nh
        if (zmhdgrid(j) .lt. zplasmin)  go to 20
        if (zmhdgrid(j) .gt. zplasmax)  go to 20
        deltaz = zmhdgrid(j)-zmagax
        do 21 i=1,nw
          if (psi(i,j) .lt. psimin)  go to 21
          if (rmhdgrid(i) .lt. rplasmin)  go to 21
          if (rmhdgrid(i) .gt. rplasmax)  go to 21
          deltar  = rmhdgrid(i)-rmagax
          d1rdpsi =  -2.0*psi(i,j)*deltar/rscale
          d2rdpsi = (-2.0/rscale)*(psi(i,j)+deltar*d1rdpsi)
          d1zdpsi =  -2.0*psi(i,j)*deltaz/zscale
          d2zdpsi = (-2.0/zscale)*(psi(i,j)+deltaz*d1zdpsi)
          sum = sum+(d2rdpsi-d1rdpsi/rmhdgrid(i)+d2zdpsi)/rmhdgrid(i)
   21 continue
   20 continue
      if (sum .eq. 0.0)  go to 1000
      anorml = -totcur*u0/(sum*drdz)
c
c --- scale psi to approximately correct total current:
c
      do 30 j=1,nh
        do 30 i=1,nw
          psi(i,j) = anorml*psi(i,j)
          kk       = (i-1)*nh+j
   30   psi1d(kk)  = psi(i,j)
c
c --- define the psi grid corresponding to rho:
c
      psimin  = anorml*psimin
      psimax  = anorml
      dpsi    = (psimin-psimax)/(nj-1)
      do 40 j=1,nj
   40 psir(j) = psimax+(j-1)*dpsi
c
c --- next we need to get the flux surface average <rmajor**2/R**2>
c --- a bicubic spline representation of psi is used here so that
c --- subroutine CNTOUR can do its thing:
c
      imslmd='4480c208'
      call my_ibcccu (psi,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,wnoperm,ier)

c
c --- trace each contour of psir(j) and do the flux surface integrals:
c
      iauto     = 1
      dang      = 5.0
      dang      = 0.1 ! 999,
      bperr     = 0.05
      bperr     = 1.e-3 ! 999,
      dr        = SQRT (drdz)
      dz        = dr
      arcl      = 0.02     ! 0.02 meters
      arcl      = 0.001    ! 0.001 meters ! 999,
      rm2inv(1) = (rmajor/rmagax)**2
      delta_psi=psir(nj-1)-psir(nj)
      do 50 j=2,nj
      psivl = psir(j)
      dx0   = 0.0
      dy0   = 0.0
      call cntour (rmagax,zmagax,PSIVL,remin,remax,zemin,zemax,
     .             zrmin,zrmax,rzmin,rzmax,dang,arcl,bperr,dx0,dy0,
     .             xlimiter(nlimiter+1),xlimiter(nlimiter+2),
     .             ylimiter(nlimiter+1),ylimiter(nlimiter+2),
     .             iauto,iautoc,rcontr,zcontr,ncontr,rmhdgrid,nw,
     .             zmhdgrid,nh,cspln,icsplnd1,icsplnd3,iounit,maxpts,
     .             ierr,wnoperm,0,delta_psi)
c
c --- now do the integrals
c
        sum  = 0.0
        sum1 = 0.0
        do i=1,ncontr-1
          bpol1 = wnoperm(i)
          bpol2 = wnoperm(i+1)
          dl    = SQRT ((rcontr(i+1)-rcontr(i))**2
     .                + (zcontr(i+1)-zcontr(i))**2)
          sum1  = sum1 + 0.5 * dl * (bpol1+bpol2)/(bpol1*bpol2)
          sum   = sum  + 0.5 * dl
     .                       * (rcontr(i)**2*bpol1+rcontr(i+1)**2*bpol2)
     .                       / (rcontr(i)**2*bpol1*rcontr(i+1)**2*bpol2)
        end do
c
   50 rm2inv(j) = rmajor**2*sum/sum1
c
c --- get the plasma pressure
c --- (neglects beams, fusion since these are not set at this point):
c
      call pressr (1, 1)
c
c --- get dp/dpsi
c

      call difydx (psir, press, pprim, nj)

c
c --- the current profile specified in inone gives values of
c --- <jtor*rmajor/r> as a function of rho. we now use this
c --- input profile to calculate ffprim,where
c ---   <jtor*rmajor/R>*rmajor = rmajor**2*pprim+
c ---                         +ffprim*<(rmajor/R)**2>/u0
c
      do j=1,nj
        ffprim(j) = u0*rmajor*(torcurnt(j)-rmajor*pprim(j))/rm2inv(j)
        ffprim(j) = -u0*rmajor*(torcurnt(j)+rmajor*pprim(j))/rm2inv(j) ! 8/25/06 HSJ
      end do
      return
c
 1000 if (iounit .ne. 0) then
        write (iounit,1010)
 1010   format (' subroutine INITMHD has detected an ERROR'           /
     .          ' mhdgrid, nlimiter points or psi guess incompatible' /
     .          ' program must stop'                                  )
        call STOP ('subroutine INITMHD: problem #1', 22)
      end if
c
 2000 if (iounit .ne. 0) then
        write  (iounit, 2010)
 2010   format (' ERROR: subroutine INITMHD has detected an ERROR'  /
     .              8x, 'RMAGAX and/or ZMAGAX are outside MHD grid' /
     .              8x, 'program must stop'                         )
        call STOP ('subroutine INITMHD: problem #2', 23)
      end if
      return
c
      end

      subroutine integ (x, y, n, sum)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- integrate y with respect to x from x(1) to x(n)
c
      dimension  x(*), y(*)
c
      sum = 0.0
      do 10 j=2,n
   10 sum = sum+(x(j)-x(j-1))*(y(j)+y(j-1))
      sum = 0.5 * sum
      return
c
      end

      subroutine limitpts (rs,zs,thets,re,ze,thete,rtan,rsl,zsl,thetsl,
     .                     rel,zel,thetel,nchord,ylim,xlim,limitr)
c
      USE constnts
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- given a straight line defined by the two points (rs,thets,zs)
c --- and (re,thete,ze),extend this line until it intersects the
c --- limiter and return the (assumed two) intersection points as
c --- (rsl,thetsl,zsl) and (rel,thetel,zel).
c --- Note that it is assumed that the interior of the limiter forms
c --- a convex set,at least with regard to the chords that are to be
c --- found! At present this routine is set up only for chords that are
c --- in poloidal or horizontal planes. A general 3-d inclination of the
c --- chord is not allowed.
c --- if rtan   is -1.0e30 then the chord is assumed to be in a poloidal plane
c --- if rtan .ne. -1.0e30 the chord is assumed to be in a horizontal plane
c --- with rtan = tangency radius.
c --- if the required intersection with the limiter is not found
c --- rsl,thetsl,zsl and rel,thetel,zel are returned as 0
c ------------------------------------------------------------------ HSJ
c
      dimension rs(*),zs(*),thets(*),re(*),ze(*),thete(*),
     .          rsl(*),zsl(*),thetsl(*),rel(*),zel(*),thetel(*),
     .          ylim(*),xlim(*),rtan(*),rr(4)
      dimension rsav(2),zsav(2)
c
c      include 'constnts.i'
c
      do 10 j=1,nchord
        rsl(j)    = rs(j)
        zsl(j)    = zs(j)
        thetsl(j) = thets(j)
        rel(j)    = re(j)
        zel(j)    = ze(j)
   10 thetel(j)   = thete(j)
c
c --- loop over chords
c
   20 do 100 j=1,nchord
        rsav(1) = 0.0
        rsav(2) = 0.0
        zsav(1) = 0.0
        zsav(2) = 0.0
c
c --- is chord in a poloidal or horizontal plane?
c --- if tangency radius is set chord is in a horizontal plane (by assumption)
c
      if (rtan(j) .ne. -1.0e30)  go to 800
      thet  = thetsl(j)*twopi/360.0
      thet1 = thetel(j)*twopi/360.0
      if (ABS (thet-thet1) .lt. twopi/360.0)  go to 120
c
c --- is chord in a horizontal plane?
c
      if (zs(j)-ze(j) .eq. 0.0)  go to 800
c
c --- chord has general 3-d inclination.  add this coding when necessary
c
      call STOP ('subroutine LIMITPTS: 3-D problem', 24)
c
c --- chord is in a horizontal plane.
c     get radius of limiter at z = zs(j)  ( = ze(j))
c
  800 ict = 0
      do 820 k=1,limitr-1
      if ((zs(j)-ylim(k))*(zs(j)-ylim(k+1)) .gt. 0.0)  go to 820
      if (zs(j) .ne. ylim(k))  go to 810
      ict = ict+1
      rsav(ict) = xlim(k)
      zsav(ict) = zs(j)
      go to 820
  810 if (zs(j) .eq. ylim(k+1))  go to 820
      ict = ict+1
      zsav(ict) = zs(j)
      if (xlim(k) .eq. xlim(k+1))  rsav(ict) = xlim(k)
      if (xlim(k) .eq. xlim(k+1))  go to 820
      a = (ylim(k+1)-ylim(k))/(xlim(k+1)-xlim(k))
      bincp = ylim(k+1)-a*xlim(k+1)
      rsav(ict) = (zs(j)-bincp)/a
  820 continue
c
      if (ict .ne. 2)
     .  call STOP ('subroutine LIMITPTS: ICT problem', 25)
      rtj = rtan(j)
      if (rtj .ne. -1.0e30)  go to 840
c
c --- get tangency radius of chord (chord may not have been specified
c --- using the rtan option)
c
      x1 = rs(j) * COS (thets(j)*twopi/360.0)
      y1 = rs(j) * SIN (thets(j)*twopi/360.0)
      x2 = re(j) * COS (thete(j)*twopi/360.0)
      y2 = re(j) * SIN (thete(j)*twopi/360.0)
c
c --- note that dx and dy can not be zero here because this would imply
c --- that the chord lies in a single poloidal plane
c
      a     = (y2-y1)/(x2-x1)
      bincp = y1-a*x1
c
c --- equation of chord is y = a*x+bincp
c --- equation of tangency radius is y = -x/a
c --- solve these two equations simultaneousl for (x,y)
c
      y1      = bincp / (1.0 + a*a)
      x1      = -a * y1
      rtj     = SQRT (x1*x1+y1*y1)
      rtan(j) = rtj
c
c --- order points in major radius
c
  840 rsvmi = MIN (rsav(1),rsav(2))
      rsvmx = MAX (rsav(1),rsav(2))
      if (j .eq. 1)  rr(1) = rsvmi
      if (j .eq. 1)  rr(4) = rsvmx
      zel(j) = zsav(1)
      zsl(j) = zsav(1)
      if (rtj .le. rsvmi)  go to 830
      if (rtj .ge. rsvmx)
     .  call STOP ('subroutine LIMITPTS: RSVMX problem', 26)
      rsl(j)    = rsvmx
      rel(j)    = rsvmx
      thetsl(j) = 0.0
      thetel(j) = 2.0 * ACOS (rtj/rsvmx) * 360.0 / twopi
      go to 100
  830 rsl(j)    = rsvmx
      rel(j)    = rsvmi
      thetsl(j) = 0.0
      ang       =      ASIN (rtj/rsl(j))
      ang2      = pi - ASIN (rtj/rel(j))
      thetel(j) = (pi-ang-ang2)*360.0/twopi
      go to 100
c
c --- chord is in a poloidal plane.
c --- z = a*r+bincp, or r = a*z+bincp depending on range
c
  120 ict = 0
      dzc = ze(j)-zs(j)
      drc = re(j)-rs(j)
c
c --- is chord vertical or horizontal?
c
      if (drc .eq. 0.0)  go to 140
      if (dzc .eq. 0.0)  go to 110
c
c --- chord is inclined
c
      thet = ATAN2 (dzc, drc)
      if (thet .lt. 0.0)  thet = twopi-thet
      if (((piov4 .lt. thet) .and. (thet .lt. tpiov4)) .or.
     .   ((fpiov4 .lt. thet) .and. (thet .lt. spiov4)))  go to 140
c
c --- get z = a*r+bincp
c
  110 a     = dzc/drc
      bincp = ze(j)-a*re(j)
      iflg  = 0
      go to 160
c
c --- get r = a*z+bincp
c
  140 a     = drc/dzc
      bincp = re(j)-a*ze(j)
      iflg  = 1
c
c --- loop over limiter (assumes xlim(limitr) = xlim(1) and
c ---                            ylim(limitr) = ylim(1))
c
  160 k = 0
      do 170 kk=1,limitr-1
      k  = k + 1
      dr = xlim(k+1)-xlim(k)
      dz = ylim(k+1)-ylim(k)
      if (dr .eq. 0.0)  go to 180
      if (dz .eq. 0.0)  go to 190
c
c --- limiter segement is inclined
c
      thet = ATAN2 (dz,dr)
      if (thet .lt. 0.0)  thet = twopi-thet
      if ((( piov4 .lt. thet) .and. (thet .lt. tpiov4)) .or.
     .    ((fpiov4 .lt. thet) .and. (thet .lt. spiov4)))  go to 200
c
c --- get z = a1*r+bincp1
c
      a1     = dz/dr
      bincp1 = ylim(k+1)-a1*xlim(k+1)
      ilflg  = 0
      go to 210
c
c --- get r = a1*z+bincp1
c
  200 a1     = dr/dz
      bincp1 = xlim(k+1)-a1*ylim(k+1)
      ilflg  = 1
c
c --- solve for intersection of limiter segement with chord
c
  210 if (iflg .ne. 0)  go to 220
      if (ilflg .ne. 0)  go to 230
c
c --- iflg = 0, ilflg = 0       z = a1*r+bincp1 and z = a*r+bincp
c --- if chords are parallel-no intersection
c
      denom = a1 - a
      if (denom .eq. 0.0)  go to 170
      rslv  = (bincp-bincp1)/denom
      zslv  = a1*rslv+bincp1
      go to 240
c
c --- iflg = 0, ilflg = 1       r = a1*z+bincp1 and z = a*r+bincp
c
  230 denom = 1.0 - a1*a
c
c --- denom = 0 for parallel (not perpendicular) chords
c
      if (denom .eq. 0.0)  go to 170
      rslv = (a1*bincp+bincp1)/denom
      zslv =  a*rslv+bincp
      go to 240
  220 if (ilflg .ne. 0)  go to 250
c
c --- iflg = 1, ilflg = 0       r = a*z+bincp and z = a1*r+bincp1
c
      denom = 1 - a * a1
      if (denom .eq. 0.0)  go to 170
      rslv  = (a*bincp1+bincp)/denom
      zslv  = a1*rslv+bincp1
      go to 240
c
c --- iflg = 1, ilflg = 1       r = a*z+bincp and r = a1*z+bincp1
c
  250 denom = a1 - a
      if (denom .eq. 0.0)  go to 170
      zslv  = (bincp-bincp1)/denom
      rslv  = a1*zslv+bincp1
      go to 240
c
c --- limiter chord segement is vertical
c
  180 rslv = xlim(k)
c
c --- if chord is also vertical
c
      if (drc .eq. 0.0)  go to 170
      if (iflg .eq. 0)  zslv = a*rslv+bincp
      if (iflg .eq. 1)  zslv = (rslv-bincp)/a
      go to 240
c
c --- limiter chord is horizontal
c
  190 zslv = ylim(k)
c
c --- if co2 chord is also horizontal,no intersection assumed
c
      if (dzc .eq. 0.0)  go to 170
      if (iflg .eq. 0.0)  rslv = (zslv-bincp)/a
      if (iflg .eq. 1  )  rslv = a*zslv+bincp
c
c --- is intersection in range?
c
  240 if ((zslv-ylim(k))*(zslv-ylim(k+1)) .gt. 0.0)  go to 170
      if ((rslv-xlim(k))*(rslv-xlim(k+1)) .gt. 0.0)  go to 170
      ict       = ict+1
      rsav(ict) = rslv
      zsav(ict) = zslv
      if (ict .eq. 2)  go to 280
c
c --- it is assumed that a chord intersects the limiter exactly twice.
c --- don't search next limiter interval. this prevents counting a single
c --- intersection twice due to the fact that the chord passes through
c --- a limiter coordinate exactly
c
      k = k + 1
      if (k .ge. limitr)  go to 280
  170 continue
c
c --- assign start and end of chord arbitrarily
c
  280 rsl(j) = rsav(1)
      rel(j) = rsav(2)
      zsl(j) = zsav(1)
      zel(j) = zsav(2)
  100 continue
      return
c
      end

      subroutine linear (fun1, fun2, fun, x1, x2, x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c linear interpolation to determine fun (x)
c ----------------------------------------------------------------------
c
      fun = ((x-x1)*fun2+(x2-x)*fun1)/(x2-x1)
      return
c
      end

      subroutine load_psincbcd
c
c
c ----------------------------------------------------------------------
c --- load psincbcd(j),j = 1,2..2*(nw+nh-2) starting at the lower
c --- left hand corner(j = 1)  and proceeding clockwise.
c --- vertical boundaries
c ----------------------------------------------------------------------
c
      USE mhdpar    
      USE mhdcom
      USE mhdbcdtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'mhdpar.i'
c      include 'mhdbcdtn.i'
c      include 'mhdcom.i'
c
          kk = nh+nw-2
          do j=1,nh
             kk = kk+1
             psincbcd(j ) = psi( 1,   j  )
             psincbcd(kk) = psi(nw,nh-j+1)
          end do
c
c --- horizontal boundaries:
c
          kk = nh
          kkk = 2*nh+nw-2
          kkk = 2*nh+nw-2
          do i=2,nw-1
              kk = kk+1
              psincbcd(kk) = psi(i,nh)
              kkk = kkk+1
              psincbcd(kkk) = psi(nw-i+1,1)
          end do
c
      return
c
      end

      subroutine magax (psi,x,y,nw,nh,isignn,iknowax,iounit,ispln,zero,
     .                  xax,yax,cspln,n2cspln,nh2,psymx,elax)

c
c ----------------------------------------------------------------------
c --- find the magnetic axis. by definition the magnetic axis is that
c --- point inside the limiter contour at which psi takes on an extremal
c --- value,either a maximum (positive current) or a minimum (negative
c --- current). convergence to the magnetic axis is considered accomplished
c --- when gradient of psi < gradtol (gradtol is defined locally below)
c --- subroutine terminates ONETWO run if magnetic axis can not be found!
c --- input
c  psi(nw,nh)       the array of psi values
c  x(nw)
c  y(nh)             the MHD grid vectors
c  isignn        =  1 search for maximum in psi
c                = -1 search for minimum in psi
c  iknowax         =1 good intial guess for magnetic axis is input
c                     in xax(1),yax(1)
c                  =0 good initial guess is not available, find intitial
c                     guess for magnetic axis by searching the MHD grid
c  iounit             fortran unit number for diagnostic output
c                     iounit = 0 suppresses output,but this routine
c                     may terminate the ONETWO run so diagnostic
c                     messages should be printed out.
c  ispln             =1 bicubic spline array cspln is input
c                    =0 bicubic spline array cspln must be calculated
c  zero(nwh)         the limiter contour indicator vector
c  xax
c  yax       intial guess of magnetic axis location if known,see iknowax
c  cspln(n2cspln,nw,nh2)   bicubic spline coeff. of psi if known
c                          if not known,set ispln = 0 to calculate it here
c  wk(ndwk)      temporary starge vector of minimum length
c                ndwk = 2*nh*nw+2*max(nh,nw) ( required if ispln=0 )
c                 alloacated below
c
c --- output
c  iknowax = 1       indicates axis is now known
c  iknowax = 0       if the axis was not found iknowax=0 is returned
c  ispln = 1         indicates that cspln is now set
c  cspln           the bicubic spline coefficients of psi
c  psymx           the value of psi on the magnetic axis
c  elax            the elongation of the plasma at the magnetic axis
c  xax
c  yax             final converged coordiantes of magnetic axis
c
c ------------------------------------------------------------------ HSJ
c

c
      USE replace_imsl,              ONLY : my_ibcccu,my_dbcevl1


      implicit  integer (i-n), real*8 (a-h, o-z)
      dimension psi(nw,nh),x(nw),y(nh)
      dimension xax(*),yax(*),psymx(*),zero(*)
      dimension cspln(n2cspln,nw,nh2),pds(6)
      real *8 ,dimension(:),allocatable :: wk
      include 'imsl.i'


      if( .not. allocated(wk))then
             ndwk = 2*nh*nw+2*max(nh,nw)
              allocate (wk(ndwk),STAT = istat)
              if(istat .ne. 0)
     .          call allocate_error("wk, magax",0,istat)
      endif


c
      gradtol = 1.0e-12
      itry = 1
      i1set = 0
      if (ispln .eq. 1)  go to 5
c
c --- set up bicubic spline coefficients:
c
         imslmd='5000c208'
      call my_ibcccu (psi,x,nw,y,nh,cspln,nw,wk,ier)
      deallocate(wk)
      ispln = 1
    5 if (iknowax .eq. 1)  go to 130   ! good initial guess available
c
c --- first find approximate location of magnetic axis by
c --- searching for max/min psi on grid inside limiter box.
c
      iknowax = 0
      i1      = 0
      j1      = 0
      if (isignn .eq. -1) then    ! find minimum psi
        psimin = +1.0d100
        do 10 j=2,nh-1
        do 10 i=2,nw-1
        kk = (i-1)*nh+j
        zzsum = zero(kk-nh-1)+zero(kk-nh)+zero(kk-nh+1)+zero(kk-1)
     .        + zero(kk)+zero(kk+1)+zero(kk-1+nh)+zero(kk+nh)
     .        + zero(kk+nh+1)
        if (zzsum .lt. 8.5)  go to 10
        if (psi(i,j) .gt. psi(i-1,j-1))  go to 10
        if (psi(i,j) .gt. psi(i-1,j  ))  go to 10
        if (psi(i,j) .gt. psi(i-1,j+1))  go to 10
        if (psi(i,j) .gt. psi(i,j-1  ))  go to 10
        if (psi(i,j) .gt. psi(i,j+1  ))  go to 10
        if (psi(i,j) .gt. psi(i+1,j-1))  go to 10
        if (psi(i,j) .gt. psi(i+1,j  ))  go to 10
        if (psi(i,j) .gt. psi(i+1,j+1))  go to 10
c
c --- the above picks smallest value of psi on the x,y grid,inside the plasma,
c --- such that the eight nearest neighbors of pt. kk are also inside the plasma
c --- thus a magnetic axis within dx or dy of limiter contour as seen by array
c --- zero will not be found using this method.
c
        psimin = MIN (psimin,psi(i,j))
        if (psimin .ne. psi(i,j))  go to 10
        i1 = i
        j1 = j
   10   continue
c
        if (i1*j1 .eq. 0) then
          if (iounit .ne. 0)  write (iounit, 125)
          iknowax = 0
          return
        end if
      else if (isignn .eq. 1) then    ! find max psi
        psimax = -1.0d100
        do 20 j=2,nh-1
        do 20 i=2,nw-1
        kk = (i-1)*nh+j
        zzsum = zero(kk-nh-1)+zero(kk-nh)+zero(kk-nh+1)+zero(kk-1)+
     .  zero(kk)+zero(kk+1)+zero(kk-1+nh)+zero(kk+nh)+zero(kk+nh+1)
        if (zzsum .lt. 8.5)  go to 20
        if (psi(i,j) .lt. psi(i-1,j-1))  go to 20
        if (psi(i,j) .lt. psi(i-1,j))  go to 20
        if (psi(i,j) .lt. psi(i-1,j+1))  go to 20
        if (psi(i,j) .lt. psi(i,j-1))  go to 20
        if (psi(i,j) .lt. psi(i,j+1))  go to 20
        if (psi(i,j) .lt. psi(i+1,j-1))  go to 20
        if (psi(i,j) .lt. psi(i+1,j))  go to 20
        if (psi(i,j) .lt. psi(i+1,j+1))  go to 20
c
c --- the above picks largest value of psi on the x,y grid,inside the plasma,
c --- such that the eight nearest neighbors of pt. kk are also inside the plasma
c --- thus a magnetic axis within dx or dy of limiter contour as seen by array
c --- zero will not be found using this method.
c
        psimax = MAX (psimax,psi(i,j))
        if (psimax .ne. psi(i,j))  go to 20
        i1 = i
        j1 = j
   20   continue
          if (i1*j1 .eq. 0) then
              if (iounit .ne. 0)  write (iounit, 125)
  125         format (' subroutine MAGAX reports:' /
     .                ' magnetic axis was not found')
              iknowax = 0
              return
          end if
      else
        if (iounit .ne. 0)  write (iounit, 25)
   25   format
     .  (' subroutine MAGAX reports incorrect setting of ISIGNN =',i5 /
     .   ' ONETWO must stop')
        call STOP ('subroutine MAGAX: unspecified problem', 27)
      end if
      i1set  = 1
      i1save = i1
      j1save = j1
      xax(1) = x(i1)
      yax(1) = y(j1)    ! estimate of magnetic axis on the grid
c
c --- final convergence by Newton's method. finds simultaneous
c --- solution of dpsi/dx = 0.0,dpsi/dy = 0.0
c --- maximum of 25 iterations is allowed
c
  130 psymx(1) = 0.0
c
      do j=1,25
        jax    = j
        if (xax(1) .lt. x(1) .or. xax(1) .gt. x(nw))  go to 115
        if (yax(1) .lt. y(1) .or. yax(1) .gt. y(nh))  go to 115
        xaxold = xax(1)
        yaxold = yax(1)
        call my_dbcevl1(x,nw,y,nh,cspln,nw,xax(1),yax(1),pds,ier,6)
        det    = pds(5)*pds(6)-pds(4)*pds(4)
        if (det .eq. 0.0)  go to 110     ! can't take any further steps
        xerr   = (-pds(2)*pds(6)+pds(4)*pds(3))/det
        yerr   = (-pds(5)*pds(3)+pds(2)*pds(4))/det
        xax(1) = xax(1)+xerr
        yax(1) = yax(1)+yerr
        if (pds(2)**2+pds(3)**2 .lt. gradtol)  go to 110
****    if (ABS (xax(1)-xaxold) + ABS (yax(1)-yaxold) .lt. 1.0e-10)  go to 110
        if (ABS (xerr) + ABS (yerr) .lt. 1.0e-6)  go to 110
      end do
c
c --- no convergence to magnetic axis.
c --- get a refined starting guess and try one more time.
c --- failure to find a magnetic axis is a fatal error
c
  115 if (itry .gt. 3) then
        if (iounit .ne. 0)
     .  write  (iounit, 120)  jax, itry, xax(1), xerr, yax(1), yerr, det
  120   format (' subroutine MAGAX reports no convergence to'    /
     .          ' magnetic axis after', i6, ' Newton iterations' /
     .          ' and ',i5,' outer passes'                       //
     .          ' current estimate is xax,xerr = ',2(2x,1pe12.4) /
     .          '                     yax,yerr = ',2(2x,1pe12.4) /
     .          ' the discriminant is ',1pe12.4                  //
     .          ' (negative discriminant means xax, yax',
     .          ' is tending toward a saddle point)')
        iknowax = 0
        return
      end if
c
      itry = itry + 1
      if (iknowax .eq. 1) then
        iknowax = 0    ! axis was not known sufficiently precisely
        go to 5        ! to allow Newton's method to converge
      end if           ! start over..
c                      ..get the guess by using the coarse grid search
c
c --- fine grid search localizes axis sufficiently that Newton's method
c --- can be applied (done only after the coarse grid search,for which
c --- iknowax = 0 was set)
c
      dxx = (x(2) - x(1)) / 20.0  ! fine step size, 1/20 of grid spacing
      dyy = (y(2) - y(1)) / 20.0
c
      if (i1set .eq. 1) then
        i1 = i1save
        j1 = j1save
      else
        iknowax = 0
        go to 5
      end if
c
      xs     = x(i1-1)
      ys     = y(j1-1)
      i1     = (x(i1+1)-x(i1-1))/dxx+1
      j1     = (y(j1+1)-y(j1-1))/dyy+1
      psimin = +1.0d100
      psimax = -1.0d100
      do 40 i=1,i1
        rma  = xs+(i-1)*dxx
        do 40 j=1,j1
          yma = ys + (j-1)*dyy
          call my_dbcevl1(x,nw,y,nh,cspln,nw,rma,yma,pds,ier,6)
          if (isignn .eq. -1) then
            psimin = MIN (pds(1),psimin)
            if (psimin .ne. pds(1))  go to 40
            xax(1) = rma
            yax(1) = yma
          else if (isignn .eq. 1) then
            psimax = MAX (psimax,pds(1))
            if (psimax .ne. pds(1))  go to 40
            xax(1) = rma
            yax(1) = yma
          end if
   40 continue
      go to 130                 ! try Newton's method with refined guess
  110 if (det .lt. 0.0)  go to 115  ! det<0 means converged to saddle pt
      iknowax  = 1                  ! converged to magnetic axis
      psymx(1) = pds(1)
c
c --- calculate elongation on axis based on general quadratic
c --- form of an ellipse;a*x**2+2b*x*y+c*y**2+d*x+e*y+f
c --- the coefficients a,b,c,d,e,f,are obtained from the taylor
c --- series of psi about the magnetic axis
c
      a = 0.5 * pds(5)
      b = 0.5 * pds(4)
      c = 0.5 * pds(6)
c
c --- rotate coordinates
c
      thet = 0.5 * ATAN2 (2*b,a-c)
      ap   = a * COS (thet)**2 + 2.0 * b * SIN (thet) * COS (thet)
     .     + c * SIN (thet)**2
      cp   = a * SIN (thet)**2 - 2.0 * b * SIN (thet) * COS (thet)
     .     + c * COS (thet)**2
      if (ap/cp .ge. 0.0)  elax = SQRT (ap/cp)
      return
c
      end

      subroutine maxpsi (xl1, yl1, xl2, yl2, x1, y1, x2, y2,
     .                   f1, f2, f3, f4, area, psimax, xtry, ytry)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c**                                                                   **
c**     MAIN PROGRAM:  MHD FITTING CODE                               **
c**                                                                   **
c**     SUBPROGRAM DESCRIPTION:                                       **
c**          maxpsi finds the largest psi value along the line        **
c**          segment y = alpha*x+beta joining the two points          **
c**          (xl1,yl1) and (xl2,yl2).                                 **
c**                                                                   **
c**     RECORD OF MODIFICATION:                                       **
c**          26/04/83..........first created                          **
c**          24/07/85..........revised                                **
c**                                                                   **
c***********************************************************************
c
      dx = xl2-xl1
      if (dx .eq. 0.0)  go to 100
      dy = yl2-yl1
      if (dy .eq. 0.0)  go to 100
      c = f1+f3-(f2+f4)
      if (c .eq. 0.0)  go to 100
      a = y2*(f2-f1)+y1*(f4-f3)
      b = x1*(f2-f3)+x2*(f4-f1)
      alpha = dy/dx
      secder = 2.0 * alpha*c
      if (secder .gt. 0.0)  go to 100
      beta = yl1-alpha*xl1
      xcrit = -(b*alpha+c*beta+a)/secder
      if ((xcrit .gt. xl2) .or. (xcrit .lt. xl1))  go to 100
      ycrit = alpha*xcrit+beta
      if ((ycrit .lt. y1) .or. (ycrit .gt. y2))  go to 200
      xl2 = xcrit
      yl2 = ycrit
      psip1 = -1.0e+35
      go to 110
  100 call fqlin(x1,y1,x2,y2,f1,f2,f3,f4,xl1,yl1,area,psip1)
  110 call fqlin(x1,y1,x2,y2,f1,f2,f3,f4,xl2,yl2,area,psip2)
      psimax = MAX (psip1,psip2)
      if (psimax .ne. psip1)  go to 120
      xtry = xl1
      ytry = yl1
      return
  120 xtch = xl2
      ytch = yl2
      return
  200 continue
      call STOP ('subroutine MAXPSI: unspecified problem', 97)
c
      end

      subroutine mhd
c
      USE param
      USE io
      USE solcon
      USE soln
      USE contour
      USE mhdpar     
      USE mhdgrid
      USE extra
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE flags
      USE constnts
      USE soln2d
      USE psig
      USE rhog
      USE mhdcom
      USE bicube
      USE shapctr
      USE flxav
      USE etc
      USE gpsi
      USE mhdbcdtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- MHD segment of 1-1/2-D transport code
c ----------------------------------------------------------------------
c



      include 'storage.i'

c
c ----------------------------------------------------------------------
c convert relevant quantitites to MKS:
c psimks   converts      from kgauss-cm**2 to volt-sec
c psikgaus converts back from volt-sec     to kgauss-cm**2
c ----------------------------------------------------------------------
c
      tocur   = 5.0 * rbp(nj)    ! tocur is in amps
      flim    = flim*1.0e-06
      rmajor  = 0.01*rmajor
      volfac  = 4.0 * pisq*rmajor
      rminor  = 0.01*rminor
      btor    = 1.0e-4*btor
      psibdry = psibdry*psimks
      psiaxis = psiaxis*psimks
      pmax    = psibdry
      pmin    = psiaxis
      rma     = rma*0.01
      zma     = zma*0.01
      xmagn1  = 0.01*xmagn1
      ymagn1  = 0.01*ymagn1
      xax(1)  = 0.01*xax(1)
      yax(1)  = 0.01*yax(1)
      rsep    = 0.01*rsep
      zsep    = 0.01*zsep
      rgeom   = 0.01*rgeom
      btgeom  = 1.0e-04*btgeom
      volume  = 1.0e-6*volume
      circum  = 0.01*circum
      cconst  = 1.0e-04
      call multpl1 (cxareanpsi,kpsi,cconst)
      call multpl1 (sfareanpsi,kpsi,cconst)
      cconst  = 1.0e-02
      call multpl1 (rmajavnpsi,kpsi,cconst)
      call multpl1 (rminavnpsi,kpsi,cconst)
      cconst  = 1.0e-06
      call multpl1 (psivolp,kpsi,cconst)
      cconst  = 0.01
      call multpl1 (rmhdgrid,nw,cconst)
      call multpl1 (zmhdgrid,nh,cconst)
      call multpl1 (rcontr,ncontr,cconst)
      call multpl1 (zcontr,ncontr,cconst)
      call multpl1 (rplasbdry,nplasbdry,cconst)
      call multpl1 (zplasbdry,nplasbdry,cconst)
      rplasmax = 0.01*rplasmax
      rplasmin = 0.01*rplasmin
      zplasmax = 0.01*zplasmax
      zplasmin = 0.01*zplasmin
      call multpl1 (psi,nwh,psimks)
      call multpl1 (p,nwh,psimks)
      call multpl1 (psival,npsi,psimks)
      call multpl1 (psir,nj,psimks)
      cconst = 1.0e7
      call multpl1 (pprim,nj,cconst)
      cconst = 1.0e-04
      call multpl1 (ffprim,nj,cconst)
      cconst = 1.0e-06
      call multpl1 (fpsi,npsi,cconst)
      cconst = 1.0e+4
      call multpl1 (curold, nj,cconst)    ! amps/m**2
      call multpl1 (curden,nj,cconst)
      call multpl1 (psiold,nj,psimks)
      cconst = 1.0e7
      call multpl1 (ppold,nj,cconst)
      cconst = 1.0e-4
      call multpl1 (ffpold,nj,cconst)
      cconst = 1.0e-2
      call multpl1 (rold,nj,cconst)
      cconst = 1.0e-1
      call multpl1 (press,nj,cconst)
      cconst = 1.0e-02
      call multpl1 (rho,npsi,cconst)
      rhomax = 0.01*rhomax
      rhoa   = 0.01*rhoa
      rhoa0  = 0.01*rhoa0
      drhoadt_geom = 0.01*drhoadt_geom
      cconst = 0.01
      call multpl1 (r,nj,cconst)
      call multpl1 (ra,nj-1,cconst)
      call multpl1 (dr,nj-1,cconst)
      cconst = 1.0e2
      call multpl1 (drr,nj,cconst)
      call multpl1 (rrm,nj,cconst)
      call multpl1 (rrp,nj-1,cconst)
      cconst = 1.0e-01
      call multpl1 (bp, nj, cconst)       ! change from kgauss to tesla
      cconst = 1.0e-02
      call multpl1 (rcap,nj,cconst)
      call multpl1 (rcap0,nj,cconst)
      call multpl1 (drcapdt,nj,cconst)
      cconst = 100.
      call multpl1 (rcapi,nj,cconst)
      call multpl1 (rcap0i,nj,cconst)
      call multpl1 (drcapidt,nj,cconst)
      cconst = 1.0e-02
      call multpl1 (ravg,npsi,cconst)
      cconst = 100.
      call multpl1 (ravgi,npsi,cconst)
      cconst = 1.0e-06
      call multpl1 (rbp, nj, cconst)      ! tesla-meters
c      print *,'c208,l5397 rbp=',rbp(nj-1),rbp(nj)
      cconst = 1.0e-01
      call multpl1 (pressb,nj,cconst)
      cconst = 1.0e-04
      call multpl1 (r2capi,nj,cconst)
      call multpl1 (r2capi0,nj,cconst)
      call multpl1 (dr2idt,nj,cconst)
      cconst = 1.0e-02
      call multpl1 (rmajorvec,nj,cconst)
      cconst = 1.0e-01
      call multpl1 (bprmaj,nj,cconst)        ! from kgauss to tesla
      call multpl1 (btotrmaj,nj,cconst)      ! from kgauss to tesla
c
c ----------------------------------------------------------------------
c get boundary values of psi. these are the values of
c psi on the psiloops for mhdmode = 'coils' and the values of
c psi on the boundary of the mhdgrid if mhdmode = 'no coils'
c get psiloop,probe,ecoil and vessel currents as appropriate
c get toroidal current,b field and loop voltage (these values may or may
c not be used below,depending on setting of irguess)
c get the desired volume if volume control is set.
c note that the Grad-Shafranov equation will have to be solved at the
c intial time (i.e., ieq = 0) if irguess .ge. 0, independent of ifixshap
c and at subsequent times (i.e., ieq>0) if ifixshap = 0, independent of irguess
c ----------------------------------------------------------------------
c
      if (((irguess .ge. 0 .and. ieq      .eq. 0)  .or.
     .     (ieq     .gt. 0 .and. ifixshap .eq. 0))
     .                     .and. mhdmethd .ne. 'tdem')
     .                                   call getmhdbc (ltest_code)
c
c ----------------------------------------------------------------------
c initialize parameters
c if ii = 1 the Grad-Shafranov equation will be solved.
c ----------------------------------------------------------------------
c
      ieqfail = 0
      dtime   = 0.0
      if (ieq .ne. 0)  dtime = time - timold
      IF(steady_state == 0.0)dtime = 1.e10 
c      dtime = dtime*steady_state  !HSJ 3/11/08
      ic = 0     ! ic counts number of passes required to converge pprim
      ii = 1
c
c     at startup (ieq = 0) solve the G.S. equation if irguess .ge. 0:
c
      if (ieq .eq. 0 .and.  irguess .lt. 0 .and.
     . mhdmethd .ne. 'toq')  ii = 0
c
c     at subsequent times (ieq > 0) solve the G.S. equation if ifixshap = 0
c
      if (ieq .gt. 0 .and. ifixshap .eq. 1 .and.
     . mhdmethd .ne. 'toq' )  ii = 0
c
c ----------------------------------------------------------------------
c calculate plasma pressure
c ----------------------------------------------------------------------
c
      ! HSJ 4/21/95 calculate pressure after TPORT
      if ( ieq .ne. 0)
     .  call pressr (1, 1)  ! even if no MHD calcs done (to write eqdsk)
c
c ----------------------------------------------------------------------
c calculate psi on the rho mesh,psir(j),j = 1,2..nj
c we get psir by integrating the equation bp0 = (dpsi/drho)/R0
c from the magnetic axis outward to the plasma edge. the constant
c of integration is called pmin and is the value of psi at the magnetic
c axis. don't call PSIRHO if ieq = 0 because gcap is not yet set !!!!!
c ----------------------------------------------------------------------
c
 2100 continue         ! HSJ 6/22/99 (see original label 2100 below)
      psiaxis = psir(1)
      psibdry = psir(nj)
****  call psiset (npsi, psiaxis, psibdry, psival, 0)
      call psiset (npsi, psiaxis, psibdry, psival, 1)
      if (ieq .eq. 0)  go to 2200
      call psirho (1)
c
c ----------------------------------------------------------------------
c calculate pprime and ffprime on psir grid,normalize the current to
c tocur if ic .ne. 0
c ----------------------------------------------------------------------
c

      call pfprim (ic)
c
c ----------------------------------------------------------------------
c skip subsequent (ieq .ne. 0) equilibrium calculations if ifixshap or istop = 1
c ----------------------------------------------------------------------
c
      if (ifixshap .eq. 1)  go to 3000
c
c ----------------------------------------------------------------------
c calculate equilibrium
c ----------------------------------------------------------------------
c
      write  (nitre, 8000)
 8000 format ('  entering outer iteration loop on current density')
c
c --- savcur copies pprim into ppold,ffprim into ffpold,
c --- curden into curold,psir into psiold,and r into rold,so that
c --- subroutine CONCUR can do its thing:
c
*2100 continue
 2200 call savcur
      call drive (ii, ic)

c
c ----------------------------------------------------------------------
c if ieqfail .ne. 0 then a catastrophic failure occurred (i.e., for some
c reason we could not find an equilibrium at all).  in this case we
c go directly to subroutine PREPAR, which resets everything and then we
c cut the time step in half in subroutine RUNTWO:
c ----------------------------------------------------------------------
c
      if (ieqfail .gt. 0)  go to 3010
c
c ----------------------------------------------------------------------
c update quantities on rho mesh
c ----------------------------------------------------------------------
c

      if (ii .eq. 0)  call savcur
      call rhoset (1)
      if (ii .eq. 0)  call curpro(ii)
c
c ----------------------------------------------------------------------
c output various quantities; iterate on pprime and ffprime
c if ic < itcur
c convergence is achieved when squared euclidean norm of pprim-ppold
c divided by squared euclidean norm of pprim is less than tolcur
c ----------------------------------------------------------------------
c
      call outgo
      if (ii .eq. 0)  go to 2500
      ic = ic + 1

      call pfprim (ic)
      call concur (ic, tolc)
      write  (nitre, 8010)  ic, itcur, tolc, tolcur
      write  (ncrt , 8010)  ic, itcur, tolc, tolcur
 8010 format (' ic,itcur,tolc,tolcur =',2(2x,i5),2(2x,1pe14.6))
      if (tolc .le. tolcur)  go to 2500
      if (ic .lt. itcur)
     .write  (ncrt, 8012)
 8012 format (/ ' subroutine MHD reports:',
     .         '  pprime not converged; MHD will be run again')
      if (ic .le. itcur)  go to 2100
      write  (nitre, 8020)
 8020 format (' done with outer current density iteration' /
     .        ' convergence not achieved')
 2500 write  (nitre, 8030)
      write  (ncrt , 8030)
 8030 format (' done with outer current density iteration' /
     .        ' convergence assumed or achieved')
c
c ----------------------------------------------------------------------
c --- load pcurrent for plotting purposes (stored in mhdcom.i)
c ----------------------------------------------------------------------
c
      call curdenmhd
c
c ----------------------------------------------------------------------
c determine gcap from Ampere's law and get initial rbp profile
c ----------------------------------------------------------------------
c

      call getg (ii, elongax)

c
c ----------------------------------------------------------------------
c calculate some MHD related parameters. These quantitites are also
c calculated in the transport section. We do them here in case transport
c calculations are not done.
c pvbar =volume average pressure
c betapmhd
c betatmhd
c ali was calculated in FLUXAV
c volfac = 4*pisq*rmajor (set here and in subroutine RHOMSH)
c here we use <bp**2> = (u0*I/circum)**2  =bpsqsurf
c ----------------------------------------------------------------------
c
      call trapv (r, press, hcap, nj, pvbar)
      pvbar    = pvbar*volfac/volume
      bpsqsurf = (u0*tocur/circum)**2
      betapmhd = 2.0 * u0*pvbar/bpsqsurf
      betatmhd = betapmhd*bpsqsurf/btor/btor
c
c ----------------------------------------------------------------------
c determine if an iteration of previous transport step is necessary
c ----------------------------------------------------------------------
c
 3000 call getrmaj1 (ncrt, nj)    ! calc rmajorvec,bprmaj
 3010 if (istop .ne. 1)  call prepar (dtime)
      timold = time

c
c ----------------------------------------------------------------------
c output MHD quantities to plot file
c ----------------------------------------------------------------------
c
      call outcoil (rmajor, btor, fcap, rmhdgrid, zmhdgrid, psir, nj)
      if (igoitr .eq. 0 .and. do_eqplot .eq. 1)  call eqplot(ii)
c
c --- note: EQPLOT destroys psi1d
c

c
c ----------------------------------------------------------------------
c convert back to gaussian units
c ----------------------------------------------------------------------
c
** 99 continue
      psibdry = psibdry*psikgaus
      psiaxis = psiaxis*psikgaus
      pmax    = psibdry
      pmin    = psiaxis
      flim    = flim*1.0e+06
      rmajor  = 100.0*rmajor
      volfac  = 4.0 * pisq*rmajor
      rminor  = 100.0*rminor
      btor    = 1.0e+4*btor
      rma     = rma * 100.0
      zma     = zma * 100.0
      xmagn1  = 100.0*xmagn1
      ymagn1  = 100.0*ymagn1
      xax(1)  = 100.0*xax(1)
      yax(1)  = 100.0*yax(1)
      rsep    = 100.0*rsep
      zsep    = 100.0*zsep
      rgeom   = 100.0*rgeom
      btgeom  = 1.0e+04*btgeom
      volume  = volume*1.0e6
      circum  = 100.0*circum
      cconst  = 1.0e4
      call multpl1 (cxareanpsi, npsi, cconst)
      call multpl1 (sfareanpsi, kpsi, cconst)
      cconst  = 100.0
      call multpl1 (rmajavnpsi, kpsi, cconst)
      call multpl1 (rminavnpsi, kpsi, cconst)
c
      cconst  = 1.0e+06
      call multpl1 (psivolp, kpsi, cconst)! used: PRENUB, POSTNUB, FREYA
      cconst  = 100.0
      call multpl1 (rmhdgrid,nw,cconst)
      call multpl1 (zmhdgrid,nh,cconst)
      call multpl1 (rcontr,ncontr,cconst)
      call multpl1 (zcontr,ncontr,cconst)
      call multpl1 (rplasbdry,nplasbdry,cconst)
      call multpl1 (zplasbdry,nplasbdry,cconst)
c
      rplasmax = 100.0*rplasmax
      rplasmin = 100.0*rplasmin
      zplasmax = 100.0*zplasmax
      zplasmin = 100.0*zplasmin
c
      call multpl1 (psi,nwh,psikgaus)
      call multpl1 (p,nwh,psikgaus)
      call multpl1 (psival,npsi,psikgaus)
      call multpl1 (psir,nj,psikgaus)
      cconst = 1.0e-7
      call multpl1 (pprim,nj,cconst)
      cconst = 1.0e+04
      call multpl1 (ffprim,nj,cconst)
      cconst = 1.0e+06
      call multpl1 (fpsi,npsi,cconst)
      cconst = 1.0e-4
      call multpl1 (curold, nj, cconst)    ! amps/cm**2
      call multpl1 (curden,nj,cconst)
      call multpl1 (psiold,nj,psikgaus)
      cconst = 1.0e-7
      call multpl1 (ppold,nj,cconst)
      cconst = 1.0e+4
      call multpl1 (ffpold,nj,cconst)
      cconst = 1.0e+2
      call multpl1 (rold,nj,cconst)
      cconst = 1.0e+1
      call multpl1 (press,nj,cconst)
      cconst = 100.0
      call multpl1 (rho,npsi,cconst)
c
      rhomax       = 100.0 * rhomax
      rhoa         = 100.0 * rhoa
      rhoa0        = 100.0 * rhoa0
      drhoadt_geom = 100.0 * drhoadt_geom
c
      call multpl1 (r,nj,cconst)
      call multpl1 (ra,nj-1,cconst)
      call multpl1 (dr,nj-1,cconst)
      cconst = 1.0e-2
      call multpl1 (drr,nj,cconst)
      call multpl1 (rrm,nj,cconst)
      call multpl1 (rrp,nj-1,cconst)
      cconst = 1.0e+01
      call multpl1 (bp, nj, cconst)
      cconst = 1.0e+02
      call multpl1 (rcap,nj,cconst)
      call multpl1 (rcap0,nj,cconst)
      call multpl1 (drcapdt,nj,cconst)
      call multpl1 (ravg,npsi,cconst)
      cconst = 1.0e-2
      call multpl1 (rcapi,nj,cconst)
      call multpl1 (rcap0i,nj,cconst)
      call multpl1 (drcapidt,nj,cconst)
      cconst = 1.0e-02
      call multpl1 (ravgi,npsi,cconst)
      cconst = 1.0e+06
      call multpl1 (rbp, nj, cconst)      ! gauss-cm
      cconst = 1.0e+1
      call multpl1 (pressb,nj,cconst)
      cconst = 1.0e+4
      call multpl1 (r2capi,nj,cconst)
      call multpl1 (r2capi0,nj,cconst)
      call multpl1 (dr2idt,nj,cconst)
      cconst = 1.0e+2
      call multpl1 (rmajorvec,nj,cconst)
      cconst = 1.0e+01
      call multpl1 (bprmaj,nj,cconst)     ! change from tesla to kgauss
      call multpl1 (btotrmaj,nj,cconst)     ! change from tesla to kgauss

      return
c
      end

      subroutine minmax (psi, nwh, nh, kn, psivl, iflag, knew)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c**                                                                   **
c**     MAIN PROGRAM:  MHD FITTING CODE                               **
c**                                                                   **
c**     SUBPROGRAM DESCRIPTION:                                       **
c**          minmax finds minimum and maximum value of psi in a       **
c**          cell.                                                    **
c**                                                                   **
c**     RECORD OF MODIFICATION:                                       **
c**          07/09/83..........first created                          **
c**                                                                   **
c******************************************************************* HSJ
c
      dimension psi(nwh)
c
      iflag = 0
      fmin  = MIN (psi(kn), psi(kn+1), psi(kn+nh), psi(kn+nh+1))
      fmax  = MAX (psi(kn), psi(kn+1), psi(kn+nh), psi(kn+nh+1))
      if ((psivl .lt. fmin) .or. (psivl .gt. fmax))  iflag = 1
      knew  = kn
      return
c
      end

      subroutine multpl1 (a, na, constant)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension a(na)
c
      do j=1,na
        a(j) = a(j) * constant
      end do
      return
c
      end

      subroutine multpl2 (a, b, na, constant)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension a(na), b(na)
c
      do j=1,na
        b(j) = a(j) * constant
      end do
      return
c
      end

