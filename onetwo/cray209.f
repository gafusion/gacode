       subroutine dif_w_Bp0 (psir, press, dpdpsi, nj_dum)
c
c ----------------------------------------------------------------------
c     subroutine dif_w_Bp0 evaluates
c       d Press /d Psi = d press / d rho * d rho / d psi
c     where  drho /d psi is just 1.0/(R0 Bp0)
c
c --- INPUT
c     psir(j)    j=1,2,3..nj  psi grid in MKS corresponding to r(j) grid
c     press(j)   j=1,2,3..nj  the corresponding pressure, mks
c
c --- OUTPUT
c     dpdpsi(j)  j=1,2,3..nj derivative of the pressure  wrt psi
c                            on the psir grid which is the psi grid that
c                            correspsonds to the r(nj) rho grid.
c ----------------------------------------------------------------------
c
      USE param
      USE extra
      USE numbrs
      USE mesh
      USE machin
      USE etc
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray209.f,v 1.84 2014/02/27 23:05:47 smithsp Exp $"/
c
      dimension dpdpsi(nj),psir(nj),press(nj)
      real *8 , dimension(:),allocatable :: dpressdrho 
c

c
      call difydx (r, press, dpdpsi, nj) ! get dpressdrho..
c                                       ..(store in dpdpsi temporarily)

      if(use_Bp0 .eq.3)then
           if( .not. allocated(dpressdrho))then
              allocate (dpressdrho(size(press)),STAT = istat)
              if(istat .ne. 0)
     .          call allocate_error("dpressdrho, dif_w_Bp0",0,istat)
           endif
           call difydx (r, press, dpressdrho, nj) ! get dpressdrho
           axis_val = dpressdrho(1)
           dpressdrho = dpdpsi                   !save dpdrho for use below
      endif

      frac = 0.5
      dpdrhoa = frac*dpdpsi(2) ! value of dpdrho at j=3/2 (since dP/drho =0 on axis)
      do j=2,nj
         dpdpsi(j)=dpdpsi(j)/(rmajor*bp(j))
      end do
c
      if (use_Bp0 .le. 1) then
c
c       get the value of the indeterminate ratio dpdrho/(R0*Bp0)
c       at the magnetic axis by moving slightly away from the axis:
c
        bpa        =  bp(2)               
c        pressa    = frac*press(1)+(1-frac)*press(2))  ! ditto
         dpdpsi(1) = dpdrhoa/(rmajor*bpa)
         dpdpsi(1) = 1.5*dpdpsi(2)
      else if (use_Bp0 .eq. 2) then ! get smooth values near axis by polynomial extrapolation
        npoints = 10  ! (3 min,limited to 10 in sub POLINT)
        do j=use_Bp0-1,1,-1 
          call POLINT (psir(use_Bp0),dpdpsi(use_Bp0),
     .                 npoints,psir(j),dpdpsi(j),DY)
        end do
      else if (use_Bp0 .eq. 3) then
         dpdpsi(1)= (q(1)/btor)*axis_val
         do j=2, size(dpdpsi)
            dpdpsi(j) = (q(j)/btor)*(dpressdrho(j)/r(j))
         enddo
      else
         dpdpsi(1) = dpdpsi(2)
      end if
      return
c
      end

      subroutine new_point (pds, arcl, xn, yn, thet, twopi, spiov4,
     .                      xmin, xmax, ymin, ymax, yaxd, xaxd, ier,
     .                      sint, cost, xns, yns, thetnew)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- subroutine used by CNTOUR to get approximation for next point on
c --- contour,given the current point and the gradient of psi .
c --- the tangent line at the point (xn,yn) is
c             df = 0=df/dx*dx+df/dy*dy
c   a circle of radius arcl,centered at (xn,yn) is
c             arcl**2 = dr**2+dz**2
c   this gives us two eq. in two unknowns, dr and dz.
c   the new point is xn +/- dx and yn +/- dy where the signs have to
c   be picked so that thet increases. Note the special treatment
c   required as thetnew crosses zero (outside this routine).
c
c ------------------------------------------------------------------ HSJ
c
      dimension pds(*)
c
      ier = 0
      dpx = pds(2)
      dpy = pds(3)
      if (ABS (dpx) .gt. ABS (dpy)) then
        if (dpx .eq. 0.0) then
          ier = 1
          return
        end if
        alpha = dpy/dpx
        dy    = arcl / SQRT (1.0+alpha*alpha)
        dx    = -alpha * dy
      else
        if (dpy .eq. 0.0) then
          ier = 1
          return
        end if
        alpha = dpx/dpy
        dx    = arcl / SQRT (1.0 + alpha*alpha)
        dy    = -alpha * dx
      end if
c
c ----------------------------------------------------------------------
c --- the sign on dx,dy must be taken so that thet increases.
c --- a unit vector in the direction of increasing thet
c --- (i.e., thet counterclockwise) is (-SIN (thet),COS (thet))
c --- the displacement vector is (dx,dy). its projection on the above
c --- vector must be positive and equals -dx*SIN (thet)+dy*COS (thet)
c ----------------------------------------------------------------------
c
      proj = -dx*sint+dy*cost
      if (proj .lt. 0.0) then
        dx = -dx
        dy = -dy
      end if
      xns = xn+dx
      yns = yn+dy
      if (xns .lt. xmin .or. xns .gt. xmax .or.
     .    yns .lt. ymin .or. yns .gt. ymax) then
        ier = 1
        return
      end if
      thetnew = ATAN2 (yns-yaxd, xns-xaxd)
      if (thetnew .lt. 0.0) then
        thetnew = thetnew+twopi
      else if (thet .gt. spiov4) then
        thetnew = thetnew+twopi
      end if
      return
c
      end

      subroutine outcoil (rmajor, btor, fcap, rmhdgrid,
     .                    zmhdgrid, psir, nj)
c
      USE param
      USE io 
      USE solcon
      USE mhdpar        
      USE soln2d
      USE mhdcom
      USE etc
      USE mhdbcdtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subroutine outputs fluxes on the loops and the coil currents
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'etc.i'
c      include 'solcon.i'
c      include 'io.i'
c      include 'mhdcom.i'
c      include 'mhdbcdtn.i'
      include 'storage.i'
c      include 'soln2d.i'
c
c --- get stark measurements, even if ifixshap = 1
c
      if (use_stark) then
c
c       first do the calculations
c
          call stark (time,nout,ncrt,rmajor,btor,fcap,
     .                rmhdgrid,zmhdgrid,psir,nj)
c
c       next do the printout
c
          tn = n
          call header (nout,time,tn)
          write  (nout, 6000)
 6000     format ('  =============== MOTIONAL STARK EFFECT RESULTS ==',
     .              '=================' / 12X,'RLOC',10X,'ZLOC',7X,
     .              'TAN(Z)EXPTL',4X,
     .              'TAN(Z)CALC',5X,'CHISQ',7x,'SIGMA')
          do j=1,nstark
            if (rstarkexptl(j) .gt. 0.0) then
               write  (nout, 6010)  j,
     .                 rstarkexptl(j), zstarkexptl(j),   tstarkexptl(j),
     .                  tstarkcalc(j),  chisqstark(j), sigstarkexptl(j)
 6010          format (2x, i3, 6(2x, 1pe12.4))
            end if
          end do
          write  (nout, 6020)  chisqstarktot
 6020     format (' total chisq, MSE only =', 1pe12.4)
      end if
c
      if (mhdmode .ne. 'coils')  return
c
      xdsumt = 0.0
      xdsum  = 0.0
c
      do i=1,nfcoil
        if (nobsrvd .ne. 0) then
          xdum(i) = ((fwtfcoil(i)*(expfcoil(i)-curfcoil(i)))**2)/nobsrvd
          xdsum   = xdsum + xdum(i)
        else
          xdum(i) = 0.0
        end if
      end do
c
      write (nout , 8100)
      write (nitre, 8100)
c
      do i=1,nfcoil
        write (nout , 8120)  i, expfcoil(i), curfcoil(i), xdum(i)
        write (nitre, 8120)  i, expfcoil(i), curfcoil(i), xdum(i)
      end do
c
      write (nout , 8130)  xdsum
      write (nitre, 8130)  xdsum
 8100 format ('1', 5x, 'coil', 3x, ' exp  cur', 4x, 'calc current',
     .           5x, 'red. chisq' /
     .           28x, 'amps'      / 5x,
     .'------------------------------------------------------------')
 8120 format (i8, 2x, 1pe12.3, 3x, 1pe12.3, 3x, 1pe12.3)
c
      xdsumt = xdsumt + xdsum
      xdsum  = 0.0
      do i=1,nsilop
        if (nobsrvd .ne. 0) then
          xdum(i) = ((fwtpsilp(i)*(psiloop(i)-psilopcl(i)))**2)/nobsrvd
          xdsum   = xdsum + xdum(i)
        else
          xdum(i) = 0.0
        end if
      end do
c
      write  (nout , 8200)
      write  (nitre, 8200)
 8200 format ('  psi loop  exptl    calc (voltsec)',
     .       '   red. chi-square' /
     . '  ------------------------------------------------------------')
c
      do i=1,nsilop
        write (nout , 8120)  i, psiloop(i), psilopcl(i), xdum(i)
        write (nitre, 8120)  i, psiloop(i), psilopcl(i), xdum(i)
      end do
c
      write (nout , 8130) xdsum
      write (nitre, 8130) xdsum
c
      xdsumt = xdsumt + xdsum
      xdsum  = 0.0
c
      do i=1,magpr2
        if (nobsrvd .ne. 0) then
          xdum(i) = ((fwtmprbe(i)*(probeval(i)-probeclc(i)))**2)/nobsrvd
          xdsum   = xdsum + xdum(i)
        else
          xdum(i) = 0.0
        end if
      end do
c
      write (nout , 8300)
      write (nitre, 8300)
 8300 format (
     .  '  mag probe  exptl    calc (tesla)   red. chisq' /
     .  ' ------------------------------------------------------------')
c
      do i=1,magpr2
        write (nout , 8120) i, probeval(i), probeclc(i), xdum(i) !BOUNDS
        write (nitre, 8120) i, probeval(i), probeclc(i), xdum(i) !BOUNDS
      end do
c
      write (nout , 8130) xdsum
      write (nitre, 8130) xdsum
 8130 format (2x, 21(1h ), 'sum red. xchisq =', 1pe12.3 ///)
      xdsumt = xdsumt + xdsum
      write (nout , 8140) xdsumt
      write (nitre, 8140) xdsumt
 8140 format ('  total reduced xchi-square =', 1pe12.3)
c
c --- predicted and calculated loop voltages:
c
      write (nout , 8190)  vlopmhdt, vloopcl(ieq+1)
      write (nitre, 8190)  vlopmhdt, vloopcl(ieq+1)
 8190 format ('  exptl and calc loop voltage :', 2(1pe12.3,2x))
      return
c
      end

      subroutine outgo
c
      USE param
      USE io 
      USE soln
      USE mhdpar 
      USE extra 
      USE numbrs      
      USE mesh
      USE geom
      USE soln2d
      USE psig
      USE rhog
      USE flxav
      USE etc
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c output from main program
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'etc.i'
c      include 'extra.i'
c      include 'flxav.i'
c      include 'geom.i'
c      include 'io.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'psig.i'
c      include 'rhog.i'
c      include 'soln.i'
c      include 'soln2d.i'
c
      do iunit=nout,nitre
        call headeq(iunit)
        write (iunit, 1010) ieq
 1010   format (/ ' subroutine OUTGO reports: equilibrium number', i3)
        write (iunit, 1020)
 1020   format (/ ' values on psi grid'//5x,'i',8x,'psi',12x,
     .            'rho',10x,'vprime',6x,'<R0**2/R**2>',5x,'q(psi)')
c
        do 30 i=1,npsi,j2prt
   30   write (iunit, 1030)
     .        i, psival(i), rho(i), vprime(i), ratave(i), qpsi(i)
c
 1030   format (' ', i5, 7(1pe15.5))
        write (iunit, 1040)
 1040   format (/ ' values on rho grid'//
     .              5x,'i',8x,'psi',12x,'rho',12x,'bp',11x,'curden',
     .              11x,'q',12x,'fcap',11x,'hcap')
        do 50 j=1,nj,jprt
   50   write (iunit, 1030)
     .        j, psir(j), r(j), bp(j), curden(j), q(j), fcap(j), hcap(j)
      end do
      return
c
      end

      subroutine pfprim (ic)
c

c
c ----------------------------------------------------------------------
c  this subroutine
c  a) renormalizes the current density curden to totcur
c  b) differentitates the pressure to get pprim
c  c) determines ffprim from curden and pprim
c --- input
c  through argument list
c  ic                    switch,ic = 0 don't renormalize curden to totcur
c                        ic  .ne. 0 means renormalize curden to tocur
c                        whenever ic .gt. 1 pprim is relaxed using omcur
c
c  through INCLUDE files:
c  INCLUDE file param:
c  kj                    dimension of vectors indexed by nj below
c
c  INCLUDE file numbrs:
c  nj                       size of vectors indexed by nj below
c
c  INCLUDE file geom
c  hcap(1..nj)           f(psilim)/f(psi)/r2cap
c  r2cap(1...nj)         <R0**2/R**2>
c
c  INCLUDE file etc:
c  tocur                 total toroidal current,amps
c  INCLUDE file constnts:
c  pi
c  u0                    4 * pi * 1.0e-07
c
c  INCLUDE file mesh:
c  r(1...nj)             transport rho grid,meters
c
c  INCLUDE file soln:
c  curden(1...nj)        <jtor*R0/R>,amps/m**2
c
c  INCLUDE file soln2d:
c  omcur                  relaxation parameter
c
c  INCLUDE file rhog:
c  psir(1..nj)            psi values corresponding to r grid
c  press(1...nj)          pressure,nt/m**2
c
c  INCLUDE file machin
c  rmajor                 =R0,location of btor,meters
c
c --- temporary storage required
c  INCLUDE file storage
c  xdum(1..nj)
c  ydum(1..nj)                 vectors min length nj
c
c --- output
c  through INCLUDE files
c  INCLUDE file soln:
c  curden(1..nj)              renormalized curden,if ic .ne. 0
c                             if ic = 0 curden is unchanged on output
c
c  INCLUDE file rhog:
c  pprim(1..nj)               d(press)/dpsi,nt/(m**2*volt*sec)
c  ffprim(1...nj)             f(psi)*d(f)/dpsi,
c
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE soln
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE constnts
      USE soln2d
      USE rhog
      USE etc
      USE mhdcom,             ONLY : rma
      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'storage.i'
c
c ----------------------------------------------------------------------
c renormalize current density if ic .ne. 0:
c note that the integral done below is exactly
c integral(jtoroidal*DR*DZ). It looks different due to the fact
c that curden =<jtoroidal*R0/R> but the result is the same.
c get the flux surface average current density,<j*R0/R>,based
c on the most recent pprim and ffprim:
c ----------------------------------------------------------------------
c
      if (ic .lt. 0) then  ! ic is never < 0 so effectivly skip this
        do j=1,nj
          curden(j) =
     .           -(rmajor*pprim(j) + ffprim(j) * r2cap(j) / (u0*rmajor))
        end do
        do j=1,nj
          xdum(j) = 2.0 * pi * hcap(j)*r(j)*curden(j)
        end do
        call trap2 (r, xdum, ydum, nj)
        fac = tocur / ydum(nj)
        call multpl1 (curden, nj, fac)
      end if
c
c ----------------------------------------------------------------------
c calculate pprim in MKS; underrelax to avoid numerical problems
c ----------------------------------------------------------------------
c
      relax = 1.0
      if (ic .gt. 1) relax = omcur
      if (use_Bp0 .eq. 0) then
        call    difydx (psir, press, ydum, nj)
      else
        call dif_w_Bp0 (psir, press, ydum, nj)
      end if
      do 60 j=1,nj
   60 pprim(j) = (1.0 - relax) * pprim(j) + relax * ydum(j)

c   
c ----------------------------------------------------------------------
c calculate ffprim in MKS. Note that  r2cap = < R0**2/R**2>
C rmajor on psi grid at z=0 is ??
c ----------------------------------------------------------------------
c
      do  j=1,nj
         ffprim(j) = -u0*rmajor*(curden(j)+rmajor*pprim(j)) / r2cap(j)
      enddo
      curaxis = -rma *pprim(1) - ffprim(1)/(u0*rma)


      return
c
      end

      subroutine prepar (dtime)
c

c
c ----------------------------------------------------------------------
c   this subroutine determines if it is necessary to iterate on
c   the last transport step (because rho at plasma boundary or the
c   geometric (i.e., ...cap) parameters are not sufficiently close
c   to their predicted values).
c
c   if solution is converged, icon = 0 is returned and
c         1) eqdsk is written
c         2) the solution vector u is updated
c         3) time derivatives of rhoa and cap quantitites are updated
c         4) save file pointers for files plotdsk, outone and quikone
c         5) if ieq = 0 print out geometric parameters
c   if solution is NOT converged, icon = 1 is returned and
c         1) read back the old solution vector (written by savit)
c         2) reset the file pointers to write over the transport
c            solution which is not acceptable
c         3) calculate the time derivative of rhoa and cap parameters
c
c   igoitr is a redundant switch here, it follows the value of icon
c ----------------------------------------------------------------------
c
      USE param
      USE io 
      USE solcon
      USE soln
      USE contour
      USE limiter
      USE mhdpar 
      USE mhdgrid
      USE rf 
      USE extra    
      USE yoka  
      USE numbrs
      USE toq_12_interface, only : toq_drive
      USE mesh
      USE machin
      USE geom
      USE tordlrot
      USE constnts
      USE soln2d
      USE bd_condtn,only : totcur
      USE rhog
      USE mhdcom
      USE shapctr
      USE flxav
      USE neo2d
      USE etc
c      USE  io,                        ONLY : ioftn, iopntr, nunits
      USE gpsi
      implicit  integer (i-n), real*8 (a-h, o-z)


c      include 'small.i'

c
      external  LENGTH
      logical   opened
      integer   pointer
      character record*132
c
****  write (6, *)  'on entry to PREPAR, nj =', nj
      igoitr  =  0
****  derwght = -0.5   ! set in cray101.f or read from inone

c
c ----------------------------------------------------------------------
c  ieq is zero only if TPORT has never been called. this can happen only
c  on the startup calculation, for the initial equilibrium. once TPORT
c  has been called (even if the MHD/transport cycle did not converge)
c  ieq .ge. 1. (ieq is incremented each time an MHD/transport cycle
c  is converged.)
c  itre counts the number of MHD/transport cycles done in order to
c  converge the transport and geometric parameters. itre starts
c  over again at 1 each time an MHD/transport cycle is converged.
c ----------------------------------------------------------------------
c



      if (ieq .eq. 0)  go to 500
c
c ----------------------------------------------------------------------
c  subroutine CHEKGM compares the (linearly) predicted cap parameters with the
c  ones just calculated using the most recent equilibrium.
c  for ieq = 1 the transport model was initially advanced assuming that the
c  derivatives of the cap quantities was zero. Hence the check in
c  chekgm for ieq = 1 is also appropriate.
c  for ieq = 1 we assume converged solution,so as to avoid startup
c  perturbations.
c ----------------------------------------------------------------------
c
c not converged due to failure
c reset everything, cut time step in half in subroutine RUNTWO
c
      call chekgm
      if (ieqfail .ne. 0)  icon = 1
c
c ----------------------------------------------------------------------
c --- if this is the initial transport/equilibrium cycle (i.e., ieq = 1)
c --- and we did not converge,even when the time step dteq was cut
c --- down, then proceed anyway
c ----------------------------------------------------------------------
c
      if (itre .eq. maxitr .and.      ieq .eq. 1  .and.
     .    icon .eq. 1      .and. 0.5 * dteq .le. dtmine) then
        write  (nout , 1) maxitr
        write  (nitre, 1) maxitr
        write  (ncrt , 1) maxitr
    1   format ('  startup calculations did not converge in', i5     /
     .          '  iterations. we will proceed assuming convergence' /
     .          '  was obtained.')
        icon = 0
      end if
      if (icon .eq. 0)  go to 500
c
c ----------------------------------------------------------------------
c  not converged, iteration necessary. This represents essentially the
c  whole ball game when MHD/transport cycles are run. that is,we are
c  just using the transport calculations as a complicated way to
c  find the fixed points of the cap functions. the right way to do
c  this problem, with Newton's method, is out of the question here since
c  it would take too many MHD/transport cycles to evaluate the Jacobian
c  or some approximation to it. Hence we do the following instead.
c  for the first two iterations of any MHD/transport cycle we use
c  fixed point iteration, with underrelaxation (i.e., derwght = -0.5).
c  however even this relaxation is not effective in many cases
c  (the gradient of the mappings is often greater than one).
c  hence after we have two previous solutions,
c  the next guess is calculated analytically,
c  assuming that the jacobian is diagonal. after we have three previous
c  solutions, we solve the problem analytically assuming a parabolic
c  shape for the cap functions.
c
c  1. get old solution vector by reading binary file "savsol"
c  2. update time derivatives drhoadt_geom, etc. (using above methods)
c  3. do not update rhoa0, gcap0, etc. (but if steady_state=0 then relax parameters)
c  4. reset disk address to write over last transport calculation
c  5. copy scratch2 file into scratch1.
c    (scratch1 was created in the most recent call to FLUXAV.
c     however since we will not keep this equilibrium we want
c     to restore the old scratch1 file which is currently stored
c     as scratch2).
c
c ----------------------------------------------------------------------
c

      call COPY_TEXT_FILE ('scratch2', 'scratch1')

c
      igoitr = 1
c
c     underrelaxed fixed point iteration
c
      do j=1,nj
        prep    = gcap0(j)+dgdt(j)*dtime
        dgdt(j) = ((1.0+derwght)*gcap(j)-derwght*prep-gcap0(j))/dtime
        dgdt(j) = dgdt(j)* cap_mult
c
        if (.not. implicit_fh) then
          prep    = fcap0(j)+dfdt(j)*dtime
          dfdt(j) = ((1.0+derwght)*fcap(j)-derwght*prep-fcap0(j))/dtime
          dfdt(j) = dfdt(j)* cap_mult
c
          prep    = hcap0(j)+dhdt(j)*dtime
          dhdt(j) = ((1.0+derwght)*hcap(j)-derwght*prep-hcap0(j))/dtime
          dhdt(j) =  dhdt(j)* cap_mult
        else
          dfdt(j) = dfdtcap(j)* cap_mult
          dhdt(j) = dhdtcap(j)* cap_mult
        end if
c
        prep     = r2cap0(j)+dr2dt(j)*dtime
        dr2dt(j) = ((1.0+derwght)*r2cap(j)-derwght*prep-r2cap0(j))/dtime
        dr2dt(j) = dr2dt(j)* cap_mult
c
        prep       = rcap0(j)+drcapdt(j)*dtime
        drcapdt(j) = ((1.0+derwght)*rcap(j)-derwght*prep-rcap0(j))/dtime
        drcapdt(j) = drcapdt(j)* cap_mult
c
        prep       = rcap0i(j)+drcapidt(j)*dtime
        drcapidt(j) = ((1.0+derwght)*rcapi(j)-derwght*prep-rcap0i(j))/
     .                                              dtime !Fixed for too many columns 2/27/2014 SPS
        drcapidt(j) = drcapidt(j)* cap_mult
c
        prep       = r2capi0(j)+dr2idt(j)*dtime
        dr2idt(j)  =
     .           ((1.0+derwght)*r2capi(j)-derwght*prep-r2capi0(j))/dtime
        dr2idt(j)  = dr2idt(j)* cap_mult
c
        if (xi_include) then
          prep       = eps0(j)+depsdt(j)*dtime
          depsdt(j)  = ((1.0+derwght)*eps(j)-derwght*prep-eps0(j))
     .               / dtime
          depsdt(j)  = depsdt(j)* cap_mult
c
          prep       = xhm20(j)+dxhm2dt(j)*dtime
          dxhm2dt(j) = ((1.0+derwght)*xhm2(j)-derwght*prep-xhm20(j))
     .               / dtime
          dxhm2dt(j) =  dxhm2dt(j)* cap_mult
c
          prep       = xi110(j)+dxi11dt(j)*dtime
          dxi11dt(j) = ((1.0+derwght)*xi11(j)-derwght*prep-xi110(j))
     .               / dtime
          dxi11dt(j) = dxi11dt(j)* cap_mult
c
          prep       = xi330(j)+dxi33dt(j)*dtime
          dxi33dt(j) = ((1.0+derwght)*xi33(j)-derwght*prep-xi330(j))
     .               / dtime
          dxi33dt(j) = dxi33dt(j)* cap_mult
c
          prep       = xips0(j)+dxipsdt(j)*dtime
          dxipsdt(j) = ((1.0+derwght)*xips(j)-derwght*prep-xips0(j))
     .               / dtime
          dxipsdt(j) = dxipsdt(j)* cap_mult
c
        end if
      end do
      IF(ABS(steady_state) .lt. 1.e-5)THEN
              rhoa0        = 0.5*(rhoa     +  rhoa0)
              gcap0(:)     = 0.5*(gcap(:)  +  gcap0(:))
              fcap0(:)     = 0.5*(fcap(:)  +  fcap0(:))
              hcap0(:)     = 0.5*(hcap(:)  +  hcap0(:))
              r2cap0(:)    = 0.5*(r2cap(:) +  r2cap0(:))
              r2capi0(:)   = 0.5*(r2capi(:)+  r2capi0(:))
              rcap0(:)     = 0.5*(rcap(:)  +  rcap0(:))
              rcap0i(:)    = 0.5*(rcapi(:) +  rcap0i(:))
      ENDIF
c
      prep   = rhoa0+drhoadt_geom*dtime
      drhoadt_geom = ((1.0+derwght)*rhoa-derwght*prep-rhoa0)/dtime
      drhoadt_geom = drhoadt_geom * cap_mult
c
c --- calc the loop voltage even though its not converged at this point.
c --- (The value is printed out for inspection in subroutine OUTCOIL.)
c --- the final converged value will be picked up under the icon = 1
c --- section of this subroutine (see further below).
c
      dpsivlop       = (psivloop-psivlop0) / dtime
      vloopcl(ieq+1) =  voltoh + dpsivlop
c
c ----------------------------------------------------------------------
c --- reset file pointers for files output, outqik and plotdsk
c --- address for next line set to iopntr(i)
c ----------------------------------------------------------------------
c
c      print *, 'prepar, rewinding units'
      do 300 i=1,nunits-1
        iunit = ioftn(i)
        call EMPTY (iunit)
        if (i .eq. 1  .and.  iprtit .eq. 1)  go to 300
        inquire (unit = iunit, opened = opened)
        if (opened) then
          rewind (unit = iunit)
        else
          call STOP ('subroutine PREPAR: file is not open', 95)
        end if
        if (iopntr(i) .gt. 0) then
          nchars = 0
          do while (nchars .lt. iopntr(i))
            read (iunit, '(a)')  record
            nchars = nchars + LENGTH (record)
          end do
        end if
****    call POSTEXT (iunit, iopntr(i), icount)
  300 continue
c      print *, 'prepar, done rewinding units'
c
c --- read transport solution, binary file "savsol", written earlier
c --- by subroutine SAVIT.
c --- the MHD section of the code assumes MKS units,
c --- so convert r and psi.
c --- restore the boundary condition vector, psincbcd,
c --- for the no coils option.
c
****  write (6, *)  'before calling READIT, nj = ', nj
      print*, 'calling READIT'
      call readit
      print*, 'done READIT'
      cconst = 1.0e-2
      call multpl1 (r, nj, cconst)
      call multpl1 (psi , nwh, psimks)
      call multpl1 (psir, nj , psimks)
      call copya (psi, p, nwh)
      psibdry = psibdry*psimks
      psiaxis = psiaxis*psimks
      call load_psincbcd
      rma = 0.01*rma
      zma = 0.01*zma
      xmagn1 = rma
      ymagn1 = zma
c
c     restore values that are not read from readit back to the values
c     they had at the beginning of this MHD/transport cycle:
c
      do j=1,nj
        hcap(j)   = hcap0(j)
        fcap(j)   = fcap0(j)
        gcap(j)   = gcap0(j)
        r2cap(j)  = r2cap0(j)
        r2capi(j) = r2capi0(j)
        rcap(j)   = rcap0(j)
        rcapi(j)  = rcap0i(j)
        eps(j)    = eps0(j)
        xhm2(j)   = xhm20(j)
        xi11(j)   = xi110(j)
        xi33(j)   = xi330(j)
        xips(j)   = xips0(j)
      end do
 
      go to 2300
c
c ----------------------------------------------------------------------
c no iteration necessary
c ----------------------------------------------------------------------
c  1. write an eqdsk for this equilibrium
c  2. update solution vector
c  3. save solution vector, etc. in binary disk file "savsol"
c     in case a later iteration is necessary
c     (this step is done in subroutine RUNTWO, where SAVIT is called)
c  4. update time derivatives
c  5. update rhoa0, gcap0, etc
c  6. for ieq = 0, print out geometric factors
c  7. record file pointers for plotdsk, output and outqik files
c  8. copy file scratch1 into file scratch2
c     file scratch2 will be copied back into file scratch1 whenever
c     the MHD/transport cycle did not converge
c ----------------------------------------------------------------------
c
  500 continue


      if (ieqdsk .ne. 0 .and. mhdmethd.ne.'toq') then

         call wrteqdsk (psi, rmhdgrid, zmhdgrid, rmajor, btor, totcur,
     .                  press(nj), nj, xsep, ysep, rcontr, zcontr,
     .                  ncontr, psibdry, psiaxis, xax, yax, versid,
     .                  xlimiter, ylimiter, nlimiter, ishot, ieqdsk,
     .                  use_Bp0,nion, nprim, nimp, ncrt, nout)
       else if(mhdmethd .eq. 'toq') then
          call using_toq_in_prepar
       endif


      !call TOQ with the eqdsk jsut created:
!      if (mhdmethd .eq. 'toq')then

         !get toq executable, write toq input file, spawn toq.
         !data passed to toq through eqdsk
!         call toq_drive(tocur)



 !     endif



c
c --- copy scratch1 into scratch2 (in case we have to restore scratch1
c --- at a later time due to non-convergence of the MHD/transport cycle)
c
c      print *,'converged section, scratch copy'
      call COPY_TEXT_FILE ('scratch1', 'scratch2')
c      print *,'converged section, done scratch copy'
c
      do 570 k=1,krf
  570 if (rfon(k) .gt. timmax)  go to 590
      jrf1min = jrf2min
      jrf1max = jrf2max
      do j=1,nrfrad
        rfrow1(j) = rfrow2(j)
      end do
c
 590  voladj1 = voladj2
      eqtim1 = eqtim2
      voladj2 = voladj
      eqtim2 = time
c
      if (ieq .eq. 0)  go to 700
c
c ----------------------------------------------------------------------
c calculate the time derivatives of the cap quantitites:
c if this is the final time then dtime = 0.0 so don't update derivatives:
c ----------------------------------------------------------------------
c
      if (dtime .ne. 0.0) then
          do 600 j=1,nj
               dgdt(j) =cap_mult * (gcap(j)-gcap0(j))/dtime
            if (.not. implicit_fh) then
                dfdt(j) = cap_mult *(fcap(j)-fcap0(j))/dtime
                dhdt(j) = cap_mult * (hcap(j)-hcap0(j))/dtime
            else
                dhdtcap(j) = cap_mult *dhdt(j)
                dfdtcap(j) = cap_mult *dfdt(j)
            end if
            drcapdt(j) = cap_mult *(rcap(j)-rcap0(j))/dtime
            drcapidt(j) = cap_mult *(rcapi(j)-rcap0i(j))/dtime
            dr2idt (j) = cap_mult *(r2capi(j)-r2capi0(j))/dtime
          if (xi_include) then
            depsdt (j) =cap_mult * (eps(j)-eps0(j))/dtime
            dxhm2dt(j) = cap_mult *  (xhm2(j)-xhm20(j))/dtime
            dxi11dt(j) =cap_mult * (xi11(j)-xi110(j))/dtime
            dxi33dt(j) =cap_mult * (xi33(j)-xi330(j))/dtime
            dxipsdt(j) =cap_mult * (xips(j)-xips0(j))/dtime
          end if
  600       dr2dt(j) =cap_mult * (r2cap(j)-r2cap0(j))/dtime
          drhoadt_geom =cap_mult * (rhoa-rhoa0)/dtime
          dpsivlop = (psivloop-psivlop0)/dtime
          vloopcl(ieq+1) = voltoh + dpsivlop   ! calculated loop voltage
      end if
c      print *,'to label 1300 in prepar'
      go to 1300
c
c ----------------------------------------------------------------------
c ieq = 0
c we don't know the time derivatives if this is an MHD/transport coupled
c run so we make an initial guess of zero for the derivatives.
c if this is a standard run (i.e., one that just uses a single eqdsk for
c all times) then we also assume zero derivaties.
c Finally if this is a multiple eqdsk run then these derivatives
c are know appriori because the mhd equilibria were pre computed
c so use this information here:
c ----------------------------------------------------------------------
c
  700 do 800 j=1,nj
      dfdt(j)    = 0.0
      dgdt(j)    = 0.0
      dhdt(j)    = 0.0
      dfdtsv(j)  = 0.0
      dhdtsv(j)  = 0.0
      dfdtcap(j) = 0.0
      dhdtcap(j) = 0.0
      dr2idt(j)  = 0.0
      drcapdt(j) = 0.0
      drcapidt(j) = 0.0
      depsdt(j)  = 0.0
      dxhm2dt(j) = 0.0
      dxi11dt(j) = 0.0
      dxi33dt(j) = 0.0
      dxipsdt(j) = 0.0
  800 dr2dt(j)   = 0.0
      drhoadt_geom = 0.0
      psivlop0   = psivloop
      vloopcl(ieq+1) = 0.0
c
c ----------------------------------------------------------------------
c print out geometric quantities
c ----------------------------------------------------------------------
c
 1300 do iunit=nout,nitre
        call headeq (iunit)
        write (iunit, 900)
  900   format (/ ' geometric quantities' //
     .            5x, 'i', 8x, 'rho', 11x, 'fcap', 11x, 'gcap',
     .           11x, 'hcap', 7x, '<R0**2/R**2>')
        do j=1,nj,jprt
          write (iunit, '(i6, 7(1pe15.5))')
     .                     j, r(j), fcap(j), gcap(j), hcap(j), r2cap(j)
        end do
      end do
c
c --- update MHD-derived quantities to current time:
c --- timcap records the time at which the cap0 are defined.
c --- timcap is used in subroutine RHOMSH to define values for the
c --- cap parameters.
c
      timcap   = time    ! defined here and nowhere else
      rhoa0    = rhoa
      psivlop0 = psivloop
      xmagn0   = xmagn1
      ymagn0   = ymagn1
c
      if (ieq .eq. 0 .or. cap_mult .gt. 0.0) then
        do j=1,nj
          gcap0(j) = gcap(j)
          if (.not. implicit_fh) then
            fcap0(j) = fcap(j)
            hcap0(j) = hcap(j)
          end if
          r2cap0(j)  = r2cap(j)
          r2capi0(j) = r2capi(j)
          rcap0(j)   = rcap(j) 
          rcap0i(j)   = rcapi(j) 
          eps0(j)    = eps(j)
          xhm20(j)   = xhm2(j)
          xi110(j)   = xi11(j)
          xi330(j)   = xi33(j)
          xips0(j)   = xips(j)
        end do
      end if
c
c ----------------------------------------------------------------------
c empty buffers and record current values of file pointers
c ----------------------------------------------------------------------
c
c      print *,'recording pointers in prepar'
       do i=1,nunits-1
        iunit = ioftn(i)
        call EMPTY_AND_RECORD_POINTER (iunit, pointer)
        iopntr(i) = pointer
      end do
c      print *,'done recording pointers in prepar'
c
      if (renormalize_rbp .eq. 0) then
        do 2200 j=1,nj
        do 2100 k=1,nk-3-iangrot
 2100     u(k           ,j) = en(j,k)
          u(nk-iangrot-2,j) = te(j)
          u(nk-iangrot-1,j) = ti(j)
          u(nk-iangrot  ,j) = rbp(j)*1.0e+6    ! gauss-cm
 2200     if (iangrot .eq. 1)  u(nk,j) = angrot(j)
      end if
c
c ----------------------------------------------------------------------
c define mesh quantities
c (r was set back to its old value in subroutine READIT above if icon=1)
c ----------------------------------------------------------------------
c
 2300 continue
****  write (6, *)  'at label 2300, nj=', nj
c
      do 2500 j=2,nj
      dr(j-1) = r(j)-r(j-1)
 2500 ra(j-1) = 0.5 * (r(j-1)+r(j))
      drr(1)  = 2.0/dr(1)
      rrp(1)  = 2.0/dr(1)
      do 2600 j=2,nj
        rrm(j) = (ra(j-1)/r(j))/dr(j-1)
        if (j .eq. nj)  go to 2600
        drr(j) = 2.0/(dr(j-1)+dr(j))
        rrp(j) = (ra(j)/r(j))/dr(j)
 2600 continue
      drr(nj) = 2.0/dr(nj-1)
c
****  write (6, *)  'on return from PREPAR, nj =', nj

      return
c
      end

      subroutine pressr (imks, iuse)
c
c
c ----------------------------------------------------------------------
c   calculates pressure from density and temperature profiles
c   results returned in MKS units (newtons/m**2) if imks = 1
c   results returned in CGS units   (dyne/cm**2) if imks = 0
c   index 1 corresponds to magnetic axis,index nj corresponds to plasma edge
c   if iuse = 1 then en,ene,te,ti,wbeam,walp are used to calculate the
c   pressure. if iuse = 2 then the corresponding values are taken from
c   the vector u(kk,kj) and ene,wbeam,walp are used.
c   if iuse = 3 then same as for iuse=2 but use usave instead of u
c   and enesav instead of ene.
c --- input:
c  ene(1...nj)            electron density,#/cm**3
c  en(1..nj,1..nion)      ion density,#/cm***3
c  wbeam(1..nj)           stored beam energy density,keV/cm**3
c                         (all energy components and all beam lines)
c  walp(1...nj)           stored alpha particle energy density,keV/cm**3
c  u or usave             may be input instead,see above,switch iuse
c --- output:
c  press(1..nj)           total pressure,j/m**3 (nt/m**2) if imks = 1
c                         total pressure,erg/cm**3 (dyne/cm**2) if imks = 0
c  pressb(1..nj)          beam pressure,j/m**3 (nt/m**2) if imks = 1
c                         beam pressure,erg/cm**3 (dyne/cm**2) if imks = 0
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE fusion
      USE io 
      USE nub2
      USE soln
      USE numbrs
      USE mesh
      USE rhog
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'numbrs.i'
c      include 'fusion.i'
c      include 'io.i'
c      include 'mesh.i'   ! r grid for zero gradient enforcement HSJ 2/1/01
c      include 'nub2.i'
c      include 'rhog.i'
c      include 'soln.i'
c
      twothirds = 2.0 / 3.0
      if (imks .eq. 0) then
        convert = 1.6022e-09
      else
        convert = 1.6022e-10
      end if
c
      do j=1,nj
          enpos = 0.0
          if (iuse .eq. 1) then
              do k=1,nion
                enpos = enpos+en(j,k)
              end do
              press(j) = ene(j)*te(j) + enpos*ti(j)
     .                 + twothirds * (wbeam(j) + walp(j))
          else if (iuse .eq. 2) then
              do k=1,nion
                enpos = enpos+u(k,j)
              end do
              press(j) = ene(j)*u(nion+1,j) + enpos*u(nion+2,j)
     .                 + twothirds * (wbeam(j) + walp(j))
          else if (iuse .eq. 3) then
              do k=1,nion
                enpos = enpos+usave(k,j)
              end do
              press(j) = enesav(j)*usave(nion+1,j)+enpos*usave(nion+2,j)
     .                 + twothirds * (wbeam(j) + walp(j))
          else
               write  (nout, 1) iuse
               write  (ncrt, 1) iuse
    1          format (' subroutine PRESSR reports:'         /
     .                 ' IUSE not set correctly, IUSE =', i5 /
     .                 ' ONETWO is terminated'                )
               call STOP ('subroutine PRESSR: unspecified problem', 30)
          end if
          press (j) = convert *            press(j)
          pressb(j) = convert * twothirds * wbeam(j)
          press_alpha(j) = convert * twothirds * walp(j)
          press_thermal(j)=press(j)-pressb(j)-press_alpha(j)
      end do
c
c     enforce zero gradient at mag axis (wbeam causes deviation from zero gradient):
       press(1)  = press(2) - (r(2)**2)*(press(3)-press(2))/
     .                                    (r(3)**2 -r(2)**2)
        
      return
c
      end

      subroutine psirho (iunit)
c ----------------------------------------------------------------------
c this subroutine calculates psir, the array of psi values on the rho mesh
c ----------------------------------------------------------------------
c
      USE param
      USE solcon
      USE soln
      USE mhdpar 
      USE numbrs       
      USE mesh
      USE machin
      USE geom
      USE rhog
      USE mhdcom
      USE etc
      USE psig,                        ONLY : psival
      USE gpsi
      implicit  integer (i-n), real*8 (a-h, o-z)


c
c ----------------------------------------------------------------------
c calculate bp in kgauss (if called from SOURCE, iunit = 0)
c calculate bp in tesla  (if called from    MHD, iunit = 1)
c ----------------------------------------------------------------------
c
      convert = 1000.0
      if (iunit .eq. 1)  convert = 1.0
c
c ----------------------------------------------------------------------
c first get bp (i.e., bp0) from rbp (which was found in transport end
c of the code). then integrate the equation bp0 = (1.0/R0)dpsi/drho to
c get psi(rho) ( = psir). the constant of integration is pmin,the value
c of psi at the magnetic axis.
c ----------------------------------------------------------------------
c
      do 10 j=2,nj
   10 bp(j) = rbp(j)/(convert*r(j)*gcap(j)*hcap(j)*fcap(j)) !bp in Tesla
c
      call trap2 (r, bp, psir, nj)        ! returns psir(1) = 0.0
      fudge = 1.0
c
      if (mhdmethd .eq.'tdem') then
c
c       if bp is evolving consistently with the eqdsk
c       poloidal flux values  then fudge=1.0
c
        delta_psi=(psir(nj)-psir(1))*rmajor
        fudge=(pmax-pmin)/delta_psi
      end if
c
      do 20 j=1,nj
   20 psir(j) = fudge*psir(j)*rmajor + pmin  ! psir(1) is set to pmin
      !keep psival(1) in sync:
      psiedge = psir(nj)
      IF(psival(2) .lt. psival(1))THEN !psival is decreasing from psilim towards mag axis
           IF(psival(2) .lt. psiedge)psival(1) = psiedge
      ELSE !if psival is increasing from psilim to mag axis
         IF(psival(2) .GT. psiedge)psival(1) = psiedge
      ENDIF
c      psival(1) = psiedge


      if(time .le. time0)then   ! save the initial map for interpolation
         call copya(r,r_at_t0,nj)
         call copya(psir,psir_at_t0,nj)
      endif
      return
c
      end

      subroutine psiset (npsi, psiaxis, psibdry, psival, iunfrm)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c calculate a vector of psi values, psival. note that
c psival(1) = psi on plasma edge,psival(npsi) = psi at magnetic axis
c ----------------------------------------------------------------------
c
      dimension psival(*)
c
      if (iunfrm .eq. 0) then
        do j=2,npsi-1
          psival(j) = psiaxis-(psiaxis-psibdry)*((npsi-j)/(npsi-1.0))**2
        end do
      else
        dpsi = (psibdry-psiaxis)/(npsi-1)
        do j=2,npsi-1
          psival(j) = psibdry-(j-1)*dpsi
        end do
      end if
      psival(1)    = psibdry
      psival(npsi) = psiaxis
      return
c
      end

      subroutine readit
c ----------------------------------------------------------------------
c --- reads solution information from previous MHD calculation
c --- written by subroutine SAVIT
c ----------------------------------------------------------------------
c
      USE param
      USE fusion
      USE io 
      USE solcon
      USE soln
      USE nub2
      USE mhdpar        
      USE mesh
      USE mixcom
      USE rhog
      USE mhdcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c     include 'param.i'
c      include 'fusion.i'
c      include 'io.i'
c      include 'mesh.i'
c      include 'nub2.i'
c      include 'solcon.i'
c      include 'soln.i'
c      include 'mhdpar.i'
c      include 'mhdcom.i'
c      include 'rhog.i'
c      include 'mixcom.i'
      logical   opened,exists
c
      write  (nitre, 1) psibdry,psiaxis,rma,zma,psi(1,1),psi(17,33),
     .                  psir(1),psir(51)
    1 format (' subroutine READIT: psibdry, psiaxis, rma, zma' /
     .          2x, 8(2x, 1pe14.7))
c
      inquire (unit = nsavsol, opened = opened,exist = exists)
      if (opened) then
          rewind (unit = nsavsol)
      else if(exists)then
          print *,'file savsol exists but is not open'
          call STOP ('subroutine READIT: file is not open', 95)
      end if
     
      read (nsavsol) u
      read (nsavsol) r
      read (nsavsol) time0, time, dt, n, eqtim0
      read (nsavsol) enbsav
      read (nsavsol) wbsav
      read (nsavsol) enasav
      read (nsavsol) wasav
      read (nsavsol) psi
      read (nsavsol) psir
      read (nsavsol) psibdry, psiaxis, rma, zma
c
c --- subroutine mix parameters
c
      read (nsavsol) dtecal, dtemix, dtical, dtimix, epste, epsti
      read (nsavsol) fuscal, fusmix, ipmix, mixpro, qmix, rmixx
      read (nsavsol) rsmixx, rsx, s3cal, s3mix, s71cal, s71mix
      read (nsavsol) s18cal, s18mix
      read (nsavsol) tem
      read (nsavsol) tep
      read (nsavsol) tim
      read (nsavsol) tip
      read (nsavsol) timmix, trcal, trmix, tsmix
      read (nsavsol) tdmix
      read (nsavsol) wmix, w0mix, w1mix, w2mix, w3mix, w4mix, w5mix
      write (nitre, 1) psibdry, psiaxis, rma, zma, psi(1,1), psi(17,33),
     .                 psir(1), psir(51)
      return
c
      end

      subroutine refinept (rcontr, zcontr, ncontr, psiwant, cspln,
     .                     n2cspln, nh2, nw, nh, rmhdgrid, zmhdgrid,
     .                     rma, zma)
c
 
c
c ----------------------------------------------------------------------
c    Given a set of points (rcontr(i),zcontr(i)),which approximately
c    describe the psi contour psi = psiwant,refine these points
c    using the bicubic spline representation of psi:
c    subroutine is identical to the guts of cntour but here
c    it is assumed that an approximate contour allready exists
c    so that the search done in cntour is eliminated.
c --- input
c  psiwant          value of psi to be traced. the contour
c                   points must be reasonably close to the desired
c                   contour so that Newton's method will work.
c  rcontr(i)
c  zcontr(i)
c  ncontr            the original contour points
c  rmhdgrid(nw)
c  zmhdgrid(nh)      the MHD grid
c  cspln(n2cspln,nw,nh2)   the bicubic spline array,must be set on input.
c  rma                     location of the magnetic axis. required so that
c  zma                     we can find a unique point on the contour.
c --- output
c  rcontr(i)
c  zcontr(i)            i = 1..ncontr,the refined contour points.
c                       if refinement of a point is not possible
c                       because Newton's method did not converge then
c                       the corresponding point of the input contour is
c                       returned.
c ------------------------------------------------------------------ HSJ
c
      USE replace_imsl,                            ONLY : my_dbcevl1
      implicit none
      include 'imsl.i'
c
      integer ncontr,iter,j,ier,n2cspln,nh2,nw,nh
      real*8  rcontr(ncontr),zcontr(ncontr),rmhdgrid(nw),zmhdgrid(nh),
     .        cspln(n2cspln,nw,nh2),psiwant,pds(6),sint,cost,xn,yn,
     .        theta,dpsids,dpsi,serr,rma,zma,dtol,drab,dzab,aserr,step,
     .        stepmax
c
      dtol    = rmhdgrid(2)-rmhdgrid(1)+zmhdgrid(2)-zmhdgrid(1)
      stepmax = 0.1 * dtol    ! limit Newton step to this magnitude
c
      imslmd = 'refinept'
c
      do 100 j=1,ncontr
        theta = ATAN2 (zcontr(j)-zma, rcontr(j)-rma)
        sint  =   SIN (theta)
        cost  =   COS (theta)
        xn = rcontr(j)
        yn = zcontr(j)
        iter = 0
  120   call my_dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,
     .  nw,xn,yn,pds,ier,3)
        if (ier .ne. 0)  go to 100
        dpsi = psiwant-pds(1)
        if (ABS (dpsi) .lt. 1.0e-10)  go to 130
        dpsids = pds(2)*cost+pds(3)*sint
        if (dpsids .eq. 0.0)  go to 100
        serr = dpsi/dpsids
        aserr = ABS (serr)
        if (aserr .lt. 1.0e-09)  go to 130
        step =  MIN (aserr, stepmax)
        serr = SIGN (step , serr   )
        xn   = xn+serr*cost
        yn   = yn+serr*sint
        iter = iter + 1
        if (iter .gt. 20)  go to 100
        go to 120
  130   drab = ABS (xn-rcontr(j))  ! required because point found can be
        dzab = ABS (yn-zcontr(j))  ! on opposite side of plasma
        if (drab+dzab .gt. dtol)  go to 100
        rcontr(j) = xn
        zcontr(j) = yn
  100 continue
      return
c
      end

      subroutine reqdsk
c
c
c ----------------------------------------------------------------------
c this subroutine reads an eqdsk as a starting guess for an equilibrium.
c and does the necessary calculations for MHD/transport startup
c it is also called from INIT to pick up limiter points if necessary
c (quantities indexed by nxeqd are read from the eqdsk in reverse
c order so that for example qpsi(1) = edge,qpsi(nxeqd) = axis value.)
c  INCLUDE file psig:
c  presspsi(1..nxeqd)    pressure (nt/m**2)
c  fpsi(1..nxeqd)        f(psi),tesla/m
c  ffppsi(1..nxeqd)      f(psi)*d(f)/dpsi
c  pppsi(1..nxeqd)       d(presspsi)/dpsi
c  qpsi(1..nxeqd)        safety factor
c  pressb(1...nxeqd)     beam pressure (if known)
c
c  on exit from this subroutine all quantities are in guassian units!
c  the switch between gaussian and MKS (and MKS back to gaussian) is
c  done in subroutine MHD.
c  The actual MHD calculations are done in MKS units.
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE fusion
      USE ions
      USE io 
      USE neut
      USE nub
      USE nub2
      USE solcon
      USE soln
      USE contour
      USE limiter
      USE mhdpar  
      USE mhdgrid
      USE rf 
      USE numbrs
      USE extra     
      USE mesh
      USE machin
      USE geom
      USE tordlrot
      USE constnts
      USE soln2d
      USE bd_condtn
      USE psig
      USE rhog
      USE mhdcom
      USE bicube
      USE shapctr
      USE flxav
      USE neo2dp
      USE etc
      USE gpsi 
      USE replace_imsl,                ONLY : my_ibcccu,my_dbcevl1
      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'imsl.i'


      include 'storage.i'

c
      character    ntitle(5)*8
c
      dimension    torflux(kstore), torcur(kstore), bpcontr(kstore)
      dimension spbb(nw)    ! also local in wrteqdsk,not currently used
      equivalence (torflux(1), xdum(1))
      equivalence                  (torcur(1), ydum(1))
      equivalence                                  (bpcontr(1), zdum(1))
c
      logical    eqtype
      integer    LENGTH, asize
      character  eqdsrce*13, tempname*64
      external   LENGTH
c
      imslmd = 'reqdsk'
c
      asize = LENGTH (eqdskin)          ! number of characters in eqdskin
      write (tempname, '(a)') eqdskin(1:asize)
      call myopen (nrguess, tempname, 2, length, iread)
      if (iread .ne. 1) then
        write  (ncrt, 8100)  eqdskin(1:asize)
        write  (nout, 8100)  eqdskin(1:asize)
 8100   format (/ ' subroutine REQDSK reports:'                   /
     .         5x, 'FATAL ERROR: could not read input eqdsk file' /
     .         5x, 'named: "', a, '"'                             /
     .         5x, 'therefore ONETWO cannot continue'              )
        call giveupus(nrguess)
        call STOP ('subroutine REQDSK: problem #1', 31)
      end if
c
      ncontr  = 0
      limread = 0
      read (nrguess, 8190) (ntitle(i), i=1,5),
     .                      dat, ipestg, nxeqd, nyeqd, eqdsrce !line1
      print *,'eqdsrce =',eqdsrce
      if (mhdmode .eq. 'no coils'       .and.
     .    eqdsrce .ne. 'FIXBDRY EQDSK') then 
        write  (nout, 8257)  eqdskin(1:asize)
        write  (ncrt, 8257)  eqdskin(1:asize)
 8257   format (/ ' subroutine REQDSK reports:'                    /
     .         5x, 'fixed boundary equilibria must use a fixed'    /
     .         5x, 'boundary eqdsk as input. The input eqdsk file' /
     .         5x, 'named: "', a, '" is not of this type,'         /
     .             'therefore ONETWO cannot continue'               )
        call giveupus(nrguess)
        call STOP ('subroutine REQDSK: problem #5', 98)
      end if
      if (eqdsrce .eq. 'FIXBDRY EQDSK')  mhdmode = 'no coils'
      if (nxeqd .gt. kpsi)  go to 5000
      read (nrguess, 8200)  xdimeqd, ydimeqd, reqd, redeqd, ymideqd !line2
      if (ABS (ymideqd) .gt. 1.0e-5) then
          write  (nout, 8251)  ymideqd
          write  (ncrt, 8251)  ymideqd
 8251     format (' subroutine REQDSK detects an error:' /
     .            ' ymideqd =', 1pe12.6                  /
     .            ' ONETWO expects ymideqd = 0.0'        /
     .            ' ONETWO must stop'                     )
          call giveupus(nrguess)
          call STOP ('subroutine REQDSK: problem #2', 32)
      end if
c
      eqtype = .false.
      if (eqdsrce .eq.  'ONETWO EQDSK '  .or.
     .    eqdsrce .eq. 'FIXBDRY EQDSK')  eqtype = .true.
c
      contrpts_set = 0    ! set to 1 below if boundary pts are available
      read (nrguess, 8200)  rma, zma, psimag, psilim, beqd       !line3
      read (nrguess, 8200)  toteqd, psimx1, psimx2, xax1, xax2   !line4
c                                   psimx1, psimx2 are not used
      read (nrguess, 8200)  zax1, zax2, psisep, rsep, zsep
      read (nrguess, 8200)  (    fpsi(i), i=nxeqd,1,-1)
c                                               fpsi(nxeqd) = axis value


      read (nrguess, 8200)  (presspsi(i), i=nxeqd,1,-1)
c                                           presspsi(nxeqd) = axis value
      read (nrguess, 8200)  (ffppsi(i), i=nxeqd,1,-1)
      read (nrguess, 8200)  (pppsi(i),  i=nxeqd,1,-1)
      read (nrguess, 8200)  ((psi(i,j), i=1,nxeqd), j=1,nyeqd)
      read (nrguess, 8200, err=8221, end=8221)  (qpsi(i), i=nxeqd,1,-1)
      read (nrguess, 8210, err=8221, end=8221)   ncontr, nlimtr
      if (ncontr .gt. nconmax)  go to 5000
      read (nrguess, 8200, err=8221, end=8221)  (rcontr(i), zcontr(i),
     .                                                i = 1,ncontr)
      if (ncontr .gt. 50)  contrpts_set = 1
      ! save eqdsk quantites:
      IF(.NOT. ALLOCATED( pppsi_eqdsk))ALLOCATE(pppsi_eqdsk(nxeqd))
      IF(.NOT. ALLOCATED( presspsi_eqdsk))
     .                    ALLOCATE(presspsi_eqdsk(nxeqd))
      IF(.NOT. ALLOCATED( fpsi_eqdsk))ALLOCATE(fpsi_eqdsk(nxeqd))
      IF(.NOT. ALLOCATED( ffppsi_eqdsk))ALLOCATE(ffppsi_eqdsk(nxeqd))
      IF(.NOT. ALLOCATED( qpsi_eqdsk))ALLOCATE(qpsi_eqdsk(nxeqd))
      pppsi_eqdsk(1:nxeqd)  = pppsi(1:nxeqd)
      presspsi_eqdsk(1:nxeqd)  = presspsi(1:nxeqd)
      fpsi_eqdsk(1:nxeqd)   = fpsi(1:nxeqd)
      ffppsi_eqdsk(1:nxeqd) = ffppsi(1:nxeqd)
      qpsi_eqdsk(1:nxeqd)   = qpsi(1:nxeqd)
      nxeqd_eqdsk           = nxeqd
      psimag_eqdsk          = psimag
      psilim_eqdsk          = psilim

c
c     to match TRANSP results:
c
      pol_flux = psilim - psimag
      if (pol_flux_lim .lt. 1.0) then
        psilim = psimag + pol_flux_lim * pol_flux
      end if
c
c --- load the plasma boundary vectors
c
      do i=1,ncontr
        rplasbdry(i) = rcontr(i)
        zplasbdry(i) = zcontr(i)
      end do
      nplasbdry = ncontr
      if (nlimiter .gt. 0) then
      call volcalc (rcontr,zcontr,ncontr,rma,zma,eqdskvol,eqdskarea)
      write (ncrt, '(1x, a, f12.6)')
     .      'volume determined from eqdsk boundary values = ', eqdskvol
      write (ncrt, '(1x, a, f12.6)')
     .      'area   determined from eqdsk boundary values = ', eqdskarea
      end if
c
c --- pick up the limiter if nlimiter is non positive. (which will only
c --- happen if REQDSK is called from INIT.) Thereafter nlimiter is a
c --- positive number so the following will not be executed. Thus, even
c --- if multiple eqdsks are read at different solution times,
c --- the limiter will not be changed.
c
 8221 if (nlimiter .le. 0) then
        if (nlimtr .gt. maxlimpt-2)  go to 5000
        nlimiter = nlimtr                        ! read from eqdsk above
        limread  = 1
        read (nrguess, 8200, err = 8220, end = 8230)
     .       (xlimiter(i), ylimiter(i), i=1,nlimiter)
        go to 4000
 8220   write  (nout, 8240)  eqdskin(1:asize)
        write  (ncrt, 8240)  eqdskin(1:asize)
 8240   format (/ ' ERROR: eqdsk file "', a, '"'                    /
     .                8x, 'has generated an error in while reading' /
     .                8x, 'the limiter points in subroutine REQDSK' /)
        nlimiter = -1
        go to 4000
 8230   write  (nout, 8250)  eqdskin(1:asize)
        write  (ncrt, 8250)  eqdskin(1:asize)
 8250   format (/ ' ERROR: eqdsk file "', a, '"'            /
     .                8x, 'does not contain limiter points' /)
        nlimiter = -1
        go to 4000
      end if
c
 8190 format (6a8, 3i4, t73, a)
 8200 format (5e16.9)
 8210 format (2i5)
 8235 format (4(2x, i5))
c
c ----------------------------------------------------------------------
c --- the following output is not part of the "official" eqdsk. It is
c --- what distinguishes eqdsk created in ONETWO from those created by EFIT:
c --- we put out enough information so that the eqdsk could be used as
c --- a restart file.
c ----------------------------------------------------------------------
c
      if (eqdsrce .eq. 'ONETWO EQDSK' .and. irguess .eq. -5) then
c
c --- skip over the limiter points (they were read above
c --- in the call to REQDSK from INIT if they were required,
c --- so don't change the logic)
c
      read (nrguess, 8200)  (wdum(1),wdum(1), i=1,nlimtr)
      read (nrguess, 8235)   nj, nprim, nimp, nti,npsil
      nion = nprim + nimp
      read (nrguess, 8200)  (atw(k), k=1,nion)
      do 760 k=1,nion
  760 read (nrguess, 8200)  (z     (j,k), j=1,nj)
      read (nrguess, 8200)  (zeff  (j  ), j=1,nj)
      read (nrguess, 8200)  (r     (j  ), j=1,nj)
      read (nrguess, 8200)  (psir  (j  ), j=1,nj)
      read (nrguess, 8200)  (curden(j  ), j=1,nj)
      read (nrguess, 8200)  (te    (j  ), j=1,nj)
      do 770 k=1,nti
  770 read (nrguess, 8200)  (ti    (j  ), j=1,nj)
      read (nrguess, 8200)  (ene   (j  ), j=1,nj)
      do 730 k=1,nprim
  730 read (nrguess, 8200)  (en    (j,k), j=1,nj)
      if (nimp .ne. 0)   then
        do 740 k=nprim+1,nion
  740   read (nrguess,8200) (en    (j,k), j=1,nj)
      end if
      read (nrguess, 8200) (spbb(i), i=nxeqd,1,-1)
c
      read (nrguess, 8200) (pressb(i), i=1,nj)
      read (nrguess, 8200) ((u(i,j), i=1,nion+4), j=1,nj)
      read (nrguess, 8200) time0,time,dt,realn,eqtim0,dtt,realibeam
      read (nrguess, 8200) (fcap(j),j=1,nj)
      read (nrguess, 8200) (gcap(j),j=1,nj)
      read (nrguess, 8200) (hcap(j),j=1,nj)
      read (nrguess, 8200) (r2cap(j),j=1,nj)
      read (nrguess, 8200) (bp(j),j=1,nj)
      read (nrguess, 8200) (pprim(j),j=1,nj)
      read (nrguess, 8200) (ffprim(j),j=1,nj)
      read (nrguess, 8200) (psival(j),j=1,npsil)
      read (nrguess, 8200) (rho(j),j=1,npsil)
      read (nrguess, 8200) (q(j),j=1,nj)
      read (nrguess, 8200) (r2capi(j),j=1,nj)
      read (nrguess, 8200) (rcap(j),j=1,nj)
      read (nrguess, 8200) (((enbsav(i,j,k),i=1,nj),j=1,ke),k=1,kb)
      read (nrguess, 8200) (((wbsav(i,j,k),i=1,nj),j=1,ke),k=1,kb)
      read (nrguess, 8200) (enasav(i),i=1,nj)
      read (nrguess, 8200) (wasav(i),i=1,nj)
      read (nrguess, 8200) (rcapi(j),j=1,nj)
c
      n     = realn           ! number of transport time steps completed
      ibeam = realibeam
      end if
c
c ----------------------------------------------------------------------
c convert eqdsk to gauss cm units
c ----------------------------------------------------------------------
c
 4000 xdimeqd = 100.0 * xdimeqd
      ydimeqd = 100.0 * ydimeqd
      reqd    = 100.0 * reqd
      redeqd  = 100.0 * redeqd
      reqdsk_box_edge = redeqd
      ymideqd = 100.0 * ymideqd
      if (limread .eq. 1)  go to 3435
      rma     = 100.0 * rma
      zma     = 100.0 * zma
      xax(1)  = rma
      yax(1)  = zma
      rsep    = 100.0 * rsep
      zsep    = 100.0 * zsep
      psimag  = psikgaus*1.0e3*psimag ! gauss-cm**2, psikgaus*1000=1.0e8
      psilim  = psikgaus*1.0e3*psilim
      call multpl1 (psi, nwh, psikgaus*1.0e3)
      beqd    = 1.0e4 * beqd                       ! gauss
      cconst  = 1.0e+6
      call multpl1 (fpsi     , nxeqd    , cconst) ! fpsi in gauss-cm
      cconst  = 1.0e+4
      call multpl1 (ffppsi   , nxeqd    , cconst) ! ffppsi in gauss
      cconst  = 10.0
      call multpl1 (presspsi , nxeqd    , cconst) ! presspsi: erg/cm**3
      cconst  = 1.0e-7
      call multpl1 (pppsi    , nxeqd    , cconst) ! pppsi in..
c                                              ..gram/gauss*cm**3*sec**2
      cconst  = 100.0
      call multpl1 (rcontr   , ncontr   , cconst) ! cm
      call multpl1 (zcontr   , ncontr   , cconst) ! cm
      call multpl1 (rplasbdry, nplasbdry, cconst) ! cm
      call multpl1 (zplasbdry, nplasbdry, cconst) ! cm
c
c ----------------------------------------------------------------------
c get values for rmajor, btor, and totcur from eqdsk if values were
c input as zero and irguess < 0 and this is initial equilibrium (iteq = 0)
c and time-dependent equilibria are not done (ifixshap = 1)
c ----------------------------------------------------------------------
c
      if (irguess .lt. 0 .and. ieq .eq. 0 .and. ifixshap .eq. 1) then
            if (rmajor .eq. 0.0)  rmajor = reqd
            reqd = rmajor
            if (btor .eq. 0.0)  btor = beqd
            beqd = btor
            flim = rmajor*btor
            if (totcur(1) .eq. 0.0)  totcur(1) = toteqd
            curfac  = totcur(1)/toteqd
            toteqd  = totcur(1) 
            rbp(nj) = 0.2 * ABS (toteqd)
            if (time .eq. bctime(1))
     .        bc(1,nk-iangrot) = 0.2 * ABS (toteqd)
            psimag  = curfac*psimag
            psilim  = curfac*psilim
            do 700 j=1,nyeqd
            do 700 i=1,nxeqd
  700         psi(i,j) = psi(i,j)*curfac
            do j=1,nxeqd
              ffppsi(j) = ffppsi(j)*curfac*curfac
              pppsi(j)  = pppsi(j)*curfac*curfac
            end do
            do j=2,nxeqd    ! j = 1 is plasma edge
              presspsi(j) = curfac * curfac * presspsi(j)
            end do
            iunfrm = 1
c
c           set up psival,with psival(1) = edge, psival(nxeqd) = axis
c
            call psiset (nxeqd,psimag,psilim,psival,iunfrm)
c
c           integrate ffprime from edge to axis
c
            call trap2 (psival, ffppsi, xdum, nxeqd)
            isgnfpsi = 1
            if (flim .lt. 0.0)  isgnfpsi = -1
            flimsq = flim * flim
            do j=1,nxeqd
              fpsi(j) = isgnfpsi * SQRT (flimsq+2.0*xdum(j))
            end do
      end if

c
      if (     btor .eq. 0.0)       btor = beqd
      if (totcur(1) .eq. 0.0)  totcur(1) = toteqd
      if (   rmajor .eq. 0.0)     rmajor = reqd
      if (rmajor*btor*totcur(1) .eq. 0.0) then
          write  (nout, 735)  btor, totcur(1), rmajor
          write  (ncrt, 735)  btor, totcur(1), rmajor
  735     format (/ ' subroutine REQDSK detected an ERROR' /
     .              ' btor      =', 1pe12.6                /
     .              ' totcur(1) =', 1pe12.6                /
     .              ' rmajor    =', 1pe12.6                /
     .              ' ONETWO will be terminated')
      end if
c
c ----------------------------------------------------------------------
c Check the eqdsk (with possible modifications done above) for
c self-consistency and adjust transport coordinate system as appropriate
c ----------------------------------------------------------------------
c

      call chekeqd (psi, nyeqd, xdimeqd, ydimeqd, redeqd, ncrt,
     .              rcontr, zcontr, ncontr, btor, rma, zma,
     .              psimag, psilim, totcur(1), rmajor, nconmax)
c
      psiaxis = psimag
      psibdry = psilim
c
c --- set up the psival grid: psival(1) = edge, psival(nw) = axis
c
      iunfrm = 1
      call psiset (nw, psiaxis, psibdry, psival, iunfrm)
c
c --- save eqdsk pressure and psi for plotting (subroutine EQPLOT)
c
      do j=1,nxeqd
        psieqdsk(j) = psival  (nxeqd-j+1) * 1.0e-08    ! volt-sec
        prseqdsk(j) = presspsi(nxeqd-j+1) * 0.1        ! nt/m**2
      end do
c
c --- interpolate the eqdsk onto the (nw,nh) grid used in ONETWO if necessary
c
      if (nxeqd .gt. nw)  go to 6000
      if (nyeqd .gt. nh)  go to 6000
      if (nxeqd .ne. nw .or. nyeqd .ne. nh) then
c
c --- be sure that new eqdsk gets written for use in rf codes,etc:
c --- ieqdsk .ne. 0 will force a call to wrteqdsk
c --- and provide an eqdskname for use in rf codes
         if(ieqdsk .eq. 0)ieqdsk = -1
       
c --- next  get the grid that corresponds to the eqdsk
c
        dxneqd      =  xdimeqd/(nxeqd-1)
        dyneqd      =  ydimeqd/(nyeqd-1)
c
        do i=1,nxeqd
          xdum(i  ) =  redeqd + (i-1)*dxneqd
        end do
c
        xdum(1    ) =  rmhdgrid(1)
        xdum(nxeqd) =  rmhdgrid(nw)
c
        do j=1,nyeqd
          ydum(j  ) = -ydimeqd * 0.5 + (j-1)*dyneqd
        end do
c
        ydum(1    ) =  zmhdgrid(1)
        ydum(nyeqd) =  zmhdgrid(nh)
c
c --- get psi on (nw,nh) grid
c
        do 30 j=1,nyeqd
        do 30 i=1,nxeqd
   30   p(i,j) = psi(i,j)
        call my_ibcccu (p,xdum,nxeqd,ydum,nyeqd,cspln,nw,wnoperm,ier)
        do 35 j=1,nh
        do 35 i=1,nw
          call my_dbcevl1 (xdum,nxeqd,ydum,nyeqd,cspln,nw,rmhdgrid(i),
     .                  zmhdgrid(j),pds,ier,1)
   35  psi(i,j) = pds(1)
c
c --- next define the eqdsk psi grid
c
        iunfrm = 1
        call psiset (nxeqd,psiaxis,psibdry,wdum,iunfrm)
c
c --- interpolate fpsi
c
        do 40 i=1,nxeqd
   40   xdum(i) = fpsi(i)
        call intrp (1, 1, wdum, xdum, nxeqd, psival, fpsi, nw)
c
c --- interpolate presspsi
c
        do 45 i=1,nxeqd
   45   xdum(i) = presspsi(i)
        call intrp (1, 1, wdum, xdum, nxeqd, psival, presspsi, nw)
c
c --- interpolate ffppsi
c
        do 50 i=1,nxeqd
   50   xdum(i) = ffppsi(i)
        call intrp (1, 1, wdum, xdum, nxeqd, psival, ffppsi, nw)
c
c --- interpolate pppsi
c
        do 55 i=1,nxeqd
   55   xdum(i) = pppsi(i)
        call intrp (1, 1, wdum, xdum, nxeqd, psival, pppsi, nw)
c
c --- interpolate qpsi
c
        do 60 i=1,nxeqd
   60   xdum(i) = qpsi(i)
        call intrp (1, 1, wdum, xdum, nxeqd, psival, qpsi, nw)
c
c --- interpolate spbb
c
        do 65 i=1,nxeqd
   65   xdum(i) = spbb(i)
        call intrp (1, 1, wdum, xdum, nxeqd, psival, spbb, nw)
      end if    ! end eqdsk interpolation
c
c --- now we have everything set on the (nw,nh) MHD grid used in ONETWO
c
c ----------------------------------------------------------------------
c --- calculate the toroidal flux
c ----------------------------------------------------------------------
c
      call trap2 (psival, qpsi, torflux, nw)
      do j=1,nw
        torflux(j) = twopi * (torflux(j) - torflux(nw))
      end do
c
c ----------------------------------------------------------------------
c --- define the rho grid
c ----------------------------------------------------------------------
c
      do 160 i=1,nw
  160   rho(i) = SQRT (torflux(i)/(pi*btor))
c
c ----------------------------------------------------------------------
c --- define the r grid
c ----------------------------------------------------------------------
c
      do 170 i=1,nj
  170 r(i) = r(i) * rho(1) / r(nj)
c
c ----------------------------------------------------------------------
c --- get q on r grid
c ----------------------------------------------------------------------
c
  605 call intrp (0, 1, rho, qpsi, nw, r, q, nj)
c
c ----------------------------------------------------------------------
c --- calculate bp0
c ----------------------------------------------------------------------
c
      do 600 j=1,nj
  600 bp(j) = r(j)*btor/(rmajor*q(j))
c
c ----------------------------------------------------------------------
c --- now define the psir grid(this is the psi grid corresponding to
c --- the r grid. it is defined by eq. 2.1-8)
c ----------------------------------------------------------------------
c
      call trap2 (r, bp, psir, nj)
      do 610 j=1,nj
  610 psir(j) = psival(nw) + (psival(1) - psival(nw))*psir(j)/psir(nj)
c
c ----------------------------------------------------------------------
c --- get pressure on psir grid
c ----------------------------------------------------------------------
c
      call intrp (1, 1, psival, presspsi, nw, psir, press, nj)
c
c ----------------------------------------------------------------------
c --- get pprim on psir grid
c ----------------------------------------------------------------------
c
      call intrp (1, 1, psival, pppsi, nw, psir, pprim, nj)
****  call difydx(psir,press,pprim,nj)
c
c ----------------------------------------------------------------------
c --- get ffprim on psir grid
c ----------------------------------------------------------------------
c
      call intrp (1, 1, psival, ffppsi, nw, psir, ffprim, nj)
c
c ----------------------------------------------------------------------
c --- get qpsir on psir grid
c ----------------------------------------------------------------------
c
      call intrp (1, 1, psival, qpsi, nw, psir, qpsir, nj)
c
c ----------------------------------------------------------------------
c --- check to make sure user has set non-DIII-D geometry correctly,
c --- if applicable
c ----------------------------------------------------------------------
c
      if ((ABS (rmhdgrid(1)-redeqd        ) .gt. 1.0e-3)  .or.
     .    (ABS (zmhdgrid(1)+ydimeqd * 0.5)) .gt. 1.0e-3) then
        write  (nout, 736) rmhdgrid(1), redeqd, zmhdgrid(1), 0.5*ydimeqd
        write  (ncrt, 736) rmhdgrid(1), redeqd, zmhdgrid(1), 0.5*ydimeqd
  736   format ('  ERROR: Non-DIII-D geometry requires that ',
     .            'XDIM, YDIM, REDGE and NLIMITER = -1'        /
     .            'be set in the third namelist of inone file' /
     .            'RMAJOR and RMINOR must be set in first namelist ',
     .            'of inone (all values are in cm)'            /
     .            'Current values are:'                        //
     .            '  rmhdgrid(1)   = ', 1pe14.6                /
     .            '  redeqd        = ', 1pe14.6                /
     .            '  zmhdgrid(1)   = ', 1pe14.6                /
     .            '  0.5 * ydimeqd = ', 1pe14.6)
        call giveupus(nrguess)
        call STOP ('subroutine REQDSK: geometry incompatibility', 99)
      end if
c
c ----------------------------------------------------------------------
c --- next get the plasma contours so we can form curden( = <jtor*R0/R>)
c --- on psir grid
c --- all contours to be found must be inside the box defined by
c --- xmin, xmax, ymin, ymax
c ----------------------------------------------------------------------
c
        xmin   = xlimiter(nlimiter+1) * 100.0  ! need cm in CNTOUR below
        xmax   = xlimiter(nlimiter+2) * 100.0
        ymin   = ylimiter(nlimiter+1) * 100.0
        ymax   = ylimiter(nlimiter+2) * 100.0
        bperr  = 0.05
c        dx     = drmhdgrd * 0.5
c        dy     = dzmhdgrd * 0.5
        arcl   =  2.0           ! 2 cm arcl length
        taxis  =  5.0
        tlim   = 30.0
        a      = (tlim-taxis) / (psir(nj)-psir(1))
        bincp  = taxis
        iounit = ncrt
c
c --- set up for bicubic spline representation of psi.
c --- CNTOUR works with psi max on axis, so temporarily invert
c
      cconst = -1.0
      call multpl1 (psi, nwh, cconst)
      call my_ibcccu  (psi,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                 wnoperm,ierr)
      call multpl1 (psi, nwh, cconst)
c
c --- for interior flux surfaces 2 to nj-1
c
      delta_psi=-psir(nj-1)+psir(nj)
      do j=2,nj-1
        iauto    = 0
        dang     = a*(psir(j)-psir(1))+bincp
        psivalue = -psir(j)
        dx0      = 0.0
        dy0      = 0.0
        call cntour (rma,zma,psivalue,rcmin,rcmax,zcmin,zcmax,
     .               zrcmin,zrcmax,rzcmin,rzcmax,dang,arcl,bperr,
     .               dx0,dy0,xmin,xmax,ymin,ymax,iauto,iautoc,rcontr,
     .               zcontr,ncontr,rmhdgrid,nw,zmhdgrid,nh,
     .               cspln,n2cspln,nh2,iounit,nconmax,ierr,bpcontr,0,
     .               delta_psi)
c
c --- get curden = <jtor*R0/R> and r2cap = <R0**2/R**2>
c
        sum  = 0.0
        sum1 = 0.0
        sum2 = 0.0
        torcur(1) =
     .       -10.0 * (rcontr(1)*pprim(j)+ffprim(j) / (fourpi*rcontr(1)))
        do i=2,ncontr 
          torcur(i) = -10.0*(rcontr(i)*pprim(j)+ffprim(j) /
     .                      (fourpi*rcontr(i)))              ! amp/cm**2
          dl   = SQRT ((rcontr(i)-rcontr(i-1))**2
     .               + (zcontr(i)-zcontr(i-1))**2)
          sum1 = sum1 + dl * (1.0 / bpcontr(i) + 1.0 / bpcontr(i-1))
          sum2 = sum2 +
     .           dl * (torcur(i  ) / (bpcontr(i  ) * rcontr(i  ))
     .               + torcur(i-1) / (bpcontr(i-1) * rcontr(i-1)))
          sum  = sum +
     .           dl * (1.0 / (rcontr(i  )**2 * bpcontr(i  ))
     .               + 1.0 / (rcontr(i-1)**2 * bpcontr(i-1)))
        end do
        if (sum1 .eq. 0.0) print *, '>> (line 10920) SUM1 = 0.0' ! DEBUG
        r2cap (j) = sum  * rmajor**2 / sum1
        curden(j) = sum2 * rmajor    / sum1
        if (j .eq. nj-1)
     .    call volcalc (rcontr,zcontr,ncontr,rma,zma,eqdskvolnjm1,
     .                  eqdskareanjm1)
      end do    ! end loop over contours j=2,....nj-1
c
c --- for outermost flux surface, j = nj, we retrace the contour only if
c --- the eqdsk is from EFIT (because EFIT does not give a sufficient
c --- number of points on the boundary)
c
      j = nj
      if (eqtype) then
        call fixedcntour (rplasbdry, zplasbdry, nplasbdry,
     .                    rcontr, zcontr, ncontr,
     .                    rcmin, rcmax, zcmin, zcmax,
     .                    rzcmin, rzcmax, zrcmin, zrcmax,
     .                    rmhdgrid, zmhdgrid, nw, nh,
     .                    bpcontr, cspln, n2cspln, nh2, pds)
      else  ! retrace the contour for EFIT-type eqdsk if use_efit_cntr=0
        if (use_efit_cntr .eq. 0) then
           iauto    = 1
           dang     = a*(psir(j)-psir(1))+bincp
           psivalue = -psir(j)
           dx0      = 0.0
           dy0      = 0.0
           call cntour (rma,zma,psivalue,rcmin,rcmax,zcmin,zcmax,
     .                  zrcmin,zrcmax,rzcmin,rzcmax,dang,arcl,bperr,
     .                  dx0,dy0,xmin,xmax,ymin,ymax,iauto,iautoc,
     .                  rcontr,zcontr,ncontr,rmhdgrid,nw,zmhdgrid,nh,
     .                  cspln,n2cspln,nh2,iounit,nconmax,ierr,bpcontr,
     .                  0,delta_psi)
           if (ierr .eq. 0) then
             call copya (rcontr,rplasbdry,ncontr)
             call copya (zcontr,zplasbdry,ncontr)
             nplasbdry = ncontr
           else                                 ! couldn't get boundary
               if (contrpts_set .eq. 1) then    ! use eqdsk boundary
                   ierr = 0                     ! if possible
                   call set_bdry(rcmin,rcmax,zcmin,zcmax,zrcmin,
     .                           zrcmax,rzcmin,rzcmax,bpcontr,
     .                           rmhdgrid,zmhdgrid,nw,nh,cspln,
     .                           n2cspln,nh2,ncrt)
                   write (nout,'(" CNTOUR did not find suitable plasma"
     .                           " boundary (called from reqdsk)." /
     .                           " Will use plasma boundary",
     .                           " directly from eqdsk")')
                   write (ncrt,'(" CNTOUR did not find suitable plasma"
     .                           " boundary (called from reqdsk)." /
     .                           " Will use plasma boundary",
     .                           " directly from eqdsk")')
               else
                   write (nout,'(" Failed to find the plasma boundary" /
     .                        " Also could not get boundary from eqdsk",
     .                        " because values were not present")')
                   write (ncrt,'(" Failed to find the plasma boundary" /
     .                        " Also could not get boundary from eqdsk",
     .                        " because values were not present")')
        call giveupus(nrguess)
        call STOP ('subroutine REQDSK: eqdsk boundary problem #1', 217)
               end if
           end if
        else    ! user requested boundary points from the eqdsk directly
               if (contrpts_set .eq. 1) then        ! use eqdsk boundary
                  call set_bdry (rcmin,rcmax,zcmin,zcmax,zrcmin,
     .                           zrcmax,rzcmin,rzcmax,bpcontr,
     .                           rmhdgrid,zmhdgrid,nw,nh,cspln,
     .                           n2cspln,nh2,ncrt)
c
c              The eqdsk contour may also be used in FLUXAV
c
               else
                 write (nout, '(" Could not get boundary from eqdsk" /
     .                          " because values were not present"   /
     .                          " try setting use_efit_cntr=0",
     .                          " in third namelist of inone")')
        call giveupus(nrguess)
        call STOP ('subroutine REQDSK: eqdsk boundary problem #2', 218)
               end if
        end if
      end if
c
      rplasmax = rcmax
      rplasmin = rcmin
      zplasmax = zcmax
      zplasmin = zcmin
c
c --- get curden = <jtor*R0/R> and r2cap = <R0**2/R**2>
c
      sum       =   0.0
      sum1      =   0.0
      sum2      =   0.0
      torcur(1) = -10.0 *
     .            (rcontr(1)*pprim(j)+ffprim(j)/(fourpi*rcontr(1)))
c
      do i=2,ncontr
        torcur(i) = -10.0 * (rcontr(i)*pprim(j)+ffprim(j) /
     .                      (fourpi*rcontr(i)))              ! amp/cm**2
        dl = SQRT ((rcontr(i)-rcontr(i-1))**2
     .           + (zcontr(i)-zcontr(i-1))**2)
        sum1 = sum1 + dl * (1.0/bpcontr(i)+1.0/bpcontr(i-1))
        sum2 = sum2 + dl * (torcur(i)/(bpcontr(i)*rcontr(i)) +
     .         torcur(i-1)/(bpcontr(i-1)*rcontr(i-1)))
        sum  = sum  + dl * (1.0 / (rcontr(i)**2*bpcontr(i)) +
     .         1.0 / (rcontr(i-1)**2*bpcontr(i-1)))
      end do
c
      if (sum1 .eq. 0.0) print *, '>> (line 10864) SUM1 = 0.0'   ! DEBUG
      if (use_efit_cntr .eq. 1) then
        dvdpsi   = (eqdskvol*1.e6-eqdskvolnjm1)/(psir(nj)-psir(nj-1))
        r2cap(j) = twopi*twopi*rmajor*qpsir(nj) / (btor*dvdpsi)
      else
        r2cap(j) = sum  * rmajor**2 / sum1
      end if
      curden(j) = sum2 * rmajor    / sum1
c
c --- finally, at the magnetic axis
c
      r2cap (1) = rmajor**2 / rma**2
      curden(1) = -10.0 * (rma * pprim(1) + ffprim(1) /
     .                    (fourpi * rma)) * (rmajor / rma)
c
c --- try a simpler approach
c
      do 185 j=1,nj
  185 curden(j) = -10.0 * (rmajor * pprim(j) + ffprim(j) * r2cap(j) /
     .                    (fourpi * rmajor))                ! amps/cm**2

            
c
c --- get f on psir grid
c
      call intrp (1, 1, psival, fpsi, nw, psir, ydum, nj)
c
c --- set boundary conditions
c
      flim = rmajor * btor
      if (time .eq. bctime(1))  bc(1,nk-iangrot) = 0.2 * totcur(1)
      rbp(nj) = 0.2*totcur(1)
      do 210 j=1,nj
  210    hcap(j) = flim/(ydum(j)*r2cap(j))
c
c --- next normalize curden
c
      do 230 j=1,nj
  230   xdum(j) = twopi*r(j)*hcap(j)*curden(j)
      call trap2 (r, xdum, ydum, nj)
      curnorml = totcur(1) / ydum(nj)
      do 240 j=1,nj
  240   curden(j) = curnorml * curden(j)
       


c
c ----------------------------------------------------------------------
c --- we have curden from eqdsk pprim and ffprim. For ABS (irguess) = 2
c --- we want to use pprim as determined from input density and temps.
c --- So we now form the pressure accordingly,differentiate to get a new
c --- pprim and then form a new ffprim consistent with given curden and pprim
c ----------------------------------------------------------------------
c
      if (IABS (irguess) .eq. 2) then
        call pressr (0, 1)
        call difydx (psir, press, pprim, nj)
        do j=1,nj
          ffprim(j) = -(0.1 * curden(j) + rmajor * pprim(j))
     .                      * fourpi * rmajor / r2cap(j)
        end do

      end if
c
c ----------------------------------------------------------------------
c --- convert psi(i,j),p(i,j),psival(i),psir(j) to kgauss cm**2
c --- convert bp(j) to kgauss
c ----------------------------------------------------------------------
c
      cconst = 1.0e-3
      call multpl1 (bp, nj, cconst)
      psibdry = psilim  * 1.0e-3
      psiaxis = psiaxis * 1.0e-3
      call multpl1 (psival, nw , cconst)
      call multpl1 (psir  , nj , cconst)
      call multpl1 (psi   , nwh, cconst)
      call copya   (psi   , p  , nwh   ) ! p must be set for irguess < 0
c
c --- SAVCUR copies pprim into ppold, ffprim into ffpold,
c --- curden into curold, psir into psiold, and r into rold,
c --- so that subroutine CONCUR can do its thing
c
      call savcur
      go to 3436
c
c --- this section of code is executed only if REQDSK was called by INIT
c --- to get the limiter points.
c --- If the mhdgrid is not set then set it from the eqdsk parameters, in cm
c
 3435 if (rmhdgrid(nw) .lt. -0.9e+30) then
      drmhdgrd = xdimeqd/(nw-1)
      do 100 i=1,nw
  100   rmhdgrid(i) = redeqd+(i-1)*drmhdgrd
      end if
      if (zmhdgrid(nh) .lt. -0.9e+30) then
      dzmhdgrd = ydimeqd/(nh-1)
      do 110 j=1,nh
  110   zmhdgrid(j) = -ydimeqd * 0.5+(j-1)*dzmhdgrd
      end if
c
c ----------------------------------------------------------------------
c close guess file
c ----------------------------------------------------------------------
c
 3436 call giveupus(nrguess)
      close (unit = nrguess)



c -------------------------------------------------------------------------
c if the total current boundary condition is not set for times
c other than the initial time then correct it now:
c --------------------------------------------------------------------------
      do j= 2,nbctim
          if(bc(j,nk-iangrot) .le. 0.0)
     .                bc(j,nk-iangrot) = bc(j-1,nk-iangrot)
      enddo




      return
c
c --- fatal errors, stop code
c
 5000 write  (ncrt, 5100)  kpsi,nconmax,maxlimpt,nxeqd,ncontr,nlimtr
      write  (nout, 5100)  kpsi,nconmax,maxlimpt,nxeqd,ncontr,nlimtr
 5100 format (/ ' ERROR detected in subroutine REQDSK'          /
     .          '   parameter settings inconsistent with eqdsk' /
     .          '   kpsi , nconmax, maxlimpt = ', 3(2x, i6)     /
     .          '   nxeqd, ncontr , nlimtr   = ', 3(2x, i6)     /
     .          '   ONETWO cannot continue')
      call giveupus(nrguess)
      call STOP ('subroutine REQDSK: eqdsk consistency problem', 33)
c
 6000 write  (ncrt, 6100)  nxeqd, nw, nyeqd, nh
      write  (nout, 6100)  nxeqd, nw, nyeqd, nh
 6100 format (' ERROR: eqdsk nxeq =', i5, '  code nw =', i5 /
     .        '        eqdsk nyeq =', i5, '  code nh =', i5 /
     .        '        We require nxeqd .le. nw and nyeqd .le. nh')
      call giveupus(nrguess)
      call STOP ('subroutine REQDSK: problem #4', 34)
c
      end

      subroutine rhoset (iscr)
c
c
c ----------------------------------------------------------------------
c  subroutine RHOSET determines rho grid to be used in transport,
c  as well as other quantitites listed below:
c  NOTE: vectors with max index npsi, are defined over the psival
c        grid, with the first value corresponding to the plasma edge
c        and the npsi'th value corresponding to the magnetic axis.
c        Corresponding vectors (obtained by spline interpolation from
c        the npsi vectors) have length nj, are defined over the rho
c        (which is called r in the transport code) grid with the first
c        value corresponding to the magnetic axis and the nj'th value
c        corresponding to the plasma edge.
c        (the rho values corresponding to the psi values in psival
c        are stored in vector rho. the rho values corresponding
c        to psi values in psir are stored in r).
c --- input
c  through argument list:
c  iscr                  Switch for writing file associated with unit
c                        nscr, see FLUXAV for explanation
c  through INCLUDE files:
c  INCLUDE file param:
c  kj                  dimension of vectors indexed by nj below
c
c  INCLUDE file mhdpar:
c  kpsi                dimension of vectors indexed by npsi below
c
c  INCLUDE file constnts:
c  twopi
c
c  INCLUDE file io:
c  nscr                 fortran unit number for file writes
c
c  INCLUDE file machin
c  flim            value of f(psi) at plasma boundary,tesla/m
c  btor            toroidal field (in vacuum at rmajor,tesla)
c  rmajor          major radius to place where btor is quoted,m
c
c  INCLUDE file numbrs:
c  nj                    size of transport grid (see vectors indexd by nj)
c
c  INCLUDE file flxav
c  npsi                  size of vectors indexed by npsi below
c
c  INCLUDE file ifs
c  rhod_psi_ifs(1..nj)       IFS model grid as function of psi
c
c  INCLUDE file neo2dp:
c  epsp(1..npsi)          inv aspect ratio  (on psival grid)
c  xhm2p(1...npsi)       see FLUXAV for reference
c  xi11p(1...npsi)                "
c  xi33p(1...npsi)                "
c  xipsp(1...npsi)                "
c
c  INCLUDE file psig
c  vprime(1..npsi)        dv/dpsi     (m**3/(volt-sec))
c  ratave(1...npsi)       <R0**2/R**2> (on psival grid)
c  ratavei(1...npsi)      <R**2>    (m**2)
c  fpsi(1...npsi)         f(psi)   (tesla/m)
c  psival(1...npsi)       psi grid vector (volt-se/rad)
c  rho(1...npsi)          rho grid vector (meters)
c  ravg(1...npsi)         <R>,meters
c  grho1npsi(1...npsi)        <ABS (grad psi)>
c  grho2npsi(1...npsi)        <(grad psi)**2>
c  elongx(1...npsi)        vertical elongation
c
c --- temporary storage used
c
c  INCLUDE file storage:
c  xdum(1...nj)               vector min length nj
c  ydum(1...nj)               vector min length nj
c  zdum(1...nj)               vector min length nj
c  udum(1...nj)               vector min length nj
c  vdum(1...nj)               vector min length nj
c  wdum(1...nj)               vector min length nj
c  sdum(1...nj)               vector min length nj
c  tdum(1...nj)               vector min length nj
c
c --- output (all through INCLUDE files):
c
c  INCLUDE file extra:
c  q(1...nj)                  safety factor (on r grid)
c  INCLUDE file flxav:
c  rhomax                 value of rho, meters, at plasma boundary
c                            (same as rhoa,remove one of them)
c
c  INCLUDE file geom
c  rhoa                r(nj) (i.e., rho at plasma boundary, meters)
c  r2cap(1...nj)       <R0**2/R**2>
c  r2capi(1...nj)      <R**2>   (m**2)
c  rcap(1...nj)        <R>   (m)
c  rcapi(1..nj)        <1/R> m-1
c  fcap(1...nj)        f(psilim)/f(psi)
c  hcap(1...nj)        fcap/<rmajor**2/R**2>
c  bsqncap(1...nj)      <(Bt0/B)**2>
c  b_avg_cap(1...nj)    <B>  ?? <B/Bt0>
c  h_factr_cap(1....nj)
c  INCLUDE file etc:
c  bp(1..nj)           poloidal b field (tesla)
c
c  INCLUDE file mesh
c  r(1..nj)            effective minor radius,meters,(the rho grid)
c  drr(1..nj)
c  rrp(1..nj)
c  rrm(1..nj)
c  dr(1...nj)
c  ra(1...nj)           these are various grid quantities used in transport
c
c  INCLUDE file neo2d:
c  eps(1..nj)          horizontal inverse aspect ratio
c  xhm2(1...nj)        see FLUXAV for reference
c  xi11(1...nj)              "
c  xi33(1...nj)              "
c  xips(1...nj)              "
c  elong_r(1...nj)   elongation on rho grid
c
c  INCLUDE file psig:
c  qpsi(1...npsi)       safety factor
c  rho(1...npsi)        rho vector(over psival grid,meters)
c  grho1npsi(1...npsi)      <ABS (grad rho)>
c  grho2npsi(1...npsi)      <(grad rhos)**2>
c
c
c  INCLUDE file ifs
c  rhod_ifs(1..nj)       IFS model grid as function of rho
c  grho1_ifs(1,..nj)      abs < grad rho > on rhod_ifs grid
c  grho2_ifs(1,...nj)    < (grad rho )**2> on rhod_ifs grid
c
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE mhdpar 
      USE io
      USE extra 
      USE numbrs 
      USE mesh
      USE machin
      USE geom
      USE constnts
      USE psig
      USE rhog
      USE ifs
      USE flxav
      USE metrics
      USE neo2d
      USE neo2dp
      USE etc
      USE cer
      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'imsl.i'

      include 'storage.i'

c
      dimension    torflux(kstore), cspline(kpsi,3)
      logical      qmonotonic, qfail
      equivalence (torflux(1)  , xdum(1))
      equivalence (cspline(1,1), vdum(1))
c
c ----------------------------------------------------------------------
c calculate safety factor (on psival grid)
c set switch if not monotonic
c ----------------------------------------------------------------------
c
      qmonotonic = .true.
      tensave    =  tensionspl
c
      do j=1,npsi
        qpsi(j) = fpsi(j) * vprime(j) * ratave(j)
     .                       / (twopi * rmajor)**2    ! eq. 2.4-3
      end do
c
      do j=2,nj
        if (q(j) - q(j-1) .ge. 0.0)  qmonotonic = .false.
      end do
c
c ----------------------------------------------------------------------
c calculate toroidal flux, torflux:
c use dvolume/dtoroidalflux = twopi / (fpsi*<R**-2>) (eq. 2.2-4)
c because vprime does not exist on the separatrix .... HSJ .... 11/24/93
c ----------------------------------------------------------------------
c
****  call trap2(psival,qpsi,torflux,npsi)
****  do 130 j=1,npsi
**130 torflux(j) = (torflux(j)-torflux(npsi))*twopi
c
      call calc_torflux (torflux, psivolp, fpsi, ratave, npsi, rmajor)
c
c ----------------------------------------------------------------------
c convert from torflux to effective radius
c if requested, use the area integration done in FLUXAV to define the
c toroidal flux of the outermost flux surface
c ----------------------------------------------------------------------
c
      if (itorfluse .eq. 1)  torflux(1) = torfluxnpsi(1)
      do j=1,npsi
        torfluxnpsi(j) = torflux(j)
        rho        (j) = SQRT (torflux(j)/(pi*btor))
      end do
      rhomax = rho(1)
c
c ----------------------------------------------------------------------
c --- flux surface avg grad rho
c     at this point we have
c                  <ABS (grad psi)> on the psi grid
c     we now convert this to
c                  <ABS (grad rho)> = ABS (q/(rho*bt0))*<ABS (grad psi)>
c     on the psi grid:
c ----------------------------------------------------------------------
c
      do j=1,npsi-1           ! j = npsi is magnetic axis
        grho1npsi(j) = ABS (grho1npsi(j) * qpsi(j) / (rho(j) * btor))
      end do
      grho1npsi(npsi) = 1.0
c
c ----------------------------------------------------------------------
c --- flux surface avg( grad rho)**2:
c     at this point we have
c                 < (grad psi)**2> on the psi grid.
c     we now convert this to
c                 <(grad rho)**2> = (q/(rho*bt0))**2*(<ABS (grad psi)>)
c     on the psi grid:
c ----------------------------------------------------------------------
c
      do j=1,npsi-1           ! j = npsi is magnetic axis
        grho2npsi(j) = grho2npsi(j)*(qpsi(j)/(rho(j)*btor))**2
      end do
      grho2npsi(npsi) = 1.0
c
c --- calculate the surface area of the flux surfaces:
c --- (this is a check on the sfareanpsi calcs done in FLUXAV):
c
      errmax = 0.0
      do j=1,npsi-1
        wdum(j) = 4.0 * (pi*rmajor)**2*btor*rho(j)*grho1npsi(j)
     .                / (fpsi(j)*ratave(j))
        errmax  = MAX (errmax, ABS (wdum(j)-sfareanpsi(j)))
      end do
c
c ----------------------------------------------------------------------
c calculate various quantities on the r mesh
c ----------------------------------------------------------------------
c

      do 160 j=1,nj
  160 r(j) = r(j)*rho(1)/r(nj)
c
      r(nj) = rho(1)                ! roundoff problems


c
c --- get q on r grid
c --- make sure monotonicity is preserved if necessary
c
      imslmd = 'rhoset'
      call intrp (0, 1, rho, qpsi, npsi, r, q, nj)
c      write(89,3123)(rho(j),qpsi(j),j=1,npsi)
c      write(89,3123)(r(j),q(j),j=1,nj)
 3123 format(2(2x,1pe14.6))
c
  100 qfail = .false.
      if (qmonotonic) then
        do j=2,nj-1
          if ((q(j+1)-q(j))*(q(j)-q(j-1)) .le. 0.0)  qfail = .true.
        end do
        if (qfail .and. tensave .eq. 0.0)  tensave = -1.0
      end if
      nvectr = kstore / nj
c
c     how many vectors size nj fit in storage vectors
c
      if (nvectr .ge. 6 .and.
     .   ((tensionspl .ne. 0.0) .or. qfail)) then
c
c         user selected tension spline (tensionspl .ne. 0)
c                          OR
c         q on psi was monotone, but q on rho
c         is not (due to poor interpolation);
c         fix the problem by forcing interpolation
c         to preserve monotonicity. (Non-monotone
c         q can be legitimate, so do this only
c         under the circumstances detected above.)
c
c         use tension spline to straighten any wiggle in q
c
          tensionspl = tensave
   90     tdum(1)    = -1.0e30
          tdum(2)    =  0.0
          tdum(3)    =  0.0
          tdum(4)    =  0.0    ! tdum sets b.c...
c                              ..zero gradient at r = 0, natural at r=a
c
          tension    = ABS (tensionspl)
          call tspline (rho,qpsi,npsi,tdum,cspline,kpsi,ier,tension,
     .           wdum(1),wdum(kpsi+1),wdum(2*kpsi+1),wdum(3*kpsi+1),
     .                  wdum(4*kpsi+1),wdum(5*kpsi+1),tenmax,r,q,nj)
          if (     tensionspl  .lt. 0.0     .and.
     .        ABS (tensionspl) .lt. tenmax) then
              do j=1,nj-1                      ! check if q is monotonic
                if (q(j+1) .lt. q(j))  go to 91
              end do
              go to 92
   91         tensionspl = 10.0 * tensionspl
              go to 90
          end if
      end if
c
   92 do 170 j=1,nj
  170 bp(j) = r(j)*btor / (rmajor*q(j))      ! eq. 2.4-4

c
      call trap2(r,bp,psir,nj)
      dpsirdif = rmajor * ABS (psir(nj))     ! PSIR(1) = 0.0 HERE
c
c --- get psir grid corresponding to r grid
c
      do j=1,nj
        psir(j) =  psival(npsi)
     .          + (psival(1) - psival(npsi)) * psir(j) / psir(nj)
      end do
c
      dpsirtrue = ABS (psir(nj) - psir(1))
c
c     modify q slightly to get correct delta psi
c
      do j=1,nj
        q (j) = dpsirdif * q(j) / dpsirtrue
        bp(j) = r(j) * btor / (rmajor * q(j))
      end do


c
      tensionspl = tensave              ! for next time RHOSET is called
c
      call intrp (0, 1, rho,    fpsi, npsi, r,    ydum, nj)
      call intrp (0, 1, rho,  ratave, npsi, r,    r2cap, nj)
      call intrp (0, 1, rho, ratavei, npsi, r,    r2capi, nj)
      call intrp (0, 1, rho,    ravg, npsi, r,    rcap, nj)
      call intrp (0, 1, rho,    ravgi, npsi, r,   rcapi, nj)
c
      call intrp (0, 1, rho, bsqinvavg, npsi, r,   bsqncap, nj)
      call intrp (0, 1, rho, bsq_avg, npsi, r,   bsq_avg_cap, nj)
      call intrp (0, 1, rho, b_avg, npsi, r,   b_avg_cap, nj)
      call intrp (0, 1, rho, h_factr, npsi, r, h_factr_cap, nj)

C      print *, 'h_factr_cap =',h_factr_cap
C      print *,'b_avg_cap =', b_avg_cap
C      call stop('209,line 2666',1)
c

      do 190 j=1,nj
      fcap(j) = flim/ydum(j)
  190 hcap(j) = fcap(j)/r2cap(j)
      rhoa = r(nj)
c
c ----------------------------------------------------------------------
c write rho for later use; add to data written by FLUXAV (iscr = 1)
c ----------------------------------------------------------------------
c

      if (iscr .eq. 1) then
        write  (nscr, 1010)  (100.0 * rho (j), j=1,npsi)
        write  (nscr, 1010)  (1.0e5 * psir(j), j=1,nj  )
 1010   format (6e12.5)
        mdum = 0
        write (nscr, '(i6)') mdum 
        call giveupus(nscr)
        close (unit = nscr)
      end if

c
c ----------------------------------------------------------------------
c update mesh quantities
c ----------------------------------------------------------------------
c
      do j=2,nj
        dr(j-1) = r(j)-r(j-1)
        ra(j-1) = 0.5 * (r(j-1)+r(j))
      end do
      drr(1) = 2.0 / dr(1)
      rrp(1) = 2.0 / dr(1)
      do j=2,nj
        rrm(j) = (ra(j-1)/r(j))/dr(j-1)
        if (j .ne. nj) then
          drr(j) = 2.0/(dr(j-1)+dr(j))
          rrp(j) = (ra(j)/r(j))/dr(j)
        end if
      end do
      drr(nj) = 2.0 / dr(nj-1)
c
c ----------------------------------------------------------------------
c quantities required for 2d neoclassical transport coeff.
c for example, given epsp, defined over the rho(npsi) grid (which
c is now in 1:1 correspondence with the psival(npsi) grid), determine
c eps, defined over the r(nj) grid,by spline interpolation with
c zero gradient condition at the magnetic axis and free (natural)
c boundary condition at the plasma edge:
c ----------------------------------------------------------------------
c
      call intrp (0, 1, rho, epsp  , npsi, r, eps    , nj)
      call intrp (0, 1, rho, xhm2p , npsi, r, xhm2   , nj)
      call intrp (0, 1, rho, xi11p , npsi, r, xi11   , nj)
      call intrp (0, 1, rho, xi33p , npsi, r, xi33   , nj)
      call intrp (0, 1, rho, xipsp , npsi, r, xips   , nj)
      call intrp (0, 1, rho, elongx, npsi, r, elong_r, nj)
      call intrp (0, 1, rho, ftnclp, npsi, r, ftncl  , nj)
      call intrp (0, 1, rho, grho1npsi, npsi, r, grho1_mesh, nj)
      call intrp (0, 1, rho, grho2npsi, npsi, r, grho2_mesh, nj)
      call intrp (0, 1, rho, rhod_psi_ifs, npsi, r,rhod_ifs, nj)
c
c     rhod_ifs is horizontal minor radius on the r grid.
c     hence grho1_mesh and grho2_mesh,both of which are defined
c     on the r grid,are also automatically valid on the rhod_ifs grid
c
      call copya (grho1_mesh, grho1_ifs, nj)
      call copya (grho2_mesh, grho2_ifs, nj)

c
c
****  call intrp (0, 1, r, grho1_mesh, nj, rhod_ifs, grho1_ifs, nj)
****  call intrp (0, 1, r, grho2_mesh, nj, rhod_ifs, grho2_ifs, nj)
c
      rhod_max_ifs = rhod_ifs(nj)*100.0 ! save norm. factor in cm
      call multpl1 (rhod_ifs,nj,1.0/rhod_ifs(nj)) ! normalize rhod_ifs
c
c ----------------------------------------------------------------------
c Quantities required for Houlberg's Bootstrap Model (jhirsh = 95 or 96)
c Transform from psival/rho grid to the r grid.
c The order of the arrays will be reversed,
c since psival(npsi) is center but r(nj) is edge.
c ----------------------------------------------------------------------
c
      if (jhirsh0 .eq. 95 .or. jhirsh0 .eq. 96) then
c
c       interpolate on to r grid (will be reversed!!)
c
        do j=1, npsi
          ydum(j) = bsq(j)
        end do
c
        call intrp (0, 1, rho, ydum, npsi, r, bsq, nj)
c
c       Save central value of bmsq (at npsi) to compute grbmsq at center
c
        bmsq0     = bmsq(npsi)
        do j=1,npsi
          ydum(j) = bmsq(j)
        end do
c
        call intrp (0, 1, rho, ydum, npsi, r, bmsq, nj)
        do j=1,npsi
          ydum(j) = grth(j)
        end do
        call intrp (0, 1, rho, ydum, npsi, r, grth, nj)
        do   m=1,3
          do j=1,npsi
            ydum(j) = gfm(m,j)
          end do
          call intrp (0, 1, rho, ydum, npsi, r, zdum, nj)
          do j=1,nj
            gfm(m,j) = zdum(j)
          end do
        end do
c
c       flux surface avg(grad rho)**2/Btot**2):
c       at this point we have:
c          <(grad psi)**2/Btot**2>
c       on the psi grid. we now convert this to:
c          <(grad rho)**2/Btot**2> = (q/(rho*bt0))**2*(<ABS (grad psi)/Btot**2>)
c       on the psi grid
c
c       Then transform onto rho grid (will be reversed!!)
c
        do j=1,npsi-1           ! j = npsi is magnetic axis
          ydum(j)  = grbmsq(j)*(qpsi(j)/(rho(j)*btor))**2
        end do
        ydum(npsi) = bmsq0
        call intrp (0, 1, rho, ydum, npsi, r, grbmsq, nj)
c
        if (j_cer .ne. 0) then
c
c         transform Bt/R along CER cord at Z=0
c
          do j=1,npsi
            ydum(j) = cer_btdr(j)
          end do
          call intrp (1, 1, rho, ydum, npsi, r, cer_btdr, nj)
c
c         transform Bz along CER cord at Z=0
c
          do j=1,npsi
            ydum(j) = cer_bp(j)
          end do
          call intrp (0, 1, rho, ydum, npsi, r, cer_bp, nj)
        end if
      end if
c
      call dump_values ('from RHOSET')
      return
c
      end

      subroutine savcur

c
c ----------------------------------------------------------------------
c     This subroutine saves the current profile and related quantities.
c ----------------------------------------------------------------------
c
c
      USE param
      USE soln
      USE numbrs
      USE mesh
      USE rhog
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'rhog.i'
c      include 'soln.i'
c
      do j=1,nj
        rold  (j) = r     (j)
        psiold(j) = psir  (j)
        ppold (j) = pprim (j)
        ffpold(j) = ffprim(j)
        curold(j) = curden(j)
      end do
      return
c
      end

      subroutine savit
c

c
c ----------------------------------------------------------------------
c saves solution info on disk in case later an iteration is required;
c information is read back by subroutine READIT
c ----------------------------------------------------------------------
c
      USE param
      USE fusion
      USE io 
      USE solcon
      USE soln
      USE nub2
      USE mhdpar  
      USE mesh
      USE mixcom
      USE rhog
      USE mhdcom
      implicit  integer (i-n), real*8 (a-h, o-z)

      
c
      rewind (unit = nsavsol)
      write (nsavsol) u
      write (nsavsol) r
      write (nsavsol) time0,time,dt,n,eqtim0
      write (nsavsol) enbsav
      write (nsavsol) wbsav
      write (nsavsol) enasav
      write (nsavsol) wasav
      write (nsavsol) psi
      write (nsavsol) psir
      write (nsavsol) psibdry,psiaxis,rma,zma
c
c --- subroutine mix parameters
c
      write (nsavsol)  dtecal,dtemix,dtical,dtimix,epste,epsti
      write (nsavsol)  fuscal,fusmix,ipmix,mixpro,qmix,rmixx
      write (nsavsol)  rsmixx,rsx,s3cal,s3mix,s71cal,s71mix
      write (nsavsol)  s18cal,s18mix
      write (nsavsol)  tem
      write (nsavsol)  tep
      write (nsavsol)  tim
      write (nsavsol)  tip
      write (nsavsol)  timmix,trcal,trmix,tsmix
      write (nsavsol)  tdmix
      write (nsavsol)  wmix,w0mix,w1mix,w2mix,w3mix,w4mix,w5mix
      return
c
****  write (nitre, 1) psibdry,psiaxis,rma,zma,psi(1,1),psi(17,33),
**** .                 psir(1),psir(51)
****1 format (' subroutine SAVIT, psibdry,psiaxis,rma,zma' /
**** .          2x,8(2x,1pe14.7))
c
      end

      subroutine set_bdry (rcmin, rcmax, zcmin, zcmax, zrcmin,
     .                     zrcmax, rzcmin, rzcmax, bpcontr,
     .                     rmhdgrid, zmhdgrid, nw, nh, cspln,
     .                     n2cspln, nh2, iounit)
c ------------------------------------------------------------------ HSJ
c set the boundary and asssociated stuff from the values on eqdsk
c ----------------------------------------------------------------------
c
      USE param
      USE contour
      USE io
      USE replace_imsl,         ONLY : my_dbcevl1
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c     include 'io.i'
c      include 'contour.i'
c
      dimension bpcontr(*),cspln(n2cspln,nw,nh2),
     .          rmhdgrid(nw),zmhdgrid(nh),pds(6)
c
      ier    = 0
      call copya (rplasbdry,rcontr,nplasbdry)
      call copya (zplasbdry,zcontr,nplasbdry)
      ncontr =   nplasbdry
      rcmax  = rplasbdry(1)
      rcmin  = rplasbdry(1)
      zcmax  = zplasbdry(1)
      zcmin  = zplasbdry(1)
      zrcmin = zcmax
      zrcmax = zcmin
      rzcmin = rcmin
      rzcmax = rcmax
      do j=1,nplasbdry
        rcmax = MAX (rcmax,rplasbdry(j))
        rcmin = MIN (rcmin,rplasbdry(j))
        zcmax = MAX (zcmax,zplasbdry(j))
        zcmin = MIN (zcmin,zplasbdry(j))
        if (rcmin .eq. rplasbdry(j))  zrcmin = zplasbdry(j)
        if (rcmax .eq. rplasbdry(j))  zrcmax = zplasbdry(j)
        if (zcmin .eq. zplasbdry(j))  rzcmin = rplasbdry(j)
        if (zcmax .eq. zplasbdry(j))  rzcmax = rplasbdry(j)
        call my_dbcevl1 (rmhdgrid, nw, zmhdgrid, nh, cspln,
     .                nw, rplasbdry(j),
     .                zplasbdry(j), pds, ier, 3)
        if (ier .ne. 0) then
          if (iounit .ne. 0)
     .    write  (iounit, 4) ier
          write  (nout, 4) ier
    4     format (/ ' subroutine SET_BDRY detects an error:' /
     .              ' DBCEVL1 returned ier =', i5 /
     .              ' while getting psi on axis')
          call STOP ('subroutine SET_BDRY: problem #1', 219)
        end if
        if (rplasbdry(j) .gt. 0) then
          bpcontr(j) = SQRT (pds(2)**2 + pds(3)**2) / rplasbdry(j)
        else
          bpcontr(j) = 0
        end if
      end do
      return
c
      end

      subroutine sgtsl (n, c, d, e, b, info)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c   tridiagonal solver with partial pivoting
c --- input
c  n          number of equations
c  d(i)       i = 1,2...n  diagonal          destroyed on output
c  c(i)                    lower diagonal    destroyed on output
c  e(i)                    upper diagonal    destroyed on output
c  b(i)                    rhs               contains solution on output
c --- output
c  b(i)       i = 1,2,...n solution vector
c  info       error flag, = 0  no error
c                         = k  if k'th pivot is zero exactly
c                              in this case no solution is possible
c ------------------------------------------------------------------ HSJ
c
      dimension  b(*), c(*), d(*), e(*)
c
      info = 0
      c(1) = d(1)
      nm1  = n - 1
      if (nm1 .lt. 1)  go to 40
      d(1) = e(1)
      e(1) = 0.0
      e(n) = 0.0
c
      do k=1,nm1
        kp1 = k + 1
c
c --- find largest of the two rows
c
        if (ABS (c(kp1)) .lt. ABS (c(k)))  go to 10
c
c --- interchange row
c
        t      = c(kp1)
        c(kp1) = c(k)
        c(k)   = t
        t      = d(kp1)
        d(kp1) = d(k)
        d(k)   = t
        t      = e(kp1)
        e(kp1) = e(k)
        e(k)   = t
        t      = b(kp1)
        b(kp1) = b(k)
        b(k)   = t
c
c --- zero elements
c
   10   if (c(k) .ne. 0.0e0)  go to 20
        info   = k
        go to 70
c
   20   t      = -c(kp1)/c(k)
        c(kp1) =  d(kp1)+t*d(k)
        d(kp1) =  e(kp1)+t*e(k)
        e(kp1) =  0.0e0
        b(kp1) =  b(kp1)+t*b(k)
      end do
c
   40 if (c(n) .ne. 0.0e0)  go to 50
      info = n
      go to 70
c
c --- back solve
c
   50 nm2  = n-2
      b(n) = b(n)/c(n)
      if (  n .eq. 1)  go to 70
      b(nm1) = (b(nm1)-d(nm1)*b(n))/c(nm1)
      if (nm2 .lt. 1)  go to 70
      do kb=1,nm2
        k    = nm2 - kb + 1
        b(k) = (b(k)-d(k)*b(k+1)-e(k)*b(k+2))/c(k)
      end do
c
   70 return
c
      end

      subroutine solvegs (iounit, ierr, kstart, kend, ifixbdry,
     .                    iavfcoil, nout, ncrt)
c
c
c ----------------------------------------------------------------------
c --- subroutine solves the grad shafranov equation as described below.
c --- The equation is treated here under the assumption of a
c --- right hand triad (r,z,phi),with psi defined as the negative of the
c --- poloidal flux between the symmetry (z) axis and the psi surface.
c --- thus we have
c                 del-star (psi) = u0*rmajor*jtoroidal
c --- with jtoroidal positive in the positive phi direction.
c --- and jtoroidal = -rmajor*dp/dpsi-(f/(u0*rmajor))*df/dpsi
c --- with this convention psi will be minimum at the magnetic
c --- axis and increase as the plasma boundary is approached,
c --- if the toroidal current is positive.
c
c --- input (through argument list):
c  iounit       fortran unit number for error messages. set iounit
c               =0 to suppress all messages from this subroutine.
c  kstart
c  kend             use these to save going over the whole grid
c  ifixbdry         =0 load psi1d on boundary
c                   =1 load not necessary (mhdmode ='no coils')
c                          or boundary value due to plasma current
c                          is fixed (mhdmode ='coils')
c  iavfcoil           =0 don't average f coil currents
c                     =1 average,see subroutine GETFCUR
c --- input (through INCLUDE files):
c --- INCLUDE file mhdpar,mhdcom and mhdbcdtn:
c  mhdmode      if mhdmode = 'no coils' it is assumed that the only
c               toroidal current present on the MHD grid is due to
c               the plasma. any other toroidal currents are outside
c               the MHD grid and their influence is manifest only
c               through the boundary conditions on psi.
c               if mhdmode ='coils' then the coils are modeled explicitly
c               as described below.
c  mhdmethd       if mhdmethd = 'green' then the Green's function is used to
c               invert del-star. if mhdmethd ='cycred' then the grad-
c               shafranov equation is solved using (single) complete
c               cyclic reduction. in the future we may change this to
c               partial cyclic reduction in combination with fast fourier
c               transforms (the facr(l) method). double cyclic reduction
c               is not used since the reduced set of nw-2 tridiagonal
c               equations is solved twice as economically using a standard
c               tri-diagonal solver.
c               note that if mhdmode = 'no coils' then mhdmethd='cycred'
c               must be used (it will be set below if set incorrectly).
c
c  pcurrent(j)   j = 1,...nwh  is the plasma current ,in amps,associated
c                with a grid filament of area dr*dz
c  gridpc(i,j)   i = 1,...nwh,j=1...nw. gives the inductive coupling
c                of the plasma current to the grid points. note:
c                the peculiar structure of gridpc is due to the fact
c                that the coupling depends on r,r0 and (z-z0)**2 only,
c                where (r0,z0) is the source point and (r,z) the field
c                point. using this symmetry gridpc is reduced from
c                an (nwh,nwh) matrix to just an (nwh,nw) matrix.)
c
c       if mhdmode = 'coils' the following coil specs are required
c
c  gridfc(i,j)  i = 1,..nwh,j=1,nfcoil is the inductive coupling of
c               f-coil j to grid point i.
c  curfcoil(j)  j = 1,...nfcoil  the fcoil currents
c  nfcoil       number of f coils
c  fixfcoil         integer,sets f coil current option (see init)
c
c  ivessel      switch,if ivessel=0 do not include vessel currents.
c               if ivessel = 1,include vessel currents.
c               if ivessel = 1 then the following vessel related quantities
c               are required:
c  gridvs(i,j)  i = 1...nwh,j=1,..nvessel. the inductive coupling of vessel
c               "filament" j to grid point i.
c  vescur(j)    j = 1,...nvessel  gives the vessel currents.
c  nvessel      number of "filaments" used to model the vessel
c
c  iecurr       switch,iecurr=0 means do not include ecoil currents
c               iecurr = 1 means include e coil currents.
c               if iecurr = 1 then the following e coil related quantities
c               must be input:
c  gridec(i,j)  i = 1,....nwh,j=1,nesum. the inductive coupling of
c               e coil bank j to grid point i.
c  ecoilcur(j)    the current flowing in e coil bank j
c  nesum        number of e coil banks.
c
c  INCLUDE file mhdgrid:
c  rmhdgrid(i)      i = 1,...nw    r coords of mhdgrid
c  zmhdgrid(j)      j = 1,...nh    z    "    "    "
c
c  INCLUDE file mhdpar:
c  nw,nh
c  nwh             give the size of the MHD grid,(rmhdgrid(nw),
c                   zmhdgrid(nh),nwh = nw*nh)
c                   note that the cyclic reduction scheme requires
c                   that nh = 2**m+1 for some positive integer m .
c                   the single reduction method used here does not
c                   place  any restriction on the radial grid so nw is
c                   arbitrary (it must of course be consistent with the
c                   values used to generateratave the Green's table).
c  INCLUDE file bicube:
c  wnoperm(j)          j = 1,2.....9(nw-2)  is a vector required if
c                   mhdmethd = 'cycred',unused otherwise. most of this
c                   storage vector can be eliminated by making
c                   appropriate changes below. It has not been done
c                   here to make the routine more easily understood.
c                   Also we will modify this routine to use the facr(l)
c                   method in the future,which will require additional
c                   storage anyway.
c  NOTE:             nw,nh,nwh,and wnoperm must be passed from INCLUDE files
c                    in order that the equivalencing below works.
c
c --- INCLUDE file mhdcom:
c  psi(i,j)              used as a temporary work array here only.
c
c  isym                 if isym = 1 then the upper half of the solution is
c                       copied into the lower half,yielding an up/down
c                       symmetric solution.
c
c --- output
c  INCLUDE file mhdcom:
c  psi1d(kk)      kk = 1,...nwh    1d version of psi (kk = (i-1)*nh+j)
c               is the computed solution. Note that this solution
c               is returned by columns (which is not the fortran default
c               convention,but is done to retain compatibility with
c               the Green's function arrays generated by older codes.
c               the fortran standard would require kk = (j-1)*nw+i) in
c               order for psi1d(kk) and psi(i,j) to be the same thing)
c  INCLUDE file mhdbcdtn:
c  curfcoil(i)     i = 1....nfcoil,the fcoil currents if ifixfcoil=0
c                  otherwise curfcoil remains unchanged on output.
c  through argument list:
c  ierr          error flag. if ierr = 0 no error ocurred. if ierr > 0
c                some error occurred and solution can not be used.
c
c ------------------------------------------------------------------ HSJ
c
      USE mhdpar 
      USE mhdgrid   
      USE constnts
      USE mhdcom
      USE bicube
      USE mhdbcdtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'mhdpar.i'
c      include 'bicube.i'
c      include 'constnts.i'
c      include 'mhdbcdtn.i'
c      include 'mhdcom.i'
c      include 'mhdgrid.i'
c
      dimension    diag(nw-2),diagl(nw-2),diagu(nw-2),diag1(nw-2),
     .             phi(nw-2),phi1(nw-2),v(nw-2),wk(nw-2),wk2(nw-2)
c      equivalence (diag(1),wnoperm(1)),(diagl(1),wnoperm(nw-1)),
c     .            (diagu(1),wnoperm(2*nw-3)),(diag1(1),wnoperm(3*nw-5)),
c     .            (phi(1),wnoperm(4*nw-7)),(phi1(1),wnoperm(5*nw-9)),
c     .            (v(1),wnoperm(6*nw-11)),(wk(1),wnoperm(7*nw-13)),
c     .            (wk2(1),wnoperm(8*nw-15))
c
      darea = (rmhdgrid(2)-rmhdgrid(1))*(zmhdgrid(2)-zmhdgrid(1))
c
c --- get the boundary condition for psi1d on the border of the MHD grid.
c --- this is required only if mhdmode= 'no coils' or mhdmode = 'coils'
c --- and mhdmethd ='cycred'.
c --- (if mhdmethd ='green' then the boundary condition
c --- is psi1d(r = 0,z) = 0 and psi1d(r=infinity,z=infinity) = 0,
c --- which is implicit in the Green's function).
c
      ierr = 0
      if (ifixbdry .eq. 1) then
c
c --- copy psincbcd into the boundary of the psi1d grid.
c --- (psincbcd was loaded on a previous call to this routine
c --- see below).
c --- left side:
c
        do 10 j=1,nh
   10   psi1d(j) = psincbcd(j)
c
c --- across the top:
c
        k0 = nh-1
        do 15 j=2,nw
           k = j*nh
   15   psi1d(k) = psincbcd(k0+j)
c
c --- down the right side
c
        k0 = nh+nw-1
        do 20 j=1,nh-1
          k = nwh-j
   20   psi1d(k) = psincbcd(k0+j)
c
c --- across the bottom from right to left
c
        k0 = 2*nh+nw-2
        do 25 j=1,nw-2
          k = nwh-(j+1)*nh+1
   25   psi1d(k) = psincbcd(k0+j)
        mhdmethd = 'cycred'        ! only method available for this case
      else if (mhdmethd .eq. 'cycred' .and. ifixbdry .eq. 0) then
c
c --- get psi1d on boundary of MHD grid due to plasma current only:
c --- vertical boundaries ((rmhdgrid(1),zmhdgrid(j)),j = 1,..nh,
c --- and (rmhdgird(nw),zmhdgrid(j)),j = 1,...nh):
c
        do 50 j=1,nh
          kk = (nw-1)*nh+j
          psi1d(j) = 0
          psi1d(kk) = 0
          do 60 ii=1,nw
            do 60 jj=1,nh
              kkkk = (ii-1)*nh+jj
              mj = IABS (j-jj)+1
              mk = (nw-1)*nh+mj
              psi1d(j) = psi1d(j)+isgngren*gridpc(mj,ii)*pcurrent(kkkk)
              psi1d(kk) = psi1d(kk)
     .                  + isgngren * gridpc(mk,ii) * pcurrent(kkkk)
   60     continue
   50   continue
c
c --- horizontal boundaries ((rmhdgrid(i),zmhdgrid(1)),i = 2,..nw-1,
c --- and (rmhdgird(i),zmhdgrid(nh)),i = 2,...nw-1):
c
        do 70 i=2,nw-1
          kk1 = (i-1)*nh
          kknh = kk1+nh
          kk1 = kk1+1
          psi1d(kk1) = 0.0
          psi1d(kknh) = 0.0
          do 80 ii=1,nw
            do 80 jj=1,nh
            kkkk = (ii-1)*nh+jj
            mj1 = IABS (jj-1)+1
            mjnh = IABS (nh-jj)+1
            mk1 = (i-1)*nh+mj1
            mknh = (i-1)*nh+mjnh
            psi1d(kk1 ) =
     .      psi1d(kk1 ) + isgngren * gridpc(mk1 ,ii) * pcurrent(kkkk)
            psi1d(kknh) =
     .      psi1d(kknh) + isgngren * gridpc(mknh,ii) * pcurrent(kkkk)
   80     continue
   70   continue
          if (ivertsbl .eq. 1) then   ! save boundary values in pisncbcd
             do 500 j=1,nh-1          ! they will be used if stability
  500        psincbcd(j) = psi1d(j)   ! logic kicks in
             k = nh-1
             do 510 i=1,nw
               kk = i*nh
               k = k+1
  510          psincbcd(k) = psi1d(kk)
              do 520 j=nh-1,1,-1
               kk = (nw-1)*nh+j
               k = k+1
  520          psincbcd(k) = psi1d(kk)
             do 530 i=nw-1,2,-1
               kk = (i-1)*nh+1
               k = k+1
  530          psincbcd(k) = psi1d(kk)
          end if
      end if
c
c --- process the fcoil currents.
c --- These are known in curfcoil if fixfcoil = 1, else we now calculate them
c
      if (fixfcoil .eq. 0)  call getfcur(ncrt,ierr,kstart,kend,iavfcoil)
      if (ierr .ne. 0) then
          write (nout  , 531) ierr
          write (ncrt  , 531) ierr
          if (iounit .ne. 0)
     .    write (iounit, 531) ierr
  531     format (' subroutine GETFCUR detected an error' /
     .            ' program must stop')
          call STOP ('subroutine SOLVEGS: unspecified problem', 100)
      end if
c
      if (mhdmethd .eq. 'green') then
c
c --- Green's function solution:
c --- note: the accumulation of psi1d below can cause difficulty on low
c --- precision machines (i.e. VAX) use double precision and/or
c --- order the terms in magnitude to avoid roundoff for such machines
c
      do 100 j=1,nh
        do 100 i=1,nw
          kk = (i-1)*nh+j
          psi1d(kk) = 0.0
c
c --- contribution from f coil:
c
          do 110 m=1,nfcoil
  110       psi1d(kk) = psi1d(kk)+isgngren*gridfc(kk,m)*curfcoil(m)
c
c --- contribution from vessel:
c
          if (ivessel .eq. 1) then
            do 120 m=1,nvessel
  120          psi1d(kk) = psi1d(kk)+isgngren*gridvs(kk,m)*vescur(m)
          end if
c
c --- contribution from e coil:
c
          if (iecurr .eq. 1) then
            do 130 m=1,nesum
  130         psi1d(kk) = psi1d(kk)+isgngren*gridec(kk,m)*ecoilcur(m)
          end if
c
c --- contribution from plasma current
c
          do 140 k0=kstart,kend
              ii = k0/nh+1
              jj = k0-(ii-1)*nh
              if (jj .eq. 0) then
                 ii = ii-1
                 jj = 65
              end if
              kkk = (ii-1)*nh+jj
              mj = IABS (j-jj)+1
              mk = (i-1)*nh+mj
  140         psi1d(kk) = psi1d(kk)+isgngren*gridpc(mk,ii)*pcurrent(kkk)
  100   continue
        if (isym .eq. 1)  call symmetrize(psi1d,nw,nh)
        return
      end if        ! end solution by Green's function
c
      if (mhdmethd .eq. 'cycred') then
c
c --- cyclic reduction solution. We solve del-star psi1d = rhs
c --- considering only the plasma current,with boundary conditions
c --- due only to the plasma current as determined above. then,if
c --- mhdmode = 'coils' we add to this solution the contribution from the coils.
c --- first set up the finite difference grid coefficients:
c --- (we don't save the results so have to recalc each time)
c --- (may want to incorparte this directly into tri diagonal solver)
c
        dr     = rmhdgrid(2) - rmhdgrid(1)
        dz     = zmhdgrid(2) - zmhdgrid(1)    ! constant grid spacing
        dzsq   = dz*dz
        dzdrsq = (dz/dr)**2
        dumy   = dzsq/(2.0*dr)
        dumy1  = 2.0 * (1.0 + dzdrsq)
        do i=1,nw-2
          diag(i)  =  dumy1
          denom    =  dumy/rmhdgrid(i+1)
          diagl(i) = -dzdrsq-denom
          diagu(i) = -dzdrsq+denom
        end do
c
c --- next set up the rhs term. Note that (-dzsq) was factored out
c --- of the finite difference expressions above and thus is included
c --- in the rhs term here. Division by darea is necessary because we
c --- have defined pcurrent as the current (in amps) of a filament of
c --- area darea. Here we require the actuall current density.
c --- the rhs is stored in the psi1d array because subroutine CMPLTCYR
c --- expects to find it there.
c
        dumy = dzsq*u0
        do 210 i=3,nw-2
          dumy1 = dumy*rmhdgrid(i)/darea
          do 210 j=2,nh-1
            kk        = (i-1)*nh+j
            psi1d(kk) = isgngren*dumy1*pcurrent(kk)
  210   continue
c
c --- finally add the finite difference expressions for i = 2,and
c --- i = nw-1 to the rhs. these terms are known boundary terms and
c --- hence do not form part of the tridiagonal matrix constructed above:
c --- for i = 2:
c
        dumy  = dzdrsq+dzsq/(2.0*rmhdgrid(2)*dr)
        dumy1 = dzsq*u0*rmhdgrid(2)/darea
        do j=2,nh-1
          kk        = nh+j
          psi1d(kk) = isgngren*dumy1*pcurrent(kk)+psi1d(j)*dumy
        end do
c
c --- for i = nw-1:
c
        dumy  = dzdrsq-dzsq/(2.0*rmhdgrid(nw-1)*dr)
        dumy1 = dzsq*u0*rmhdgrid(nw-1)/darea
        k  = (nw-2)*nh
        kk = k+nh
        do 230 j=2,nh-1
          kkk        = k+j
          kkkk       = kk+j
          psi1d(kkk) = isgngren*dumy1*pcurrent(kkk)+psi1d(kkkk)*dumy
  230   continue
c
c --- some additional input required in cmpltcyr:
c
        nw1 = nw
        i   = 1
        do 240 j=1,15
          i = i*2
          if (i .eq. nh-1) then
            nhpwr = j
            go to 250
          end if
  240   continue
        if (iounit .ne. 0)  write (iounit, 260) nh
        write  (nout, 260) nh
        write  (ncrt, 260) nh
  260   format (' subroutine SOLVEGS detects error,nh =', i5 /
     .          ' nh must be 2**m+1 for some integer m')
        ierr = 1      ! wrong nh
        return
c
  250   do j=1,nh
          do i=1,nw
            kk       = (i-1)*nh+j
            psi(i,j) = psi1d(kk)
          end do
        end do
c
c --- now call the solution routine.
c --- subroutine cmpltcyr requires that the rhs is
c --- initially stored in psi(i,j),i = 2,...nw-1,j=2,..nh-1.
c --- the boundary values are stored in psi(1,j),psi(nw,j),j = 1,..nh
c --- and psi(i,1),psi(i,nh),i = 2,....nw-1.  on output the boundary
c --- psi remains unchanged and the interior psi values are the
c --- solution found with the given boundary condition.
c
        call cmpltcyr (diagl,diag,diagu,psi,nw,nw1,nh,nhpwr,phi,
     .                 phi1,diag1,v,wk,wk2)
c
c --- temporary, check solution:
c
c     errmax = 0.0
c     errsum = 0.0
c     dzsq = (zmhdgrid(2)-zmhdgrid(1))**2
c     dr = rmhdgrid(2)-rmhdgrid(1)
c     dzdrsq = dzsq/(dr*dr)
c     kct = 0
c     do i=3,nw-2
c     t1 = -dzdrsq-dzsq/(2.0*dr*rmhdgrid(i))
c     t2 = 2.0 * (dzdrsq+1.0)
c     t3 = -dzdrsq+dzsq/(2.0*dr*rmhdgrid(i))
c     do j=2,nh-1
c     kct = kct+1
c     aleft = psi(i-1,j)*t1+psi(i,j)*t2+psi(i+1,j)*t3-psi(i,j-1)
**** .     -psi(i,j+1)
c     k = (i-1)*nh+j
c     err=ABS (aleft-psi1d(k))
c     errmax = MAX (err,errmax)
****  if (err .eq. errmax) then
c       ierror = i
c       jerror = j
****  end if
c     errsum = errsum+err
c     end do
c     end do
c     erravg = errsum/kct
c     write (6,26000)errmax,ierror,jerror,erravg
c     write (8,26000)errmax,ierror,jerror,erravg
*6000 format ('  errmax,i,j =',1pe12.3,2x,i4,i4,
**** .        '  erravg = ',1pe12.6)
c
c --- end of temporary solution check
c
        do j=1,nh
          do i=1,nw
            kk = (i-1)*nh+j
            psi1d(kk) = psi(i,j)
          end do
        end do
c
        if (mhdmode .eq. 'coils') then
c
c --- finally add psi1d contributions from current sources other than plasma.
c
          do 280 kk=1,nwh
            if (ivessel .ne. 0) then
              do 290 m=1,nvessel
  290           psi1d(kk) = psi1d(kk)+isgngren*gridvs(kk,m)*vescur(m)
            end if
            if (iecurr .ne. 0) then
              do 300 m=1,nesum
  300           psi1d(kk) = psi1d(kk)+isgngren*gridec(kk,m)*ecoilcur(m)
            end if
            do 310 m=1,nfcoil
  310         psi1d(kk) = psi1d(kk)+isgngren*gridfc(kk,m)*curfcoil(m)
  280     continue
        end if
        if (isym .eq. 1)  call symmetrize(psi1d,nw,nh)
        return
       end if
c
c --- error:
c
      if (iounit .ne. 0)
     .write  (iounit, 90)  mhdmode, mhdmethd
      write  (nout  , 90)  mhdmode, mhdmethd
      write  (ncrt  , 90)  mhdmode, mhdmethd
   90 format (' subroputine SOLVEGS reports ERROR: mhdmode, mhdmeth = ',
     .          a, 1x, a /
     .        ' program must stop')
      ierr = 1
      return
c
      end

      subroutine stark (time, nout, ncrt, rmajor, btor, fcap,
     .                  rmhdgrid, zmhdgrid, psir, nj)
c
c
c ----------------------------------------------------------------------
c --- subroutine gets exptl mse data by interpolation in time of input
c --- data and calculates corresponding theoretical values.
c
c --- input:
c  time         current time
c  rmajor       major radius at which btor is defined
c  btor         toroidal field in vac. at rmajor
c  fcap(j)          j = 1,2,..nj
c  rmhdgrid(i)      i = 1,2..nw
c  zmhdgrid(j)      j = 1,2..nh
c  psir(j)          j = 1,2,..nj  psi on rho grid (corresponding to fcap)
c
c  from include files:
c  file mhdpar.i
c          nh
c          nwh
c          nstark
c
c  file bicube.i
c          cspln
c          wnoperm
c          pds
c
c  file mhdbcdtn.i
c    in these arrays j indicates channel no. (up to nstark of them)
c    and i indicates time index (tied to timestark(i))
c    timestark(i)        i = 1,2..istark
c    tstark(j,i)         j = 1,2...?
c    sigstark(j,i)
c    fwtstark(j,i)
c    rstark(j,i)
c    zstark(j,i)
c    a1stark(j,i)
c    a2stark(j,i)
c    a3stark(j,i)
c    a4stark(j,i)
c
c  file mhdcom.i:
c  psi(i,j)         i = 1,2..nw,j=1,2..nh
c
c --- output:
c  file mhdbcdtn.i:
c     tstarkexptl(j)        j = 1,nstark current values (at t=time)
c     sigstarkexptl(j)
c     rstarkexptl(j)   this value is returned as zero for invalid pts.
c     zstarkexptl(j)
c     tstarkcalc(j)
c     chisqstark(j)
c     chisqstarktot
c     chisqstark(j)
c
c ------------------------------------------------------------------ HSJ
c
      USE mhdpar 
      USE mhdcom
      USE bicube
      USE mhdbcdtn
      USE replace_imsl,          ONLY : my_ibcccu, my_dbcevl1
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'mhdpar.i'
c      include 'bicube.i'
c      include 'mhdcom.i'
c      include 'mhdbcdtn.i'
      include 'storage.i'
c
      integer  j, jj, jjm1, jlo, nj, nout, ncrt, ier
c
c      real*8   time,timestark,dtime,dt,dval,rmhdgrid(nw),zmhdgrid(nh),
c     .         psir(nj),fcap(nj) , fwtstark,fwtstarkexptl(nstark),
c     .         tstark,tstarkexptl, sigstark,sigstarkexptl,
c     .         rstark,rstarkexptl, zstark,zstarkexptl,
c     .         a1stark,a1(nstark), a2stark,a2(nstark),
c     .         a3stark,a3(nstark), a4stark,a4(nstark)
      real*8   time,dtime,dt,dval,rmhdgrid(nw),zmhdgrid(nh),
     .         psir(nj),fcap(nj) , fwtstarkexptl(nstark),
     .         a1(nstark), a2(nstark),a3(nstark),a4(nstark)
c
      equivalence (xdum(1)         , a1(1))
      equivalence (xdum(nstark+1)  , a2(1))
      equivalence (xdum(2*nstark+1), a3(1))
      equivalence (xdum(3*nstark+1), a4(1))
      equivalence (xdum(4*nstark+1), fwtstarkexptl(1))
c
      do j=1,istark
        jj = j
        if (timestark(j) .ge. time )  go to 10
      end do
c
          if (istark .le. 0) then
            write (nout, 90)  istark, (timestark(j), j=1,mxtbcmhd)
            write (ncrt, 90)  istark, (timestark(j), j=1,mxtbcmhd)
   90       format (/ ' subroutine STARK reports:' /
     .                ' ERROR: istark =', i5       /
     .                ' timestark = '              /
     .                  2x, (5(2x, 1pe12.4)))
            call STOP ('subroutine STARK: unspecified problem', 101)
          end if
          write  (nout, 100)  time, timestark(1), timestark(istark)
          write  (ncrt, 100)  time, timestark(1), timestark(istark)
  100     format (' subroutine STARK reports:'                        /
     .            ' time is outside range of MSE values'              /
     .            ' calculations will not be done, time =',   1pe12.6 /
     .            ' timestark(1), timestark(istark) =', 2(2x, 1pe12.6))
          go to 20
   10 if (jj .eq. 1) then
        write (nout, 100)  time, timestark(1), timestark(istark)
        write (ncrt, 100)  time, timestark(1), timestark(istark)
        go to 20
      end if
c
c timestark(jj-1) .lt. time .le. timestark(jj)
c
      jjm1  = jj-1
      dt    = time-timestark(jjm1)
      dtime = timestark(jj)-timestark(jjm1)
      do j=1,nstark
       if (fwtstark(j,jj) .ne. 0.0 .and. fwtstark(j,jjm1) .ne. 0.0) then
              dval = fwtstark(j,jj)-fwtstark(j,jjm1)
              fwtstarkexptl(j) = fwtstark(j,jjm1)+dt*dval/dtime
c
              dval = tstark(j,jj)-tstark(j,jjm1)
              tstarkexptl(j) = tstark(j,jjm1)+dt*dval/dtime
c
              dval = sigstark(j,jj)-sigstark(j,jjm1)
              sigstarkexptl(j) = sigstark(j,jjm1)+dt*dval/dtime
c
              dval = rstark(j,jj)-rstark(j,jjm1)
              rstarkexptl(j) = rstark(j,jjm1)+dt*dval/dtime
c
              dval = zstark(j,jj)-zstark(j,jjm1)
              zstarkexptl(j) = zstark(j,jjm1)+dt*dval/dtime
c
              dval  = a1stark(j,jj)-a1stark(j,jjm1)
              a1(j) = a1stark(j,jjm1)+dt*dval/dtime
c
              dval  = a2stark(j,jj)-a2stark(j,jjm1)
              a2(j) = a2stark(j,jjm1)+dt*dval/dtime
c
              dval  = a3stark(j,jj)-a3stark(j,jjm1)
              a3(j) = a3stark(j,jjm1)+dt*dval/dtime
c
              dval  = a4stark(j,jj)-a4stark(j,jjm1)
              a4(j) = a4stark(j,jjm1)+dt*dval/dtime
c
       else
         tstarkexptl  (j) = 0.0
         sigstarkexptl(j) = 0.0
         rstarkexptl  (j) = 0.0
       end if
      end do
c
c the corresponding theoretical values are now found:
c
c --- first set up bicubic spline array for psi:
c
      call my_ibcccu(psi,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,wnoperm,ier)
c
c --- next get psi at the (R,Z) values from mse:
c
      jlo = 0
      do j=1,nstark
          if (rstarkexptl(j) .ne. 0.0) then
              call dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                     rstarkexptl(j),zstarkexptl(j),pds,ier,3)
              psiloc = pds(1)
              br     =  pds(3) / rstarkexptl(j)       ! radial   B field
              bz     = -pds(2) / rstarkexptl(j)       ! vertical B field
c
c             get toroidal b field:
c
              call tableintrp (psir, nj, psiloc, jlo) ! find psiloc
c                                                       in psir
              if (jlo .eq. 0 .or. jlo .eq. nj) then
                write (nout, 110) rstarkexptl(j),psiloc,psir(1),psir(nj)
                write (ncrt, 110) rstarkexptl(j),psiloc,psir(1),psir(nj)
  110           format (' subroutine STARK reports:'               /
     .                  ' psivalue for stark at r =', 1pe14.6      /
     .                  ' with psival =',1pe14.6                   /
     .                  ' is outside the range psir(1)...psir(nj)' /
     .                  ' psir(1) =', 1pe14.6, '  psir(nj) =', 1pe14.6)
                tstarkcalc(j) = 0.0
              else
                dval          = fcap(jlo+1)-fcap(jlo)
                dpsi          = psir(jlo+1)-psir(jlo)
                fcapa         = fcap(jlo)+(psiloc-psir(jlo))*dval/dpsi
                bt            = rmajor*btor/(fcapa*rstarkexptl(j))
                tstarkcalc(j) = a1(j)*bz/(a2(j)*bt+a3(j)*br+a4(j)*bz)
              end if
          else
              tstarkcalc(j)   = 0.0
          end if
      end do
c
c --- next get chi-square:
c
      chisqstarktot = 0.0
      do j=1,nstark
          if (rstarkexptl(j) .ne. 0.0) then
              chisqstark(j) = fwtstarkexptl(j)*((tstarkexptl(j)
     .                      - tstarkcalc(j))/sigstarkexptl(j))**2
              chisqstarktot = chisqstarktot+chisqstark(j)
          end if
      end do
c
      return
c
   20 do j=1,nstark
        rstarkexptl(j) = 0.0    ! signifies no data at this time
      end do
      return
c
      end

      subroutine symmetrize (psi1d, nw, nh)
c
      implicit none
c
c ----------------------------------------------------------------------
c     this subroutine makes the solution up/down symmetric about the
c     point jmid = nh/2 +1
c     it is assumed that nh = 2**n+1 for some integer n
c ------------------------------------------------------------------ HSJ
c
      integer  nw, nh, jmid, jmidm1, klow, kup, jcon, jp, jmidp1, j, i
      real*8   psi1d(*)
c
      jmidm1 = nh/2
      jmid   = jmidm1+1
      jmidp1 = jmid+1
      jcon   = jmidm1/(1-jmidm1)
      do   i=1,nw
        do j=jmidp1,nh
          jp          = jcon*(j-nh)+1
          klow        = (i-1)*nh+jp
          kup         = (i-1)*nh+j
          psi1d(klow) = psi1d(kup)
        end do
      end do
      return
c
      end

      subroutine tridiag (a, b, c, r, u, wk, n)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- solve the tridiagonal system
c                   c*u = r
c --- where c has lower diagonal a
c               diagonal       b
c               upper diagonal c
c --- pivoting is not implemented here so c should be diagonally
c --- dominant. this is NOT checked for in this routine.
c --- this routine requires approximately 5*n operations
c --- (multiply, divide).
c --- use LINPACK routine sgtsl if matrix is not diagonally dominant
c --- input
c  a(i)
c  b(i)
c  c(i)            i = 1,2....n
c                  are the diagonal elements
c  r(i)            is the rhs
c  n               number of equations
c  wk              work vector length n
c --- output:
c  u(i)            solution vector,i = 1,2,...n
c ------------------------------------------------------------------ HSJ
c
      dimension  a(n), b(n), c(n), r(n), wk(n), u(n)
c
      bet = b(1)
      u(1) = r(1)/bet
      do j=2,n
        wk(j) = c(j-1)/bet
        bet = b(j)-a(j)*wk(j)
        u(j) = (r(j)-a(j)*u(j-1))/bet
      end do
c
c --- now back sustitute
c
      do j=n-1,1,-1
        u(j) = u(j)-wk(j+1)*u(j+1)
      end do
      return
c
      end

      subroutine tridiagHSJ (a, b, c, r, u, wk, n, ierr)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- solve the tridiagonal system
c                   c*u = r
c --- where c has lower diagonal a
c               diagonal       b
c               upper diagonal c
c --- pivoting is not implemented here so c should be diagonally
c --- dominant. this is NOT checked for in this routine.
c --- this routine requires approximately 5*n operations
c --- (multiply, divide).
c --- use LINPACK routine SGTSL if matrix is not diagonally dominant
c --- input
c  a(i)
c  b(i)
c  c(i)            i = 1,2....n
c                  are the diagonal elements
c  r(i)            is the rhs
c  n               number of equations
c  wk              work vector length n
c --- output:
c  u(i)            solution vector,i = 1,2,...n
c  ierr             =0 no problerm
c                   =1 encountered zero pivot,solution is hossed
c ------------------------------------------------------------------ HSJ
c
      dimension  a(n), b(n), c(n), r(n), wk(n), u(n)
c
      ierr = 0
      bet  = b(1)
      if (bet .eq. 0.0) then
        ierr = 1
        return
      end if
      u(1) = r(1)/bet
      do j=2,n
        wk(j)  = c(j-1)/bet
        bet    = b(j)-a(j)*wk(j)
        if (bet .eq. 0.0) then
          ierr = 1
          return
        end if
        u(j)   = (r(j)-a(j)*u(j-1))/bet
      end do
c
c --- now back sustitute
c
      do j=n-1,1,-1
        u(j) = u(j)-wk(j+1)*u(j+1)
      end do
      return
c
      end

      subroutine volcalc (rcontr, zcontr, ncontr, rma,zma, volume, area)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c get the plasma volume, using the contour given by
c rcontr(j),zcontr(j), j = 1,2..ncontr
c use  vol = 2 * pi * contour integral(r*z*dr)
c use area = contour integral(z*dr)
c where the contour is given in rcontr, zcontr
c ----------------------------------------------------------------------
c
      dimension  rcontr(*), zcontr(*)
      data       twopi /6.28318530718 /
c
      xym    = rcontr(1)*zcontr(1)
      area   = 0.0
      xyma   = zcontr(1)
      zzm    = zma-zcontr(1)
      volsum = 0.0
      do j=2,ncontr
        xyp    = rcontr(j)*zcontr(j)
        xypa   = zcontr(j)
        dx     = rcontr(j)-rcontr(j-1)
        volsum = volsum+(xyp+xym)/2.0*dx
        area   = area+(xyma+xypa)/2.0*dx
        xym    = xyp
        xyma   = xypa
        zzp    = zma-zcontr(j)
        if (zzp*zzm .le. 0.0) then
          slope = (rcontr(j)-rcontr(j-1))/(zcontr(j)-zcontr(j-1))
          if (rcontr(j) .lt. rma)  rminzm = rcontr(j) + zzp*slope
          if (rcontr(j) .gt. rma)  rmaxzm = rcontr(j) + zzp*slope
        end if
        zzm = zzp
      end do
      volume = ABS (volsum)*twopi
      area   = ABS (area)
      return
c
      end

      subroutine wrteqdsk (psi, rmhdgrid, zmhdgrid, rmajor, btor,
     .                     torcur, pedge, nj, xsep, ysep, xcontr,
     .                     ycontr, ncontr, psilim, psiaxis, xax, yax,
     .                     vid, xlimiter, ylimiter, nlimiter, ishot,
     .                  ieqdsk, use_Bp0,nion, nprim, nimp, ncrt, nout)
c

c
c ----------------------------------------------------------------------
c  this routine creates an eqdsk file called g0.time.
c  after the standard eqdsk output,we append the transport information
c  necessary to use this eqdsk as a startup file.
c --- input (through argument list)
c  psi(i,j)                   psi in volt-sec
c  rmhdgrid(1..nw)
c  zmhdgrid(1..nh)            mhdgrid,meters
c  rmajor                     R0,meters
c  btor                       b field at R0,tesla
c  torcur                     total toroidal current,amps
c  nj                         size of vectors on transport grid
c                             (see vectors indexed by nj in this list)
c  pedge                      pressure at plasma boundary,nt/m**2
c  psilim                     psi on plasma boundary,volt-sec/rad
c  psiaxis                    psi on magnetic axis,volt-sec/rad
c  xax,yax                    R,z of magnetic axis,meters
c  xsep,ysep                  R,z of xpoint,meters
c  vid                        8 character version id
c
c  xlimiter(1..nlimiter)
c  ylimiter(1..nlimiter)
c  nlimiter                   the limiter description
c
c  xcontr(1..ncontr)
c  ycontr(1..ncontr)
c  ncontr                     The R,z points defining plasma boundary
c
c  ishot                      shot number
c                             ishot and time are used to generate an
c                             eqdsk name of the form gssssss.ttttt
c                             the name is always exactly 13 characters
c                             long. ssssss is the shot number with as
c                             many leading zeros as necessary and ttttt
c                             is the time in msec, with as many leading
c                             zeros as necessary (not true on CRAY; see below)
c
c  nprim          number of primary species
c  nimp           number of impurity species ( .ge. 0)
c  nion           (nion = nprim+nimp)
c  nti            number of ion temperatures (1 only in ONETWO)
c
c  ieqdsk                 =  1 create a new eqdsk file
c                         = -1 overwrite the old eqdsk file
c  use_Bp0                controls  pprim on axis 
c
c  through INCLUDE files:
c  INCLUDE file etc
c  bp(kj)             poloidal b field
c
c  INCLUDE file extra:
c  q(kj)
c
c  INCLUDE file fusion:
c  enasav(kj)                 alpha particle density
c  wasav(kj)                  alpha particle energy density
c
c  INCLUDE file soln:
c  u(kk,kj)                   the transport solution array
c  te             electron temp in keV.
c  ti             temp of ions (keV).
c  ene            electron density
c  en(nj,nion)    primary and impurity ion densities
c
c  INCLUDE file geom:
c  fcap(kj)         various geometric factors
c  gcap(kj)
c  hcap(kj)
c  r2capi(kj)
c  rcap(kj)
c  rcapi(kj)
c
c  INCLUDE file nub2:
c  enbsav(kj,ke,kb)           fast ion density
c  wbsav(kj,ke,kb)            fast ion stored energy density
c
c  INCLUDE file psig:
c  psival(kpsi)
c  rho(kpsi)
c
c  INCLUDE file rhog:
c  pprim(1..nj)               dp/dpsi, nt/(m**2-volt-sec) ( = amp/m**3)
c  ffprim(1...nj)             f*fprim, tesla,( = kg/(coul-sec))
c  pressb(1...nj)             beam pressure,nt/m**2
c  psir           poloidal flux function (flux/(2*pi)) on r-mesh.
c
c  INCLUDE file solcon:
c  time0                      start time of the analysis
c  time                       current transport time
c  dt                         current transport time step
c  dtt                        theta*dt
c  n                          number of transport steps completed
c  eqtim0                     time of last equilibrium calculation
c                             (eqtim0 is updated to time after a
c                              successful transport/mhd step)
c  INCLUDE file ions:
c  atw            mass number of ions (units of proton mass).
c  z              charge state of ions (units of electronic charge).
c  zeff
c
c  INCLUDE file nub
c  ibeam                    status flag for beam
c
c  INCLUDE file mesh:
c  r              rho, SQRT (toroidal flux/(pi*btor)), not normalized.
c  temporary storage required:
c  INCLUDE file bicube:
c  wnoperm                    vector of min length 10*nw
c
c  EVERYTHING IS IN MKS UNITS  (EXCEPT TEMPERATURES, WHICH ARE IN keV).
c  to conform to the eqdsk specification we must interpolate the parameters
c  from the transport grid to the MHD grid.
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE fusion
      USE ions
      USE solcon
      USE soln
      USE mhdpar
      USE nub 
      USE ename
      USE nub2
      USE extra
      USE mesh
      USE geom
      USE psig
      USE rhog
      USE soln2d,                     ONLY : ieq,ifixshap
      USE kinetic_efit,               ONLY : wrt_kin_efit_naml
      USE bicube
      USE etc
      USE flxav,                      ONLY : npsi
      USE  eq_pfprim,                 ONLY : nwg_pfprim



      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'imsl.i'
c      include 'etc.i'

c
      integer use_Bp0,set_grid,kine_stat
      dimension    psi(nw,nh), rmhdgrid(*), zmhdgrid(*),
     .             xlimiter(*), ylimiter(*), xcontr(*), ycontr(*)
      dimension    psi_eqdsk(nw), xpp(nw), xffp(nw), qpsinw(nw),
     .             xfprime(nw)
      dimension    sf(nw), sp(nw), spbb(nw), xax(2), yax(2), dum(kj)
      character*(*)    vid, ntitle*40
      logical      monotonic
      INTERFACE 
      SUBROUTINE check_monotonic (array, n, monotonic, incr)
        USE nrtype,   ONLY : I4B,DP
        INTEGER(I4B), INTENT(IN) :: n,incr
        REAL(DP),INTENT(IN),DIMENSION(:) :: array(n)
        LOGICAL, INTENT(OUT) ::    monotonic
        END SUBROUTINE check_monotonic
      END INTERFACE 
c
c
c --- check to see if psir is monotonic
c --- because psir is determined from bp, sometimes
c --- psir becomes non-monotonic near the magnetic axis;
c --- in this case we can not define the psi grid for
c --- the eqdsk, so diagnostic action is taken: 
c

      call check_monotonic (psir, nj, monotonic, 1)
     
c
c --- check enbsav, wbsav, enasav, wasav, for numbers that are too
c --- small in single precision: causes problems with some programs
c --- that try to read these eqdsks
c
      too_small = 1.0d-38
      do    ib=1,kb
        do  ie=1,ke
          do j=1,kj
             wbsav(j,ie,ib) = MAX ( wbsav(j,ie,ib), too_small)
            enbsav(j,ie,ib) = MAX (enbsav(j,ie,ib), too_small)
          end do
        end do
      end do
c
      kine_stat = 1
      psir_mono :    if (monotonic) then         ! psir is monotonic so write the eqdsk
      eqdskoldname = eqdskfilename
      xdim         =  rmhdgrid(nw) - rmhdgrid(1)
      ydim         =  zmhdgrid(nh) - zmhdgrid(1)
      ymideqd      = (zmhdgrid( 1) + zmhdgrid(nh))/2.0
c
c --- construct uniform, nw in length, psi grid for eqdsk:
      set_grid =1      ! only psi_eqdsk is determined in nwg_pfprim if = 1
      CALL nwg_pfprim(psi_eqdsk,sf,sp,xpp,xffp,nw,psiaxis,psilim, 
     .                      set_grid)
      set_grid =0

c
c --- interpolate pprim, ffprim, pressb and q onto eqdsk psi grid
c
      imslmd = 'wrteqdsk'
c
      call intrp (1, 1, psir,  pprim, nj, psi_eqdsk, xpp   , nw)
      call intrp (1, 1, psir, ffprim, nj, psi_eqdsk, xffp  , nw)
      call intrp (1, 1, psir, pressb, nj, psi_eqdsk, spbb  , nw)
      call intrp (1, 1, psir,      q, nj, psi_eqdsk, qpsinw, nw)

c
      rsign = 1.0
      if (qpsinw(1) .lt. 0.0)  rsign = -1.0
      call multpl1 (qpsinw, nw, rsign) ! output q as a positive quantity
c
c ----------------------------------------------------------------------
c integrate xpp and xffp to get pressure and f
c unfortunately the uniform psi grid spacing required by the eqdsk
c is typically too coarse (if nw=33 or even 65) to get the levelling
c off of the pressure profile typically found near the magnetic axis
c when the pressure profile is determined as a function of RHO .
c typically the first psi point away from the magnetic axis accounts for
c as many as 10 rho grid points in the same interval so the structure
c in press(rho) is lost in press(psi). The following calcuations either
c get the pressure and f (i.e., sp and sf) self-consistently by integrating
c the interpolated pprim and ffprim (i.e., xpp and xffp) to form sp and sf
c if comp_methd_eqdsk =0
c                                OR
c sp and sf are found by interpolation of press and fpsi from the rho
c to the psi grid and then xpp and xffp are determined by differentiating
c sp and sf, if comp_methd_eqdsk=1
c
c for fine enough grids the two methods should agree quite well but
c for the grids typically used in mhd calculations there is some
c difference between the methods.
c
c   HSJ
c ----------------------------------------------------------------------
c
      if (comp_methd_eqdsk .eq. 0) then
         !now determine sf,sp, because set_grid =0 in nwg_pfprim
         CALL nwg_pfprim(psi_eqdsk,sf,sp,xpp,xffp,nw,psiaxis,psilim, 
     .                      set_grid)
 
      elseif(comp_methd_eqdsk .gt. 0 ) then  
c
c         get pressure (sp) and f (sf) directly on the eqdsk psi grid
c         (i.e., psi_eqdsk) by interpolation from the transport psi grid
c         (i.e., psir) recall that psir is the psi grid corresponding to
c         the rho (i.e., r) grid.
c
          call intrp (1, 1, psir,  press, nj, psi_eqdsk, sp, nw)
c          call intrp (1, 1, psir,   fpsi, nj, psi_eqdsk, sf, nw)
c          replaced above with following 8/24/06 HSJ:
          call intrp (1, 1, psival, fpsi, npsi, psi_eqdsk, sf, nw)

c
c         now differentiate to get xpp
c
          call difydx (psi_eqdsk, sp, xpp, nw)
c
c         differentiate to get fprime
c
          call difydx (psi_eqdsk, sf, xfprime, nw)

          if (use_Bp0 .eq. 2) then ! get smooth values near axis by polynomial extrapolation
            npoints = 3  ! use lowest degree (limited to 10 in sub POLINT)
            call POLINT (psi_eqdsk(use_Bp0),xpp(use_Bp0),
     .                 npoints,psi_eqdsk(1),xpp(1),DY)

            call POLINT (psi_eqdsk(use_Bp0),xfprime(use_Bp0),
     .                 npoints,psi_eqdsk(1),xfprime(1),dy)
          endif


c
c         form f*fprime:
c
          do j=1,nw
            xffp(j)=sf(j)*xfprime(j)
          end do

      end if
c
      xax(2) = 0.0
      yax(2) = 0.0
      xdum1  = 0.0
      ipestg = 4
c
c ----------------------------------------------------------------------
c construct an eqdsk file name
c ----------------------------------------------------------------------
c

       call eqdsk_name(eqdskfilename)

c
c --- form a title
c

      itms            = ABS( 1000*time )

      IF(ifixshap == 1)THEN
         write (intfl, 9710)  ishot, itms
         read  (intfl, 8025)  ntitle
         ntitle = adjustl(ntitle)
      ELSE
         ntitle ='FIXBDRY EQDSK  99999 '       ! note chekgm also has this
      ENDIF
c
c --- get an available unit number for output
c
      call getioun(iounit,66)
c
c --- destroy the old eqdsk file if necessary
c
      if (ieqdsk .eq. -1 .and. ieq .ne. 0 ) call DESTROY (eqdskoldname)
c
c --- create the file
c
      open (unit = iounit, file = eqdskfilename, status = 'UNKNOWN',
     .       err = 500)
c
      dpsimag = psi_eqdsk( 1)
      dpsilim = psi_eqdsk(nw)
c
c --- finally do the eqdsk writes
c
      IF(ifixshap == 1)THEN
         write (iounit, 8000)  ntitle(1:len_trim(ntitle)), 
     .             vid(1:len_trim(vid)), ipestg, nw, nh
      ELSE
         write (iounit, 8001)  ntitle(1:len_trim(ntitle)), 
     .             vid(1:len_trim(vid)), ipestg, nw, nh
      ENDIF
      write (iounit, 8010)  xdim,ydim,rmajor,rmhdgrid(1     ),
     .                                       zmhdgrid(nh/2+1)
      write (iounit, 8010)  xax(1),yax(1),dpsimag,dpsilim,btor
      write (iounit, 8010)  torcur,psi_eqdsk(1),psi_eqdsk(nw),
     .                                           xax(1),xax(2)
      write (iounit, 8010)  yax(1),yax(2),psi_eqdsk(nw),xsep,ysep
      write (iounit, 8010) (sf(i),i=1,nw)
      write (iounit, 8010) (sp(i),i=1,nw)
      write (iounit, 8010) (xffp(i),i=1,nw)
      write (iounit, 8010) (xpp(i),i=1,nw)
      write (iounit, 8010) ((psi(i,j),i=1,nw),j=1,nh)
      write (iounit, 8010) (qpsinw(i),i=1,nw)
      write (iounit, 8020)  ncontr,nlimiter
      write (iounit, 8010) (xcontr(i),ycontr(i),i=1,ncontr)
      write (iounit, 8010) (xlimiter(i),ylimiter(i),i=1,nlimiter)



c
c ----------------------------------------------------------------------
c --- the following output is not part of the "official" eqdsk. It is
c --- what distinguishes eqdsk created in ONETWO from those created by EFIT:
c --- we put out enough information so that the eqdsk could be used as
c --- a restart file.
c --- use nti as a switch to indicate that
c ----------------------------------------------------------------------
c
      nti = -1
      write (iounit, 8030)  nj, nprim, nimp, nti,npsi
      nti = +1
      write (iounit, 8010) (atw(k), k=1,nion)
c
      do 700 k=1,nion
 700  write (iounit, 8010) (z(j,k),j=1,nj)
      write (iounit, 8010) (zeff(j),j=1,nj)
      write (iounit, 8010) (r(j),j=1,nj)
      write (iounit, 8010) (psir(j),j=1,nj)
      write (iounit, 8010) (curden(j),j=1,nj)
      write (iounit, 8010) (te(j),j=1,nj)
      do 710 k=1,nti
 710  write (iounit, 8010) (ti(j),j=1,nj)
      do 720 j=1,nj
 720  dum(j) = 1.0e6*ene(j)
      write (iounit, 8010)  (dum(j),j=1,nj)
      do 730 k=1,nprim
      do 731 j=1,nj
 731  dum(j) = 1.0e6*en(j,k)
 730  write (iounit, 8010) (dum(j),j=1,nj)
      if (nimp .ne. 0)   then
      do 740 k=nprim+1,nion
      do 741 j=1,nj
 741  dum(j) = 1.0e6*en(j,k)
 740  write (iounit, 8010) (dum(j),j=1,nj)
      end if
      write (iounit, 8010) (spbb(i),i=1,nw)
c
c --- write out beam pressure on transport grid.
c --- signal that this write was done by printing out nti = -1 above.
c
      write (iounit, 8010)  (pressb(j),j=1,nj)
      write (iounit, 8010)  ((u(i,j),i=1,nion+4),j=1,nj)
      realn     = n
      realibeam = ibeam
      write (iounit, 8010)  time0,time,dt,realn,eqtim0,dtt,realibeam
      write (iounit, 8010)  (fcap(j),j=1,nj)
      write (iounit, 8010)  (gcap(j),j=1,nj)
      write (iounit, 8010)  (hcap(j),j=1,nj)
      write (iounit, 8010)  (r2cap(j),j=1,nj)
      write (iounit, 8010)  (bp(j),j=1,nj)
      write (iounit, 8010)  (pprim(j),j=1,nj)
      write (iounit, 8010)  (ffprim(j),j=1,nj)
      write (iounit, 8010)  (psival(j),j=1,npsi)
      write (iounit, 8010)  (rho(j),j=1,npsi)
      write (iounit, 8010)  (q(j),j=1,nj)
      write (iounit, 8010)  (r2capi(j),j=1,nj)
      write (iounit, 8010)  (rcap(j),j=1,nj)
      write (iounit, 8010)  (((enbsav(i,j,k),i=1,nj),j=1,ke),k=1,kb)
      write (iounit, 8010)  (((wbsav(i,j,k),i=1,nj),j=1,ke),k=1,kb)
      write (iounit, 8010)  (enasav(i),i=1,nj)
      write (iounit, 8010)  (wasav(i),i=1,nj)
      write (iounit, 8010)  (press_alpha(i),i=1,nj)
      write (iounit, 8010)  (press_thermal(i),i=1,nj)
      write (iounit, 8010)  (rcapi(j),j=1,nj)
      write (iounit,FMT ='(" done with eqdsk write")')
c
 500  close (unit = iounit)
      call giveupus(iounit)


      end if psir_mono   ! psir monotonic branch

c --- restore enbsav, wbsav, enasav, wasav to zero if necessary
c --- required for Houlberg calc of temperatures in nclboot HSJ
c
c
      do    ib=1,kb
        do  ie=1,ke
          do j=1,kj
             if ( wbsav(j,ie,ib) .lt. 1.1*too_small)
     .                wbsav(j,ie,ib) =0.0
             if ( enbsav(j,ie,ib) .lt. 1.1*too_small)
     .                enbsav(j,ie,ib) =0.0
          end do
        end do
      end do
c

c    write input file for kinetic efit runs:
c      kine_message = eqdskfilename

      call wrt_kin_efit_naml(0,kine_stat)

      

      return
c
 8000 format (a, a,t49, 3i4, t73, 'ONETWO EQDSK')
 8001 format (a, a,t49, 3i4, t73, 'FIXBDRY EQDSK')

 8010 format (5e16.9)
 8020 format (2i5)
 8025 format (a)
 8030 format (4(2x, i5))
 9000 format ('.0000' , i1)
 9020 format ('.000'  , i2)
 9030 format ('.00'   , i3)
 9040 format ('.0'    , i4)
 9050 format ('.'     , i5)
 9500 format ('g00000', i1)
 9520 format ('g0000' , i2)
 9530 format ('g000'  , i3)
 9540 format ('g00'   , i4)
 9550 format ('g0'    , i5)
 9560 format ('g'     , i6)
 9700 format ('g0'    ,  a)
 9710 format (i12,i12, ' transport/mhd ')
c
      end








      subroutine wrt_tdem_eqdsk(time_eqdsk,eqdsk_name) 
c
c ----------------------------------------------------------------------
c  this routine creates an eqdsk file called g0.time
c  from a tdem netcdf eqdsk file.
c --- input (through argument list)
c time_eqdsk              Time at which eqdsk is to be written

c --- output through argment list:
c eqdsk_name              name of eqdsk created gu this subroutine
c ------------------------------------------------------------------ HSJ
c
      USE param,   only :  kj
      USE contour, only :  nrplas,nconmax,rsep,zsep
      USE mhdpar,  only  : nw,nh
      USE mhdgrid, only  : zmhdgrid
      USE ename ,  only  : intfl,eqdskfilename
      USE extra,   only  : toteqd
      USE tdem,    only  : qpsi_tdem,pprimpsi_tdem,ffprimpsi_tdem,
     .                     presspsi_tdem,fpsi_tdem
      USE yoka,    only  : ishot
      USE io,      only  : versid
      implicit  none


      character (len =*) eqdsk_name

      real * 8, allocatable,dimension(:,:) :: psi_t           !temp psi points
      real * 8, allocatable,dimension(:) :: rlim_eqd_t,zlim_eqd_t  !limiter points
      real * 8, allocatable,dimension(:) :: rp_eqd_t,zp_eqd_t !plasma limiting surface

      real * 8 time_eqdsk                  !in sec
      real * 8 btor_t                      !in tesla
      real * 8 xmaxis_t                    !mag axis in m
      real * 8 zmaxis_t                    !mag axis in m
      real * 8 psimagax_t                  !mag axis ,voltsec/rad
      real * 8 psilimitr_t                 !limiting lux surface,voltsec/rad
      real * 8 xdim_t                      !r  length  of grid,m
      real * 8 zdim_t                      !z  length of grid, m
      real * 8 rmajor_t                    !major radius,m
      real * 8 redgem_t                    !inside ede of R grid,m
      real * 8 zmid_t                      !midpoint of z grid,m
      real * 8 xsep_t,zsep_t               !(R,Z) of x point

      integer *4  nr_t,nz_t                !size of (R,Z) grid
      integer *4  nlim_eqd_t               !no limiter points
      integer *4  np_eqd_t                 !no points on plasma bdy


      integer *4 itms,iashot,ipestg,iounit,i,j
      character (len =6) itchar
      character (len =7) ischar
      character (len =40) ntitle


      allocate(psi_t(nw,nh),rlim_eqd_t(nconmax),zlim_eqd_t(nconmax),
     .             rp_eqd_t(nrplas),zp_eqd_t(nrplas) )

            call  ech_netcdf_interface (psi_t,nw,nh,
     .                         rlim_eqd_t,zlim_eqd_t,nlim_eqd_t,
     .                         rp_eqd_t,zp_eqd_t,np_eqd_t,
     .                         xdim_t,zdim_t,rmajor_t, redgem_t,zmid_t,  
     .                         xmaxis_t,zmaxis_t,psimagax_t,psilimitr_t,
     .                         btor_t,time_eqdsk)




c
c ----------------------------------------------------------------------
c construct an eqdsk file name
c ----------------------------------------------------------------------
c
      itms = ABS (1000 * time_eqdsk)
      if (itms .lt. 10) then
            write (intfl, 9000) itms
            read  (intfl, 8025) itchar
          else if (itms .lt. 100   ) then
            write (intfl, 9020) itms
            read  (intfl, 8025) itchar
          else if (itms .lt. 1000  ) then
            write (intfl, 9030) itms
            read  (intfl, 8025) itchar
          else if (itms .lt. 10000 ) then
            write (intfl, 9040) itms
            read  (intfl, 8025) itchar
          else if (itms .lt. 100000) then
            write (intfl, 9050) itms
            read  (intfl, 8025) itchar
          else
            itms = itms / 10
  200       if (itms .ge. 100000) then
              itms = itms / 10
              go to 200
            end if
            write (intfl, 9050) itms
            read  (intfl, 8025) itchar
      end if
c
c --- itchar is now a 6-character symbol; do the same with the shot number
c
      iashot = IABS (ishot)
      if (iashot .lt. 10) then
            write (intfl, 9500) iashot
            read  (intfl, 8025) ischar
          else if (iashot .lt. 100    ) then
            write (intfl, 9520) iashot
            read  (intfl, 8025) ischar
          else if (iashot .lt. 1000   ) then
            write (intfl, 9530) iashot
            read  (intfl, 8025) ischar
          else if (iashot .lt. 10000  ) then
            write (intfl, 9540) iashot
            read  (intfl, 8025) ischar
          else if (iashot .lt. 100000 ) then
            write (intfl, 9550) iashot
            read  (intfl, 8025) ischar
          else if (iashot .lt. 1000000) then
            write (intfl, 9560) iashot
            read  (intfl, 8025) ischar
          else
            iashot = iashot / 10
  300       if (iashot .ge. 1000000) then
              iashot = iashot / 10
              go to 300
            end if
            write (intfl, 9550) iashot
            read  (intfl, 8025) ischar
      end if
c
c --- ischar is now a 7-character symbol
c --- form the file name (on CRAY we must truncate to 8 characters)
c
      write (intfl, 9700)  itchar
      read  (intfl, 8025)  eqdskfilename
c
c --- form a title
c

      write (intfl, 9710)  ischar, itchar
      read  (intfl, 8025)  ntitle
c
c --- get an available unit number for output
c
      call getioun(iounit,66)
c
c --- destroy the old eqdsk file if necessary
c
c
c --- create the file
c
      if (eqdsk_name .eq. 'g0.restart' ) then !JMP START
      	eqdskfilename='g0.restart'
      end if !JMP END

      eqdsk_name = eqdskfilename
      open (unit = iounit, file = eqdskfilename, status = 'UNKNOWN',
     .       err = 500)
c

c 
c --- finally do the eqdsk writes
c
      ipestg = 4
      xsep_t = rsep/100.
      zsep_t = zsep/100.
      rp_eqd_t(:) = rp_eqd_t(:)*0.01
      zp_eqd_t(:) = zp_eqd_t(:)*0.01

      write (iounit, 8000)  ntitle, versid, ipestg, nw, nh
      write (iounit, 8010)  xdim_t,zdim_t,rmajor_t, redgem_t,
     .                                       zmhdgrid(nh/2+1)
      write (iounit, 8010)  xmaxis_t,zmaxis_t,psimagax_t,
     .                      psilimitr_t,btor_t
      write (iounit, 8010)  toteqd,psimagax_t,psilimitr_t,xmaxis_t,
     .                                                    xmaxis_t
      write (iounit, 8010)  zmaxis_t,zmaxis_t,psilimitr_t,xsep_t,zsep_t
      write (iounit, 8010) (fpsi_tdem(i),i=nw,1,-1)
      write (iounit, 8010) (presspsi_tdem(i),i=nw,1,-1)
      write (iounit, 8010) (ffprimpsi_tdem(i),i=nw,1,-1)
      write (iounit, 8010) (pprimpsi_tdem(i),i=nw,1,-1)
      write (iounit, 8010) ((psi_t(i,j),i=1,nw),j=1,nh)
      write (iounit, 8010) (qpsi_tdem(i),i=nw,1,-1)
      write (iounit, 8020)  np_eqd_t,nlim_eqd_t
      write (iounit, 8010) (rp_eqd_t(i),zp_eqd_t(i),i=1,np_eqd_t)
      write (iounit, 8010) (rlim_eqd_t(i),zlim_eqd_t(i),i=1,nlim_eqd_t)

 
      close (unit = iounit)
      call giveupus(iounit)
      deallocate(psi_t,rlim_eqd_t,zlim_eqd_t,
     .             rp_eqd_t,zp_eqd_t)


      return
c
 8000 format (a, a8, 3i4, t73, 'TDEM  EQDSK')
 8010 format (5e16.9)
 8020 format (2i5)
 8025 format (a)
 8030 format (4(2x, i5))
 9000 format ('.0000' , i1)
 9020 format ('.000'  , i2)
 9030 format ('.00'   , i3)
 9040 format ('.0'    , i4)
 9050 format ('.'     , i5)
 9500 format ('g00000', i1)
 9520 format ('g0000' , i2)
 9530 format ('g000'  , i3)
 9540 format ('g00'   , i4)
 9550 format ('g0'    , i5)
 9560 format ('g'     , i6)
 9700 format ('g0'    ,  a)
 9710 format (a, a, ' transport/mhd')
c
 500  call STOP('unable t open eqdsk file in sub wrt_tdem_eqdsk',1)
      end













      subroutine xpoint (psi,x,y,nw,nh,xsep,ysep,iknowxpt,cspln,n2cspln,
     .                   nh2,wk,ndwk,ispln,xlimiter,ylimiter,nlimiter,
     .                   iounit,ndimxc,ichklim,xcontr,ycontr,ncontr,
     .                   xax,yax,xc,yc,psisep,xmin,xmax,ymin,ymax,
     .                   xymin,xymax,yxmin,yxmax,wdum)

c
c ----------------------------------------------------------------------
c --- subroutine is similar to MAGAX except that a stationary point of psi
c --- is found rather than an extremum. The basic difference is that
c --- xpoint searches along a contour rather than over the MHD grid.
c --- the object is to find a saddle point in psi on or "near" the input
c --- contour. This saddle point defines the xpoint of the plasma if the
c --- contour that passes through the saddle point does not pass outside
c --- the limiter.
c --- Given a contour (usually found by bound) we search
c --- this contour for the minimum poloidal b field. This point is used
c --- as a starting guess to converge to dpsi/dx = 0 and
c ---                                    dpsi/dy = 0 simultaneously.
c --- convergence is assumed achieved if square of gradient is less
c --- than gradtol (defined locally below).
c --- if such a point is found within 2dx and 2dy of the starting guess
c --- then,on option (see ichklim),a contour passing through this point is
c --- generated (using subroutine CNTOUR) and the new contour points are checked
c --- to see if all of them are inside the limiter,using subroutine LIMITPTS.
c --- if so the new contour takes the place of the one found by bound.
c --- if the new contour passes outside the limiter it is rejected
c --- and the original contour found by bound is returned.
c
c --- input
c
c  psi(nw,nh)       the array of psi values
c  x(nw)
c  y(nh)             the MHD grid vectors
c  xsep
c  ysep            initial guess for xpoint location if known, see iknowxpt
c  iknowxpt       = 1 good intial guess for xpoint is input
c                     in xsep,ysep
c                 = 0 good initial guess is not available, find initial
c                     guess for xpoint by searching the
c                     input contour for min b poloidal
c  cspln(n2cspln,nw,nh2)   bicubic spline coeff. of psi if known
c                          if not known,set ispln = 0 to calculate it here
c  wk(ndwk)          temporary starge vector of minimum length
c                    ndwk = 2*nh*nw+2*max(nh,nw) ( required if ispln=0 )
c  ispln            = 1 bicubic spline array cspln is input
c                   = 0 bicubic spline array cspln must be calculated
c  xlimiter
c  ylimiter
c  nlimiter           the limiter points vectors
c  iounit             Fortran unit number for diagnostic output
c                     iounit = 0 suppresses output.
c  ndimxc             dimension of xcontr, ycontr, exactly as in dimension
c                     statement where these vectors are defined.
c                     this parameter is required in case cntour finds
c                     more points than can be stored,in which case an
c                     error exit must be taken.
c  ichklim            if an xpoint is found ichklim = 1 means trace
c                     a new contour corresponding the psi value at
c                     the x point and check this contour to make sure
c                     it is inside the limiter. As implemented here
c                     this method is way too slow (about 10 sec on 33 by
c                     65  grid) to use routinely. It exists primarily as
c                     a result of some testing done previously.
c                     if ichklim = 0 this check is not done,a new contour
c                     is not generated. instead the single point in the
c                     input contour which was closest to the x point is
c                     replaced by the xpoint and psisep is updated.
c                     this is the fast way to incorporate an x point in the
c                     contour. Note however that only one single point on
c                     the contour is changed. thus the remaining points may
c                     or may not adequately represent the new contour
c                     psi = psisep.this method takes about 20 msec.
c                     The final possibility is ichklim = -1. In this case
c                     it is assumed that the original input contour points
c                     are close to the required contour. We update these
c                     contour points to give psi = psisep exactly,BUT DO NOT
c                     CHECK FOR THE LIMITER. Hence using ichklim = -1 means
c                     we will get an accurate contour but for marginally
c                     diverted plasmas the contour found may pass slightly
c                     outside the limiter.this method is reasonably fast
c                     (about 50 msec).
c  xax
c  yax                location of magnetic axis (required if ichklim .ne. 0)
c  wdum(ndimxc)          work storage
c
c --- if ichklim = 1 the following temporary storage is required
c  xc(ndimxc)
c  yc(ndimxc)
c
c --- output
c
c  ispln = 1         indicates that cspln is now set
c  cspln           the bicubic spline coefficients of psi
c  psisep          the value of psi on the xpoint
c  iknowxpt      =1 indicates (acceptable) xpoint was found
c                =0 the xpoint was not found,because there is none,
c                  (within the required distance from the staring guess)
c                   or the xpoint contour passes outside the limiter.
c  xsep            final converged coordiantes of xpoint if found
c  ysep            returned as 0.0 if acceptable xpoint was not found.
c --- if iknowxpt =1 then the following parameters are reset
c  psisep
c  xmin
c  xmax
c  ymin
c  ymax
c  xymin
c  xymax
c  yxmax
c  yxmin
c  xcontr          (depends on setting of ichklim,see above)
c  ycontr
c  ncontr
c
c ------------------------------------------------------------------ HSJ
c
c
      USE replace_imsl,                   ONLY : my_ibcccu,my_dbcevl1
      implicit  integer (i-n), real*8 (a-h, o-z)
      include 'imsl.i'
c
      dimension x(nw),y(nh),psi(nw,nh),xc(*),yc(*),pds(6),
     .          ii(2),jj(2),cspln(2,nw,nh2),wk(ndwk),wdum(ndimxc),
     .          xcontr(ncontr),ycontr(ncontr),xlimiter(*),ylimiter(*)
c
      imslmd  = 'xpoint'
      gradtol = 1.0e-12
      psisave = psisep
      itry    = 1
      dx      = x(2) - x(1)
      dy      = y(2) - y(1)
      gradm   = 1.0d+100
      if (ispln .eq. 1)  go to 10
c
c     fit two-d spline to psi
c
      call my_ibcccu (psi, x, nw, y, nh, cspln, nw, wk, ier)
      ispln = 1
c
c     first find approximate location of x point by
c     finding minimum grad psi on contour points
c
   10 if ((    xsep .eq. 0.0) .and. (ysep .eq. 0.0))  go to 20
      if ( iknowxpt .eq. 1                         )  go to 30
c
c --- good guess not available; search the contour for starting guess
c
   20 do i=1,ncontr
        call my_dbcevl1(x,nw,y,nh,cspln,nw,xcontr(i),ycontr(i),
     .                  pds,ier,3)
        grad  = pds(2)**2+pds(3)**2
        gradm = MIN (grad,gradm)
        if (gradm .eq. grad)  k = i
      end do
c
      xsep = xcontr(k)
      ysep = ycontr(k)
      iknowxpt = 0
c
c --- final convergence to x point by Newton's method. find solution of
c --- dpsi/dx = 0 and dpsi/dy=0 using xsep,ysep as starting values.
c --- allow 15 iterations.
c
   30 do j=1,15
        call my_dbcevl1(x,nw,y,nh,cspln,nw,xsep,ysep,pds,ier,6)
        det = pds(5)*pds(6)-pds(4)*pds(4)
        if (det .eq. 0.0)  go to 60   ! can generate no more corrections
        xerr = (-pds(2)*pds(6)+pds(4)*pds(3))/det
        yerr = (-pds(5)*pds(3)+pds(2)*pds(4))/det
        xsep = xsep+xerr
        ysep = ysep+yerr
        erxpt = xerr**2+yerr**2
        if (              erxpt .lt. 1.0e-16)  go to 60
        if (pds(2)**2+pds(3)**2 .lt. gradtol)  go to 60
      end do
c
c --- no convergence
c --- Newton's method may have failed because start was too far away.
c --- get refined start and try one more time.
c --- NOTE: failure to find an x point is not a fatal error.
c
  120 if (itry .lt. 2)  go to 70
      if (iounit .ne. 0)  write (iounit, 1000)
 1000 format (' no convergence to xpoint in subroutine XPOINT')
      psisep   = psisave
      xsep     = 0.0
      ysep     = 0.0
      iknowxpt = 0
      return
c
   70 itry = itry+1
      if (iknowxpt .eq. 1) then
        itry = 1
        go to 20
      end if
c
c --- do a fine grid search near the starting point
c
      xsep  = xcontr(k)
      ysep  = ycontr(k)
   80 x1    = xcontr(k)-dx
      x2    = x1+dx
      y1    = ycontr(k)-dy
      y2    = y1+dy
      ii(1) = 1+(x1-x(1))/dx
      ii(2) = 1+(x2-x(1))/dx
      jj(1) = 1+(y1-y(1))/dy
      jj(2) = 1+(y2-y(1))/dy
      gradm = 1.0d+100
      dxx   = dx/20.0
      dyy   = dy/20.0
      nxx   = dx/dxx+1
      nyy   = dy/dyy+1
c
      do i=1,2
      do j=1,2
        do 110 ik=1,nxx
          xl = x(ii(i))+(ik-1)*dxx
          xl = MIN (xl,x(ii(i)+1))
          do 110 jk=1,nyy
            yl = y(jj(j))+(jk-1)*dyy
            yl = MIN (yl,y(jj(j)+1))
            call my_dbcevl1(x,nw,y,nh,cspln,nw,xl,yl,pds,ier,3)
            grad = pds(2)**2+pds(3)**2
            gradm = MIN (grad,gradm)
            if (grad .ne. gradm)  go to 110
            xsep = xl
            ysep = yl
  110   continue
      end do
      end do
c
      go to 30    ! try again with refined starting point
c
c --- an x point was found at (xsep,ysep). Now check its validity
c --- first check
c
   60 if (ABS (xsep-xcontr(k)) .gt. 2.0*dx .or.
     .    ABS (ysep-ycontr(k)) .gt. 2.0*dy) then
        if (itry .ge. 2)  go to 120
        itry = 2
        go to 80
      end if
      if (ichklim .eq. 1) then
c
c --- second check, is contour inside limiter?
c --- generate the contour corresponding to psisep
c
         delta_psi=0.0
        dang    = 5
        arcl    = 0.02
        bperr   = 0.03
        iaouto  = 1
        dxx     = dx/4.0
        dyy     = dy/4.0
        ierr    = 0
        psicopy = psisep      ! CNTOUR may change this value (iauto = 1)
        dx0     = 0.0
        dy0     = 0.0
        call cntour (xax,yax,psicopy,xemin,xemax,yemin,yemax,yxemin,
     .               yxemax,xyemin,xyemax,dang,arcl,bperr,dx0,dy0,
     .               xlimiter(nlimiter+1),xlimiter(nlimiter+2),
     .               ylimiter(nlimiter+1),ylimiter(nlimiter+2),
     .               iaouto,iautoc,xc,yc,ipts,x,nw,y,nh,cspln,n2cspln,
     .               nh2,iounit,ndimxc,ierr,wdum,0,delta_psi)
        if (ierr .ne. 0)  go to 450
c
c --- cntour has returned the contour points (xc(i),yc(i)),i = 1..ipts,
c --- which correspond to the value psi = psicopy.
c --- now check for relationship to limiter:
c
        thets   =  0.0
        thete   =  0.0
        rtan    = -1.0e30
        nchord  =  1
        xlimmin = xlimiter(nlimiter+1)
        xlimmax = xlimiter(nlimiter+2)
        ylimmin = ylimiter(nlimiter+1)
        ylimmax = ylimiter(nlimiter+2)
        do 300 j=1,ipts
          if (xc(j) .lt. xlimmin)  go to 350
          if (xc(j) .gt. xlimmax)  go to 350
          if (yc(j) .lt. ylimmin)  go to 350
          if (yc(j) .gt. ylimmax)  go to 350
          rs = xc(j)
          re = xc(j)+0.01
          zs = yc(j)
          ze = yc(j)+0.01
          call limitpts(rs,zs,thets,re,ze,thete,rtan,rsl,zsl,thetsl,
     .         rel,zel,thetel,nchord,ylimiter,
     .         xlimiter,nlimiter)
          rs = MIN (rsl,rel)
          re = MAX (rsl,rel)
          if (xc(j) .lt. rs .or. xc(j) .gt. re)  go to 350
          zs = MIN (zsl,zel)
          ze = MAX (zsl,zel)
          if (yc(j) .lt. zs .or. yc(j) .gt. ze)  go to 350
          go to 300
  350     if (itry .eq. 2)  go to 120
          itry = 2
          go to 80
  300   continue
c
c --- we fell through the loop so the new contour is ok.
c --- overwrite old points
c
      ncontr = ipts
      do 400 j=1,ncontr
        xcontr(j) = xc(j)
  400   ycontr(j) = yc(j)
        xmin      = xemin
        xmax      = xemax
        ymin      = yemin
        ymax      = yemax
        xymin     = xyemin
        xymax     = xyemax
        yxmin     = yxemin
        yxmax     = yxemax
        psisep    = psicopy
      end if
c
  450 if (ichklim .ne. 1) then
c
c --- replace xcontr(k),ycontr(k) with xsep,ysep and store the value of
c --- psi on the seperatrix in silimit
c
      xcontr(k) = xsep
      ycontr(k) = ysep
      psisep    = pds(1)
c
c --- modify the contour if requested
c
      if (ichklim .eq. -1) then
c
c       refine the k-1 contour points which come before the x point:
c
        ncontr1 = k-1
        if (ncontr1 .lt. 1)  go to 500
        call refinept (xcontr,ycontr,ncontr1,psisep,cspln,2,nh2,
     .                 nw,nh,x,y,xax,yax)
c
c       refine the k+1 to ncontr contour points after the xpoint:
c
  500   ncontr1 = ncontr-k
        call refinept (xcontr(k+1),ycontr(k+1),
     .                 ncontr1,psisep,cspln,2,nh2,nw,nh,x,y,xax,yax)
      end if
c
c --- search contour for new extremes. a general search is done rather
c --- than simple replacement for completely general geometries
c
      xmin =  1.0e10
      xmax = -1.0e10
      ymin =  1.0e10
      ymax = -1.0e10
      do j=1,ncontr
        xmax = MAX (xmax,xcontr(j))
        xmin = MIN (xmin,xcontr(j))
        ymax = MAX (ymax,ycontr(j))
        ymin = MIN (ymin,ycontr(j))
        if (ymin .eq. ycontr(j))  xymin = xcontr(j)
        if (ymax .eq. ycontr(j))  xymax = xcontr(j)
        if (xmin .eq. xcontr(j))  yxmin = ycontr(j)
        if (xmax .eq. xcontr(j))  yxmax = ycontr(j)
      end do
      end if
      iknowxpt = 1
      return
c
      end

      subroutine zlim (zero, nw, nh, limitr, xlim, ylim, x, y, iflag)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c**                                                                   **
c**     MAIN PROGRAM:  MHD FITTING CODE                               **
c**                                                                   **
c**     SUBPROGRAM DESCRIPTION:                                       **
c**          zlim determines whether points on the (x,y) grid are     **
c**          inside or outside of the boundary set by the limiters.   **
c**                                                                   **
c**     CALLING ARGUMENTS:                                            **
c**       zero............1 if inside and 0 otherwise                 **
c**       nw..............dimension of x                              **
c**       nh..............dimension of y                              **
c**       limitr..........number of limiter points                    **
c**       xlim............r coordinates of limiter                    **
c**       ylim............z coordinates of limiter                    **
c**       x...............r grid                                      **
c**       y...............z grid                                      **
c**       iflag...........1 convex geometry                           **
c**                       2 general geometry                          **
c**                                                                   **
c**     RECORD OF MODIFICATION:                                       **
c**          26/04/83..........first created                          **
c**          13/08/85..........iflag added                            **
c**                                                                   **
c***********************************************************************
c
c      dimension  zero(*),xlim(*),ylim(*),x(*),y(*),zzsum
      logical b,c,d,inside,bold
      INTEGER inout,nw,nh,limitr
      REAL*8 zero(*),zzsum,x(*),y(*),xlim(*),ylim(*)
c
      zero(1:nw*nh) = 0.0
      IF(limitr .lt. 4)THEN
 
         DO j=1,limitr
            print *,'j,xlim,ylim =',j,xlim(j),ylim(j)
         ENDDO
         CALL STOP('limitr in zlim error',1)
      ENDIF
      zzsum =0.0
      DO j=1,nh
         DO i =1,nw
            k=(i-1)*nh+j
            call pnpoly(x(i),y(j),xlim,ylim,limitr,inout)
            if(inout == 1)zero(k) = 1.0
            zzsum = zzsum + zero(k)
         ENDDO
      ENDDO


      RETURN
      
      go to (10, 200) iflag
c
   10 kk = 0
      do 100 i=1,nw
      do 100 j=1,nh
      kk = kk + 1
      zero(kk) = 1.0
      ncross = 0
c
      do 20 k=1,limitr-1
        if ((ylim(k) .lt. y(j)) .and. (ylim(k+1) .lt. y(j)))  go to 20
        if (x(i) .eq. xlim(k))  go to 20
        t = x(i) - xlim(k)
        s = xlim(k+1) - x(i)
        if ((t*s) .lt. 0.0)  go to 20
        di = (ylim(k+1)-ylim(k)) / (xlim(k+1)-xlim(k))
        f = ylim(k) + di*(x(i)-xlim(k))
        if (f .lt. y(j))  go to 20
        ncross = ncross + 1
   20 continue
c
      mcross = 0.5 * ncross
      mcross =   2 * mcross
      if (ncross .eq. mcross)  zero(kk) = 0.0
  100 continue
      return
c
  200 kk = 0
      do 2000 i=1,nw
      do 2000 j=1,nh
        kk = kk+1
        d = .false.
        b = .true.
        n = 0
        inside = .false.
        bold = b
        do k=1,limitr-1
          c = .false.
          if (y(j) .le. ylim(k) .and. y(j) .ge. ylim(k+1)
     .         .or. y(j) .ge. ylim(k) .and. y(j) .le. ylim(k+1)) then
              c = .true.
              d = .true.
              n = n+1
          end if
          if (c .and.
     .      (y(j)-ylim(k))*(xlim(k+1)-xlim(k))-
     .      (ylim(k+1)-ylim(k))*(x(i)-xlim(k)) .gt. 0.0)
     .       b = .not.b
          if (n .eq. 2) then
              n = 0
              if (bold .eqv. b)  inside = .true.
              bold = b
          end if
        end do
        zero(kk) = 0.0
        if (inside .and. d)  zero(kk) = 1.0
 2000 continue
      return
c
      end
      
