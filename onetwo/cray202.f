
      subroutine curdencalc2 (curden, nw, nh, rmhdgrid, wzero, psi,
     .                        psivl, npsivl, u0, totc, pprime,
     .                        ffprime, jsymetric)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c calculate a new current density:
c    curden(i,j) = -R*dp/dpsi-f*dfdpsi/(u0*R)
c ----------------------------------------------------------------------
c
      dimension  curden(nw,nh), wzero(*), psi(nw,nh), rmhdgrid(*),
     .           pprime(*), ffprime(*), psivl(*)
c
      jlo  = 0
      totc = 0.0
c
      do j=jsymetric,nh
          do i=1,nw
            k0 = (i-1)*nh + j
            if (wzero(k0) .ne. 0.0 ) then
                call tableintrp(psivl,npsivl,psi(i,j),jlo)
                if (jlo .eq. 0 .or. jlo .eq. npsivl) then
                   ppri = 0.0
                   ffpri = 0.0
                   itest = 0
                   if (psivl(1) .le. psi(i,j) .and. psi(i,j)  .le.
     .                 psivl(npsivl)) itest = -1
                   if (itest .ne. 0)
     .             call STOP ('subroutine CURDENCALC2: problem #1', 35)
                else
                   itest = 0
                   if (psivl(jlo) .le. psi(i,j) .and. psi(i,j)  .le.
     .                 psivl(jlo+1))  itest = 1
                   if (itest .eq. 0)
     .             call STOP ('subroutine CURDENCALC2: problem #2', 126)
                   dpsi = psi(i,j)-psivl(jlo)
                   dpdpsi = (pprime(jlo+1)-pprime(jlo))/
     .                      (psivl(jlo+1)-psivl(jlo))
                   dfdpsi = (ffprime(jlo+1)-ffprime(jlo))/
     .                      (psivl(jlo+1)-psivl(jlo))
                   ppri = pprime(jlo)+dpdpsi*dpsi
                   ffpri = ffprime(jlo)+dfdpsi*dpsi
                end if
                curden(i,j) = -rmhdgrid(i)*ppri-ffpri/(u0*rmhdgrid(i))
            else
                curden(i,j) = 0.0
            end if
            totc = totc+curden(i,j)
          end do
      end do
      return
c
      end

      subroutine curtotcalc (rmhdgrid, nw, zmhdgrid, nh, cspln, n2cspln,
     .                       nh2, psi2d, pds, wnoperm, rplasbdry,
     .                       zplasbdry, nplasbdry, u0, pcuriter, psilow,
     .                       psiavg, psihigh)
c

      USE replace_imsl,                     ONLY : my_ibcccu,my_dbcevl1
      implicit none
      include 'imsl.i'
c
c ----------------------------------------------------------------------
c subroutine calculates the total plasma current by integrating the
c poloidal B field over the plasma boundary
c ----------------------------------------------------------------------
c


      integer  nw, nh, j, nplasbdry, n2cspln, nh2, ier, icalc
c
      real*8   rmhdgrid(*),zmhdgrid(*),rplasbdry(*),zplasbdry(*),
     .         pds(*),cspln(n2cspln,nw,nh2),psi2d(nw,nh),wnoperm(*),
     .         pcuriter,sumdlbar,psisum,psilow,psihigh,rbar,zbar,
     .         dlbar,bpbar,u0,psiavg
c
c --- get the bicubic spline coefficients:
c
      imslmd ='82 c202'
      call my_ibcccu (psi2d,rmhdgrid,nw,zmhdgrid,nh,cspln,
     .                nw,wnoperm,ier)
c
c --- do the integral over the plasma boundary:
c
      pcuriter =  0.0
      sumdlbar =  0.0
      icalc    =  3
      psisum   =  0.0
      psilow   =  1.0d100
      psihigh  = -1.0d100
c
      do j=1,nplasbdry-1
        rbar     = 0.5*(rplasbdry(j)+rplasbdry(j+1))
        zbar     = 0.5*(zplasbdry(j)+zplasbdry(j+1))
        call my_dbcevl1 (rmhdgrid,nw,zmhdgrid,nh,cspln,nw,rbar,zbar,
     .                pds,ier,icalc)
        dlbar    = SQRT ((rplasbdry(j+1)-rplasbdry(j))**2 +
     .                   (zplasbdry(j+1)-zplasbdry(j))**2)
        bpbar    = SQRT (pds(2)**2+pds(3)**2)/rbar
        pcuriter = pcuriter+bpbar*dlbar
        sumdlbar = sumdlbar+dlbar
        psisum   = psisum+pds(1)
        psilow   = MIN (psilow,pds(1))
        psihigh  = MAX (psihigh,pds(1))
      end do
c
      pcuriter = pcuriter / u0
      psiavg   =   psisum / (nplasbdry-1)
      return
c
      end

      subroutine cycsoltn (nw,nh,nwh,psi,p,psi1d,nj,pprim,ffprim,
     .                     psir,iounit,isgngren,ieq,
     .                     kside,rmhdgrid,zmhdgrid,wzero,zero,
     .                     omeq,curden2d,u0,xax,yax,xlimiter,
     .                     ylimiter,nlimiter,rcontr,zcontr,ncontr,
     .                     rplasbdry,zplasbdry,nplasbdry,cspln,
     .                     n2cspln,nh2,wnoperm,pcurrent,psincbcd,
     .                     boundrcntr,boundzcntr,nboundpts,psibdry,
     .                     toleq,iteq,outer_iterations,bdindex,
     .                     psidif,residout,pds,elongax,ideriv,
     .                     ieqfail,navgbd,nconmax,nwork,rplasmin,
     .                     anormal,tocur,psinctemp,psiaxis)
c
      USE io,only : n66,n77
      implicit none
c
c ----------------------------------------------------------------------
c --- driver for single cyclic reduction solution of fixed boundary calculations
c
c --- input
c
c  psinctemp(i)       i = 1,2..kside temp work storage
c --- output
c ------------------------------------------------------------------ HSJ
c
      logical  converged
c
      integer  i,j,nw,nh,iter,ipass,k0,nwh,nj,bdindex(*),
     .         ispln,jsymetric,limfag,nerr,ix,ncontr,nlimiter,
     .         iounit,isgngren,kside,nplasbdry,navgbd,
     .         ncrt,isignpsi,n2cspln,nh2,nboundpts,iteq,
     .         outer_iterations,ideriv,ieqfail,
     .         iknowax,imax,jmax,ind,nconmax,nwork,ieq,map
c
      real*8   psi(nw,nh),p(nw,nh),psi1d(nwh),pprim(nj),ffprim(nj),
     .         psir(nj),rmhdgrid(nw),zmhdgrid(nh),wzero(nwh),zero(nwh),
     .         omeq,curden2d(nw,nh),u0,darea,xax(*),yax(*),
     .         xlimiter(*),ylimiter(*),rcontr(*),zcontr(*),
     .         rplasbdry(*),zplasbdry(*),cspln(n2cspln,nw,nh2),
     .         wnoperm(nwork),pcurrent(nwh),psincbcd(kside),
     .         boundrcntr(*),boundzcntr(*),psibdry,toleq,
     .         psidif(*),xmin,xmax,ymin,ymax,rymin,rymax,zxmin,
     .         zxmax,boundpsilim,psilow,psiavg,psihigh,pds(*),
     .         dpsi,elongax,pcuriter,psiaxis,dumy,ptrace,
     .         radold,relerr,rma,zma,residout,rplasmin,
     .         rmin,totc,anormal,tocur,psinctemp(*),
     .         psimaxbd,psilocbd,omeqloc
c

c
      jsymetric = 1
      darea     = (rmhdgrid(2)-rmhdgrid(1))*(zmhdgrid(2)-zmhdgrid(1))
      ieqfail   = 0
c
c     start the outer iterations (ipass counts them)
c     (required because we change the boundary values of psi on
c     the rectangular grid iteratively):
c
      ipass     = 0
 1010 iter      = 0
      ipass     = ipass + 1
      converged = .false.
c
c    start the inner iteration here (iter counts them);
c    (required because j toroidal depends on psi)
c
 1000 iter = iter + 1
c
c --- save the previous psi for use in relaxation:
c
      do i=1,nw
        do j=1,nh
          p(i,j)    = psi(i,j)
          k0        = (i-1)*nh+j
          psi1d(k0) = psi(i,j)
        end do
      end do
c
c ----------------------------------------------------------------------
c   calculate the plasma current in AMPS
c   NOTE THAT MULTIPLICATION BY DAREA IS DONE HERE TO AVOID HAVING
c   TO DO IT IN SEVERAL OTHER PLACES. THUS PCURRENT IS NOT A CURRENT
c   DENSITY. IT IS THE TOTAL CURRENT ASSOCIATED WITH A FILAMENT OF
c   CROSS-SECTIONAL AREA DAREA
c   pprim and ffprim are defined on the psir(1...nj) psi grid.
c ----------------------------------------------------------------------
c
      call curdencalc2 (curden2d, nw, nh, rmhdgrid, wzero, psi,
     .                  psir, nj, u0, totc, pprim, ffprim, jsymetric)
      do j=1,nh
        do i=1,nw
          k0           = (i-1)*nh+j
          pcurrent(k0) = curden2d(i,j)*wzero(k0)*darea
        end do
      end do
c
c --- load psi1d with boundary values of psi stored in psincbcd:
c --- psincbcd itself is loaded in subroutine getmhdbc (i.e.,
c --- get MHD boundary conditions)
c
      call setbdry (psi1d, nw, nh, nwh, kside, psincbcd)
c
c --- solve the GS equation with the given current density
c
      call solvefbd (psi, psi1d, pcurrent, u0, isgngren, iounit)
c
c --- relax the solution:
c
      omeqloc = omeq
      if (ipass .gt. outer_iterations/2)  omeqloc = omeq * 0.5
      if ( iter .gt. 3) then
        do   j=1,nh
          do i=1,nw
            psi(i,j) = omeq*psi(i,j)+(1.0-omeq)*p(i,j)
          end do
        end do
      end if
c
c --- find the magnetic axis
c
      if (iter .eq. 1)  iknowax = 0
      ispln    =  0
      isignpsi = -1
      ncrt     =  6
      call magax (psi,rmhdgrid,zmhdgrid,nw,nh,isignpsi,iknowax,
     .            ncrt,ispln,zero,rma,zma,cspln,n2cspln,nh2,
     .            psiaxis,elongax)
      if (iknowax .eq. 0) then
        call dump_psi_values (psi, nw, nh, rmhdgrid, zmhdgrid, n77,
     .              isignpsi, zero, 0, map, .false., dumy, 0, ptrace)
        write (6, *) 'ERROR: the magnetic axis was not found'
        write (6, *) '       at inner, outer iterations =', iter, ipass
        write (6, *) '       this most likely was caused by unrealistic'
        write (6, *) '       total toroidal current or bt0'
        call STOP ('subroutine CYCSOLTN: problem #2', 106)
      end if
c
c --- get average value of psi on plasma boundary, get total current.
c
      pcuriter = 0.0
      call curtotcalc (rmhdgrid, nw, zmhdgrid, nh, cspln, n2cspln, nh2,
     .                 psi, pds, wnoperm, rplasbdry, zplasbdry,
     .                 nplasbdry, u0, pcuriter, psilow, psiavg, psihigh)
      if (pcuriter .le. 0.0) then
        ieqfail = 1
        return
      end if
c
c --- rescale the psir array to new psiaxis value
c --- here psibdry remains unchanged
c
      call psiscale (psiaxis, psibdry, psir, nj, isignpsi)
c
c --- get the relative error, relerr and check for convergence
c
      call concek (psi, p, psiaxis, psibdry, nw, nh, toleq, ind, relerr,
     .             imax, jmax)
      if (ind .eq. 1)  converged = .true.
      if (converged)  go to 2000
c
c --- go to next inner iteration if max iterations allows
c
      if (iter .le. iteq)  go to 1000
c
 2000 write (6, *)  'ipass, iter, residual ', ipass, iter, relerr
c
c --- for the initial equilibrium (ieq = 0) we have an eqdsk
c --- we use the boundary values of psi from this eqdsk so the
c --- following is not necessary
c
      if (ipass .lt. outer_iterations .and. ieq .ne. 0) then
c
c --- temporarily save psincbcd in psinctemp
c
      call copya (psincbcd, psinctemp, kside)
c
c --- set up new boundary condition on border of MHD grid based on
c --- current solution
c
      call psibd (psi, rmhdgrid, zmhdgrid, rplasbdry, zplasbdry,
     .            nplasbdry, psibdry, psi1d, bdindex, ideriv,
     .            psincbcd, navgbd, jsymetric)
c
c --- relax the boundary
c
      do j=1,kside
        psincbcd(j) = 0.5*(psincbcd(j)+psinctemp(j))
      end do
      if (.not. converged)  go to 1010  ! inner iterations not converged
c
c       inner iterations converged; now check if the values of psi on
c       the rectangular boundary of the MHD grid are also converged
c
        psimaxbd = -1.0
        do j=1,kside
          psilocbd = ABS ((psincbcd(j)-psinctemp(j))/psiaxis)
          psimaxbd = MAX (psimaxbd,psilocbd)
        end do
        if (psimaxbd .gt. residout)  go to 1010   ! boundary not conv.
      end if
c
      if (.not. converged) then
        write (6, *)  ' did not reach required residual of ', toleq
        write (6, *)  ' in allowed inner iterations of =', iteq
        write (6, *)  ' will continue with best solution found'
      end if
      if (psimaxbd .gt. residout .and. ieq .ne. 0) then
        write (6, *)  ' did not reach required residual of ', residout
        write (6, *)  ' on the rectangular boundary of the MHD grid'
        write (6, *)  ' in allowed number of passes =',outer_iterations
        write (6, *)  ' will continue with best solution found'
      end if
c
c --- done with the solution,(either converged or maximum outer iterations
c --- were done). next do some processing of the solution:
c
      write (6, *)  'converged at outer iteration number ', ipass
      write (6, *)  'with relative error ', relerr
      write (6, *)  'psilow, psiavg, psihigh ', psilow, psiavg, psihigh
      do j=1,nh
        do i=1,nw
          k0        = (i-1)*nh+j
          psi1d(k0) = psi(i,j)
        end do
      end do
      radold = rplasmin
      limfag = 1
      ix     = 1
      rmin   = rplasmin
      nerr   = -1
      call bound (psi1d,nw,nh,nwh,zero,rmhdgrid,zmhdgrid,rma,zma,rmin,
     .     ix,nlimiter,xlimiter,ylimiter,nerr,limfag,radold,nconmax,
     .     rcontr,zcontr,ncontr,dpsi,xmin,xmax,ymin,ymax,
     .     rymin,rymax,zxmin,zxmax,boundpsilim)
c
      if (nerr .ne. 0) then
        write (6, *)  'subroutine BOUND, called by FIXDBDRY, reports:'
        write (6, *)  '  boundary not found, nerr =', nerr
        write (6, *)  '  code must stop'
        call STOP ('subroutine CYCSOLTN: problem #1', 36)
      end if
      write (6, *)  'psibdry, psilim from BOUND ', psibdry, boundpsilim
c
c --- the boundary values found by bound may not be
c --- exactly the same as the given fixed plasma contour.
c --- save the points found by bound for comparison elsewhere:
c
      nboundpts = ncontr
      do j=1,nboundpts
        boundrcntr(j) = rcontr(j)
        boundzcntr(j) = zcontr(j)
      end do
      xax(1)  = rma    ! needed in rest of code
      yax(1)  = zma
      anormal = tocur/pcuriter
      return
c
      end

      subroutine extrplt (psi, coeff, map, nw, nh, jsymetric)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- this subroutine extrapolates the values of psi to the grid
c --- boundary, using only values of psi interior to the plasma.
c ----------------------------------------------------------------------
c
      integer map(*)
      real*8  psi(nw,nh), coeff(9,*)
c
c --- zero the solution outside the plasma
c
      do j=jsymetric,nh
        do i=1,nw
          k = (i-1)*nh+j
          if (map(k) .eq. 0)  psi(i,j) = 0
        end do
      end do
c
c     vertical values:
c
      do j=jsymetric,nh
          do i=2,nw-1
              k = (i-1)*nh+j
              if ( map(k) .gt. 0 .and. (map(k-nh) .eq. 0  .or.
     .             map(k+nh) .eq. 0) ) then
                  l = map(k)
                  all = coeff(6,l)
                  alr = coeff(8,l)
                  allr = all*alr*(all+alr)
                  psill = psi(i-1,j)
                  psic = psi(i,j)
                  psilr = psi(i+1,j)
                  dpsidr = -psilr*all*all+psic*(all*all-alr*alr)
     .                             +psill*alr*alr
                  dpsidr = -dpsidr/allr
                  dsqpsidrsq = -psilr*all+psic*(all+alr)-psill*alr
                  dsqpsidrsq = -dsqpsidrsq/allr
                  if (map(k-nh) .eq. 0) then
                        psi(1,j) = psi(i,j)-dpsidr*(i-1)
     .                                    +dsqpsidrsq*(i-1)*(i-1)
                  else
                        psi(nw,j) = psi(i,j)+dpsidr*(nw-i)
     .                          +dsqpsidrsq*(nw-i)*(nw-i)
                  end if
              end if
          end do
      end do
c
c --- horizontal values
c
      jstart = 2
      if (jsymetric .ne. 1)  jstart = jsymetric
      do j=jstart,nh-1
          do i=1,nw
              k = (i-1)*nh+j
              if ( map(k) .gt. 0 .and. (map(k-1) .eq. 0  .or.
     .             map(k+1) .eq. 0) ) then
                  l = map(k)
                  alt = coeff(7,l)
                  alb = coeff(9,l)
                  allz = alb*alt*(alb+alt)
                  psilt = psi(i,j+1)
                  psic = psi(i,j)
                  psilb = psi(i,j-1)
                  dpsidz = psilt*alb*alb+psic*(alt*alt-alb*alb)
     .                      -psilb*alt*alt
                  dpsidz = dpsidz/allz
                  dsqpsidzsq = -psilt*alb+psic*(alt+alb)-psilb*alt
                  dsqpsidzsq = -dsqpsidzsq/allz
                  if (map(k-1) .eq. 0) then
                      psi(i,1) = psi(i,j)-dpsidz*(j-1)
     .                           +dsqpsidzsq*(j-1)*(j-1)
                  else
                     psi(i,nh) = psi(i,j)+dpsidz*(nh-j)
     .                               +dsqpsidzsq*(nh-j)*(nh-j)
                  end if
              end if
          end do
      end do
      return
c
      end

      subroutine findomegae (residin,maxinit,iterations,updownsym,
     .                       jsymetric,nw,nh,map,psi1d,coeff,deltar,
     .                       deltarzsq,deltarsq,interior,optomegae,
     .                       pds,rmhdgrid,psibdry,residmax,psi2d,
     .                       zmhdgrid)
c
      implicit none
c
c ----------------------------------------------------------------------
c find the optimum relaxation parameter for SOR solution method,
c --- in the exterior region between the plasma and the edge of
c --- the rectangular grid.
c
c --- input (most of this input is passed to subroutine outer_solution)
c
c  residin               maximum residual allowed. convergence
c                        will be assumed achieved if the largest
c                        residual is less than residin.
c  maxinit               maximum number of SOR iterations allowed
c                        for the right omega we expect the iterations
c                        to converge in about 2*(nw+nh) iterations
c                        (but also depends on value of residin)
c  psi1d(k)              for k corresponding to grid points on the
c                        boundary  of the rectangular grid,psi1d(k)
c                        is input (in volt*sec/rad)
c                        For these values of k psi1d remains unchanged
c                        on output.
c  pds(6)                used  as work space.
c  psibdry               value of psi on plasma boundary
c  rmhdgrid(i)           i = 1,2..nw the major radius values
c  zmhdgrid(i)           i = 1,2..nh the vertical grid values
c  updownsym             logical,ture if up/down sym solution is to be found
c  nw,nh                 rectangular grid size (for this method of
c                        solution no restriction is placed on nw,nh.
c                        However if updown symetric solution is requested
c                        then it is assumed that nh is odd !!!!!!!!!
c  jsymetric             defined
c  map(k)                defined
c  coeff(9,*)            defined
c  deltar                defined
c  deltarzsq             defined
c  deltarsq              defined
c  interior         logical,set to false for this routine
c
c --- output
c  optomegae           optimum relaxation factor
c
c ----------------------------------------------------------------------
c
      logical  interior, updownsym
      integer  nw, nh, map(*), iterations, maxinit, jsymetric,
     .         i, j, k, ih, ic, itermaxpf
      real*8   residmax,residin,psi1d(*),pds(*),coeff(9,*),optomegae,
     .         deltar,deltarzsq,deltarsq,pi,rmhdgrid(*),psibdry,
     .         psi2d(nw,nh),omegamax,omegamin,delomega,delomegamin,
     .         omega,omegah,omegac,rhojac,rhoest,optest,rhosor,
     .         drdzsq,zmhdgrid(nh), two
c
      data     pi /3.14159265358979/
c
c --- we know omega decreases linearly toward the optimum omega
c --- starting from 2.0. hence the following grid search is done:
c
      two = 2.0
      write (6, *)
     .      'optimal relaxation parameter search for exterior problem:'
      itermaxpf   = nw*nh  ! max number of inner iters done in sorpicard
      if (updownsym)  itermaxpf = itermaxpf/2
      delomega    = 0.10   ! starting coarse grid search interval
      delomegamin = 0.005  ! optimal omega within omega +/- delomegamin
      omegamax    = 2.0
      omegamin    = 1.0    ! accelerated conv. only between these limits
      write (6, *)  'initial omegamax = ', omegamax
      write (6, *)  'initial omegamin = ', omegamin
  200 ih    = 0
      omega = omegamax
      do while (omega .ge. omegamin-1.0e-7)
          do j=jsymetric,nh
            do i=1,nw
              k        = (i-1)*nh+j
              psi1d(k) = psi2d(i,j)  ! start with same boundary values
            end do
          end do
          call outer_solution (residin,maxinit,iterations,
     .                     updownsym,jsymetric,nw,nh,map,psi1d,coeff,
     .                     deltar,deltarzsq,deltarsq,interior,omega,pds,
     .                     rmhdgrid,psibdry,residmax)
          ic     = iterations
          omegac = omega    ! current values
          write (6, *)  'iters, omega =', ic, omegac
          if (ih  .gt. 0 .and. ic .gt. ih)  go to 100
          omegah = omegac
          ih     = ic
          omega  = omega - delomega
      end do
c
  100 write (6, *)  ' '
      write (6, *)  ' '
c
c --- optimal omega is between omegac and omegah+delomega
c --- (it is omegah+delomega because we may have stepped past the minimum).
c --- more refined search is now done
c
      optomegae = omegah
      omegamax  = omegah+delomega
      omegamax  = MIN (omegamax, two   )
      omegamin  = MAX (omegamin, omegac)
      delomega  = delomega * 0.5
      write (6, *)  'new omegamax =', omegamax
      write (6, *)  'new omegamin =', omegamin
      if (delomega .gt. delomegamin)  go to 200
c
c --- done, now print some informative output
c
      rhojac = (2.0 / optomegae - 1.0)**2
      if (rhojac .le. 1.0) then
        rhojac = SQRT (1.0 - rhojac)    ! Jacobi spectral radius
      else
        rhojac = 1.0
      end if
c
      drdzsq = (rmhdgrid(2)-rmhdgrid(1))/(zmhdgrid(2)-zmhdgrid(1))
      drdzsq = drdzsq*drdzsq
      rhoest = COS (pi/nw) + drdzsq * COS (pi/nh)
      rhoest = rhoest / (1.0 + drdzsq)
      if (rhoest .gt. 1.0)
     .rhoest = 1.0
      optest =     2.0 / (1.0 + SQRT (1.0-rhoest**2))
      rhosor = (rhojac / (1.0 + SQRT (1.0-rhojac**2)))**2
c
      write (6, *)  ' done with optimal relaxation parameter search'
      write (6, *)  ' spectral radius of SOR iteration matrix =', rhosor
      write (6, *)  ' optomegae =', optomegae
      write (6, *)
     .        '    spectral radius of Jacobi iteration matrix =', rhojac
      write (6, *)
     .        ' approximate radius of Jacobi iteration matrix =', rhoest
      write (6, *)  ' approximate optomegae based on approximate'
      write (6, *)  ' radius of Jacobi matrix =',optest
      return
c
      end

      subroutine findomegai (residin,maxinit,iterations,
     .                       updownsym,jsymetric,nw,nh,map,psi1d,coeff,
     .                       optomegai,rmhdgrid,psibdry,residmax,psi2d,
     .                       zmhdgrid,source)
c
      implicit none
c
c ----------------------------------------------------------------------
c find  the optimum relaxation parameter for SOR solution method,
c --- in the interior region of the plasma.
c
c --- input (most of this input is passed to subroutine INNER_SOLUTION)
c
c  residin               maximum residual allowed. convergence
c                        will be assumed achieved if the largest
c                        residual is less than residin.
c  maxinit               maximum number of SOR iterations allowed
c                        for the right omega we expect the iterations
c                        to converge in about 2*(nw+nh) iterations
c                        (but also depends on value of residin)
c  psibdry               value of psi on plasma boundary
c  rmhdgrid(i)           i = 1,2..nw the major radius values
c  zmhdgrid(i)           i = 1,2..nh the vertical grid values
c  updownsym             logical,ture if up/down sym solution is to be found
c  nw,nh                 rectangular grid size (for this method of
c                        solution no restriction is placed on nw,nh.
c                        However if updown symetric solution is requested
c                        then it is assumed that nh is odd !!!!!!!!!
c  jsymetric             defined
c  map(k)                defined
c  coeff(9,*)            defined
c  source(k)             rhs
c
c --- output
c  optomegai           optimum relaxation factor
c
c ----------------------------------------------------------------------
c
      logical updownsym
      integer nw, nh, map(*), iterations, maxinit, jsymetric,
     .        i, j, k, ih, ic, itermaxpf
      real*8  residmax, residin, psi1d(*), coeff(9,*), optomegai,
     .        pi, source(*), rmhdgrid(*), psibdry, psi2d(nw,nh),
     .        omegamax, omegamin, delomega, delomegamin, omega,
     .        omegah, omegac, rhojac, rhoest, optest, rhosor, drdzsq,
     .        zmhdgrid(nh), two
c
      data    pi /3.14159265358979/
c
c --- we know omega decreases linearly toward the optimum omega
c --- starting from 2.0. hence the following grid search is done.
c
      two = 2.0
      write (6, *)
     .      'optimal relaxation parameter search for interior problem:'
      itermaxpf = nw*nh    ! max number of inner iters done in sorpicard
      if (updownsym)  itermaxpf = itermaxpf / 2
      delomega    = 0.10   ! starting coarse grid search interval
      delomegamin = 0.005  ! optimal omega within omega +/- delomegamin
      omegamax    = 2.0
      omegamin    = 1.0    ! accelerated conv. only between these limits
      write (6, *)  'initial omegamax = ', omegamax
      write (6, *)  'initial omegamin = ', omegamin
  200 ih    = 0
      omega = omegamax
      do while (omega .ge. omegamin-1.0e-7)
          do j=jsymetric,nh
            do i=1,nw
              k          = (i-1)*nh + j
              psi1d(k)   = 0.0           ! start with same initial guess
              psi2d(i,j) = 0.0
            end do
          end do
          call inner_solution (residin,maxinit,iterations,
     .                         updownsym,jsymetric,nw,nh,map,psi1d,
     .                         coeff,omega,psibdry,residmax,source)
          ic     = iterations
          omegac = omega    ! current values
          write (6, *)  'iters, omega =', ic, omegac
          if (ih  .gt. 0 .and. ic .gt. ih)  go to 100
          omegah = omegac
          ih     = ic
          omega  = omega-delomega
      end do
c
  100 write (6, *)  ' '
      write (6, *)  ' '
c
c --- optimal omega is between omegac and omegah+delomega
c --- (it is omegah+delomega because we may have steped past the
c --- minimum).
c --- more refined search is now done.
c
      optomegai = omegah
      omegamax  = omegah + delomega
      omegamax  = MIN (omegamax, two   )
      omegamin  = MAX (omegamin, omegac)
      delomega  = delomega * 0.5
      write (6, *)  'new omegamax =',omegamax
      write (6, *)  'new omegamin =',omegamin
      if (delomega .gt. delomegamin)  go to 200
c
c --- done, print some informative output:
c
      rhojac = (2.0 / optomegai-1.0)**2
      if (rhojac .le. 1.0) then
          rhojac = SQRT (1.0 - rhojac)    ! Jacobi spectral radius
      else
          rhojac = 1.0
      end if
      drdzsq = (rmhdgrid(2)-rmhdgrid(1))/(zmhdgrid(2)-zmhdgrid(1))
      drdzsq = drdzsq * drdzsq
      rhoest = COS (pi/nw) + drdzsq * COS (pi/nh)
      rhoest = rhoest/(1.0 + drdzsq)
      if (rhoest .gt. 1.0)  rhoest = 1.0
      optest = 2.0 / (1.0 + SQRT (1.0-rhoest**2))
      rhosor = (rhojac/(1.0 + SQRT (1.0-rhojac**2)))**2
      write (6, *)  ' done with optimal relaxation parameter search '
      write (6, *)  ' spectral radius of SOR iteration matrix =',rhosor
      write (6, *)  ' optomegai =',optomegai
      write (6, *)
     .        '    spectral radius of Jacobi iteration matrix =', rhojac
      write (6, *)
     .        ' approximate radius of JACOBI iteration matrix =', rhoest
      write (6, *)  ' approximate optomegai based on approximate'
      write (6, *)  ' radius of Jacobi matrix =', optest
      return
c
      end
      subroutine fixdbdry (iteq, toleq, iounit, indeq, relerr, tocur,
     .                     ieq, omeq,ic)
c
c
c ----------------------------------------------------------------------
c --- fixed boundary equilibrium calculation module
c --- this subroutine is set up to solve  del-star (psi) = u0*r*jphi
c --- where jphi( = curden(i,j)) = -R*dp/dpsi-f*dfdpsi/(u0*r)
c --- psilim is greater than psiaxis for positive current flow.
c
c --- description of the input:
c --- units are MKS
c
c  iteq                     max number of inner iterations
c  toleq                    relative error for convergence
c  omeq                     relaxation parameter, 0.0 .le. omeq .le. 1.0
c  tocur                    total current, amps
c  bt0                      toroidal b field, tesla
c  psidifmax                relative error for convergence of outer iterations
c  toleqmin                 relative error for convergence of inner iterations
c  INCLUDE file contour.i:
c  rplasbdry(i)
c  zplasbdry(i)               i = 1,2...nplasbdry defines the fixed plasma
c  nplasbdry                  boundary.
c  ic                        pprim iteration counter passed to using_toq
c  INCLUDE file mhdcom.i:
c  psibdry                    value of psi on plasma boundary
c
c ------------------------------------------------------------------ HSJ
c

      USE param
      USE io
      USE contour
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE numbrs
      USE toq_12_interface, only : toq_drive,npsi_toq,
     .                              toqfile,psinarray,ffparray,pparray
      USE constnts
      USE rhog
      USE mhdcom
      USE bicube
      USE gpsi
      USE mhdbcdtn
      implicit  integer (i-n), real*8 (a-h, o-z)
c      parameter (ipassmax = 60, itermaxp = 200)
      parameter (ipassmax = 120, itermaxp = 1000)
c
c --- ipassmax is max no of outer iterations
c --- itermaxp is max no of inner iterations
c

c      include 'mhdbcdtn.i'

c      include 'small.i'
c      include 'zerocom.i'

      include 'fixbdry.i'
c
      logical    converged, updownsym
      integer    outer_iterations,outeritused
      dimension  psidif(ipassmax), map(nwh), coeff(nwh,9),
     .           curden2d(nw,nh), source(nwh), psiprevious(nw,nh),
     .           psinocoil(kside), psidifbdext(ipassmax),
     .           psidifbd(itermaxp), psinctemp(kside)

c
c --- this module is called only if the fixed boundary MHD calculations
c --- are being used. In this case the Green's table is not required
c --- and we use some of the Green's table arrays as storage for pointers
c --- required in the variable finite difference scheme:
c --- some large arrays available are:
c --- gridfc(nwh,nfcoil),gridpc(nwh,nw)
c --- rsilpc(nsilop,nwh),rmp2pc(magpr2,nwh),rgowpc(nrogow,nwh)
c
c      equivalence (gridpc(1,1), coeff(1,1))    ! assumes nw .ge. 9
c      equivalence (gridfc(1,1), map(1))
c      equivalence (gridfc(1,2), curden2d(1,1))
c      equivalence (gridfc(1,3), source(1))
c      equivalence (gridfc(1,4), psiprevious(1,1))
c      equivalence (gridfc(1,5), psinocoil(1))
c      equivalence (gridfc(1,6), psidifbdext(1))
c      equivalence (gridfc(1,7), psidifbd(1))
c      equivalence (gridfc(1,8), psinctemp(1))
c
      data icall /0/
c
      if (nw .lt. 9) then               ! restriction due to equivalence
        write  (nout, 515) nw
        write  (ncrt, 515) nw
  515   format (' ERROR: subroutine FIXDBDRY cannot run, nw =', i5)
        call STOP ('subroutine FIXDBDRY: problem #1', 102)
      end if
c
      psidifmax        = 1.0e-4
      toleqmin         = toleq
      outer_iterations = ipassmax
      inner_iterations = iteq
      ideriv           =  1
      isgngren         = -1
      if (omeq .le. 0.0)
     .omeq             = 0.5
      toleqloc         = toleqmin
      converged        = .false.
      itermax          = inner_iterations
      ipass            = 0
      darea = (rmhdgrid(2)-rmhdgrid(1))*(zmhdgrid(2)-zmhdgrid(1))
c
      if (icall .eq. 0) then  ! do the following only on first call to..
c                             ..this routine. save the results so they
c                             ..are available on subsequent calls
            icall = icall + 1
            if (nlimiter .le. 0) then
              write  (iounit, 10) nlimiter
   10         format (' subroutine FIXDBDRY reports:'        /
     .                ' nlimiter is not set, nlimiter =', i5 /
     .                ' ONETWO must stop')
              call STOP ('subroutine FIXDBDRY: problem #2', 103)
            end if
c
c --- set up array zero ( = 1 inside limiter, =0 outside limiter)
c
            iflag = 1
            call zlim (zero,nw,nh,nlimiter,xlimiter,ylimiter,rmhdgrid,
     .                 zmhdgrid,iflag)
c
c --- set up array wzero (= 1 inside plasma, = 0 outside plasma)
c
            if (nplasbdry .le. 0) then
                write  (iounit, 11) nplasbdry
   11           format (' subroutine FIXDBDRY reports:'          /
     .                  ' nplasbdry is not set, nplasbdry =', i5 /
     .                  ' ONETWO must stop')
                call STOP ('subroutine FIXDBDRY: problem #3', 104)
            end if
            call zlim (wzero,nw,nh,nplasbdry,rplasbdry,zplasbdry,
     .                 rmhdgrid,zmhdgrid,iflag)
            rplasmin = rplasbdry(1)
            rplasmax = rplasbdry(1)
            zplasmin = zplasbdry(1)
            zplasmax = zplasbdry(1)
            do j=1,nplasbdry
                rplasmin = MIN (rplasmin,rplasbdry(j))
                rplasmax = MAX (rplasmax,rplasbdry(j))
                zplasmin = MIN (zplasmin,zplasbdry(j))
                zplasmax = MAX (zplasmax,zplasbdry(j))
            end do
c
c --- set up the pointer and coefficient arrays if necessary:
c
             updownsym = .false.
             if (mhdmethd .eq. 'sorpicrd' .or. mhdmethd .eq. 'toq')then
                if (isym .eq. 1)  updownsym = .true.
                call mapvar(nw,nh,wzero,map,nplas,rplasbdry,
     .            zplasbdry,nplasbdry,rmhdgrid,zmhdgrid,coeff,updownsym)
             end if
c
c --- set up vector bdindex (used in determination of boundary condition
c --- on MHD grid below)
c
            jsymetric = 1
            if (updownsym .and. (mhdmethd .eq. 'sorpicrd'  .or.
     .          mhdmethd .eq. 'toq')) jsymetric = nh / 2 + 1
            call load_index (nw, nh, rmhdgrid, zmhdgrid,
     .                       rplasbdry, zplasbdry, nplasbdry,
     .                       bdindex, jsymetric)
      end if
c
c ----------------------------------------------------------------------
c decide on solution method, either cyc reduction or successive overrelaxation
c or inverse solver (TOQ):
c ----------------------------------------------------------------------
c
      if (mhdmethd .eq. 'sorpicrd' ) then    ! variable finite diff. grid
        !need to create eqdsk even if toq is run.
        ideriv   = 1
        navgbd   = 1
        residin  = 0.1*toleq
        residout =     toleq
        call sorpicard (nw,nh,map,coeff,rmhdgrid,psi1d,u0 ,
     .                  psiaxis,psilim,curden2d,wzero,psi,
     .                  rplasbdry,zplasbdry,nplasbdry,
     .                  cspln,n2cspln,nh2,pds,zmhdgrid,wnoperm,
     .                  tocur,psidifbd,zero,source,
     .                  psiprevious,psinocoil,ieq,
     .                  bdindex,ideriv,navgbd,psibdry,
     .                  itermaxp,ipassmax,nj,psir,
     .                  outeritused,xax,yax,psidifbdext,
     .                  iextiters,updownsym,optwi,optomegai,optwe,
     .                  optomegae,anormal,residin,residout,
     .                  pprim,ffprim,elongax,psincbcd,kside,omeq)
c
c       optimum relaxation parameters need be found only
c       on first call to SORPICARD. turn them off now:
c
        optwi = .false.
        optwe = .false.
      else if (mhdmethd .eq. 'cycred') then
        navgbd   = 1
        residin  = 0.1 * toleq
        residout =       toleq
        call cycsoltn (nw,nh,nwh,psi,p,psi1d,nj,pprim,ffprim,
     .                 psir,iounit,isgngren,ieq,
     .                 kside,rmhdgrid,zmhdgrid,wzero,zero,
     .                 omeq,curden2d,u0,xax,yax,xlimiter,
     .                 ylimiter,nlimiter,rcontr,zcontr,ncontr,
     .                 rplasbdry,zplasbdry,nplasbdry,cspln,
     .                 n2cspln,nh2,wnoperm,pcurrent,psincbcd,
     .                 boundrcntr,boundzcntr,nboundpts,psibdry,
     .                 toleq,iteq,outer_iterations,bdindex,
     .                 psidif,residout,pds,elongax,ideriv,
     .                 ieqfail,navgbd,nconmax,nwork,rplasmin,
     .                 anormal,tocur,psinctemp,psiaxis)
        else if(mhdmethd.eq.'toq') then
           call using_toq(ic)
           return
      else
        write  (nout, 510)  mhdmethd
        write  (ncrt, 510)  mhdmethd
  510   format (' ERROR: mhdmethd =', a, ' is not valid here')
        call STOP ('subroutine FIXDBDRY: problem #4', 105)

      endif



      if (mhdmethd .eq. 'toq_dumy')then  !skip this branch for now

         !get toq executable, write toq input file, spawn toq.
         !data passed to toq though eqdsk called fneqdsk.
         !on very first call fneqdsk is set to eqdskin in sub init(cray102.f).
         !On subsequqnt calls we do what ??
         call toq_drive(tocur)


         !read toq output file, dskonetwo which supplies flux surface
         !average quantities (need to work on this)
         !call read_toq_output(npsi_toq,psiaxis_toq,psibdry_toq)
         !in subsequent calls to toq use the latest toq guess file
         !ieqdtoq =0,ieqdsk = 2 in toq namelist  must be set
         !need to load psi, psi1d,p (in small.i,loaded below)

         if( .not. allocated(psinarray))then
              allocate (psinarray(npsi_toq),STAT = istat)
              if(istat .ne. 0)
     .          call allocate_error("psinarray",0,istat)
          endif
         if( .not. allocated(ffparray))then
              allocate (ffparray(npsi_toq),STAT = istat)
              if(istat .ne. 0)
     .          call allocate_error("ffparray",0,istat)
          endif
         if( .not. allocated(pparray))then
              allocate (pparray(npsi_toq),STAT = istat)
              if(istat .ne. 0)
     .          call allocate_error("pparray",0,istat)
          endif
         !get information for running fixed boundary code,
         !from toq output file dskfixb :
         call toqfile(rcfixb,raxtoq,zaxtoq,btorfixb,totcurtoq,
     &   psinarray,ffparray,pparray,npsitoq,rplasbdry,zplasbdry,
     &                                                nplasbdry)
        ideriv   = 1
        navgbd   = 1
        residin  = 0.1*toleq
        residout =     toleq

        call sorpicard (nw,nh,map,coeff,rmhdgrid,psi1d,u0 ,
     .                  psiaxis,psilim,curden2d,wzero,psi,
     .                  rplasbdry,zplasbdry,nplasbdry,
     .                  cspln,n2cspln,nh2,pds,zmhdgrid,wnoperm,
     .                  tocur,psidifbd,zero,source,
     .                  psiprevious,psinocoil,ieq,
     .                  bdindex,ideriv,navgbd,psibdry,
     .                  itermaxp,ipassmax,nj,psir,
     .                  outeritused,xax,yax,psidifbdext,
     .                  iextiters,updownsym,optwi,optomegai,optwe,
     .                  optomegae,anormal,residin,residout,
     .                  pprim,ffprim,elongax,psincbcd,kside,omeq)


c        print *,'raxtoq,xax(1),zaxtoq,yax(1)= ',
c     .           raxtoq,xax(1),zaxtoq,yax(1)
c
c       optimum relaxation parameters need be found only
c       on first call to SORPICARD. turn them off now:
c
        optwi = .false.
        optwe = .false.

      end if
c
      rma = xax(1)
      zma = yax(1)
c
c --- now scale the solution to the desired current
c
      write (6, *)  'anormal norm factor:', anormal
      do j=1,nh
        do i=1,nw
          k        = (i-1) * nh+j
          psi1d(k) = anormal * (psi1d(k)-psibdry)+psibdry
          psi(i,j) = psi1d(k)
          p(i,j)   = psi(i,j)  ! p (in small.i) is used in cray321.f
        end do
      end do
      psilow   = anormal * (psilow  - psibdry) + psibdry
      psihigh  = anormal * (psihigh - psibdry) + psibdry
      psiavg   = anormal * (psiavg  - psibdry) + psibdry
      psiaxis  = anormal * (psiaxis - psibdry) + psibdry
      pmin     = psiaxis     ! pmin is saved in small.i, used as an..
c                            ..integration constant in subroutine PSIRHO
      isignpsi = -1
      call psiscale (psiaxis, psibdry, psir, nj, isignpsi)
c
c --- check the new total current
c
      call curtotcalc (rmhdgrid, nw, zmhdgrid, nh, cspln, n2cspln, nh2,
     .                 psi, pds, wnoperm, rplasbdry, zplasbdry,
     .                 nplasbdry, u0, pcuriter, psilow, psiavg, psihigh)
      write (6, *)  'total current after normalization ', pcuriter
      write (6, *)  'psilow, psiavg, psihigh ', psilow, psiavg, psihigh
      write (6, *)                  'psibdry ', psibdry
      return
c
      end

      subroutine fixedcntour (rplasbdry,zplasbdry,nplasbdry,
     .                        rcontr,zcontr,ncontr,
     .                        rcmin,rcmax,zcmin,zcmax,
     .                        rzcmin,rzcmax,zrcmin,zrcmax,
     .                        rmhdgrid,zmhdgrid,nw,nh,
     .                        bpcontr,cspln,n2cspln,nh2,pds)
c
      USE replace_imsl,       ONLY : my_dbcevl1

      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- subroutine handles the fixed boundary plasma contour case
c ------------------------------------------------------------------ HSJ
c
      dimension rplasbdry(*),zplasbdry(*),rcontr(*),zcontr(*),
     .          pds(*),cspln(n2cspln,nw,nh2),rmhdgrid(nw),
     .          zmhdgrid(nh),bpcontr(*)
c
          rcmin = rplasbdry(1)
          rcmax = rplasbdry(1)
          zcmin = zplasbdry(1)
          zcmax = zplasbdry(1)
          do j=1,nplasbdry
              rcontr(j) = rplasbdry(j)
              zcontr(j) = zplasbdry(j)
              rcmin = MIN (rcmin,rcontr(j))
              if (rcmin .eq. rcontr(j))  zrcmin = zcontr(j)
              rcmax = MAX (rcmax,rcontr(j))
              if (rcmax .eq. rcontr(j))  zrcmax = zcontr(j)
              zcmin = MIN (zcmin,zcontr(j))
              if (zcmin .eq. zcontr(j))  rzcmin = rcontr(j)
              zcmax = MAX (zcmax,zcontr(j))
              if (zcmax .eq. zcontr(j))  rzcmax = rcontr(j)
              call my_dbcevl1 (rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                      rcontr(j),zcontr(j),pds,ier,3)
              bpcontr(j) = (1.0/rcontr(j)) * SQRT (pds(2)**2+pds(3)**2)
          end do
          ncontr = nplasbdry
c
      return
c
      end

      subroutine inner_solution (residin,maxinit,iterations,
     .                           updownsym,jsymetric,nw,nh,map,psi1d,
     .                           coeff,omega,psibdry,residmax,source)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- find the solution in the plasma interior.
c --- the inner iterations are standard SOR with relaxation parameter
c --- omega. odd/even (red/black) ordering is not done because the
c --- structure of the matrix is such that this ordering cannot be
c --- maintained. similarly ADI methods offer no advantage here.
c --- chebychev acceleration is not done because we are dealing
c --- mostly with the asymptotic solution where the optimal omega
c --- should be best.  Algebraic multigird is not used because the setup
c --- time is not acceptable once we have a good approximation to the
c --- solution.  We could use AMG for the first few outer iterations
c --- but this seems like too much trouble to be worth the effort.
c
c --- input
c
c  residin               maximum residual allowed. convergence
c                        will be assumed achieved if the largest
c                        residual is less than residin.
c  maxinit               maximum number of SOR iterations allowed
c                        for the right omega we expect the iterations
c                        to converge in about 2*(nw+nh) iterations
c                        (but also depends on value of residin)
c  omega                 relaxation parameter to use in SOR. only
c                        values between 1.0 and 2.0 make sense.
c  psibdry               value of psi on plasma boundary
c  updownsym             logical,true if up/down sym solution is to be found
c  nw,nh                 rectangular grid size (for this method of
c                        solution no restriction is placed on nw,nh.
c                        However if updown symetric solution is requested
c                        then it is assumed that nh is odd !!!!!!!!!
c  jsymetric             defined
c  map(k)                defined
c  coeff(9,*)            defined
c  source(k)             the rhs of the finite difference
c                        form of the GS equation. only values of k
c                        corresponding to the plasma interior will be used.
c
c --- output
c
c  iterations            the number of iterations required for convergence
c                       (If the residual does not drop below residmax
c                        in less than maxinit iterations then iterations=
c                        maxint will be returned)
c  psi1d(k)             solution for k corresponding to grid points
c                       inside the plasma.
c  residmax              maximum residual (if iterations is less than
c                        maxinit then residmax will be less than
c                        the requested value of residin.
c                        Otherwise convergence was not achieved in
c                        the allowed maxinit iterations. )
c ------------------------------------------------------------------ HSJ
c
      logical  updownsym
      integer  nw,nh,map(*),iterations,maxinit,jsymetric,
     .         i,j,k,l,ll,lt,lr,lb
      real*8   residmax,residin,psi1d(*),coeff(9,*),omega,
     .         residual,psibdry,psill,psilt,psilr,psilb,cl,source(*)
c
      iterations = 0
      residmax = residin + 1.0e10
      do while (residmax .gt. residin .and. iterations .lt. maxinit)
        iterations = iterations+1
        residmax = 0.0
        do j=jsymetric,nh
            do i=1,nw
                k = (i-1)*nh+j
                if (map(k) .ne. 0) then   ! do only if pt l is inside
c
                    ll = map(k-nh)        ! ll to left of l
                    if (ll .eq. 0) then
                      psill = psibdry     ! ll on left boundary
                    else
                      psill = psi1d(k-nh) ! ll left of l, inside plasma
                    end if
c
                    lt = map(k+1)        ! lt above l
                    if (lt .eq. 0) then
                      psilt = psibdry    ! lt on top boundary
                    else
                      psilt = psi1d(k+1) ! lt on top of l, inside plasma
                    end if
c
                    lr = map(k+nh)        ! lr to right of l
                    if (lr .eq. 0) then
                      psilr = psibdry     ! lr on right boundary
                    else
                      psilr = psi1d(k+nh) ! lr right of l, inside plasma
                    end if
c
                    if (updownsym .and. j .eq. jsymetric) then
c
c                         to force up/down symmetry the boundary condition
c                         at the mid plane (zmhdgrid(jsymetric) = 0.0 ) is
c                         dpsi/dz = 0.0.
c                         this is enforced by setting psilb = psilt:
c
                        psilb = psilt
                    else
                        lb = map(k-1)        ! lb below l
                        if (lb .eq. 0) then
                          psilb = psibdry    ! lb on bottom boundary
                        else
                          psilb = psi1d(k-1) ! lb below l, inside plasma
                        end if
                    end if
c
                    l = map(k)               ! l the current unknown
c
c                     solve the linear equation
c                       psil +cl-sl*coeff(5,l) = 0.0
c                     for the unknown psil.
c
                    cl = coeff(1,l)*psill+coeff(2,l)*psilt
     .                         +coeff(3,l)*psilr+coeff(4,l)*psilb
                    residual = psi1d(k)+cl-source(k)*coeff(5,l)
                    psi1d(k) = psi1d(k)-omega*residual
                     residmax = MAX (ABS (residual),residmax)
c
                end if          ! map(k) inside
            end do              ! end on row    index i
        end do                  ! end on column index j
      end do                    ! end do-while loop
c
      return
c
      end
      subroutine limiter_check(rcminm,rcmaxm,zcminm,zcmaxm,xlimiter,
     .                    ylimiter,nlimiter)
c-------------------------------------------------------------------
c     check to make sure that limiter remains consistent with
c     fixed boundary calcualtions when the boundary is allowed
c     to move slightly ( limiter may have been set tightly up
c     against the plasma,not allowing for any variation)
c     Return the limiter points modified to accomodate the new
c     extremes in the boundary. Dont allow limiter to ,move off
c     the computational grid however. (such an adjustment will require
c     a new starting eqdsk).
c     All quantities are in meters  
c     
c---------------------------------------------------------HSJ-1-20-01
      USE mhdpar,only :  nw,nh
      implicit none
      real *8 rcminm,rcmaxm,zcminm,zcmaxm
      real *8 ylimiter(*),xlimiter(*)
      real *8 rmhdgrid, zmhdgrid,xdim, ydim, ymid, redge
      integer *4 nlimiter
c      include 'mhdpar.i'  ! nw,nh
c      include 'mhdgrid.i' ! rmhdgrid,zmhdgrid
      
      xlimiter(nlimiter+1) = MIN(xlimiter(nlimiter+1),rcminm)
c      xlimiter(nlimiter+1) = MAX(xlimiter(nlimiter+1),rmhdgrid(1))
      xlimiter(nlimiter+2) = MAX(xlimiter(nlimiter+2),rcmaxm)
c      xlimiter(nlimiter+2) = MIN(xlimiter(nlimiter+2),rmhdgrid(nw))
      ylimiter(nlimiter+1) = MIN(ylimiter(nlimiter+1),zcminm)
c      ylimiter(nlimiter+1) = MAX(ylimiter(nlimiter+1),zmhdgrid(1))
      ylimiter(nlimiter+2) = MAX(ylimiter(nlimiter+2),zcmaxm)
c      ylimiter(nlimiter+2) = MIN(ylimiter(nlimiter+2),zmhdgrid(nh))
      return 
      end



 

      subroutine load_index (nw,nh,rmhdgrid,zmhdgrid,
     .                       rplasbdry,zplasbdry,nplasbdry,
     .                       bdindex,jsymetric)
c
      implicit none
c
c ----------------------------------------------------------------------
c  calculate the index of the point on the plasma boundary
c  which is closest to point j on the boundary of the MHD grid.
c  store the results in integer vector bdindex(j),j = 1,2,3...2(nw+nh)-4 .
c  bindex(i) = k means the point on the plasma boundary closest to
c  point i on the boundary of the mhdgrid is given by
c  rplasbdry(k),zplasbdry(k)
c  bindex is arranged to start with j = 1 in the lower left hand
c  corner of the grid and proceeds clockwise. hence we have
c
c    bindex(1) =                 (1,1)
c    bindex(nh) =                  (1,nh)
c    bindex(nh+nw-1) =           (nw,nh)
c    bindex(nh+nw-1+nh-1)      = (nw,1)
c    bindex(2*(nh+nw)-4)       = (2,1)
c
c  if the plasma is up/down symmetric then only the upper
c  half of the boundary is used. we then have
c
c    bindex(1:jsymetric-1)                        not used
c    bindex(jsymetric)                            (1,jsymetric)
c    bindex(nh)                                   (1,nh)
c    bindex(nh+nw-1) =                            (nw,nh)
c    bindex(nh+nw-1+jsymetric-1) =                (nw,jsymetric)
c    bindex(nh+nw-1+jsymetric: 2*(nh+nw)-4)       not used
c
c --- input:
c  nw             size of radial MHD grid
c  nh             size of vertical MHD grid
c  rmhdgrid(j)    j = 1,2,..nw radial MHD grid
c  zmhdgrid(j)    j = 1,2...nh vertical MHD grid
c  rplasbdry(j)   j = 1,2..nplasbdry  radial plasma boundary points
c  zplasbdry(j)   j = 1,2..nplasbdry  vertical plasam boundary points
c  nplasbdry      number of points on plasma boundary
c  jsymetric      set to 1 if plasma is not to be treated as up/down
c                 symetric
c                 set to nh/2+1  if plasma is to be treated as up/down
c                 symetric
c
c --- output
c  bdindex(i)     i = 1,2,3,...2*(nw+nh)-4 (see above for description )
c                          or
c  bdindex(i)     i = jsymetric,.....nh+nw+jsymetric-2
c
c ------------------------------------------------------------------ HSJ
c
      integer nw,nh,nplasbdry,bdindex(*),ib,kk,kk1,j,ir,il,l,i,
     .        jsymetric
      real*8 rmhdgrid(nw),zmhdgrid(nh),
     .       rplasbdry(nplasbdry),zplasbdry(nplasbdry),
     .       drr,drl,rrm,rlm,dzr,dzl,rl,rr
      ib = 2*(nw+nh) - 4                ! length of bdindex
c
      if (jsymetric .eq. 1) then        ! non-symmetric case
c
c --- the vertical boundaries(rmhdgrid(1),zmhdgrid(j),j = 1,2..nh)
c --- and rmhdgrid(nw),zmhdgrid(j),j = 1,2..nh)
c
            do j=1,nh
              kk = 2*nh+nw-1-j
              rrm = 1.0d+100
              rlm = 1.0d+100
              do l=1,nplasbdry
                drr = -rplasbdry(l)+rmhdgrid(1)
                drl = -rplasbdry(l)+rmhdgrid(nw)
                dzr = -zplasbdry(l)+zmhdgrid(j)
                dzl = -zplasbdry(l)+zmhdgrid(j)
                rr = SQRT (drr**2+dzr**2)
                rl = SQRT (drl**2+dzl**2)
                rrm = MIN (rrm,rr)
                rlm = MIN (rlm,rl)
                if (rrm .eq. rr)  ir = l
                if (rlm .eq. rl)  il = l
              end do
              bdindex(j) = ir
              bdindex(kk) = il
            end do
c
c --- horizontal boundaries ((rmhdgrid(i),zmhdgrid(1)),i = 2,..nw-1,
c --- and (rmhdgird(i),zmhdgrid(nh)),i = 2,...nw-1):
              kk1 = nh
              kk = ib+1
              do i=2,nw-1
                kk1 = kk1+1
                kk = kk-1
                rrm = 1.0d+100
                rlm = 1.0d+100
                do l=1,nplasbdry
                    drr = -rplasbdry(l)+rmhdgrid(i)
                    drl = -rplasbdry(l)+rmhdgrid(i)
                    dzr = -zplasbdry(l)+zmhdgrid(1)
                    dzl = -zplasbdry(l)+zmhdgrid(nh)
                    rr = SQRT (drr**2+dzr**2)
                    rl = SQRT (drl**2+dzl**2)
                    rrm = MIN (rrm,rr)
                    rlm = MIN (rlm,rl)
                    if (rrm .eq. rr)  ir = l
                    if (rlm .eq. rl)  il = l
                end do
                bdindex(kk) = ir
                bdindex(kk1) = il
              end do
      else    ! symmetric case
c
c --- the vertical boundaries(rmhdgrid(1),zmhdgrid(j),j = jsymetric..nh)
c --- and rmhdgrid(nw),zmhdgrid(j),j = jsymetric,..nh)
c
            do j=jsymetric,nh
              kk = nh+nw+2*jsymetric-2-j
              rrm = 1.0d+100
              rlm = 1.0d+100
              do l=1,nplasbdry
                   if (zplasbdry(l) .ge. 0.0) then
                      drr = -rplasbdry(l)+rmhdgrid(1)
                      drl = -rplasbdry(l)+rmhdgrid(nw)
                      dzr = -zplasbdry(l)+zmhdgrid(j)
                      dzl = -zplasbdry(l)+zmhdgrid(j)
                      rr = SQRT (drr**2+dzr**2)
                      rl = SQRT (drl**2+dzl**2)
                      rrm = MIN (rrm,rr)
                      rlm = MIN (rlm,rl)
                      if (rrm .eq. rr)  ir = l
                      if (rlm .eq. rl)  il = l
                   end if
              end do
              bdindex(j) = ir
              bdindex(kk) = il
            end do
c
c --- horizontal boundary
c --- (rmhdgird(i),zmhdgrid(nh)),i = 2,...nw-1):
              kk1 = nh
              kk = ib+1
              do i=2,nw-1
                kk1 = kk1+1
                kk = kk-1
                rrm = 1.0d+100
                rlm = 1.0d+100
                do l=1,nplasbdry
                   if (zplasbdry(l) .ge. 0.0) then
                          drl = -rplasbdry(l)+rmhdgrid(i)
                          dzl = -zplasbdry(l)+zmhdgrid(nh)
                          rl = SQRT (drl**2+dzl**2)
                          rlm = MIN (rlm,rl)
                          if (rlm .eq. rl)  il = l
                   end if
                end do
                bdindex(kk1) = il
              end do
      end if
c
      return
c
      end

      subroutine magax2 (psi,x,y,nw,nh, isignn, iknowax,iounit,ispln,
     .                   zero,xax,yax,cspln,n2cspln,nh2,
     .                   psymx,elax,map)

c
c ----------------------------------------------------------------------
c --- find the magnetic axis. by definition the magnetic axis is that
c --- point inside the limiter contour at which psi takes on an extremal
c --- value,either a maximum (positive current) or a minimum (negative
c --- current). convergence to the magnetic axis is considered accomplished
c --- when gradient of psi is .lt. gradtol (gradtol is defined locally below)
c --- subroutine terminates ONETWO run if magnetic axis can not be found!
c --- this routine is similar to subroutine MAGAX but limits the search to
c --- a local region inside the plasma. this is done so that psi
c --- values external to the plasma can't influence the results.
c --- input
c  psi(nw,nh)         the array of psi values
c  x(nw)
c  y(nh)              the MHD grid vectors
c  isignn        =  1 search for maximum in psi
c                = -1 search for minimum in psi
c  iknowax         =1 good intial guess for magnetic axis is input
c                     in xax(1),yax(1)
c                  =0 good initial guess is not available, find initial
c                     guess for magnetic axis by searching the MHD grid
c  iounit             Fortran unit number for diagnostic output
c                     iounit = 0 suppresses output, but this routine
c                     may terminate the ONETWO run so diagnostic
c                     messages should be printed out.
c  ispln             =1 bicubic spline array cspln is input
c                    =0 bicubic spline array cspln must be calculated
c  zero(nwh)         the limiter contour indicator vector
c  xax
c  yax       intial guess of magnetic axis location if known, see iknowax
c  cspln(n2cspln,nw,nh2)   bicubic spline coeff. of psi if known
c                          if not known,set ispln = 0 to calculate it here
c  wk(ndwk)      temporary starge vector of minimum length
c                ndwk = 2*nh*nw+2*max(nh,nw) ( required if ispln=0 )
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

      USE replace_imsl,                 ONLY : my_ibcccu,my_dbcevl1

      implicit  integer (i-n), real*8 (a-h, o-z)
      dimension psi(nw,nh),x(nw),y(nh)
      dimension xax(*),yax(*),psymx(*),zero(*)
      dimension cspln(n2cspln,nw,nh2),pds(6)
      dimension xx(5),yy(5),psiloc(5,5),map(*),csloc(2,5,10)
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
      iknowax = 0
      i1 = 0
      j1 = 0
      if (isignn .eq. -1) then     ! find minimum psi
        psimin = +1.0d100
        do 10 j=2,nh-1
        do 10 i=2,nw-1
        kk = (i-1)*nh+j
        zzsum = zero(kk-nh-1)+zero(kk-nh)+zero(kk-nh+1)+zero(kk-1)+
     .  zero(kk)+zero(kk+1)+zero(kk-1+nh)+zero(kk+nh)+zero(kk+nh+1)
        if (zzsum .lt. 8.5)  go to 10
        if (psi(i,j) .gt. psi(i-1,j-1))  go to 10
        if (psi(i,j) .gt. psi(i-1,j))  go to 10
        if (psi(i,j) .gt. psi(i-1,j+1))  go to 10
        if (psi(i,j) .gt. psi(i,j-1))  go to 10
        if (psi(i,j) .gt. psi(i,j+1))  go to 10
        if (psi(i,j) .gt. psi(i+1,j-1))  go to 10
        if (psi(i,j) .gt. psi(i+1,j))  go to 10
        if (psi(i,j) .gt. psi(i+1,j+1))  go to 10
c
c --- the above picks smallest value of psi on the x,y grid, inside the plasma,
c --- such that the eight nearest neighbors of pt. kk are also inside the plasma
c --- thus a magnetic axis within dx or dy of limiter contour as seen by array
c --- zero will not be found using this method.
c
        psimin = MIN (psimin,psi(i,j))
        if (psimin .ne. psi(i,j))  go to 10
        i1 = i
        j1 = j
   10   continue
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
c --- the above picks largest value of psi on the x,y grid, inside the plasma,
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
  125         format (' subroutine MAGAX2 reports:' /
     .                ' magnetic axis was not found')
              iknowax = 0
              return
          end if
      else
        if (iounit .ne. 0)  write (iounit, 25) isignn
   25   format (' subroutine MAGAX2 reports'         /
     .          ' incorrect setting of ISIGNN =', i5 /
     .          ' ONETWO must stop'                  /)
        call STOP ('subroutine MAGAX2: problem #1', 37)
      end if
      xax(1) = x(i1)
      yax(1) = y(j1)    ! estimate of magnetic axis on the grid
c
c --- set up grid interior to plasma,centered on approximate
c --- magnetic axis, nx = ny=4 is minimum allowed by ibcccu:
c
          xx(1) = x(i1-2)
          xx(2) = x(i1-1)
          xx(3) = x(i1)
          xx(4) = x(i1+1)
          xx(5) = x(i1+2)
          yy(1) = y(j1-2)
          yy(2) = y(j1-1)
          yy(3) = y(j1)
          yy(4) = y(j1+1)
          yy(5) = y(j1+2)
          do i=1,5
            do j=1,5
                psiloc(i,j) = psi(i1-3+i,j1-3+j)
              k = (i1-4+i)*nh+j1-3+j
              if (map(k) .eq. 0) then
                write (6, *)  'MAGAX2 reports grid is outside plasma'
                call STOP ('subroutine MAGAX: problem #2', 107)
              end if
            end do
          end do
c
c --- fit spline over local grid:
c
      imslmd ='1622c202'
      call my_ibcccu (psiloc,xx,5,yy,5,csloc,5,wk,ier)
      deallocate(wk)
      ispln = 0     ! because we changed definition of cspln
c
c --- final convergence by Newton's method. finds simultaneous
c --- solution of dpsi/dx = 0.0,dpsi/dy = 0.0
c --- maximum of 25 iterations are allowed.
c
      
  130 psymx(1) = 0.0
      icalc = 6
      do 100 j=1,25
      jax = j
      if (xax(1) .lt. xx(1) .or. xax(1) .gt. xx(5))  go to 115
      if (yax(1) .lt. yy(1) .or. yax(1) .gt. yy(5))  go to 115
      xaxold = xax(1)
      yaxold = yax(1)
      call my_dbcevl1 (xx,5,yy,5,csloc,5,xax(1),yax(1),pds,ier,6)
      pds1 = pds(1)
      det = pds(5)*pds(6)-pds(4)*pds(4)
      if (det .eq. 0.0)  go to 110           ! can take no further steps
      xerr = (-pds(2)*pds(6)+pds(4)*pds(3))/det
      yerr = (-pds(5)*pds(3)+pds(2)*pds(4))/det
   21 xnew = xaxold+xerr
      if (xnew .gt. xx(5) .or. xnew .lt. xx(1)) then
          xerr = xerr*0.5
          if (ABS (xerr) .gt.  1.0e-14)  go to 21
      end if
   22 ynew = yaxold+yerr
      if (ynew .gt. yy(5) .or. ynew .lt. yy(1)) then
          yerr = yerr*0.5
          if (ABS (yerr) .gt. 1.0e-14)  go to 22
      end if
c
c --- limit the Newton step to ensure convergence
c --- by cutting back xerr,yerr until a smaller (isignn = -1)
c --- or larger (isignn = +1) psi is obtained.
c
   23 xnew = xaxold+xerr
      ynew = yaxold+yerr
      call my_dbcevl1 (xx,5,yy,5,csloc,5,xnew,ynew,pds,ier,6)
      if (isignn .eq. -1) then
          if (pds(1) .gt. pds1) then
               xerr = xerr*0.5
               yerr = yerr*0.5
               if (ABS (xerr) .gt. 1.0e-14 .or. ABS (yerr) .gt. 1.0e-14)
     .                        go to 23
          end if
      else
          if (pds(1) .lt. pds1) then
              xerr = xerr*0.5
              yerr = yerr*0.5
               if (ABS (xerr) .gt. 1.0e-14 .or. ABS (yerr) .gt. 1.0e-14)
     .                        go to 23
          end if
      end if
      xax(1) = xax(1)+xerr
      yax(1) = yax(1)+yerr
      if (pds(2)**2+pds(3)**2 .lt. gradtol)  go to 110
****  if (ABS (xax(1)-xaxold) + ABS (yax(1)-yaxold) .lt. 1.0e-10)  go to 110
      if (ABS (xerr) + ABS (yerr) .lt. 1.0e-6)  go to 110
  100 continue
c --- failure to find a magnetic axis is a fatal error
  115 continue
      if (iounit .ne. 0)
     .write (iounit,120)jax,itry,xax(1),xerr,yax(1),yerr,det
  120 format (' subroutine MAGAX2 reports no convergence to'   /
     .        ' magnetic axis after ',i5,' Newton iterations'  /
     .        ' and ',i5,' outer passes'                       //
     .        ' current estimate is xax,xerr = ',2(2x,1pe12.4) /
     .        '                     yax,yerr = ',2(2x,1pe12.4) /
     .        ' the discriminant is ',1pe12.4                  //
     .        ' (negative discriminant means xax,yax',
     .        ' is tending toward a saddle point)')
      iknowax = 0
      return
  110 if (det .lt. 0.0) then  ! det < 0  means converged to saddle point
         xax(1) = x(i1)
          yax(1) = y(j1)
          pds(1) = psimax
          if (isignn .eq. -1)  pds(1) = psimin
      end if
      iknowax  = 1            ! converged to magnetic axis
      psymx(1) = pds(1)
c
c --- calculate elongation on axis based on general quadratic
c --- form of an ellipse;a*x**2+2b*x*y+c*y**2+d*x+e*y+f
c --- the coefficients a,b,c,d,e,f,are obtained from the taylor
c --- series of psi about the magnetic axis
c
      a = 0.5*pds(5)
      b = 0.5*pds(4)
      c = 0.5*pds(6)
c
c --- rotate coordinates
c
      thet = 0.5 * ATAN2 (2*b, a-c)
      ap   = a * COS (thet)**2 + 2.0 * b * SIN (thet) * COS (thet)
     .     + c * SIN (thet)**2
      cp   = a * SIN (thet)**2 - 2.0 * b * SIN (thet) * COS (thet)
     .     + c * COS (thet)**2
c      elax = SQRT (ap/cp)
      if (cp  .eq. 0.0) then
          write (6, *)  'cp=0.0 in MAGAX2'
          write (6, *)  'xax ,yax =',xax(1),yax(1)
          write (6, *)  'psiaxis =',psymx(1)
          write (6, *)  'thet,a,b,c=',thet,a,b,c
          write (6, *)  'dsqsdxsq,dsqsdysq=',pds(5),pds(6)
          elax = 1.0
      else if (ap/cp .lt. 0.0) then
          write (6, *)  'ap =',ap
          write (6, *)  'cp=',cp
          write (6, *)  'xax ,yax =',xax(1),yax(1)
          write (6, *)  'psiaxis =',psymx(1)
          write (6, *)  'thet,a,b,c=',thet,a,b,c
          write (6, *)  'dsqsdxsq,dsqsdysq=',pds(5),pds(6)
          elax = 1.0
      else
          elax = SQRT (ap/cp)
      end if
      return
c
      end

      subroutine magax3 (psi, x, y, nw, nh, isignn, iknowax, iounit,
     .                   xax, yax, c, psymx, elax, map)
c

      USE replace_imsl,                ONLY : my_icsicu,my_icsevu
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- find the magnetic axis. by definition the magnetic axis is that
c --- point inside the limiter contour at which psi takes on an extremal
c --- value, either a maximum or a minimum
c --- convergence to the magnetic axis is considered accomplished
c --- when gradient of psi < gradtol (gradtol is defined locally below)
c --- subroutine terminates ONETWO run if magnetic axis can not be found!
c --- this routine is accomplishes the same task as MAGAX and MAGAX2 but
c --- is limited to a line search along the up/down symmetry axis
c --- so the search is one-dimensional rather than two-dimensional.
c --- input
c  psi(nw,nh)         the array of psi values
c  x(nw)
c  y(nh)              the MHD grid vectors
c  isignn        =  1 search for maximum in psi
c                = -1 search for minimum in psi
c  iknowax       =  1 good initial guess for magnetic axis is input
c                     in xax(1), yax(1)
c                =  0 good initial guess is not available, find initial
c                     guess for magnetic axis by searching the MHD grid
c  iounit             Fortran unit number for diagnostic output
c                     iounit = 0 suppresses output.
c  xax
c  yax       intial guess of magnetic axis location if known,see iknowax
c  c(nw,3)            work array (not usable outside this subroutine
c                     because of dimension statements).
c
c  map(k)             k = 1,2,..nw*nh, used as indicator vector to
c                     search only grid points inside the plasma.
c                     map(k) is set up in subroutine MAPVAR.
c
c --- output
c  iknowax = 1       indicates axis is now known
c  iknowax = 0       if the axis was not found iknowax = 0 is returned
c  psymx           the value of psi on the magnetic axis
c  elax            the elongation of the plasma at the magnetic axis
c  xax
c  yax             final converged coordiantes of magnetic axis
c
c ------------------------------------------------------------------ HSJ
c
      dimension  psi(nw,nh), x(nw), y(nh)
      dimension  xax(*), yax(*), psymx(*), map(*)
      dimension  c(nw,3), bpar(4)
      logical    icoarse
c
      icoarse = .false.    ! coarse grid search not yet done
      gradtol = 1.0e-12
      steptol = 1.0e-08
      jmid    = nh/2+1
      dx      = x(2)-x(1)
      dxx     = dx
      dysq    = (y(2)-y(1))**2
      itry    = 1
      if (iknowax .eq. 1) then
        is    = (xax(1)-x(1))/dx+1
        isave = is
        go to 400
      end if
c
c --- rough search
c
  300 if (isignn .eq. 1) then    ! search for maximum in psi
          psimax = -1.0d100
          do i=1,nw
              k = (i-1)*nh+jmid
              if (map(k) .gt. 0) then
                  psimax = MAX (psi(i,jmid),psimax)
                  if (psimax .eq. psi(i,jmid))  is = i
              end if
          end do
      else                      ! search for minimum in psi
          psimin = 1.0d+100
          do i=1,nw
              k = (i-1)*nh+jmid
              if (map(k) .gt. 0) then
                  psimin = MIN (psi(i,jmid),psimin)
                  if (psimin .eq. psi(i,jmid))  is = i
              end if
          end do
      end if
      icoarse = .true.     ! coarse grid search is done
      isave   = is
      xax(1)  = x(is)
      yax(1)  = y(jmid)    ! estimate of magnetic axis on the grid
c
c --- final convergence, use cubic spline to get derivatives off the grid:
c --- subroutine icsicu (really, icsicu1) gets the coefficient array c.
c
c       the value of the spline approximation at  xx is
c             s = ((c(j,3)*d+c(j,2))*d+c(j,1))*d+y(j)
c        where x(j) .le. xx .lt. x(j+1) and
c              d = xx-x(j).
c       we also have ds/dx at xx:
c             ds/dx = (3*c(j,3)*d+2.0*c(j,2))*d+c(j,1)
c       and dsq(s)/dxsq at xx:
c             dsq(s)/dxsq = 6.0 * c(j,3)*d+2.0*c(j,2)
c
  400 do i=1,4
        bpar(i) = 0.0    ! selects natural cubic spline
      end do
      ier = 0
c      call icsicu1 (x, psi(1,jmid), nw, bpar, c, nw, ier)
      call my_icsicu (x, psi(1,jmid), nw, bpar, c, nw, ier)
      if (ier .ne. 0) then
          if (iounit .ne. 0)  write (iounit, 410)
  410     format (' subroutine MAGAX3 detects ERROR' /
     .            ' ier from ICSICU1 =', i5)
          iknowax = 0
          return
      end if
  500 il     = is - 2
      il     = MAX0 (il, 1)
      ip     = is + 2
      ip     = MIN0 (ip, nw-1)
      deltax = 0.0
      xx     = xax(1)
  450 do j=1,20       ! allow 20 Newton steps to find xx where dpsi/dx=0
          do i=il,ip  ! dpsi/dy  = 0.0 is enforced by boundary condition
c                       on the solution psi
              if (xx .ge. x(i) .and. xx .lt. x(i+1)) then
                  d        = xx-x(i)
                  dsdx     = (3.0*c(i,3)*d+2.0*c(i,2))*d+c(i,1)
                  dsqsdxsq = 6.0 * c(i,3)*d+2.0*c(i,2)
                  if (dsqsdxsq .eq. 0.0)  go to 1000
                  deltax   = -dsdx/dsqsdxsq
                  go to 50
              end if
          end do
              go to 200    ! drop through the loop, error exit
   50         xx = xx+deltax
              if (xx .gt. x(nw) .or. xx .lt. x(1))  go to 200
              if (ABS (dsdx) .lt. gradtol)  go to 100
              if (ABS (deltax) .lt. steptol)  go to 100
      end do
      go to 200    ! dropped through the Newton loop, error exit
c
c --- Newton converged.
c --- set axis and psi value at axis:
c
  100 xax(1)   = xx
      yax(1)   = y(jmid)
      psymx(1) = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+psi(i,jmid)
      iknowax  = 1                ! converged to magnetic axis
c
c --- calculate elongation on axis based on general quadratic
c --- form of an ellipse;a*x**2+2b*x*y+c*y**2+d*x+e*y+f
c --- the coefficients a,b,c,d,e,f,are obtained from the taylor
c --- series of psi about the magnetic axis
c
      a = 0.5*dsqsdxsq
c
c --- mixed derivative, dsqpsi/dxdy, = 0 for up/down symmetric case
c
      b = 0.0
c
c --- fit psi on the line one grid space up in z:
c
c      call icsicu1 (x, psi(1,jmid+1), nw, bpar, c, nw, ier)
      call my_icsicu (x, psi(1,jmid+1), nw, bpar, c, nw, ier)
      if (ier .ne. 0) then
        write (6, *)  'ier from ICSICU1 in MAGAX3 =', ier
      end if
c
c --- get the value of psi,psiz, at xax(1),y(jmid+1):
c
      i        = (xax(1)-x(1))/dx+1
      d        = xax(1)-x(i)
      psiz     = ((c(i,3)*d+c(i,2))*d+c(i,1))*d+psi(i,jmid+1)
      dsqsdysq = 2.0 * (psiz-psymx(1))/dysq
      cc       = 0.5 * dsqsdysq
c
c --- rotate coordinates
c
      thet = 0.5 * ATAN2 (2*b, a-cc)
      ap   =  a * COS (thet)**2 + 2.0 * b * SIN (thet) * COS (thet)
     .     + cc * SIN (thet)**2
      cp   =  a * SIN (thet)**2 - 2.0 * b * SIN (thet) * COS (thet)
     .     + cc * COS (thet)**2
      if (cp  .eq. 0.0) then
          write (6, *)  'cp = 0.0 in MAGAX3'
          write (6, *)  'xax ,yax =',xax(1),yax(1)
          write (6, *)  'psiaxis =',psymx(1)
          write (6, *)  'thet,a,b,cc=',thet,a,b,cc
          write (6, *)  'dsqsdxsq,dsqsdysq=',dsqsdxsq,dsqsdysq
          elax = 1.0
      else if (ap/cp .lt. 0.0) then
          write (6, *)  'ap =',ap
          write (6, *)  'cp=',cp
          write (6, *)  'xax ,yax =',xax(1),yax(1)
          write (6, *)  'psiaxis =',psymx(1)
          write (6, *)  'thet,a,b,cc=',thet,a,b,cc
          write (6, *)  'dsqsdxsq,dsqsdysq=',dsqsdxsq,dsqsdysq
          elax = 1.0
      else
          elax = SQRT (ap/cp)
      end if
      return
c
c --- no convergence to magnetic axis. get a refined starting guess
c --- and try one more time.
c --- failure to find a magnetic axis is a fatal error
c
  200 if (itry .gt. 3) then
          if (iounit .ne. 0)  write (iounit, 120)  xax(1), yax(1)
  120     format (' subroutine MAGAX3 reports no convergence to' /
     .            ' magnetic axis,xax,yax =',2(2x,1pe12.4))
          iknowax = 0
          return
      end if
      itry = itry + 1
      iknowax = 0
c
c --- Newton method failed to converge. If coarse grid search has
c --- not been done (because iknowax was set on entry) then go back
c --- and do it now:
c
      if (.not. icoarse)  go to 300
c
c --- Newton method failed to converge.  coarse grid search has been done.
c --- try a fine grid search:
c
      ism    = isave
      ism1   = ism-1
      isp1   = ism+1
      km     = ism1*nh+jmid
      km1    = (ism1-1)*nh+jmid
      kp1    = ism*nh+jmid
      psimax = -1.0d+100
      psimin =  1.0d+100
      if (map(km)*map(km1)*map(kp1) .eq. 0) then
c
c --- error,extrema in psi is too cloase to plasma boundary
c --- (or outside plasma).
c
          iknowax = 0
          if (iounit .gt. 0)  write (iounit, 1500)
 1500     format (' subroutine MAGAX3 reports that magnetic' /
     .            ' axis could not be found inside plasma'   /
     .            ' km1,km,kp1 =',3(2x,i5)                   /
     .            ' map(km1),map(km),map(kp1) =',3(2x,i5))
          return
      end if
c
      dxf = 0.05*(x(2)-x(1))    ! fine grid spacing
      xx  = x(ism1)
      do while (xx .le. x(isp1))
          xx = xx+dxf
          k  = ism
          if (xx .lt. x(ism))  k = ism1
          dxk    = xx-x(k)
          psival = ((c(k,3)*dxk+c(k,2))*dxk+c(k,1))*dxk+psi(k,jmid)
          if (isignn .eq. 1) then
            psimax = MAX (psimax,psival)
            if (psimax .eq. psival)  xsave = xx
          else    ! search for minimum in psi
            psimin = MIN (psimin,psival)
            if (psimin .eq. psival)  xsave = xx
          end if
      end do
      if (xsave .lt. x(ism)) then
          il = ism1
      else
          il = ism
      end if
      ip = il+1
      xx = xsave
      go to 450
c
 1000 if (iounit .ne. 0)  write (iounit, 1010)  xax(1), y(jmid)
 1010 format (' ERROR in subroutine MAGAX3: converged to saddle point' /
     .        ' xax,yax = ', 2(2x, 1pe12.4))
      iknowax = 0
      return
c
      end

      subroutine outer_solution (residin,maxinit,iterations,
     .                     updownsym,jsymetric,nw,nh,map,psi1d,coeff,
     .                     deltar,deltarzsq,deltarsq,interior,omega,pds,
     .                     rmhdgrid,psibdry,residmax)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- subroutine uses SOR to iterate the solution for psi in the
c --- region between the plasma and the edge of the rectangular
c --- grid to convergence.
c
c --- input
c
c  residin               maximum residual allowed. convergence
c                        will be assumed achieved if the largest
c                        residual is less than residin.
c  maxinit               maximum number of SOR iterations allowed
c                        for the right omega we expect the iterations
c                        to converge in about 2*(nw+nh) iterations
c                        (but also depends on value of residin)
c  psi1d(k)              for k corresponding to grid points on the
c                        boundary  of the rectangular grid,psi1d(k)
c                        is input (in volt*sec/rad)
c                        For these values of k psi1d remains unchanged
c                        on output.
c  omega                 relaxation parameter to use in SOR. only
c                        values between 1.0 and 2.0 make sense.
c  pds(6)                used  as work space.
c  psibdry               value of psi on plasma boundary
c  rmhdgrid(i)           i = 1,2..nw the major radius values
c  updownsym             logical,ture if up/down sym solution is to be found
c  nw,nh                 rectangular grid size (for this method of
c                        solution no restriction is placed on nw,nh.
c                        However if updown symetric solution is requested
c                        then it is assumed that nh is odd !!!!!!!!!
c  jsymetric             defined
c  map(k)                defined
c  coeff(9,*)            defined
c  deltar                defined
c  deltarzsq             defined
c  deltarsq              defined
c  interior         logical,set to false for this routine
c
c --- output
c
c  iterations            the number of iterations required for convergence
c                       (If the residual does not drop below residmax
c                        in less than maxinit iterations then iterations=
c                        maxint will be returned)
c  psi1d(k)             solution for k corresponding to grid points
c                       outside the plasma and not on the boundary
c                       of the rectangular grid
c  residmax              maximum residual (if iterations is less than
c                        maxinit then residmax will be less than
c                        the requested value of residin.
c                        Otherwise convergence was not achieved in
c                        the allowed maxinit iterations.)
c ------------------------------------------------------------------ HSJ
c
      logical  interior,updownsym
      integer  nw,nh,map(*),iterations,maxinit,jsymetric,
     .         i,j,k,kk,l,jstart
      real*8   residmax,residin,psi1d(*),pds(*),coeff(9,*),omega,
     .         deltar,deltarzsq,deltarsq,psill,all,psilt,alt,
     .         psilr,alr,psilb,alb,rl,residual,rmhdgrid(*),psibdry
c
      iterations = 0
      residmax = residin+1.0e10
      do while (residmax .gt. residin .and. iterations .lt. maxinit)
        iterations = iterations+1
        residmax = 0.0
        jstart = 2
        if (updownsym)  jstart = jsymetric
        do j=jstart,nh-1
            do i=2,nw-1
                k = (i-1)*nh+j
                if (map(k) .eq. 0) then    ! do only for exterior points
                          kk = k - nh
c
c                         pt to left is outside
c
                          if (map(kk) .eq. 0) then
                               psill = psi1d(kk)
                               all = 1.0
                          else
                              l = map(kk)    ! l is index of interior pt
                              all = 1.0-coeff(8,l)
                              psill = psibdry
                          end if
                          kk = k+nh
c
c                         pt to right is outside
c
                          if (map(kk) .eq. 0) then
                               psilr = psi1d(kk)
                               alr = 1.0
                          else
                              l = map(kk)    ! l is index of interior pt
                              alr = 1.0-coeff(6,l)
                              psilr = psibdry
                          end if
                          kk = k+1
c
c                         pt on top is outside
c
                          if (map(kk) .eq. 0) then
                               psilt = psi1d(kk)
                               alt   = 1.0
                          else
                              l = map(kk)    ! l is index of interior pt
                              alt   = 1.0-coeff(9,l)
                              psilt = psibdry
                          end if
                          kk = k-1
                          if (.not.updownsym .or. j .gt. jsymetric) then
c
c                             pt on bot is outside
c
                              if (map(kk) .eq. 0) then
                                   psilb = psi1d(kk)
                                   alb   = 1.0
                              else
                                l = map(kk)  ! l is index of interior pt
                                alb   = 1.0-coeff(7,l)
                                psilb = psibdry
                              end if
                          else          ! up/down sym with j = jsymetric
                              alb   = alt
                              psilb = psilt
                          end if
                          rl = rmhdgrid(i)
                      call coefcalc(pds,all,alt,alr,alb,
     .                            deltar,deltarzsq,deltarsq,rl,interior)
                        residual = psi1d(k)+psill*pds(1)+psilt*pds(2)
     .                                    +psilr*pds(3)+psilb*pds(4)
                        psi1d(k) = psi1d(k)-omega*residual
                        residmax = MAX (residmax, ABS (residual))
                end if           ! end point outside plasma branch
            end do               ! end i (row   ) loop
        end do                   ! end j (column) loop
      end do                     ! end do-while loop
c
      return
c
      end

      subroutine psibd (psi, rmhdgrid, zmhdgrid, rplasbdry, zplasbdry,
     .                  nplasbdry, psibdry, psi1d, bdindex, ideriv,
     .                  psinocoil, navgbd, jsymetric)
c
      USE mhdpar
      USE bicube
      USE replace_imsl,                     ONLY : my_ibcccu,my_dbcevl1
      implicit  integer (i-n), real*8 (a-h, o-z)
c
       include 'imsl.i'
c      include 'mhdpar.i'
c      include 'bicube.i'
c
      integer    bdindex(*)
      dimension  psi(nw,nh), rmhdgrid(nw), zmhdgrid(nh), psi1d(*)
      dimension  rplasbdry(*), zplasbdry(*), psinocoil(*)

c
      ibd = 2 * (nw + nh) - 4
      nhh = nh
      if (jsymetric .ne. 1)  nhh = jsymetric
      imslmd ='2194c202'
      call my_ibcccu (psi(1,jsymetric), rmhdgrid, nw,
     .             zmhdgrid(jsymetric), nhh, cspln, nw, wnoperm, ier)
      icalc = 6
      if (ideriv .eq. 1)
     .icalc = 3
      kk1   = 2 * nh + nw - 1
      if (jsymetric .ne. 1)  kk1 = nh + nw + jsymetric - 1
      do j=jsymetric,nh
        i = bdindex(j)    ! left vertical boundary
        call my_dbcevl1 (rmhdgrid, nw, zmhdgrid(jsymetric), nhh, 
     .                   cspln, nw,rplasbdry(i), zplasbdry(i), pds,
     .                   ier, 6)


        dr = rmhdgrid(1) - rplasbdry(i)
        dz = zmhdgrid(j) - zplasbdry(i)
        if (ideriv .eq. 2) then
            psi1d(j) = psibdry+(pds(2)+0.5*pds(5)*dr)*dr
     .               + (pds(3)+0.5*pds(6)*dz)*dz
     .               +  pds(4)*dz*dr
        else
            psi1d(j) = psibdry+pds(2)*dr
     .               + pds(3)*dz
        end if
c
c --- average around a number of plasma boundary points:
c
          if (navgbd .gt. 0) then
              ntotl  =  2 * navgbd + 1
              isignn = -1
   10         do jl=1,navgbd
                 k = i + jl * isignn
                 if (k .lt. 1        )  k =     nplasbdry - k
                 if (k .gt. nplasbdry)  k = k - nplasbdry
            call my_dbcevl1(rmhdgrid,nw,zmhdgrid(jsymetric),nhh,cspln,
     .           nw,rplasbdry(k),zplasbdry(k),pds,ier,6)
                 dr = rmhdgrid(1) - rplasbdry(k)
                 dz = zmhdgrid(j) - zplasbdry(k)
                 if (ideriv .eq. 2) then
                      psiav = psibdry+(pds(2)+0.5*pds(5)*dr)*dr
     .                         +(pds(3)+0.5*pds(6)*dz)*dz
     .                         +pds(4)*dz*dr
                  else
                      psiav = psibdry + pds(2)*dr + pds(3)*dz
                  end if
                  psi1d(j) = psi1d(j)+psiav
              end do
              isignn   = -isignn
              if (isignn .eq. 1)  go to 10
              psi1d(j) = psi1d(j)/ntotl
          end if
c
c --- end of average
c
        psinocoil(j) = psi1d(j)
c
        kk1 = kk1 - 1
        i   = bdindex(kk1)    ! right vertical boundary
        call dbcevl1 (rmhdgrid,nw,zmhdgrid(jsymetric),nhh,cspln,nw,
     .                rplasbdry(i),zplasbdry(i),pds,ier,6)
        dr  = rmhdgrid(nw) - rplasbdry(i)
        dz  = zmhdgrid(j ) - zplasbdry(i)
        kk  = (nw - 1) * nh + j
        if (ideriv .eq. 2) then
          psi1d(kk) = psibdry + (pds(2)+0.5*pds(5)*dr)*dr
     .              + (pds(3) + 0.5*pds(6)*dz) * dz + pds(4)*dz*dr
        else
          psi1d(kk) = psibdry + pds(2)*dr + pds(3)*dz
        end if
c
c --- average around a number of plasma boundary points:
c
          if (navgbd .gt. 0) then
              ntotl  =  2 * navgbd + 1
              isignn = -1
   20         do jl=1,navgbd
                 k = i + jl * isignn
                 if (k .lt. 1        )  k =     nplasbdry - k
                 if (k .gt. nplasbdry)  k = k - nplasbdry
                 call my_dbcevl1(rmhdgrid,nw,zmhdgrid(jsymetric),
     .                           nhh,cspln,
     .                       nw,rplasbdry(k),zplasbdry(k),pds,ier,6)
                 dr = rmhdgrid(nw) - rplasbdry(k)
                 dz = zmhdgrid(j ) - zplasbdry(k)
                 if (ideriv .eq. 2) then
                    psiav = psibdry+(pds(2)+0.5*pds(5)*dr)*dr
     .                    + (pds(3)+0.5*pds(6)*dz)*dz + pds(4)*dz*dr
                  else
                    psiav = psibdry + pds(2)*dr + pds(3)*dz
                  end if
                  psi1d(kk) = psi1d(kk) + psiav
              end do
              isignn    = -isignn
              if (isignn .eq. 1)  go to 20
              psi1d(kk) = psi1d(kk)/ntotl
          end if
c
c --- end of average
c
        psinocoil(kk1) = psi1d(kk)
      end do
c
c --- horizontal boundaries ((rmhdgrid(i),zmhdgrid(1)),i = 2,..nw-1,
c --- and (rmhdgrid(i),zmhdgrid(nh)),i = 2,...nw-1):
c
        kkk   = nh             ! top    boundary
        kkkk  = 2*(nh+nw) - 3  ! bottom boundary
        do i=2,nw-1
          kk  = i*nh           ! top    boundary
          kkk = kkk + 1        ! top    boundary
          j   = bdindex(kkk)
          call my_dbcevl1(rmhdgrid,nw,zmhdgrid(jsymetric),nhh,cspln,nw,
     .                 rplasbdry(j),zplasbdry(j),pds,ier,6)
          dr  = rmhdgrid(i ) - rplasbdry(j)
          dz  = zmhdgrid(nh) - zplasbdry(j)
          if (ideriv .eq. 2) then
              psi1d(kk) = psibdry+(pds(2)+0.5*pds(5)*dr)*dr
     .                   +(pds(3)+0.5*pds(6)*dz)*dz
     .                   +pds(4)*dz*dr
          else
              psi1d(kk) = psibdry+pds(2)*dr
     .                   +pds(3)*dz
          end if
c
c --- average around a number of plasma boundary points:
c
          if (navgbd .gt. 0) then
              ntotl = 2*navgbd+1
              isignn = -1
   30         do jl=1,navgbd
                 k = j + jl * isignn
                 if (k .lt. 1        )  k =     nplasbdry - k
                 if (k .gt. nplasbdry)  k = k - nplasbdry
                 call my_dbcevl1(rmhdgrid,nw,zmhdgrid(jsymetric),
     .                           nhh,cspln,
     .                       nw,rplasbdry(k),zplasbdry(k),pds,ier,6)
                 dr = rmhdgrid(i ) - rplasbdry(k)
                 dz = zmhdgrid(nh) - zplasbdry(k)
                 if (ideriv .eq. 2) then
                   psiav = psibdry + (pds(2) + 0.5*pds(5)*dr)*dr
     .                   + (pds(3) + 0.5*pds(6)*dz)*dz + pds(4)*dz*dr
                 else
                   psiav = psibdry + pds(2)*dr + pds(3)*dz
                 end if
                 psi1d(kk) = psi1d(kk) + psiav
              end do
              isignn    = -isignn
              if (isignn .eq. 1)  go to 30
              psi1d(kk) = psi1d(kk)/ntotl
          end if
c
c --- end of average
c
          psinocoil(kkk) = psi1d(kk)
c
          if (jsymetric .eq. 1) then
          kk1  = kk   - nh + 1    ! bottom boundary
          kkkk = kkkk - 1         ! bottom boundary
          j    = bdindex(kkkk)
          call my_dbcevl1(rmhdgrid,nw,zmhdgrid(jsymetric),nhh,cspln,nw,
     .                 rplasbdry(j),zplasbdry(j),pds,ier,6)
          dr = rmhdgrid(i) - rplasbdry(j)
          dz = zmhdgrid(1) - zplasbdry(j)
          if (ideriv .eq. 2) then
              psi1d(kk1) = psibdry+(pds(2)+0.5*pds(5)*dr)*dr
     .                   +(pds(3)+0.5*pds(6)*dz)*dz
     .                   +pds(4)*dz*dr
          else
              psi1d(kk1) = psibdry+pds(2)*dr + pds(3)*dz
          end if
c
c --- average around a number of plasma boundary points:
c
          if (navgbd .gt. 0) then
              ntotl  =  2 * navgbd + 1
              isignn = -1
   40         do jl=1,navgbd
                 k = j + jl * isignn
                 if (k .lt. 1        )  k =     nplasbdry - k
                 if (k .gt. nplasbdry)  k = k - nplasbdry
                 call my_dbcevl1(rmhdgrid,nw,zmhdgrid(jsymetric),
     .                           nhh,cspln,
     .                       nw,rplasbdry(k),zplasbdry(k),pds,ier,6)
                 dr = rmhdgrid(i) - rplasbdry(k)
                 dz = zmhdgrid(1) - zplasbdry(k)
                 if (ideriv .eq. 2) then
                      psiav = psibdry+(pds(2)+0.5*pds(5)*dr)*dr
     .                         +(pds(3)+0.5*pds(6)*dz)*dz
     .                         +pds(4)*dz*dr
                  else
                      psiav = psibdry+pds(2)*dr + pds(3)*dz
                  end if
                  psi1d(kk1) = psi1d(kk1)+psiav
              end do
              isignn     = -isignn
              if (isignn .eq. 1)  go to 40
              psi1d(kk1) = psi1d(kk1)/ntotl
          end if
c
c --- end of average
c
          psinocoil(kkkk) = psi1d(kk1)
          end if     ! end jsymetric = 1 branch
        end do       ! end do loop over index i
c
      return
c
      end

      subroutine psiscale (psiaxis, psibdry, psir, nj, isignpsi)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- rescale the psi array,psir(1,..nj),to the new axis and edge values,
c --- keeping the same fraction of the total poloidal flux inside the plasma
c ------------------------------------------------------------------ HSJ
c
      dimension psir(nj)
c
          psi1 = psir(1)
          if (isignpsi .lt. 0) then
              pmin = psiaxis
              pmax = psibdry
          else
              pmin = psibdry
              pmax = psiaxis
          end if
          do j=1,nj
              percnf = (psir(j)-psi1)/(psir(nj)-psi1)
              psir(j) = psiaxis+percnf*isignpsi*(pmin-pmax)
          end do
c
      return
c
      end

      subroutine setbdry (psi1d, nw, nh, nwh, ibd, psincbcd)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension  psi1d(*), psincbcd(*)
c
c --- vertical boundaries
c
      kk1 = 2*nh+nw-1
      do j=1,nh
        psi1d(j) = psincbcd(j)
c
        kk1 = kk1-1
        kk = (nw-1)*nh+j
        psi1d(kk) = psincbcd(kk1)
      end do
c
c --- horizontal boundaries ((rmhdgrid(i),zmhdgrid(1)),i = 2,..nw-1,
c --- and (rmhdgrid(i),zmhdgrid(nh)),i = 2,...nw-1):
c
        kkk  = nh               ! top    boundary
        kkkk = 2*(nh+nw) - 3    ! bottom boundary
        do i=2,nw-1
          kk  = i*nh            ! top    boundary
          kkk = kkk + 1         ! top    boundary
          psi1d(kk) = psincbcd(kkk)
c
          kk1  = kk - nh + 1    ! bottom boundary
          kkkk = kkkk - 1       ! bottom boundary
          psi1d(kk1) = psincbcd(kkkk)
        end do
c
      return
c
      end

      subroutine solvefbd (psi, psi1d, pcurrent, u0, isgngren, iounit)
c
      USE mhdpar
      USE mhdgrid
      USE bicube 
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- cyclic reduction solution. We solve del-star psi1d = rhs
c --- considering only the plasma current,with boundary conditions
c --- due only to the plasma current.
c --- first set up the finite difference grid coefficients.
c --- (we don't save the results so have to recalc each time)
c --- (may want to incorparte this directly into tri diagonal solver)
c ------------------------------------------------------------------ HSJ
c
c      include 'mhdpar.i'
c      include 'mhdgrid.i'
c      include 'bicube.i'
c
      dimension psi(nw,nh),psi1d(*),
     .           pcurrent(*)
      dimension diag(nw-2),diagl(nw-2),diagu(nw-2),diag1(nw-2),
     .          phi(nw-2),phi1(nw-2),v(nw-2),wk(nw-2),wk2(nw-2)
c      equivalence (diag(1),wnoperm(1)),(diagl(1),wnoperm(nw-1)),
c     .        (diagu(1),wnoperm(2*nw-3)),(diag1(1),wnoperm(3*nw-5)),
c     .        (phi(1),wnoperm(4*nw-7)),(phi1(1),wnoperm(5*nw-9)),
c     .        (v(1),wnoperm(6*nw-11)),(wk(1),wnoperm(7*nw-13)),
c     .        (wk2(1),wnoperm(8*nw-15))
c
      darea  = (rmhdgrid(2)-rmhdgrid(1))*(zmhdgrid(2)-zmhdgrid(1))
      dr     = rmhdgrid(2) - rmhdgrid(1)
      dz     = zmhdgrid(2) - zmhdgrid(1)    ! constant grid spacing
      dzsq   = dz*dz
      dzdrsq = (dz/dr)**2
      dumy   = dzsq/(2.0*dr)
      dumy1  = 2.0 * (1.0 + dzdrsq)
      do 200 i=1,nw-2
        diag(i) = dumy1
        denom = dumy/rmhdgrid(i+1)
        diagl(i) = -dzdrsq-denom
  200   diagu(i) = -dzdrsq+denom
c
c --- next set up the rhs term. Note that (-dzsq) was factored out
c --- of the finite difference expressions above and thus is included
c --- in the rhs term here. Division by darea is necessary because we
c --- have defined pcurrent as the current (in amps) of a filament of
c --- area darea. Here we require the actuall current density.
c --- the rhs is stored in the psi1d array because subroutine CMPLTCYR expects
c --- to find it there.
c
        dumy = dzsq*u0
        do 210 i=3,nw-2
          dumy1 = dumy*rmhdgrid(i)/darea
          do 210 j=2,nh-1
            kk = (i-1)*nh+j
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
        do 220 j=2,nh-1
          kk = nh+j
          psi1d(kk) = isgngren*dumy1*pcurrent(kk)+psi1d(j)*dumy
  220   continue
c
c --- for i = nw-1:
c
        dumy  = dzdrsq-dzsq/(2.0*rmhdgrid(nw-1)*dr)
        dumy1 = dzsq*u0*rmhdgrid(nw-1)/darea
        k = (nw-2)*nh
        kk = k+nh
        do 230 j=2,nh-1
          kkk = k+j
          kkkk = kk+j
          psi1d(kkk) = isgngren*dumy1*pcurrent(kkk)+psi1d(kkkk)*dumy
  230   continue
c
c --- some additional input required in cmpltcyr:
c
        nw1 = nw
        i   = 1
c
        do j=1,15
          i = i*2
          if (i .eq. nh-1) then
            nhpwr = j
            go to 250
          end if
        end do
c
        if (iounit .ne. 0)  write (iounit, 260) nh
  260   format (' ERROR: subroutine SOLVEFBD detects error, nh =', i5 /
     .          '        nh must be 2**m+1 for some integer m')
        ierr = 1
        return
c
  250   do 255 j=1,nh
         do 255 i=1,nw
           kk = (i-1)*nh+j
           psi(i,j) = psi1d(kk)
  255   continue
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
        do 265 j=1,nh
          do 265 i=1,nw
            kk = (i-1)*nh+j
            psi1d(kk) = psi(i,j)
  265   continue
c
        return
c
      end

      subroutine sorpicard (nw,nh,map,coeff,rmhdgrid,psi1d,u0,
     .                      psiaxis,psilim,curden,wzero,psi2d,
     .                      rplasbdry,zplasbdry,nplasbdry,
     .                      cspln,n2cspln,nh2,pds,zmhdgrid,wnoperm,
     .                      totcur,psidifbd,zero,source,
     .                      psiprevious,psinocoil,ieq,
     .                      bdindex,ideriv,navgbd,psibdry,
     .                      itermaxp,ipassmax,nj,psir,
     .                      outeritused,rma,zma,psidifbdext,
     .                      iextiters,updownsym,optwi,optomegai,optwe,
     .                      optomegae,anormal,residin,residout,
     .                      pprim,ffprim,elongax,psincbcd,kside,omeq)
c
      USE io,only : n66,n77
      implicit none
c
c ----------------------------------------------------------------------
c Sucessive over relaxation, Picard method of solving the set of equations
c using the exact plasma boundary with a rectangular finite difference
c grid. Points (interior to the plasma) adjacent to the boundary have
c the usual rectangular five point difference scheme modified to
c account for the presence of the boundary.
c this routine first solves the (nonlinear) Grad-Shafranov equation
c in the plasma interior. then psi is found in the exterior region
c between the plasma surface and the rectangular boundary of the
c computational grid.
c
c --- input
c   nw          number of grid points in radial direction
c   nh          number of grid points in axial (vertical)  direction
c               (nw and nh are not restricted to particular grid sizes)
c   map(i)      i = 1,2...nw*nh,a pointer vector,see subroutine mapvar
c rmhdgrid(i)   i = 1,2..nw horizontal grid
c zmhdgrid(j)   j = 1,2..nh vertical grid
c updownsym         =true  if case is up/down symmetric
c                   =false otherwise
c totcur       desired total current
c
c  psibdry     desired value of psi on plasma boundary.
c
c    coeff(1,l)        coefficient of psi ll (see subroutine mapvar)
c         (2,l)                           lt
c         (3,l)                           lr
c         (4,l)                           lb
c         (5,l)                           source term
c         (6,l)
c         (7,l)
c         (8,l)
c         (9,l)
c
c   optwi
c   optwe          optwi is set to true if optimum relaxation
c                  parameter for interations inside the plasma
c                  is to be found. optwe controls the same thing
c                  for iterations outside the plasma.
c   optomegai        value of relaxation parameter to use for
c                    iterations inside the plasma
c   optomegae        value of relaxation parameter to use for
c                    iterations outside the plasma
c    omeq              relaxation for outer,current iterations
c  residin          maximum allowed residual in inner iterations
c  residout         maximum allowed residual in outer iterations
c  ipassmax         maximum allowed number of outer iterations
c                   (for both the interior and exterior problem).
c  itermaxp         maximum allowed number of inner iterations
c                   (for both the interior and exterior problem).
c  u0               4.0 * pi*10e-7
c  ieq
c  psincbcd(kside)
c  rplasbdry(i)     i = 1,2..nplasbdry
c  zplasbdry(i)     i = 1,2..nplasbdry
c  nplasbdry    (rplasbdry(i),zplasbdry(i)) specifies the fixed plasma
c               boundary
c  cspln(n2cspln,nw,nh2)  used to store the bicubic spline coefficients
c  pds(6)          used to store the evaluation of bicubic spline
c  wnoperm(i)      i = 1,2,..    work vector for bicubic spline
c  wzero(i)        i = 1,2,..nw*nh  wzero(i) = 1.  inside plasma
c                  wzero(i) = 0.  outside plasma
c   zero(i)      i = 1,2,..nw*nh  zero(i) = 1.  inside limiter
c                 zero(i) = 0.  outside limiter
c   psiprevious(nw,nh)   a work storage array
c   psinocoil(i)         i = 1,2,..2*nw+2*nh-4 stores values of
c                        psi on boundary of rectangular grid
c   source(k)            k = 1,2,..nw*nh work storage used to hold
c                        rhs of Grad-Shafranov equation
c                        psi on boundary of rectangular grid
c   curden(nw,nh)        used to hold current density
c   ideriv        set to 1
c   navgbd        set to 1
c  pprim(j)          j = 1,2..nj
c  ffprim(j)         j = 1,2..nj the p prime and ffprime functions
c                    defined over the rho grid,psir(j).
c
c --- output:
c
c  anormal          current normalization factor
c
c  psi2d(nw,nh)     psi normalized to the desired input current totcur
c  psi1d(i)         i = 1,2...nw*nh,used to hold 1d version of psi
c  iextiters        number of outer iterations required to solve the
c                   exterior problem.
c  psidifbdext(i)   i = 1,2...iextiters the maximum residual on the
c                   rectangular boundary of the MHD grid as a
c                   function of outer iteration number
c  outeritused      number of outer iterations required to solve the
c                   plasma interior  problem.
c  psidifbd(i)      i = 1,2...outeritused the maximum outer residual for the
c                   plasma interior problem as a
c                   function of outer iteration number
c  bdindex(i)       i = 1,2..2*nw+2*nh-4 (integer vector)
c                   see subroutine load index for explanation
c  rma,zma          r,z,coordinates of magnetic axis
c  psiaxis          value of psi on magetic axis
c  psilim           value of psi on plasma boundary
c  psir(j)             j = 1,2..nj,new psi vector on transport rho grid
c  elongax          elongation on axis
c  other routines called by this subroutine:
c   findomegai   (calls inner_solution)
c   inner_solution
c   magax2    (calls ibcccu,dbcevl)
c   magax3    (calls icsicu1)
c   curdencalc  (calls psicalc,ffprimec,pprimec)
c                   (pprimec,ffprimec call find)
c   extrplt
c   findomegae (calls outer_solution)
c   outer_solution (calls coefcalc)
c   psibd
c   symcopy
c ---------------------------------------------------------- HSJ 3/26/93
c
      logical  interior, updownsym, optwe, optwi
c
      integer  k,i,j,nw,nh,map(*),outer_iterations,outermaxit,
     .         n2cspln,nh2,nplasbdry,ieqfail,outeritused,
     .         iknowax,ncrt,ispln,isignpsi,inner_iterations,
     .         maxinit,nwh,ideriv,navgbd,bdindex(*),iextiters,
     .         itermaxp,ipassmax,ieq,kside,jsymetric,nj,isymflg,iounit
c
      real*8   rmhdgrid(*),coeff(9,*),rl,psi1d(*),
     .         psibdry,psiaxis,psilim,u0,psirel,omega,omegasave,
     .         curden(nw,nh),psi2d(nw,nh),wzero(*),totc,psincbcd(kside),
     .         cspln(n2cspln,nw,nh2),wnoperm(*),pds(*),zmhdgrid(*),
     .         rplasbdry(*),zplasbdry(*),pcuriter,
     .         anormal,totcur,psidifbd(*),zero(*),dumy,ptrace,
     .         rma,zma,elongax,source(*),psiprevious(nw,nh),
     .         psinocoil(*),drdzgrid,psir(*),omext,omeq,
     .         deltar,deltaz,optomegai,optomegae,residiter,
     .         deltarsq,deltarzsq,residmax,residin,
     .         residout,psidifbdext(*),residbdry,psirelbd,
     .         pi,rhoest,drdzsq,pprim(*),ffprim(*),
     .         one
c
      include 'fitparms.i'

c
      data pi /3.14159265358979/
c
      one        = 1.0
      omext      = omeq
      if (omeq .le. 0.0)  omext = 0.5
      jsymetric  = 1
      if (updownsym)  jsymetric = nh/2+1
      deltar     = rmhdgrid(2)-rmhdgrid(1)
      deltaz     = zmhdgrid(2)-zmhdgrid(1)
      drdzgrid   = deltar*deltaz
      nwh        = nw*nh
      maxinit    = itermaxp
      psilim     = psibdry
      outermaxit = ipassmax
      iknowax    =  0
      isignpsi   = -1
      ncrt       =  6
      outer_iterations = 0
c
c --- initialize psi1d
c
      do i=1,nw
        do j=1,nh
          k        = (i-1)*nh+j
          psi1d(k) = psi2d(i,j)
        end do
      end do
c
c --- calculate initial current density
c
      call curdencalc2 (curden,nw,nh,rmhdgrid,wzero,psi2d,
     .                  psir,nj,u0,totc,pprim,ffprim,jsymetric)
      totc = totc*drdzgrid
c
c --- get initial source vector
c
      do j=jsymetric,nh
        do i=1,nw
          k = (i-1)*nh+j
          if (map(k) .ne. 0) then
            rl        = rmhdgrid(i)
            source(k) = u0*rl*curden(i,j)
          end if
        end do
      end do
c
c --- find optimal relaxation parameter for interior solution if called for:
c --- (finds optomegai such that the number of iterations required
c ---  for residmax to drop below the input value residiter is minimized).
c
      residiter = 1.0e-6
      if (optwi)       ! find optimal relaxation parameter, optomegai
     .      call findomegai (residiter,maxinit,inner_iterations,
     .                       updownsym,jsymetric,nw,nh,map,psi1d,coeff,
     .                       optomegai,rmhdgrid,psibdry,residmax,psi2d,
     .                       zmhdgrid,source)
c
c        at this point optomegai is a good approximation to the optimal
c        omega if optwi was set to true. otherwise optomegai is the input
c        value from the namelist. if the input value was not set then
c        use the Jacobi approximation to it:
c
      if (optomegai .eq. 0.0) then
        drdzsq = (rmhdgrid(2)-rmhdgrid(1)) / (zmhdgrid(2) - zmhdgrid(1))
        drdzsq = drdzsq*drdzsq
        rhoest = COS (pi/nw) + drdzsq * COS (pi/nh)
        rhoest = rhoest/(1.0 + drdzsq)
        if (rhoest .gt. 1.0)  rhoest = 1.0
        optomegai = 2.0/(1.0 + SQRT (1.0-rhoest**2))
        if (optomegai .lt. 1.0)  optomegai = 1.00
        if (optomegai .ge. 2.0)  optomegai = 1.99
        write (6, *) ' Jacobi estimation of omega =', optomegai
      end if
c
      omega = optomegai
c
c --- start the outer iterations
c
  110 outer_iterations = outer_iterations+1
      do j=jsymetric,nh
        do i=1,nw
          psiprevious(i,j) = psi2d(i,j)
        end do
      end do
c
c --- for the initial outer iterations we don't converge as tightly
c --- as for the final ones
c
      if (outer_iterations .gt. 5) then
        residiter = residin
        residiter = MIN (residiter, 0.01*residout)
      else
        residiter = residout
      end if
c
c --- do the inner iterations required to solve the problem:
c --- (does sor until residmax drops below residiter,using the fixed
c --- value omega.)
c
      call inner_solution (residiter,maxinit,inner_iterations,
     .                     updownsym,jsymetric,nw,nh,map,psi1d,coeff,
     .                     omega,psibdry,residmax,source)
c
c --- done with  inner iterations. Set up for next outer iteration:
c
c --- get psi2d:
c
      residmax = 0.0
      do j=jsymetric,nh
        do i=1,nw
           k = (i-1)*nh+j
           if (outer_iterations .gt. 1) then
              psi2d(i,j) = omext*psi1d(k) +(1.0-omext)*psiprevious(i,j)
              psirel     = ABS ((psi2d(i,j)-psiprevious(i,j))/psiaxis)
           else
              psi2d(i,j) = psi1d(k)
              psirel     = 1.0
           end if
           psi1d(k)      = psi2d(i,j)    ! for next pass
           if (map(k) .ne. 0)  residmax = MAX (residmax, psirel)
        end do
      end do
      psidifbd(outer_iterations) = residmax
c
c --- get the location of the magnetic axis,rma,zma, and the value of
c --- psi there,psiaxis(return with iknowax =1 if axis is found):
c
      if (updownsym) then
        call magax3 (psi2d,rmhdgrid,zmhdgrid,nw,nh,isignpsi,
     .               iknowax,ncrt,rma,zma,cspln,psiaxis,elongax,map)
      else
        call magax2 (psi2d,rmhdgrid,zmhdgrid,nw,nh,isignpsi,iknowax,
     .               ncrt,ispln,zero,rma,zma,cspln,n2cspln,nh2,
     .               psiaxis,elongax,map)
      end if
c
      if (iknowax .eq. 0) then        ! must always have a magnetic axis
c
c       we have a fatal error condition. copy upper half of psi
c       into lower half if necessary then dump the psi values and quit.
c
        isymflg = 0
        if (updownsym) then
          call symcopy (psi2d,psi1d,nw,nh,jsymetric)
          isymflg = 1
        end if
        write  (ncrt, 30)
   30   format (' ERROR: magnetic axis not found')
        call dump_psi_values (psi2d, nw, nh, rmhdgrid, zmhdgrid, n77,
     .           isignpsi, zero, isymflg, map, .true., dumy, 0, ptrace)
        call STOP ('subroutine SORPICARD: problem #2', 108)
      end if
c
c --- get new psi vector, psir
c
      call psiscale (psiaxis,psibdry,psir,nj,isignpsi)
c
c --- calculate new current density
c
      call curdencalc2 (curden,nw,nh,rmhdgrid,wzero,psi2d,
     .                  psir,nj,u0,totc,pprim,ffprim,jsymetric)
      totc     = totc*drdzgrid
      pcuriter = totc
      if (pcuriter .gt. 0.0)  go to 144
      ieqfail  = 1
      call STOP ('subroutine SORPICARD: problem #1', 39)
c
  144 do j=jsymetric,nh
          do i=1,nw
            k = (i-1)*nh+j
            if (map(k) .ne. 0) then
              rl = rmhdgrid(i)
c
c             evaluate the source term required for picard iteration:
c
              source(k) = u0*rl*curden(i,j)
            else
              source(k) = 0.0
            end if
          end do
      end do
c
****      write (6, *)  'iter-out,iter-in,residmax ', outer_iterations,
**** .                                      inner_iterations, residmax
****      write (6, *)  'psiaxis,rma,zma  ',psiaxis,rma,zma
****      write (6, *)  'totc,pcuriter,totcur ',totc,pcuriter,totcur
c
          if (outer_iterations .lt. outermaxit   ! do another outer iter
     .          .and. residmax .gt. residout)  go to 110
          outeritused = outer_iterations    ! done with outer iterations
          write (6, *)  'results for solution of interior problem:'
          write (6, *)  'iter-out,iter-in,residmax ', outer_iterations,
     .                                      inner_iterations, residmax
c
c ----------------------------------------------------------------------
c  the current normalization factor:
c
           if (updownsym)  pcuriter = 2.0 * pcuriter
           anormal = totcur / pcuriter
c
c ----------------------------------------------------------------------
c --- we are done with the interior solution. now,to be compatible
c --- with other methods of calculation,we need to get the solution
c --- in the exterior region between the plasma and the boundary
c --- of the rectangular grid.  To generate this exterior solution
c --- we need values of psi on the boundary of the rectangular grid.
c --- Inititialy we obtain these boundary values by analytic
c --- continuation,using taylor series,with derivatives evaluated
c --- from the interior solution only.
c --- load psi1d with boundary values of psi stored in psincbcd:
c --- psincbcd itself is loaded in subroutine GETMHDBC (i.e., get MHD boundary
c --- conditions):
c
c      if (ieq .eq. 0) then !HSJ 
      if (ieq .ge. 0) then   !always use the initial bc values
          call setbdry(psi1d,nw,nh,nwh,kside,psincbcd)
          do j=1,nh
            psi2d(1,j) = psi1d(j)
          end do
          do i=2,nw
            k = i*nh
            psi2d(i,nh) = psi1d(k)
          end do
          do j=nh-1,1,-1
            k = (nw-1)*nh+j
            psi2d(nw,j) = psi1d(k)
          end do
          do i=nw-1,2,-1
            k          = (i-1)*nh+1
            psi2d(i,1) = psi1d(k)
          end do
      else
        call extrplt(psi2d,coeff,map,nw,nh,jsymetric)
      end if
c
c ----------------------------------------------------------------------
c
c --- set up psi1d,which will be used to find the exterior solution.
c --- set up psiprevious,which will monitor the changes in the
c --- boundary values:
c
      do j=jsymetric,nh
        do i=1,nw
          k        = (i-1)*nh+j
          psi1d(k) = psi2d(i,j)
          psiprevious(i,j) = psi2d(i,j)
        end do
      end do
c
c --- now solve the exterior problem
c
      interior  = .false.
      deltarsq  = deltar*deltar
      deltaz    = (zmhdgrid(2)-zmhdgrid(1))
      deltarzsq = deltarsq/(deltaz*deltaz)
      omegasave = omega           ! save omega used in interior solution
c
c --- first get the optimum relaxation parameter for the exterior
c --- problem if requested:
c
      residiter = 1.0e-6
      if (optwe) then
        call findomegae (residiter,maxinit,inner_iterations,
     .                   updownsym,jsymetric,nw,nh,map,psi1d,coeff,
     .                   deltar,deltarzsq,deltarsq,interior,optomegae,
     .                   pds,rmhdgrid,psibdry,residmax,psi2d,zmhdgrid)
      end if
c
c --- if optwe is true then the omega we want to use is the
c --- value of optomegae just found above.
c --- if optwe is false then the value we want to use is either
c         a) the input value,if optomegae is input as non zero
c             or
c         b) the default value,if optomegae  = 0.0
c
      if (optomegae .gt. 0.0) then
        omega = optomegae ! from input, or from subroutine FINDOPTOMEGAE
      else
        omega = omegasave - 0.1     ! exterior omega slightly < interior
        omega = MAX (omega, one)
      end if
c
      outer_iterations = 0
  305 outer_iterations = outer_iterations + 1
      if (outer_iterations .gt. 5) then
        residiter = residin
      else
        residiter = residout
      end if
      call outer_solution(residiter,maxinit,inner_iterations,
     .            updownsym,jsymetric,nw,nh,map,psi1d,coeff,
     .            deltar,deltarzsq,deltarsq,interior,omega,pds,
     .            rmhdgrid,psibdry,residmax)
      residbdry =residmax  ! for ieq =0 case
c
c --- we are done with the inner iterations of the exterior solution.
c --- we now have psi over the entire rectangular grid.
c --- (or over the upper half if updownsym  = .true.)
c
       do j=jsymetric,nh
           do i=1,nw
               k = (i-1)*nh+j
               psi2d(i,j) = psi1d(k)
           end do
       end do
c
      if (ieq .ne. 0) then
c
c       now we can get new boundary values of psi on the rectangular
c       border of the grid (PSIBD calls the bicubic spline routines):
c
        ideriv = 1        ! ideriv = 2 does not work !
c
        call psibd (psi2d, rmhdgrid, zmhdgrid, rplasbdry, zplasbdry,
     .              nplasbdry, psibdry, psi1d, bdindex, ideriv,
     .              psinocoil, navgbd, jsymetric)
c
c       new boundary values are now stored in psinocoil and psi1d
c
        residbdry = 0.0
        do j=jsymetric,nh     ! up lhs
          psirelbd  = ABS ((psi1d(j)-psiprevious(1,j))/psiaxis)
          residbdry = MAX (residbdry,psirelbd)
          psiprevious(1,j) = psi1d(j)
        end do
        do j=nh,jsymetric,-1   ! down rhs
          psirelbd  = ABS ((psi1d(j)-psiprevious(nw,j))/psiaxis)
          residbdry = MAX (residbdry,psirelbd)
          psiprevious(nw,j) = psi1d(j)
        end do
        if (.not. updownsym) then
          do i=2,nw-1         ! across bottom
            k         = (i-1)*nh+1
            psirelbd  = ABS ((psi1d(k)-psiprevious(i,1))/psiaxis)
            residbdry = MAX (residbdry,psirelbd)
            psiprevious(i,1) = psi1d(k)
          end do
        end if
        do i=2,nw-1             ! across top
          k = i*nh
          psirelbd  = ABS ((psi1d(k)-psiprevious(i,nh))/psiaxis)
          residbdry = MAX (residbdry,psirelbd)
          psiprevious(i,nh) = psi1d(k)
        end do
        psidifbdext(outer_iterations) = residbdry  ! save for plotting
****    write (6, *)  'iter-in,residbdry =',inner_iterations,residbdry
c
c       now check for convergence of exterior solution:
c
        if (residbdry .gt. residout .and. outer_iterations .lt.
     .      outermaxit)  go to 305        ! do another outer iteration
c
      end if    ! ieq .ne. 0  branch
c
      write (6, *)  'for solution of exterior problem we have:'
      write (6, *)  'iter-in, residbdry=', inner_iterations, residbdry
c
c -------------------- RECTANGULAR BOUNDARY FITTING ----------------------------
c
      if (updownsym .and. fitboundary) then
        iextiters = 0
        iounit    = 6
        call fitdriver (iounit,rplasbdry,zplasbdry,nplasbdry,
     .                  residin,maxinit,inner_iterations,updownsym,
     .                  jsymetric,nw,nh,map,psi1d,coeff,deltar,
     .                  deltarzsq,deltarsq,interior,omega,psi2d,
     .                  rmhdgrid,zmhdgrid,psibdry,residmax,
     .                  psidifbdext,iextiters)
      end if
c
c --- finally, we are done with both inner and outer solutions.
c --- if this was an up/down symmetric case then expand the solution
c --- to cover the lower half of the grid. symcopy copies upper half of
c --- psi2d into lower half and also copies the resulting psi2d into psi1d
c
      if (updownsym)  call symcopy (psi2d, psi1d, nw, nh, jsymetric)
      return
c
      end

      subroutine symcopy (psi2d, psi1d, nw, nh, jsymetric)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- copy the upper half of the symmetric solution into the
c --- lower half. (this subroutine is similar to subroutine SYMMETRIZE,
c --- but does both 2d and 1d psi arrays).
c ----------------------------------------------------------------------
c
      dimension  psi2d(nw,nh), psi1d(*)
c
      jj = jsymetric
c
      do j=jsymetric-1,1,-1
        jj = jj + 1
        do i=1,nw
          psi2d(i,j) = psi2d(i,jj)
        end do
      end do
c
      do j=1,nh
        do i=1,nw
          k        = (i-1)*nh+j
          psi1d(k) = psi2d(i,j)
        end do
      end do
c
      return
c
      end

      subroutine tableintrp (xx, n, x, jlo)
c
      implicit none
c
c ----------------------------------------------------------------------
c --- correlated table search routine. use jlo from previous call to
c --- get jlo for current value of x. if jlo from previous call is
c --- no good(jlo=0 or jlo=n) then do a binary search.
c --- this routine returns jlo so that
c   xx(jlo) .le. x .le. xx(jlo+1) if jlo = 1
c   xx(jlo) .lt. x .le. xx(jlo+1) if jlo =2,3...n-1
c   it is assumed that xx(j),j = 1,2..n is monotonically
c   increasing; this is NOT checked for.
c   this is a modified version of the Numerical Recipes subroutine HUNT
c ------------------------------------------------------------------ HSJ
c
      integer  n, jlo, jhi, jmid, inc
      real*8   x, xx(n)
c
      if      (x .lt. xx(1)) then
        jlo = 0                     ! indicates out of range below xx(1)
      else if (x .le. xx(2)) then
        jlo = 1                     ! xx(1) .le. x .le. xx(2)
      else if (x .le. xx(n)) then   ! xx(2) .lt. x .le. xx(n)
c
c       check if jlo from previous call is usable:
c
        if (jlo .le. 0 .or. jlo .ge. n) then
          jlo = 2
          jhi = n
          go to 15           ! no correlation go directly to bisection
        else                 ! 1 .le. jlo .lt. n
c
c         bracket x, then use bisection:
c         start with  jlo  from previous call:
c
          inc = 1  
          if (x .gt. xx(jlo)) then    ! search up
    4       jhi = jlo + inc
            if (jhi .ge. n) then
              jhi = n
            else if (x .gt. xx(jhi)) then
              inc = inc + inc
              jlo = jhi
              go to 4
            end if
          else
    5       jhi = jlo
            jlo = jlo-inc
            if (jlo .le. 1) then
              jlo = 1
            else if (x .le. xx(jlo)) then
              inc = inc + inc
              go to 5
            end if
          end if
        end if
c
c         bisection:
c
   10     if (jhi-jlo .eq. 1)  return
   15     jmid = (jhi + jlo) / 2
          if (x .gt. xx(jmid)) then
            jlo = jmid
          else           ! x .le. xx(jmid)
            jhi = jmid
          end if
          go to 10
c
      else               ! x .gt. xx(n)
        jlo = n
      end if
c
      return
c
      end


