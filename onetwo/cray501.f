
      subroutine adaptive_func (x, njdum, dzeta, f)
c
c ----------------------------------------------------------------------
c     sets up the weight function required by adaptive_grid_gen.
c
c INPUT
c     x(j)      j=1,..nj the  variable for
c               which a new grid is to be found by
c               adaptive_grid_gen.
c     njdum        number of grid pts
c     dzeta     uniform increment in zeta grid
c
c INPUT from include files:
c       see list below next to include statements
c
c OUTPUT
c      f(j)      j=1,..nj-1 is the value of 1./w  on the half grid
c
c ----------------------------------------------------------3/7/97--HSJ
c
      USE param
      USE     soln
      USE numbrs
      USE adaptive
      implicit none
c
      integer njdum, i, j         ! njdum=nj in numbrs.i
      real*8  dzeta,dxdzeta,dfdzeta,d2xdzeta2,d2fdzeta2,
     .        dzetasq,d2fdzeta,dfdzetau,curve_weight,curve,
     .        x(*),f(*)
c
c      include 'param.i'      ! kk,kj,etc
c      include 'adaptive.i'   ! eps_adaptive
c      include 'numbrs.i'     ! nion
c      include 'soln.i'       ! u(kk,kj)
c
      data  curve_weight / 5.0 /   ! set in cray102.f or inone
****  include_curvature=.true.     ! set in cray102.f or inone
c
      dzetasq=dzeta*dzeta
      do j=1,nj-1       ! on the half grid
         dxdzeta=(x(j+1)-x(j))/dzeta
         if (include_curvature .gt. 0) then
            if (j .gt. 1 .and. j .lt. nj-1) then
                 d2xdzeta2=(x(j+1)-2.*x(j)+x(j-1))/dzetasq
            else if (j .eq. 1) then
                 d2xdzeta2=(x(3)-2.*x(2)+x(1))/dzetasq
            else ! j=nj-1
                 d2xdzeta2=(x(nj)-2.*x(nj-1)+x(nj-2))/dzetasq
            end if
         end if
         dfdzeta=0.0
         d2fdzeta=0.0
         do i=1,nion+4
             dfdzetau=(u(i,j+1)-u(i,j))/dzeta
             dfdzeta=dfdzeta+eps_adaptive(i)*dfdzetau**2
             if (include_curvature .gt. 0) then
                 if (j .gt. 1 .and. j .lt. nj-1) then
                   d2fdzeta2=(u(i,j+1)-2.*u(i,j)+u(i,j-1))/dzetasq
                 else if (j .eq. 1) then
                   d2fdzeta2=(u(i,3)-2.*u(i,2)+u(i,1))/dzetasq
                 else ! j=nj-1
                   d2fdzeta2=(u(i,nj)-2.*u(i,nj-1)+u(i,nj-2))/dzetasq
                 end if
                 d2fdzeta = d2fdzeta+d2fdzeta2*curve_eps_adaptive(i)
              end if
         end do
         f(j) = SQRT (dxdzeta**2+dfdzeta) / ABS (dxdzeta) ! reciprocal..
c                                                         ..of weight
         if (include_curvature .gt. 0) then
           curve = ABS (dxdzeta*d2fdzeta - SQRT (dfdzeta) * d2xdzeta2)
     .           / ((dxdzeta**2+dfdzeta)**1.5)
           curve = (1.0 + curve_weight*curve)
           f(j)  = curve*f(j)
         else
           continue
         end if
c
      end do
      return
c
      end
      subroutine adaptive_grid_gen (sol, nj, rhoa)
c ------------------------------------------------------------------
c --- subroutine generates a grid adapted to some aspect
c --- of the solution as determined by subroutine adaptive_func
c ---
c  INPUT
c  sol(j)         j=1,2,...nj the old sol
c  nj             number of grid points Note that nj includes
c                 sol(1) and sol(nj),which are actually boundary
c                 conditions on the solution.
c  rhoa           the value of rho  at the plasma boundary in cm
c                 at this time. The time dependence of rhoa
c                 is determined by the mhd calcs in the usual way
c  From include files:
c  dzeta          the normalized [0,1] grid has uniform
c                 grid spacing given by dzeta
c
c  OUPUT
c  sol(j)         j=1,2,...nj the new solution (i.e., grid) addapted to
c                 the weight function given by adaptive_func
c  sol(1)         on output is the same as sol(1) on input
c  sol(nj)        on output is the same as sol(1) on input
c
c  EXTERNAL REFERENCES
c  tridiag
c  sgtsl           tridiagonal system solvers
c  adaptive_func  routine that will supply the
c  values F_i
c
c  Method:
c        the non linear dif eq.
c             d((1./w)(dx/dz))/dz = 0         (see notes)
c          with bc x(1)=a,x(nj)=b (see below,a=0.0,b=rhoa)
c        is finite differenced and repeatedly solved until
c        the nonlinearity due to w dependence on x is
c        converged. At grid point i we call (1/w) = F_i
c        F_i is linearized by using the values of x_i initially from
c        the input guess and then from the previous iteration
c        after the intial solution.
c        For a given fixed zeta grid spacing the tridiagonal
c        matrix used is diagonally dominant only if
c
c --------------------------------------------------- 3/4/97 -- HSJ ----
c
      USE param
      USE verbose
      USE adaptive
      implicit none
c
      include 'storage.i'
c
      logical    not_converged
      integer    iloop_max
      parameter (iloop_max = 100)    ! max loops to resolve nonlinearity
      integer    nj,j,ierr,iloop,valid_found,jm1,
     .           njm1,njm2,jmin,jmax,jmax_diffl,jmax_diff(iloop_max)
      real*8     grid(kj),sub(kj),diag(kj),supr(kj),rhs(kj),
     .           sol(*),old_sol(kj),work(kj),f(kj),min_dr,max_dr,
     .           dzeta,dzetaov2,max_diff(iloop_max),tolerance,a,b,
     .           dsol,max_diffl,diffl,dr_adaptive,dr_adaptive_crit,
     .           r_input(kj),omega,omega1,rhoa,ovwl,ovwr
      data       tolerance, dr_adaptive_crit /1.0e-4, 0.005/ ! max error
      data       omega /0.5/
c
c     note dr_adaptive_crit is in cm (unfortunately),so we are saying that
c     the grid should not become more dense than 0.01 cm between grid points
c     this routine will save the best valid solution and use it if the
c     solution of the nonlinear ode doesnt converge.
c
      equivalence (tdum(1),grid(1))
      equivalence (tdum(kj+1),sub(1))
      equivalence (tdum(2*kj+1),diag(1))
      equivalence (tdum(3*kj+1),supr(1))
      equivalence (tdum(4*kj+1),rhs(1))
      equivalence (tdum(5*kj+1),old_sol(1))
      equivalence (tdum(6*kj+1),work(1))
      equivalence (tdum(7*kj+1),f(1))
      equivalence (tdum(8*kj+1),r_input(1))
c
      if (include_adaptive .le. 0)  return
      omega1        = 1.0 - omega
      not_converged = .true.
      call copya (sol, r_input, nj) ! save in case of failure
      iloop         = 0
      dzeta         = dzeta_adaptive
      dzetaov2      = 0.5 * dzeta
      a             = sol(1) ! solution will be from sol(2) to sol(nj-1)
      b             = rhoa
      sol(nj)       = b      ! rhoa may be different from soln(nj)
      njm2          = nj-2   ! number of equations to be solved
      njm1          = nj-1   ! from j = 2 to j = nj-1
      valid_found   = 0
c
      do while (not_converged)  ! ======================================
c
        iloop = iloop+1
        if (iloop .gt. 1) then    ! sol = omega*sol+(1-omega)*old_sol
          call multpl1 (    sol, nj, omega )
          call multpl1 (old_sol, nj, omega1)
          call    adda (old_sol, sol, nj)
        end if
        call adaptive_func(sol,nj,dzeta,f)
c
c     set up the tridiagonal matrix:
c        first equation:
c
        diag(1) =  f(1)
        rhs (1) = -a*diag(1)
        supr(1) =  f(2)
        diag(1) = -supr(1)-diag(1)
c
c     second through penultimate equations:
c
        do j=3,njm2
          jm1=j-1
          ovwl=f(jm1)
          ovwr=f(j)              ! w(j) is on the half grid
          sub(jm1)=ovwl          ! sub diagonal,not used for j=1
          diag(jm1)=-(ovwl+ovwr) ! diagonal
          supr(jm1)=ovwr         ! supper diagonal,not used for j=njm2
          rhs(jm1)=0.0
        end do
c
c       last equation:
c
        diag(njm2) =  f(njm1)
        rhs(njm2)  = -b*diag(njm2)
        sub(njm2)  =  f(njm2)
        diag(njm2) = -sub(njm2)-diag(njm2)
c
        call copya (sol, old_sol, nj) ! copy sol into old_sol
c
c       solve the nj-2 equations,store solution in sol(2) to sol(nj-1)
c       sol(1) and sol(nj) are fixed by boundary conditions:
c
        call tridiagHSJ (sub, diag, supr, rhs, sol(2), work, njm2,ierr)
        if (ierr .gt. 0) then ! encountered zero pivot..
c                             ..use partial pivoting
          ierr = 0
          call sgtsl (njm2, sub, diag, supr, rhs, ierr)
          if (ierr .ne. 0)
     .      call STOP ('subroutine GRIDGEN: IERR is non-zero', 244)
          call copya (rhs,sol(2),njm2) ! copy rhs into vector sol
        end if
c
c       check convergence and validity of solution
c
        max_diffl = 0.0
        do j=2,njm1
          dr_adaptive = sol(j)-sol(j-1)
          if (dr_adaptive .lt. dr_adaptive_crit) then ! trouble, maybe
            iloop = iloop_max + 1                     ! not monotonic
            go to 25                                  ! stop searching
          end if
          diffl = abs ((old_sol(j)-sol(j))/sol(j))
          max_diffl = MAX (max_diffl, diffl)
          if (max_diffl .eq. diffl)  jmax_diffl = j
        end do
        valid_found = 1
        call copya (sol, r_adaptive, nj) ! save best so far
        if (max_diffl .lt. tolerance)  not_converged=.false.
         max_diff(iloop) =  max_diffl
        jmax_diff(iloop) = jmax_diffl
   25   if (iloop .ge. iloop_max) then ! something is wrong
            if (gridgenvb .eq. 1 .and. iloop .ne. iloop_max+1) then
              do j=1,iloop_max
                write (*,'("j,   max_diff in rho =",i5,3x,1pe12.4)')
     .          jmax_diff(j),max_diff(j)
              end do
              write (*,'("grid calcs did not converge, tolerance =",
     .                    1pe12.4)') tolerance
              write (*,'("last rho grid estimate:" / (5(2x,1pe14.6)))')
     .                                               (sol(j),j=1,nj)
            end if
            if (valid_found .eq. 1) then
              call copya(r_adaptive,sol,nj) ! restore best solution
            else
              call copya (r_input, sol, nj) ! restore previous solution
            end if
            not_converged = .false.         ! allow exit
         end if
c
      end do                ! ==========================================
c
      if ((sol(1) .ne. a) .or. (sol(nj) .ne. b))
     .  call STOP ('subroutine GRIDGEN: unspecified problem', 245)
      if (gridgenvb .eq. 1) then
         min_dr =  sol(nj)
         max_dr = -sol(nj)
         do j=2,njm1
            dsol   = sol(j)-sol(j-1)
            min_dr = MIN (min_dr, dsol)
            max_dr = MAX (max_dr, dsol)
            if (max_dr .eq. dsol)  jmax = j
            if (min_dr .eq. dsol)  jmin = j
c            print *,'j,sol(j)',j,sol(j)
         end do
         write (*, '("min_dr,j =",1pe12.2,1x,i5)')  min_dr, jmin
         write (*, '("max_dr,j =",1pe12.2,1x,i5)')  max_dr, jmax
      end if
      
      return
c
      end

      subroutine intrp_adp (r_in, f_in, n_in, r_out, f_out, n_out)
c
c -------------------------------------------------------------------
c     interpolate f_in onto f_out using splines or ?
c -------------------------------------------------------------------
c
      USE param
      USE replace_imsl,            ONLY : my_icsicu,my_icsevu
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'   ! kj
c
      dimension  r_in(*), f_in(*), r_out(*), f_out(*)
      dimension  bpar(4), c(kj,3)
c
      bpar(1) = 1.0       ! zero gradient at rho = 0
      bpar(2) = 6.0 *(f_in(2)-f_in(1))/((r_in(2)-r_in(1))**2)
      bpar(3) = 0.0       ! natural at rho = 1
      bpar(4) = 0.0
      ier     = 0
      call icsicu1 (r_in, f_in, n_in, bpar, c, kj, ier)
      if (ier .ne. 0)  go to 10
c      call icsevu1 (r_in, f_in, n_in, c, kj, r_out, f_out, n_out, ier)
      call my_icsevu (r_in, f_in, n_in, c, kj, r_out, f_out, n_out, ier)
   10 if (ier .eq. 0)  return
      call STOP ('subroutine INTRP_ADP: non-zero IER from IMSL', 247)
c
      end

      subroutine meshgen (r, roa, ra, dr, rrp, drr, rrm, nj, codeid,
     .                    time, rhod_ifs, rhod_max_ifs)
c
c ----------------------------------------------------------------------
c     isolate the mesh calcs to this subroutine
c ------------------------------------------------------- 3/10/97 HSJ --
c
      implicit none
c
      integer   j, nj
      real*8    r(*), roa(*), ra(*), drr(*), rrp(*), dr(*), rrm(*),
     .          time, rhod_max_ifs, rhod_ifs(*)
      character codeid*(*)
c
      do j=1,nj-1
        dr(j)  = r(j+1)-r(j)
        roa(j) = r(j)/r(nj)
        ra (j) = 0.5*(r(j)+r(j+1))
      end do
c
      roa(nj) = 1.0
      if (codeid .eq. 'onedee')  call rhomsh (time)
c
c     MESHGEN is called from RHOMSH only if codeid .ne. 'onedee'
c
      drr(1) = 2.0 / dr(1)
      rrp(1) = 2.0 / dr(1)
c
      do  j=2,nj
        rrm(j)  = (ra(j-1)/r(j))/dr(j-1)
        if (j .ne. nj) then
          drr(j)  = 2.0/(dr(j-1)+dr(j))
          rrp(j)  = (ra(j)/r(j))/dr(j)
        end if
      end do
c
      drr(nj) = 2.0 / dr(nj-1)
c
      if (codeid .eq. 'onedee') then
        call copya(r,rhod_ifs,nj)
        rhod_max_ifs = rhod_ifs(nj)
        call multpl1 (rhod_ifs, nj, 1.0/rhod_max_ifs)
      end if
      return
c
      end

      subroutine meshgeninitial (imesh, reqsp, kappa, rminor, reff,
     .                           r, dr, widths, delr, dzeta_adaptive,
     .                           include_adaptive, nout, ncrt, kj,
     .                           nj, njs, nps, ierr)
c
c ----------------------------------------------------------------------
c compute mesh quantities
c ----------------------------------------------------------------------
c
      USE geom,                     ONLY : codeid
      implicit none
c
      logical reqsp
      integer kj,nj,njs,nps,j,imesh,nout,ncrt,ierr,include_adaptive
      real*8  kappa,rminor,reff,delr,dzeta_adaptive,dr(*),r(*),
     .        dr2tmp,r2tmp,r2tmpn,widths,tol,sum,xmin
c
      reqsp = .false.
      reff =rminor
      if(codeid == 'onedee') 
     .      reff  = SQRT (kappa) * rminor
      if (imesh .eq. 3)  go to 3265
      if (imesh .eq. 2)  go to 3270
      if (imesh .eq. 1)  go to 3290
c
c ----------------------------------------------------------------------
c default is uniform mesh
c ----------------------------------------------------------------------
c
      delr  = reff/(nj-1)
      dzeta_adaptive = 1./(nj-1)
      do 3260 j=1,nj-1
      dr(j) = delr
 3260 r(j)  = (j-1)*delr
      reqsp = .true.
      go to 3310
c
c ----------------------------------------------------------------------
c equally spaced in r**2
c ----------------------------------------------------------------------
c
 3265 r(1)   = 0.0
      dr2tmp = reff**2/(nj-1)
      r2tmp  = 0.0
      do 3266 j=1,nj-1
      r2tmpn = r2tmp+dr2tmp
      dr(j)  = SQRT (r2tmpn) - SQRT (r2tmp)
      r2tmp  = r2tmpn
 3266 r(j+1) = r(j)+dr(j)
      go to 3310
c
c ----------------------------------------------------------------------
c nonuniform mesh with dr(j) specified
c ----------------------------------------------------------------------
c
 3270 r(1)   = 0.0
      do 3280 j=1,nj-1
 3280 r(j+1) = r(j)+dr(j)
      go to 3310
c
c ----------------------------------------------------------------------
c nonuniform mesh with r(j) specified
c ----------------------------------------------------------------------
c
 3290 if(r(nj) .le. 1.0)then  !input grid is on [0,1.0] instead of [0.0,rminor]
         do j=1,nj
            r(j) = rminor*r(j)
         enddo
      endif
      r(1)  = 0.0
      do 3300 j=1,nj-1
 3300 dr(j) = r(j+1)-r(j)
c
c ----------------------------------------------------------------------
c make sure r(nj) = reff
c ----------------------------------------------------------------------
c
 3310 r(nj)    = reff
      dr(nj-1) = reff-r(nj-1)
c
c ----------------------------------------------------------------------
c  set up mesh for scrape-off layer
c ----------------------------------------------------------------------
c
      njs     = nj
      if (widths .le. 0.0)  go to 3318
      njs     = nj+nps
      delr    = 0.9*widths/(nps-1)
      do 3315 j=1,nps
 3315 r(nj+j) = r(nj) + 0.1*widths + delr*(j-1)
      reqsp   = .false.
c
c ----------------------------------------------------------------------
c check number of mesh points
c ----------------------------------------------------------------------
c
 3318 if (njs .le. kj)  go to 1970
      write (nout, 7970) kj
      write (ncrt, 7970) kj
 7970 format (' *** ERROR: njs must be <= ',i8)
      ierr = 1
c
c ----------------------------------------------------------------------
c check mesh
c ----------------------------------------------------------------------
c
 1970 sum  = 0.0
      tol  = rminor * 1.0e-6
      xmin = rminor / 10000.0
      xmin = rminor / 100000.0     !changed 07/02, hsj, seems arbitrary at any rate ??
c
      do 1980 j=2,nj
        delr = r(j) - r(j-1)
        sum  =  sum + delr
        if (delr .gt. xmin)  go to 1980
        write  (nout, 7980)  j, r(j-1), r(j)
        write  (ncrt, 7980)  j, r(j-1), r(j)
 7980   format (' *** ERROR: r(j)-r(j-1) is too small' /
     .          ' j,r(j-1),r(j)',i8,2x,2e12.3)
        ierr = 1
 1980 continue
c
      if (ABS (sum - SQRT (kappa)*rminor) .lt. tol)  go to 1990
      write  (nout, 7990  )
      write  (ncrt, 7990  )
 7990 format (' *** ERROR: sum of delr s # rminor ')
      ierr = 1
c
c ----------------------------------------------------------------------
c calculate some mesh-related quantities.
c ----------------------------------------------------------------------
c
 1990 if (include_adaptive .gt. 0)  reqsp = .false.


      return
c
      end

      subroutine set_rho (times)
c
c -------------------------------------------------------------------
c
      USE param
      USE solcon
      USE     soln
      USE mhdpar
      USE numbrs
      USE mesh
      USE adaptive
      USE geom
      USE soln2d
      USE ifs
      USE neo2d
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'adaptive.i'
c      include 'mhdpar.i'   ! kpsi
c      include 'ifs.i'      ! rhod_ifs (needs kpsi from mhdpar.i)
c      include 'mesh.i'
c      include 'geom.i'     ! rhoa, rhoa_save
c      include 'numbrs.i'
c      include 'neo2d.i'    ! eps, etc
c      include 'solcon.i'
c      include 'soln2d.i'   ! implicit_fh
      include 'storage.i'
c      include 'soln.i'
c
      dimension roa_old(kj), rhod_ifs_old(kj)
      real*8 zero
c
      equivalence (udum(1   ), roa_old     (1))
      equivalence (udum(kj+1), rhod_ifs_old(1))
c
      data speed_adaptive_min /1.0e-10/
c
      if (times .gt. freeze_adaptive) then
        zero             = 0.0
        include_adaptive = 0
        call multpl1 (drhodt_adaptive, nj, zero)
      end if
      if (include_adaptive .eq. 0)  return
c
c     given the current normalized r grid in roa,which
c     corresponds to the mid plane grid rhod_ifs,
c     get a new rhod_ifs grid that corresponds to the
c     new r grid
c
      call copya(roa,roa_old,nj)
      call copya(rhod_ifs,rhod_ifs_old,nj)
      call adaptive_grid_gen(r,nj,rhoa)  ! redefines r
 100  drmin = r(nj)
c
      do j=2,nj-1
c          r(j)  = 0.25*r(j)+ 0.75*r_mesh(j) 
           r(j)  = dt*1000.0*speed_adaptive*(r(j)-r_mesh(j))+r_mesh(j)
           drmin = MIN (r(j)-r(j-1), drmin)
      end do
c
      drmin_crit=1.e-5 ! just a guess
      if (drmin .lt. drmin_crit) then
cHSJ        speed_adaptive = speed_adaptive*0.5
        if (speed_adaptive .gt. speed_adaptive_min)  go to 100
        call STOP ('subroutine SET_RHO: speed_adaptive too small', 246)
      end if
c
      call meshgen (r,roa,ra,dr,rrp,drr,rrm,nj,codeid,times,
     .              rhod_ifs,rhod_max_ifs)
c
c        get other profiles on adaptive grid:
c
****     call intrp_adp (roa_old, rhod_ifs_old, nj, roa, rhod_ifs, nj)
c
         call intrp_adp (roa_old, gcap0, nj, roa, gcap, nj)
         if (.not. implicit_fh) then
            call copya(hcap,hcap0,nj)
            call copya(fcap,fcap0,nj)
            call intrp_adp (roa_old, fcap0, nj, roa, fcap , nj)
            call intrp_adp ( roa_old, hcap0, nj, roa, hcap , nj)
         end if
         call copya(r2cap,r2cap0,nj)
         call copya(rcap,rcap0,nj)
         call copya(rcapi,rcap0i,nj)
         call copya(r2capi,r2capi0,nj)
         call intrp_adp (roa_old, r2cap0,  nj, roa, r2cap  ,nj)
         call intrp_adp (roa_old, r2capi0, nj, roa, r2capi ,nj)
         call intrp_adp (roa_old, rcap0,   nj, roa, rcap   ,nj)
         call intrp_adp (roa_old, rcap0i,   nj, roa, rcapi ,nj)
         if (xi_include) then
             call copya(eps,eps0,nj)
             call copya(xhm2,xhm20,nj)
             call copya(xi11,xi110,nj)
             call copya(xi33,xi330,nj)
             call copya(xips,xips0,nj)
             call intrp_adp (roa_old, eps0,  nj, roa, eps  ,nj)
             call intrp_adp (roa_old, xhm20, nj, roa, xhm2 ,nj)
             call intrp_adp(roa_old, xi110, nj, roa, xi11 ,nj)
             call intrp_adp (roa_old, xi330, nj, roa, xi33 ,nj)
             call intrp_adp (roa_old, xips0, nj, roa, xips ,nj)
         end if
c
      drtmax = 0.0
      drmin  = r(nj)
            do j=1,nj
               drhodt_adaptive(j) = (r(j)-r_mesh(j))/dt
               drtabs = ABS (drhodt_adaptive(j))
               drtmax = MAX (drtmax, drtabs)
               if (drtmax .eq. drtabs)  jmax = j
               if (j .lt. nj) then
                  drjd  = r(j+1) - r(j)
                  drmin = MIN (drmin, drjd)
                  if (drmin .eq. drjd)  jmind = j
               end if
            end do
      write (6, '("drmin,j =",1pe12.2,1x,i5 /
     .            "max drhodt,j",1pe12.2,1x,i5)')drmin,jmind,drtmax,jmax
          if (drtmax .gt. 100.0) then
cHSJ         speed_adaptive = speed_adaptive * 0.5
         write (6, '(" ********************* drtmax too large "     /
     .               " ********************* speed_adaptive halved" /
     .               " ********************* new value = ", 1pe12.4)')
     .                 speed_adaptive
      end if
      return
c
      end
