!------------------------------------------------------------
! gyro_collision_setup.f90 
!
! PURPOSE:
!  Make time-advance matrix, combining RBF expansion in 
!  space with Crank-Nicholson scheme in time.  
! 
! NOTES:
!  The parameter ord_rbf=1,3,5 controls the order of the 
!  RBF expansion.  
!
!  The irregular gridset in the (xi,theta) plane is expanded 
!  in an RBF spline series.
!
!  Near the boundaries at xi=-1 and xi=1, the SNaK method is 
!  used to improve accuracy.
!
!  The pitch-angle scattering operator is
!
!                d       2   d
!  C = nu_total --- (1-xi ) ---
!               dxi         dxi
!
!
!  C[v_par] = -2 nu_total v_par
!
!-------------------------------------------------------------

subroutine gyro_collision_setup

  use gyro_globals
  use gyro_collision_private
  use gyro_pointers
  use math_constants

  !-----------------------------------------------------------
  implicit none
  !
  integer :: isp
  integer :: p
  integer :: pp
  integer :: ic
  !
  real :: h_coll
  real :: x 
  real :: r1
  real :: r2
  real, dimension(n_rbf) :: xp
  real, dimension(n_rbf) :: yp
  real, dimension(n_rbf) :: xc
  real, dimension(4*n_rbf) :: work
  integer, dimension(n_rbf) :: iwork
  real, dimension(n_rbf,n_rbf) :: a1
  real, dimension(n_rbf,n_rbf) :: a2
  real, dimension(n_rbf,n_rbf) :: d1
  real, dimension(n_rbf,n_rbf) :: d2
  real, dimension(n_rbf,n_rbf) :: d3
  real :: rcond
  real :: anorm
  !
  real, external :: rbf,irbf
  real, external :: DLANGE
  !-----------------------------------------------------------

  nu_total(:,:,:) = 0.0

  condition_number = 0.0

  !-----------------------------------------------------------
  ! Electrons
  !
  do ie=1,n_energy

     x = sqrt(energy(ie))

     h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*erf(x)

     nu_total(:,ie,indx_e) = &
          0.5*(nu_s(indx_e,:)/x**3)*(z_eff_s(:)+h_coll)

  enddo
  !
  ! Ions
  !
  do is=1,indx_e-1
     do isp=1,indx_e-1

        do ie=1,n_energy

           x = sqrt(energy(ie))*mu(is)/mu(isp)

           h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*erf(x)

           x = sqrt(energy(ie))

           nu_total(:,ie,is) = nu_total(:,ie,is)+&
                0.5*(nu_s(is,:)/x**3)*h_coll*&
                den_s(isp,:)/den_s(is,:)*z(isp)**2

        enddo
     enddo
  enddo

  !-----------------------------------------------------------

  !-----------------------------------------------------------------
  ! Remove n=0 collisions?
  !
  if (kill_coll_flag == 1 .and. n_1(in_1) == 0) then
     nu_total(:,:,:) = 0.0
  endif
  !-----------------------------------------------------------------

  !------------------------------------------------
  ! Define xi = v_parallel/v, including sign of 
  ! velocity.
  !
  do k=1,n_lambda
     ck = class(k)
     do m=1,n_stack
        if (m <= n_theta(ck)) then
           xi(:,k,m) = -sqrt(abs(1.0-lambda(:,k)*b0_t(:,k,m))) 
        else
           xi(:,k,m) = sqrt(abs(1.0-lambda(:,k)*b0_t(:,k,m))) 
        endif
     enddo
  enddo
  !
  !--------------------------------------------------

  allocate(i_piv(n_rbf))

  p_ine_loc = 0

  do p_ine = 1+i_proc_1,n_ine_1,n_proc_1

     p_ine_loc = p_ine_loc+1

     i  = ine_i(p_ine)
     ie = ine_e(p_ine)

     ! Pack single-index RBF array
     ! -1 < xp < 1
     ! -1 <= yp < 1
     p = 0
     do k=1,n_lambda
        do m=1,n_stack

           p = p+1
           xp(p) = xi(i,k,m)
           xc(p) = xi(i,k,m)
           yp(p) = theta_t(i,k,m)/pi

           ! Shift centers near boundary according to SNaK 
           ! method for improved boundary accuracy and 
           ! significantly better time-advance stability. 
           if (k <= 2) then
              if (xp(p) < 0.0) then
                 xc(p) = -1.0-0.25*k
              else
                 xc(p) = 1.0+0.25*k
              endif
           endif

        enddo
     enddo

     ! Build RBF expansion matrix
     do p=1,n_rbf
        do pp=1,n_rbf
           a1(p,pp) = rbf(pi*(yp(p)-yp(pp)),xp(p)-xc(pp),ord_rbf)
        enddo
     enddo

     ! Dense (LAPACK) factorization
     anorm = DLANGE('I',n_rbf,n_rbf,a1,n_rbf,work)
     call DGETRF(n_rbf,n_rbf,a1,n_rbf,i_piv,info)
     if (info /= 0) then
        call catch_error('ERROR: Factorization in gyro_collision_setup')
     endif
     call DGECON('I',n_rbf,a1,n_rbf,anorm,rcond,work,iwork,info)

     if (1.0/rcond > condition_number) condition_number = 1.0/rcond

     ! Inverse (LAPACK)
     call DGETRI(n_rbf,a1,n_rbf,i_piv,work,n_rbf,info)
     if (info /= 0) then
        call catch_error('ERROR: Inversion in gyro_collision_setup')
     endif

     ! Original method
     do p=1,n_rbf
        do pp=1,n_rbf
           call drbf(pi*(yp(p)-yp(pp)),xp(p)-xc(pp),ord_rbf,r1,r2)
           ! Lorentz derivative  
           a2(p,pp) = (1-xp(p)**2)*r2-2*xp(p)*r1
        enddo
     enddo

     ! Perform matrix-matrix multiply with BLAS 
     ! to obtain derivative (Legendre) matrix: 
     call DGEMM('N',&
          'N',&
          n_rbf,&
          n_rbf,&
          n_rbf,&
          1.0,&
          a2,& 
          n_rbf,&
          a1,&
          n_rbf,&
          0.0,&
          d1,&
          n_rbf)

     if (linsolve_method == 3) then

        d1_rbf(:,:) = transpose(d1(:,:))
        goto 100

     else

        do ic=1,n_coll

           is = c_map(ic)

           d2(:,:) = 0.0
           d3(:,:) = 0.0
           ! Identity matrix
           do p=1,n_rbf
              d2(p,p) = 1.0
              d3(p,p) = 1.0
           enddo

           ! "Implicit" part (Crank-Nicholson)
           d2(:,:) = d2(:,:)-0.5*dt*nu_total(i,ie,is)*d1(:,:)
           ! "Explicit" part (Crank-Nicholson)
           d3(:,:) = d3(:,:)+0.5*dt*nu_total(i,ie,is)*d1(:,:)

           call DGETRF(n_rbf,n_rbf,d2,n_rbf,i_piv,info)
           call DGETRI(n_rbf,d2,n_rbf,i_piv,work,4*n_rbf,info)

           ! Perform matrix-matrix multiply with BLAS 
           ! to obtain collisional time-advance matrix:
           call DGEMM('N',&
                'N',&
                n_rbf,&
                n_rbf,&
                n_rbf,&
                1.0,&
                d2,&
                n_rbf,&
                d3,&
                n_rbf,&
                0.0,&
                a1,&
                n_rbf)

           ! Take a transpose so collision advance has
           ! optimal index order.
           d_rbf(:,:,p_ine_loc,ic) = transpose(a1(:,:))

        enddo

     endif

  enddo

100 continue

  deallocate(i_piv)

  if (i_proc == 0) call gyro_collision_grid_write(trim(path)//'out.gyro.collision_grid',1)

  if (i_proc == 0 .and. debug_flag == 1) &
       print *,'[gyro_collision_setup]'

end subroutine gyro_collision_setup

real function rbf(dt,dy,s)

  implicit none

  real, intent(in) :: dt,dy
  integer, intent(in) :: s
  real :: rl

  rl  = sqrt(2.0-2.0*cos(dt)+dy**2)
  rbf = rl**s

end function rbf

subroutine drbf(dt,dy,s,r1,r2)

  implicit none

  real, intent(in) :: dt,dy
  integer, intent(in) :: s
  real, intent(inout) :: r1,r2
  real :: rl

  rl = sqrt(2.0-2.0*cos(dt)+dy**2)

  r1 = s*rl**(s-2)*dy

  if (rl > 0.0) then
     r2 = (s*(s-2)*dy**2+s*rl**2)*rl**(s-4)
  else
     r2 = 0.0
  endif

end subroutine drbf
