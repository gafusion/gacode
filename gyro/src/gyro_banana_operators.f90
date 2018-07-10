!-----------------------------------------------------
! gyro_banana_operators.f90
!
! PURPOSE:
!  Compute arrays connected with solution in theta 
!  direction. 
! 
! NOTES:
!
! Pointers are allocated with enough tau's for
! trapped particles.  The alternative is to 
! have separate theta_t's and tau's for trapped 
! and passing particles.  We don't want that.
!
!  - m_map(ck,j,p) 
!
!   returns index along orbit (tau-direction) for a 
!   given real position in theta (j) and sign of 
!   velocity.
!
!  - theta_t(k,m)
!
!   the physical poloidal angle, theta, corresponding 
!   to tau(k,m) (below).
!
!  - tau(k,m)
!  
!   normalized time along orbit.
!
!  - m_cyc(n_class,<wider than n_stack>,2)
!
!   Cyclic version of m-index.
!
!  - p_cyc(n_class,<wider than n_stack>,2)
!
!   Cyclic version of phase multiplier.
!-----------------------------------------------------

subroutine gyro_banana_operators

  use gyro_globals
  use math_constants

  !-------------------------------------------
  implicit none
  !
  integer :: p
  !
  real :: period
  real :: banana_theta(n_tau(1)+1)
  real :: banana_tau(n_tau(1)+1)
  !-------------------------------------------

  m_cyc = 0

  ! Physical signs of velocity (needed only for finite-beta coding).

  ! ** Passing **

  do i=1,n_x

     call gyro_to_geo(i)

     call gyro_banana_init(nint_ORB_s)

     do k=1,n_pass

        class(k) = 1

        ck = class(k)

        sigma(ck,1) = -1.0
        sigma(ck,2) = 1.0

        do m=1,n_stack
           if (m <= n_theta(ck)) then
              m_phys(ck,m) = m
              p_phys(ck,m) = 1
           else
              m_phys(ck,m) = m-n_theta(ck)
              p_phys(ck,m) = 2
           endif
        enddo

        do j=1,n_theta(ck)

           m_cyc(ck,j+n_theta(ck),1) = j
           p_cyc(ck,:,j+n_theta(ck),1) = phase(in_1,:)

           m_cyc(ck,j-n_theta(ck),1) = j
           p_cyc(ck,:,j-n_theta(ck),1) = 1.0/phase(in_1,:)

           ! Primary (v<0)
           m_cyc(ck,j,1) = j
           p_cyc(ck,:,j,1) = 1.0

           m = j+n_theta(ck)

           m_cyc(ck,m+n_theta(ck),2) = m
           p_cyc(ck,:,m+n_theta(ck),2) = phase(in_1,:)

           m_cyc(ck,m-n_theta(ck),2) = m
           p_cyc(ck,:,m-n_theta(ck),2) = 1.0/phase(in_1,:)

           ! Primary (v>0)
           m_cyc(ck,m,2) = m
           p_cyc(ck,:,m,2) = 1.0

        enddo

        do j=1,n_theta(ck)
           m_map(ck,j,1) = j
           m_map(ck,j,2) = j+n_theta(ck)
        enddo

        call gyro_banana_uniform_taugrid(lambda(i,k),&
             n_tau(1)+1,&
             nint_ORB_do,&
             banana_theta,&
             banana_tau)

        ! Parameterization of entire orbit, given 
        ! the passing orbit computed by ORB library.

        period = banana_tau(n_tau(1)+1)

        do j=1,n_tau(1)

           tau(i,k,j)           = banana_tau(j)/period
           tau(i,k,j+n_stack/2) = tau(i,k,j)

           theta_t(i,k,j)           = banana_theta(j)
           theta_t(i,k,j+n_stack/2) = theta_t(i,k,j)

        enddo

        ! omega(k) => coefficient of d/d_tau

        omega(i,k) = 1.0/period

     enddo

     do k=n_pass+1,n_lambda

        ! ** Trapped **

        class(k) = 2

        ck = class(k)

        sigma(ck,1) = -1.0
        sigma(ck,2) = 1.0

        do m=1,n_stack
           if (m < n_theta(ck)) then
              m_phys(ck,m) = m
              p_phys(ck,m) = 1
           else
              m_phys(ck,m) = n_stack-m+2
              p_phys(ck,m) = 2
           endif
        enddo

        do j=1,n_stack
           m_cyc(ck,j+n_stack,:) = j
           p_cyc(ck,:,j+n_stack,:) = 1.0
           m_cyc(ck,j,:) = j
           p_cyc(ck,:,j,:) = 1.0
           m_cyc(ck,j-n_stack,:) = j
           p_cyc(ck,:,j-n_stack,:) = 1.0
        enddo

        do j=1,n_theta(ck)
           m_map(ck,j,1) = j
           m_map(ck,j,2) = n_tau(ck)-j+2
        enddo
        m_map(ck,1,2) = 1

        call gyro_banana_uniform_taugrid(lambda(i,k),&
             n_tau(1)+1,&
             nint_ORB_do,&
             banana_theta,&
             banana_tau)

        ! Parameterization of entire orbit, given
        ! the trapped orbit computed by ORB library.

        ! Period here is really a half-period

        period = banana_tau(n_tau(1)+1)

        do j=1,n_tau(1)+1
           tau(i,k,j)     = banana_tau(j)/period
           theta_t(i,k,j) = banana_theta(j)
        enddo
        do j=2,n_tau(1)
           tau(i,k,j+n_tau(1)) = 2*tau(i,k,n_tau(1)+1)-tau(i,k,n_tau(1)+2-j)
           theta_t(i,k,j+n_tau(1)) = banana_theta(n_tau(1)+2-j)
        enddo

        ! omega(k) => coefficient of d/d_tau

        omega(i,k) = 1.0/period

     enddo

  enddo ! i

  !----------------------------------------------------
  ! Theta grid for solution of Maxwell equations:
  ! (evenly-spaced on [-pi,pi])
  !
  do j=1,n_theta_int
     theta_int(j) = -pi+(j-1)*pi_2/n_theta_int
  enddo

  if (n_theta_plot > 1) then

     ! Span [-pi,pi]

     do j=1,n_theta_plot
        theta_plot(j) = -pi+(j-1)*pi_2/n_theta_plot
     enddo

  else

     ! Special case for one point: 
     !  only theta=0 (outboard midplane)

     theta_plot(1) = 0.0

  endif

  !
  ! Print the grid to stdout:
  !
  if (i_proc == 0 .and. verbose_flag == 1) then

     i = ir_norm

     print *
     print *,' Interpolation Grid'
     print *,'-------------------'
     do j_int=1,n_theta_int
        print 10,j_int,theta_int(j_int)
     enddo

     print *
     print *,' Passing tau grid'
     print *,'------------------'
     do m=1,n_stack
        print 10,m,tau(i,1,m),theta_t(i,1:n_pass,m)
     enddo

     print *
     print *,' Trapped tau grid'
     print *,'------------------'
     do m=1,n_stack
        print 10,m,tau(i,n_lambda,m),theta_t(i,n_pass+1:n_lambda,m)
     enddo

     print *
     print *,' tau_hat grid'
     print *,'------------------'
     do k=1,n_lambda
        print 10,k,lambda(i,k),1/omega(i,k)
     enddo
  endif
  !--------------------------------------------------------

  !--------------------------------------------------------
  ! sigma-related operator for gyro_tau_derivative vectorization:
  !
  do ck=1,2
     do p=1,2
        if (ck == 2) then
           sigma_tau(ck,p) = 1.0
        else
           sigma_tau(ck,p) = -sigma(ck,p)
        endif
     enddo ! p
  enddo ! ck
  !--------------------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_banana_operators done]'
  endif

10 format(i2,2x,6(f12.9,1x))

end subroutine gyro_banana_operators
