!-----------------------------------------------------------------
! gyro_bessel_operator.f90
!
! PURPOSE:
! Compute all the Bessel-function kernel (averaging) stencils
! to the required bandwidth (m_gyro). Truncation is accomplished
! by assuming functions are constant in the truncated region.
!
! In the limit of maximal m_gyro (n_x/2), these are the
! pseudospectral operators.
!
! NOTES:
!
! itype=1: J_0(x)
! 2: J_0^2(x)
! 3: -(i/2)*k_x*rho*[ J_0(x)+J_2(x) ]
! 4: G_perp = (1/2)*[ J_0(x)+J_2(x) ]
! 5: G_perp^2(x)
! 6: G_perp(x)*J_0(x)
! 7: I_0(x^2)
! 8: i*k_x*rho*[ J_0(x)-J_1(x)/x ] / x^2
!
! x = rho*sqrt((2*pi*p*a+u)^2+v^2)
!
! kx -> pi_2*p*a + u
! ky -> v
!
!
! WARNING:
!  RJBESL requires non-negative argument, so we need to operate on 
!  argument with abs().
!
! INPUT:
!
! a = |grad(r)|/L
! u = k_theta*Gq*Theta
! v = k_theta*Gq
!
! rho = v_perp/Omega_c (itype=1-6,8)
! rho = rhos_unit*sqrt(T)/(mu*z*b0) (itype=7)
!------------------------------------------------------------------

subroutine gyro_bessel_operator(rho,a,u,v,g,itype)

  use gyro_globals, only : z_gyro, i_gyro, m_gyro, n_x, truncation_method
  use math_constants

  !-----------------------------------------
  implicit none
  !
  real, intent(in) :: a
  real, intent(in) :: u
  real, intent(in) :: v
  real, intent(in) :: rho
  integer, intent(in) :: itype
  complex, intent(inout) :: g(-m_gyro:m_gyro-i_gyro)
  !
  integer :: p
  integer :: p0
  integer :: m
  !
  real :: x
  real :: func(-n_x/2:n_x/2-1)
  complex :: g0
  complex :: gp
  !
  real, external :: besei0
  !-----------------------------------------

  p0 = n_x/2

  select case (itype)

  case (1)

     ! J_0
     do p=-p0,p0-1
        x = abs(rho)*sqrt((pi_2*p*a+u)**2+v**2)
        func(p) = bessel_j0(x)/n_x
     enddo

  case (2)

     ! J_0^2
     do p=-p0,p0-1
        x = abs(rho)*sqrt((pi_2*p*a+u)**2+v**2)
        func(p) = bessel_j0(x)**2/n_x
     enddo

  case (3)

     ! -(i/2)*k_x*rho*[ J_0(z)+J_2(z) ]

     ! The factor -(i/2) will be applied outside this loop
     do p=-p0,p0-1
        x = abs(rho)*sqrt((pi_2*p*a+u)**2+v**2)
        func(p) = (pi_2*p*a+u)*rho*(bessel_j0(x)+bessel_jn(2,x))/n_x
     enddo

  case (4)

     ! G = (1/2)*[ J_0(z)+J_2(z) ]
     do p=-p0,p0-1
        x = abs(rho)*sqrt((pi_2*p*a+u)**2+v**2)
        func(p) = 0.5*(bessel_j0(x)+bessel_jn(2,x))/n_x
     enddo

  case (5)

     ! G^2
     do p=-p0,p0-1
        x = abs(rho)*sqrt((pi_2*p*a+u)**2+v**2)
        func(p) = (0.5*(bessel_j0(x)+bessel_jn(2,x)))**2/n_x
     enddo

  case (6)

     ! G * J_0
     do p=-p0,p0-1
        x = abs(rho)*sqrt((pi_2*p*a+u)**2+v**2)
        func(p) = 0.5*(bessel_j0(x)+bessel_jn(2,x))*bessel_j0(x)/n_x
     enddo

  case (7)

     ! I_0
     do p=-p0,p0-1
        x = abs(rho)*sqrt((pi_2*p*a+u)**2+v**2)
        func(p) = besei0(x*x)/n_x
     enddo

  case (8)

     ! i*k_x*rho*[ J_0(x)-J_1(x)/x ] / x^2

     ! The factor i will be applied outside this loop
     if (u**2+v**2 > 1e-8) then
        do p=-p0,p0-1
           x = abs(rho)*sqrt((pi_2*p*a+u)**2+v**2)
           func(p) = (pi_2*p*a+u)*rho*(bessel_j0(x)-bessel_j1(x)/x)/x**2/n_x
        enddo
     else
        ! Case is n=0, which will ultimately give zero contribution so we 
        ! zero it here explicitly to avoid 0/0.
        func(:) = 0.0
     endif

  end select

  !-----------------------------------------------------------
  ! Real-space operators:
  !
  ! Construct real-space forms by summing over Fourier modes.
  ! Real-space gridpoints in truncated region are ignored.
  !
  do m=-m_gyro,m_gyro-i_gyro
     g0 = (0.0,0.0)
     do p=-p0,p0-1
        g0 = g0 + z_gyro(m,p)*func(p)
     enddo
     g(m) = g0
  enddo
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  ! Final factors
  !
  ! Multiply by -i/2 for case 3
  if (itype == 3) g = -(i_c/2.0)*g

  ! Multiply by i for case 8
  if (itype == 8) g = i_c*g
  !-----------------------------------------------------------


  if (truncation_method == 1) then

     ! ORIGINAL METHOD

     ! Enforce reality if n=0:
     if (abs(u) <= epsilon(0.0)) then
        if (i_gyro /= 1) then

           ! Correct truncated gyroaverage
           g0 = sum(g(:))-1.0
           g(0) = g(0)-g0

        endif
        g = real(g)
     endif

  else

     ! ALTERNATIVE METHOD

     !-----------------------------------------------------------
     ! Make end-coefficients equal to sum of neglected ones.
     ! This is equivalent to assuming that functions are
     ! constant over the truncated region.
     !
     if (m_gyro < p0) then
        do m=-p0,-m_gyro-1
           g0 = (0.0,0.0)
           do p=-p0,p0-1
              g0 = g0 + z_gyro(m,p)*func(p)
           enddo
           ! Periodic point (give 1/2 this bit to other end)
           if (m == -p0) then
              gp = g0/2.0
              g(-m_gyro) = g(-m_gyro)+gp
           else
              g(-m_gyro) = g(-m_gyro)+g0
           endif
        enddo
        do m=m_gyro+1,p0-1
           g0 = (0.0,0.0)
           do p=-p0,p0-1
              g0 = g0 + z_gyro(m,p)*func(p)
           enddo
           g(m_gyro) = g(m_gyro)+g0
        enddo
        ! Ficticious m=p0 gets half the contribution from m=-p0:
        g(m_gyro) = g(m_gyro)+gp
     endif
     !-----------------------------------------------------------

     ! Enforce reality if n=0:
     if (abs(u) <= epsilon(0.0)) g = real(g)

  endif

end subroutine gyro_bessel_operator

