!-----------------------------------------------------
! gyro_banana_init.f90
!
! PURPOSE:
!  Simple routine to initialize internal variables.
!-----------------------------------------------------

subroutine gyro_banana_init(n_IN)

  use gyro_banana_private
  use geo

  !-----------------------------------------------------------
  implicit none
  !
  integer :: i
  integer, intent(in) :: n_IN 
  !
  real, dimension(:), allocatable :: f
  real :: d_theta
  !-----------------------------------------------------------

  n = n_IN

  !----------------------------------------------
  ! Determine maximum lambda, and lambda at
  ! the trapped-passing boundary:
  !
  allocate(ttmp(1))
  ttmp(1) = 0.0
  call geo_interp(1,ttmp,.true.)
  !
  lambda_max = 1.0/GEO_b(1)
  !
  ttmp(1) = pi
  call geo_interp(1,ttmp,.false.)
  !
  lambda_tp = 1.0/GEO_b(1)
  deallocate(ttmp)
  !-----------------------------------------------
  
  allocate(f(n))
  allocate(ttmp(n))

  d_theta = 2.0*pi/(n-1)

  do i=1,n
     ttmp(i) = -pi+(i-1)*d_theta
  enddo
  call geo_interp(n,ttmp,.false.)
  
  f(:) = GEO_g_theta(:)/GEO_b(:)

  fluxave = (0.5*(f(1)+f(n))+sum(f(2:n-1)))*d_theta

  deallocate(f,ttmp)

end subroutine gyro_banana_init
