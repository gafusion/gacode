!-----------------------------------------------------------------
! cgyro_step.f90
!
! PURPOSE:
!  Common work variables for adaptive RK routines
!-----------------------------------------------------------------

module cgyro_step
  
  integer :: itrk,conv

  real :: tol
  real :: delta_t_last,delta_t_last_step
  real :: deltah2
  real :: deltah2_max,deltah2_min
  real :: delta_x,delta_x_min,delta_x_max
  real :: local_max_error
  real :: rel_error,var_error,scale_x
  real :: delta_t_tot

  !---------------------------------
  ! GPU-only
  real :: error_rhs,error_hx
  !---------------------------------

  real, dimension(2) :: error_sum,error_x
  
  integer, parameter :: itrk_max = 1000
  real, parameter :: eps = 2e-12
  
end module cgyro_step

  
