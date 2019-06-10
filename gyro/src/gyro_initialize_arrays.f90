!----------------------------------------------------
! gyro_initialize_arrays.f90
!
! PURPOSE:
!  Post-allocation initialization of arrays.
!------------------------------------------------

subroutine gyro_initialize_arrays

  use gyro_globals

  !---------------------------------------------------
  implicit none
  !---------------------------------------------------

  !---------------------------------------------------
  h     = 0.0
  h_old = 0.0
  h_0   = 0.0
  h_err = 0.0
  !---------------------------------------------------

  !---------------------------------------------------
  omega_linear = 0.0
  time_error   = 0.0

  ! No field_blend
  field_blend_dot  = 0.0
  field_blend_old  = 0.0
  field_blend_old2 = 0.0

  ! No field_tau
  field_tau_old  = 0.0
  field_tau_old2 = 0.0

  ! No cap_h 
  h_cap      = 0.0
  h_cap_dot  = 0.0
  h_cap_old  = 0.0
  h_cap_old2 = 0.0

  gyro_uv_old  = 0.0
  gyro_uv_old2 = 0.0

  entropy(:,:) = 0.0
  !---------------------------------------------------

  !---------------------------------------------------
  ! Initialization to prevent MPI-grid errors
  if (collision_flag == 1) then
     f_coll  = 0.0
     fb_coll = 0.0
  endif
  !---------------------------------------------------

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_initialize_arrays done]' 
  endif

end subroutine gyro_initialize_arrays
