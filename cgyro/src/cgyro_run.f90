!---------------------------------------------------------------
! cgyro_run.f90
!
! PURPOSE:
!  This is the actual system call to run CGYRO via subroutine 
!  interface (i.e, from TGYRO).
!
! NOTES:
! - To run CGYRO as a subroutine, do this:
!   call cgyro_init
!   call cgyro_init_kernel
!   call cgyro_run
!   call cgyro_final_kernel
!-------------------------------------------------------------

subroutine cgyro_run(test_flag_in,var_in,n_species_out,flux_tave_out,tave_min_out,tave_max_out,status_out)

  use cgyro_globals, only: test_flag,max_time,error_status, &
                           nonlinear_flag,signal,n_species
  use cgyro_flux_mod, only: cgyro_flux_sum

  implicit none

  ! Input parameters (IN) - REQUIRED
  integer, intent(in)   :: test_flag_in
  integer, intent(in)   :: var_in           
  integer, intent(out)  :: n_species_out
  real, intent(out)     :: tave_min_out, tave_max_out
  real, intent(out)     :: flux_tave_out(11,3)
  integer, intent (out) :: status_out  ! 0 is good to continue; > 0 means stop

  ! Set corresponding global variables
  test_flag = test_flag_in

  ! Re-set max time
  if (var_in > 0.0) then
     max_time = var_in
  endif
  
  ! Run GYRO
  call cgyro_kernel

  if (error_status == 0) then
     if (nonlinear_flag == 0 .and. signal > 0) then
        ! linear converged or underflow
        status_out = 1
     endif
     ! no errors
     status_out = 0
  else
     ! something wrong
     status_out = 2
  endif


  n_species_out = n_species
  call cgyro_flux_sum(tave_min_out, tave_max_out, flux_tave_out)
  
end subroutine cgyro_run
