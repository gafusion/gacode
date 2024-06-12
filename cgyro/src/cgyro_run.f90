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

  use mpi
  use cgyro_globals
  use timer_lib

  implicit none

  ! Input parameters (IN) - REQUIRED
  integer, intent(in)   :: test_flag_in
  integer, intent(in)   :: var_in           
  integer, intent(out)  :: n_species_out
  real, intent(out)     :: tave_min_out, tave_max_out
  real, intent(out)     :: flux_tave_out(11,3)
  integer, intent (out) :: status_out  ! 0 is good to continue; > 0 means stop
  real, dimension(:,:), allocatable :: sum_out

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


  ! Return time-averaged flux data (need to reduce across n first)
  flux_tave_out(:,:) = 0.0
  allocate(sum_out(n_species,3))
  if (abs(gamma_e) > 1e-10) then
     call MPI_ALLREDUCE(cflux_tave(:,:), &
       sum_out(:,:), &
       size(sum_out), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)
  else
     call MPI_ALLREDUCE(gflux_tave(:,:), &
       sum_out(:,:), &
       size(sum_out), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)
  endif
  
  n_species_out = n_species
  tave_min_out = tave_min
  tave_max_out = tave_max
  flux_tave_out(1:n_species,:) = sum_out(1:n_species,:)/tave_step

  deallocate(sum_out)
  
end subroutine cgyro_run
