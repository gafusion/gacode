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
!   call cgyro_run
!-------------------------------------------------------------

subroutine cgyro_run(test_flag_in)

  use mpi
  use cgyro_globals
  use timer_lib

  implicit none

  ! Input parameters (IN) - REQUIRED
  integer, intent(in) :: test_flag_in

  ! Set corresponding global variables
  test_flag = test_flag_in

  ! Run GYRO
  call cgyro_kernel

end subroutine cgyro_run
