!-----------------------------------------------------------
! tgyro_multi_driver.f90
!
! PURPOSE:
!  Main driver for multi-job utility.
!----------------------------------------------------------

subroutine tgyro_multi_driver

  use mpi
  use tgyro_globals
  use gyro_interface

  implicit none

  gyro_restart_method = 1

  ! See gyro/src/gyro_globals.f90 for definition of transport_method
  transport_method = 1

  if (lcode == 'gyro') then

     ! Initialize GYRO
     call gyro_init(lpath,gyro_comm)

     ! Run GYRO
     call gyro_run(gyrotest_flag,gyro_restart_method,transport_method)

  else 

     ! Initialize CGYRO
     call cgyro_init(lpath,gyro_comm)

     ! Run CGYRO
     call cgyro_run(gyrotest_flag)

  endif

  ! Artificially trigger error to print message
  gyro_error_status_out = 1
  call tgyro_trap_component_error(gyro_error_status_out,gyro_error_message_out)

end subroutine tgyro_multi_driver

