!--------------------------------------------------------------
! EXPRO_read.f90
!
! PURPOSE:
!  Read experimental profiles (serial).
!
! NOTES:
!  Variable naming corresponds to input.profiles 
!  documentation at 
!
!    http://fusion.gat.com/theory/input.profiles
!--------------------------------------------------------------

subroutine EXPRO_read

  use EXPRO_interface

  implicit none

  call EXPRO_read_driver
  call EXPRO_compute_derived

end subroutine EXPRO_read

