!---------------------------------------------------------
! cgyro_shear.f90
!
! PURPOSE:
!  Manage ExB shear algorithms.
!---------------------------------------------------------

subroutine cgyro_shear

  use cgyro_globals
  use timer_lib

  implicit none

  ! Spectral ExB shear
  call timer_lib_in('shear')
  select case(shear_method)
  case (1)
     ! Discrete shift (Hammett) 
     call cgyro_shear_hammett
  case (3)
     ! Linear shift (more accurate than discrete shift)
     call cgyro_shear_linear
  case (4)
     ! Linear shift (more accurate than discrete shift)
     call cgyro_shear_crank
  end select
  call timer_lib_out('shear')

end subroutine cgyro_shear
