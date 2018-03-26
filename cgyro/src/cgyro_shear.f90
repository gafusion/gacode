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
  if (shear_method == 1) then
     ! Discrete shift (Hammett) 
     call timer_lib_in('shear')
     call cgyro_shear_hammett
     call timer_lib_out('shear')
  endif

end subroutine cgyro_shear
