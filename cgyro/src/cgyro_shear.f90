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
  if (shear_method == 1) then
     ! Discrete shift (Hammett) 
     call cgyro_shear_hammett
  endif
  call timer_lib_out('shear')

end subroutine cgyro_shear
