subroutine create_write

  use create_globals
  use EXPRO_interface

  implicit none

  open(unit=1,file='input.profiles.create',status='replace')

  !---------------------------------------------------------------
  ! Basic information
  !
  write(1,20) '# * This file was generated with create'

  write(1,20) '#'
  write(1,20) '# See https://fusion.gat.com/theory/input.profiles for complete documentation.'
  write(1,20) '#'
  write(1,20) '# * Data description:'
  write(1,20) '#'
  write(1,30) '#    RADIAL GRIDPOINTS : ',EXPRO_n_exp
  !---------------------------------------------------------------

  write(1,20) '# '
  write(1,24) 'N_ION=',EXPRO_n_ion
  if (EXPRO_n_exp > 99) then
     write(1,30) 'N_EXP=',EXPRO_n_exp
  else
     write(1,25) 'N_EXP=',EXPRO_n_exp
  endif
  write(1,'(a,sp,1pe14.7)') 'BT_EXP=',EXPRO_b_ref
  write(1,60) 'ARHO_EXP=',EXPRO_arho

  call EXPRO_write(1)

  close(1)

20 format(30(a))
24 format(a,i1)
25 format(a,i2)
30 format(a,i3)
60 format(a,1pe13.7)

end subroutine create_write
