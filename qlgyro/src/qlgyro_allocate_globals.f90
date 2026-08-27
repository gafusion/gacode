subroutine qlgyro_allocate_globals

  use qlgyro_globals

  implicit none

  ! Collision frequencies
  if (.not. allocated(nui)) allocate(nui(n_ion_max))
  if (.not. allocated(ti)) allocate(ti(n_ion_max))
  if (.not. allocated(dlntidr)) allocate(dlntidr(n_ion_max))
  if (.not. allocated(ni)) allocate(ni(n_ion_max))
  if (.not. allocated(dlnnidr)) allocate(dlnnidr(n_ion_max))
  
  pi = 4.0*atan(1.0)
  !
  mp  = 1

end subroutine qlgyro_allocate_globals
