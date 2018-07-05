subroutine EXPRO_alloc(path_in,flag)

  use EXPRO_globals
  use EXPRO_interface

  implicit none

  character(len=*), intent(in) :: path_in
  integer, intent(in) :: flag

  path = path_in
  comm = -1

  call EXPRO_alloc_control(0,flag)

end subroutine EXPRO_alloc
