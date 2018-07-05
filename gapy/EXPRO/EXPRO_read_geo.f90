!--------------------------------------------------------------
! EXPRO_read_geo.f90
!
! PURPOSE:
!  Read EXPRO_geo from input.profiles.geo.
!--------------------------------------------------------------

subroutine EXPRO_read_geo

  use EXPRO_globals
  use EXPRO_interface

  implicit none

  integer, parameter :: io=1
  integer :: i

  open(unit=io,&
       file=trim(path)//'input.profiles.geo',&
       status='old')

  call EXPRO_skip_header(io)
  read(io,*) EXPRO_nfourier
  do i=1,EXPRO_n_exp
     read(io,*) EXPRO_geo(:,:,i)
  enddo

  close(io)

end subroutine EXPRO_read_geo
