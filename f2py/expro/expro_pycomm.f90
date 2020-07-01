! Simplified (non-MPI) version of expro_comm routines suitable for f2py

subroutine expro_icomm(p)

  implicit none

  integer, intent(inout) :: p

  read(1,*) p

end subroutine expro_icomm

subroutine expro_rcomm(x)

  implicit none

  double precision, intent(inout) :: x

  read(1,10) x

10 format(1pe14.7)

end subroutine expro_rcomm

subroutine expro_lcomm(x,n)
  
  implicit none

  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x

  read(1,10) x

10 format(10(1pe14.7))

end subroutine expro_lcomm

subroutine expro_scomm(x,n)
  
  implicit none

  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x

  read(1,*) x

end subroutine expro_scomm

subroutine expro_vcomm(x,n)

  implicit none

  integer :: idum,i
  integer, intent(in) :: n
  double precision, intent(inout), dimension(n) :: x

  do i=1,n
     read(1,10) idum,x(i)
  enddo

10 format(i3,1x,1pe14.7)

end subroutine expro_vcomm

subroutine expro_tcomm(x,n)

  implicit none

  integer, intent(in) :: n
  character*10, intent(inout) :: x(n)
  
  read(1,*) x

end subroutine expro_tcomm

subroutine expro_acomm(x,m,n)

  implicit none

  integer :: idum,i
  integer, intent(in) :: m,n
  double precision, intent(inout), dimension(m,n) :: x

  do i=1,n
     read(1,10) idum,x(:,i)
  enddo

10 format(i3,1x,10(1pe14.7,1x))

end subroutine expro_acomm


