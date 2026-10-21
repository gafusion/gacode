! Print NEO's F(m,n,lambda) and FB(m,n,lambda) from neo_compute_fcoll.
program fcoll_dump
  implicit none
  integer :: m0, m, n
  real :: lambda
  real, allocatable :: f(:,:), fb(:,:)
  character(len=64) :: arg
  call get_command_argument(1,arg); read(arg,*) m0
  call get_command_argument(2,arg); read(arg,*) lambda
  allocate(f(-m0:m0,-m0:m0), fb(-m0:m0,-m0:m0))
  call neo_compute_fcoll(m0,lambda,f,fb)
  do m=-m0,m0
     do n=-m0,m0
        if (m+n >= -1) write(*,'(2i5,2es45.35)') m, n, f(m,n), fb(m,n)
     enddo
  enddo
end program fcoll_dump
