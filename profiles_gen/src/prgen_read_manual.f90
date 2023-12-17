!--------------------------------------------------------------
! prgen_read_genf.f90
!
! PURPOSE:
!  Read profiles from General Fusion Statefile
!--------------------------------------------------------------

subroutine prgen_read_manual

  use prgen_globals

  implicit none

  integer :: i,n0
  
  character(len=1) :: a
  real, dimension(:,:), allocatable :: x
  real, dimension(:), allocatable :: rho0

  ! Number of resampled points
  nx = 16
  call prgen_allocate('')
  
  call getlen('te.txt',n0)
  
  allocate(x(2,n0))
  allocate(rho0(n0))

  ! te
  open(unit=1,file='te.txt',status='old')
  read(1,*) x(:,:)
  close(1)
  rho0 = x(1,:)
  do i=1,nx
     rho(i) = rho0(1)+(i-1.0)/(nx-1.0)*(rho0(n0)-rho0(1))
  enddo
  call bound_interp(rho0,x(2,:)*1e-3,n0,rho,te_kev,nx)

  ! ti
  open(unit=1,file='ti.txt',status='old')
  read(1,*) x(:,:)
  close(1)
  call bound_interp(rho0,x(2,:)*1e-3,n0,rho,ti_kev,nx)

  ! ne
  open(unit=1,file='ne.txt',status='old')
  read(1,*) x(:,:)
  close(1)
  call bound_interp(rho0,x(2,:)*1e-19,n0,rho,ne_e19m3,nx)

  ! ni
  open(unit=1,file='ni.txt',status='old')
  read(1,*) x(:,:)
  close(1)
  call bound_interp(rho0,x(2,:)*1e-19,n0,rho,ni_e19m3,nx)

end subroutine prgen_read_manual

subroutine getlen(infile,n)

  implicit none
  
  character(len=*), intent(in) :: infile
  integer, intent(inout) :: n
  character(len=1) :: a
  integer :: ios
  
  open(unit=1,file=trim(infile),status='old')
  n = -1
  do
     read(1,*,iostat=ios) a
     if (ios < 0) exit
     n = n+1
  enddo
  close(1)
  
end subroutine getlen

