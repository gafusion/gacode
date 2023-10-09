!--------------------------------------------------------
! Read optional Ti file
!
! col 1: rho (root of normalized toroidal flux)
! col 2: Ti (keV)
!--------------------------------------------------------

subroutine prgen_read_ti

  use prgen_globals

  implicit none

  ! number of columns in Ti file

  character (len=1) :: a
  integer :: n_in
  integer :: i
  integer :: ierr
  integer :: count
  real :: xa,xb
  real :: ta,tb
  real, dimension(2) :: x
  real, dimension(:), allocatable :: rho_in
  real, dimension(:), allocatable :: t_in

  !-------------------------------------------------------
  ! Determine vector length (assume no header)
  !
  count = 0 
  open(unit=1,file=file_ti,status='old')
  do 
     read(1,*,iostat=ierr) x(1)
     if (ierr < 0) exit
     count = count+1
  enddo
  close(1)
  !
  n_in = count
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Read data (assume no header)
  !
  allocate(rho_in(n_in))
  allocate(t_in(n_in))
  !
  open(unit=1,file=file_ti,status='old')
  do i=1,n_in 
     read(1,*,iostat=ierr) x(:)
     if (ierr < 0) then
        print '(t2,a,a)','ERROR: Ti format error.'
        stop
     endif
     rho_in(i) = x(1)
     t_in(i) = x(2) ! Ti ()
  enddo
  print '(a,a)','INFO: Assuming 2-column format (root_torflux,Ti)',trim(file_ti)

  close(1)
  !-------------------------------------------------------

  print *,rho_in(:)
  print *,rho
  !-------------------------------------------------------
  ! Spline interpolation from CER grid to ITERDB rho grid.
  !
  call cub_spline(rho_in,t_in,n_in,rho,ti_kev,nx)
  !-------------------------------------------------------

end subroutine prgen_read_ti
