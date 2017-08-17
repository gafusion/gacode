!--------------------------------------------------------
! Read optional CER file
!
! col 1: rho (root of normalized toroidal flux)
! col 2: Ti (keV)
! col 3: nc (10^20/m^3) 
! col 4: vpolc (km/s)
! col 5: vtorc (km/s)
! col 6: Er (kV/m)
! col 7: omega0 (krad/s) [(c*Er)/(R*Bp)]
!
! NOTE also that the CER toroidal angle sign convention 
! is fixed and opposite the GACODE convention.  The CER 
! angle is phi in the (R,phi,Z) system defined here 
!
! https://fusion.gat.com/theory/Gyrofieldorient 
!--------------------------------------------------------

subroutine prgen_read_cer

  use prgen_globals

  implicit none

  ! number of columns in CER file
  integer, parameter :: n_col = 7

  character (len=1) :: a
  integer :: n_in
  integer :: i
  integer :: ierr
  integer :: count
  real :: xa,xb
  real :: fa(3),fb(3)
  real, dimension(n_col) :: x
  real, dimension(:), allocatable :: rho_in
  real, dimension(:,:), allocatable :: f_in

  !-------------------------------------------------------
  ! Determine vector length
  !
  count = 0 
  open(unit=1,file=cer_file,status='old')
  read(1,*) a
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
  ! Read data
  !
  allocate(rho_in(n_in))
  allocate(f_in(n_in,3))
  !
  open(unit=1,file=cer_file,status='old')
  read(1,*) a
  do i=1,n_in 
     read(1,*,iostat=ierr) x(:)
     if (ierr < 0) then
        print '(t2,a,a)','ERROR: CER format error.'
        stop
     endif
     rho_in(i) = x(1)
     f_in(i,1) = x(4) ! vpol (km/s)
     f_in(i,2) = x(5) ! vtor (km/s)
     f_in(i,3) = x(7) ! omega0 (krad/s)
  enddo
  print '(a,a)','INFO: Assuming 7-column Solomon CER format in ',trim(cer_file)

  close(1)
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Map last gridpoint so that rho_in(n_er) = 1.
  !
  xa = rho_in(n_in-1)
  xb = rho_in(n_in)
  fa(:) = f_in(n_in-1,:)
  fb(:) = f_in(n_in,:)
  !
  rho_in(n_in) = 1.0
  f_in(n_in,:) = fa(:)*(xb-1.0)/(xb-xa)+fb(:)*(1.0-xa)/(xb-xa)
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Spline interpolation from CER grid to ITERDB rho grid.
  !
  ! COORDINATES: convert vpolc to m/s (- sign for CER phi to GACODE varphi)
  call cub_spline(rho_in,-f_in(:,1)*1e3,n_in,rho,vpolc_exp,nx)
  !
  ! COORDINATES: convert vtorc to m/s (- sign for CER phi to GACODE varphi)
  call cub_spline(rho_in,-f_in(:,2)*1e3,n_in,rho,vtorc_exp,nx)
  !
  ! COORDINATES: convert omega to (1/s) (- sign for CER phi to GACODE varphi)
  call cub_spline(rho_in,-f_in(:,3)*1e3,n_in,rho,omega0,nx)
  !-------------------------------------------------------

end subroutine prgen_read_cer
