!---------------------------------------------------------
! prgen_read_gato.f90
!
! PURPOSE:
!  Read data from GATO (eqgta, grid.dat, o1gta):
!
!  - psi (poloidal flux over 2pi)
!  - V(psi) (volume)
!  - dV(psi)/dpsi
!  - R(psi,arclength)
!  - Z(psi,arclength)
!  - q
!----------------------------------------------------------

subroutine prgen_read_gato

  use prgen_globals

  implicit none

  integer :: i
  integer :: ip
  real :: fa,fb
  real, dimension(:), allocatable :: gato_psi  
  real, dimension(:), allocatable :: gato_dummy
  real, dimension(:), allocatable :: gato_q
  real, dimension(:), allocatable :: gato_p
  real, dimension(:,:), allocatable :: gato_bigr
  real, dimension(:,:), allocatable :: gato_bigz
  real, dimension(:,:), allocatable :: gvec
  real, dimension(:,:,:), allocatable :: g3vec
  real, dimension(:,:,:), allocatable :: g3rho

  character (len=70) :: cdum

  !--------------------------------------------------
  ! Read gfile header info
  !
  open(unit=1,file='eqgta',status='old')
  read(1,'(a)') efit_header
  close(1) 
  !--------------------------------------------------
  
  !---------------------------------------------------
  ! Read flux-surface data in GATO's grid.dat file
  !
  ! nsurf = number of finite flux surfaces
  ! 
  open(unit=1,file='out.gato.fluxsurf',status='old')
  read(1,'(a)') cdum
  read(1,*) nsurf,narc

  allocate(gato_psi(0:nsurf))
  allocate(gato_dummy(0:nsurf))
  allocate(gato_q(0:nsurf))
  allocate(gato_p(0:nsurf))
  allocate(gato_bigr(narc,nsurf))
  allocate(gato_bigz(narc,nsurf))

  read(1,'(a)') cdum

  ! psi mesh
  read(1,'(a)') cdum
  read(1,*) gato_psi(:) 

  ! Volume
  read(1,'(a)') cdum
  read(1,*) gato_dummy(:)

  ! dVol/dpsi
  read(1,'(a)') cdum
  read(1,*) gato_dummy(:)

  ! Pressure
  read(1,'(a)') cdum
  read(1,*) gato_p(:)

  ! Safety Factor
  read(1,'(a)') cdum
  read(1,*) gato_q(:)

  ! 2D (R,Z)
  read(1,'(a)') cdum
  read(1,*) gato_bigr(:,:)
  read(1,'(a)') cdum
  read(1,*) gato_bigz(:,:)

  close(1)
  !----------------------------------------------------

  !----------------------------------------------------------------
  ! Correct flux variation in profile data (ITERDB, etc) 
  !
  if (format_type == 0 .or. format_type == 7) then
     ! Case 1: raw gfile mode ; dpsi is undefined
     dpsi_data = gato_psi(nsurf)
     dpsi_efit = gato_psi(nsurf)
     do i=1,nx
        ! Shrink slightly for safety
        dpsi(i) = (i-1)*dpsi_efit/(nx-1+1e-12)
     enddo
  else
     ! Case 2: typical case 
     dpsi_data = dpsi(nx)
     dpsi_efit = gato_psi(nsurf)
     
     ! Ensure max(dpsi) = max(gato_psi) 
     dpsi(:)  = dpsi(:)*dpsi_efit/dpsi_data
     ! Extra insurance against roundoff
     dpsi(nx) = dpsi_efit
  endif
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Get Miller-style (model shape) coefficients using 
  ! fluxfit routines:
  !
  ! 1  r (cm)
  ! 2  z_mag (cm)
  ! 3  R0 (cm)
  ! 4  kappa
  ! 5  delta
  ! 6  zeta
  !
  call fluxfit_driver(1,1,nsurf,narc,gato_bigr,gato_bigz,verbose_flag)
  !
  allocate(gvec(6,0:nsurf))

  open(unit=1,file='fluxfit.profile',status='old')
  read(1,*) cdum
  read(1,*) gvec(:,1:nsurf)
  close(1)

  ! Explicitly set rmin=0 at origin
  gvec(1,0) = 0.0

  ! Use extrapolation to get values of shape parameters at origin
  do i=2,6
     call bound_extrap(fa,fb,gvec(i,:),gato_psi,nsurf+1)
     gvec(i,0) = fa
  enddo

  ! Map shape coefficients onto poloidal flux (dpsi) grid:
  ! NOTE: need sqrt here to get sensible behaviour as r -> 0.
  call cub_spline(sqrt(gato_psi),gvec(1,:),nsurf+1,sqrt(dpsi),rmin,nx)
  call cub_spline(gato_psi,gvec(2,:),nsurf+1,dpsi,zmag,nx)
  call cub_spline(gato_psi,gvec(3,:),nsurf+1,dpsi,rmaj,nx)
  call cub_spline(gato_psi,gvec(4,:),nsurf+1,dpsi,kappa,nx)
  call cub_spline(gato_psi,gvec(5,:),nsurf+1,dpsi,delta,nx)
  call cub_spline(gato_psi,gvec(6,:),nsurf+1,dpsi,zeta,nx)

  ! Total pressure and q from GATO-EFIT 
  ! (possible overwrite from raw data input in statefile)
  if (nop_flag == 0) then
     print '(a)','INFO: (prgen) Using total pressure from gfile.'
     call cub_spline(gato_psi,gato_p,nsurf+1,dpsi,p_tot,nx)
  endif
  if (noq_flag == 0) then
     print '(a)','INFO: (prgen) Using safety factor (q) from gfile.'
     call cub_spline(gato_psi,gato_q,nsurf+1,dpsi,q,nx)
  endif

  ! Explicitly set rmin=0 at origin
  rmin(1) = 0.0

  deallocate(gvec)
  !----------------------------------------------------------------

  !----------------------------------------------------
  ! Get general geometry coefficients
  !
  call fluxfit_driver(2,nfourier,nsurf,narc,gato_bigr,gato_bigz,verbose_flag)

  allocate(g3vec(4,0:nfourier,0:nsurf))
  allocate(g3rho(4,0:nfourier,nx))

  open(unit=1,file='fluxfit.geo',status='old')
  read(1,*) ip
  read(1,*) g3vec(:,:,1:nsurf)
  close(1)

  ! Explicitly set rmin=0 at origin
  g3vec(:,1:nfourier,0) = 0.0

  ! Extrapolate centers (R0,Z0) to origin
  do i=1,4
     call bound_extrap(fa,fb,g3vec(i,0,:),gato_psi,nsurf+1)
     g3vec(i,0,0) = fa
  enddo

  ! Map results onto poloidal flux (dpsi) grid:
  ! NOTE: need sqrt here to get sensible behaviour as r -> 0.
  do i=1,4
     do ip=0,nfourier
        call cub_spline(sqrt(gato_psi),g3vec(i,ip,:),nsurf+1,sqrt(dpsi),g3rho(i,ip,:),nx)
     enddo
  enddo

  ! Explicitly set rmin=0 at origin
  g3rho(:,1:nfourier,1) = 0.0

  open(unit=1,file='input.profiles.geo',status='replace')
  write(1,'(a)') '# input.profiles.geo'
  write(1,'(a)') '#'
  write(1,'(a)') '# See https://fusion.gat.com/theory/input.profiles.geo for complete documentation.'
  write(1,'(a)') '#'
  write(1,'(a)') '# File format:'
  write(1,'(a)') '#-------------------'
  write(1,'(a)') '# nfourier'
  write(1,'(a)') '# a[4,0:nfourier,nx]'  
  write(1,'(a)') '#-------------------'
  write(1,'(a)') '#'
  write(1,'(a)') '# NOTE: nx=EXPRO_n_exp is defined in input.profiles'
  write(1,*) nfourier
  write(1,'(1pe20.13)') g3rho(:,:,:)
  close(1)
  !----------------------------------------------------

  ! Clean up
  deallocate(g3vec)
  deallocate(g3rho)
  deallocate(gato_psi)
  deallocate(gato_q)
  deallocate(gato_p)
  deallocate(gato_dummy)
  deallocate(gato_bigr)
  deallocate(gato_bigz)

end subroutine prgen_read_gato
