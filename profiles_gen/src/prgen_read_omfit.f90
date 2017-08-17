!----------------------------------------------------------------
! prgen_read_omfit.f90
!
! PURPOSE:
!  Read data from OMFIT flux-surface mapper and fit to 
!  geometry coefficients.
!----------------------------------------------------------------

subroutine prgen_read_omfit

  use prgen_globals

  implicit none

  integer :: i,i1,i2
  integer :: ip
  integer :: ntot
  real :: fa,fb
  character (len=1) :: a
  integer, dimension(:), allocatable :: narcv
  real, dimension(:), allocatable :: r_raw,z_raw,l_raw
  real, dimension(:,:), allocatable :: r2,z2
  real, dimension(:), allocatable :: sqdpsi
  real, dimension(:), allocatable :: psi
  real, dimension(:), allocatable :: q_omfit
  real, dimension(:), allocatable :: p_omfit
  real, dimension(:,:), allocatable :: gvec
  real, dimension(:,:,:), allocatable :: g3vec
  real, dimension(:,:,:), allocatable :: g3rho

  !----------------------------------------------------------------
  ! Read gfile header info
  !
  open(unit=1,file='gfile',status='old')
  read(1,'(a)') efit_header
  close(1) 
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Read the OMFIT mapper file
  open(unit=1,file='out.omfit.fluxsurf',status='old')

  do i=1,8
     read(1,*) a
  enddo
  read(1,*) nsurf

  allocate(narcv(nsurf))
  read(1,*) narcv(:)

  ntot = sum(narcv)
  narc = narcv(1)

  allocate(psi(nsurf))
  allocate(q_omfit(nsurf))
  allocate(p_omfit(nsurf))
  allocate(r_raw(ntot))
  allocate(z_raw(ntot))
  allocate(l_raw(ntot))

  read(1,*) psi(:)
  read(1,*) q_omfit(:)
  read(1,*) p_omfit(:)
  read(1,*) r_raw(:)
  read(1,*) z_raw(:)
  read(1,*) l_raw(:)
  close(1)

  psi(:) = psi(:)-psi(1)

  allocate(r2(narc,nsurf))
  allocate(z2(narc,nsurf))

  ! Repack data into 2D array

  do i=1,nsurf
     if (i == 1) then
        i1 = 1
     else
        i1 = sum(narcv(1:(i-1)))+1
     endif
     i2 = sum(narcv(1:i))
     ! Reverse order
     r2(:,i) = r_raw(i2:i1:-1)
     z2(:,i) = z_raw(i2:i1:-1)
  enddo

  deallocate(r_raw)
  deallocate(z_raw)
  deallocate(l_raw)
  deallocate(narcv)
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Correct flux variation in profile data (ITERDB, etc) 
  !
  if (format_type == 0 .or. format_type == 7) then
     ! Case 1: raw gfile mode ; dpsi is undefined
     dpsi_data = psi(nsurf)
     dpsi_efit = psi(nsurf)
     do i=1,nx
        dpsi(i) = (i-1)*dpsi_efit/(nx-1)
     enddo
  else
     ! Case 2: typical case 
     dpsi_data = dpsi(nx)
     dpsi_efit = psi(nsurf)

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
  call fluxfit_driver(1,1,nsurf,narc,r2,z2,verbose_flag)
  !
  allocate(gvec(6,nsurf))
  !
  open(unit=1,file='fluxfit.profile',status='old')
  read(1,*) a
  read(1,*) gvec(:,:)
  close(1)

  allocate(sqdpsi(nx))
  sqdpsi = sqrt(dpsi)

  ! Explicitly set rmin=0 at origin
  gvec(1,1) = 0.0

  ! Use extrapolation to get values of shape parameters at origin
  do i=2,6
     call bound_extrap(fa,fb,gvec(i,:),psi,nsurf)
     gvec(i,1) = fa
  enddo

  ! Map shape coefficients onto poloidal flux (dpsi) grid:
  call cub_spline(sqrt(psi),gvec(1,:),nsurf,sqdpsi,rmin,nx)
  call cub_spline(sqrt(psi),gvec(2,:),nsurf,sqdpsi,zmag,nx)
  call cub_spline(sqrt(psi),gvec(3,:),nsurf,sqdpsi,rmaj,nx)
  call cub_spline(sqrt(psi),gvec(4,:),nsurf,sqdpsi,kappa,nx)
  call cub_spline(sqrt(psi),gvec(5,:),nsurf,sqdpsi,delta,nx)
  call cub_spline(sqrt(psi),gvec(6,:),nsurf,sqdpsi,zeta,nx)

  ! Map q and p into poloidal flux grid
  if (nop_flag == 0) then
     print '(a)','INFO: (prgen) Using total pressure from gfile.'
     call cub_spline(sqrt(psi),p_omfit,nsurf,sqdpsi,p_tot,nx)
  endif
  if (noq_flag == 0) then 
     print '(a)','INFO: (prgen) Using safety factor (q) from gfile.'
     call cub_spline(sqrt(psi),q_omfit,nsurf,sqdpsi,q,nx)
  endif

  ! Explicitly set rmin=0 at origin
  rmin(1) = 0.0

  deallocate(gvec)
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Get general geometry coefficients
  !
  call fluxfit_driver(2,nfourier,nsurf,narc,r2,z2,verbose_flag)

  allocate(g3vec(4,0:nfourier,nsurf))
  allocate(g3rho(4,0:nfourier,nx))

  open(unit=1,file='fluxfit.geo',status='old')
  read(1,*) ip
  read(1,*) g3vec(:,:,:)
  close(1)

  ! Explicitly set rmin=0 at origin
  g3vec(:,1:nfourier,1) = 0.0

  ! Extrapolate centers (R0,Z0) to origin
  do i=1,4
     call bound_extrap(fa,fb,g3vec(i,0,:),psi,nsurf)
     g3vec(i,0,1) = fa
  enddo

  ! Map Fourier coefficients onto poloidal flux (dpsi) grid 
  do i=1,4
     do ip=0,nfourier
        call cub_spline(sqrt(psi),g3vec(i,ip,:),nsurf,sqdpsi,g3rho(i,ip,:),nx)
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
  !----------------------------------------------------------------

  ! Clean up
  deallocate(g3vec)
  deallocate(g3rho)
  deallocate(r2)
  deallocate(z2)
  deallocate(psi)
  deallocate(q_omfit)
  deallocate(p_omfit)

end subroutine prgen_read_omfit
