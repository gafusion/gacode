!---------------------------------------------------------
! prgen_read_dskgato.f90
!
! PURPOSE:
!  Read DSKGATO format.
!  efit_method=4: old dskgato
!  efit_method=5: new dskgato
!----------------------------------------------------------

subroutine prgen_read_dskgato

  use prgen_globals

  implicit none

  integer :: i,ip
  integer :: neqtyp
  integer :: ntht,neqsym
  real :: fa,fb
  real :: dummy(4)
  real, dimension(:), allocatable :: psi
  real, dimension(:), allocatable :: q_dsk
  real, dimension(:), allocatable :: pdum
  real, dimension(:), allocatable :: tdum
  real, dimension(:,:), allocatable :: xs
  real, dimension(:,:), allocatable :: zs

  real, dimension(:,:), allocatable :: gvec
  real, dimension(:,:,:), allocatable :: g3vec
  real, dimension(:,:,:), allocatable :: g3rho

  character (len=70) :: cdum

  !---------------------------------------------------
  ! 
  if (efit_method == 4) then
     neqtyp=1
  else
     neqtyp=0
  endif
  open(unit=1,file='grid.dat',status='old')
  if (neqtyp == 0) then
     read(1,'(a)') cdum
     read(1,'(a)') cdum
     read(1,*) nsurf,ntht,neqsym
     print '(a)','INFO: (prgen) Assuming new-type dskgato flux-surface format.'
  else
     read(1,*) nsurf,ntht
     neqsym = 1
     print '(a)','INFO: (prgen) Assuming old-type dskgato flux-surface format.' 
  endif
  ! Accounting for magnetic axis
  nsurf = nsurf-1

  if (neqsym == 0) then
     narc = ntht
  else
     ! Repeat point
     narc = 2*(ntht-1)+1
  endif

  read (1,10) dummy(1:4)

  if (neqtyp == 0) read(1,10) dummy(1:2) ! totcur,axddxz
  if (neqtyp == 1) read(1,10) dummy(1:3) ! totcur,axddxz,dnnorm

  allocate(tdum(ntht))
  allocate(psi(0:nsurf))
  allocate(q_dsk(0:nsurf))
  allocate(pdum(0:nsurf))
  allocate(xs(narc,0:nsurf))
  allocate(zs(narc,0:nsurf))

  read(1,10) psi(:)
  read(1,10) pdum(:) ! fval
  read(1,10) pdum(:) ! ffprime
  read(1,10) pdum(:) ! sp
  read(1,10) pdum(:) ! pprime
  read(1,10) q_dsk(:) ! qsfin
  read(1,10) pdum(:) ! nel

  read(1,10) tdum(:) ! seqdpdr
  read(1,10) tdum(:) ! seqdpdz

  ! These arrays should repeat endpoint
  read(1,10) ((xs(i,ip),ip=0,nsurf),i=1,ntht)
  read(1,10) ((zs(i,ip),ip=0,nsurf),i=1,ntht)
  close(1)

  if (neqsym == 1) then
     zs(1,:)    = 0.0
     zs(ntht,:) = 0.0
     do i=ntht+1,narc
        ip = narc+1-i
        xs(i,:) = +xs(ip,:)
        zs(i,:) = -zs(ip,:)
     enddo
  endif

  !do i=1,narc
  !   print *,i,xs(i,nsurf/2),zs(i,nsurf/2)
  !enddo

  ! Set psi(1) = 0
  psi(:) = psi(:)-psi(0)

  ! COORDINATES: set sign of poloidal flux
  if (psi(nsurf) < 0.0) then
     ! Flux negative
     psi(:) = -psi(:)
     ipccw  = 1
  else
     ipccw = -1
  endif
  !----------------------------------------------------

  !----------------------------------------------------------------
  ! Correct flux variation in profile data (ITERDB, etc) 
  !
  if (format_type == 0 .or. format_type == 7) then
     ! Case 1: raw gfile mode ; dpsi is undefined
     dpsi_data = psi(nsurf)
     dpsi_efit = psi(nsurf)
     do i=1,nx
        ! Shrink slightly for safety
        dpsi(i) = (i-1)*dpsi_efit/(nx-1+1e-12)
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
  call fluxfit_driver(1,1,nsurf,narc,xs(:,1:nsurf),zs(:,1:nsurf),verbose_flag)
  !
  allocate(gvec(6,0:nsurf))
  !
  open(unit=1,file='fluxfit.profile',status='old')
  read(1,*) cdum
  read(1,*) gvec(:,1:nsurf)
  close(1)

  ! Explicitly set rmin=0 at origin
  gvec(1,0) = 0.0

  ! Use extrapolation to get values of shape parameters at origin
  do i=2,6
     call bound_extrap(fa,fb,gvec(i,:),psi,nsurf+1)
     gvec(i,0) = fa
  enddo

  ! Map shape coefficients onto poloidal flux (dpsi) grid:
  call cub_spline(sqrt(psi),gvec(1,:),nsurf+1,sqrt(dpsi),rmin,nx)
  call cub_spline(psi,gvec(2,:),nsurf+1,dpsi,zmag,nx)
  call cub_spline(psi,gvec(3,:),nsurf+1,dpsi,rmaj,nx)
  call cub_spline(psi,gvec(4,:),nsurf+1,dpsi,kappa,nx)
  call cub_spline(psi,gvec(5,:),nsurf+1,dpsi,delta,nx)
  call cub_spline(psi,gvec(6,:),nsurf+1,dpsi,zeta,nx)

  ! Explicitly set rmin=0 at origin
  rmin(1) = 0.0

  deallocate(gvec)
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Get general geometry coefficients
  !
  call fluxfit_driver(2,nfourier,nsurf,narc,xs(:,1:nsurf),zs(:,1:nsurf),verbose_flag)

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
     call bound_extrap(fa,fb,g3vec(i,0,:),psi,nsurf+1)
     g3vec(i,0,0) = fa
  enddo

  ! Map Fourier coefficients onto poloidal flux (dpsi) grid 
  do i=1,4
     do ip=0,nfourier
        call cub_spline(psi,g3vec(i,ip,:),nsurf+1,dpsi,g3rho(i,ip,:),nx)
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

  if (noq_flag == 0) then
     call cub_spline(psi,q_dsk,nsurf+1,dpsi,q,nx)
  endif

  ! Cleanup
  deallocate(g3vec)
  deallocate(g3rho)
  deallocate(xs)
  deallocate(zs)
  deallocate(psi)
  deallocate(tdum)
  deallocate(pdum)
  deallocate(q_dsk)
 
10 format(1p4e19.12)

end subroutine prgen_read_dskgato
