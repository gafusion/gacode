!----------------------------------------------------------------
! prgen_geometry.f90
!
! PURPOSE:
!  Read data from flux-surface mapper
!----------------------------------------------------------------

subroutine prgen_geometry

  use prgen_globals
  use expro

  implicit none

  integer :: npsi,nf,i,j,ip
  real, dimension(:,:), allocatable :: efit_si,efit_ci
  real, dimension(:), allocatable :: efit_rho,efit_psi,efit_q,efit_p,efit_fpol
  real, dimension(:), allocatable :: efit_rmin,efit_rmaj,efit_kappa,efit_zmaj
  real, dimension(:,:,:), allocatable :: g3vec
  real, dimension(:,:,:), allocatable :: g3rho
  real :: psi_sep,pratio

  !----------------------------------------------------------------
  open(unit=1,file='out.dim',status='old')
  read(1,*) npsi
  read(1,*) nf
  read(1,*) rcentr
  read(1,*) bcentr
  read(1,*) current
  read(1,*) psi0
  read(1,*) psi1
  read(1,*) psi_sep
  close(1)

  allocate(efit_si(npsi,0:nf-1))
  allocate(efit_ci(npsi,0:nf-1))
  allocate(efit_rmin(npsi))
  allocate(efit_rmaj(npsi))
  allocate(efit_kappa(npsi))
  allocate(efit_zmaj(npsi))
  allocate(efit_psi(npsi))
  allocate(efit_q(npsi))
  allocate(efit_p(npsi))
  allocate(efit_fpol(npsi))
  allocate(efit_rho(npsi))

  open(unit=1,file='out.data',status='old',access='stream')
  read(1) efit_psi
  read(1) efit_q
  read(1) efit_p
  read(1) efit_fpol
  read(1) efit_si
  read(1) efit_ci
  read(1) efit_rmin
  read(1) efit_rmaj
  read(1) efit_kappa
  read(1) efit_zmaj
  close(1)

  ! PFILE only contains normalized flux
  if (format_type == 3) then
     dpsi = (psi1-psi0)*dpsi
  endif
  
  efit_psi = efit_psi-efit_psi(1)

  ! Deal with null statefile
  if (format_type == 0) then
     nx = npsi
     call prgen_allocate('')
     dpsi = efit_psi
  endif

  !--------------------------------------------------------------------------------------
  ! Flux diagnostics 
  print 10,'INFO: (prgen_geometry)    GACODE/EFIT dpsi:',(psi_sep-psi0)/(psi1-psi0)
  print 10,'INFO: (prgen_geometry)    MAPPED/EFIT dpsi:',efit_psi(npsi)/(psi1-psi0),' (adjust with -psinorm)'
  print 10,'INFO: (prgen_geometry) STATEFILE/EFIT dpsi:',dpsi(nx)/(psi1-psi0)

  ! Problem condition: dpsi(nx) > efit_psi(npsi) 
  pratio = dpsi(nx)/efit_psi(npsi)
  if (pratio > 1.0) then
     if (pratio > 1.05) then
        print '(a)','WARNING: (prgen_geometry) Detected statefile dpsi(nx) >> mapped dpsi_max [BAD]'
     else
        print '(a)','WARNING: (prgen_geometry) Detected statefile dpsi(nx) > mapped dpsi_max'
     endif
     print '(a,f5.3,a)','WARNING: (prgen_geometry) Compressing statefile flux by ',pratio,' to match mapped flux'
     dpsi(:) = dpsi(:)/(pratio+1e-10)
  endif
  !--------------------------------------------------------------------------------------

  ! Get rho on efit mesh (so have psi,rho,q)
  call prgen_get_chi(npsi,efit_q,efit_psi,efit_rho,torfluxa)
  ! Get q
  call cub_spline(efit_psi,efit_q,npsi,dpsi,q,nx)
  ! Get rho
  call prgen_get_chi(nx,q,dpsi,rho,torfluxa)

  call bound_interp(efit_rho,efit_rmin,npsi,rho,rmin,nx)
  call bound_interp(efit_psi,efit_q,npsi,dpsi,q,nx)
  call bound_interp(efit_psi,efit_p,npsi,dpsi,p_tot,nx)
  call bound_interp(efit_psi,efit_fpol,npsi,dpsi,fpol,nx)
  call bound_interp(efit_rho,efit_zmaj,npsi,rho,zmag,nx)
  call bound_interp(efit_rho,efit_rmaj,npsi,rho,rmaj,nx)
  call bound_interp(efit_rho,efit_kappa,npsi,rho,kappa,nx)

  ! Miller extended harmonic shape coefficients
  call bound_interp(efit_rho,efit_ci(:,0),npsi,rho,shape_cos(0,:),nx)

  if (nf-1 > 6) then
     print '(a)','ERROR: (prgen_geometry) max nharm is 6'
     stop
  endif

  do i=1,nf-1
     call bound_interp(efit_rho,efit_si(:,i),npsi,rho,shape_sin(i,:),nx)
     call bound_interp(efit_rho,efit_ci(:,i),npsi,rho,shape_cos(i,:),nx)
  enddo

  delta = sin(shape_sin(1,:))
  zeta  = -shape_sin(2,:)
  !==============================================================================

  if (nfourier > 0) then

     ! Legacy direct Fourier representation
     allocate(g3vec(npsi,0:nfourier,4))
     allocate(g3rho(nx,0:nfourier,4))

     open(unit=1,file='fluxfit.geo',status='old',access='stream')
     do i=1,4
        read(1) g3vec(:,:,i)
     enddo
     close(1)

     ! Map results onto poloidal flux (dpsi) grid:
     ! NOTE: need sqrt here to get sensible behaviour as r -> 0.
     do i=1,4
        do ip=0,nfourier
           call bound_interp(efit_rho,g3vec(:,ip,i),npsi,rho,g3rho(:,ip,i),nx)
        enddo
     enddo

     open(unit=1,file='input.gacode.geo',status='replace')
     write(1,'(a)') '# input.gacode.geo'
     write(1,'(a)') '#'
     write(1,'(a)') '# File format:'
     write(1,'(a)') '#-------------------'
     write(1,'(a)') '# nfourier'
     write(1,'(a)') '# a[4,0:nfourier,nx]'
     write(1,'(a)') '#-------------------'
     write(1,'(a)') '#'
     write(1,'(a)') '# NOTE: nx=EXPRO_n_exp is defined in input.gacode'
     write(1,*) nfourier
     do j=1,nx
        do ip=0,nfourier
           do i=1,4
              write(1,'(1pe20.13)') g3rho(j,ip,i)
           enddo
        enddo
     enddo
     close(1)
     !----------------------------------------------------

     ! Clean up
     deallocate(g3vec)
     deallocate(g3rho)

  endif

  deallocate(efit_si)
  deallocate(efit_ci)
  deallocate(efit_rmin)
  deallocate(efit_rmaj)
  deallocate(efit_kappa)
  deallocate(efit_zmaj)
  deallocate(efit_psi)
  deallocate(efit_q)
  deallocate(efit_p)
  deallocate(efit_fpol)
  deallocate(efit_rho)

10 format(t1,a,f9.6,a)

end subroutine prgen_geometry
