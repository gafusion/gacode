!--------------------------------------------------------------
! prgen_read_peqdsk.f90
!
! PURPOSE:
!  Extract peqdsk data from native TEXT file.
!--------------------------------------------------------------

subroutine prgen_read_peqdsk

  use prgen_globals

  implicit none

  integer, parameter :: ncol=3
  integer :: i, num
  real    :: y1, y2, y3
  logical :: ierr
  real, dimension(:,:), allocatable :: xv

  !----------------------------------------------------
  ! Read the peqdsk file
  !

  open(unit=1,file='pfile.ne',status='old')
  read(1,*) i
  nx = i
  call allocate_internals
  call allocate_peqdsk_vars
  allocate(xv(ncol,i))
  read(1,*) xv
  peqdsk_psi(:) = xv(1,:)
  deallocate(xv)
  close(1)
  ! Use ne grid as peqdsk_psi grid for all peqdsk quantities
  ! Except if ne grid psi > 1.0, then use fixed, uniform mesh
  if(peqdsk_psi(nx) > 1.0) then
     print '(a)','INFO: (prgen) Using fixed, uniform psi grid.'
     do i=1,nx
        peqdsk_psi(i) = (i-1.0)/(nx-1)
     enddo
  else
     print '(a)','INFO: (prgen) Using pfile.ne psi grid.'
  endif

  ! psi_norm and ne(10^20/m^3)
  open(unit=1,file='pfile.ne',status='old')
  read(1,*) i
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_ne,nx)
  deallocate(xv)
  close(1)

  ! te(KeV)
  inquire(file='pfile.te',exist=ierr)
  if (ierr) then
     open(unit=1,file='pfile.te',status='old')
  else
     open(unit=1,file='pfile.Te',status='old')
  endif
  read(1,*) i
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_te,nx)
  deallocate(xv)
  close(1)

  ! ni(10^20/m^3)
  open(unit=1,file='pfile.ni',status='old')
  read(1,*) i
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_ni,nx)
  deallocate(xv)
  close(1)

  ! ti(KeV)
  inquire(file='pfile.ti',exist=ierr)
  if (ierr) then
     open(unit=1,file='pfile.ti',status='old')
  else
     open(unit=1,file='pfile.Ti',status='old')
  endif
  read(1,*) i
  allocate(xv(ncol,i))
  read(1,*) xv
  call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_ti,nx)
  deallocate(xv)
  close(1)

  ! ptot (KPa)
  inquire(file='pfile.ptot',exist=ierr)
  if (ierr) then
     open(unit=1,file='pfile.ptot',status='old')
     read(1,*) i
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,p_tot,nx)
     deallocate(xv)
     close(1)
  else
     p_tot(:) = 0.0
  endif

  ! nb(10^20/m^3)
  inquire(file='pfile.nb',exist=ierr)
  if (ierr) then
     open(unit=1,file='pfile.nb',status='old')
     read(1,*) i
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_nb,nx)
     deallocate(xv)
     close(1)
     peqdsk_nbeams = 1
  else
     peqdsk_nb(:) = 0.0
     peqdsk_nbeams = 0
  endif

  ! pb(KPa)
  inquire(file='pfile.pb',exist=ierr)
  if (ierr) then
     open(unit=1,file='pfile.pb',status='old')
     read(1,*) i
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_pb,nx)
     deallocate(xv)
     close(1)
  else
     peqdsk_pb(:) = 0.0
  endif

  peqdsk_nimp = 0
  peqdsk_nz(:,:) = 0.0
  inquire(file='pfile.nz1',exist=ierr)
  if (ierr) then
     open(unit=1,file='pfile.nz1',status='old')
     read(1,*) i
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_nz(1,:),nx)
     deallocate(xv)
     close(1)
     peqdsk_nimp = peqdsk_nimp + 1
     inquire(file='pfile.nz2',exist=ierr)
     if (ierr) then
        open(unit=1,file='pfile.nz2',status='old')
        read(1,*) i
        allocate(xv(ncol,i))
        read(1,*) xv
        call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_nz(2,:),nx)
        deallocate(xv)
        close(1)
        peqdsk_nimp = peqdsk_nimp + 1
        inquire(file='pfile.nz3',exist=ierr)
        if (ierr) then
           open(unit=1,file='pfile.nz3',status='old')
           read(1,*) i
           allocate(xv(ncol,i))
           read(1,*) xv
           call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_nz(3,:),nx)
           deallocate(xv)
           close(1)
           peqdsk_nimp = peqdsk_nimp + 1 
           inquire(file='pfile.nz4',exist=ierr)
           if (ierr) then
              print '(a)','ERROR: (prgen) Too many impurity species in pfile'
              stop
           endif
        endif
     endif
  endif

  peqdsk_type(:) = type_therm
  peqdsk_z(:)    = 0.0
  peqdsk_m(:)    = 0.0

  inquire(file='pfile.species',exist=ierr)
  if (ierr) then
     peqdsk_fmt = 1
     open(unit=1,file='pfile.species',status='old')
     read(1,*) num
     if (num > (1 + peqdsk_nimp + peqdsk_nbeams)) then
        print '(a)','ERROR: (prgen) Species number in pfile incorrect'
        stop
     endif
     do i=1,peqdsk_nimp
        read(1,*) y1, y2, y3
        peqdsk_z(1+i) = y2
        peqdsk_m(1+i) = y3
     enddo
     read(1,*) y1, y2, y3
     peqdsk_z(1) = y2
     peqdsk_m(1) = y3
     if (peqdsk_nbeams == 1) then
        read(1,*) y1, y2, y3
        peqdsk_z(1+peqdsk_nimp+1) = y2
        peqdsk_m(1+peqdsk_nimp+1) = y3
        peqdsk_type(1+peqdsk_nimp+1) = type_fast
     endif
  else
     peqdsk_fmt = 0
  endif

  ! omeg(kRad/s)
  inquire(file='pfile.omeg',exist=ierr)
  if (ierr) then
     open(unit=1,file='pfile.omeg',status='old')
     read(1,*) i
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_omegat,nx)
     deallocate(xv)
     close(1)
  else
     peqdsk_omegat(:) = 0.0
  endif

  ! omgeb(kRad/s) [not a typo]
  inquire(file='pfile.omgeb',exist=ierr)
  if (ierr) then
     open(unit=1,file='pfile.omgeb',status='old')
     read(1,*) i
     allocate(xv(ncol,i))
     read(1,*) xv
     call cub_spline(xv(1,:),xv(2,:),i,peqdsk_psi,peqdsk_omgeb,nx)
     deallocate(xv)
     close(1)
  else
     peqdsk_omgeb(:) = 0.0
  endif

  dpsi(:)   = peqdsk_psi(:)-peqdsk_psi(1)
  rmin(:)   = 0.0
  rmaj(:)   = 0.0
  rho(:)    = 0.0
  kappa(:)  = 0.0
  delta(:)  = 0.0
  zmag(:)   = 0.0
  zeta(:)   = 0.0
  omega0(:) = 0.0
  
end subroutine prgen_read_peqdsk
