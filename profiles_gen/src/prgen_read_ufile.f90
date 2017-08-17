!--------------------------------------------------------------
! prgen_read_ufile.f90
!
! PURPOSE:
!  Read profiles from post-processed (by gacode_ufile_tool) 
!  UFILEs. 
!--------------------------------------------------------------

subroutine prgen_read_ufile

  use prgen_globals

  implicit none

  integer :: i,ix,j
  real :: chi_ta
  real, dimension(:), allocatable :: chi_t
  character(len=16) :: a

  open(unit=1,file='out.BT.ave',status='old') 
  read(1,*) ufile_bref
  close(1)

  ! Chi_t(a), where Chi_t is the total toroidal flux
  open(unit=1,file='out.PHIA.ave',status='old') 
  read(1,*) chi_ta
  chi_ta = chi_ta/(2*pi)
  close(1)

  open(unit=1,file='out.com',status='old') 
  read(1,*) raw_data_file
  read(1,*) ufile_tok
  read(1,*) a
  read(1,*) a
  read(1,*) ufile_shot
  read(1,*) ufile_time
  read(1,*) a
  do i=1,8
     read(1,*) ufile_z(i)
     read(1,*) ufile_m(i)
  enddo
  close(1)
 
  ! NOTE: this is a hardwired dimension
  nx = 51

  call allocate_internals
  call allocate_ufile_vars

  allocate(chi_t(nx))

  ! Classic normalized rho
  do i=1,nx
     rho(i) = (i-1)/(nx-1.0)
  enddo

  call ufile_mapper('out.NE.ave',rho,ufile_ne,nx,1)
  call ufile_mapper('out.TE.ave',rho,ufile_te,nx,1)
  call ufile_mapper('out.TI.ave',rho,ufile_ti(:,1),nx,1)
  call ufile_mapper('out.ZEFFR.ave',rho,ufile_zeff,nx,0)
  call ufile_mapper('out.RMAJOR.ave',rho,rmaj,nx,0)
  call ufile_mapper('out.RMINOR.ave',rho,rmin,nx,0)
  ! Correct rmin on magnetic axis
  rmin(1) = 0.0
  call ufile_mapper('out.Q.ave',rho,q,nx,0)
  call ufile_mapper('out.KAPPAR.ave',rho,kappa,nx,0)
  call ufile_mapper('out.DELTAR.ave',rho,delta,nx,0)
  call ufile_mapper('out.PRES.ave',rho,ufile_pres,nx,1)
  call ufile_mapper('out.VROT.ave',rho,ufile_vrot,nx,0)
  call ufile_mapper('out.VROTM.ave',rho,ufile_vrotm,nx,0)
  call ufile_mapper('out.VOLUME.ave',rho,ufile_volume,nx,1)
  call ufile_mapper('out.QFUSI.ave',rho,ufile_qfusi,nx,0)
  call ufile_mapper('out.QFUSE.ave',rho,ufile_qfuse,nx,0)
  call ufile_mapper('out.QNBII.ave',rho,ufile_qnbii,nx,0)
  call ufile_mapper('out.QNBIE.ave',rho,ufile_qnbie,nx,0)
  call ufile_mapper('out.QICRHI.ave',rho,ufile_qicrhi,nx,0)
  call ufile_mapper('out.QICRHE.ave',rho,ufile_qicrhe,nx,0)
  call ufile_mapper('out.QEI.ave',rho,ufile_qei,nx,0)
  call ufile_mapper('out.QRAD.ave',rho,ufile_qrad,nx,0)
  call ufile_mapper('out.QECHE.ave',rho,ufile_qeche,nx,0)
  call ufile_mapper('out.QECHI.ave',rho,ufile_qechi,nx,0)
  call ufile_mapper('out.QLHE.ave',rho,ufile_qlhe,nx,0)
  call ufile_mapper('out.QLHI.ave',rho,ufile_qlhi,nx,0)
  call ufile_mapper('out.QOHM.ave',rho,ufile_qohm,nx,0)
  call ufile_mapper('out.QWALLI.ave',rho,ufile_qwalli,nx,0)
  call ufile_mapper('out.QWALLE.ave',rho,ufile_qwalle,nx,0)
  call ufile_mapper('out.TORQ.ave',rho,ufile_torq,nx,0)
  call ufile_mapper('out.SWALL.ave',rho,ufile_swall,nx,0)
  call ufile_mapper('out.SNBII.ave',rho,ufile_snbii,nx,0)
  call ufile_mapper('out.SNBIE.ave',rho,ufile_snbie,nx,0)

  i=0

  if (ufile_m(1) > 0.0) then
     i = i+1
     call ufile_mapper('out.NM1.ave',rho,ufile_ni(:,i),nx,1)
     ufile_type(i)='therm'
  endif
  if (ufile_m(2) > 0.0) then
     i = i+1
     call ufile_mapper('out.NM2.ave',rho,ufile_ni(:,i),nx,1)
     ufile_type(i)='therm'
  endif
  if (ufile_m(3) > 0.0) then
     i = i+1
     call ufile_mapper('out.NM3.ave',rho,ufile_ni(:,i),nx,1)
     ufile_type(i)='therm'
  endif
  if (ufile_m(4) > 0.0) then
     i = i+1
     call ufile_mapper('out.NM4.ave',rho,ufile_ni(:,i),nx,1)
     ufile_type(i)='therm'
  endif
  if (ufile_m(5) > 0.0) then
     i = i+1
     call ufile_mapper('out.NM5.ave',rho,ufile_ni(:,i),nx,1)
     ufile_type(i)='therm'
  endif
  if (ufile_m(6) > 0.0) then
     i = i+1
     call ufile_mapper('out.NFAST1.ave',rho,ufile_ni(:,i),nx,1)
     ufile_type(i)='fast'
  endif
  if (ufile_m(7) > 0.0) then
     i = i+1
     call ufile_mapper('out.NFAST2.ave',rho,ufile_ni(:,i),nx,1)
     ufile_type(i)='fast'
  endif
  if (ufile_m(8) > 0.0) then
     i = i+1
     call ufile_mapper('out.NFAST3.ave',rho,ufile_ni(:,i),nx,1)
     ufile_type(i)='fast'
  endif

  ufile_nion=i

  ! Now compress m/z
  i=1
  do ix=1,8
     if (int(ufile_m(i)) == 0) then
        do j=i+1,8
           ufile_m(j-1) = ufile_m(j)
           ufile_z(j-1) = ufile_z(j)
        enddo
     else
        i = i+1
     endif
  enddo

  ! Set ion temperatures
  do i=2,ufile_nion
     if (ufile_type(i) == 'therm') then
        ufile_ti(:,i) = ufile_ti(:,1) 
     else
        ufile_ti(:,i) = 5*ufile_ti(:,1)
        print '(a)','WARNING: (prgen_read_ufile) Setting bogus fast-ion temperature.'
     endif
  enddo

  ! Compute the quasineutrality error with all ions:

  quasi_err = 0.0
  ix = min(ufile_nion,n_ion_max)
  do i=1,nx
     quasi_err = quasi_err+sum(ufile_ni(i,1:ix)*ufile_z(1:ix))
  enddo
  quasi_err = abs(quasi_err/sum(ufile_ne(:))-1.0)

  !------------------------------------------------------------------------
  ! Use classic parameterization chi_t = B_ref/2 rho^2
  ! where chi_t is toroidal flux over 2pi.
  !
  ufile_arho = sqrt(2*chi_ta/abs(ufile_bref))
  chi_t(:)   = chi_ta*rho**2
  !
  ! Compute psi_p by integrating d(chi_t)/q = d(psi_p)
  ! using the trapezoidal rule
  !
  dpsi(1) = 0.0
  do i=2,nx
     dpsi(i) = dpsi(i-1) + 0.5*(1.0/q(i-1)+1.0/q(i))*(chi_t(i)-chi_t(i-1))
  enddo
  !------------------------------------------------------------------------

end subroutine prgen_read_ufile


!---------------------------------------------------------------------
! Routine to read split 2D UFILE data.
!---------------------------------------------------------------------

subroutine ufile_mapper(file,x,y,nx,neg_protect)

  character(len=*), intent(in) :: file
  integer, intent(in) :: nx
  integer, intent(in) :: neg_protect
  real, intent(in) :: x(nx)
  real, intent(inout) :: y(nx)

  real, dimension(:), allocatable :: x0
  real, dimension(:), allocatable :: y0
  real :: ya,yb
  integer :: ierr

  open(unit=1,file='out.dim',status='old')
  read(1,*) nx0
  close(1)

  allocate(x0(0:nx0+1))
  allocate(y0(0:nx0+1))

  open(unit=1,file=trim(file),status='old',iostat=ierr)
  if (ierr == 0) then
     ! File exists
     do i=1,nx0
        read(1,*) x0(i),y0(i)
     enddo
  else
     ! File does NOT exist: pad output with zeros
     y(:) = 0.0
     return
  endif

  ! Check for zero profile (1.0 is effectively zero).
  if (sum(y0(1:nx0)) < 1.0) then 
     y(:) = 0.0
     return
  endif
  
  if (x0(1) > 0.0 .and. x0(nx0) < 1.0) then
     ! Extrapolate to 0 and 1
     x0(0)     = 0.0
     x0(nx0+1) = 1.0
     if (neg_protect==0) then
        call bound_extrap(ya,yb,y0,x0,nx0+2)
     else
        call bound_exp(ya,yb,y0,x0,nx0+2)
     endif
     y0(0) = ya
     y0(nx0+1) = yb
     call cub_spline(x0,y0,nx0+2,x,y,nx)
  else if (x0(1) > 0.0 .and. abs(x0(nx0)-1.0) <= epsilon(0.0)) then
     ! Extrapolate to 0 only
     x0(0) = 0.0
     call bound_extrap(ya,yb,y0(0:nx0),x0(0:nx0),nx0+1)
     y0(0) = ya
     call cub_spline(x0(0:nx0),y0(0:nx0),nx0+1,x,y,nx)
  else
     print '(a)','WARNING: (prgen) Bad data in prgen_read_ufile: '//file
     y(:) = 0.0
     return
  endif

  deallocate(x0)
  deallocate(y0)

end subroutine ufile_mapper


!---------------------------------------------------------------------
! Routine, based on math/bound_extrap.f90, to use extrapolation based
! on an inferred scale length at the right end of the data.  This is 
! appropriate for profiles (n,T for example) that must remain positive.
!---------------------------------------------------------------------

subroutine bound_exp(fa,fb,f,r,n)

  implicit none

  integer, intent(in) :: n

  real, intent(inout) :: fa
  real, intent(inout) :: fb
  real, intent(in), dimension(n) :: f
  real, intent(in), dimension(n) :: r

  real :: r1,r2,r3
  real :: ra,rb
  real :: f1,f2,f3
  real :: l,c

  ! Left boundary

  ra = r(1)

  r2 = r(2)
  r3 = r(3)
  
  f2 = f(2)
  f3 = f(3)

  fa = (ra-r2)/(r3-r2)*f3+(r3-ra)/(r3-r2)*f2

  ! Right boundary (fit to exp(-r/l) with scale length l)

  rb = r(n)

  r1 = r(n-2)
  r2 = r(n-1)
  
  f1 = f(n-2)
  f2 = f(n-1)

  l = (r2-r1)/log(f1/f2)
  c = f2*exp(r2/l)

  fb = c*exp(-rb/l)

end subroutine bound_exp
