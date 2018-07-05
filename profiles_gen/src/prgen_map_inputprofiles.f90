!------------------------------------------------------------
! prgen_map_inputprofiles.f90
!
! PURPOSE:
!  Map input profiles plus gfile data.
!------------------------------------------------------------

subroutine prgen_map_inputprofiles

  use prgen_globals
  use EXPRO_interface

  implicit none

  integer :: i
  real :: xnew(nx)

  EXPRO_n_exp = nx
  call EXPRO_alloc('./',1) 
  call EXPRO_read

  if (efit_method > 1) then

     ! Compute rho, bref and arho based on GATO dpsi and q.
     call prgen_get_chi(nx,q,kappa,rmin,dpsi,rho,EXPRO_b_ref,EXPRO_arho)

     ! Align with new (GATO) rho grid, not old rho grid, EXPRO_rho:
     call cub_spline(EXPRO_rho,EXPRO_te,nx,rho,xnew,nx)
     EXPRO_te = xnew
     call cub_spline(EXPRO_rho,EXPRO_ne,nx,rho,xnew,nx)
     EXPRO_ne = xnew
     call cub_spline(EXPRO_rho,EXPRO_z_eff,nx,rho,xnew,nx)
     EXPRO_z_eff = xnew
     call cub_spline(EXPRO_rho,EXPRO_w0,nx,rho,xnew,nx)
     EXPRO_w0 = xnew
     call cub_spline(EXPRO_rho,EXPRO_flow_mom,nx,rho,xnew,nx)
     EXPRO_flow_mom = xnew
     call cub_spline(EXPRO_rho,EXPRO_pow_e,nx,rho,xnew,nx)
     EXPRO_pow_e = xnew
     call cub_spline(EXPRO_rho,EXPRO_pow_i,nx,rho,xnew,nx)
     EXPRO_pow_i = xnew
     call cub_spline(EXPRO_rho,EXPRO_pow_ei,nx,rho,xnew,nx)
     EXPRO_pow_ei = xnew
     call cub_spline(EXPRO_rho,EXPRO_flow_beam,nx,rho,xnew,nx)
     EXPRO_flow_beam = xnew
     call cub_spline(EXPRO_rho,EXPRO_flow_wall,nx,rho,xnew,nx)
     EXPRO_flow_wall = xnew
     call cub_spline(EXPRO_rho,EXPRO_ptot,nx,rho,xnew,nx)
     EXPRO_ptot = xnew
     do i=1,10
        call cub_spline(EXPRO_rho,EXPRO_ni(i,:),nx,rho,xnew,nx)
        EXPRO_ni(i,:) = xnew(:)
        call cub_spline(EXPRO_rho,EXPRO_ti(i,:),nx,rho,xnew,nx)
        EXPRO_ti(i,:) = xnew(:)
        call cub_spline(EXPRO_rho,EXPRO_vtor(i,:),nx,rho,xnew,nx)
        EXPRO_vtor(i,:) = xnew(:)
        call cub_spline(EXPRO_rho,EXPRO_vpol(i,:),nx,rho,xnew,nx)
        EXPRO_vpol(i,:) = xnew(:)
     enddo

     !---------------------------------------------------------
     ! Map profile data into EXPRO interface.
     ! NOTE: EXPRO_alloc already called in prgen_read_inputprofiles
     !
     EXPRO_rho(:)  = rho(:)
     EXPRO_rmin(:) = rmin
     EXPRO_rmaj(:) = rmaj(:)
     ! COORDINATES: set sign of q
     EXPRO_q(:)     = abs(q(:))*ipccw*btccw
     EXPRO_kappa(:) = kappa(:)
     EXPRO_delta(:) = delta(:)
     EXPRO_zeta(:)  = zeta(:)
     EXPRO_zmag(:)  = zmag(:)
     ! COORDINATES: set sign of poloidal flux
     EXPRO_polflux = abs(dpsi(:))*(-ipccw)

  endif

  call swap(EXPRO_ni(:,1:10))  
  call swap(EXPRO_ti(:,1:10))  
  call swap(EXPRO_vpol(:,1:10))  
  call swap(EXPRO_vtor(:,1:10))  

  !---------------------------------------------------------
  ! Read the cer file and overlay
  !
  if (cer_file /= "null") then
    print '(a)','INFO: (prgen_map_inputprofiles) IGNORING cer file'
  endif
  !---------------------------------------------------------

end subroutine prgen_map_inputprofiles

subroutine swap(f)

  use prgen_globals

  implicit none

  real, intent(inout), dimension(10,nx) :: f

  integer :: i,ip
  real, dimension(10,nx) :: ftmp

  ftmp = f
  do i=1,10
     ip     = reorder_vec(i)
     f(i,:) = ftmp(ip,:)
  enddo
  
end subroutine swap
  
