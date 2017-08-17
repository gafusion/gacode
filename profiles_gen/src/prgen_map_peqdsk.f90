!------------------------------------------------------------
! prgen_map_peqdsk.f90
!
! PURPOSE:
!  Map native peqdsk data onto input.profiles standard.  
!------------------------------------------------------------

subroutine prgen_map_peqdsk

  use prgen_globals
  use EXPRO_interface
  
  implicit none
  
  integer :: i
  real, dimension(nx) :: ni_d
  real, dimension(3,nx) :: ni_imp
  real, dimension(nx) :: ni_b
  real, dimension(nx) :: z_eff

  ! Compute rho, bref and arho:
  call prgen_get_chi(nx,q,kappa,rmin,dpsi,rho,peqdsk_bref,peqdsk_arho)

  ni_d(:) = 10*peqdsk_ni(:)
  if(peqdsk_fmt == 0) then
     ! old p-file, assume missing density is carbon
     ni_imp(1,:) = 10*(peqdsk_ne(:)-peqdsk_ni(:)-peqdsk_nb(:))/6.0
     peqdsk_nimp  = 1
     peqdsk_z(1)  = 1.0
     peqdsk_m(1)  = 2.0
     peqdsk_z(2)  = 6.0
     peqdsk_m(2)  = 12.0
     if(peqdsk_nbeams == 1) then
        peqdsk_z(3)  = 1.0
        peqdsk_m(3)  = 2.0
     endif
  else
     ! new p-file
     do i=1,peqdsk_nimp
        ni_imp(i,:) = 10*peqdsk_nz(i,:)
     enddo
  endif
  ni_b(:) = 10*peqdsk_nb(:)

  ! Compute Z_eff
  z_eff(:) = peqdsk_z(1)**2 * ni_d(:)
  do i=1,peqdsk_nimp
     z_eff(:) = z_eff(:)+peqdsk_z(1+i)**2 * ni_imp(i,:)
  enddo
  z_eff(:) = z_eff(:)+peqdsk_z(1+peqdsk_nimp+1)**2 * ni_b(:)
  z_eff(:) = z_eff(:)/(10*peqdsk_ne(:))

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  EXPRO_n_exp = nx
  call EXPRO_alloc('./',1)
  !
  EXPRO_rho(:)  = rho(:)
  EXPRO_rmin(:) = rmin(:)
  EXPRO_rmaj(:)  = rmaj(:)
  ! COORDINATES: set sign of q
  EXPRO_q(:)      = abs(q(:))*ipccw*btccw
  EXPRO_kappa(:)  = kappa(:)
  EXPRO_delta(:)  = delta(:)
  EXPRO_te(:)     = peqdsk_te(:)
  EXPRO_ne(:)     = peqdsk_ne(:)*10
  EXPRO_z_eff(:)  = z_eff(:)
  ! COORDINATES: -ipccw accounts for DIII-D toroidal angle convention
  EXPRO_w0(:)        = -ipccw*1e3*peqdsk_omgeb(:) 
  EXPRO_flow_mom(:)  = 0.0      ! flow_mom
  EXPRO_pow_e(:)     = 0.0      ! pow_e
  EXPRO_pow_i(:)     = 0.0      ! pow_i 
  EXPRO_pow_ei(:)    = 0.0      ! pow_ei_exp
  EXPRO_zeta(:)      = zeta(:)
  EXPRO_flow_beam(:) = 0.0      ! flow_beam
  EXPRO_flow_wall(:) = 0.0      ! flow_wall_exp
  EXPRO_zmag(:)      = zmag(:)  
  EXPRO_ptot(:)      = p_tot(:)      
  ! COORDINATES: set sign of poloidal flux
  EXPRO_polflux = abs(dpsi(:))*(-ipccw)

  ! ni, nc, nb
  EXPRO_ni(1,:) = ni_d(:)
  do i=1,peqdsk_nimp
     EXPRO_ni(1+i,:) = ni_imp(i,:)
  enddo
  EXPRO_ni(1+peqdsk_nimp+1,:) = ni_b(:)

  ! ti, tc, tb
  EXPRO_ti(1,:) = peqdsk_ti(:)
  do i=1,peqdsk_nimp
     EXPRO_ti(1+i,:) = peqdsk_ti(:)
  enddo
  do i=1,nx
     if (peqdsk_nb(i) > epsilon(0.)) then
        EXPRO_ti(1+peqdsk_nimp+1,i) = peqdsk_pb(i)/(peqdsk_nb(i)*10)/1.602
     else
        EXPRO_ti(1+peqdsk_nimp+1,i) = 0.0
     endif
  enddo

  ! vphi
  ! COORDINATES: -ipccw accounts for DIII-D toroidal angle convention
  EXPRO_vtor(:,:) = 0.0
  EXPRO_vtor(2,:) = -ipccw*1e3*peqdsk_omegat(:)*(rmaj(:)+rmin(:))

  ! vpol
  EXPRO_vpol(:,:) = 0.0

end subroutine prgen_map_peqdsk
