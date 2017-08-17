!------------------------------------------------------------
! prgen_map_corsica.f90
!
! PURPOSE:
!  Map native corsica data onto input.profiles standard.  
!------------------------------------------------------------

subroutine prgen_map_corsica

  use prgen_globals
  use EXPRO_interface
  
  implicit none

  ! Compute rho, bref and arho:
  ! Ignore corsica input "rho", "rmin", and "q"; use gato and corsica poloidal
  ! flux
  call prgen_get_chi(nx,q,kappa,rmin,dpsi,rho,corsica_bref,corsica_arho)

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  EXPRO_n_exp = corsica_nvals
  call EXPRO_alloc('./',1)
  !
  EXPRO_rho(:)       = rho(:)
  EXPRO_rmin(:)      = rmin(:)
  EXPRO_rmaj(:)      = rmaj(:)
  EXPRO_q(:)         = q(:)
  EXPRO_kappa(:)     = kappa(:)
  EXPRO_delta(:)     = delta(:)
  EXPRO_te(:)        = corsica_te(:)
  EXPRO_ne(:)        = corsica_ne(:)*10.
  EXPRO_z_eff(:)     = corsica_zeff(:)
  EXPRO_w0(:)        = 0.0      ! omega
  EXPRO_flow_mom(:)  = 0.0      ! flow_mom
  EXPRO_pow_e(:)     = 0.0      ! pow_e
  EXPRO_pow_i(:)     = 0.0      ! pow_i 
  EXPRO_pow_ei(:)    = 0.0      ! pow_ei_exp
  EXPRO_zeta(:)      = zeta(:)
  EXPRO_flow_beam(:) = 0.0      ! flow_beam
  EXPRO_flow_wall(:) = 0.0      ! flow_wall_exp
  EXPRO_zmag(:)      = zmag(:)  
  EXPRO_ptot(:)      = p_tot(:)
  EXPRO_polflux(:)   = dpsi(:)

  ! Construct ion densities and temperatures assuming corsica ion species
  ! (n D+T) is two species, each with 1/2 n_corsica and same temperature

  ! ni
  EXPRO_ni(1,:) = 0.5 * corsica_ndt(:)*10.0
  EXPRO_ni(2,:) = 0.5 * corsica_ndt(:)*10.0
  
  ! ti
  EXPRO_ti(1,:) = corsica_ti(:)
  EXPRO_ti(2,:) = corsica_ti(:)

  ! vphi
  EXPRO_vtor(:,:) = 0.0

  ! vpol
  EXPRO_vpol(:,:) = 0.0
  
end subroutine prgen_map_corsica
