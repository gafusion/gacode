!------------------------------------------------------------
! prgen_map_corsica.f90
!
! PURPOSE:
!  Map native corsica data onto input.profiles standard.  
!------------------------------------------------------------

subroutine prgen_map_corsica

  use prgen_globals
  use expro
  
  implicit none

  ! Ignore corsica input "rho", "rmin", and "q"; use gato and corsica poloidal flux
  call prgen_get_chi(nx,q,dpsi,rho,torfluxa)

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  expro_n_exp = corsica_nvals
  expro_n_ion = 2
  call expro_init(1)
  !
  expro_rho(:)       = rho(:)
  expro_rmin(:)      = rmin(:)
  expro_rmaj(:)      = rmaj(:)
  expro_te(:)        = corsica_te(:)
  expro_ne(:)        = corsica_ne(:)*10.0
  expro_z_eff(:)     = corsica_zeff(:)
  expro_w0(:)        = 0.0      ! omega
  expro_ptot(:)      = p_tot(:)

  ! Construct ion densities and temperatures assuming corsica ion species
  ! (n D+T) is two species, each with 1/2 n_corsica and same temperature

  expro_mass(:) = 2.0
  expro_z(:) = 1.0
  expro_name(:) = 'D'
  expro_type(:) = type_therm
  
  ! ni
  expro_ni(1,:) = 0.5 * corsica_ndt(:)*10.0
  expro_ni(2,:) = 0.5 * corsica_ndt(:)*10.0
  
  ! ti
  expro_ti(1,:) = corsica_ti(:)
  expro_ti(2,:) = corsica_ti(:)

  expro_vtor(:,:) = 0.0
  expro_vpol(:,:) = 0.0
  
end subroutine prgen_map_corsica
