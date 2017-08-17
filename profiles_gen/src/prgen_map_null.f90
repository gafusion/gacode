!------------------------------------------------------------
! prgen_map_null.f90
!
! PURPOSE:
!  Map shape data only.  
!------------------------------------------------------------

subroutine prgen_map_null

  use prgen_globals
  use EXPRO_interface

  implicit none

  ! Compute rho, bref and arho:
  call prgen_get_chi(nx,q,kappa,rmin,dpsi,rho,null_bref,null_arho)

  !---------------------------------------------------------
  ! Map profile data into EXPRO interface names:
  !
  EXPRO_n_exp = nx
  call EXPRO_alloc('./',1)
  
  EXPRO_rho  = rho
  EXPRO_rmin = rmin(:)
  EXPRO_rmaj = rmaj(:)
  ! COORDINATES: set sign of q
  EXPRO_q = abs(q(:))*ipccw*btccw
  EXPRO_ptot = p_tot(:)
  EXPRO_kappa = kappa(:)
  EXPRO_delta = delta(:)
  EXPRO_zeta = zeta(:)
  EXPRO_zmag = zmag(:)
  EXPRO_te = 1.0
  EXPRO_ne = 1.0

  ! COORDINATES: set sign of poloidal flux
  EXPRO_polflux = abs(dpsi(:))*(-ipccw)

  ! Ion temperatures and densities
  EXPRO_ni(1,:) = 1.0
  EXPRO_ti(1,:) = 1.0

 end subroutine prgen_map_null

