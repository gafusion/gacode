!------------------------------------------------------------
! prgen_map_null.f90
!
! PURPOSE:
!  Map shape data only.  
!------------------------------------------------------------

subroutine prgen_map_null

  use prgen_globals
  use expro

  implicit none

  call prgen_get_chi(nx,q,dpsi,rho,torfluxa)

  !---------------------------------------------------------
  ! Map profile data into expro interface names:
  !
  expro_n_exp = nx
  expro_n_ion = 1
  call expro_init(1)
  
  expro_rho  = rho
  expro_rmin = rmin(:)
  expro_rmaj = rmaj(:)
  expro_ptot = p_tot(:)
  expro_te = 1.0
  expro_ne = 1.0

  ! Ion temperatures and densities
  expro_ni(1,:) = 1.0
  expro_ti(1,:) = 1.0

 end subroutine prgen_map_null

