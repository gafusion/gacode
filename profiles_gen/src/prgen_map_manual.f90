!------------------------------------------------------------
! prgen_map_manual.f90
!
! PURPOSE:
!  Map MANUAL data
!------------------------------------------------------------

subroutine prgen_map_manual

  use prgen_globals
  use expro
  
  implicit none

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  expro_n_exp = nx
  expro_n_ion = 2
  call expro_init(1)
  !
  expro_rho(:)   = rho(:)
  expro_rmin(:)  = rmin(:)
  expro_rmaj(:)  = rmaj(:)
  expro_te(:)    = te_kev(:)
  expro_ne(:)    = ne_e19m3(:)
  expro_w0(:)    = omega0(:)

  expro_mass(1) = 2.0
  expro_z(1) = 1.0
  expro_name(1) = 'D'
  expro_type(1) = type_therm
  
  expro_ni(1,:) = ni_e19m3(:)
  expro_ti(1,:) = ti_kev(:)

  expro_mass(2) = 12.0
  expro_z(2) = 6.0
  expro_name(2) = 'C'
  expro_type(2) = type_therm
  
  expro_ni(2,:) = nz_e19m3(:)
  expro_ti(2,:) = ti_kev(:)
  
end subroutine prgen_map_manual
