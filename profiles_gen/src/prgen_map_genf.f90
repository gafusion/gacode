!------------------------------------------------------------
! prgen_map_genf.f90
!
! PURPOSE:
!  Map native GENF data onto input.gacode standard
!------------------------------------------------------------

subroutine prgen_map_genf

  use prgen_globals
  use expro
  
  implicit none

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  expro_n_exp = nx
  expro_n_ion = 1
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

  expro_qohme(:) = qohm(:)
  
end subroutine prgen_map_genf
