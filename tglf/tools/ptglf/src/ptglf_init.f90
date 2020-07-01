subroutine ptglf_init

  use ptglf_globals
  use tglf_interface

  implicit none

  ! Want fluxes from TGLF
  tglf_use_transport_model_in = .true.

  tglf_ns_in = 2

  tglf_zs_in(1) = -1.0
  tglf_zs_in(2) = 1.0

  tglf_mass_in(1) = 1e-4
  tglf_mass_in(2) = 1.0

  ! Density ratios: ne/ne,ni(1)/ne,ni(2)/ne,
  tglf_as_in(:) = 1.0
  !
  ! Density gradients (e,i,z)
  tglf_rlns_in(:) = 1.0
  !
  ! Temperature gradients (e,i,z)
  tglf_rlts_in(:) = 3.0
  !
  ! Temperature ratios: Te/Te,Ti(1)/Te,Ti(2)/Te
  tglf_taus_in(:) = 1.0

end subroutine ptglf_init
