subroutine tgyro_etgcrit(eflux_e_loc,eflux_i_loc,pflux_e_loc,pflux_i_loc)

  use tgyro_globals

  implicit none

  real, intent(inout) :: eflux_e_loc
  real, intent(inout) :: eflux_i_loc
  real, intent(inout) :: pflux_e_loc
  real, intent(inout) :: pflux_i_loc
  
  real :: etae_in
  real :: alte_in
  real :: alne_in
  real :: chi_e_out
  real :: d_e_out

  ! Input parameters
  etae_in = dlntedr(i_r)/dlnnedr(i_r)
  alte_in = r_min*dlntedr(i_r)
  alne_in = r_min*dlnnedr(i_r)

  !----------------------------------------------
  ! Diffusivity model
  chi_e_out = 1.5*(etae_in-1.4)*(alte_in/60.0)
  d_e_out   = 0.02*chi_e_out
  !----------------------------------------------

  ! Conversion from diffusivity to flux
  eflux_e_loc = chi_e_out * alte_in
  pflux_e_loc =   d_e_out * alne_in
  eflux_i_loc = 0.0
  pflux_i_loc = pflux_e_loc

end subroutine tgyro_etgcrit
