module EXPRO_locsim_interface

  integer :: n_species_exp
  
  real, parameter :: mass_deuterium  = 3.3452
  real, parameter :: temp_norm_fac   = 1602.2
  real, parameter :: charge_norm_fac = 1.6022

  real, dimension(:), allocatable :: rmin_exp

  real, dimension(:,:), allocatable :: temp_exp
  real, dimension(:,:), allocatable :: dens_exp
  real, dimension(:,:), allocatable :: dlntdr_exp
  real, dimension(:,:), allocatable :: dlnndr_exp
  real, dimension(:,:), allocatable :: sdlntdr_exp
  real, dimension(:,:), allocatable :: sdlnndr_exp

  real, dimension(:), allocatable :: gamma_e_exp
  real, dimension(:), allocatable :: gamma_p_exp
  real, dimension(:), allocatable :: mach_exp
  
  real, dimension(:,:,:), allocatable :: geo_yin_exp

  ! Local values

  real :: shift_loc
  real :: q_loc
  real :: s_loc
  real :: kappa_loc
  real :: delta_loc
  real :: zeta_loc
  real :: s_kappa_loc
  real :: s_delta_loc
  real :: s_zeta_loc
  real :: zmag_loc
  real :: dzmag_loc
  real :: gamma_e_loc
  real :: gamma_p_loc
  real :: mach_loc
  real :: rmaj_loc
  real :: rhos_loc
  real :: z_eff_loc
  real :: b_unit_loc
  real :: rho_norm_loc
  real :: psi_norm_loc
  real :: psi_a_loc
  
  real, dimension(9) :: dens_loc
  real, dimension(9) :: temp_loc
  real, dimension(9) :: dlnndr_loc
  real, dimension(9) :: dlntdr_loc
  real, dimension(9) :: sdlnndr_loc
  real, dimension(9) :: sdlntdr_loc

  integer :: geo_ny_loc
  real, dimension(:,:), allocatable :: geo_yin_loc

end module EXPRO_locsim_interface
