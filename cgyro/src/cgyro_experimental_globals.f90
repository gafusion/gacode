module cgyro_experimental_globals

  integer :: n_grid_exp
  integer :: n_species_exp

  real, parameter :: mass_deuterium  = 3.3452
  real, parameter :: temp_norm_fac   = 1602.2
  real, parameter :: charge_norm_fac = 1.6022

  real :: dens_norm, temp_norm, vth_norm

  real    :: a_meters
  real    :: b_unit

  real, dimension(:), allocatable :: rmin_exp
  real, dimension(:), allocatable :: rmaj_exp
  real, dimension(:), allocatable :: q_exp
  real, dimension(:), allocatable :: s_exp
  real, dimension(:), allocatable :: shift_exp
  real, dimension(:), allocatable :: kappa_exp
  real, dimension(:), allocatable :: s_kappa_exp
  real, dimension(:), allocatable :: delta_exp
  real, dimension(:), allocatable :: s_delta_exp
  real, dimension(:), allocatable :: zeta_exp
  real, dimension(:), allocatable :: s_zeta_exp
  real, dimension(:), allocatable :: zmag_exp
  real, dimension(:), allocatable :: dzmag_exp

  real, dimension(:), allocatable :: te_ade_exp
  real, dimension(:), allocatable :: ne_ade_exp
  real, dimension(:), allocatable :: dlntdre_ade_exp
  real, dimension(:), allocatable :: dlnndre_ade_exp

  real, dimension(:,:), allocatable :: temp_exp
  real, dimension(:,:), allocatable :: dens_exp
  real, dimension(:,:), allocatable :: dlntdr_exp
  real, dimension(:,:), allocatable :: dlnndr_exp
  real, dimension(:,:), allocatable :: sdlntdr_exp
  real, dimension(:,:), allocatable :: sdlnndr_exp

  real, dimension(:), allocatable :: z_eff_exp
  real, dimension(:), allocatable :: b_unit_exp

  real, dimension(:), allocatable :: gamma_e_exp
  real, dimension(:), allocatable :: gamma_p_exp
  real, dimension(:), allocatable :: mach_exp
  real, dimension(:), allocatable :: rhos_exp

  real, dimension(:,:,:), allocatable :: geo_yin_exp

end module cgyro_experimental_globals


