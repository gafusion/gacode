module gyro_profile_exp

  !--------------------------------------------------
  ! Experimental grid range: 
  !
  !          i_exp = (0,1,...,n_grid_exp)
  !
  integer :: i_exp
  !--------------------------------------------------

  !--------------------------------------------------
  ! (keV/sec)/MW:
  !
  real, parameter :: kevdsecpmw = 1.6022e-19*1e3*1e-6
  !--------------------------------------------------

  real, dimension(:), allocatable :: rhogrid_exp
  real, dimension(:), allocatable :: rmin_exp
  real, dimension(:), allocatable :: rmaj_exp
  real, dimension(:), allocatable :: q_exp
  real, dimension(:), allocatable :: kappa_exp
  real, dimension(:), allocatable :: delta_exp
  real, dimension(:), allocatable :: zeta_exp

  real, dimension(:,:), allocatable :: tem_exp
  real, dimension(:,:), allocatable :: den_exp

  real, dimension(:), allocatable :: w0_exp
  real, dimension(:), allocatable :: w0p_exp
  real, dimension(:), allocatable :: gamma_e_exp
  real, dimension(:), allocatable :: gamma_p_exp
  real, dimension(:), allocatable :: mach_exp
  real, dimension(:), allocatable :: z_eff_exp
  real, dimension(:), allocatable :: zmag_exp
  real, dimension(:), allocatable :: ptot_exp

  real, dimension(:), allocatable :: r_p
  real, dimension(:), allocatable :: b_unit_p
  real, dimension(:), allocatable :: rhosda_p
  real, dimension(:), allocatable :: csda_p
  real, dimension(:), allocatable :: drmaj_p
  real, dimension(:), allocatable :: dzmag_p
  real, dimension(:), allocatable :: s_delta_p
  real, dimension(:), allocatable :: s_zeta_p
  real, dimension(:), allocatable :: s_kappa_p
  real, dimension(:), allocatable :: q_p
  real, dimension(:), allocatable :: shat_p

  real, dimension(:,:), allocatable :: dlnndr_p
  real, dimension(:,:), allocatable :: dlntdr_p
  real, dimension(:), allocatable :: dlnptotdr_p

  real, dimension(:), allocatable :: beta_unit_p
  real, dimension(:), allocatable :: beta_unit_ptot_p

  real, dimension(:,:,:), allocatable :: geo_p

end module gyro_profile_exp


