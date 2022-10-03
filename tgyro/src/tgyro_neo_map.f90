subroutine tgyro_neo_map

  use tgyro_globals
  use neo_interface

  implicit none
  real :: gamma_p0,u000
  integer :: is,i0
  real :: m_norm,t_norm,n_norm
  real :: sum_nz
  
  ! Initialize NEO
  call neo_init(paths(i_r-1),gyro_comm)

  neo_test_flag_in = gyrotest_flag
  
  ! Simulation mode (dke solve vs. analytic)
  if (tgyro_neo_method == 2) then
     ! Kinetic NEO and theory (no NCLASS)
     neo_sim_model_in = 2
  else
     ! Analytic theory
     neo_sim_model_in = 0
  endif

  ! Resolution 
  neo_n_energy_in = tgyro_neo_n_energy
  neo_n_xi_in     = tgyro_neo_n_xi
  neo_n_theta_in  = tgyro_neo_n_theta

  ! Geometry
  neo_equilibrium_model_in = 2
  neo_rmin_over_a_in       = r(i_r)/r_min
  neo_rmaj_over_a_in       = r_maj(i_r)/r_min
  neo_q_in                 = abs(q(i_r))
  neo_shear_in             = s(i_r)
  neo_shift_in             = shift(i_r)
  neo_kappa_in             = kappa(i_r)
  neo_s_kappa_in           = s_kappa(i_r)
  neo_delta_in             = delta(i_r)
  neo_s_delta_in           = s_delta(i_r)
  neo_zeta_in              = zeta(i_r)
  neo_s_zeta_in            = s_zeta(i_r)
  neo_zmag_over_a_in       = zmag(i_r)/r_min
  neo_s_zmag_in            = dzmag(i_r)

  ! New HAM shape parameters
  neo_shape_sin3_in        = shape_sin3(i_r)
  neo_shape_s_sin3_in      = shape_ssin3(i_r)
  neo_shape_cos0_in        = shape_cos0(i_r)
  neo_shape_s_cos0_in      = shape_scos0(i_r)
  neo_shape_cos1_in        = shape_cos1(i_r)
  neo_shape_s_cos1_in      = shape_scos1(i_r)
  neo_shape_cos2_in        = shape_cos2(i_r)
  neo_shape_s_cos2_in      = shape_scos2(i_r)
  neo_shape_cos3_in        = shape_cos3(i_r)
  neo_shape_s_cos3_in      = shape_scos3(i_r)
  
  neo_ipccw_in = -signb*signq
  neo_btccw_in = -signb

  neo_rho_star_in  = 0.001

  neo_n_species_in = sum(calc_flag(1:loc_n_ion))+1
 
  if (neo_n_species_in > 9) then
     call tgyro_catch_error('ERROR: (TGYRO) n_species > 9 not supported in NEO interface.') 
  endif

  ! Assuming NEO mass norm is m_deuterium and temp norm is based on first ion
  ! dens norm based on electrons
  m_norm = md
  t_norm = ti(1,i_r)
  n_norm = ne(i_r)
  
  ! Electrons
  neo_z_in(1)      = -1
  neo_mass_in(1)   = me*loc_me_multiplier/m_norm
  neo_dens_in(1)   = 1.0
  neo_temp_in(1)   = te(i_r)/t_norm
  neo_dlnndr_in(1) = r_min*dlnnedr(i_r)
  neo_dlntdr_in(1) = r_min*dlntedr(i_r)

  neo_nu_1_in      = nui(1,i_r)*r_min/c_s(i_r)/sqrt(t_norm/te(i_r)) &
       * ne(i_r)/ni(1,i_r) * sqrt(mi(1)/(me*loc_me_multiplier)) &
       * (ti(1,i_r)/te(i_r))**1.5 / (1.0*zi_vec(1))**4

  ! Ions
  i0 = 1
  do is=1,loc_n_ion
     if (calc_flag(is) == 0) cycle
     i0 = i0 + 1
     neo_z_in(i0)      = zi_vec(is)
     neo_mass_in(i0)   = mi(is)/m_norm
     neo_dens_in(i0)   = ni(is,i_r)/n_norm
     neo_temp_in(i0)   = ti(is,i_r)/t_norm
     neo_dlnndr_in(i0) = r_min*dlnnidr(is,i_r)
     neo_dlntdr_in(i0) = r_min*dlntidr(is,i_r)
  enddo

  ! Force quasineutrality
  if (loc_n_ion == 1) then
     ! Quasineutrality offset
     sum_nz = sum(neo_dens_in(:)*neo_z_in(:))
     ! Modify main ion to force quasineutrality
     neo_dens_in(2) = neo_dens_in(2)-sum_nz/neo_z_in(2)
  endif

  ! Setting density gradient artificially to zero to compute D and v
  if (tgyro_zero_dens_grad_flag /= 0) then
     neo_dlnndr_in(tgyro_zero_dens_grad_flag) = 0
  endif

  ! Rotation is always *on* in NEO.
  !
  ! COORDINATES: The signs of all rotation-related quantities below are 
  ! inherited (unchanged) from input.profiles.  In general NEO expects 
  ! these to be correctly signed/oriented.

  gamma_p0 = -r_maj(i_r)*w0p(i_r)
  u000     = r_maj(i_r)*w0(i_r)

  neo_rotation_model_in = 2
  neo_omega_rot_in = u000 * r_min / r_maj(i_r) &
       / (c_s(i_r) * sqrt(t_norm/te(i_r)))
  neo_omega_rot_deriv_in = -gamma_p0 * r_min**2 / r_maj(i_r) &
       / (c_s(i_r) * sqrt(t_norm/te(i_r)))

  ! Parameter only used for global runs
  neo_rmin_over_a_2_in = neo_rmin_over_a_in

end subroutine tgyro_neo_map
