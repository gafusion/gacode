subroutine tgyro_gyro_map

  use tgyro_globals
  use gyro_interface

  implicit none

  real :: gamma_eb0
  real :: gamma_p0
  real :: u000

  if (loc_n_ion > 5) then
     call tgyro_catch_error('ERROR: (TGYRO) n_ion > 5 not supported in GYRO interface.')
  endif

  ! Initialize GYRO
  call gyro_init(paths(i_r-1), gyro_comm)

  if (worker == 0) then
     gyro_output_flag_in = 1
  else
     gyro_output_flag_in = 0
  endif

  ! Force TGYRO to handle profile parameter calculation
  gyro_radial_profile_method_in = 5

  gyro_safety_factor_in = abs(q(i_r))
  gyro_shear_in = s(i_r)
  gyro_delta_in = delta(i_r)
  gyro_s_delta_in = s_delta(i_r)
  gyro_kappa_in = kappa(i_r)
  gyro_s_kappa_in = s_kappa(i_r)
  gyro_shift_in = shift(i_r)
  gyro_aspect_ratio_in = r_maj(i_r)/r_min
  gyro_radius_in = r(i_r)/r_min
  gyro_betae_unit_in = betae_unit(i_r)*loc_betae_scale
  if (tgyro_ptot_flag == 1) then
     error_flag = 1
     error_msg  = 'Error: GYRO cannot use ptot via TGYRO.'
  endif

  gyro_mu_electron_in = sqrt(mi(1)/(me*loc_me_multiplier))
  gyro_nu_ei_in = nue(i_r)*r_min/c_s(i_r)*loc_nu_scale
  gyro_dlnndr_electron_in = r_min*dlnnedr(i_r)
  gyro_dlntdr_electron_in = r_min*dlntedr(i_r)

  ! Note that GYRO normalizing mass is mD = 2*mp.

  gyro_z_in          = zi_vec(1)
  gyro_ni_over_ne_in = ni(1,i_r)/ne(i_r)
  gyro_ti_over_te_in = ti(1,i_r)/te(i_r)
  gyro_mu_1_in       = sqrt(2*mp/mi(1))
  gyro_dlnndr_in     = r_min*dlnnidr(1,i_r)
  gyro_dlntdr_in     = r_min*dlntidr(1,i_r)

  if (loc_n_ion > 1) then
     gyro_z_2_in          = zi_vec(2)
     gyro_ni_over_ne_2_in = ni(2,i_r)/ne(i_r)
     gyro_ti_over_te_2_in = ti(2,i_r)/te(i_r)
     gyro_mu_2_in         = sqrt(2*mp/mi(2))
     gyro_dlnndr_2_in     = r_min*dlnnidr(2,i_r)
     gyro_dlntdr_2_in     = r_min*dlntidr(2,i_r)
  endif

  if (loc_n_ion > 2) then
     gyro_z_3_in          = zi_vec(3)
     gyro_ni_over_ne_3_in = ni(3,i_r)/ne(i_r)
     gyro_ti_over_te_3_in = ti(3,i_r)/te(i_r)
     gyro_mu_3_in         = sqrt(2*mp/mi(3))
     gyro_dlnndr_3_in     = r_min*dlnnidr(3,i_r)
     gyro_dlntdr_3_in     = r_min*dlntidr(3,i_r)
  endif

  if (loc_n_ion > 3) then
     gyro_z_4_in          = zi_vec(4)
     gyro_ni_over_ne_4_in = ni(4,i_r)/ne(i_r)
     gyro_ti_over_te_4_in = ti(4,i_r)/te(i_r)
     gyro_mu_4_in         = sqrt(2*mp/mi(4))
     gyro_dlnndr_4_in     = r_min*dlnnidr(4,i_r)
     gyro_dlntdr_4_in     = r_min*dlntidr(4,i_r)
  endif

  if (loc_n_ion > 4) then
     gyro_z_5_in          = zi_vec(5)
     gyro_ni_over_ne_5_in = ni(5,i_r)/ne(i_r)
     gyro_ti_over_te_5_in = ti(5,i_r)/te(i_r)
     gyro_mu_5_in         = sqrt(2*mp/mi(5))
     gyro_dlnndr_5_in     = r_min*dlnnidr(5,i_r)
     gyro_dlntdr_5_in     = r_min*dlntidr(5,i_r)
  endif

  gyro_z_eff_in  = z_eff(i_r)
  gyro_zmag_in   = zmag(i_r)
  gyro_dzmag_in  = dzmag(i_r)
  gyro_zeta_in   = zeta(i_r)
  gyro_s_zeta_in = s_zeta(i_r)

  ! General geometry Fourier coefficients

  gyro_num_equil_flag_in = loc_num_equil_flag

  gyro_n_fourier_geo_in = n_fourier_geo
  gyro_a_fourier_geo_in(:,:) = a_fourier_geo(:,:,i_r)


  if (tgyro_rotation_flag == 1) then

     ! COORDINATES: The signs of all rotation-related quantities below are 
     ! inherited (unchanged) from input.profiles.  In general GYRO expects 
     ! these to be correctly signed/oriented.

     u000      = r_maj(i_r)*w0(i_r)
     gamma_p0  = -r_maj(i_r)*w0p(i_r)
     gamma_eb0 = gamma_p0*r(i_r)/(q(i_r)*r_maj(i_r)) 

     gyro_gamma_e_in = gamma_eb0*r_min/c_s(i_r)
     gyro_gamma_p_in = gamma_p0*r_min/c_s(i_r)
     gyro_mach_in    = u000/c_s(i_r)

  else

     gyro_gamma_e_in = 0.0
     gyro_gamma_p_in = 0.0
     gyro_mach_in = 0.0

  endif

  gyro_ipccw_in = -real(signb*signq)
  gyro_btccw_in = -real(signb)

end subroutine tgyro_gyro_map
