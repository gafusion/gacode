subroutine qlgyro_gyro_map

  use qlgyro_globals
  use gyro_interface

  implicit none

  if (n_ion > 5) then
     call qlgyro_catch_error('ERROR: (QLGYRO) n_ion > 5 not supported in GYRO interface.')
  endif

  a_min = 1
  c_s = 1
  ne = 1.0
  te = 1.0

  !Find no. of ions
  if (gyro_ni_over_ne_2_in == 0) then
     n_ion = 1
  else if (gyro_ni_over_ne_3_in == 0) then
     n_ion = 2
  else if (gyro_ni_over_ne_4_in == 0) then
     n_ion = 3
  else if (gyro_ni_over_ne_5_in == 0) then
     n_ion = 4
  else
     n_ion = 5
  end if

  n_species = n_ion+1
  
  q = abs(gyro_safety_factor_in)
  s = gyro_shear_in
  delta = gyro_delta_in
  s_delta = gyro_s_delta_in
  kappa = gyro_kappa_in
  s_kappa = gyro_s_kappa_in
  shift = gyro_shift_in
  r_maj = gyro_aspect_ratio_in * a_min
  r_min = gyro_radius_in * a_min
  betae_unit = gyro_betae_unit_in

  ! Note that GYRO normalizing mass is mD = 2*mp.

  zi_vec(1)  = gyro_z_in
  ni(1) = gyro_ni_over_ne_in * ne
  ti(1) = gyro_ti_over_te_in * te
  mi(1) = 2*mp/gyro_mu_1_in**2
  dlnnidr(1) = gyro_dlnndr_in / a_min
  dlntidr(1) = gyro_dlntdr_in / a_min
  
  if (n_ion > 1) then
     zi_vec(2)  = gyro_z_2_in
     ni(2) = gyro_ni_over_ne_2_in * ne
     ti(2) = gyro_ti_over_te_2_in * te
     mi(2) = 2*mp/gyro_mu_2_in**2
     dlnnidr(2) = gyro_dlnndr_2_in / a_min
     dlntidr(2) = gyro_dlntdr_2_in / a_min
  endif

  if (n_ion > 2) then
     zi_vec(3)  = gyro_z_3_in
     ni(3) = gyro_ni_over_ne_3_in * ne
     ti(3) = gyro_ti_over_te_3_in * te
     mi(3) = 2*mp/gyro_mu_3_in**2
     dlnnidr(3) = gyro_dlnndr_3_in / a_min
     dlntidr(3) = gyro_dlntdr_3_in / a_min
  endif

  if (n_ion > 3) then
     zi_vec(4)  = gyro_z_4_in
     ni(4) = gyro_ni_over_ne_4_in * ne
     ti(4) = gyro_ti_over_te_4_in * te
     mi(4) = 2*mp/gyro_mu_4_in**2
     dlnnidr(4) = gyro_dlnndr_4_in / a_min
     dlntidr(4) = gyro_dlntdr_4_in / a_min
  endif

  if (n_ion > 4) then
     zi_vec(5)  = gyro_z_5_in
     ni(5) = gyro_ni_over_ne_5_in * ne
     ti(5) = gyro_ti_over_te_5_in * te
     mi(5) = 2*mp/gyro_mu_5_in**2
     dlnnidr(5) = gyro_dlnndr_5_in / a_min
     dlntidr(5) = gyro_dlntdr_5_in / a_min
  endif

  me =  mi(1)/ gyro_mu_electron_in**2
  nue = gyro_nu_ei_in * c_s/a_min
  dlnnedr = gyro_dlnndr_electron_in / a_min
  dlntedr = gyro_dlntdr_electron_in / a_min

  z_eff  = gyro_z_eff_in
  zmag   = gyro_zmag_in
  dzmag  = gyro_dzmag_in
  zeta   = gyro_zeta_in
  s_zeta = gyro_s_zeta_in

  ! General geometry Fourier coefficients

  qlgyro_num_equil_flag = gyro_num_equil_flag_in

  n_fourier_geo = gyro_n_fourier_geo_in
  a_fourier_geo(:,:) = gyro_a_fourier_geo_in(:,:)

  !rho_s = c_s* 

  if (qlgyro_rotation_flag == 1) then

     ! COORDINATES: The signs of all rotation-related quantities below are 
     ! inherited (unchanged) from input.profiles.  In general GYRO expects 
     ! these to be correctly signed/oriented.
     ! NEED TO DOUBLE CHECK ALL OF THIS

     !gamma_eb0 = gyro_gamma_e_in * c_s / a_min
     !gamma_p0  = gyro_gamma_p_in * c_s / a_min

     ! Not normalisedd correcttly as a_min set to 1
     !w0        = gyro_mach_in*c_s/r_maj

     ! Set flow terms in input.qlgyro

     !gamma_eb0 = 0.0
     gamma_p0 = 0.0
     w0 = 0.0
     

  else

     !gamma_eb0 = 0.0
     gamma_p0 = 0.0
     w0 = 0.0

  endif

  signb = -int(gyro_btccw_in)
  signq = int(gyro_btccw_in*gyro_ipccw_in)

  bunit = gyro_b_unit_norm_in
  a_meters = gyro_a_meters_in
  temp_norm = gyro_tem_norm_in
  dens_norm = gyro_den_norm_in
  rho_star_norm = gyro_rhos_norm_in
  csda_norm = gyro_csda_norm_in

end subroutine qlgyro_gyro_map


