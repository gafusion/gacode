subroutine qlgyro_cgyro_map

  use qlgyro_globals
  use qlgyro_cgyro_interface

  implicit none

  integer :: i

  n_species = cgyro_n_species_in
  n_ion = n_species - 1

  if (n_ion .eq. 0) n_ion = 1

  do i=1, n_species
     if (cgyro_z_in(i) .lt. 0) then
        elec_ind = i
        if (i .eq. 1) then
           ion_start = 2
           ion_end = n_species
        else
           ion_start = 1
           ion_end = n_species-1
        end if
     endif
  end do
  
  q = abs(cgyro_q_in)
  s = cgyro_s_in
  delta = cgyro_delta_in
  s_delta = cgyro_s_delta_in
  kappa = cgyro_kappa_in
  s_kappa = cgyro_s_kappa_in
  shift = cgyro_shift_in
  r_maj = cgyro_rmaj_in
  r_min = cgyro_rmin_in
  betae_unit = cgyro_betae_unit_in

  if (cgyro_ae_flag_in .eq. 0) then
     me = cgyro_mass_in(elec_ind)
     ne = cgyro_dens_in(elec_ind)
     te = cgyro_temp_in(elec_ind)
     dlnnedr = cgyro_dlnndr_in(elec_ind)
     dlntedr = cgyro_dlntdr_in(elec_ind)
  else
     me = cgyro_mass_ae_in
     ne = cgyro_dens_ae_in
     te = cgyro_temp_ae_in
     dlnnedr = cgyro_dlnndr_ae_in
     dlntedr = cgyro_dlntdr_ae_in
  end if
  ! Note that CGYRO normalizing mass is mD = 2*mp.

  zi_vec(1:n_ion) = cgyro_z_in(ion_start:ion_end)
  ni(1:n_ion) = cgyro_dens_in(ion_start:ion_end)
  ti(1:n_ion) = cgyro_temp_in(ion_start:ion_end)
  mi(1:n_ion) = cgyro_mass_in(ion_start:ion_end)
  dlnnidr(1:n_ion) = cgyro_dlnndr_in(ion_start:ion_end)
  dlntidr(1:n_ion) = cgyro_dlntdr_in(ion_start:ion_end)

  nue = cgyro_nu_ee_in

  lambda_star = cgyro_lambda_star_in
  
  z_eff  = cgyro_z_eff_in
  zmag   = cgyro_zmag_in
  dzmag  = cgyro_dzmag_in
  zeta   = cgyro_zeta_in
  s_zeta = cgyro_s_zeta_in
  shape_sin3 = cgyro_shape_sin3_in
  shape_ssin3 = cgyro_shape_s_sin3_in
  shape_cos0 = cgyro_shape_cos0_in
  shape_scos0 = cgyro_shape_s_cos0_in
  shape_cos1 = cgyro_shape_cos1_in
  shape_scos1 = cgyro_shape_s_cos1_in
  shape_cos2 = cgyro_shape_cos2_in
  shape_scos2 = cgyro_shape_s_cos2_in
  shape_cos3 = cgyro_shape_cos3_in
  shape_scos3 = cgyro_shape_s_cos3_in

  ! General geometry Fourier coefficients

  qlgyro_num_equil_flag = cgyro_equilibrium_model_in

  !rho_s = c_s* 

  if (qlgyro_rotation_flag == 1) then

     gamma_e = -cgyro_gamma_e_in
     gamma_p = -cgyro_gamma_p_in
     w0 = cgyro_mach_in

  else

     gamma_e = 0.0
     gamma_p = 0.0
     w0 = 0.0

  endif

  signb = -int(cgyro_btccw_in)
  signq = int(cgyro_btccw_in*cgyro_ipccw_in)

  bunit = cgyro_b_unit_in
  bgs2 = cgyro_b_gs2_in
  a_meters = cgyro_a_meters_in
  temp_norm = cgyro_temp_norm_in
  dens_norm = cgyro_dens_norm_in
  rho_star_norm = cgyro_rho_star_norm_in
  csda_norm = cgyro_vth_norm_in

end subroutine qlgyro_cgyro_map


!------------------------------------------------------------
! cgyro_qlgyro_map.f90
!
! PURPOSE:
!  Mapping from QLGYRO internal variables to CGYRO interface.
!------------------------------------------------------------

subroutine cgyro_qlgyro_map

  use qlgyro_globals
  use qlgyro_cgyro_interface

  implicit none

  integer :: i0
  real :: gamma_eb0

  !----------------------------------------------------------------
  ! Signs of toroidal magnetic field and current
  cgyro_btccw_in = -1.0*signb
  cgyro_btccw_in = -1.0*signb*signq
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Number of species (max=11)
  cgyro_n_species_in = n_species

  if (cgyro_n_species_in > 11) then
     call qlgyro_catch_error('ERROR: (qlgyro_cgyro_map) Too many ions in CGYRO.')
  endif
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Species loop:
  !
  ! Charges: e,i,z
  cgyro_z_in(elec_ind) = -1.0
  !
  ! Mass ratios: me/md,m(1)/md,m(2)/md, ...
  cgyro_mass_in(elec_ind) = me
  !
  ! Density ratios: ne/ne,ni(1)/ne,ni(2)/ne, ...
  cgyro_dens_in(elec_ind) = ne
  !
  ! Density gradients (e,i,z)
  cgyro_dlnndr_in(elec_ind) = dlnnedr
  !
  ! Temperature gradients (e,i,z)
  cgyro_dlntdr_in(elec_ind) = dlntedr
  !
  ! Temperature ratios: Te/Te,Ti(1)/Te,Ti(2)/Te
  cgyro_temp_in(elec_ind) = te

  cgyro_z_in(ion_start:ion_end) = zi_vec(1:n_ion)
  cgyro_dens_in(ion_start:ion_end) = ni(1:n_ion)
  cgyro_temp_in(ion_start:ion_end) = ti(1:n_ion)
  cgyro_mass_in(ion_start:ion_end) = mi(1:n_ion)
  cgyro_dlnndr_in(ion_start:ion_end) = dlnnidr(1:n_ion)
  cgyro_dlntdr_in(ion_start:ion_end) = dlntidr(1:n_ion)

  !----------------------------------------------------------------
  !   debye length/rhos   te in ev, rho_s in cm ne in 10^13/cm^3
  cgyro_lambda_star_in = lambda_star

  ! Model (Miller) shape
  cgyro_equilibrium_model_in = 2

  cgyro_rmin_in    = r_min
  cgyro_rmaj_in    = r_maj
  cgyro_zmag_in    = zmag
  cgyro_shift_in = shift
  cgyro_dzmag_in = dzmag
  cgyro_kappa_in   = kappa
  cgyro_s_kappa_in = s_kappa
  cgyro_delta_in   = delta
  cgyro_s_delta_in = s_delta
  cgyro_zeta_in    = zeta
  cgyro_s_zeta_in  = s_zeta
  cgyro_q_in       = q
  cgyro_s_in = s
  cgyro_shape_sin3_in = shape_sin3
  cgyro_shape_s_sin3_in = shape_ssin3
  cgyro_shape_cos0_in = shape_cos0
  cgyro_shape_s_cos0_in = shape_scos0
  cgyro_shape_cos1_in = shape_cos1
  cgyro_shape_s_cos1_in = shape_scos1
  cgyro_shape_cos2_in = shape_cos2
  cgyro_shape_s_cos2_in = shape_scos2
  cgyro_shape_cos3_in = shape_cos3
  cgyro_shape_s_cos3_in = shape_scos3

  ! What should this be?
  cgyro_beta_star_scale_in = 1.0

  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Electron beta used for electromagnetic calculations
  cgyro_betae_unit_in = betae_unit
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! Collisions:
  !
  ! Electron collision frequency
  cgyro_nu_ee_in = nue
  !
  ! Zeff
  cgyro_z_eff_in = z_eff
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Gamma_ExB (ExB shearing rate, units of a/cs)
  if (qlgyro_rotation_flag == 1) then

     cgyro_gamma_e_in  = -gamma_e
     cgyro_gamma_p_in = -gamma_p
     cgyro_mach_in = -w0
  else

     cgyro_gamma_e_in = 0.0
     cgyro_gamma_p_in = 0.0
     cgyro_mach_in = 0.0

  endif

  cgyro_delta_t_method_in = 1

  cgyro_b_unit_in = bunit
  cgyro_b_gs2_in = bgs2
  cgyro_a_meters_in = a_meters
  cgyro_temp_norm_in = temp_norm
  cgyro_dens_norm_in = dens_norm
  cgyro_rho_star_norm_in = rho_star_norm
  cgyro_vth_norm_in = csda_norm


end subroutine cgyro_qlgyro_map
