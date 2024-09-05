
!------------------------------------------------------------
! tgyro_cgyro_map.f90
!
! PURPOSE:
!  Mapping from TGYRO internal variables to CGYRO interface.
!------------------------------------------------------------

subroutine tgyro_cgyro_map

  use tgyro_globals
  use qlgyro_cgyro_interface

  implicit none

  integer :: i0, i_ion
  real :: gamma_eb0
  real :: gamma_p0

  integer :: qn_spec
  real :: delta_n, delta_dn
  ! call cgyro_init(lpath, gyro_comm)

  cgyro_profile_model_in = 1

  !----------------------------------------------------------------
  ! Signs of toroidal magnetic field and current 
  cgyro_btccw_in = -1.0*signb
  cgyro_btccw_in = -1.0*signb*signq
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Number of species (max=11)
  cgyro_n_species_in = sum(calc_flag(1:loc_n_ion))+1

  if (cgyro_n_species_in > 11) then
     call tgyro_catch_error('ERROR: (tgyro_cgyro_map) Too many ions in CGYRO.')
  endif
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Species loop:
  !
  ! Quasieutrality enforcing
  delta_n = 0.0
  delta_dn = 0.0
  qn_spec = -1

  ! Charges: e,i,z
  cgyro_z_in(1) = -1.0
  !
  ! Mass ratios: me/md,m(1)/md,m(2)/md, ... 
  !              [assume md is normalizing mass]
  cgyro_mass_in(1) = (me*loc_me_multiplier)/md
  !
  ! Density ratios: ne/ne,ni(1)/ne,ni(2)/ne, ...
  cgyro_dens_in(1) = 1.0
  !
  ! Density gradients (e,i,z)
  cgyro_dlnndr_in(1) = r_min*dlnnedr(i_r)
  !
  ! Temperature gradients (e,i,z)
  cgyro_dlntdr_in(1) = r_min*dlntedr(i_r)
  !
  ! Temperature ratios: Te/Te,Ti(1)/Te,Ti(2)/Te
  cgyro_temp_in(1) = 1.0
  if (evo_e(0) /= -1) then
     delta_n = delta_n - 1.0
     delta_dn = delta_dn - 1.0*r_min*dlnnedr(i_r)
  else
     qn_spec = 1
  endif

  i0 = 1

  do i_ion=1,loc_n_ion
     if (calc_flag(i_ion) == 0) cycle
     i0 = i0+1 
     cgyro_z_in(i0)   = zi_vec(i_ion)
     cgyro_mass_in(i0) = mi_vec(i_ion)
     cgyro_dens_in(i0)   = ni(i_ion,i_r)/ne(i_r) 
     cgyro_dlnndr_in(i0) = r_min*dlnnidr(i_ion,i_r)
     cgyro_dlntdr_in(i0) = r_min*dlntidr(i_ion,i_r)
     cgyro_temp_in(i0) = ti(i_ion,i_r)/te(i_r)

     if (evo_e(i_ion) /= -1) then
        delta_n = delta_n + zi_vec(i_ion) * ni(i_ion, i_r) / ne(i_r)
        delta_dn = delta_dn + zi_vec(i_ion) * ni(i_ion, i_r) / ne(i_r) * r_min * dlnnidr(i_ion,i_r)
     else
        qn_spec = i0
     end if
  end do

  ! Force quasineutrality
  if (qn_spec .gt. 0) then
     cgyro_dens_in(qn_spec) = -delta_n
     cgyro_dlnndr_in(qn_spec) = delta_dn / delta_n
  end if

  ! Setting density gradient artificially to zero to compute D and v
  if (tgyro_zero_dens_grad_flag /= 0) then
     cgyro_dlnndr_in(tgyro_zero_dens_grad_flag) = 0
  endif

  ! how to handle scaling of den, ti, te
  ! cgyro_dlntdr_scale_in = tgyro_input_dlntdr_scale

  !----------------------------------------------------------------
  !   debye length/rhos   te in ev, rho_s in cm ne in 10^13/cm^3
  cgyro_lambda_star_in = 7.43e2*sqrt(te(i_r)/(ne(i_r)))/abs(rho_s(i_r))

  ! Model (Miller) shape
  cgyro_equilibrium_model_in = 2

  cgyro_rmin_in    = r(i_r)/r_min
  cgyro_rmaj_in    = r_maj(i_r)/r_min
  cgyro_zmag_in    = zmag(i_r)/r_min
  cgyro_shift_in = shift(i_r)
  cgyro_dzmag_in = dzmag(i_r)
  cgyro_kappa_in   = kappa(i_r)
  cgyro_s_kappa_in = s_kappa(i_r)
  cgyro_delta_in   = delta(i_r)
  cgyro_s_delta_in = s_delta(i_r)
  cgyro_zeta_in    = zeta(i_r)
  cgyro_s_zeta_in  = s_zeta(i_r)
  cgyro_q_in       = q(i_r)
  cgyro_s_in = s(i_r)
  cgyro_shape_sin3_in = shape_sin3(i_r)
  cgyro_shape_s_sin3_in = shape_ssin3(i_r)
  cgyro_shape_cos0_in = shape_cos0(i_r)
  cgyro_shape_s_cos0_in = shape_scos0(i_r)
  cgyro_shape_cos1_in = shape_cos1(i_r)
  cgyro_shape_s_cos1_in = shape_scos1(i_r)
  cgyro_shape_cos2_in = shape_cos2(i_r)
  cgyro_shape_s_cos2_in = shape_scos2(i_r)
  cgyro_shape_cos3_in = shape_cos3(i_r)
  cgyro_shape_s_cos3_in = shape_scos3(i_r)
  
  ! What should this be?
  cgyro_beta_star_scale_in = 1.0

  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Electron beta used for electromagnetic calculations
  cgyro_betae_unit_in = betae_unit(i_r)
  cgyro_betae_unit_scale_in = loc_betae_scale
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! Collisions:
  !
  ! Electron collision frequency
  cgyro_nu_ee_in = nue(i_r)*r_min/c_s(i_r)
  cgyro_nu_ee_scale_in = loc_nu_scale
  !
  ! Zeff
  cgyro_z_eff_in = z_eff(i_r)
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Gamma_ExB (ExB shearing rate, units of a/cs)
  if (tgyro_rotation_flag == 1) then
     gamma_p0  = -r_maj(i_r)*f_rot(i_r)*w0_norm
     gamma_eb0 = gamma_p0*r(i_r)/(q(i_r)*r_maj(i_r)) 
     ! Currently TGLF uses toroidal current as reference direction
     ! Overall minus sign is due to TGYRO toroidal angle in CW direction
     cgyro_gamma_e_in    = -gamma_eb0*r_min/c_s(i_r)  
     cgyro_gamma_p_in = -gamma_p0*r_min/c_s(i_r)
     cgyro_mach_in       = -r_maj(i_r)*w0(i_r)/c_s(i_r)
  else

     cgyro_gamma_e_in = 0.0
     cgyro_gamma_p_in = 0.0
     cgyro_mach_in = 0.0

  endif

  ! cgyro_mach_scale_in = tgyro_w0_scale


  !----------------------------------------------------------------
  ! CONTROL PARAMETERS
  !
  if (loc_betae_scale <= 0.0) then
     cgyro_n_field_in = 1
  endif

  cgyro_delta_t_method_in = 1

  cgyro_b_unit_in = b_unit(i_r)
  cgyro_a_meters_in = r_min
  cgyro_temp_norm_in = te(i_r)
  cgyro_dens_norm_in = ne(i_r)
  cgyro_rho_star_norm_in = rho_s(i_r) / r_min
  cgyro_vth_norm_in = c_s(i_r)

  cgyro_worker_index_in = worker_index
  cgyro_iteration_method_in = tgyro_iteration_method

end subroutine tgyro_cgyro_map
