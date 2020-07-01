subroutine tgyro_allocate_globals

  use tgyro_globals

  implicit none

  ! particle fluxes
  allocate(pflux_e_neo(n_r))
  allocate(pflux_i_neo(loc_n_ion,n_r))
  allocate(pflux_e_tur(n_r))
  allocate(pflux_i_tur(loc_n_ion,n_r))
  allocate(pflux_e_tot(n_r))
  allocate(pflux_i_tot(loc_n_ion,n_r))
  allocate(pflux_e_target(n_r))
  allocate(pflux_he_target(n_r))

  ! momentum fluxes
  allocate(mflux_e_neo(n_r))
  allocate(mflux_i_neo(loc_n_ion,n_r))
  allocate(mflux_e_tur(n_r))
  allocate(mflux_i_tur(loc_n_ion,n_r))
  allocate(mflux_tot(n_r))
  allocate(mflux_target(n_r))

  ! energy fluxes
  allocate(eflux_e_neo(n_r))
  allocate(eflux_i_neo(loc_n_ion,n_r))
  allocate(eflux_e_tur(n_r))
  allocate(eflux_i_tur(loc_n_ion,n_r))
  allocate(eflux_e_tot(n_r))
  allocate(eflux_i_tot(n_r))
  allocate(eflux_e_target(n_r))
  allocate(eflux_i_target(n_r))

  ! exchange power densities
  allocate(expwd_e_tur(n_r))
  allocate(expwd_i_tur(loc_n_ion,n_r))

  ! Moving gyroBohm diffusivity
  allocate(chi_gb(n_r))

  ! Moving gyroBohm fluxes
  allocate(gamma_gb(n_r))
  allocate(pi_gb(n_r))
  allocate(q_gb(n_r))
  allocate(s_gb(n_r))

  ! Collision frequencies
  allocate(nue(n_r))
  allocate(nui(loc_n_ion,n_r))
  allocate(nui_HH(n_r))
  allocate(nue_HH(n_r))
  allocate(nue_star(n_r))

  allocate(z_eff(n_r))

  ! Formulary exchange rate
  allocate(nu_exch(n_r))

  ! Alpha heating parameters
  allocate(frac_ai(n_r))
  allocate(frac_ae(n_r))
  allocate(e_cross(n_r))

  ! Electron and ion temperatures
  allocate(te(n_r))
  allocate(dlntedr(n_r))
  allocate(ti(loc_n_ion,n_r))
  allocate(dlntidr(loc_n_ion,n_r))

  ! Electron and ion densities
  allocate(ne(n_r))
  allocate(dlnnedr(n_r))
  allocate(ni(loc_n_ion,n_r))
  allocate(dlnnidr(loc_n_ion,n_r))

  ! Rotation quantities
  allocate(w0(n_r))
  allocate(w0p(n_r))
  allocate(gamma_eb(n_r))
  allocate(gamma_p(n_r))
  allocate(u00(n_r))

  allocate(pr(n_r))
  allocate(ptot(n_r))
  allocate(pext(n_r))
  allocate(dpext(n_r))
  allocate(dlnpdr(n_r))
  allocate(dlnptotdr(n_r))
  allocate(beta_unit(n_r))
  allocate(betae_unit(n_r))
  allocate(c_s(n_r))
  allocate(v_i(n_r))
  allocate(er(n_r))
  allocate(f_rot(n_r))

  allocate(rho(n_r))
  allocate(r(n_r))
  allocate(r_maj(n_r))
  allocate(q(n_r))
  allocate(s(n_r))
  allocate(kappa(n_r))
  allocate(s_kappa(n_r))
  allocate(delta(n_r))
  allocate(s_delta(n_r))
  allocate(shift(n_r))
  allocate(zmag(n_r))
  allocate(dzmag(n_r))
  allocate(zeta(n_r))
  allocate(s_zeta(n_r))
  allocate(b_ref(n_r))
  allocate(b_unit(n_r))
  allocate(volp(n_r))
  allocate(vol(n_r))
  allocate(ave_grad_r(n_r))
  allocate(rho_i(n_r))
  allocate(rho_s(n_r))

  allocate(p_i_in(n_r))
  allocate(p_e_in(n_r))
  allocate(p_exch_in(n_r))
  allocate(p_i_fus(n_r))
  allocate(p_e_fus(n_r))
  allocate(p_i_fus_in(n_r))
  allocate(p_e_fus_in(n_r))
  allocate(p_i_aux_in(n_r))
  allocate(p_e_aux_in(n_r))
  allocate(p_i(n_r))
  allocate(p_e(n_r))
  allocate(s_alpha_i(n_r))
  allocate(s_alpha_e(n_r))
  allocate(sn_alpha(n_r))
  allocate(s_brem(n_r))
  allocate(p_brem(n_r))
  allocate(s_sync(n_r))
  allocate(p_sync(n_r))
  allocate(p_line(n_r))
  allocate(p_sync_in(n_r))
  allocate(p_brem_in(n_r))
  allocate(p_line_in(n_r))
  allocate(s_exch(n_r))
  allocate(p_exch(n_r))
  allocate(s_expwd(n_r))
  allocate(p_expwd(n_r))
  allocate(f_b_in(n_r))
  allocate(f_w_in(n_r))
  allocate(f_he_fus(n_r))
  allocate(mf_in(n_r))

  allocate(res(p_max))
  allocate(res0(p_max))
  allocate(relax(p_max))
  allocate(therm_vec(sum(therm_flag(1:loc_n_ion))))

  allocate(a_fourier_geo(8,0:32,n_r))

  allocate(b_flag(p_max))

  allocate(pmap(2:n_r,n_evolve))
  allocate(mask(p_max,4+loc_n_ion))

end subroutine tgyro_allocate_globals
