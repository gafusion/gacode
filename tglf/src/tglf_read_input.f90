!--------------------------------------------------------------
! tglf_read_input.f90
!
! PURPOSE:
!  Complete read of tglf input parameters as specified by 
!  tglf_interface (but ignoring Fourier and ELITE geometry 
!  coefficients).
!--------------------------------------------------------------

subroutine tglf_read_input

  use tglf_interface

  implicit none

  open(unit=1,file=trim(tglf_path_in)//'input.tglf.gen',status='old')

  read(1,*) tglf_units_in
  read(1,*) tglf_use_transport_model_in
  read(1,*) tglf_geometry_flag_in
  read(1,*) tglf_write_wavefunction_flag_in

  ! Data passed to: put_signs
  read(1,*) tglf_sign_bt_in
  read(1,*) tglf_sign_it_in

  ! Data passed to: put_rare_switches
  read(1,*) tglf_theta_trapped_in
  read(1,*) tglf_park_in
  read(1,*) tglf_ghat_in
  read(1,*) tglf_gchat_in
  read(1,*) tglf_wd_zero_in
  read(1,*) tglf_linsker_factor_in
  read(1,*) tglf_gradb_factor_in
  read(1,*) tglf_filter_in
  read(1,*) tglf_damp_psi_in
  read(1,*) tglf_damp_sig_in

  ! Data passed to: put_switches
  read(1,*) tglf_iflux_in
  read(1,*) tglf_use_bper_in
  read(1,*) tglf_use_bpar_in
  read(1,*) tglf_use_mhd_rule_in
  read(1,*) tglf_use_bisection_in
  read(1,*) tglf_use_inboard_detrapped_in
  read(1,*) tglf_ibranch_in
  read(1,*) tglf_nmodes_in
  read(1,*) tglf_nbasis_max_in
  read(1,*) tglf_nbasis_min_in
  read(1,*) tglf_nxgrid_in
  read(1,*) tglf_nky_in

  ! Data passed to: put_model_parameters
  read(1,*) tglf_adiabatic_elec_in
  read(1,*) tglf_alpha_mach_in
  read(1,*) tglf_alpha_e_in
  read(1,*) tglf_alpha_p_in
  read(1,*) tglf_alpha_quench_in
  read(1,*) tglf_alpha_zf_in
  read(1,*) tglf_xnu_factor_in
  read(1,*) tglf_debye_factor_in
  read(1,*) tglf_etg_factor_in
  read(1,*) tglf_sat_rule_in
  read(1,*) tglf_kygrid_model_in
  read(1,*) tglf_xnu_model_in
  read(1,*) tglf_vpar_model_in
  read(1,*) tglf_vpar_shear_model_in

  ! Data passed to: put_species
  read(1,*) tglf_ns_in
  read(1,*) tglf_mass_in(1)
  read(1,*) tglf_mass_in(2)
  read(1,*) tglf_mass_in(3)
  read(1,*) tglf_mass_in(4)
  read(1,*) tglf_mass_in(5)
  read(1,*) tglf_mass_in(6)
  read(1,*) tglf_zs_in(1)
  read(1,*) tglf_zs_in(2)
  read(1,*) tglf_zs_in(3)
  read(1,*) tglf_zs_in(4)
  read(1,*) tglf_zs_in(5)
  read(1,*) tglf_zs_in(6)

  ! Data passed to: put_kys
  read(1,*) tglf_ky_in

  ! Data passed to: put_gaussian_width
  read(1,*) tglf_width_in
  read(1,*) tglf_width_min_in
  read(1,*) tglf_nwidth_in
  read(1,*) tglf_find_width_in

  ! Data passed to: put_gradients
  read(1,*) tglf_rlns_in(1)
  read(1,*) tglf_rlns_in(2)
  read(1,*) tglf_rlns_in(3)
  read(1,*) tglf_rlns_in(4)
  read(1,*) tglf_rlns_in(5)
  read(1,*) tglf_rlns_in(6)
  read(1,*) tglf_rlts_in(1)
  read(1,*) tglf_rlts_in(2)
  read(1,*) tglf_rlts_in(3)
  read(1,*) tglf_rlts_in(4)
  read(1,*) tglf_rlts_in(5)
  read(1,*) tglf_rlts_in(6)
  read(1,*) tglf_vpar_shear_in(1) 
  read(1,*) tglf_vpar_shear_in(2) 
  read(1,*) tglf_vpar_shear_in(3) 
  read(1,*) tglf_vpar_shear_in(4) 
  read(1,*) tglf_vpar_shear_in(5) 
  read(1,*) tglf_vpar_shear_in(6) 
  read(1,*) tglf_vexb_shear_in

  ! Data passed to: put_profile_shear
  read(1,*) tglf_vns_shear_in(1)
  read(1,*) tglf_vns_shear_in(2)
  read(1,*) tglf_vns_shear_in(3)
  read(1,*) tglf_vns_shear_in(4)
  read(1,*) tglf_vns_shear_in(5)
  read(1,*) tglf_vns_shear_in(6)
  read(1,*) tglf_vts_shear_in(1)
  read(1,*) tglf_vts_shear_in(2)
  read(1,*) tglf_vts_shear_in(3)
  read(1,*) tglf_vts_shear_in(4)
  read(1,*) tglf_vts_shear_in(5)
  read(1,*) tglf_vts_shear_in(6)

  ! Data passed to: put_averages
  read(1,*) tglf_taus_in(1)
  read(1,*) tglf_taus_in(2)
  read(1,*) tglf_taus_in(3)
  read(1,*) tglf_taus_in(4)
  read(1,*) tglf_taus_in(5)
  read(1,*) tglf_taus_in(6)
  read(1,*) tglf_as_in(1)
  read(1,*) tglf_as_in(2)
  read(1,*) tglf_as_in(3)
  read(1,*) tglf_as_in(4)
  read(1,*) tglf_as_in(5)
  read(1,*) tglf_as_in(6)
  read(1,*) tglf_vpar_in(1)
  read(1,*) tglf_vpar_in(2)
  read(1,*) tglf_vpar_in(3)
  read(1,*) tglf_vpar_in(4)
  read(1,*) tglf_vpar_in(5)
  read(1,*) tglf_vpar_in(6)
  read(1,*) tglf_vexb_in
  read(1,*) tglf_betae_in
  read(1,*) tglf_xnue_in
  read(1,*) tglf_zeff_in
  read(1,*) tglf_debye_in

  ! Data passed to: put_eikonal
  read(1,*) tglf_new_eikonal_in

  ! Data passed to: put_s_alpha_geometry
  read(1,*) tglf_rmin_sa_in
  read(1,*) tglf_rmaj_sa_in
  read(1,*) tglf_q_sa_in
  read(1,*) tglf_shat_sa_in
  read(1,*) tglf_alpha_sa_in
  read(1,*) tglf_xwell_sa_in
  read(1,*) tglf_theta0_sa_in
  read(1,*) tglf_b_model_sa_in
  read(1,*) tglf_ft_model_sa_in

  ! Data passed to: put_Miller_geometry
  read(1,*) tglf_rmin_loc_in
  read(1,*) tglf_rmaj_loc_in
  read(1,*) tglf_zmaj_loc_in
  read(1,*) tglf_drmindx_loc_in
  read(1,*) tglf_drmajdx_loc_in
  read(1,*) tglf_dzmajdx_loc_in
  read(1,*) tglf_kappa_loc_in
  read(1,*) tglf_s_kappa_loc_in
  read(1,*) tglf_delta_loc_in
  read(1,*) tglf_s_delta_loc_in
  read(1,*) tglf_zeta_loc_in
  read(1,*) tglf_s_zeta_loc_in
  read(1,*) tglf_q_loc_in
  read(1,*) tglf_q_prime_loc_in
  read(1,*) tglf_p_prime_loc_in
  read(1,*) tglf_kx0_loc_in

  ! Set threshold for TGLF-NN execution versus full TGLF calculation
  read(1,*) tglf_nn_max_error_in

  close(1)

end subroutine tglf_read_input
