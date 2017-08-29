!-----------------------------------------------------------
! neo_harvest.f90
!
! PURPOSE:
!  subroutines to send neo input and output parameters
!  to harvest
!---------------------------------------------------------

module neo_harvest

  implicit none

  integer :: ierr
  character(len=65507) :: harvest_sendline
  character(len=255) :: harvest_tag
  character NUL
  parameter(NUL = char(0))
  character(len=2) :: NUM

  include 'harvest_lib.inc'

contains

  subroutine neo_harvest_input

    use neo_globals
    use neo_interface
    use EXPRO_interface

    implicit none

    ierr=init_harvest('Neo_jbs'//NUL,harvest_sendline,len(harvest_sendline))
    ierr=set_harvest_protocol('UDP'//NUL)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_n_energy_in"//NUL,neo_n_energy_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_n_xi_in"//NUL,neo_n_xi_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_n_theta_in"//NUL,neo_n_theta_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_n_radial_in"//NUL,neo_n_radial_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_matsz_scalefac_in"//NUL,neo_matsz_scalefac_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_rmin_over_a_in"//NUL,neo_rmin_over_a_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_rmin_over_a_2_in"//NUL,neo_rmin_over_a_2_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_rmaj_over_a_in"//NUL,neo_rmaj_over_a_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_silent_flag_in"//NUL,neo_silent_flag_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_sim_model_in"//NUL,neo_sim_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_equilibrium_model_in"//NUL,neo_equilibrium_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_collision_model_in"//NUL,neo_collision_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_profile_model_in"//NUL,neo_profile_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_profile_erad0_model_in"//NUL,neo_profile_erad0_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_profile_equilibrium_model_in"//NUL,neo_profile_equilibrium_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_ipccw_in"//NUL,neo_ipccw_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_btccw_in"//NUL,neo_btccw_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_te_ade_in"//NUL,neo_te_ade_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_ne_ade_in"//NUL,neo_ne_ade_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_dlntdre_ade_in"//NUL,neo_dlntdre_ade_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_dlnndre_ade_in"//NUL,neo_dlnndre_ade_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_rotation_model_in"//NUL,neo_rotation_model_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_omega_rot_in"//NUL,neo_omega_rot_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_omega_rot_deriv_in"//NUL,neo_omega_rot_deriv_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_spitzer_model_in"//NUL,neo_spitzer_model_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_epar0_spitzer_in"//NUL,neo_epar0_spitzer_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_coll_uncoupledei_model_in"//NUL,neo_coll_uncoupledei_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_coll_uncoupledaniso_model_in"//NUL,neo_coll_uncoupledaniso_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_n_species_in"//NUL,neo_n_species_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_nu_1_in"//NUL,neo_nu_1_in)
    ierr=set_harvest_payload_int_array(harvest_sendline,"neo_z_in"//NUL,neo_z_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_mass_in"//NUL,neo_mass_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_dens_in"//NUL,neo_dens_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_temp_in"//NUL,neo_temp_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_dlnndr_in"//NUL,neo_dlnndr_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_dlntdr_in"//NUL,neo_dlntdr_in,11)
    ierr=set_harvest_payload_int_array(harvest_sendline,"neo_aniso_model_in"//NUL,neo_aniso_model_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_temp_para_in"//NUL,neo_temp_para_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_dlntdr_para_in"//NUL,neo_dlntdr_para_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_temp_perp_in"//NUL,neo_temp_perp_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_dlntdr_perp_in"//NUL,neo_dlntdr_perp_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_profile_dlnndr_scale_in"//NUL,neo_profile_dlnndr_scale_in,11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_profile_dlntdr_scale_in"//NUL,neo_profile_dlntdr_scale_in,11)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_dphi0dr_in"//NUL,neo_dphi0dr_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_epar0_in"//NUL,neo_epar0_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_q_in"//NUL,neo_q_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_rho_star_in"//NUL,neo_rho_star_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_shear_in"//NUL,neo_shear_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_shift_in"//NUL,neo_shift_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_kappa_in"//NUL,neo_kappa_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_s_kappa_in"//NUL,neo_s_kappa_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_delta_in"//NUL,neo_delta_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_s_delta_in"//NUL,neo_s_delta_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_zeta_in"//NUL,neo_zeta_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_s_zeta_in"//NUL,neo_s_zeta_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_zmag_over_a_in"//NUL,neo_zmag_over_a_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_s_zmag_in"//NUL,neo_s_zmag_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_beta_star_in"//NUL,neo_beta_star_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_profile_delta_scale_in"//NUL,neo_profile_delta_scale_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_profile_zeta_scale_in"//NUL,neo_profile_zeta_scale_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_profile_zmag_scale_in"//NUL,neo_profile_zmag_scale_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_threed_model_in"//NUL,neo_threed_model_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_threed_exb_model_in"//NUL,neo_threed_exb_model_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_threed_exb_dphi0dr_in"//NUL,neo_threed_exb_dphi0dr_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_threed_drift_model_in"//NUL,neo_threed_drift_model_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_threed_hyperxi_in"//NUL,neo_threed_hyperxi_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_laguerre_method_in"//NUL,neo_laguerre_method_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_write_cmoments_flag_in"//NUL,neo_write_cmoments_flag_in)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_geo_ny_in"//NUL,neo_geo_ny_in)
  ! ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_geo_yin_in"//NUL,neo_geo_yin_in,)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_subroutine_flag"//NUL,neo_subroutine_flag)
    ierr=set_harvest_payload_int(harvest_sendline,"neo_test_flag_in"//NUL,neo_test_flag_in)

  end subroutine neo_harvest_input

  subroutine neo_harvest_output

    use neo_globals
    use neo_interface
    use EXPRO_interface

    implicit none

    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_pflux_thHH_out"//NUL,neo_th_out(1))
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_eflux_thHHi_out"//NUL,neo_th_out(2))
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_eflux_thHHe_out"//NUL,neo_th_out(3))
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_eflux_thCHi_out"//NUL,neo_th_out(4))
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_pflux_thHS_out"//NUL,neo_thHS_out(1:11,1),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_eflux_thHS_out"//NUL,neo_thHS_out(1:11,2),11)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_jpar_thS_out"//NUL,neo_th_out(5))
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_jpar_thK_out"//NUL,neo_th_out(6))
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_jpar_thN_out"//NUL,neo_th_out(7))
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_jtor_thS_out"//NUL,neo_th_out(8))
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_pflux_dke_out"//NUL,neo_dke_out(1:11,1),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_efluxtot_dke_out"//NUL,neo_dke_out(1:11,2),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_mflux_dke_out"//NUL,neo_dke_out(1:11,3),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_efluxncv_dke_out"//NUL,neo_dke_out(1:11,4),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_vpol_dke_out"//NUL,neo_dke_out(1:11,5),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_vtor_dke_out"//NUL,neo_dke_out(1:11,6),11)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_jpar_dke_out"//NUL,neo_dke_1d_out(1))
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_jtor_dke_out"//NUL,neo_dke_1d_out(2))
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_pflux_gv_out"//NUL,neo_gv_out(1:11,1),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_efluxtot_gv_out"//NUL,neo_gv_out(1:11,2),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_mflux_gv_out"//NUL,neo_gv_out(1:11,3),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_efluxncv_gv_out"//NUL,neo_gv_out(1:11,4),11)
    
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_pflux_nclass_out"//NUL,neo_nclass_out(1:11,1),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_efluxtot_nclass_out"//NUL,neo_nclass_out(1:11,2),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_vpol_nclass_out"//NUL,neo_nclass_out(1:11,3),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_vtor_nclass_out"//NUL,neo_nclass_out(1:11,4),11)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,"neo_nclassvis_out"//NUL,neo_nclass_out(1:11,5),11)
    ierr=set_harvest_payload_dbl(harvest_sendline,"neo_jpar_nclass_out"//NUL,neo_nclass_1d_out)
    
    ierr=set_harvest_payload_int(harvest_sendline,"neo_error_status_out"//NUL,neo_error_status_out)

    ierr=harvest_send(harvest_sendline)

  end subroutine neo_harvest_output

end module neo_harvest
