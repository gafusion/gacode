!-------------------------------------------------------------------------
! cgyro_interface.f90
!
! PURPOSE:
!  Provides interface description for CGYRO.
!
! CALLING SEQUENCE:
!  call cgyro_init(...)
!  set cgyro_*_in variables
!  call cgyro_kernel(...)
!  get cgyro_*_out variables
!-------------------------------------------------------------------------

module qlgyro_cgyro_interface

  implicit none

  ! Input parameters (set to default values from cgyro/bin/cgyro_parse.py)
  integer :: cgyro_n_energy_in = 8
  integer :: cgyro_n_xi_in = 16
  integer :: cgyro_n_theta_in = 24
  integer :: cgyro_n_radial_in = 4
  integer :: cgyro_n_toroidal_in = 1
  integer :: cgyro_n_field_in = 1
  real :: cgyro_e_max_in = 8.0
  real :: cgyro_alpha_poly_in = 0.0
  integer :: cgyro_e_fix_in = 2
  integer :: cgyro_delta_t_method_in = 0
  real :: cgyro_delta_t_in = 0.01
  real :: cgyro_error_tol_in = 6e-5
  real :: cgyro_max_time_in = 100.0
  integer :: cgyro_print_step_in = 100
  integer :: cgyro_restart_step_in = 10
  real :: cgyro_freq_tol_in = 0.001
  integer :: cgyro_restart_flag_in = 0
  real :: cgyro_up_radial_in = 1.0
  real :: cgyro_up_theta_in = 1.0
  real :: cgyro_up_alpha_in = 0.0
  integer :: cgyro_nup_radial_in = 3
  integer :: cgyro_nup_theta_in = 3
  integer :: cgyro_nup_alpha_in = 3
  integer :: cgyro_n_wave_in = 2
  integer :: cgyro_constant_stream_flag_in = 1
  integer :: cgyro_explicit_trap_flag_in = 0
  real :: cgyro_ky_in = 0.3
  integer :: cgyro_box_size_in = 1
  real :: cgyro_ipccw_in = -1.0
  real :: cgyro_btccw_in = -1.0
  integer :: cgyro_silent_flag_in = 0
  integer :: cgyro_profile_model_in = 1
  integer :: cgyro_equilibrium_model_in = 2
  integer :: cgyro_collision_model_in = 4
  integer :: cgyro_collision_mom_restore_in = 1
  integer :: cgyro_collision_ene_restore_in = 1
  integer :: cgyro_collision_ene_diffusion_in = 1
  integer :: cgyro_collision_kperp_in = 1
  integer :: cgyro_collision_field_model_in = 1
  integer :: cgyro_collision_ion_model_in = 0
  integer :: cgyro_collision_precision_mode_in = 0
  real :: cgyro_z_eff_in = 1.0
  integer :: cgyro_z_eff_method_in = 2
  integer :: cgyro_zf_test_mode_in = 0
  integer :: cgyro_nonlinear_flag_in = 0
  integer :: cgyro_ae_flag_in = 0
  real :: cgyro_temp_ae_in = 1.0
  real :: cgyro_dens_ae_in = 1.0
  real :: cgyro_mass_ae_in = 0.0002724486
  real :: cgyro_dlntdr_ae_in = 1.0
  real :: cgyro_dlnndr_ae_in = 1.0
  real :: cgyro_lambda_star_in = 0.0
  integer :: cgyro_h_print_flag_in = 0
  integer :: cgyro_moment_print_flag_in = 0
  integer :: cgyro_gflux_print_flag_in = 0
  integer :: cgyro_field_print_flag_in = 0
  real :: cgyro_amp0_in = 0
  real :: cgyro_amp_in = 0.1
  real :: cgyro_gamma_e_in = 0.0
  real :: cgyro_gamma_p_in = 0.0
  real :: cgyro_mach_in = 0.0
  integer :: cgyro_rotation_model_in = 1
  integer :: cgyro_toroidals_per_proc_in = 1
  integer :: cgyro_mpi_rank_order_in = 2
  integer :: cgyro_velocity_order_in = 1
  integer :: cgyro_hiprec_flag_in = 0
  integer :: cgyro_udsymmetry_flag_in = 0
  integer :: cgyro_shear_method_in = 2
  integer :: cgyro_n_global_in = 4
  real :: cgyro_nu_global_in = 15.0
  integer :: cgyro_psym_flag_in = 0
  integer :: cgyro_profile_shear_flag_in = 0
  integer :: cgyro_theta_plot_in = 1
  integer :: cgyro_gpu_bigmem_flag_in = 1
  integer :: cgyro_upwind_single_flag_in = 0
  real :: cgyro_px0_in = 0.0
  integer :: cgyro_stream_term_in = 0
  real :: cgyro_stream_factor_in = 1.0
  integer :: cgyro_exch_flag_in = 0
  real :: cgyro_res_weight_power_in = 1.0

  real :: cgyro_rmin_in = 0.5
  real :: cgyro_rmaj_in = 3.0
  real :: cgyro_q_in = 2.0
  real :: cgyro_s_in = 1.0
  real :: cgyro_shift_in = 0.0
  real :: cgyro_kappa_in = 1.0
  real :: cgyro_s_kappa_in = 0.0
  real :: cgyro_delta_in = 0.0  
  real :: cgyro_s_delta_in = 0.0
  real :: cgyro_zeta_in = 0.0
  real :: cgyro_s_zeta_in = 0.0
  real :: cgyro_zmag_in = 0.0  
  real :: cgyro_dzmag_in = 0.0
  integer, parameter :: n_shape=6
  real :: cgyro_shape_sin3_in = 0.0
  real :: cgyro_shape_s_sin3_in = 0.0
  real :: cgyro_shape_sin4_in = 0.0
  real :: cgyro_shape_s_sin4_in = 0.0
  real :: cgyro_shape_sin5_in = 0.0
  real :: cgyro_shape_s_sin5_in = 0.0
  real :: cgyro_shape_sin6_in = 0.0
  real :: cgyro_shape_s_sin6_in = 0.0
  real :: cgyro_shape_cos0_in = 0.0
  real :: cgyro_shape_s_cos0_in = 0.0
  real :: cgyro_shape_cos1_in = 0.0
  real :: cgyro_shape_s_cos1_in = 0.0
  real :: cgyro_shape_cos2_in = 0.0
  real :: cgyro_shape_s_cos2_in = 0.0
  real :: cgyro_shape_cos3_in = 0.0
  real :: cgyro_shape_s_cos3_in = 0.0
  real :: cgyro_shape_cos4_in = 0.0
  real :: cgyro_shape_s_cos4_in = 0.0
  real :: cgyro_shape_cos5_in = 0.0
  real :: cgyro_shape_s_cos5_in = 0.0
  real :: cgyro_shape_cos6_in = 0.0
  real :: cgyro_shape_s_cos6_in = 0.0
  real :: cgyro_betae_unit_in = 0.0
  integer :: cgyro_n_species_in = 1

  real :: cgyro_nu_ee_in = 0.1
  
  real, dimension(11) :: cgyro_z_in = 1.0
  real, dimension(11) :: cgyro_mass_in = 0.0
  real, dimension(11) :: cgyro_dens_in = (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  real, dimension(11) :: cgyro_temp_in = 1.0
  real, dimension(11) :: cgyro_dlnndr_in = 1.0
  real, dimension(11) :: cgyro_dlntdr_in = 1.0
  real, dimension(11) :: cgyro_sdlnndr_in = 0.0
  real, dimension(11) :: cgyro_sdlntdr_in = 0.0

  integer :: cgyro_subroutine_flag_in = 1 
  integer :: cgyro_quasineutral_flag_in = 1 
  real :: cgyro_lambda_star_scale_in = 1.0
  real :: cgyro_gamma_e_scale_in = 1.0
  real :: cgyro_gamma_p_scale_in = 1.0
  real :: cgyro_mach_scale_in = 1.0
  real :: cgyro_beta_star_scale_in = 1.0
  real :: cgyro_betae_unit_scale_in = 1.0
  real :: cgyro_nu_ee_scale_in = 1.0
  
  real, dimension(11) :: cgyro_dlnndr_scale_in = 1.0
  real, dimension(11) :: cgyro_dlntdr_scale_in = 1.0

  ! Normalisation values
  real :: cgyro_b_unit_in=1
  real :: cgyro_b_gs2_in=1.0
  real :: cgyro_a_meters_in=1
  real :: cgyro_temp_norm_in=1.0
  real :: cgyro_dens_norm_in=1.0
  real :: cgyro_rho_star_norm_in=8.3333333e-4
  real :: cgyro_vth_norm_in=2.1850387e5
  
  ! Output parameters
  real, dimension(:), allocatable :: cgyro_elec_pflux_out 
  real, dimension(:), allocatable :: cgyro_elec_mflux_out 
  real, dimension(:), allocatable :: cgyro_elec_eflux_out
  real, dimension(:), allocatable :: cgyro_elec_expwd_out 

  real, dimension(:,:), allocatable :: cgyro_ion_pflux_out
  real, dimension(:,:), allocatable :: cgyro_ion_mflux_out
  real, dimension(:,:), allocatable :: cgyro_ion_eflux_out
  real, dimension(:,:), allocatable :: cgyro_ion_expwd_out


  complex, dimension(:, :, :), allocatable :: cgyro_wavefunction_out
  real, dimension(:), allocatable :: cgyro_thetab_out
  
  real, dimension(:), allocatable :: cgyro_r_out

  real, dimension(:, :, :, :), allocatable :: cgyro_gbflux_out
  real, dimension(:, :), allocatable :: cgyro_k_perp_out
  real, dimension(:), allocatable :: cgyro_jacobian_out

  complex, dimension(:), allocatable :: cgyro_omega_out
  complex, dimension(:), allocatable :: cgyro_omega_error_out

  integer :: cgyro_signal_out
  integer :: cgyro_error_status_out
  character(len=80) :: cgyro_error_message_out

  character(len=80) :: cgyro_path_in
  real, dimension(2) :: cgyro_integration_error_in

  logical :: cgyro_printout_in

contains


  ! Map GLOBAL variables to INTERFACE parameters
  subroutine map_global2interface()

    use cgyro_globals

    implicit none
    
    cgyro_n_energy_in = n_energy
    cgyro_n_xi_in = n_xi
    cgyro_n_theta_in = n_theta
    cgyro_n_radial_in = n_radial
    cgyro_n_toroidal_in = n_toroidal
    cgyro_n_field_in = n_field
    cgyro_e_max_in = e_max
    cgyro_alpha_poly_in = alpha_poly
    cgyro_e_fix_in = e_fix
    cgyro_delta_t_method_in = delta_t_method
    cgyro_delta_t_in = delta_t
    cgyro_error_tol_in = error_tol
    cgyro_max_time_in = max_time
    cgyro_print_step_in = print_step
    cgyro_restart_step_in = restart_step
    cgyro_freq_tol_in = freq_tol
    cgyro_restart_flag_in = restart_flag
    cgyro_up_radial_in = up_radial
    cgyro_up_theta_in = up_theta
    cgyro_up_alpha_in = up_alpha
    cgyro_nup_radial_in = nup_radial
    cgyro_nup_theta_in = nup_theta
    cgyro_nup_alpha_in = nup_alpha
    cgyro_n_wave_in = n_wave
    cgyro_constant_stream_flag_in = constant_stream_flag
    cgyro_explicit_trap_flag_in = explicit_trap_flag
    cgyro_ky_in = ky
    cgyro_box_size_in = box_size
    cgyro_ipccw_in = ipccw
    cgyro_btccw_in = btccw
    cgyro_silent_flag_in = silent_flag
    cgyro_profile_model_in = profile_model
    cgyro_equilibrium_model_in = equilibrium_model
    cgyro_collision_model_in = collision_model
    cgyro_collision_mom_restore_in = collision_mom_restore
    cgyro_collision_ene_restore_in = collision_ene_restore
    cgyro_collision_ene_diffusion_in = collision_ene_diffusion
    cgyro_collision_kperp_in = collision_kperp
    cgyro_collision_field_model_in = collision_field_model
    cgyro_collision_ion_model_in = collision_ion_model
    cgyro_collision_precision_mode_in = collision_precision_mode
    cgyro_z_eff_in = z_eff
    cgyro_z_eff_method_in = z_eff_method
    cgyro_zf_test_mode_in = zf_test_mode
    cgyro_nonlinear_flag_in = nonlinear_flag
    cgyro_ae_flag_in = ae_flag
    cgyro_temp_ae_in = temp_ae
    cgyro_dens_ae_in = dens_ae
    cgyro_mass_ae_in = mass_ae
    cgyro_dlntdr_ae_in = dlntdr_ae
    cgyro_dlnndr_ae_in = dlnndr_ae
    cgyro_lambda_star_in = lambda_star
    cgyro_h_print_flag_in = h_print_flag
    cgyro_moment_print_flag_in = moment_print_flag
    cgyro_gflux_print_flag_in = gflux_print_flag
    cgyro_field_print_flag_in = field_print_flag
    cgyro_amp0_in = amp0
    cgyro_amp_in = amp
    cgyro_gamma_e_in = gamma_e
    cgyro_gamma_p_in = gamma_p
    cgyro_mach_in = mach
    cgyro_rotation_model_in = rotation_model
    cgyro_toroidals_per_proc_in = nt_loc
    cgyro_mpi_rank_order_in = mpi_rank_order
    cgyro_velocity_order_in = velocity_order
    cgyro_hiprec_flag_in = hiprec_flag
    cgyro_udsymmetry_flag_in = udsymmetry_flag
    cgyro_shear_method_in = shear_method
    cgyro_n_global_in = n_global
    cgyro_nu_global_in = nu_global
    cgyro_psym_flag_in = psym_flag
    cgyro_profile_shear_flag_in = profile_shear_flag
    cgyro_theta_plot_in = theta_plot
    cgyro_gpu_bigmem_flag_in = gpu_bigmem_flag
    cgyro_upwind_single_flag_in = upwind_single_flag
    cgyro_px0_in = px0
    cgyro_stream_term_in = stream_term
    cgyro_stream_factor_in = stream_factor
    cgyro_exch_flag_in = exch_flag
    cgyro_res_weight_power_in = res_weight_power

    cgyro_rmin_in = rmin
    cgyro_rmaj_in = rmaj
    cgyro_q_in = q
    cgyro_s_in = s
    cgyro_shift_in = shift
    cgyro_kappa_in = kappa
    cgyro_s_kappa_in = s_kappa
    cgyro_delta_in = delta
    cgyro_s_delta_in = s_delta
    cgyro_zeta_in = zeta
    cgyro_s_zeta_in = s_zeta
    cgyro_zmag_in = zmag
    cgyro_dzmag_in = dzmag
    cgyro_shape_sin3_in = shape_sin(3)
    cgyro_shape_s_sin3_in = shape_s_sin(3)
    cgyro_shape_sin4_in = shape_sin(4)
    cgyro_shape_s_sin4_in = shape_s_sin(4)
    cgyro_shape_sin5_in = shape_sin(5)
    cgyro_shape_s_sin5_in = shape_s_sin(5)
    cgyro_shape_sin6_in = shape_sin(6)
    cgyro_shape_s_sin6_in = shape_s_sin(6)
    cgyro_shape_cos0_in = shape_cos(0)
    cgyro_shape_s_cos0_in = shape_s_cos(0)
    cgyro_shape_cos1_in = shape_cos(1)
    cgyro_shape_s_cos1_in = shape_s_cos(1)
    cgyro_shape_cos2_in = shape_cos(2)
    cgyro_shape_s_cos2_in = shape_s_cos(2)
    cgyro_shape_cos3_in = shape_cos(3)
    cgyro_shape_s_cos3_in = shape_s_cos(3)
    cgyro_shape_cos4_in = shape_cos(4)
    cgyro_shape_s_cos4_in = shape_s_cos(4)
    cgyro_shape_cos5_in = shape_cos(5)
    cgyro_shape_s_cos5_in = shape_s_cos(5)
    cgyro_shape_cos6_in = shape_cos(6)
    cgyro_shape_s_cos6_in = shape_s_cos(6)
    cgyro_betae_unit_in = betae_unit
    cgyro_n_species_in = n_species
    
    cgyro_nu_ee_in = nu_ee
    
    cgyro_z_in = z
    cgyro_mass_in = mass
    cgyro_dens_in = dens
    cgyro_temp_in = temp
    cgyro_dlnndr_in = dlnndr
    cgyro_dlntdr_in = dlntdr
    cgyro_sdlnndr_in = sdlnndr
    cgyro_sdlntdr_in = sdlntdr

    cgyro_subroutine_flag_in = subroutine_flag
    cgyro_quasineutral_flag_in = quasineutral_flag
    cgyro_lambda_star_scale_in = lambda_star
    cgyro_gamma_e_scale_in = gamma_e_scale
    cgyro_gamma_p_scale_in = gamma_p_scale
    cgyro_mach_scale_in = mach_scale
    cgyro_beta_star_scale_in = beta_star_scale
    cgyro_betae_unit_scale_in = betae_unit_scale
    cgyro_nu_ee_scale_in = nu_ee_scale
    
    cgyro_dlnndr_scale_in = dlnndr_scale
    cgyro_dlntdr_scale_in = dlntdr_scale

    cgyro_path_in = path
    cgyro_integration_error_in = integration_error

    ! Normalisation Values
    if (b_unit .ne. 0) cgyro_b_unit_in = b_unit
    if (b_gs2 .ne. 0) cgyro_b_gs2_in = b_gs2
    if (a_meters .ne. 0) cgyro_a_meters_in = a_meters
    if (temp_norm .ne. 0) cgyro_temp_norm_in = temp_norm
    if (dens_norm .ne. 0) cgyro_dens_norm_in = dens_norm
    if (rho_star_norm .ne. 0) cgyro_rho_star_norm_in = rho_star_norm
    if (vth_norm .ne. 0) cgyro_vth_norm_in = vth_norm
    
    ! Outputs
    cgyro_error_status_out = error_status
    cgyro_error_message_out = error_message

    cgyro_printout_in = printout

  end subroutine map_global2interface

  subroutine allocate_cgyro_interface

    use cgyro_globals

    implicit none
    
    ! Output parameters
        ! Allocate output arrays (deallocated in cgyro_cleanup)
    if (.not.allocated(cgyro_elec_pflux_out)) allocate(cgyro_elec_pflux_out(n_radial))
    if (.not.allocated(cgyro_elec_mflux_out)) allocate(cgyro_elec_mflux_out(n_radial))
    if (.not.allocated(cgyro_elec_eflux_out)) allocate(cgyro_elec_eflux_out(n_radial))
    if (.not.allocated(cgyro_elec_expwd_out)) allocate(cgyro_elec_expwd_out(n_radial))
    if (.not.allocated(cgyro_ion_pflux_out)) allocate(cgyro_ion_pflux_out(n_radial,5))
    if (.not.allocated(cgyro_ion_mflux_out)) allocate(cgyro_ion_mflux_out(n_radial,5))
    if (.not.allocated(cgyro_ion_eflux_out)) allocate(cgyro_ion_eflux_out(n_radial,5))
    if (.not.allocated(cgyro_ion_expwd_out)) allocate(cgyro_ion_expwd_out(n_radial,5))

    if(.not.allocated(cgyro_gbflux_out)) allocate(cgyro_gbflux_out(3, n_species, n_field, n_toroidal))

    if (.not.allocated(cgyro_wavefunction_out)) allocate(cgyro_wavefunction_out(n_toroidal, n_field, n_theta*n_radial))
    if (.not.allocated(cgyro_thetab_out)) allocate(cgyro_thetab_out(n_theta*n_radial))
    if (.not.allocated(cgyro_k_perp_out)) allocate(cgyro_k_perp_out(n_theta*n_radial, n_toroidal))
    if (.not.allocated(cgyro_jacobian_out)) allocate(cgyro_jacobian_out(n_theta*n_radial))

    if (.not.allocated(cgyro_omega_out)) allocate(cgyro_omega_out(n_toroidal))
    if (.not.allocated(cgyro_omega_error_out)) allocate(cgyro_omega_error_out(n_toroidal))

  end subroutine allocate_cgyro_interface

  ! Map INTERFACE parameters to GLOBAL variables
  subroutine map_interface2global()

    use cgyro_globals

    implicit none

    n_energy = cgyro_n_energy_in
    n_xi = cgyro_n_xi_in
    n_theta = cgyro_n_theta_in
    n_radial = cgyro_n_radial_in
    n_toroidal = cgyro_n_toroidal_in
    n_field = cgyro_n_field_in
    e_max = cgyro_e_max_in
    alpha_poly = cgyro_alpha_poly_in
    e_fix = cgyro_e_fix_in
    delta_t_method = cgyro_delta_t_method_in
    delta_t = cgyro_delta_t_in
    error_tol = cgyro_error_tol_in
    max_time = cgyro_max_time_in
    print_step = cgyro_print_step_in
    restart_step = cgyro_restart_step_in
    freq_tol = cgyro_freq_tol_in
    restart_flag = cgyro_restart_flag_in
    up_radial = cgyro_up_radial_in
    up_theta = cgyro_up_theta_in
    up_alpha = cgyro_up_alpha_in
    nup_radial = cgyro_nup_radial_in
    nup_theta = cgyro_nup_theta_in
    nup_alpha = cgyro_nup_alpha_in
    n_wave = cgyro_n_wave_in
    constant_stream_flag = cgyro_constant_stream_flag_in
    explicit_trap_flag = cgyro_explicit_trap_flag_in
    ky = cgyro_ky_in
    box_size = cgyro_box_size_in
    ipccw = cgyro_ipccw_in
    btccw = cgyro_btccw_in
    silent_flag = cgyro_silent_flag_in
    profile_model = cgyro_profile_model_in
    equilibrium_model = cgyro_equilibrium_model_in
    collision_model = cgyro_collision_model_in
    collision_mom_restore = cgyro_collision_mom_restore_in
    collision_ene_restore = cgyro_collision_ene_restore_in
    collision_ene_diffusion = cgyro_collision_ene_diffusion_in
    collision_kperp = cgyro_collision_kperp_in
    collision_field_model = cgyro_collision_field_model_in
    collision_ion_model = cgyro_collision_ion_model_in
    collision_precision_mode = cgyro_collision_precision_mode_in
    z_eff = cgyro_z_eff_in
    z_eff_method = cgyro_z_eff_method_in
    zf_test_mode = cgyro_zf_test_mode_in
    nonlinear_flag = cgyro_nonlinear_flag_in
    ae_flag = cgyro_ae_flag_in
    temp_ae = cgyro_temp_ae_in
    dens_ae = cgyro_dens_ae_in
    mass_ae = cgyro_mass_ae_in
    dlntdr_ae = cgyro_dlntdr_ae_in
    dlnndr_ae = cgyro_dlnndr_ae_in
    lambda_star = cgyro_lambda_star_in
    h_print_flag = cgyro_h_print_flag_in
    moment_print_flag = cgyro_moment_print_flag_in
    gflux_print_flag = cgyro_gflux_print_flag_in
    field_print_flag = cgyro_field_print_flag_in
    amp0 = cgyro_amp0_in
    amp = cgyro_amp_in
    gamma_e = cgyro_gamma_e_in
    gamma_p = cgyro_gamma_p_in
    mach = cgyro_mach_in
    rotation_model = cgyro_rotation_model_in
    nt_loc = cgyro_toroidals_per_proc_in
    mpi_rank_order = cgyro_mpi_rank_order_in
    velocity_order = cgyro_velocity_order_in
    hiprec_flag = cgyro_hiprec_flag_in
    udsymmetry_flag = cgyro_udsymmetry_flag_in
    shear_method = cgyro_shear_method_in
    n_global = cgyro_n_global_in
    nu_global = cgyro_nu_global_in
    psym_flag = cgyro_psym_flag_in
    profile_shear_flag = cgyro_profile_shear_flag_in
    theta_plot = cgyro_theta_plot_in
    gpu_bigmem_flag = cgyro_gpu_bigmem_flag_in
    px0 = cgyro_px0_in
    stream_term = cgyro_stream_term_in
    stream_factor = cgyro_stream_factor_in
    exch_flag = cgyro_exch_flag_in
    res_weight_power = cgyro_res_weight_power_in

    rmin = cgyro_rmin_in
    rmaj = cgyro_rmaj_in
    q = cgyro_q_in
    s = cgyro_s_in
    shift = cgyro_shift_in
    kappa = cgyro_kappa_in
    s_kappa = cgyro_s_kappa_in
    delta = cgyro_delta_in
    s_delta = cgyro_s_delta_in
    zeta = cgyro_zeta_in
    s_zeta = cgyro_s_zeta_in
    zmag = cgyro_zmag_in
    dzmag = cgyro_dzmag_in

    shape_sin(3) = cgyro_shape_sin3_in
    shape_s_sin(3) = cgyro_shape_s_sin3_in
    shape_sin(4) = cgyro_shape_sin4_in
    shape_s_sin(4) = cgyro_shape_s_sin4_in
    shape_sin(5) = cgyro_shape_sin5_in
    shape_s_sin(5) = cgyro_shape_s_sin5_in
    shape_sin(6) = cgyro_shape_sin6_in
    shape_s_sin(6) = cgyro_shape_s_sin6_in
    shape_cos(0) = cgyro_shape_cos0_in
    shape_s_cos(0) = cgyro_shape_s_cos0_in
    shape_cos(1) = cgyro_shape_cos1_in
    shape_s_cos(1) = cgyro_shape_s_cos1_in
    shape_cos(2) = cgyro_shape_cos2_in
    shape_s_cos(2) = cgyro_shape_s_cos2_in
    shape_cos(3) = cgyro_shape_cos3_in
    shape_s_cos(3) = cgyro_shape_s_cos3_in
    shape_cos(4) = cgyro_shape_cos4_in
    shape_s_cos(4) = cgyro_shape_s_cos4_in
    shape_cos(5) = cgyro_shape_cos5_in
    shape_s_cos(5) = cgyro_shape_s_cos5_in
    shape_cos(6) = cgyro_shape_cos6_in
    shape_s_cos(6) = cgyro_shape_s_cos6_in
    betae_unit = cgyro_betae_unit_in
    n_species = cgyro_n_species_in
    
    nu_ee = cgyro_nu_ee_in
    
    z = cgyro_z_in
    mass = cgyro_mass_in
    dens = cgyro_dens_in
    temp = cgyro_temp_in
    dlnndr = cgyro_dlnndr_in
    dlntdr = cgyro_dlntdr_in
    sdlnndr = cgyro_sdlnndr_in
    sdlntdr = cgyro_sdlntdr_in
    
    subroutine_flag = cgyro_subroutine_flag_in

    quasineutral_flag = cgyro_quasineutral_flag_in
    lambda_star = cgyro_lambda_star_scale_in
    gamma_e_scale = cgyro_gamma_e_scale_in
    gamma_p_scale = cgyro_gamma_p_scale_in
    mach_scale = cgyro_mach_scale_in
    beta_star_scale = cgyro_beta_star_scale_in
    betae_unit_scale = cgyro_betae_unit_scale_in
    nu_ee_scale = cgyro_nu_ee_scale_in
    
    dlnndr_scale = cgyro_dlnndr_scale_in
    dlntdr_scale = cgyro_dlntdr_scale_in

    ! Normalisation Values
    b_unit = cgyro_b_unit_in
    b_gs2 = cgyro_b_gs2_in
    a_meters = cgyro_a_meters_in
    temp_norm = cgyro_temp_norm_in
    dens_norm = cgyro_dens_norm_in
    rho_star_norm = cgyro_rho_star_norm_in
    vth_norm = cgyro_vth_norm_in

    path = cgyro_path_in
    integration_error = cgyro_integration_error_in
    printout = cgyro_printout_in

    signal = cgyro_signal_out
    error_status = cgyro_error_status_out

  end subroutine map_interface2global


  subroutine interfacelocaldump

    use cgyro_globals

    implicit none

    if (i_proc > 0) return

    open(unit=1, file=trim(path)//'out.cgyro.localdump', status='replace')
    write(1,20) 'N_ENERGY', cgyro_n_energy_in
    write(1,20) 'N_XI', cgyro_n_xi_in
    write(1,20) 'N_THETA', cgyro_n_theta_in
    write(1,20) 'N_RADIAL', cgyro_n_radial_in
    write(1,20) 'N_TOROIDAL', cgyro_n_toroidal_in
    write(1,20) 'N_FIELD', cgyro_n_field_in
    write(1,30) 'E_MAX', cgyro_e_max_in
    write(1,30) 'ALPHA_POLY', cgyro_alpha_poly_in
    write(1,20) 'E_FIX', cgyro_e_fix_in
    write(1,20) 'DELTA_T_METHOD', cgyro_delta_t_method_in
    write(1,30) 'DELTA_T', cgyro_delta_t_in
    write(1,30) 'ERROR_TOL',cgyro_error_tol_in
    write(1,30) 'MAX_TIME', cgyro_max_time_in
    write(1,20) 'PRINT_STEP',cgyro_print_step_in
    write(1,20) 'RESTART_STEP',cgyro_restart_step_in
    write(1,30) 'FREQ_TOL',cgyro_freq_tol_in
    write(1,30) 'UP_RADIAL',cgyro_up_radial_in
    write(1,30) 'UP_THETA',cgyro_up_theta_in
    write(1,30) 'UP_ALPHA',cgyro_up_alpha_in
    write(1,20) 'NUP_RADIAL',cgyro_nup_radial_in
    write(1,20) 'NUP_THETA',cgyro_nup_theta_in
    write(1,20) 'NUP_ALPHA',cgyro_nup_alpha_in
    write(1,20) 'N_WAVE',cgyro_n_wave_in
    write(1,20) 'CONSTANT_STREAM_FLAG',cgyro_constant_stream_flag_in
    write(1,20) 'EXPLICIT_TRAP_FLAG',cgyro_explicit_trap_flag_in
    write(1,30) 'KY',cgyro_ky_in
    write(1,20) 'BOX_SIZE',cgyro_box_size_in
    write(1,30) 'IPCCW',cgyro_ipccw_in
    write(1,30) 'BTCCW',cgyro_btccw_in
    write(1,20) 'SILENT_FLAG',cgyro_silent_flag_in
    write(1,20) 'PROFILE_MODEL',cgyro_profile_model_in
    write(1,20) 'EQUILIBRIUM_MODEL',cgyro_equilibrium_model_in
    write(1,20) 'COLLISION_MODEL',cgyro_collision_model_in
    write(1,20) 'COLLISION_MOM_RESTORE',cgyro_collision_mom_restore_in
    write(1,20) 'COLLISION_ENE_RESTORE',cgyro_collision_ene_restore_in
    write(1,20) 'COLLISION_ENE_DIFFUSION',cgyro_collision_ene_diffusion_in
    write(1,20) 'COLLISION_KPERP',cgyro_collision_kperp_in
    write(1,20) 'COLLISION_FIELD_MODEL',cgyro_collision_field_model_in
    write(1,20) 'COLLISION_ION_MODEL',cgyro_collision_ion_model_in
    write(1,20) 'COLLISION_PRECISION_MODE',cgyro_collision_precision_mode_in
    write(1,30) 'Z_EFF',cgyro_z_eff_in
    write(1,20) 'Z_EFF_METHOD',cgyro_z_eff_method_in
    write(1,20) 'ZF_TEST_MODE',cgyro_zf_test_mode_in
    write(1,20) 'NONLINEAR_FLAG',cgyro_nonlinear_flag_in
    write(1,20) 'AE_FLAG',cgyro_ae_flag_in
    write(1,30) 'TEMP_AE',cgyro_temp_ae_in
    write(1,30) 'DENS_AE',cgyro_dens_ae_in
    write(1,30) 'MASS_AE',cgyro_mass_ae_in
    write(1,30) 'DLNTDR_AE',cgyro_dlntdr_ae_in
    write(1,30) 'DLNNDR_AE',cgyro_dlnndr_ae_in
    write(1,30) 'LAMBDA_STAR',cgyro_lambda_star_in
    write(1,20) 'H_PRINT_FLAG',cgyro_h_print_flag_in
    write(1,20) 'MOMENT_PRINT_FLAG',cgyro_moment_print_flag_in
    write(1,20) 'GFLUX_PRINT_FLAG',cgyro_gflux_print_flag_in
    write(1,20) 'FIELD_PRINT_FLAG',cgyro_field_print_flag_in
    write(1,30) 'AMP0',cgyro_amp0_in
    write(1,30) 'AMP',cgyro_amp_in
    write(1,30) 'GAMMA_E',cgyro_gamma_e_in
    write(1,30) 'GAMMA_P',cgyro_gamma_p_in
    write(1,30) 'MACH',cgyro_mach_in
    write(1,20) 'ROTATION_MODEL',cgyro_rotation_model_in
    write(1,20) 'TOROIDALS_PER_PROC',cgyro_toroidals_per_proc_in
    write(1,20) 'MPI_RANK_ORDER',cgyro_mpi_rank_order_in
    write(1,20) 'VELOCITY_ORDER',cgyro_velocity_order_in
    write(1,20) 'HIPREC_FLAG',cgyro_hiprec_flag_in
    write(1,20) 'UDSYMMETRY_FLAG',cgyro_udsymmetry_flag_in
    write(1,20) 'SHEAR_METHOD',cgyro_shear_method_in
    write(1,20) 'N_GLOBAL',cgyro_n_global_in
    write(1,30) 'NU_GLOBAL',cgyro_nu_global_in
    write(1,20) 'PSYM_FLAG',cgyro_psym_flag_in
    write(1,20) 'PROFILE_SHEAR_FLAG',cgyro_profile_shear_flag_in
    write(1,20) 'THETA_PLOT', cgyro_theta_plot_in
    write(1,30) 'PX0',cgyro_px0_in
    write(1,20) 'STREAM_TERM',cgyro_stream_term_in
    write(1,30) 'STREAM_FACTOR',cgyro_stream_factor_in
    write(1,20) 'EXCH_FLAG',cgyro_exch_flag_in
    write(1,30) 'RES_WEIGHT_POWER',cgyro_res_weight_power_in

    write(1,30) 'RMIN',cgyro_rmin_in
    write(1,30) 'RMAJ',cgyro_rmaj_in
    write(1,30) 'Q',cgyro_q_in
    write(1,30) 'S',cgyro_s_in
    write(1,30) 'SHIFT',cgyro_shift_in
    write(1,30) 'KAPPA',cgyro_kappa_in
    write(1,30) 'S_KAPPA',cgyro_s_kappa_in
    write(1,30) 'DELTA',cgyro_delta_in
    write(1,30) 'S_DELTA',cgyro_s_delta_in
    write(1,30) 'ZETA',cgyro_zeta_in
    write(1,30) 'S_ZETA',cgyro_s_zeta_in
    write(1,30) 'ZMAG',cgyro_zmag_in
    write(1,30) 'DZMAG',cgyro_dzmag_in
    write(1,30) 'SHAPE_SIN3',cgyro_shape_sin3_in
    write(1,30) 'SHAPE_S_SIN3',cgyro_shape_s_sin3_in
    write(1,30) 'SHAPE_COS0',cgyro_shape_cos0_in
    write(1,30) 'SHAPE_S_COS0',cgyro_shape_s_cos0_in
    write(1,30) 'SHAPE_COS1',cgyro_shape_cos1_in
    write(1,30) 'SHAPE_S_COS1',cgyro_shape_s_cos1_in
    write(1,30) 'SHAPE_COS2',cgyro_shape_cos2_in
    write(1,30) 'SHAPE_S_COS2',cgyro_shape_s_cos2_in
    write(1,30) 'SHAPE_COS3',cgyro_shape_cos3_in
    write(1,30) 'SHAPE_S_COS3',cgyro_shape_s_cos3_in
    write(1,30) 'BETAE_UNIT',cgyro_betae_unit_in
    write(1,20) 'SUBROUTINE_FLAG',cgyro_subroutine_flag_in

    write(1,20) 'N_SPECIES',cgyro_n_species_in
    write(1,30) 'NU_EE',cgyro_nu_ee_in
    write(1,30) 'Z_1',cgyro_z_in(1)
    write(1,30) 'MASS_1',cgyro_mass_in(1)
    write(1,30) 'DENS_1',cgyro_dens_in(1)
    write(1,30) 'TEMP_1',cgyro_temp_in(1)
    write(1,30) 'DLNNDR_1',cgyro_dlnndr_in(1)
    write(1,30) 'DLNTDR_1',cgyro_dlntdr_in(1)
    write(1,30) 'SDLNNDR_1',cgyro_sdlnndr_in(1)
    write(1,30) 'SDLNTDR_1',cgyro_sdlntdr_in(1)
    write(1,30) 'Z_2',cgyro_z_in(2)
    write(1,30) 'MASS_2',cgyro_mass_in(2)
    write(1,30) 'DENS_2',cgyro_dens_in(2)
    write(1,30) 'TEMP_2',cgyro_temp_in(2)
    write(1,30) 'DLNNDR_2',cgyro_dlnndr_in(2)
    write(1,30) 'DLNTDR_2',cgyro_dlntdr_in(2)
    write(1,30) 'SDLNNDR_2',cgyro_sdlnndr_in(2)
    write(1,30) 'SDLNTDR_2',cgyro_sdlntdr_in(2)
    write(1,30) 'Z_3',cgyro_z_in(3)
    write(1,30) 'MASS_3',cgyro_mass_in(3)
    write(1,30) 'DENS_3',cgyro_dens_in(3)
    write(1,30) 'TEMP_3',cgyro_temp_in(3)
    write(1,30) 'DLNNDR_3',cgyro_dlnndr_in(3)
    write(1,30) 'DLNTDR_3',cgyro_dlntdr_in(3)
    write(1,30) 'SDLNNDR_3',cgyro_sdlnndr_in(3)
    write(1,30) 'SDLNTDR_3',cgyro_sdlntdr_in(3)
    write(1,30) 'Z_4',cgyro_z_in(4)
    write(1,30) 'MASS_4',cgyro_mass_in(4)
    write(1,30) 'DENS_4',cgyro_dens_in(4)
    write(1,30) 'TEMP_4',cgyro_temp_in(4)
    write(1,30) 'DLNNDR_4',cgyro_dlnndr_in(4)
    write(1,30) 'DLNTDR_4',cgyro_dlntdr_in(4)
    write(1,30) 'SDLNNDR_4',cgyro_sdlnndr_in(4)
    write(1,30) 'SDLNTDR_4',cgyro_sdlntdr_in(4)
    write(1,30) 'Z_5',cgyro_z_in(5)
    write(1,30) 'MASS_5',cgyro_mass_in(5)
    write(1,30) 'DENS_5',cgyro_dens_in(5)
    write(1,30) 'TEMP_5',cgyro_temp_in(5)
    write(1,30) 'DLNNDR_5',cgyro_dlnndr_in(5)
    write(1,30) 'DLNTDR_5',cgyro_dlntdr_in(5)
    write(1,30) 'SDLNNDR_5',cgyro_sdlnndr_in(5)
    write(1,30) 'SDLNTDR_5',cgyro_sdlntdr_in(5)
    write(1,30) 'Z_6',cgyro_z_in(6)
    write(1,30) 'MASS_6',cgyro_mass_in(6)
    write(1,30) 'DENS_6',cgyro_dens_in(6)
    write(1,30) 'TEMP_6',cgyro_temp_in(6)
    write(1,30) 'DLNNDR_6',cgyro_dlnndr_in(6)
    write(1,30) 'DLNTDR_6',cgyro_dlntdr_in(6)
    write(1,30) 'SDLNNDR_6',cgyro_sdlnndr_in(6)
    write(1,30) 'SDLNTDR_6',cgyro_sdlntdr_in(6)
    
    write(1,20) 'QUASINEUTRAL_FLAG',cgyro_quasineutral_flag_in
    write(1,30) 'LAMBDA_STAR_SCALE',cgyro_lambda_star_scale_in
    write(1,30) 'GAMMA_E_SCALE',cgyro_gamma_e_scale_in
    write(1,30) 'GAMMA_P_SCALE',cgyro_gamma_p_scale_in
    write(1,30) 'MACH_SCALE',cgyro_mach_scale_in
    write(1,30) 'BETA_STAR_SCALE',cgyro_beta_star_scale_in
    write(1,30) 'BETAE_UNIT_SCALE',cgyro_betae_unit_scale_in
    write(1,30) 'NU_EE_SCALE',cgyro_nu_ee_scale_in
    write(1,30) 'DLNNDR_1_SCALE',cgyro_dlnndr_in(1)
    write(1,30) 'DLNTDR_1_SCALE',cgyro_dlntdr_in(1)
    write(1,30) 'DLNNDR_2_SCALE',cgyro_dlnndr_in(2)
    write(1,30) 'DLNTDR_2_SCALE',cgyro_dlntdr_in(2)
    write(1,30) 'DLNNDR_3_SCALE',cgyro_dlnndr_in(3)
    write(1,30) 'DLNTDR_3_SCALE',cgyro_dlntdr_in(3)
    write(1,30) 'DLNNDR_4_SCALE',cgyro_dlnndr_in(4)
    write(1,30) 'DLNTDR_4_SCALE',cgyro_dlntdr_in(4)
    write(1,30) 'DLNNDR_5_SCALE',cgyro_dlnndr_in(5)
    write(1,30) 'DLNTDR_5_SCALE',cgyro_dlntdr_in(5)
    write(1,30) 'DLNNDR_6_SCALE',cgyro_dlnndr_in(6)
    write(1,30) 'DLNTDR_6_SCALE',cgyro_dlntdr_in(6)
    close(1)
    
20  format(a,'=',i3)
30  format(a,'=',1pe12.5)
    
  end subroutine interfacelocaldump
    
end module qlgyro_cgyro_interface
