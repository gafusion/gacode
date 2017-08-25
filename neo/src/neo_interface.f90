!-------------------------------------------------------------------------
! neo_interface.f90
!
! PURPOSE:
!  Provides interface description for NEO.
!
! CALLING SEQUENCE:
!  call neo_init(...)
!  set neo_*_in variables
!  call neo_run(...)
!  get neo_*_out variables
!-------------------------------------------------------------------------

module neo_interface

  implicit none

  ! Input parameters (set to default values from python/neo_dict.py)
  integer :: neo_n_energy_in = 6
  integer :: neo_n_xi_in = 17
  integer :: neo_n_theta_in = 17
  integer :: neo_n_radial_in = 1
  integer :: neo_matsz_scalefac_in = 50
  real    :: neo_rmin_over_a_in = 0.5
  real    :: neo_rmin_over_a_2_in = 0.6
  real    :: neo_rmaj_over_a_in = 3.0
  integer :: neo_silent_flag_in = 0
  integer :: neo_sim_model_in = 2
  integer :: neo_equilibrium_model_in = 0
  integer :: neo_collision_model_in = 4
  integer :: neo_profile_model_in = 1
  integer :: neo_profile_erad0_model_in = 1
  integer :: neo_profile_equilibrium_model_in = 1
  integer :: neo_ipccw_in = -1
  integer :: neo_btccw_in = -1
  real    :: neo_te_ade_in = 1.0
  real    :: neo_ne_ade_in = 1.0
  real    :: neo_dlntdre_ade_in  = 1.0
  real    :: neo_dlnndre_ade_in  = 1.0
  integer :: neo_rotation_model_in = 1
  real    :: neo_omega_rot_in = 0.0
  real    :: neo_omega_rot_deriv_in = 0.0
  integer :: neo_spitzer_model_in = 0
  real    :: neo_epar0_spitzer_in = 1.0
  integer :: neo_coll_uncoupledei_model_in = 0
  integer :: neo_coll_uncoupledaniso_model_in = 0
  integer :: neo_n_species_in = 1
  real                   :: neo_nu_1_in = 0.1
  real,    dimension(11) :: neo_z_in = (/1,1,1,1,1,1,1,1,1,1,1/)
  real,    dimension(11) :: neo_mass_in = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  real,    dimension(11) :: neo_dens_in = (/1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/)
  real,    dimension(11) :: neo_temp_in = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  real,    dimension(11) :: neo_dlnndr_in = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  real,    dimension(11) :: neo_dlntdr_in = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  integer, dimension(11) :: neo_aniso_model_in  = (/1,1,1,1,1,1,1,1,1,1,1/)
  real,    dimension(11) :: neo_temp_para_in    = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  real,    dimension(11) :: neo_dlntdr_para_in  = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  real,    dimension(11) :: neo_temp_perp_in    = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  real,    dimension(11) :: neo_dlntdr_perp_in  = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  real,    dimension(11) :: neo_profile_dlnndr_scale_in = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  real,    dimension(11) :: neo_profile_dlntdr_scale_in = (/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
  real    :: neo_dphi0dr_in = 0.0
  real    :: neo_epar0_in = 0.0
  real    :: neo_q_in = 2.0
  real    :: neo_rho_star_in = 0.001
  real    :: neo_shear_in = 1.0
  real    :: neo_shift_in = 0.0
  real    :: neo_kappa_in = 1.0
  real    :: neo_s_kappa_in = 0.0
  real    :: neo_delta_in = 0.0
  real    :: neo_s_delta_in = 0.0
  real    :: neo_zeta_in = 0.0
  real    :: neo_s_zeta_in = 0.0
  real    :: neo_zmag_over_a_in = 0.0
  real    :: neo_s_zmag_in = 0.0
  real    :: neo_beta_star_in = 0.0
  real    :: neo_profile_delta_scale_in = 1.0
  real    :: neo_profile_zeta_scale_in = 1.0
  real    :: neo_profile_zmag_scale_in = 1.0
  integer :: neo_threed_model_in = 0
  integer :: neo_threed_exb_model_in = 0
  real    :: neo_threed_exb_dphi0dr_in = 0.0
  integer :: neo_threed_drift_model_in = 0
  real    :: neo_threed_hyperxi_in = 0.0
  integer :: neo_laguerre_method_in = 1
  integer :: neo_write_cmoments_flag_in = 0
  integer :: neo_geo_ny_in = 0
  real, dimension(8,0:32) :: neo_geo_yin_in = 0.0
  ! the exception of the default is subroutine_flag
  integer :: neo_subroutine_flag = 1
  integer :: neo_test_flag_in = 0
  
  ! Output parameters
  ! theory
  real    :: neo_pflux_thHH_out  = 0.0
  real    :: neo_eflux_thHHi_out = 0.0
  real    :: neo_eflux_thHHe_out = 0.0
  real    :: neo_eflux_thCHi_out = 0.0
  real, dimension(11) :: neo_pflux_thHS_out = 0.0
  real, dimension(11) :: neo_eflux_thHS_out = 0.0
  real    :: neo_jpar_thS_out    = 0.0
  real    :: neo_jpar_thK_out    = 0.0
  real    :: neo_jpar_thN_out    = 0.0
  real    :: neo_jtor_thS_out    = 0.0
  ! drift-kinetic soln
  real, dimension(11) :: neo_pflux_dke_out    = 0.0
  real, dimension(11) :: neo_efluxtot_dke_out = 0.0
  real, dimension(11) :: neo_efluxncv_dke_out = 0.0
  real, dimension(11) :: neo_mflux_dke_out    = 0.0
  real, dimension(11) :: neo_vpol_dke_out     = 0.0
  real, dimension(11) :: neo_vtor_dke_out     = 0.0
  real                :: neo_jpar_dke_out     = 0.0
  real                :: neo_jtor_dke_out     = 0.0
  ! gyro-viscosity
  real, dimension(11) :: neo_pflux_gv_out  = 0.0
  real, dimension(11) :: neo_efluxtot_gv_out  = 0.0
  real, dimension(11) :: neo_efluxncv_gv_out  = 0.0
  real, dimension(11) :: neo_mflux_gv_out  = 0.0
  ! nclass
  real, dimension(11) :: neo_nclassvis_out       = 0.0
  real, dimension(11) :: neo_pflux_nclass_out    = 0.0
  real, dimension(11) :: neo_efluxtot_nclass_out = 0.0
  real, dimension(11) :: neo_vpol_nclass_out     = 0.0
  real, dimension(11) :: neo_vtor_nclass_out     = 0.0
  real                :: neo_jpar_nclass_out     = 0.0
  ! error checking
  integer :: neo_error_status_out=0
  character(len=80) :: neo_error_message_out=''

contains

  ! Map GLOBAL variables to INTERFACE parameters
  subroutine map_global2interface()

    use neo_globals

    implicit none

    neo_n_energy_in = n_energy
    neo_n_xi_in = n_xi
    neo_n_theta_in = n_theta
    neo_n_radial_in = n_radial
    neo_matsz_scalefac_in = matsz_scalefac
    neo_rmin_over_a_in = rmin_1_in
    neo_rmin_over_a_2_in = rmin_2_in
    neo_rmaj_over_a_in = rmaj_in
    neo_silent_flag_in = silent_flag
    neo_sim_model_in = sim_model
    neo_equilibrium_model_in = equilibrium_model
    neo_collision_model_in = collision_model
    neo_profile_model_in = profile_model
    neo_profile_erad0_model_in = profile_erad0_model
    neo_profile_equilibrium_model_in = profile_equilibrium_model
    neo_ipccw_in  = ipccw_in
    neo_btccw_in  = btccw_in
    neo_te_ade_in = te_ade_in
    neo_ne_ade_in = ne_ade_in
    neo_dlntdre_ade_in = dlntdre_ade_in
    neo_dlnndre_ade_in = dlnndre_ade_in
    neo_rotation_model_in = rotation_model
    neo_omega_rot_in = omega_rot_in
    neo_omega_rot_deriv_in = omega_rot_deriv_in
    neo_spitzer_model_in = spitzer_model
    neo_epar0_spitzer_in = epar0_spitzer
    neo_coll_uncoupledei_model_in = coll_uncoupledei_model
    neo_coll_uncoupledaniso_model_in = coll_uncoupledaniso_model
    neo_n_species_in = n_species
    neo_nu_1_in = nu_1_in  
    neo_z_in = z_in
    neo_mass_in = mass_in
    neo_dens_in = dens_in
    neo_temp_in = temp_in
    neo_dlnndr_in = dlnndr_in
    neo_dlntdr_in = dlntdr_in
    neo_aniso_model_in = aniso_model_in
    neo_temp_para_in   = temp_para_in
    neo_dlntdr_para_in = dlntdr_para_in
    neo_temp_perp_in   = temp_perp_in
    neo_dlntdr_perp_in = dlntdr_perp_in
    neo_profile_dlnndr_scale_in = profile_dlnndr_scale
    neo_profile_dlntdr_scale_in = profile_dlntdr_scale
    neo_dphi0dr_in = dphi0dr_in
    neo_epar0_in = epar0_in
    neo_q_in = q_in
    neo_rho_star_in = rho_in
    neo_shear_in = shear_in
    neo_shift_in = shift_in
    neo_kappa_in = kappa_in
    neo_s_kappa_in = s_kappa_in
    neo_delta_in = delta_in
    neo_s_delta_in = s_delta_in
    neo_zeta_in = zeta_in
    neo_s_zeta_in = s_zeta_in
    neo_zmag_over_a_in = zmag_in
    neo_s_zmag_in = s_zmag_in
    neo_beta_star_in = beta_star_in
    neo_profile_delta_scale_in = profile_delta_scale
    neo_profile_zeta_scale_in = profile_zeta_scale
    neo_profile_zmag_scale_in = profile_zmag_scale
    neo_threed_model_in = threed_model
    neo_threed_exb_model_in = threed_exb_model
    neo_threed_exb_dphi0dr_in = threed_exb_dphi0dr
    neo_threed_drift_model_in = threed_drift_model
    neo_threed_hyperxi_in = threed_hyperxi
    neo_laguerre_method_in = laguerre_method 
    neo_write_cmoments_flag_in = write_cmoments_flag 
    neo_geo_ny_in = geo_ny_in
    neo_geo_yin_in(:,:) = geo_yin_in(:,:)
    neo_subroutine_flag = subroutine_flag

  end subroutine map_global2interface

  ! Map INTERFACE parameters to GLOBAL variables
  subroutine map_interface2global()

    use neo_globals

    implicit none

    n_energy = neo_n_energy_in
    n_xi = neo_n_xi_in
    n_theta = neo_n_theta_in
    n_radial = neo_n_radial_in
    matsz_scalefac = neo_matsz_scalefac_in
    rmin_1_in = neo_rmin_over_a_in
    rmin_2_in = neo_rmin_over_a_2_in
    rmaj_in = neo_rmaj_over_a_in
    silent_flag = neo_silent_flag_in
    sim_model = neo_sim_model_in
    equilibrium_model = neo_equilibrium_model_in
    collision_model = neo_collision_model_in
    profile_model = neo_profile_model_in
    profile_erad0_model = neo_profile_erad0_model_in
    profile_equilibrium_model = neo_profile_equilibrium_model_in
    ipccw_in  = neo_ipccw_in
    btccw_in  = neo_btccw_in
    te_ade_in = neo_te_ade_in
    ne_ade_in = neo_ne_ade_in
    dlntdre_ade_in = neo_dlntdre_ade_in
    dlnndre_ade_in = neo_dlnndre_ade_in   
    rotation_model = neo_rotation_model_in
    omega_rot_in = neo_omega_rot_in
    omega_rot_deriv_in = neo_omega_rot_deriv_in
    spitzer_model = neo_spitzer_model_in
    epar0_spitzer = neo_epar0_spitzer_in
    coll_uncoupledei_model = neo_coll_uncoupledei_model_in
    coll_uncoupledaniso_model = neo_coll_uncoupledaniso_model_in
    n_species = neo_n_species_in
    nu_1_in = neo_nu_1_in
    z_in = neo_z_in
    mass_in = neo_mass_in
    dens_in = neo_dens_in
    temp_in = neo_temp_in
    dlnndr_in = neo_dlnndr_in
    dlntdr_in = neo_dlntdr_in
    aniso_model_in = neo_aniso_model_in 
    temp_para_in   = neo_temp_para_in   
    dlntdr_para_in = neo_dlntdr_para_in 
    temp_perp_in   = neo_temp_perp_in   
    dlntdr_perp_in = neo_dlntdr_perp_in 
    profile_dlnndr_scale = neo_profile_dlnndr_scale_in
    profile_dlntdr_scale = neo_profile_dlntdr_scale_in
    dphi0dr_in = neo_dphi0dr_in
    epar0_in = neo_epar0_in
    q_in = neo_q_in
    rho_in = neo_rho_star_in
    shear_in = neo_shear_in
    shift_in = neo_shift_in
    kappa_in = neo_kappa_in
    s_kappa_in = neo_s_kappa_in
    delta_in = neo_delta_in
    s_delta_in = neo_s_delta_in
    zeta_in = neo_zeta_in
    s_zeta_in = neo_s_zeta_in
    zmag_in = neo_zmag_over_a_in
    s_zmag_in = neo_s_zmag_in
    beta_star_in = neo_beta_star_in
    profile_delta_scale = neo_profile_delta_scale_in
    profile_zeta_scale = neo_profile_zeta_scale_in
    profile_zmag_scale = neo_profile_zmag_scale_in
    threed_model = neo_threed_model_in
    threed_exb_model = neo_threed_exb_model_in
    threed_exb_dphi0dr = neo_threed_exb_dphi0dr_in
    threed_drift_model = neo_threed_drift_model_in
    threed_hyperxi = neo_threed_hyperxi_in
    laguerre_method = neo_laguerre_method_in
    write_cmoments_flag = neo_write_cmoments_flag_in 
    geo_ny_in = neo_geo_ny_in
    geo_yin_in(:,:) = neo_geo_yin_in(:,:)
    subroutine_flag = neo_subroutine_flag

    call interfacelocaldump

  end subroutine map_interface2global

  subroutine interfacelocaldump
    use neo_globals
    implicit none

    if(neo_silent_flag_in > 0 .or. i_proc > 0) return

    open(unit=1,file=trim(path)//'out.neo.localdump',status='replace')

    if(neo_n_radial_in > 1) then
       write(1,*) 'Localdump only works with n_radial=1'
       close(1)
       return
    endif
    if(neo_profile_model_in /= 1) then
       write(1,*) 'Localdump only works with profile_model=1'
       close(1)
       return
    endif

    write(1,20) 'N_ENERGY=',neo_n_energy_in
    write(1,20) 'N_XI=',neo_n_xi_in
    write(1,20) 'N_THETA=',neo_n_theta_in
    write(1,20) 'N_RADIAL=',neo_n_radial_in
    write(1,20) 'MATSZ_SCALEFAC=',neo_matsz_scalefac_in
    write(1,30) 'RMIN_OVER_A=',neo_rmin_over_a_in
    write(1,30) 'RMAJ_OVER_A=',neo_rmaj_over_a_in
    write(1,20) 'SILENT_FLAG=',neo_silent_flag_in
    write(1,20) 'SIM_MODEL=',neo_sim_model_in
    write(1,20) 'EQUILIBRIUM_MODEL=',neo_equilibrium_model_in
    write(1,20) 'COLLISION_MODEL=',neo_collision_model_in
    write(1,20) 'IPCCW=',neo_ipccw_in
    write(1,20) 'BTCCW=',neo_btccw_in
    write(1,30) 'TE_ADE=',neo_te_ade_in
    write(1,30) 'NE_ADE=',neo_ne_ade_in
    write(1,30) 'DLNTDRE_ADE=',neo_dlntdre_ade_in
    write(1,30) 'DLNNDRE_ADE=',neo_dlnndre_ade_in
    write(1,20) 'ROTATION_MODEL=',neo_rotation_model_in
    write(1,30) 'OMEGA_ROT=',neo_omega_rot_in
    write(1,30) 'OMEGA_ROT_DERIV=',neo_omega_rot_deriv_in
    write(1,20) 'SPITZER_MODEL=',neo_spitzer_model_in
    write(1,30) 'EPAR0_SPITZER=',neo_epar0_spitzer_in
    write(1,20) 'COLL_UNCOUPLEDEI_MODEL=',neo_coll_uncoupledei_model_in
    write(1,20) 'COLL_UNCOUPLEDANISO_MODEL=',neo_coll_uncoupledaniso_model_in
    write(1,20) 'N_SPECIES=',neo_n_species_in
    write(1,30) 'NU_1=',neo_nu_1_in

    write(1,30) 'Z_1=',neo_z_in(1)
    write(1,30) 'MASS_1=',neo_mass_in(1)
    write(1,30) 'DENS_1=',neo_dens_in(1)
    write(1,30) 'TEMP_1=',neo_temp_in(1)
    write(1,30) 'DLNNDR_1=',neo_dlnndr_in(1)
    write(1,30) 'DLNTDR_1=',neo_dlntdr_in(1)
    write(1,20) 'ANISO_MODEL_1=', neo_aniso_model_in(1)  
    if(neo_aniso_model_in(1) == 2) then
       write(1,30) 'TEMP_PARA_1=', neo_temp_para_in(1)  
       write(1,30) 'DLNTDR_PARA_1=', neo_dlntdr_para_in(1)  
       write(1,30) 'TEMP_PERP_1=', neo_temp_perp_in(1)  
       write(1,30) 'DLNTDR_PERP_1=', neo_dlntdr_perp_in(1)  
    endif
    if (neo_n_species_in > 1) then
       write(1,30) 'Z_2=',neo_z_in(2)
       write(1,30) 'MASS_2=',neo_mass_in(2)
       write(1,30) 'DENS_2=',neo_dens_in(2)
       write(1,30) 'TEMP_2=',neo_temp_in(2)
       write(1,30) 'DLNNDR_2=',neo_dlnndr_in(2)
       write(1,30) 'DLNTDR_2=',neo_dlntdr_in(2)
       write(1,20) 'ANISO_MODEL_2=', neo_aniso_model_in(2)  
       if(neo_aniso_model_in(2) == 2) then 
          write(1,30) 'TEMP_PARA_2=', neo_temp_para_in(2) 
          write(1,30) 'DLNTDR_PARA_2=', neo_dlntdr_para_in(2)
          write(1,30) 'TEMP_PERP_2=', neo_temp_perp_in(2) 
          write(1,30) 'DLNTDR_PERP_2=', neo_dlntdr_perp_in(2) 
       endif
    endif
    if (neo_n_species_in > 2) then
       write(1,30) 'Z_3=',neo_z_in(3)
       write(1,30) 'MASS_3=',neo_mass_in(3)
       write(1,30) 'DENS_3=',neo_dens_in(3)
       write(1,30) 'TEMP_3=',neo_temp_in(3)
       write(1,30) 'DLNNDR_3=',neo_dlnndr_in(3)
       write(1,30) 'DLNTDR_3=',neo_dlntdr_in(3)
       write(1,20) 'ANISO_MODEL_3=', neo_aniso_model_in(3)  
       if(neo_aniso_model_in(3) == 2) then              
          write(1,30) 'TEMP_PARA_3=', neo_temp_para_in(3) 
          write(1,30) 'DLNTDR_PARA_3=', neo_dlntdr_para_in(3)
          write(1,30) 'TEMP_PERP_3=', neo_temp_perp_in(3) 
          write(1,30) 'DLNTDR_PERP_3=', neo_dlntdr_perp_in(3) 
       endif
    endif
    if (neo_n_species_in > 3) then
       write(1,30) 'Z_4=',neo_z_in(4)
       write(1,30) 'MASS_4=',neo_mass_in(4)
       write(1,30) 'DENS_4=',neo_dens_in(4)
       write(1,30) 'TEMP_4=',neo_temp_in(4)
       write(1,30) 'DLNNDR_4=',neo_dlnndr_in(4)
       write(1,30) 'DLNTDR_4=',neo_dlntdr_in(4)
       write(1,20) 'ANISO_MODEL_4=', neo_aniso_model_in(4)  
       if(neo_aniso_model_in(4) == 2) then              
          write(1,30) 'TEMP_PARA_4=', neo_temp_para_in(4) 
          write(1,30) 'DLNTDR_PARA_4=', neo_dlntdr_para_in(4)
          write(1,30) 'TEMP_PERP_4=', neo_temp_perp_in(4) 
          write(1,30) 'DLNTDR_PERP_4=', neo_dlntdr_perp_in(4) 
       endif
    endif
    if (neo_n_species_in > 4) then
       write(1,30) 'Z_5=',neo_z_in(5)
       write(1,30) 'MASS_5=',neo_mass_in(5)
       write(1,30) 'DENS_5=',neo_dens_in(5)
       write(1,30) 'TEMP_5=',neo_temp_in(5)
       write(1,30) 'DLNNDR_5=',neo_dlnndr_in(5)
       write(1,30) 'DLNTDR_5=',neo_dlntdr_in(5)
       write(1,20) 'ANISO_MODEL_5=', neo_aniso_model_in(5)  
       if(neo_aniso_model_in(5) == 2) then              
          write(1,30) 'TEMP_PARA_5=', neo_temp_para_in(5) 
          write(1,30) 'DLNTDR_PARA_5=', neo_dlntdr_para_in(5)
          write(1,30) 'TEMP_PERP_5=', neo_temp_perp_in(5) 
          write(1,30) 'DLNTDR_PERP_5=', neo_dlntdr_perp_in(5) 
       endif
    endif
    if (neo_n_species_in > 5) then
       write(1,30) 'Z_6=',neo_z_in(6)
       write(1,30) 'MASS_6=',neo_mass_in(6)
       write(1,30) 'DENS_6=',neo_dens_in(6)
       write(1,30) 'TEMP_6=',neo_temp_in(6)
       write(1,30) 'DLNNDR_6=',neo_dlnndr_in(6)
       write(1,30) 'DLNTDR_6=',neo_dlntdr_in(6)
       write(1,20) 'ANISO_MODEL_6=', neo_aniso_model_in(6)  
       if(neo_aniso_model_in(6) == 2) then              
          write(1,30) 'TEMP_PARA_6=', neo_temp_para_in(6) 
          write(1,30) 'DLNTDR_PARA_6=', neo_dlntdr_para_in(6)
          write(1,30) 'TEMP_PERP_6=', neo_temp_perp_in(6) 
          write(1,30) 'DLNTDR_PERP_6=', neo_dlntdr_perp_in(6) 
       endif
    endif
    if (neo_n_species_in > 6) then
       write(1,30) 'Z_7=',neo_z_in(7)
       write(1,30) 'MASS_7=',neo_mass_in(7)
       write(1,30) 'DENS_7=',neo_dens_in(7)
       write(1,30) 'TEMP_7=',neo_temp_in(7)
       write(1,30) 'DLNNDR_7=',neo_dlnndr_in(7)
       write(1,30) 'DLNTDR_7=',neo_dlntdr_in(7)
       write(1,20) 'ANISO_MODEL_7=', neo_aniso_model_in(7)  
       if(neo_aniso_model_in(7) == 2) then              
          write(1,30) 'TEMP_PARA_7=', neo_temp_para_in(7) 
          write(1,30) 'DLNTDR_PARA_7=', neo_dlntdr_para_in(7)
          write(1,30) 'TEMP_PERP_7=', neo_temp_perp_in(7)
          write(1,30) 'DLNTDR_PERP_7=', neo_dlntdr_perp_in(7) 
       endif
    endif
    if (neo_n_species_in > 7) then
       write(1,30) 'Z_8=',neo_z_in(8)
       write(1,30) 'MASS_8=',neo_mass_in(8)
       write(1,30) 'DENS_8=',neo_dens_in(8)
       write(1,30) 'TEMP_8=',neo_temp_in(8)
       write(1,30) 'DLNNDR_8=',neo_dlnndr_in(8)
       write(1,30) 'DLNTDR_8=',neo_dlntdr_in(8)
       write(1,20) 'ANISO_MODEL_8=', neo_aniso_model_in(8)  
       if(neo_aniso_model_in(8) == 2) then              
          write(1,30) 'TEMP_PARA_8=', neo_temp_para_in(8) 
          write(1,30) 'DLNTDR_PARA_8=', neo_dlntdr_para_in(8)
          write(1,30) 'TEMP_PERP_8=', neo_temp_perp_in(8) 
          write(1,30) 'DLNTDR_PERP_8=', neo_dlntdr_perp_in(8) 
       endif
    endif
    if (neo_n_species_in > 8) then
       write(1,30) 'Z_9=',neo_z_in(9)
       write(1,30) 'MASS_9=',neo_mass_in(9)
       write(1,30) 'DENS_9=',neo_dens_in(9)
       write(1,30) 'TEMP_9=',neo_temp_in(9)
       write(1,30) 'DLNNDR_9=',neo_dlnndr_in(9)
       write(1,30) 'DLNTDR_9=',neo_dlntdr_in(9)
       write(1,20) 'ANISO_MODEL_9=', neo_aniso_model_in(9)  
       if(neo_aniso_model_in(9) == 2) then              
          write(1,30) 'TEMP_PARA_9=', neo_temp_para_in(9) 
          write(1,30) 'DLNTDR_PARA_9=', neo_dlntdr_para_in(9)
          write(1,30) 'TEMP_PERP_9=', neo_temp_perp_in(9) 
          write(1,30) 'DLNTDR_PERP_9=', neo_dlntdr_perp_in(9) 
       endif
    endif
    if (neo_n_species_in > 9) then
       write(1,30) 'Z_10=',neo_z_in(10)
       write(1,30) 'MASS_10=',neo_mass_in(10)
       write(1,30) 'DENS_10=',neo_dens_in(10)
       write(1,30) 'TEMP_10=',neo_temp_in(10)
       write(1,30) 'DLNNDR_10=',neo_dlnndr_in(10)
       write(1,30) 'DLNTDR_10=',neo_dlntdr_in(10)
       write(1,20) 'ANISO_MODEL_10=', neo_aniso_model_in(10)  
       if(neo_aniso_model_in(10) == 2) then              
          write(1,30) 'TEMP_PARA_10=', neo_temp_para_in(10) 
          write(1,30) 'DLNTDR_PARA_10=', neo_dlntdr_para_in(10)
          write(1,30) 'TEMP_PERP_10=', neo_temp_perp_in(10) 
          write(1,30) 'DLNTDR_PERP_10=', neo_dlntdr_perp_in(10) 
       endif
    endif
    if (neo_n_species_in > 10) then
       write(1,30) 'Z_11=',neo_z_in(11)
       write(1,30) 'MASS_11=',neo_mass_in(11)
       write(1,30) 'DENS_11=',neo_dens_in(11)
       write(1,30) 'TEMP_11=',neo_temp_in(11)
       write(1,30) 'DLNNDR_11=',neo_dlnndr_in(11)
       write(1,30) 'DLNTDR_11=',neo_dlntdr_in(11)
       write(1,20) 'ANISO_MODEL_11=', neo_aniso_model_in(11)  
       if(neo_aniso_model_in(11) == 2) then              
          write(1,30) 'TEMP_PARA_11=', neo_temp_para_in(11) 
          write(1,30) 'DLNTDR_PARA_11=', neo_dlntdr_para_in(11)
          write(1,30) 'TEMP_PERP_11=', neo_temp_perp_in(11) 
          write(1,30) 'DLNTDR_PERP_11=', neo_dlntdr_perp_in(11) 
       endif
    endif

    write(1,30) 'DPHI0DR=',neo_dphi0dr_in
    write(1,30) 'EPAR0=',neo_epar0_in
    write(1,30) 'Q=',neo_q_in
    write(1,30) 'RHO_STAR=',neo_rho_star_in
    write(1,30) 'SHEAR=',neo_shear_in
    write(1,30) 'SHIFT=',neo_shift_in
    write(1,30) 'KAPPA=',neo_kappa_in
    write(1,30) 'S_KAPPA=',neo_s_kappa_in
    write(1,30) 'DELTA=',neo_delta_in
    write(1,30) 'S_DELTA=',neo_s_delta_in
    write(1,30) 'ZETA=',neo_zeta_in
    write(1,30) 'S_ZETA=',neo_s_zeta_in
    write(1,30) 'ZMAG_OVER_A=',neo_zmag_over_a_in
    write(1,30) 'S_ZMAG=',neo_s_zmag_in
    write(1,30) 'BETA_STAR=',neo_beta_star_in

    write(1,20) 'LAGUERRE_METHOD=',neo_laguerre_method_in
    write(1,20) 'WRITE_CMOMENTS_FLAG=',neo_write_cmoments_flag_in

    close(1)

    if(neo_equilibrium_model_in == 3) then
       open(unit=1,file=trim(path)//'out.neo.localdump_geo',status='replace')
       write(1,10) neo_geo_ny_in
       write(1,40) neo_geo_yin_in(:,0:neo_geo_ny_in)
       close(1)
    endif

10  format(i3)
40  format(1pe12.5)
20  format(t2,a,i3)
30  format(t2,a,1pe12.5)

  end subroutine interfacelocaldump
  
  
end module neo_interface
