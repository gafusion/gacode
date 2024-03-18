!------------------------------------------------------------
! tgyro_tglf_map.f90
!
! PURPOSE:
!  Mapping from TGYRO internal variables to TGLF interface.
!------------------------------------------------------------

subroutine qlgyro_tglf_map

  use qlgyro_globals
  use tglf_interface
  USE tglf_kyspectrum
  
  implicit none

  ! Local variables
  integer :: i0
  real :: q_abs
  real :: q_prime
  real :: p_prime

  q_abs = abs(q)

  ! Initialize TGLF
  call tglf_init('.', qlgyro_comm)

  ! Are we just testing (to get dump file, for example)?
  tglf_test_flag_in = gyrotest_flag

  ! Want fluxes from TGLF
  tglf_use_transport_model_in = .true.

  tglf_path_in='./'
  !----------------------------------------------------------------
  ! Signs of toroidal magnetic field and current 
  tglf_sign_bt_in = -1.0*signb
  tglf_sign_it_in = -1.0*signb*signq
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! No. of modes
  tglf_nmodes_in = n_modes

  !----------------------------------------------------------------
  ! Sat rule
  tglf_sat_rule_in = sat_rule

  !----------------------------------------------------------------
  ! Ky grid model

  tglf_kygrid_model_in = kygrid_model
  tglf_ky_in = ky_in
  tglf_nky_in = n_ky
  
  !----------------------------------------------------------------
  ! No. of fields
  if (n_field .eq. 2) tglf_use_bper_in = .true.
  if (n_field .eq. 3) tglf_use_bpar_in = .true.
  
  !----------------------------------------------------------------
  ! Number of species (max=6)
  tglf_ns_in = n_ion+1

  if (tglf_ns_in > nsm) then
     call qlgyro_catch_error('ERROR: (qlgyro_tglf_map) Too many ions in TGLF.')
  endif
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Species loop:
  !
  ! Charges: e,i,z
  tglf_zs_in(1) = -1.0
  !
  ! Mass ratios: me/mi(1),mi(1)/mi(1),m(2)/mi(1), ... 
  !              [assume mi(1) is normalizing mass]
  tglf_mass_in(1) = me
  !
  ! Density ratios: ne/ne,ni(1)/ne,ni(2)/ne, ...
  tglf_as_in(1) = ne
  !
  ! Density gradients (e,i,z)
  tglf_rlns_in(1) = dlnnedr
  !
  ! Temperature gradients (e,i,z)
  tglf_rlts_in(1) = dlntedr
  !
  ! Temperature ratios: Te/Te,Ti(1)/Te,Ti(2)/Te
  tglf_taus_in(1) = te

  i0 = 1
  do i_ion=1,n_ion
     i0 = i0+1 
     tglf_zs_in(i0)   = zi_vec(i_ion)
     tglf_mass_in(i0) = mi(i_ion)
     tglf_as_in(i0)   = ni(i_ion)
     tglf_rlns_in(i0) = dlnnidr(i_ion)
     tglf_rlts_in(i0) = dlntidr(i_ion)
     tglf_taus_in(i0) = ti(i_ion)
  enddo
  
 
  !----------------------------------------------------------------
  !   Debye length/rhos   te in ev, rho_s in cm ne in 10^13/cm^3
 ! tglf_debye_in = 7.43e2*sqrt(te/(ne))/abs(rho_s)

  !----------------------------------------------------------------
  ! TGLF-specific quantities
  q_prime = (q_abs/r_min)**2*s
  p_prime = (q_abs/r_min)*(beta_unit/(8*pi))*(-dlnpdr)
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Geometry parameters:
  !
  ! s-alpha (not really needed)
  tglf_rmin_sa_in     = r_min
  tglf_rmaj_sa_in     = r_maj
  tglf_q_sa_in        = q_abs
  tglf_shat_sa_in     = s
  tglf_alpha_sa_in    = r_maj*beta_unit*dlnpdr*q_abs**2
  tglf_xwell_sa_in    = 0.0
  tglf_theta0_sa_in   = 0.0
  tglf_b_model_sa_in  = 0
  tglf_ft_model_sa_in = 1

  if (qlgyro_num_equil_flag == 1) then

     ! Numerical (Fourier) shape
     tglf_geometry_flag_in = 2

     tglf_q_fourier_in       = q_abs
     tglf_q_prime_fourier_in = q_prime
     tglf_p_prime_fourier_in = p_prime
     tglf_nfourier_in        = n_fourier_geo
     tglf_fourier_in(:,0:n_fourier_geo) = a_fourier_geo(:,0:n_fourier_geo)

  else

     ! Model (Miller) shape
     tglf_geometry_flag_in = 1 

     tglf_rmin_loc_in    = r_min
     tglf_rmaj_loc_in    = r_maj
     tglf_zmaj_loc_in    = zmag
     tglf_drmajdx_loc_in = shift
     tglf_dzmajdx_loc_in = dzmag
     tglf_kappa_loc_in   = kappa
     tglf_s_kappa_loc_in = s_kappa
     tglf_delta_loc_in   = delta
     tglf_s_delta_loc_in = s_delta
     tglf_zeta_loc_in    = zeta
     tglf_s_zeta_loc_in  = s_zeta
     tglf_q_loc_in       = q_abs
     tglf_q_prime_loc_in = q_prime
     tglf_p_prime_loc_in = p_prime
     tglf_beta_loc_in = betae_unit

  endif
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Electron beta used for electromagnetic calculations
  tglf_betae_in = betae_unit
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! Collisions:
  !
  ! Electron collision frequency
  tglf_xnue_in = nue
  tglf_xnu_model_in = xnu_model
  !
  ! Zeff
  tglf_zeff_in = z_eff
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Gamma_ExB (ExB shearing rate, units of a/cs)
  if (qlgyro_rotation_flag == 1) then
     ! Currently TGLF uses toroidal current as reference direction
     ! Overall minus sign is due to QLGYRO toroidal angle in CW direction
     !tglf_vexb_shear_in    = -tglf_sign_it_in*gamma_eb0*a_min/c_s  
     !tglf_vpar_shear_in(:) = -tglf_sign_it_in*gamma_p0*a_min/c_s
     !tglf_vpar_in(:)       = -tglf_sign_it_in*r_maj*w0/c_s
     ! Currently GAMMA_EB in input.qlgyro reads in TGLF_VEB_SHEAR
     tglf_vexb_shear_in = -tglf_sign_it_in*gamma_exb
     tglf_alpha_quench_in = 1.0
  endif
  !----------------------------------------------------------------
  ! set ky_in to value for n=1
  if (tglf_kygrid_model_in == 3) tglf_ky_in = abs(rho_s*q/r_min)
  !-----------------------------------
  ! Number of high-k modes
  !  nky=12 (default to include ETG)
  !  nky=0  (low-k only)
  !  tglf_nky_in = 12
  !-----------------------------------

  !----------------------------------------------------------------
  ! Linear mode selection
  !
  ! i_branch_tg = -1: use nmodes_tg modes (default for transport)
  ! i_branch_tg =  0: most unstable electron mode (1) and ion mode (2)
  !
  tglf_ibranch_in = -1
  !
  ! Number of modes in transport calculation;
  ! should be 2-4.
  !tglf_nmodes_in = 4
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Width selection parameters
  !
  ! Use bisection method to determine Gaussian width;
  ! bisection is better, faster.
  !  tglf_use_bisection_in = .true.
  !
  ! Number of Hermite functions to determine Gaussian Width (1-2)
  !  tglf_nbasis_min_in = 1
  !
  ! Number of Hermite function in full expansion (4 or more, even)
  !  tglf_nbasis_max_in = 4
  !
  ! Number of Hermite quadrature nodes
  ! tglf_nxgrid_in = 16
  !
  ! Maximum number of widths sampled (could be higher than 21)
  !  tglf_nwidth_in = 21
  !
  ! Bisection search interval; must increase nwidth_tg if
  ! increasing this interval to maintain accuracy:
  !
  !  accuracy ~ (width_max_tg-width_min_tg)/nwidth_tg
  !
  !  tglf_width_min_in = 0.3
  !  tglf_width_in = 1.65
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! CONTROL PARAMETERS
  !
  ! Use adiabatic electrons
  tglf_adiabatic_elec_in = .false.
  !
  ! Compute eikonal (use stored value if false)
  tglf_new_eikonal_in = .true.
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! OTHER DEFAULT CONTROL PARAMETERS
  !
  !  tglf_find_width_in    = .true.
  !  tglf_iflux_in         = .true.
  !  tglf_theta_trapped_in = 0.7
  !  tglf_xnu_factor_in    = 1.0
  !  tglf_debye_factor_in  = 1.0
  !  tglf_filter_in        = 2.0
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! NEW TGLF SETTINGS
  !
  ! select case (qlgyro_tglf_revision)

  ! case(0)

  !    ! use defaults and overwrites

  ! case (1)

  !    ! APS07 

  !    tglf_alpha_quench_in = 1.0
  !    tglf_xnu_model_in    = 1

  ! case (2)

  !    ! Summer 2009

  !    tglf_alpha_quench_in = 1.0
  !    tglf_xnu_model_in    = 2

  ! case (3)

  !    ! IAEA-2012 with spectral shift ExB shear model

  !    tglf_alpha_quench_in = 0.0
  !    tglf_xnu_model_in    = 2
  !    tglf_alpha_e_in = 1.0
  !    tglf_alpha_p_in = 1.0

  ! end select

  
  !----------------------------------------------------------------
  ! DUMP PARAMETERS
  !
  if (qlgyro_tglf_dump_flag == 0) then
     tglf_dump_flag_in   = .false.
  else
     tglf_dump_flag_in = .true.
  endif

  !----------------------------------------------------------------
  ! VERBOSITY
  !
  tglf_quiet_flag_in = .true.

  !----------------------------------------------------------------
  ! TGLFNN ACTIVATION THRESHOLD
  !
  tglf_nn_max_error_in = qlgyro_tglf_nn_max_error

  if (flux_method == 6) then
     
  end if

  call put_units(tglf_units_in)

  call put_signs(tglf_sign_Bt_in,tglf_sign_It_in)

  call put_species(tglf_ns_in, &
       tglf_zs_in, &
       tglf_mass_in)

  call put_kys(tglf_ky_in)

  call put_gaussian_width(tglf_width_in, &
       tglf_width_min_in, &
       tglf_nwidth_in, &
       tglf_find_width_in)

  call put_eikonal(tglf_new_eikonal_in)

  call put_gradients(tglf_rlns_in, &
       tglf_rlts_in, &
       tglf_vpar_shear_in, &
       tglf_vexb_shear_in)

  call put_averages(tglf_taus_in, &
       tglf_as_in, &
       tglf_vpar_in, &
       tglf_vexb_in, &
       tglf_betae_in, &
       tglf_xnue_in, &
       tglf_zeff_in, &
       tglf_debye_in)

  call put_switches(tglf_iflux_in, &
       tglf_use_bper_in, &
       tglf_use_bpar_in, &
       tglf_use_mhd_rule_in, &
       tglf_use_bisection_in, &
       tglf_use_inboard_detrapped_in, &
       tglf_ibranch_in, &
       tglf_nmodes_in, &
       tglf_nbasis_max_in, &
       tglf_nbasis_min_in, &
       tglf_nxgrid_in, &
       tglf_nky_in, &
       tglf_use_ave_ion_grid_in)  

  call put_rare_switches(tglf_theta_trapped_in, &
       tglf_wdia_trapped_in, &
       tglf_park_in, &
       tglf_ghat_in, &
       tglf_gchat_in, &
       tglf_wd_zero_in, &
       tglf_linsker_factor_in, &
       tglf_gradB_factor_in, &
       tglf_filter_in, &
       tglf_damp_psi_in, &
       tglf_damp_sig_in)

  call put_model_parameters(tglf_adiabatic_elec_in, &
       tglf_alpha_e_in, &
       tglf_alpha_p_in, &
       tglf_alpha_mach_in, &
       tglf_alpha_quench_in, &
       tglf_alpha_zf_in, &
       tglf_xnu_factor_in, &
       tglf_debye_factor_in, &
       tglf_etg_factor_in, &
       tglf_rlnp_cutoff_in, &
       tglf_sat_rule_in, &
       tglf_kygrid_model_in, &
       tglf_xnu_model_in, &
       tglf_vpar_model_in, &
       tglf_vpar_shear_model_in)

  tglf_kygrid_model_in = kygrid_model

  if (tglf_geometry_flag_in == 1 ) then

     call put_miller_geometry(tglf_rmin_loc_in, &
          tglf_rmaj_loc_in, &
          tglf_zmaj_loc_in, &
          tglf_drmindx_loc_in, &
          tglf_drmajdx_loc_in, &
          tglf_dzmajdx_loc_in, &
          tglf_kappa_loc_in, &
          tglf_s_kappa_loc_in, &
          tglf_delta_loc_in, &
          tglf_s_delta_loc_in, &
          tglf_zeta_loc_in, &
          tglf_s_zeta_loc_in, &
          tglf_q_loc_in, &
          tglf_q_prime_loc_in, &
          tglf_p_prime_loc_in, &
          tglf_beta_loc_in,  &
          tglf_kx0_loc_in)

  elseif (tglf_geometry_flag_in == 2)then

     call put_fourier_geometry(tglf_q_fourier_in,  &
          tglf_q_prime_fourier_in, &
          tglf_p_prime_fourier_in, &
          tglf_nfourier_in, &
          tglf_fourier_in)

  else

     call put_s_alpha_geometry(tglf_rmin_sa_in, &
          tglf_rmaj_sa_in, &
          tglf_q_sa_in, &
          tglf_shat_sa_in, &
          tglf_alpha_sa_in, &
          tglf_xwell_sa_in, &
          tglf_theta0_sa_in, &
          tglf_b_model_sa_in, &
          tglf_ft_model_sa_in)

  endif

  call tglf_startup

end subroutine qlgyro_tglf_map
