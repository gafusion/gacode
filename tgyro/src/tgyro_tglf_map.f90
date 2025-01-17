
!------------------------------------------------------------
! tgyro_tglf_map.f90
!
! PURPOSE:
!  Mapping from TGYRO internal variables to TGLF interface.
!------------------------------------------------------------

subroutine tgyro_tglf_map

  use tgyro_globals
  use tglf_interface

  implicit none

  ! Local variables
  integer :: i_ion,i0
  real :: q_abs
  real :: q_prime
  real :: p_prime
  real :: beta
  real :: gamma_eb0
  real :: gamma_p0

  real :: shape_fac
  
  q_abs = abs(q(i_r))

  ! Initialize TGLF
  call tglf_init(paths(i_r-1), gyro_comm)

  ! Are we just testing (to get dump file, for example)?
  tglf_test_flag_in = gyrotest_flag

  ! Want fluxes from TGLF
  tglf_use_transport_model_in = .true.

  !----------------------------------------------------------------
  ! Signs of toroidal magnetic field and current 
  tglf_sign_bt_in = -1.0*signb
  tglf_sign_it_in = -1.0*signb*signq
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Number of species (max=6)
  tglf_ns_in = sum(calc_flag(1:loc_n_ion))+1

  if (tglf_ns_in > nsm) then
     call tgyro_catch_error('ERROR: (tgyro_tglf_map) Too many ions in TGLF.')
  endif
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Species loop:
  !
  ! Charges: e,i,z
  tglf_zs_in(1) = -1.0
  !
  ! Mass ratios: me/md,m(1)/md,m(2)/md, ... 
  !              [assume md is normalizing mass]
  tglf_mass_in(1) = (me*loc_me_multiplier)/md
  !
  ! Density ratios: ne/ne,ni(1)/ne,ni(2)/ne, ...
  tglf_as_in(1) = 1.0
  !
  ! Density gradients (e,i,z)
  tglf_rlns_in(1) = r_min*dlnnedr(i_r)
  !
  ! Temperature gradients (e,i,z)
  tglf_rlts_in(1) = r_min*dlntedr(i_r)
  !
  ! Temperature ratios: Te/Te,Ti(1)/Te,Ti(2)/Te
  tglf_taus_in(1) = 1.0

  i0 = 1
  do i_ion=1,loc_n_ion
     if (calc_flag(i_ion) == 0) cycle
     i0 = i0+1 
     tglf_zs_in(i0)   = zi_vec(i_ion)
     tglf_mass_in(i0) = mi(i_ion)/md
     tglf_as_in(i0)   = ni(i_ion,i_r)/ne(i_r) 
     tglf_rlns_in(i0) = r_min*dlnnidr(i_ion,i_r)
     tglf_rlts_in(i0) = r_min*dlntidr(i_ion,i_r)
     tglf_taus_in(i0) = ti(i_ion,i_r)/te(i_r)
  enddo

  ! Setting density gradient artificially to zero to compute D and v
  if (tgyro_zero_dens_grad_flag /= 0) then
     tglf_rlns_in(tgyro_zero_dens_grad_flag) = 0
  endif

  !----------------------------------------------------------------
  !   debye length/rhos   te in ev, rho_s in cm ne in 10^13/cm^3
  tglf_debye_in = 7.43e2*sqrt(te(i_r)/(ne(i_r)))/abs(rho_s(i_r))

  !----------------------------------------------------------------
  ! TGLF-specific quantities
  q_prime = (q_abs/(r(i_r)/r_min))**2*s(i_r)
  if(tgyro_tglf_ptot_flag == 1)then
    beta = beta_unit(i_r)*ptot(i_r)/pr(i_r)
    p_prime = (q_abs/(r(i_r)/r_min))*(beta/(8*pi))*(-r_min*dlnptotdr(i_r))
  else
    beta = 0.0 ! tglf will compute its own dlnpdr
    p_prime = (q_abs/(r(i_r)/r_min))*(beta_unit(i_r)/(8*pi))*(-r_min*dlnpdr(i_r))
  endif
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Geometry parameters:
  !
  ! s-alpha (not really needed)
  tglf_rmin_sa_in     = r(i_r)/r_min
  tglf_rmaj_sa_in     = r_maj(i_r)/r_min
  tglf_q_sa_in        = q_abs
  tglf_shat_sa_in     = s(i_r)
  tglf_alpha_sa_in    = r_maj(i_r)*beta_unit(i_r)*dlnpdr(i_r)*q_abs**2
  tglf_xwell_sa_in    = 0.0
  tglf_theta0_sa_in   = 0.0
  tglf_b_model_sa_in  = 0
  tglf_ft_model_sa_in = 1

  ! Model (Miller) shape
  ! Force Miller geometry in TGLF for now
  ! In subsequent edits, maybe propagate the XMH flag to TGLF
  ! Match the XMH flag to an available tglf_geometry_flag
  tglf_geometry_flag_in = 1 
  select case( tgyro_tglf_mxh_flag )
  case ( 1 )
     shape_fac = 1.0 !Strictly Miller geometry for default geometry flag == 1
  case default
     shape_fac = 0.0
  end select

  tglf_rmin_loc_in    = r(i_r)/r_min
  tglf_rmaj_loc_in    = r_maj(i_r)/r_min
  tglf_zmaj_loc_in    = zmag(i_r)/r_min
  tglf_drmajdx_loc_in = shift(i_r)
  tglf_dzmajdx_loc_in = dzmag(i_r)
  tglf_kappa_loc_in   = kappa(i_r)
  tglf_s_kappa_loc_in = s_kappa(i_r)
  tglf_delta_loc_in   = delta(i_r)
  tglf_s_delta_loc_in = s_delta(i_r)
  tglf_zeta_loc_in    = zeta(i_r)
  tglf_s_zeta_loc_in  = s_zeta(i_r)
  ! eXtended Miller Harmonic shape coefficients

  tglf_shape_sin3_loc_in   = shape_fac * shape_sin3(i_r)
  tglf_shape_s_sin3_loc_in = shape_fac * shape_ssin3(i_r)
!  tglf_shape_sin4_loc_in   = 0.0 ! XMH coefficients defined up to n=3 in tgyro_globals.f90
!  tglf_shape_s_sin4_loc_in = 0.0 ! leave additional coefficients available within TGLF
!  tglf_shape_sin5_loc_in   = 0.0 ! with their TGLF default values
!  tglf_shape_s_sin5_loc_in = 0.0
!  tglf_shape_sin6_loc_in   = 0.0
!  tglf_shape_s_sin6_loc_in = 0.0
  tglf_shape_cos0_loc_in   = shape_fac * shape_cos0(i_r)
  tglf_shape_s_cos0_loc_in = shape_fac * shape_scos0(i_r)
  tglf_shape_cos1_loc_in   = shape_fac * shape_cos1(i_r)
  tglf_shape_s_cos1_loc_in = shape_fac * shape_scos1(i_r)
  tglf_shape_cos2_loc_in   = shape_fac * shape_cos2(i_r)
  tglf_shape_s_cos2_loc_in = shape_fac * shape_scos2(i_r)
  tglf_shape_cos3_loc_in   = shape_fac * shape_cos3(i_r)
  tglf_shape_s_cos3_loc_in = shape_fac * shape_scos3(i_r)
!  tglf_shape_cos4_loc_in   = 0.0  ! XMH coefficients defined up to n=3 in tgyro_globals.f90
!  tglf_shape_s_cos4_loc_in = 0.0
!  tglf_shape_cos5_loc_in   = 0.0
!  tglf_shape_s_cos5_loc_in = 0.0
!  tglf_shape_cos6_loc_in   = 0.0
!  tglf_shape_s_cos6_loc_in = 0.0
  !
  tglf_q_loc_in       = q_abs
  tglf_q_prime_loc_in = q_prime
  tglf_p_prime_loc_in = p_prime
  tglf_beta_loc_in = beta
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Electron beta used for electromagnetic calculations
  tglf_betae_in = betae_unit(i_r)*loc_betae_scale
  !----------------------------------------------------------------
  !----------------------------------------------------------------
  ! Collisions:
  !
  ! Electron collision frequency
  tglf_xnue_in = nue(i_r)*r_min/c_s(i_r)*loc_nu_scale
  !
  ! Zeff
  tglf_zeff_in = z_eff(i_r)
  !----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Gamma_ExB (ExB shearing rate, units of a/cs)
  if (tgyro_rotation_flag == 1) then
     gamma_p0  = -r_maj(i_r)*f_rot(i_r)*w0_norm
     gamma_eb0 = gamma_p0*r(i_r)/(q_abs*r_maj(i_r)) 
     ! Currently TGLF uses toroidal current as reference direction
     ! Overall minus sign is due to TGYRO toroidal angle in CW direction
     tglf_vexb_shear_in    = -tglf_sign_it_in*gamma_eb0*r_min/c_s(i_r)  
     tglf_vpar_shear_in(:) = -tglf_sign_it_in*gamma_p0*r_min/c_s(i_r)
     tglf_vpar_in(:)       = -tglf_sign_it_in*r_maj(i_r)*w0(i_r)/c_s(i_r)
  endif
  !----------------------------------------------------------------
  ! set ky_in to value for n=1
  if (tglf_kygrid_model_in == 3) tglf_ky_in = abs(rho_s(i_r)*q(i_r)/r(i_r))
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
  if (loc_betae_scale <= 0.0) then
     tglf_use_bper_in = .false.
     tglf_use_bpar_in = .false.
  endif
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
  select case (tgyro_tglf_revision)

  case(0)

     ! use defaults and overwrites

  case (1)

     ! Original 1-D saturation model sat_rule = 0 (SAT0) with the Waltz quench rule for ExB velocity shear and the innacurate 
     ! electron collision model
           !   J. E. Kinsey, G. M. Staebler, and R. E. Waltz, “The First Transport Code Simulations using the Trapped Gyro-Landau Fluid Transport Model”, 
           !    Phys. Plasmas 15, 055908 (2008).

     tglf_sat_rule_in = 0
     tglf_units_in = 'GYRO'
     tglf_alpha_quench_in = 1.0
     tglf_xnu_model_in    = 1

  case (2)

           !  Improved 2-D (kx,ky)saturation model with multi-scale ETG zonal flow mixing: sat_rule = 1  (SAT1)
           !   G. M. Staebler, N. T. Howard, J. Candy, and C. Holland, "A model of the saturation of coupled electron and ion scale gyrokinetic turbulence", Nucl. Fusion 57 (2017) 066046.
     ! Improved electron collision model (now the default) 
           !   G. M. Staebler and J. E. Kinsey, “Electron Collisions in the Trapped Gyro-Landau Fluid Transport Model”, Phys. Plasmas 17, 122309 (2010) 
     ! Turn on the "spectral shift" ExB shear model and turn off the Waltz quench rule
     ! G. Staebler, R. Waltz, J. Candy, and J. Kinsey, "A new paradigm for suppression of gyrokinetic turbulence by velocity shear" ,Phys. Rev. Lett. 110, 055003 (2013)

      tglf_sat_rule_in = 1
      tglf_xnu_model_in = 2
      tglf_units_in = 'GYRO'
      tglf_alpha_quench_in = 0.0
      tglf_alpha_e_in = 1.0

  case (3)

     ! SAT_RULE=2 (SAT2) settings determined for JET DTE2 experiments.
     ! These are the recommended setting for using SAT_RULE =2
     ! Note that SAT2 is a 3-D fit (kx,ky,theta) to the saturated potential fluctuation intensity from 64 CGYRO runs
     ! SAT2 was fit to CGYRO using the most unstable CGYRO linear eigenmode spectrum not TGLF. 
           !   G. Staebler, J. Candy, E. Belli, J. Kinsey, N. Bonanomi, and B. Patel. Geometry dependence of the fluctuation intensity in gyrokinetic turbulence. 
           !   Plasma Phys. Control. Fusion 63, 015013 (2020). doi.org/10.1088/1361-6587/abc861
     ! G. Staebler, E. Belli, J. Candy, J. Kinsey, H. Dudding, and B. Patel. "Verification of a quasi-linear model for gyrokinetic turbulent transport."
           !   Nucl. Fusion 61, 116007, (2021). doi.org/10.1088/1361-6587/abc861

     tglf_sat_rule_in = 2
     tglf_units_in = 'CGYRO'
     tglf_alpha_quench_in = 0.0
     tglf_alpha_e_in = 1.0
     tglf_alpha_p_in = 1.0
     tglf_alpha_mach_in = 0.0
     tglf_use_bper_in = .true.
     tglf_use_bpar_in = .false.
     tglf_use_ave_ion_grid_in = .true.
     tglf_kygrid_model_in = 4
     tglf_nbasis_max_in = 6
     tglf_nmodes_in = 8
     tglf_geometry_flag_in = 1
     tglf_use_mhd_rule_in = .false.
     tglf_nky_in = 18
 
  case (4)

     ! Momentum transport without EM terms.
     ! due to a bug in TGLF that shows up when you have both parallel flows (alpha_mach=1) and EM flutter use_bper = T 
     ! the parallel flow should not be included for EM runs. However, parallel flows give an important pinch for momentum transport 
     ! so for low-beta plasmas parallel flow should be included. These setting should be used for such cases with all saturation rules

     tglf_alpha_quench_in = 0.0
     tglf_alpha_e_in = 1.0
     tglf_alpha_p_in = 1.0
     tglf_alpha_mach_in = 1.0
     tglf_use_bper_in = .false.
     tglf_use_bpar_in = .false.

 
  end select

  !----------------------------------------------------------------
  ! DUMP PARAMETERS
  !
  if (tgyro_tglf_dump_flag == 0) then
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
  tglf_nn_max_error_in = tgyro_tglf_nn_max_error

end subroutine tgyro_tglf_map
