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
  integer :: harvest_err
  real :: q_abs
  real :: q_prime
  real :: p_prime
  real :: gamma_eb0
  real :: gamma_p0

  include 'harvest_lib.inc'

  CHARACTER NUL
  PARAMETER(NUL = CHAR(0))

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
  if (tgyro_quickfast_flag == 1) then
     tglf_ns_in = sum(therm_flag(1:loc_n_ion))+1
  else
     tglf_ns_in = loc_n_ion+1
  endif
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
  ! Mass ratios: me/mi(1),mi(1)/mi(1),m(2)/mi(1), ... 
  !              [assume mi(1) is normalizing mass]
  tglf_mass_in(1) = (me*loc_me_multiplier)/mi(1)
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
     if (therm_flag(i_ion) == 0 .and. tgyro_quickfast_flag == 1) cycle
     i0 = i0+1 
     tglf_zs_in(i0)   = zi_vec(i_ion)
     tglf_mass_in(i0) = mi(i_ion)/mi(1)
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
  tglf_debye_in = 7.43D2*SQRT(te(i_r)/(ne(i_r)*1.D13))/ABS(rho_s(i_r))

  !----------------------------------------------------------------
  ! TGLF-specific quantities
  q_prime = (q_abs/(r(i_r)/r_min))**2*s(i_r)
  p_prime = (q_abs/(r(i_r)/r_min))*(beta_unit(i_r)/(8*pi))*(-r_min*dlnpdr(i_r))
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

  if (loc_num_equil_flag == 1) then

     ! Numerical (Fourier) shape
     tglf_geometry_flag_in = 2

     tglf_q_fourier_in       = q_abs
     tglf_q_prime_fourier_in = q_prime
     tglf_p_prime_fourier_in = p_prime
     tglf_nfourier_in        = n_fourier_geo
     tglf_fourier_in(:,0:n_fourier_geo) = a_fourier_geo(:,0:n_fourier_geo,i_r)

  else

     ! Model (Miller) shape
     tglf_geometry_flag_in = 1 

     tglf_rmin_loc_in    = r(i_r)/r_min
     tglf_rmaj_loc_in    = r_maj(i_r)/r_min
     tglf_zmaj_loc_in    = zmag(i_r)
     tglf_drmajdx_loc_in = shift(i_r)
     tglf_dzmajdx_loc_in = dzmag(i_r)
     tglf_kappa_loc_in   = kappa(i_r)
     tglf_s_kappa_loc_in = s_kappa(i_r)
     tglf_delta_loc_in   = delta(i_r)
     tglf_s_delta_loc_in = s_delta(i_r)
     tglf_zeta_loc_in    = zeta(i_r)
     tglf_s_zeta_loc_in  = s_zeta(i_r)
     tglf_q_loc_in       = q_abs
     tglf_q_prime_loc_in = q_prime
     tglf_p_prime_loc_in = p_prime

  endif
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
  tglf_nmodes_in = 2
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
  if (loc_betae_scale == 0.0) then
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

     ! APS07 

     tglf_alpha_quench_in = 1.0
     tglf_xnu_model_in    = 1

  case (2)

     ! Summer 2009

     tglf_alpha_quench_in = 1.0
     tglf_xnu_model_in    = 2

  case (3)

     ! IAEA-2012 with spectral shift ExB shear model

     tglf_alpha_quench_in = 0.0
     tglf_xnu_model_in    = 2
     tglf_alpha_e_in = 1.0
     tglf_alpha_p_in = 1.0

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

  !----------------------------------------------------------------
  ! HARVEST: NEO AND TARGET FLUXES, GYRO-BOHM NORMALIZATIONS, SHOT
  !
  if (i_tran == 0 .and. 0==1) then
     ! Initialization
     tglf_harvest_extra_in = NUL
     harvest_err=set_harvest_verbose(0)

     ! Target fluxes
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_eflux_e_target'//NUL,eflux_e_target(i_r))
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_eflux_i_target'//NUL,eflux_i_target(i_r))
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_pflux_e_target'//NUL,pflux_e_target(i_r))
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_mflux_target'//NUL,mflux_target(i_r))

     ! Neoclassical fluxes
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_eflux_e_neo'//NUL,eflux_e_neo(i_r))
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_sum_eflux_i_neo'//NUL,&
          & sum(eflux_i_neo(therm_vec(:),i_r)))
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_pflux_e_neo'//NUL,pflux_e_neo(i_r))
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_sum_mflux_i_neo'//NUL,&
          & sum(mflux_i_neo(therm_vec(:),i_r)))

     ! Turbulent fluxes interface
     CALL tglf_harvest_local

     ! Gyrobohm normalizations
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_q_gb'//NUL,q_gb(i_r))
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_pi_gb'//NUL,pi_gb(i_r))
     harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_gamma_gb'//NUL,gamma_gb(i_r))

     ! Indication of thermal ions
     harvest_err=set_harvest_payload_int_array(tglf_harvest_extra_in,'tgyro_therm_vec'//NUL,therm_vec(:),size(therm_vec))

     ! Experimental shot
     harvest_err=set_harvest_payload_int(tglf_harvest_extra_in,'shot'//NUL,shot)
  endif

end subroutine tgyro_tglf_map
