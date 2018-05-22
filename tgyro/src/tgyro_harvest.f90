subroutine tgyro_harvest

  use tglf_pkg
  use tgyro_globals
  use tglf_interface
  
  implicit none
  
  integer :: harvest_err

  character NUL
  parameter(NUL = char(0))

  include 'harvest_lib.inc'

  !----------------------------------------------------------------
  ! HARVEST: NEO AND TARGET FLUXES, GYRO-BOHM NORMALIZATIONS, SHOT
  !
  ! Initialization
  tglf_harvest_extra_in = NUL
  harvest_err=set_harvest_verbose(0)

  ! Target fluxes
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,&
       'tgyro_eflux_e_target'//NUL,eflux_e_target(i_r))
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,&
       'tgyro_eflux_i_target'//NUL,eflux_i_target(i_r))
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,&
       'tgyro_pflux_e_target'//NUL,pflux_e_target(i_r))
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,&
       'tgyro_mflux_target'//NUL,mflux_target(i_r))

  ! Neoclassical fluxes
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,&
       'tgyro_eflux_e_neo'//NUL,eflux_e_neo(i_r))
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,&
       'tgyro_sum_eflux_i_neo'//NUL,&
       & sum(eflux_i_neo(therm_vec(:),i_r)))
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,&
       'tgyro_pflux_e_neo'//NUL,pflux_e_neo(i_r))
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,&
       'tgyro_sum_mflux_i_neo'//NUL,&
       & sum(mflux_i_neo(therm_vec(:),i_r)))

  ! Turbulent fluxes interface
  call tglf_harvest_local

  ! Gyrobohm normalizations
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_q_gb'//NUL,q_gb(i_r))
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_pi_gb'//NUL,pi_gb(i_r))
  harvest_err=set_harvest_payload_dbl(tglf_harvest_extra_in,'tgyro_gamma_gb'//NUL,gamma_gb(i_r))

  ! Indication of thermal ions
  harvest_err=set_harvest_payload_int_array(tglf_harvest_extra_in,&
       'tgyro_therm_vec'//NUL,therm_vec(:),size(therm_vec))

  ! Experimental shot
  harvest_err=set_harvest_payload_int(tglf_harvest_extra_in,'shot'//NUL,shot)

end subroutine tgyro_harvest


subroutine tglf_harvest_local()

  ! Harvest LOCAL INTERFACE variables
  !
  ! _bol : boolean
  ! _int : integer
  ! _dbl : double

  use tglf_interface
  use tglf_pkg

  integer :: ierr, i, j, nky
  character(LEN=65507) :: harvest_sendline
  character(LEN=2) :: NUM
  character NUL
  parameter(NUL = char(0))

  real, dimension(10) :: tmp
  integer, dimension(10) :: ions_order

  real, dimension(:), allocatable :: spectrum

  if (.not.tglf_use_transport_model_in) then
     write(1,*) 'HARVEST ONLY WHEN `TRANSPORT_MODEL=.TRUE.`' !only when computing fluxes
     return
  endif

  if (.not.tglf_iflux_in) then
     write(1,*) 'HARVEST ONLY AVAILABLE WHEN `IFLUX=.TRUE.' !only when computing fluxes
     return
  endif

  if (tglf_geometry_flag_in /= 1) then
     write(1,*) 'HARVEST ONLY SUPPORTS MILLER `GEOMETRY_FLAG=1`' !only when using miller
     return
  endif

  !table is set according to HARVEST_TABLE environmental variable
  ierr=init_harvest(NUL,harvest_sendline,len(harvest_sendline))

  !no underscore to allow different versions of the same run
  ierr=set_harvest_payload_str(harvest_sendline,'VERSION'//NUL,'db9432992b7d457b91bebb0d5237f5cc5551dde0'//NUL) 

  !   ---------------------------------------------------
  !    Plasma parameters
  !   ---------------------------------------------------
  ierr=set_harvest_payload_int(harvest_sendline,'SIGN_BT'//NUL,tglf_sign_bt_in)
  ierr=set_harvest_payload_int(harvest_sendline,'SIGN_IT'//NUL,tglf_sign_it_in)
  if (tglf_kygrid_model_in.ne.1) then
     ierr=set_harvest_payload_dbl(harvest_sendline,'KY'//NUL,tglf_ky_in)
  endif
  ierr=set_harvest_payload_dbl(harvest_sendline,'VEXB_SHEAR'//NUL,tglf_vexb_shear_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'BETAE'//NUL,tglf_betae_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'XNUE'//NUL,tglf_xnue_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'ZEFF'//NUL,tglf_zeff_in)
  if (tglf_debye_in.ne.0) then
     ierr=set_harvest_payload_dbl(harvest_sendline,'DEBYE'//NUL,tglf_debye_in)
  endif

  !   ---------------------------------------------------
  !    Sort ions by A, Z, Te/Ti
  !    electrons always first (as required by TGLF)
  !    !main ion always second (affects GB normalization)
  !   ---------------------------------------------------
  tmp=0.0
  tmp(1)=1e10
  tmp(2)=1e9
  ions_order=0
  do i=3,tglf_ns_in
     tmp(i)=tglf_mass_in(i)*100000+tglf_zs_in(i)*1000+1./(1+tglf_rlts_in(i))
  enddo
  do i=1,tglf_ns_in
     j=maxloc(tmp,dim=1)
     ions_order(i)=j
     tmp(j)=0
  enddo

  !   ---------------------------------------------------
  !    Species vectors
  !   ---------------------------------------------------
  do i = 1,tglf_ns_in
     if (i < 10) then
        write (NUM, "(I01,A1)") i,NUL
     else
        write (NUM, "(I02,A1)") i,NUL
     endif

     j=ions_order(i)

     ierr=set_harvest_payload_dbl(harvest_sendline,&
          '+ZS_'//NUM,tglf_zs_in(j))     !make new tables for different charges
     ierr=set_harvest_payload_dbl(harvest_sendline,&
          '+MASS_'//NUM,tglf_mass_in(j)) !make new tables for different masses
     ierr=set_harvest_payload_dbl(harvest_sendline,'RLNS_'//NUM,tglf_rlns_in(j))
     ierr=set_harvest_payload_dbl(harvest_sendline,'RLTS_'//NUM,tglf_rlts_in(j))
     if (i /= 1) then !these are always 1, by definition
        ierr=set_harvest_payload_dbl(harvest_sendline,'TAUS_'//NUM,tglf_taus_in(j))
        ierr=set_harvest_payload_dbl(harvest_sendline,'AS_'//NUM,tglf_as_in(j))
     endif
     ierr=set_harvest_payload_dbl(harvest_sendline,'VPAR_'//NUM,tglf_vpar_in(j))
     ierr=set_harvest_payload_dbl(harvest_sendline,'VPAR_SHEAR_'//NUM,tglf_vpar_shear_in(j))
  enddo

  !   ---------------------------------------------------
  !    Miller geometry parameters
  !   ---------------------------------------------------
  ierr=set_harvest_payload_dbl(harvest_sendline,'RMIN_LOC'//NUL,tglf_rmin_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'RMAJ_LOC'//NUL,tglf_rmaj_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'ZMAJ_LOC'//NUL,tglf_zmaj_loc_in)
  if (tglf_drmindx_loc_in /= 1) then
     ierr=set_harvest_payload_dbl(harvest_sendline,&
          'DRMINDX_LOC'//NUL,tglf_drmindx_loc_in) ! different derivatives affect everything!
  endif
  ierr=set_harvest_payload_dbl(harvest_sendline,'DRMAJDX_LOC'//NUL,tglf_drmajdx_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'DZMAJDX_LOC'//NUL,tglf_dzmajdx_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'Q_LOC'//NUL,tglf_q_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'KAPPA_LOC'//NUL,tglf_kappa_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'S_KAPPA_LOC'//NUL,tglf_s_kappa_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'DELTA_LOC'//NUL,tglf_delta_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'S_DELTA_LOC'//NUL,tglf_s_delta_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'ZETA_LOC'//NUL,tglf_zeta_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'S_ZETA_LOC'//NUL,tglf_s_zeta_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'P_PRIME_LOC'//NUL,tglf_p_prime_loc_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'Q_PRIME_LOC'//NUL,tglf_q_prime_loc_in)
  if (tglf_kx0_loc_in /= 0) then
     ierr=set_harvest_payload_dbl(harvest_sendline,&
          'KX0_LOC'//NUL,tglf_kx0_loc_in) !should be always 0.0
  endif

  !   ---------------------------------------------------
  !    Gaussian width parameters
  !   ---------------------------------------------------
  ierr=set_harvest_payload_dbl(harvest_sendline,'+WIDTH'//NUL,tglf_width_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+WIDTH_MIN'//NUL,tglf_width_min_in)
  ierr=set_harvest_payload_bol(harvest_sendline,'+FIND_WIDTH'//NUL,tglf_find_width_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+NWIDTH'//NUL,tglf_nwidth_in)

  !   ---------------------------------------------------
  !    Control parameters
  !   ---------------------------------------------------
  ierr=set_harvest_payload_bol(harvest_sendline,'+USE_TRANSPORT_MODEL'//NUL,tglf_use_transport_model_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+GEOMETRY_FLAG'//NUL,tglf_geometry_flag_in)
  ierr=set_harvest_payload_bol(harvest_sendline,'+IFLUX'//NUL,tglf_iflux_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+NS'//NUL,tglf_ns_in)
  ierr=set_harvest_payload_bol(harvest_sendline,'+ADIABATIC_ELEC'//NUL,tglf_adiabatic_elec_in)
  ierr=set_harvest_payload_bol(harvest_sendline,'+USE_BPER'//NUL,tglf_use_bper_in)
  ierr=set_harvest_payload_bol(harvest_sendline,'+USE_BPAR'//NUL,tglf_use_bpar_in)
  ierr=set_harvest_payload_bol(harvest_sendline,'+USE_MHD_RULE'//NUL,tglf_use_mhd_rule_in)
  ierr=set_harvest_payload_bol(harvest_sendline,'+USE_BISECTION'//NUL,tglf_use_bisection_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+SAT_RULE'//NUL,tglf_sat_rule_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+KYGRID_MODEL'//NUL,tglf_kygrid_model_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+XNU_MODEL'//NUL,tglf_xnu_model_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+VPAR_MODEL'//NUL,tglf_vpar_model_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+VPAR_SHEAR_MODEL'//NUL,tglf_vpar_shear_model_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+IBRANCH'//NUL,tglf_ibranch_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+NMODES'//NUL,tglf_nmodes_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+NBASIS_MAX'//NUL,tglf_nbasis_max_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+NBASIS_MIN'//NUL,tglf_nbasis_min_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+NXGRID'//NUL,tglf_nxgrid_in)
  ierr=set_harvest_payload_int(harvest_sendline,'+NKY'//NUL,tglf_nky_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_MACH'//NUL,tglf_alpha_mach_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_E'//NUL,tglf_alpha_e_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_P'//NUL,tglf_alpha_p_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_QUENCH'//NUL,tglf_alpha_quench_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+ALPHA_ZF'//NUL,tglf_alpha_zf_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+XNU_FACTOR'//NUL,tglf_xnu_factor_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+DEBYE_FACTOR'//NUL,tglf_debye_factor_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+ETG_FACTOR'//NUL,tglf_etg_factor_in)

  !   ---------------------------------------------------
  !    Expert parameters
  !   ---------------------------------------------------
  ierr=set_harvest_payload_dbl(harvest_sendline,'+DAMP_PSI'//NUL,tglf_damp_psi_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+DAMP_SIG'//NUL,tglf_damp_sig_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+PARK'//NUL,tglf_park_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+GHAT'//NUL,tglf_ghat_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+GCHAT'//NUL,tglf_gchat_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+WD_ZERO'//NUL,tglf_wd_zero_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+LINSKER_FACTOR'//NUL,tglf_linsker_factor_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+GRADB_FACTOR'//NUL,tglf_gradB_factor_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+FILTER'//NUL,tglf_filter_in)
  ierr=set_harvest_payload_dbl(harvest_sendline,'+THETA_TRAPPED'//NUL,tglf_theta_trapped_in)

  !   ---------------------------------------------------
  !    Disabled harvest
  !   ---------------------------------------------------
  !   ierr=set_harvest_payload_bol(harvest_sendline,'+NEW_EIKONAL'//NUL,tglf_new_eikonal_in) !DISABLED because not used anymore and should not make difference anyways
  !   ierr=set_harvest_payload_dbl(harvest_sendline,'+VEXB'//NUL,tglf_vexb_in) !DISABLED because not used anymore
  !   ierr=set_harvest_payload_int(harvest_sendline,'+WRITE_WAVEFUNCTION_FLAG'//NUL,tglf_write_wavefunction_flag_in)  ! DISABLED because it does not matter for us
  !   ierr=set_harvest_payload_dbl(harvest_sendline,'VNS_SHEAR_'//NUM,tglf_vns_shear_in(j)) !DISABLED because not used anymore
  !   ierr=set_harvest_payload_dbl(harvest_sendline,'VTS_SHEAR_'//NUM,tglf_vts_shear_in(j)) !DISABLED because not used anymore

  !   ---------------------------------------------------
  !    Output results
  !   ---------------------------------------------------
  do i = 1,tglf_ns_in
     if (i < 10) then
        write (NUM, "(I01,A1)") i,NUL
     else
        write (NUM, "(I02,A1)") i,NUL
     endif

     j=ions_order(i)

     ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_PARTICLE_FLUX_ESTATIC_'//NUM,get_particle_flux(j,1))
     ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_ENERGY_FLUX_ESTATIC_'//NUM,get_energy_flux(j,1))
     ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_PAR_ESTATIC_'//NUM,get_stress_par(j,1))
     ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_TOR_ESTATIC_'//NUM,get_stress_tor(j,1))
     ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_EXCHANGE_ESTATIC_'//NUM,get_exchange(j,1))

     if (tglf_use_bper_in) then
        ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_PARTICLE_FLUX_EMPER_'//NUM,get_particle_flux(j,2))
        ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_ENERGY_FLUX_EMPER_'//NUM,get_energy_flux(j,2))
        ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_PAR_EMPER_'//NUM,get_stress_par(j,2))
        ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_TOR_EMPER_'//NUM,get_stress_tor(j,2))
        ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_EXCHANGE_EMPER_'//NUM,get_exchange(j,2))
     endif

     if (tglf_use_bpar_in) then
        ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_PARTICLE_FLUX_EMPAR_'//NUM,get_particle_flux(j,3))
        ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_ENERGY_FLUX_EMPAR_'//NUM,get_energy_flux(j,3))
        ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_PAR_EMPAR_'//NUM,get_stress_par(j,3))
        ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_TOR_EMPAR_'//NUM,get_stress_tor(j,3))
        ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_EXCHANGE_EMPAR_'//NUM,get_exchange(j,3))
     endif

  enddo

  !   ---------------------------------------------------
  !    Spectra
  !   ---------------------------------------------------
  nky = get_nky_out()
  allocate(spectrum(nky))

  do i=1,nky
     spectrum(i) = get_ky_spectrum_out(i)
  enddo

  ierr=set_harvest_payload_dbl_array(harvest_sendline,'KY_SPECTRUM'//NUL,spectrum,nky)

  do i=1,tglf_nmodes_in
     if (i < 10) then
        write (NUM, "(I01,A1)") i,NUL
     else
        write (NUM, "(I02,A1)") i,NUL
     endif

     do j=1,nky
        spectrum(j) = get_eigenvalue_spectrum_out(1,j,i)
     enddo
     ierr=set_harvest_payload_dbl_array(harvest_sendline,'OUT_EIGENVALUE_SPECTRUM_GAMMA'//NUM,spectrum,nky)
     do j=1,nky
        spectrum(j) = get_eigenvalue_spectrum_out(2,j,i)
     enddo
     ierr=set_harvest_payload_dbl_array(harvest_sendline,'OUT_EIGENVALUE_SPECTRUM_OMEGA'//NUM,spectrum,nky)
  enddo
  deallocate(spectrum)

  !   ---------------------------------------------------
  !    Additional entries
  !   ---------------------------------------------------
  ierr=set_harvest_payload_raw(harvest_sendline,trim(tglf_harvest_extra_in)//NUL)

  !   ---------------------------------------------------
  !    Send data
  !   ---------------------------------------------------
  ierr=harvest_send(harvest_sendline)

end subroutine tglf_harvest_local
