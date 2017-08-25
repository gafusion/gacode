  ! Harvest LOCAL INTERFACE variables
  !
  ! _bol : boolean
  ! _int : integer
  ! _dbl : double

  SUBROUTINE tglf_harvest_local()

    USE tglf_interface

    INTEGER :: ierr, i, j, nky, get_nky_out
    CHARACTER(LEN=65507) :: harvest_sendline
    CHARACTER(LEN=2) :: NUM
    CHARACTER NUL
    PARAMETER(NUL = CHAR(0))

    REAL, DIMENSION(10) :: tmp
    INTEGER, DIMENSION(10) :: ions_order
    
    !REAL, DIMENSION(:), ALLOCATABLE :: spectrum

    EXTERNAL :: get_nky_out

    IF (.NOT.tglf_use_transport_model_in) THEN
        WRITE(1,*) 'HARVEST ONLY WHEN `TRANSPORT_MODEL=.TRUE.`' !only when computing fluxes
        RETURN
    ENDIF

    IF (.NOT.tglf_iflux_in) THEN
        WRITE(1,*) 'HARVEST ONLY AVAILABLE WHEN `IFLUX=.TRUE.' !only when computing fluxes
        RETURN
    ENDIF

    IF (tglf_geometry_flag_in .NE. 1 ) THEN
       WRITE(1,*) 'HARVEST ONLY SUPPORTS MILLER `GEOMETRY_FLAG=1`' !only when using miller
       RETURN
    ENDIF

    ierr=init_harvest('TGYRO_TGLF_NEO'//NUL,harvest_sendline,LEN(harvest_sendline))
    ierr=set_harvest_protocol('UDP'//NUL)

    ierr=set_harvest_payload_str(harvest_sendline,'VERSION'//NUL,'APS15_1'//NUL) !no underscore to allow different versions of the same run

!   ---------------------------------------------------
!    Plasma parameters
!   ---------------------------------------------------
    ierr=set_harvest_payload_int(harvest_sendline,'SIGN_BT'//NUL,tglf_sign_bt_in)
    ierr=set_harvest_payload_int(harvest_sendline,'SIGN_IT'//NUL,tglf_sign_it_in)
    IF (tglf_kygrid_model_in.NE.1) THEN
        ierr=set_harvest_payload_dbl(harvest_sendline,'KY'//NUL,tglf_ky_in)
    ENDIF
    ierr=set_harvest_payload_dbl(harvest_sendline,'VEXB_SHEAR'//NUL,tglf_vexb_shear_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'BETAE'//NUL,tglf_betae_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'XNUE'//NUL,tglf_xnue_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'ZEFF'//NUL,tglf_zeff_in)
    IF (tglf_debye_in.NE.0) THEN
        ierr=set_harvest_payload_dbl(harvest_sendline,'DEBYE'//NUL,tglf_debye_in)
    ENDIF

!   ---------------------------------------------------
!    Sort ions by A, Z, Te/Ti
!    electrons always first (as required by TGLF)
!    !main ion always second (affects GB normalization)
!   ---------------------------------------------------
    tmp=0.0
    tmp(1)=1E10
    tmp(2)=1E9
    ions_order=0
    DO i = 3,tglf_ns_in
        tmp(i)=tglf_mass_in(i)*100000+tglf_zs_in(i)*1000+1./(1+tglf_rlts_in(i))
    ENDDO
    DO i = 1,tglf_ns_in
        j=MAXLOC(tmp, DIM=1)
        ions_order(i)=j
        tmp(j)=0
    ENDDO

!   ---------------------------------------------------
!    Species vectors
!   ---------------------------------------------------
    DO i = 1,tglf_ns_in
      IF (i < 10) THEN
         write (NUM, "(I01,A1)") i,NUL
      ELSE
         write (NUM, "(I02,A1)") i,NUL
      ENDIF

      j=ions_order(i)

      ierr=set_harvest_payload_dbl(harvest_sendline,'+ZS_'//NUM,tglf_zs_in(j))     !make new tables for different charges
      ierr=set_harvest_payload_dbl(harvest_sendline,'+MASS_'//NUM,tglf_mass_in(j)) !make new tables for different masses
      ierr=set_harvest_payload_dbl(harvest_sendline,'RLNS_'//NUM,tglf_rlns_in(j))
      ierr=set_harvest_payload_dbl(harvest_sendline,'RLTS_'//NUM,tglf_rlts_in(j))
      if (i.NE.1) THEN !these are always 1, by definition
        ierr=set_harvest_payload_dbl(harvest_sendline,'TAUS_'//NUM,tglf_taus_in(j))
        ierr=set_harvest_payload_dbl(harvest_sendline,'AS_'//NUM,tglf_as_in(j))
      ENDIF
      ierr=set_harvest_payload_dbl(harvest_sendline,'VPAR_'//NUM,tglf_vpar_in(j))
      ierr=set_harvest_payload_dbl(harvest_sendline,'VPAR_SHEAR_'//NUM,tglf_vpar_shear_in(j))
    ENDDO

!   ---------------------------------------------------
!    Miller geometry parameters
!   ---------------------------------------------------
    ierr=set_harvest_payload_dbl(harvest_sendline,'RMIN_LOC'//NUL,tglf_rmin_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'RMAJ_LOC'//NUL,tglf_rmaj_loc_in)
    ierr=set_harvest_payload_dbl(harvest_sendline,'ZMAJ_LOC'//NUL,tglf_zmaj_loc_in)
    IF (tglf_drmindx_loc_in.NE.1) THEN
        ierr=set_harvest_payload_dbl(harvest_sendline,'DRMINDX_LOC'//NUL,tglf_drmindx_loc_in) ! different derivatives affect everything!
    ENDIF
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
    IF (tglf_kx0_loc_in.NE.0) THEN
        ierr=set_harvest_payload_dbl(harvest_sendline,'KX0_LOC'//NUL,tglf_kx0_loc_in) !should be always 0.0
    ENDIF

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
    DO i = 1,tglf_ns_in
      IF (i < 10) THEN
         write (NUM, "(I01,A1)") i,NUL
      ELSE
         write (NUM, "(I02,A1)") i,NUL
      ENDIF

      j=ions_order(i)

      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_PARTICLE_FLUX_ESTATIC_'//NUM,get_particle_flux(j,1))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_ENERGY_FLUX_ESTATIC_'//NUM,get_energy_flux(j,1))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_PAR_ESTATIC_'//NUM,get_stress_par(j,1))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_TOR_ESTATIC_'//NUM,get_stress_tor(j,1))
      ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_EXCHANGE_ESTATIC_'//NUM,get_exchange(j,1))

      if (tglf_use_bper_in) THEN
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_PARTICLE_FLUX_EMPER_'//NUM,get_particle_flux(j,2))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_ENERGY_FLUX_EMPER_'//NUM,get_energy_flux(j,2))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_PAR_EMPER_'//NUM,get_stress_par(j,2))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_TOR_EMPER_'//NUM,get_stress_tor(j,2))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_EXCHANGE_EMPER_'//NUM,get_exchange(j,2))
      ENDIF

      if (tglf_use_bpar_in) THEN
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_PARTICLE_FLUX_EMPAR_'//NUM,get_particle_flux(j,3))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_ENERGY_FLUX_EMPAR_'//NUM,get_energy_flux(j,3))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_PAR_EMPAR_'//NUM,get_stress_par(j,3))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_TOR_EMPAR_'//NUM,get_stress_tor(j,3))
          ierr=set_harvest_payload_dbl(harvest_sendline,'OUT_STRESS_EXCHANGE_EMPAR_'//NUM,get_exchange(j,3))
      ENDIF

   ENDDO

! !   ---------------------------------------------------
! !    Spectra (to gather the spectra this function must be in TGLF itself)
! !   ---------------------------------------------------
!    nky = get_nky_out()

!    ALLOCATE(spectrum(nky))
   
!    DO i = 1, nky
!       spectrum(i) = get_ky_spectrum_out(i)
!    ENDDO

!    ierr=set_harvest_payload_dbl_array(harvest_sendline,'KY_SPECTRUM'//NUL,spectrum,nky)
   
!    DO i = 1, tglf_nmodes_in
!       IF (i < 10) THEN
!          write (NUM, "(I01,A1)") i,NUL
!       ELSE
!          write (NUM, "(I02,A1)") i,NUL
!       ENDIF
      
!       DO j = 1, nky
!          spectrum(j) = get_eigenvalue_spectrum_out(1,j,i)
!       ENDDO 
!       ierr=set_harvest_payload_dbl_array(harvest_sendline,'OUT_EIGENVALUE_SPECTRUM_GAMMA'//NUM,spectrum,nky)
!       DO j = 1, nky
!          spectrum(j) = get_eigenvalue_spectrum_out(2,j,i)
!       ENDDO
!       ierr=set_harvest_payload_dbl_array(harvest_sendline,'OUT_EIGENVALUE_SPECTRUM_OMEGA'//NUM,spectrum,nky)
!    ENDDO
!    DEALLOCATE(spectrum)

!   ---------------------------------------------------
!   Additional entries
!   ---------------------------------------------------
    ierr=set_harvest_payload_raw(harvest_sendline,TRIM(tglf_harvest_extra_in)//NUL)

!   ---------------------------------------------------
!    Send data
!   ---------------------------------------------------
   ierr=harvest_send(harvest_sendline)

  END SUBROUTINE tglf_harvest_local
