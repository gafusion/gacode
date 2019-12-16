!---------------------------------------------------------
! tglf_run_mpi.f90
!
! PURPOSE:
!  Manage call to TGLF simulation for both standalone and
!  TGYRO usage.
!---------------------------------------------------------

#ifdef MPI_TGLF
subroutine tglf_run_mpi()
#else
subroutine tglf_run()
#endif

  use tglf_pkg
  use tglf_interface
  use tglf_global
  
#ifdef MPI_TGLF
  use tglf_mpi
#endif

  implicit none

  integer :: i_ion,n

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
       tglf_nky_in)
       
  call put_rare_switches(tglf_theta_trapped_in, &
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
       tglf_sat_rule_in, &
       tglf_kygrid_model_in, &
       tglf_xnu_model_in, &
       tglf_vpar_model_in, &
       tglf_vpar_shear_model_in, &
       tglf_write_wavefunction_flag_in)

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

  ! Create parameter dump file
  if (tglf_dump_flag_in .eqv. .true.) then
     call tglf_dump_local
     if (tglf_test_flag_in == 1) return
  endif

  if (tglf_use_transport_model_in) then

     ! Call TGLFNN and later TGLF if TGLFNN was not accurate
     if (tglf_nn_max_error_in .gt. 0) then
        call tglf_nn_tm
     endif
     if (.not. valid_nn) then
#ifdef MPI_TGLF
        call tglf_tm_mpi
#else
        call tglf_tm
#endif
     endif
     
     !---------------------------------------------
     ! Output (normalized to Q_GB)
     ! 
     ! Electrons

     ! Gammae/Gamma_GB
     tglf_elec_pflux_out = get_particle_flux(1,1)  &
               + get_particle_flux(1,2)            &
               + get_particle_flux(1,3)

     ! Qe/Q_GB
     tglf_elec_eflux_low_out = get_q_low(1)
     tglf_elec_eflux_out     = get_energy_flux(1,1) &
               + get_energy_flux(1,2)               &
               + get_energy_flux(1,3)

     ! Pi_e/Pi_GB
     tglf_elec_mflux_out = get_stress_tor(1,1)   &
               + get_stress_tor(1,2)             &
               + get_stress_tor(1,3)

     ! S_e/S_GB
     tglf_elec_expwd_out = get_exchange(1,1)  &
               + get_exchange(1,2)            &
               + get_exchange(1,3)

     ! Ions

     do i_ion=1,ns_in

             ! Gammai/Gamma_GB
             tglf_ion_pflux_out(i_ion) = get_particle_flux(i_ion+1,1) &
                  + get_particle_flux(i_ion+1,2)                      &
                  + get_particle_flux(i_ion+1,3)

             ! Qi/Q_GB
             tglf_ion_eflux_low_out(i_ion) = get_q_low(i_ion+1)
             ! _low already included EM terms
             tglf_ion_eflux_out(i_ion)     = get_energy_flux(i_ion+1,1) &
                  + get_energy_flux(i_ion+1,2)                          &
                  + get_energy_flux(i_ion+1,3)

             ! Pi_i/Pi_GB
             tglf_ion_mflux_out(i_ion) = get_stress_tor(i_ion+1,1)  &
                  + get_stress_tor(i_ion+1,2)                       &
                  + get_stress_tor(i_ion+1,3)

             ! S_i/S_GB
             tglf_ion_expwd_out(i_ion) = get_exchange(i_ion+1,1)    &
                  + get_exchange(i_ion+1,2)                         &
                  + get_exchange(i_ion+1,3)

     enddo

     call get_error_status(tglf_error_message,tglf_error_status)

  else

     ! Run single-ky linear stability
     call tglf_ky

     ! Collect linear eigenvalues
     do n=1,tglf_nmodes_in
        tglf_eigenvalue_out(n) = get_frequency(n) + xi*get_growthrate(n)
     enddo

     ! Print eigenfunction if flag set
#ifdef MPI_TGLF
     if (iProcTglf == iProc0Tglf .and. tglf_write_wavefunction_flag_in == 1) then
        call write_wavefunction_out(trim(tglf_path_in)//'out.tglf.wavefunction')
     endif     
#else
     if (tglf_write_wavefunction_flag_in == 1) then
        call write_wavefunction_out(trim(tglf_path_in)//'out.tglf.wavefunction')
     endif
#endif

  endif

  ! Diagnostic output
  interchange_DR = get_DR()
  interchange_DM = get_DM()

#ifdef MPI_TGLF
end subroutine tglf_run_mpi
#else
end subroutine tglf_run
#endif
