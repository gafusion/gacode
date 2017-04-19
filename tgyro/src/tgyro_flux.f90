!------------------------------------------------------------
! tgyro_flux.f90
!
! PURPOSE:
!  Manage calls to obtain neoclassical and turbulent fluxes.
!  
! NOTE:
!  Component errors/warnings are trapped and managed here.
!------------------------------------------------------------

subroutine tgyro_flux

  use mpi
  use tgyro_globals
  use gyro_interface
  use neo_interface
  use tglf_interface
  use glf23_interface

  implicit none

  integer :: i_ion,i0, is, is_ion
  integer :: i1,i2
  integer :: n_12
  real :: d_max
  real, dimension(8) :: x_out

  !-------------------------------------------
  ! IFS-PPPL parameters
  real :: rltcrit
  real :: rltcritz
  real :: chii0
  real :: chie0
  real :: dummy1
  real :: dummy2
  !-------------------------------------------
  !-------------------------------------------
  ! NEO parameters
  real :: Gamma_neo_GB
  real :: Q_neo_GB
  real :: Pi_neo_GB
  !-------------------------------------------

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),position='append')
     close(1)
  endif

  !-----------------------------
  ! Counter
  flux_counter = flux_counter+1
  !-----------------------------

  !-----------------------------------------------------------
  ! Neoclassical fluxes
  !-----------------------------------------------------------

  pflux_i_neo(:,:) = 0.0
  pflux_e_neo(:)   = 0.0
  eflux_i_neo(:,:) = 0.0
  eflux_e_neo(:)   = 0.0
  mflux_i_neo(:,:) = 0.0
  mflux_e_neo(:)   = 0.0

  ! Map TGYRO parameters to NEO
  call tgyro_neo_map  

  ! NEO normalization to GB output conversions
  ! Gamma_norm = n0_norm vt_norm
  ! Q_norm     = n0_norm vt_norm T_norm
  ! Pi_norm    = n0_norm T0_norm a_norm
  ! where
  ! vt_norm = sqrt(T_norm / mass_norm)
  !
  ! Gamma_GB = n_0e*c_s*(rho_s/a)^2
  ! Q_GB     = n_0e*c_s*(rho_s/a)^2*T_0e
  ! Pi_GB    = n_0e*T_0e*a*(rho_s/a)^2
  ! where
  ! c_s    = sqrt(T_0e/m_i)
  ! rho_se = c_s/Omega_ci
  Gamma_neo_GB = neo_dens_in(1) * neo_temp_in(1)**1.5 * neo_rho_star_in**2
  Q_neo_GB     = neo_dens_in(1) * neo_temp_in(1)**2.5 * neo_rho_star_in**2
  Pi_neo_GB    = neo_dens_in(1) * neo_temp_in(1)**2   * neo_rho_star_in**2

  select case (loc_neo_method) 

  case (1)

     ! NEO analytic theory: Hinton-Hazeltine
     !  
     ! NOTE: impurities have zero flux in this case

     call neo_run

     pflux_e_neo(i_r)   = neo_pflux_thHH_out /Gamma_neo_GB 
     eflux_e_neo(i_r)   = neo_eflux_thHHe_out /Q_neo_GB

     is_ion = -1
     d_max = -1.0
     do is=2, neo_n_species_in
        if(neo_z_in(is) > 0 .and. neo_dens_in(is) > d_max) then
           is_ion = is
           d_max  = neo_dens_in(is)
        endif
     enddo
     i0 = 1
     do i_ion=1,loc_n_ion
        if (therm_flag(i_ion) == 0 .and. tgyro_quickfast_flag == 1) cycle
        i0 = i0+1
        if(i0 == is_ion) exit
     enddo
     
     pflux_i_neo(i_ion,i_r) = neo_pflux_thHH_out/Gamma_neo_GB
     if (loc_chang_hinton == 1) then
        eflux_i_neo(i_ion,i_r) = neo_eflux_thCHi_out/Q_neo_GB
     else
        eflux_i_neo(i_ion,i_r) = neo_eflux_thHHi_out/Q_neo_GB
     endif
     
  case (2)

     ! Call NEO kinetic calculation
     if (gyrotest_flag == 0) call neo_run

     pflux_e_neo(i_r) = (neo_pflux_dke_out(1)+tgyro_neo_gv_flag*neo_pflux_gv_out(1)) &
          /Gamma_neo_GB 
     eflux_e_neo(i_r) = (neo_efluxncv_dke_out(1)+tgyro_neo_gv_flag*neo_efluxncv_gv_out(1)) &
          /Q_neo_GB
     mflux_e_neo(i_r) = (neo_mflux_dke_out(1)+tgyro_neo_gv_flag*neo_mflux_gv_out(1)) &
          /Pi_neo_GB

     i0 = 1
     do i_ion=1,loc_n_ion
        if (therm_flag(i_ion) == 0 .and. tgyro_quickfast_flag == 1) cycle
        i0 = i0+1 
        pflux_i_neo(i_ion,i_r) = (neo_pflux_dke_out(i0) + &
             tgyro_neo_gv_flag*neo_pflux_gv_out(i0))/Gamma_neo_GB 
        eflux_i_neo(i_ion,i_r) = (neo_efluxncv_dke_out(i0) + &
             tgyro_neo_gv_flag*neo_efluxncv_gv_out(i0))/Q_neo_GB
        mflux_i_neo(i_ion,i_r) = (neo_mflux_dke_out(i0) + &
             tgyro_neo_gv_flag*neo_mflux_gv_out(i0))/Pi_neo_GB
     enddo

  case (3)

     ! NEO analytic theory: Hirshman-Sigmar

     call neo_run

     pflux_e_neo(i_r)   = neo_pflux_thHS_out(1)/Gamma_neo_GB 
     eflux_e_neo(i_r)   = neo_eflux_thHS_out(1)/Q_neo_GB

     i0 = 1
     do i_ion=1,loc_n_ion
        if (therm_flag(i_ion) == 0 .and. tgyro_quickfast_flag == 1) cycle
        i0 = i0+1 
        pflux_i_neo(i_ion,i_r) = neo_pflux_thHS_out(i0)/Gamma_neo_GB
        eflux_i_neo(i_ion,i_r) = neo_eflux_thHS_out(i0)/Q_neo_GB
     enddo

  end select

  call tgyro_trap_component_error(neo_error_status_out,neo_error_message_out)

  !-----------------------------------------------------------
  ! Turbulent fluxes
  !-----------------------------------------------------------

  pflux_i_tur(:,:) = 0.0
  pflux_e_tur(:)   = 0.0
  eflux_i_tur(:,:) = 0.0
  eflux_e_tur(:)   = 0.0
  mflux_i_tur(:,:) = 0.0
  mflux_e_tur(:)   = 0.0
  expwd_i_tur(:,:) = 0.0
  expwd_e_tur(:)   = 0.0

  select case (flux_method)

  case (0)

     ! No fluxes (tgyro_noturb_flag=1)

  case (1) 

     call ifs_pppl(r_maj(i_r)*dlntidr(1,i_r),&
          r_maj(i_r)*dlnnidr(1,i_r),&
          r_maj(i_r)*dlnnedr(i_r),&
          abs(q(i_r)),&
          kappa(i_r),&
          s(i_r),&
          1.0,&
          0.0,&
          ti(1,i_r)/te(i_r),&
          r(i_r)/r_maj(i_r),&
          2.5e-7*ne(i_r)/(te(i_r)**1.5*ti(1,i_r)**0.5)*(r_maj(i_r)/100), &
          r_maj(i_r),&
          rho_i(i_r),&
          v_i(i_r), &
          rltcrit,&
          rltcritz,&
          dummy1,&
          dummy2, &
          chii0,&
          chie0)

     ! Normalize to chi_GB = rho_s^2 c_s/a
     x_out(2) = chie0/chi_gb(i_r)
     x_out(4) = chii0/chi_gb(i_r)

     ! Convert to flux
     eflux_e_tur(i_r) = x_out(2)*r_min*dlntedr(i_r)

     eflux_i_tur(1,i_r) = x_out(4)*r_min*dlntidr(1,i_r)*&
          ni(1,i_r)/ne(i_r)*ti(1,i_r)/te(i_r)

     call tgyro_trap_component_error(0,'null')

  case (2)

     ! Map TGYRO parameters to TGLF
     call tgyro_tglf_map

     call tglf_run_mpi

     call tgyro_trap_component_error(tglf_error_status,tglf_error_message)

     pflux_e_tur(i_r) = tglf_elec_pflux_out
     if (tglf_q_low_flag == 1) then
        eflux_e_tur(i_r) = tglf_elec_eflux_low_out
     else
        eflux_e_tur(i_r) = tglf_elec_eflux_out
     endif
     mflux_e_tur(i_r) = -tglf_sign_It_in*tglf_elec_mflux_out
     expwd_e_tur(i_r) = tglf_elec_expwd_out

     i0 = 0
     do i_ion=1,loc_n_ion
        if (therm_flag(i_ion) == 0 .and. tgyro_quickfast_flag == 1) cycle
        i0 = i0+1 
        pflux_i_tur(i_ion,i_r) = tglf_ion_pflux_out(i0)
        if (tglf_q_low_flag == 1) then
           eflux_i_tur(i_ion,i_r) = tglf_ion_eflux_low_out(i0)
        else
           eflux_i_tur(i_ion,i_r) = tglf_ion_eflux_out(i0)
        endif
        mflux_i_tur(i_ion,i_r) = -tglf_sign_It_in*tglf_ion_mflux_out(i0)
        expwd_i_tur(i_ion,i_r) = tglf_ion_expwd_out(i0)
     enddo

  case (3)  ! Map TGYRO parameters to GLF23

     call tgyro_glf23_map

     if (gyrotest_flag == 0) call glf23_run

     !     call tgyro_trap_component_error(glf23_error_status,glf23_error_message)

     pflux_e_tur(i_r) = glf23_elec_pflux_out
     eflux_e_tur(i_r) = glf23_elec_eflux_out
     mflux_e_tur(i_r) = glf23_elec_mflux_out
     expwd_e_tur(i_r) = glf23_elec_expwd_out

     pflux_i_tur(1:loc_n_ion,i_r) = glf23_ion_pflux_out(1:loc_n_ion)
     eflux_i_tur(1:loc_n_ion,i_r) = glf23_ion_eflux_out(1:loc_n_ion)
     mflux_i_tur(1:loc_n_ion,i_r) = glf23_ion_mflux_out(1:loc_n_ion)
     expwd_i_tur(1:loc_n_ion,i_r) = glf23_ion_expwd_out(1:loc_n_ion)

  case (4)

     ! Map TGYRO parameters to GYRO
     call tgyro_gyro_map

     call gyro_run(gyrotest_flag,gyro_restart_method,transport_method)

     call tgyro_trap_component_error(gyro_error_status_out,gyro_error_message_out)

     if (tgyro_mode == 1 .and. gyrotest_flag == 0) then

        i1 = 1+gyro_explicit_damp_grid_in
        i2 = gyro_radial_grid_in-gyro_explicit_damp_grid_in
        n_12 = gyro_radial_grid_in-2*gyro_explicit_damp_grid_in

        ! Map GYRO (local simulation) output to TGYRO
        ! Need to convert from GBD units to GB1 units:

        pflux_e_tur(i_r) = sum(gyro_elec_pflux_out(i1:i2))/n_12*(2*mp/mi(1))**1.5
        eflux_e_tur(i_r) = sum(gyro_elec_eflux_out(i1:i2))/n_12*(2*mp/mi(1))**1.5
        mflux_e_tur(i_r) = sum(gyro_elec_mflux_out(i1:i2))/n_12*(2*mp/mi(1))**2
        expwd_e_tur(i_r) = sum(gyro_elec_expwd_out(i1:i2))/n_12*(2*mp/mi(1))**1.5

        do i_ion=1,loc_n_ion
           pflux_i_tur(i_ion,i_r) = sum(gyro_ion_pflux_out(i1:i2,i_ion))/n_12*(2*mp/mi(1))**1.5
           eflux_i_tur(i_ion,i_r) = sum(gyro_ion_eflux_out(i1:i2,i_ion))/n_12*(2*mp/mi(1))**1.5
           mflux_i_tur(i_ion,i_r) = sum(gyro_ion_mflux_out(i1:i2,i_ion))/n_12*(2*mp/mi(1))**2
           expwd_i_tur(i_ion,i_r) = sum(gyro_ion_expwd_out(i1:i2,i_ion))/n_12*(2*mp/mi(1))**1.5
        enddo

     endif

  case default

     call tgyro_catch_error('ERROR: (TGYRO) No matching flux method in tgyro_flux.')

  end select

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !------------------------------------------------------------------
  ! Compute total fluxes given neoclassical and turbulent components:
  !
  ! NOTE: therm_vec here ensures ion sum is only over thermal species
  !
  ! Particle flux for every species
  pflux_e_tot(i_r) = pflux_e_neo(i_r)+pflux_e_tur(i_r)
  pflux_i_tot(:,i_r) = pflux_i_neo(:,i_r)+pflux_i_tur(:,i_r)

  ! Energy fluxes for electrons and summed thermal ions
  eflux_e_tot(i_r) = eflux_e_neo(i_r)+eflux_e_tur(i_r)
  eflux_i_tot(i_r) = sum(eflux_i_neo(therm_vec(:),i_r)+eflux_i_tur(therm_vec(:),i_r))

  ! Momentum flux
  mflux_tot(i_r) = mflux_e_neo(i_r)+mflux_e_tur(i_r)+&
       sum(mflux_i_neo(therm_vec(:),i_r)+mflux_i_tur(therm_vec(:),i_r))
  !-------------------------------------------------------------------

end subroutine tgyro_flux
