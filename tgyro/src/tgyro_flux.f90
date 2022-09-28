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
  use neo_interface
  use tglf_interface

  implicit none

  integer :: i_ion,i0
  real, dimension(8) :: x_out
  real :: a1,a2,a3,a4

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

  select case (tgyro_neo_method) 

  case (0)

     ! Zero neoclassical flux

     pflux_e_neo(i_r) = 0.0
     eflux_e_neo(i_r) = 0.0

     pflux_i_neo(:,i_r) = 0.0
     eflux_i_neo(:,i_r) = 0.0

  case (1)

     ! NEO analytic theory: Hirshman-Sigmar

     call neo_run

     pflux_e_neo(i_r)   = neo_pflux_thHS_out(1)/Gamma_neo_GB 
     eflux_e_neo(i_r)   = neo_eflux_thHS_out(1)/Q_neo_GB

     i0 = 1
     do i_ion=1,loc_n_ion
        if (calc_flag(i_ion) == 0) cycle
        i0 = i0+1 
        pflux_i_neo(i_ion,i_r) = neo_pflux_thHS_out(i0)/Gamma_neo_GB
        eflux_i_neo(i_ion,i_r) = neo_eflux_thHS_out(i0)/Q_neo_GB
     enddo

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
        if (calc_flag(i_ion) == 0) cycle
        i0 = i0+1 
        pflux_i_neo(i_ion,i_r) = (neo_pflux_dke_out(i0) + &
             tgyro_neo_gv_flag*neo_pflux_gv_out(i0))/Gamma_neo_GB 
        eflux_i_neo(i_ion,i_r) = (neo_efluxncv_dke_out(i0) + &
             tgyro_neo_gv_flag*neo_efluxncv_gv_out(i0))/Q_neo_GB
        mflux_i_neo(i_ion,i_r) = (neo_mflux_dke_out(i0) + &
             tgyro_neo_gv_flag*neo_mflux_gv_out(i0))/Pi_neo_GB
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
     eflux_e_tur(i_r) = tglf_elec_eflux_out
     mflux_e_tur(i_r) = -tglf_sign_It_in*tglf_elec_mflux_out
     expwd_e_tur(i_r) = tglf_elec_expwd_out

     i0 = 0
     do i_ion=1,loc_n_ion
        if (calc_flag(i_ion) == 0) cycle
        i0 = i0+1 
        pflux_i_tur(i_ion,i_r) = tglf_ion_pflux_out(i0)
        eflux_i_tur(i_ion,i_r) = tglf_ion_eflux_out(i0)
        mflux_i_tur(i_ion,i_r) = -tglf_sign_It_in*tglf_ion_mflux_out(i0)
        expwd_i_tur(i_ion,i_r) = tglf_ion_expwd_out(i0)
     enddo

  case (5) 

     call tgyro_etgcrit(a1,a2,a3,a4)
     eflux_e_tur(i_r)   = a1
     eflux_i_tur(1,i_r) = a2
     pflux_e_tur(i_r)   = a3
     pflux_i_tur(1,i_r) = a4

     call tgyro_trap_component_error(0,'null')

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
