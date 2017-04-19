subroutine vgen_compute_neo(i,vtor_diff, rotation_model, er0, &
     omega, omega_deriv, simntheta)

  use vgen_globals
  use neo_interface
  use EXPRO_interface
  use mpi

  implicit none

  integer, intent(in) :: i
  integer, intent(in) :: rotation_model
  real, intent(in)    :: er0            ! kV
  real, intent(in)    :: omega          ! 1/s
  real, intent(in)    :: omega_deriv    ! 1/(m s)
  real, intent(out)   :: vtor_diff      ! vtor_exp - vtor_neo (m/s) 
  integer, intent(out) :: simntheta

  integer :: j, n, is
  real :: cc, loglam

  integer :: nmin, nmax, nth
  real :: cpu_in, cpu_out

  ! Set the local NEO input parameters
  neo_silent_flag_in = 1

  ! Normalizations
  
  dens_norm = EXPRO_ne(i)
  temp_norm = EXPRO_te(i)
  mass_norm = neo_mass_in(1)
  vth_norm  = sqrt(temp_norm * temp_norm_fac &
       / (mass_norm * mass_deuterium)) &
       * 1.0e4 / EXPRO_rmin(EXPRO_n_exp)

  cc = sqrt(2.0) * pi * charge_norm_fac**4 &
       * 1.0 / (4.0 * pi * 8.8542)**2 &
       * 1.0 / (sqrt(mass_deuterium) * temp_norm_fac**1.5) &
       * 1e9
  
  ! Rotation paramters
  
  neo_rotation_model_in  = rotation_model
  neo_dphi0dr_in         = -er0 / EXPRO_grad_r0(i) &
       * EXPRO_rmin(EXPRO_n_exp) / temp_norm
  neo_omega_rot_in       = omega / vth_norm 
  neo_omega_rot_deriv_in = omega_deriv / vth_norm * EXPRO_rmin(EXPRO_n_exp) 
  
  ! Geometry parameters
  
  neo_rmin_over_a_in = EXPRO_rmin(i)/EXPRO_rmin(EXPRO_n_exp)
  neo_rmaj_over_a_in = EXPRO_rmaj(i)/EXPRO_rmin(EXPRO_n_exp)
  neo_q_in           = abs(EXPRO_q(i))
  neo_shear_in       = EXPRO_s(i)
  neo_shift_in       = EXPRO_drmaj(i)
  neo_kappa_in       = EXPRO_kappa(i)
  neo_s_kappa_in     = EXPRO_skappa(i)
  neo_delta_in       = EXPRO_delta(i)
  neo_s_delta_in     = EXPRO_sdelta(i)
  neo_zeta_in        = EXPRO_zeta(i) 
  neo_s_zeta_in      = EXPRO_szeta(i) 
  neo_zmag_over_a_in = EXPRO_zmag(i)/EXPRO_rmin(EXPRO_n_exp)
  neo_s_zmag_in      = EXPRO_dzmag(i) 
  
  neo_rho_star_in = sqrt(temp_norm * temp_norm_fac &
       * mass_norm * mass_deuterium) &
       / (charge_norm_fac * abs(EXPRO_bunit(i))) &
       * 1.0e-4 / EXPRO_rmin(EXPRO_n_exp)
  
  if(neo_equilibrium_model_in == 3) then
     neo_geo_ny_in = EXPRO_nfourier
     neo_geo_yin_in(:,:) = 0.0
     do n=0,EXPRO_nfourier
        do j=1,4  
           neo_geo_yin_in(j,n)   = EXPRO_geo(j,n,i)/EXPRO_rmin(EXPRO_n_exp)
           neo_geo_yin_in(j+4,n) = EXPRO_dgeo(j,n,i)
        enddo
     enddo
  endif
  
  ! Species-dependent parameters
  
  neo_ne_ade_in = EXPRO_ne(i) / dens_norm
  neo_te_ade_in = EXPRO_te(i) / temp_norm
  
  loglam = 24.0 - log(sqrt(EXPRO_ne(i)*1e13)/(EXPRO_te(i)*1000))
  
  ! Species 1 -- must be an ion species
  neo_dens_in(1)   = EXPRO_ni_new(i) / dens_norm
  neo_dlnndr_in(1) = EXPRO_dlnnidr_new(i) * EXPRO_rmin(EXPRO_n_exp)
  neo_temp_in(1)   = EXPRO_ti(1,i) / temp_norm
  neo_dlntdr_in(1) = EXPRO_dlntidr(1,i) * EXPRO_rmin(EXPRO_n_exp)
  neo_nu_1_in      = cc * loglam * neo_dens_in(1) * dens_norm &
       * neo_z_in(1)**4 / (sqrt(neo_mass_in(1)) &
       * (neo_temp_in(1) * temp_norm)**1.5) / vth_norm
  
  ! Species 2
  do is=2,neo_n_species_in
     if(neo_z_in(is) /= -1) then
        neo_dens_in(is)   = EXPRO_ni(is,i) / dens_norm
        neo_dlnndr_in(is) = EXPRO_dlnnidr(is,i) * EXPRO_rmin(EXPRO_n_exp)
        neo_temp_in(is)   = EXPRO_ti(is,i) / temp_norm 
        neo_dlntdr_in(is) = EXPRO_dlntidr(is,i) * EXPRO_rmin(EXPRO_n_exp)
     else
        neo_dens_in(is)   = EXPRO_ne(i) / dens_norm
        neo_temp_in(is)   = EXPRO_te(i) / temp_norm
        neo_dlnndr_in(is) = EXPRO_dlnnedr(i) * EXPRO_rmin(EXPRO_n_exp)
        neo_dlntdr_in(is) = EXPRO_dlntedr(i) * EXPRO_rmin(EXPRO_n_exp)
     endif
  enddo

  nmin = (nth_min - 1)/2
  nmax = (nth_max - 1)/2
  nth = nint(nmin*sqrt(EXPRO_thetascale(i)/EXPRO_thetascale(2)))
  if(nth < nmax) then
     if(nth < nmin) then
        neo_n_theta_in = 2*nmin + 1
     else
        neo_n_theta_in = 2*nth + 1
     endif
  else
     neo_n_theta_in = 2*nmax + 1
  endif

  if (i == 2+i_proc) then
     open(unit=1,file='out.vgen.neontheta'//tag(i_proc+1),status='replace')
     close(1)
  endif
  open(unit=1,file='out.vgen.neontheta'//tag(i_proc+1),status='old',position='append')
  write(1,'(e16.8)',advance='no') EXPRO_rho(i)
  write(1,'(i3)',advance='no') neo_n_theta_in
  write(1,*)
  close(1)
  simntheta=neo_n_theta_in

  if(i == 2+i_proc) then
     open(unit=1,file='out.vgen.neoexpnorm'//tag(i_proc+1),status='replace')
     close(1)
  endif
  open(unit=1,file='out.vgen.neoexpnorm'//tag(i_proc+1),status='old',position='append')
  write(1,'(e16.8)',advance='no') EXPRO_rho(i)
  write(1,'(e16.8)',advance='no') EXPRO_rmin(EXPRO_n_exp)
  write(1,'(e16.8)',advance='no') mass_norm
  write(1,'(e16.8)',advance='no') dens_norm
  write(1,'(e16.8)',advance='no') temp_norm
  write(1,'(e16.8)',advance='no') vth_norm*EXPRO_rmin(EXPRO_n_exp)
  write(1,'(e16.8)',advance='no') EXPRO_bunit(i)
  write(1,*)
  close(1)

  if(i == 2+i_proc) then
     open(unit=1,file='out.vgen.neoequil'//tag(i_proc+1),status='replace')
     close(1)
  endif
  open(unit=1,file='out.vgen.neoequil'//tag(i_proc+1),status='old',position='append')
  write(1,'(e16.8)',advance='no') EXPRO_rho(i)
  write(1,'(e16.8)',advance='no') neo_rmin_over_a_in
  write(1,'(e16.8)',advance='no') neo_q_in
  write(1,'(e16.8)',advance='no') neo_rho_star_in
  write(1,'(e16.8)',advance='no') neo_rmaj_over_a_in
  write(1,'(e16.8)',advance='no') neo_nu_1_in
  do is=1,neo_n_species_in
     write(1,'(e16.8)',advance='no') neo_dens_in(is)
     write(1,'(e16.8)',advance='no') neo_temp_in(is)
     write(1,'(e16.8)',advance='no') neo_dlnndr_in(is)
     write(1,'(e16.8)',advance='no') neo_dlntdr_in(is)
  enddo
  close(1)

  ! Run NEO
  cpu_in = MPI_Wtime()
  call neo_run()
  cpu_out = MPI_Wtime()

  if(timing_flag == 1) then
     if (i == 2+i_proc) then
        open(unit=1,file='out.vgen.neotime'//tag(i_proc+1),status='replace')
        close(1)
     endif
     open(unit=1,file='out.vgen.neotime'//tag(i_proc+1),status='old',position='append')
     write(1,'(e16.8)',advance='no') EXPRO_rho(i)
     write(1,'(e16.8)',advance='no') cpu_out-cpu_in
     write(1,*)
     close(1)
  endif

  if (neo_error_status_out > 0) then
     print *,neo_error_message_out
     stop
  endif

  ! Assign output flows and bootstrap current

  vtor_diff = vtor_measured(i) &
       - neo_vtor_dke_out(erspecies_indx) * vth_norm * EXPRO_rmin(EXPRO_n_exp)
  
  do j=1,n_ions
     EXPRO_vpol(j,i) = neo_vpol_dke_out(j) &
          * vth_norm * EXPRO_rmin(EXPRO_n_exp)
     EXPRO_vtor(j,i) = neo_vtor_dke_out(j) &
          * vth_norm * EXPRO_rmin(EXPRO_n_exp)
  enddo
  jbs_norm = charge_norm_fac*dens_norm*vth_norm &
       *EXPRO_rmin(EXPRO_n_exp)/1e6
  jbs_neo(i)    = neo_jpar_dke_out*jbs_norm
  jbs_sauter(i) = neo_jpar_thS_out*jbs_norm
  jbs_koh(i)    = neo_jpar_thK_out*jbs_norm
  jbs_nclass(i) = neo_jpar_thN_out*jbs_norm
  pflux_sum(i)  = 0.0
  do j=1,neo_n_species_in
     pflux_sum(i) = pflux_sum(i) + zfac(j)*neo_pflux_dke_out(j)
  enddo
  pflux_sum(i) = pflux_sum(i) / neo_rho_star_in**2
  
end subroutine vgen_compute_neo
