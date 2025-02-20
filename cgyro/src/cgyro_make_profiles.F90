subroutine cgyro_make_profiles

  use cgyro_globals
  use cgyro_io
  use expro_locsim_interface

  implicit none

  integer :: is,ir
  integer :: j
  integer :: num_ele

  real :: cc,loglam

  !-------------------------------------------------------------
  ! Adiabatic-electron logic
  if (n_species == 1) then
     ae_flag = 1
  endif
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Local geometry treatment
  !
  if (equilibrium_model == 1) then
     ! s-alpha
     geo_numeq_flag = -1
     geo_ny = 0      
     allocate(geo_yin(8,0:geo_ny))
     geo_yin(:,:) = 0.0
  else if (equilibrium_model == 2) then
     ! MXH
     geo_numeq_flag = 0
     geo_ny = 0
     allocate(geo_yin(8,0:geo_ny))
     geo_yin(:,:) = 0.0
  else
     ! Fourier
     geo_numeq_flag = 1
     geo_ny = geo_ny_in  
     allocate(geo_yin(8,0:geo_ny))
     do j=0,geo_ny
        geo_yin(:,j) = geo_yin_in(:,j)
     enddo
  endif
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Plasma radial profiles (n,T,etc)
  !
  ! FIELD ORIENTATION NOTES:
  !  Field orientation is accomplished by giving signs to a minimal 
  !  set of quantities:
  !  
  !  1. sign(b_unit)   = -btccw
  !  2. sign(q)        = ipccw*btccw
  !  3. sign(rho_star) = -btccw
  !-----------------------------------------------------------------------

  if (profile_model == 2) then

     ! Experimental profiles

     call expro_locsim_profiles(&
          geo_numeq_flag,&
          udsymmetry_flag,&
          quasineutral_flag,&
          n_species+ae_flag,&
          rmin,&
          btccw,&
          ipccw,&
          a_meters,&
          path,&
          CGYRO_COMM_WORLD)

     do is=1,n_species
        z(is)    = z_loc(is)
        mass(is) = mass_loc(is)/2.0
     enddo

     shift   = shift_loc
     kappa   = kappa_loc
     delta   = delta_loc
     zeta    = zeta_loc
     s_kappa = s_kappa_loc
     s_delta = s_delta_loc
     s_zeta  = s_zeta_loc

     ! MXH (cos's will be reset to 0.0 if udsymmetry_flag=1)
     shape_sin(3)   = shape_sin3_loc
     shape_s_sin(3) = shape_s_sin3_loc
     shape_sin(4)   = shape_sin4_loc
     shape_s_sin(4) = shape_s_sin4_loc
     shape_sin(5)   = shape_sin5_loc
     shape_s_sin(5) = shape_s_sin5_loc
     shape_sin(6)   = shape_sin6_loc
     shape_s_sin(6) = shape_s_sin6_loc
     shape_cos(0)   = shape_cos0_loc
     shape_s_cos(0) = shape_s_cos0_loc
     shape_cos(1)   = shape_cos1_loc
     shape_s_cos(1) = shape_s_cos1_loc
     shape_cos(2)   = shape_cos2_loc
     shape_s_cos(2) = shape_s_cos2_loc
     shape_cos(3)   = shape_cos3_loc
     shape_s_cos(3) = shape_s_cos3_loc
     shape_cos(4)   = shape_cos4_loc
     shape_s_cos(4) = shape_s_cos4_loc
     shape_cos(5)   = shape_cos5_loc
     shape_s_cos(5) = shape_s_cos5_loc
     shape_cos(6)   = shape_cos6_loc
     shape_s_cos(6) = shape_s_cos6_loc
 
     q       = q_loc
     s       = s_loc
     zmag    = zmag_loc
     dzmag   = dzmag_loc
     gamma_e = gamma_e_loc
     gamma_p = gamma_p_loc
     mach    = mach_loc
     rmaj    = rmaj_loc
     rhos    = rhos_loc
     z_eff   = z_eff_loc
     b_unit  = b_unit_loc

     dens(1:n_species)    = dens_loc(1:n_species)     
     temp(1:n_species)    = temp_loc(1:n_species)     
     dlnndr(1:n_species)  = dlnndr_loc(1:n_species)     
     dlntdr(1:n_species)  = dlntdr_loc(1:n_species)
     sdlnndr(1:n_species) = sdlnndr_loc(1:n_species)     
     sdlntdr(1:n_species) = sdlntdr_loc(1:n_species)     
     sbeta                = sbeta_loc
 
     if (ae_flag == 1) then
        is_ele = n_species+1
     else
        is_ele = n_species
     endif

     dens_ele = dens_loc(is_ele)
     temp_ele = temp_loc(is_ele)
     mass_ele = mass(is_ele)
     dlnndr_ele = dlnndr_loc(is_ele)
     dlntdr_ele = dlntdr_loc(is_ele)

     ! Get interpolated geometry coefficients
     if (geo_numeq_flag == 1) then
        geo_ny = geo_ny_loc
        deallocate(geo_yin)
        allocate(geo_yin(8,0:geo_ny))
        geo_yin = geo_yin_loc
     endif

     ! Normalizing quantities
     dens_norm = dens_ele       ! ne e19/m3
     temp_norm = temp_ele       ! Te keV
     mass_norm = mass_deuterium ! mD e-27 kg

     ! Compute vth (m/s) using dimensional quantities.  
     ! mass(i) is thus measured in units of deuterium mass.
     do is=1,n_species
        vth(is) = sqrt(temp(is) * temp_norm_fac / (mass(is) * mass_norm)) &
             * 1.0e4 
     enddo
     vth_norm  = sqrt(temp_norm * temp_norm_fac / mass_norm) * 1.0e4 ! c_s (m/s)

     ! Normalizing rho_star and GB flux factors
     rho_star_norm = sqrt(temp_norm * temp_norm_fac * mass_deuterium) &
          / (charge_norm_fac * b_unit) * 1.0e-4 / a_meters   ! rho_s/a
     gamma_gb_norm = dens_norm * vth_norm * rho_star_norm**2 ! e19 m-2 s-1
     q_gb_norm     = gamma_gb_norm * temp_norm * temp_norm_fac / 1.0e6 ! MW/m2
     pi_gb_norm    = dens_norm * temp_norm * temp_norm_fac * a_meters &
          * rho_star_norm**2 ! N/m

     ! Compute collision frequency

     cc = sqrt(2.0) * pi * charge_norm_fac**4 &
          * 1.0 / (4.0 * pi * 8.8542)**2 &
          * 1.0 / (sqrt(mass_deuterium) * temp_norm_fac**1.5) &
          * 1e9

     loglam = 24.0 - log(sqrt(dens_ele*1e13)/(temp_ele*1e3))
     nu_ee  = cc * loglam * dens_ele / (sqrt(mass_ele)*temp_ele**1.5) &
          / (vth_norm/a_meters) 
     
     ! Electron beta
     betae_unit = betae_loc
     
     ! Debye length (from NRL plasma formulary):
     ! Use input lambda_debye as scaling parameter

     lambda_star = 7.434*sqrt((1e3*temp_norm)/(1e13*dens_norm))/rhos

     ! Normalize
     do is=1,n_species
        dens(is) = dens(is)/dens_norm
        temp(is) = temp(is)/temp_norm
        vth(is)  = vth(is)/vth_norm
     enddo
     dens_ae   = dens_ae/dens_norm
     temp_ae   = temp_ae/temp_norm
     dens_ele  = dens_ele/dens_norm
     temp_ele  = temp_ele/temp_norm
     gamma_e   = gamma_e/(vth_norm/a_meters)
     gamma_p   = gamma_p/(vth_norm/a_meters)
     mach      = mach/vth_norm

     do is=1,n_species
        nu(is) = nu_ee *z(is)**4 &
             * dens(is) / dens_ele &
             * sqrt(mass_ele/mass(is)) * (temp_ele/temp(is))**1.5
     enddo

     ! Re-scaling
     lambda_star      = lambda_star * lambda_star_scale
     gamma_e          = gamma_e     * gamma_e_scale
     gamma_p          = gamma_p     * gamma_p_scale
     mach             = mach        * mach_scale    
     betae_unit       = betae_unit  * betae_unit_scale
     do is=1,n_species
        dlnndr(is) = dlnndr(is) * dlnndr_scale(is) 
        dlntdr(is) = dlntdr(is) * dlntdr_scale(is)  
        nu(is)     = nu(is)     * nu_ee_scale
     enddo

     ! Set beta_* consistent with re-scaled beta and gradients and then re-scale
     ! note: beta_star(0) will be over-written with sonic rotation
     call set_betastar
     beta_star(0)     = beta_star(0)  * beta_star_scale
     beta_star_fac    = beta_star_fac * beta_star_scale
     
  else

     a_meters      = 0.0
     b_unit        = 0.0
     dens_norm     = 0.0
     temp_norm     = 0.0
     vth_norm      = 0.0
     mass_norm     = 0.0
     rho_star_norm = 0.0
     gamma_gb_norm = 0.0
     q_gb_norm     = 0.0
     pi_gb_norm    = 0.0

     ! Set sign of q (assumed positive in input)
     q = abs(q)*(ipccw)*(btccw)

     if (ae_flag == 1) then
        is_ele = -1
        dens_ele = dens_ae
        temp_ele = temp_ae
        mass_ele = mass_ae
        dlnndr_ele = dlnndr_ae
        dlntdr_ele = dlntdr_ae
     else
        is_ele = -1
        do is=1,n_species
           if (z(is) < 0.0) then
              is_ele = is
              exit
           endif
        enddo
        if (is_ele == -1) then
           call cgyro_error('No electron species specified')
           return
        endif
        dens_ele = dens(is_ele)
        temp_ele = temp(is_ele)
        mass_ele = mass(is_ele)
        dlnndr_ele = dlnndr(is_ele)
        dlntdr_ele = dlntdr(is_ele)
     endif

     do is=1,n_species

        ! thermal velocity
        vth(is) = sqrt(temp(is)/mass(is))

        ! collision frequency
        nu(is) = nu_ee *z(is)**4 &
             * dens(is) / dens_ele &
             * sqrt(mass_ele/mass(is)) * (temp_ele/temp(is))**1.5

     enddo

     ! Always compute beta_* consistently with parameters in input.cgyro and then re-scale
     ! note: beta_star(0) will be over-written with sonic rotation
     call set_betastar
     beta_star(0)  = beta_star(0)*beta_star_scale
     beta_star_fac = beta_star_fac*beta_star_scale

  endif

  ! z_eff -- only use value from input.cgyro or input.gacode
  ! if simple electron Lorentz collisions (1,5) and z_eff_method=1
  ! else compute z_eff from the input ions' densities and charges
  if (.not. (z_eff_method == 1 .and. &
       (collision_model==1 .or. collision_model==5))) then
     z_eff = 0.0
     do is=1,n_species
        if (z(is) > 0.0) then 
           z_eff = z_eff+dens(is)*z(is)**2/dens_ele
        endif
     enddo
  endif

  ! Check electron species consistency
  num_ele = 0
  do is=1,n_species
     if (z(is) < 0) then
        num_ele = num_ele + 1
     endif
  enddo
  if(num_ele == 0) then
     if(ae_flag == 0) then
        call cgyro_error('No electron species specified')
        return
     endif
  else if(num_ele == 1) then
     if(ae_flag == 1) then
        call cgyro_error('Electron species specified with adiabatic electron flag')
        return
     endif
  else
     call cgyro_error('Only one electron species allowed')
     return
  endif

  if (udsymmetry_flag == 1) then
     zmag = 0.0 ; dzmag = 0.0
     shape_cos(:)   = 0.0
     shape_s_cos(:) = 0.0
  endif

  !-------------------------------------------------------------
  ! Manage simulation type (n=0,linear,nonlinear)
  !
  ! Note: nt1,nt2,nt_loc properly initialized in cgyro_mpi_grid
  !
  if (zf_test_mode > 0 .or. collision_test_mode>0) then

     if (zf_test_mode > 2) then
        ! The Apar and Bpar initial conditions could later be
        ! ZF_TEST_MODE = 3 and 4.  Currently these are not available.
        call cgyro_error('ZF_TEST_MODE > 2 is not available')
        return
     endif

     ! Zonal flow (n=0) test

     k_theta_base = q/rmin
     rho     = abs(ky/k_theta_base)*(-btccw)
     length  = abs(box_size/(s*k_theta_base))

     ! my_toroidal == 0
     ! k_theta == 0

     if (n_radial == 1) then
        call cgyro_info('ZF test with n_radial = 1 ; setting UP_RADIAL=0')
     else
        call cgyro_info('ZF test with n_radial > 1 ; setting UP_RADIAL=0')
     endif

     up_radial = 0.0
     
  else if (n_toroidal == 1) then

     ! Single linear mode (assume n=1, compute rho)

     k_theta_base = q/rmin
     rho     = abs(ky/k_theta_base)*(-btccw)
     length  = abs(box_size/(s*k_theta_base))

     ! my_toroidal == 1
     ! k_theta == k_theta_base

     call cgyro_info('Single-mode linear analysis')

  else

     ! Multiple modes (n=0,1,2,...,n_toroidal-1)

     k_theta_base = q/rmin
     rho     = abs(ky/k_theta_base)*(-btccw)
     length  = abs(box_size/(s*k_theta_base))

     ! Now define individual k_thetas

     ! my_toroidal == i_group_1
     ! k_theta == my_toroidal*k_theta_base

     call cgyro_info('Multiple toroidal harmonics')

  endif
  !
  ! To account for abs() operations above, we must define a variable
  ! with sign of (qs) 
  !
  if (ipccw*btccw*sign(1.0,s) < 0.0) then
     sign_qs = -1
  else
     sign_qs = 1
  endif

  if (q*rho < 0.0) then
     call cgyro_info('Ion direction: omega > 0') 
  else
     call cgyro_info('Ion direction: omega < 0') 
  endif
  !-------------------------------------------------------------

  !------------------------------------------------------------------------
  ! ExB and profile shear
  !
  source_flag = 0
  if (abs(gamma_e) > 1e-10 .and. nonlinear_flag > 0) then
     omega_eb_base = k_theta_base*length*gamma_e/(2*pi)
     select case (shear_method)
     case (1)
        call cgyro_info('ExB shear: Hammett discrete shift (WARNING: TESTING PURPOSES ONLY)') 
     case (2)
        call cgyro_info('ExB shear: Wavenumber advection') 
        source_flag = 1
     case default
        call cgyro_error('Unknown ExB shear method') 
     end select
  else
     omega_eb_base = 0.0
     call cgyro_info('ExB shear: OFF') 
  endif

  if (global_flag == 1) then
     call cgyro_info('Global profile variation: ON') 
     source_flag = 1
  else
     sdlnndr(1:n_species) = 0.0
     sdlntdr(1:n_species) = 0.0
     sbeta                = 0.0
  endif


  !------------------------------------------------------------------------

  lambda_debye = lambda_star*rho

  !-------------------------------------------------------------
  ! Fourier index mapping
  !
  allocate(px(n_radial))
  if (collision_test_mode>0) then
     do ir=1,n_radial
        px(ir) = ir-1
     enddo
  else if (zf_test_mode > 0) then
     ! Need positive k_r Fourier coefficients only
     do ir=1,n_radial
        px(ir) = ir
     enddo
  else
     do ir=1,n_radial
        px(ir) = -n_radial/2 + (ir-1)
     enddo
  endif

  !-------------------------------------------------------------
#if defined(OMPGPU)
!$omp target enter data map(to:px,z,temp)
#elif defined(_OPENACC)
!$acc enter data copyin(px,z,temp)
#endif

end subroutine cgyro_make_profiles

subroutine set_betastar

  use cgyro_globals

  implicit none

  integer :: is

  ! This is dp/dr at theta=0
  ! Note that in rotation, this will be overwritten as dp_eff/dr (theta=0)
  
  beta_star(:) = 0.0
  do is=1,n_species
     beta_star(0) = beta_star(0) &
          + dens(is)*temp(is)/(dens_ele*temp_ele) &
          *(dlnndr(is)+dlntdr(is))
  enddo
  if (ae_flag == 1) then
     beta_star(0) = beta_star(0) + (dlnndr_ele + dlntdr_ele)
  endif
  beta_star(0)  = beta_star(0)*betae_unit
  ! 8*pi/Bunit^2 * scaling factor
  beta_star_fac = -betae_unit/(dens_ele*temp_ele)
  
end subroutine set_betastar
