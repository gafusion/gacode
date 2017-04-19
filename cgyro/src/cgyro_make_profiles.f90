subroutine cgyro_make_profiles

  use cgyro_globals
  use cgyro_io
  use cgyro_experimental_globals

  implicit none

  integer :: is,ir,ix 
  integer :: j
  integer :: num_ele

  real :: cc,loglam

  !-------------------------------------------------------------
  ! Manage electrons
  !
  num_ele = 0
  do is=1,n_species
     if (z(is) < 0) then
        num_ele = num_ele + 1
        is_ele = is
     endif
  enddo

  if (num_ele == 0) then
     ! Adiabatic electrons
     ae_flag = 1
     call cgyro_info('Using adiabatic electrons')
  else if (num_ele == 1) then
     ! GK electrons
     ae_flag = 0
     call cgyro_info('Using gyrokinetic electrons')
  else
     call cgyro_error('Only one electron species allowed')
     return
  endif

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
     ! Miller
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

     call cgyro_experimental_profiles
     call cgyro_experimental_map
     call cgyro_experimental_alloc(0)

     if (ae_flag == 1) then
        dens_ele = ne_ade
        temp_ele = te_ade
        mass_ele = masse_ade
        dlnndr_ele = dlnndre_ade
        dlntdr_ele = dlntdre_ade
     else 
        dens_ele = dens(is_ele)
        temp_ele = temp(is_ele)
        mass_ele = mass(is_ele)
        dlnndr_ele = dlnndr(is_ele)
        dlntdr_ele = dlntdr(is_ele)
     endif

     ! Normalizing quantities
     dens_norm = dens_ele
     temp_norm = temp_ele

     ! Compute vth (m/s) using dimensional quantities.  
     ! mass(i) is thus measured in units of deuterium mass.
     do is=1,n_species
        vth(is) = sqrt(temp(is) * temp_norm_fac &
             / (mass(is) * mass_deuterium)) * 1.0e4 
     enddo
     vth_norm  = sqrt(temp_ele * temp_norm_fac &
          / (mass_deuterium)) * 1.0e4

     ! Compute collision frequency

     cc = sqrt(2.0) * pi * charge_norm_fac**4 &
          * 1.0 / (4.0 * pi * 8.8542)**2 &
          * 1.0 / (sqrt(mass_deuterium) * temp_norm_fac**1.5) &
          * 1e9

     loglam = 24.0 - log(sqrt(dens_ele*1e13)/(temp_ele*1e3))
     nu_ee  = cc * loglam * dens_ele / (sqrt(mass_ele)*temp_ele**1.5) &
          / (vth_norm/a_meters) 

     ! beta calculation in CGS:
     !
     !         8*pi ( n[1e19/m^3]*1e-6*1e19 )( T[keV]*1.6022*1e-9 )
     ! beta = ------------------------------------------------------
     !                           ( 1e4*B[T] )^2
     !
     !      = 4.027e-3 n[1e19/m^3]*T[keV]/B[T]^2

     betae_unit = 4.027e-3 * dens_ele * temp_ele / b_unit**2

     ! Debye length (from NRL plasma formulary):
     ! Use input lambda_debye as scaling parameter

     lambda_star = 7.43 * sqrt((1e3*temp_norm)/(1e13*dens_norm))/rhos

     ! Normalize
     do is=1,n_species
        dens(is) = dens(is)/dens_norm
        temp(is) = temp(is)/temp_norm
        vth(is)  = vth(is)/vth_norm
     enddo
     ne_ade   = ne_ade/dens_norm
     te_ade   = te_ade/temp_norm
     dens_ele = dens_ele/dens_norm
     temp_ele = temp_ele/temp_norm
     gamma_e  = gamma_e/(vth_norm/a_meters)
     gamma_p  = gamma_p/(vth_norm/a_meters)
     mach     = mach/vth_norm

     do is=1,n_species
        nu(is) = nu_ee *z(is)**4 &
             * dens(is) / dens_ele &
             * sqrt(mass_ele/mass(is)) * (temp_ele/temp(is))**1.5
     enddo

     ! Always compute beta_* consistently
     call set_betastar
     
     ! Re-scaling
     lambda_star      = lambda_star * lambda_star_scale
     gamma_e          = gamma_e      * gamma_e_scale
     gamma_p          = gamma_p      * gamma_p_scale
     mach             = mach         * mach_scale
     q                = q            * q_scale
     s                = s            * s_scale
     shift            = shift        * shift_scale
     kappa            = kappa        * kappa_scale
     delta            = delta        * delta_scale
     zeta             = zeta         * zeta_scale
     s_kappa          = s_kappa      * s_kappa_scale
     s_delta          = s_delta      * s_delta_scale
     s_zeta           = s_zeta       * s_zeta_scale
     beta_star(0)     = beta_star(0) * beta_star_scale
     betae_unit       = betae_unit   * betae_unit_scale
     do is=1,n_species
        dlnndr(is) = dlnndr(is)  * dlnndr_scale(is) 
        dlntdr(is) = dlntdr(is)  * dlntdr_scale(is)  
        nu(is)     = nu(is)      * nu_ee_scale
     enddo

  else

     a_meters  = 1.0
     b_unit    = 1.0
     dens_norm = 1.0
     temp_norm = 1.0
     vth_norm  = 1.0
     
     q = abs(q)*(ipccw)*(btccw)

     if (ae_flag == 1) then
        dens_ele = ne_ade
        temp_ele = te_ade
        mass_ele = masse_ade
        dlnndr_ele = dlnndre_ade
        dlntdr_ele = dlntdre_ade
     else 
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

     ! Always compute beta_* consistently
     call set_betastar
     beta_star(0) = beta_star(0)*beta_star_scale
     
  endif

  
  !-------------------------------------------------------------
  ! Manage simulation type (n=0,linear,nonlinear)
  !
  if (zf_test_flag == 1) then

     ! Zonal flow (n=0) test

     k_theta = q/rmin
     rho     = abs(ky/k_theta)*(-btccw)
     length  = abs(box_size/(s*k_theta))

     k_theta = 0

     call cgyro_info('Triggered zonal flow test')

     if (n_radial /= 1) then
        call cgyro_info('Zonal flow test with n_radial>1')
     endif

     n = 0

  else if (n_toroidal == 1) then

     ! Single linear mode (assume n=1, compute rho)

     k_theta = q/rmin
     rho     = abs(ky/k_theta)*(-btccw)
     length  = abs(box_size/(s*k_theta))

     n = 1

     call cgyro_info('Single-mode linear analysis')

  else

     ! Multiple modes (n=0,1,2,...,n_toroidal-1)

     k_theta = q/rmin
     rho     = abs(ky/k_theta)*(-btccw)
     length  = abs(box_size/(s*k_theta))

     ! Now define individual k_thetas

     n = i_group_1

     k_theta = n*k_theta

     call cgyro_info('Multiple toroidal harmonics')

  endif
  !-------------------------------------------------------------

  !------------------------------------------------------------------------
  ! ExB and profile shear
  !
  if (abs(gamma_e) > 1e-10 .and. nonlinear_flag > 0) then
     omega_eb = k_theta*length*gamma_e/(2*pi)
     select case (shear_method)
     case (1)
        call cgyro_info('ExB shear: Hammett discrete shift') 
     case (2)
        call cgyro_info('ExB shear: Explicit wavenumber advection') 
     case (3)
        call cgyro_info('ExB shear: Linearized Hammett shift') 
     case (4)
        call cgyro_info('ExB shear: Implicit wavenumber advection') 
     case (5)
        call cgyro_info('ExB shear: Compact-difference wavenumber advection') 
     end select
  else
     omega_eb = 0.0
     shear_method = 0
     call cgyro_info('ExB shear: OFF') 
  endif

  if (profile_shear_flag == 1) then
     call cgyro_info('Profile shear: Continuous wavenumber advection') 
  endif
  !------------------------------------------------------------------------

  lambda_debye = lambda_star*rho

  !-------------------------------------------------------------
  ! Fourier index mapping
  !
  allocate(indx_xi(n_xi))
  do ix=1,n_xi
     indx_xi(ix) = ix-1
  enddo
  allocate(px(n_radial))
  if (zf_test_flag == 1) then
     do ir=1,n_radial
        px(ir) = ir
        ! only need positive k_r Fourier coefficients.
     enddo
  else
     do ir=1,n_radial
        px(ir) = -n_radial/2 + (ir-1)
     enddo
  endif

  !-------------------------------------------------------------

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
     beta_star(0) = beta_star(0) + (dlnndre_ade + dlntdre_ade)
  endif
  beta_star(0)  = beta_star(0)*betae_unit
  ! 8pi/Bunit^2 * scaling factor
  beta_star_fac = -betae_unit/(dens_ele*temp_ele)*beta_star_scale  
  
end subroutine set_betastar
