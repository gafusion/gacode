subroutine cgyro_make_profiles

  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: is,ir
  integer :: j
  integer :: num_ele

  !-------------------------------------------------------------
  ! Adiabatic-electron logic
  if (n_species == 1) then
     ae_flag = 1
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

  ! Compute effective charge (diagnostic)
  z_eff = 0.0
  do is=1,n_species
     if (z(is) > 0.0) then 
        z_eff = z_eff+dens(is)*z(is)**2/dens_ele
     endif
  enddo

  ! Check electron species consistency
  num_ele = 0
  do is=1,n_species
     if (z(is) < 0) then
        num_ele = num_ele + 1
     endif
  enddo
  if (num_ele == 0) then
     if (ae_flag == 0) then
        call cgyro_error('No electron species specified')
        return
     endif
  else if (num_ele == 1) then
     if (ae_flag == 1) then
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
  if (zf_test_mode > 0 .or. collision_test_mode > 0) then
 
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
  if (collision_test_mode>0) then
     px_zero = -(-1)
  else if (zf_test_mode > 0) then
     ! Need positive k_r Fourier coefficients only
     px_zero = 0
  else
     px_zero = -(-n_radial/2 - 1)
  endif
  ! Historical note
  ! px(ir) := ir-px_zero

  !-------------------------------------------------------------
#if defined(OMPGPU)
  !$omp target enter data map(to:z,temp)
#elif defined(_OPENACC)
  !$acc enter data copyin(z,temp)
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
