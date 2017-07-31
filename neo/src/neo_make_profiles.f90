subroutine neo_make_profiles

  use neo_globals
  use neo_profile_exp, only : &
       rmaj_p,&
       b_unit_p,& 
       q_exp,&
       n_grid_exp
  use neo_allocate_profile
  
  implicit none
  
  integer :: ir, is, ip, num_ele, j
  integer, parameter :: io=20
  
  real :: cc
  real :: loglam
  
  call PROFILE_SIM_alloc(1)
  
  ! equally-spaced radial grid if only first-order (local) problem
  if(n_radial == 1) then
     r(1) = rmin_1_in
  else
     do ir=1,n_radial
        r(ir) = rmin_1_in+(rmin_2_in-rmin_1_in)*(ir-1.0)/(n_radial-1)
     enddo
  endif
  
  num_ele = 0
  do is=1, n_species
     if(z_in(is) < 0.0) then
        num_ele = num_ele + 1
     endif
  enddo
  if(num_ele == 0) then
     adiabatic_ele_model = 1
  else if(num_ele == 1) then
     adiabatic_ele_model = 0
  else
     call neo_error('ERROR: (NEO) Only one electron species allowed')
     return
  end if
  
  select case (profile_model) 
     
  case (1)
     
     ! Standard local simulation (one point)
     
     if(btccw_in > 0) then
        sign_bunit = -1.0
     else
        sign_bunit =  1.0
     endif
     
     if(ipccw_in > 0) then
        sign_q = -sign_bunit
     else
        sign_q =  sign_bunit
     endif
     
     ir = 1
     
     rmaj(ir)      = rmaj_in
     q(ir)         = abs(q_in) * sign_q
     rho(ir)       = abs(rho_in) * sign_bunit
     shear(ir)     = shear_in      
     shift(ir)     = shift_in     
     kappa(ir)     = kappa_in     
     s_kappa(ir)   = s_kappa_in   
     delta(ir)     = delta_in    
     s_delta(ir)   = s_delta_in
     zeta(ir)      = zeta_in    
     s_zeta(ir)    = s_zeta_in
     zmag(ir)      = zmag_in    
     s_zmag(ir)    = s_zmag_in
     beta_star(ir) = beta_star_in
     
     ! general geometry -- accessible only from interface 
     ! via parameters geo_ny_in and geo_yin_in
     if(equilibrium_model == 3) then
        geo_numeq_flag = 1
        deallocate(geo_yin)
        geo_ny = geo_ny_in   
        allocate(geo_yin(8,0:geo_ny,1))
        do j=0,geo_ny
           geo_yin(:,j,1) = geo_yin_in(:,j)
        enddo
     endif

     dphi0dr(ir)   = dphi0dr_in
     epar0(ir)     = epar0_in

     te_ade(ir)      = te_ade_in
     ne_ade(ir)      = ne_ade_in
     dlntdre_ade(ir) = dlntdre_ade_in
     dlnndre_ade(ir) = dlnndre_ade_in

     omega_rot(ir)       = omega_rot_in 
     omega_rot_deriv(ir) = omega_rot_deriv_in 

     do is=1,n_species
        z(is)           = z_in(is)
        mass(is)        = mass_in(is)
        dens(is,ir)     = dens_in(is)
        dlnndr(is,ir)   = dlnndr_in(is)
        aniso_model(is) = aniso_model_in(is)
        if(aniso_model(is) == 2) then
           !temp(is,ir)        = temp_in(is)
           !dlntdr(is,ir)      = dlntdr_in(is)
           temp_para(is,ir)   = temp_para_in(is)
           temp_perp(is,ir)   = temp_perp_in(is)
           dlntdr_para(is,ir) =  dlntdr_para_in(is)
           dlntdr_perp(is,ir) =  dlntdr_perp_in(is)
           ! for anisotropic species, define a Teff = 1/3 Tpar + 2/3 Tperp
           temp(is,ir)   = 1.0/3.0*temp_para(is,ir) + 2.0/3.0*temp_perp(is,ir)
           dlntdr(is,ir) = 1.0/3.0*temp_para(is,ir)/temp(is,ir) &
                *dlntdr_para(is,ir) + 2.0/3.0*temp_perp(is,ir) &
                /temp(is,ir) *dlntdr_perp(is,ir)
        else
           temp(is,ir)        = temp_in(is)
           dlntdr(is,ir)      = dlntdr_in(is)
           temp_para(is,ir)   = temp(is,ir)
           temp_perp(is,ir)   = temp(is,ir)
           dlntdr_para(is,ir) = dlntdr(is,ir)
           dlntdr_perp(is,ir) = dlntdr(is,ir)
        endif
        nu(is,ir)     = nu_1_in * z(is)**4 / z(1)**4 &
             * dens(is,ir) / dens(1,ir) &
             * sqrt(mass(1)/mass(is)) * (temp(1,ir)/temp(is,ir))**1.5
     enddo

     do is=1, n_species
        vth(is,ir) = sqrt(temp(is,ir)/mass(is))
        vth_para(is,ir) = sqrt(temp_para(is,ir)/mass(is))
     enddo

     ! These normalizations are arbitrary for local profiles
     temp_norm_fac = 1.0
     charge_norm_fac = 1.0
     dens_norm(:) = 1.0
     vth_norm(:) = 1.0
     b_unit(:) = 1.0
     a_meters        = 1.0
     temp_norm(:)    = 1.0
     psiN_polflux(:) = 0.0
     psiN_polflux_a  = 0.0


  case (2)
     
     ! Standard simulation with experimental profiles
     
     ! EAB: Epar is not in INPUT_profiles -- use INPUT val and assume
     ! radially constant
     do ir=1, n_radial
        epar0(ir) = epar0_in
     enddo
     
     do is=1,n_species
        z(is)    = z_in(is)
        mass(is) = mass_in(is)
        aniso_model(is) = aniso_model_in(is)
     enddo
     
     call neo_experimental_profiles
     if(error_status > 0) return
     call neo_map_experimental_profiles
           
     ! ** Normalizing quantities **
     temp_norm_fac   = 1.6022*1000
     charge_norm_fac = 1.6022
     dens_norm(:) = dens(1,:)
     temp_norm(:) = temp(1,:)
     
     ! Compute vth/a (1/s) using dimensional quantities.  
     ! mass(i) is thus measured in units of deuterium mass.
     do is=1,n_species
        vth(is,:) = sqrt(temp(is,:) * temp_norm_fac &
             / (mass(is) * mass_deuterium)) &
             * 1.0e4 / a_meters
     enddo
     vth_norm(:)  = vth(1,:) * sqrt(mass(1))
     
     ! Determine the equilibrium parameters
     select case (profile_equilibrium_model) 
     case (0)
        ! input equilibrium not used -- re-set for s-alpha geometry
        ! use the input profile-averaged rmaj, q, and b_unit
        equilibrium_model = 0
        do ir=1, n_radial
           rmaj(ir)   = 0.0
           q(ir)      = 0.0
           b_unit(ir) = 0.0
           do ip=1,n_grid_exp
              rmaj(ir)   = rmaj(ir)   + rmaj_p(ip)
              q(ir)      = q(ir)      + q_exp(ir)
              b_unit(ir) = b_unit(ir) + b_unit_p(ip)
           enddo
           rmaj(ir)   = rmaj(ir)   / (1.0 * n_grid_exp)
           q(ir)      = q(ir)      / (1.0 * n_grid_exp)
           b_unit(ir) = b_unit(ir) / (1.0 * n_grid_exp)
        enddo
        
     case(1)
        ! use the input equilibrium parameters (miller)
        equilibrium_model = 2
        
     case(2)
        ! use the input equilibrium parameters (general)
        equilibrium_model = 3
        
     end select

     ! EAB: note about beta_star
     ! This is presently only used for anisotropic species, which currently
     ! only works in local profile mode
     beta_star(:) = 0.0
     
     ! Compute the rotation parameters
     select case (rotation_model)     
     case (1)
        ! no rotation effects
        omega_rot(:) = 0.0
        omega_rot_deriv(:) = 0.0
     case (2)
        ! normalize
        omega_rot(:) = omega_rot(:) / vth_norm(:)
        omega_rot_deriv(:) = omega_rot_deriv(:) / vth_norm(:) * a_meters
        
     end select
     
     ! Compute the equilibrium radial electric field
     select case (profile_erad0_model) 
     case (0)
        ! er0 not included -- reset to be zero
        dphi0dr(:) = 0.0
     case (1)
        ! input erad0 is included
        ! Normalize erad0 from exp profiles
        dphi0dr(:) = dphi0dr(:) * a_meters / temp_norm(:) / 1000
     end select
     
     ! Compute rho/a for using dimensional quantities
     ! mass(i) is thus measured in units of deuterium mass.
     ! EAB: Note on 05/09/16 fixed error in which this had
     ! an extra factor of sqrt(mass(1))
     rho(:) = sqrt(temp(1,:) * temp_norm_fac &
          * mass_deuterium) &
          / (charge_norm_fac * b_unit(:)) &
          * 1.0e-4 / a_meters
     
     ! Compute collision frequency
     
     ! Numerical coefficient (relative to M_D)
     cc = sqrt(2.0) * pi * charge_norm_fac**4 &
          * 1.0 / (4.0 * pi * 8.8542)**2 &
          * 1.0 / (sqrt(mass_deuterium) * temp_norm_fac**1.5) &
          * 1e9

     do is=1,n_species
        do ir=1, n_radial
           
           ! Coulomb logarithm
           ! EAB: 03/22/09 redefined this wrt electron species
           ! (was previously defined wrt species 1)
           loglam = 24.0 - log(sqrt(ne_ade(ir)*1e13)/(te_ade(ir)*1000))

           ! Collision rate (1/sec)
           nu(is,ir) = cc * loglam * dens(is,ir) * z(is)**4 &
                / (sqrt(mass(is)) * temp(is,ir)**1.5)
           
           ! Express in local dimensionless NEO units:
           nu(is,ir) = nu(is,ir)/vth_norm(ir) 
           
        enddo
     enddo
     
     do is=1,n_species
        do ir=1, n_radial
           ! Normalize N and T to value at r for species 1.
           dens(is,ir) = dens(is,ir)/dens_norm(ir)
           temp(is,ir) = temp(is,ir)/temp_norm(ir)
           if(is == 1) then
              te_ade(ir) = te_ade(ir) / temp_norm(ir)
              ne_ade(ir) = ne_ade(ir) / dens_norm(ir)
           endif
           ! Normalize vth/a
           vth(is,ir) = vth(is,ir)/vth_norm(ir)
        enddo
     enddo
     
     ! Anisotropic mode is currently not available for global profiles
     aniso_model(:)   = 1
     temp_perp(:,:)   = temp(:,:)
     temp_para(:,:)   = temp(:,:)
     dlntdr_perp(:,:) = dlntdr(:,:)
     dlntdr_para(:,:) = dlntdr(:,:)
     vth_para(:,:)    = vth(:,:)

     call PROFILE_EXP_alloc(0)
     
  end select
  
  if(rotation_model == 2) then
     ! In the strong rotation limit, only use phi_(-1) (set by omega_rot)
     ! phi_(0) is used only in the weak rotation limit
     dphi0dr(:) = 0.0
  endif
  
  
  ! Print the re-mapped equilibrium data
  if(silent_flag == 0 .and. i_proc == 0) then
     
     open(unit=io,file=trim(path)//'out.neo.equil',status='replace')
     do ir=1,n_radial
        write (io,'(e16.8)',advance='no') r(ir)
        write (io,'(e16.8)',advance='no') dphi0dr(ir)
        write (io,'(e16.8)',advance='no') q(ir)
        write (io,'(e16.8)',advance='no') rho(ir)
        write (io,'(e16.8)',advance='no') rmaj(ir)
        write (io,'(e16.8)',advance='no') omega_rot(ir)
        write (io,'(e16.8)',advance='no') omega_rot_deriv(ir)
        do is=1,n_species
           write (io,'(e16.8)',advance='no') dens(is,ir)
           write (io,'(e16.8)',advance='no') temp(is,ir)
           write (io,'(e16.8)',advance='no') dlnndr(is,ir)
           write (io,'(e16.8)',advance='no') dlntdr(is,ir)
           write (io,'(e16.8)',advance='no') nu(is,ir)
        enddo
        write (io,*)
     enddo
     close(io)
     
     if(profile_model >= 2) then
        open(unit=io,file=trim(path)//'out.neo.expnorm',status='replace')
        do ir=1,n_radial
           write (io,'(e16.8)',advance='no') r(ir)
           write (io,'(e16.8)',advance='no') a_meters
           write (io,'(e16.8)',advance='no') mass_deuterium
           write (io,'(e16.8)',advance='no') dens_norm(ir)
           write (io,'(e16.8)',advance='no') temp_norm(ir)
           write (io,'(e16.8)',advance='no') vth_norm(ir)*a_meters
           write (io,'(e16.8)',advance='no') b_unit(ir)
           write (io,*)
        enddo
        close(io)
     end if

  endif
  
end subroutine neo_make_profiles
