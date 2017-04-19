subroutine cgyro_check

  use cgyro_globals
  use cgyro_io

  implicit none

  integer :: is
  logical :: lfe
  character(len=1), dimension(7) :: ctag

  !-----------------------------------------------------------------------
  ! Grid parameter checks
  !
  if (modulo(n_xi,2) /= 0) then 
     call cgyro_error('n_xi must be even.')
     return
  endif

  if (zf_test_flag == 0 .and. modulo(n_radial,box_size) /= 0) then 
     call cgyro_info('RESOLUTION WARNING -- n_radial not a multiple of box_size.')
  endif

  if (n_radial < (n_toroidal-1)*box_size .and. zf_test_flag == 0) then
     call cgyro_info('RESOLUTION WARNING -- n_radial < n*box_size.')
  endif

  if (n_radial < box_size .and. zf_test_flag == 0) then
     call cgyro_info('SEVERE RESOLUTION WARNING -- n_radial < box_size.')
     return
  endif

  if (n_species > 6) then
     call cgyro_error('n_species <= 6.')
     return
  endif

  if (box_size > 1 .and. h_print_flag == 1) then
     call cgyro_error('Distribution output not available for box_size > 1')
     return
  endif
  !------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Profile checks
  !
  select case(profile_model)

  case (1)
     call cgyro_info('Profile model: local input (input.cgyro)')

  case (2)
     call cgyro_info('Profile model: experimental (input.profiles)')

  case default
     call cgyro_error('Invalid value for profile_model')
     return

  end select
  !-----------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Equilibrium model
  !
  select case (equilibrium_model)  

  case (1) 
     call cgyro_info('Equilibrium: s-alpha')
     if (profile_model == 2) then
        call cgyro_error('s-alpha equilibrium model not valid with experimental profiles')
        return
     endif

  case (2) 
     call cgyro_info('Equilibrium: Miller')

  case (3) 
     call cgyro_info('Equilibrium: General (Fourier)')

     if (geo_ny <= 0) then
        call cgyro_error('Fourier geometry coefficients missing.')
        return
     endif
     if (udsymmetry_flag == 1) then
        call cgyro_error('Cannot have UDSYMMETRY_FLAG=1 with general geometry.')
        return
     endif

  case default

     call cgyro_error('Invalid value for equilibrium_model')
     return

  end select
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Field consistency checks
  if (n_field > 1) then
     if (abs(betae_unit) < epsilon(0.0)) then
        call cgyro_error('BETAE_UNIT must be nonzero for electromagnetic simulation.')
        return
     endif
     if (ae_flag == 1) then
        call cgyro_error('Electrons must be gyrokinetic for electromagnetic simulation.')
        return
     endif
  endif
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Fields
  select case(n_field)

  case (1)
     call cgyro_info('Electrostatic fluctuations (Phi)')

  case (2)
     call cgyro_info('Transverse EM fluctuations (Phi,A_par)')

  case (3)
     call cgyro_info('Transverse and compressional EM fluctuations (Phi,A_par,B_par)')

  case default
     call cgyro_error('Invalid value for n_field')
     return
  end select
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Rotation model
  !
  select case (rotation_model)

  case(1)
     call cgyro_info('Rotation terms: O(mach) only (traditional GYRO)')
  case(2)
     call cgyro_info('Rotation terms: O(mach + mach^2) (full Sugama rotation)')
  case(3)
     call cgyro_info('Rotation terms: O(mach^2) only (isolated centrifugal for testing)')
  case(4)
     call cgyro_info('Rotation terms: O(mach^2) GKW CF TRAP term only (for testing)') 
  case(5)
     call cgyro_info('Rotation terms: O(mach^2) GKW CF DRIFT term only (for testing)')   

  case default
     call cgyro_error('Invalid value for rotation_model')
     return

  end select
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Collision model and settings
  !
  ctag(:) = ' '
  ctag(1) = 'x' 

  select case (collision_model)  

  case (1) 
     call cgyro_info('Collision model: Lorentz ee+ei')
  case (2) 
     call cgyro_info('Collision model: Connor')
  case (4) 
     call cgyro_info('Collision model: Sugama')
     ctag(2) = 'x'
  case(5)
     call cgyro_info('Collision model: Simple Lorentz ee+ei')
  case default
     call cgyro_error('Invalid value for collision_model')
     return

  end select

  if (collision_model /= 5) then

     if(collision_model /= 1) then
        select case (collision_mom_restore)
        case(0)
           ! Collision momentum restoring: off       
        case(1)
           ! Collision momentum restoring: on
           ctag(3)= 'x'
        case default
           call cgyro_error('Invalid value for collision_mom_restore')
           return
        end select
     endif

     select case (collision_field_model)
     case(0)
        ! Collision field corrections : off
     case (1)
        ! Collision field corrections : on
        ctag(7) = 'x'
     case default
        call cgyro_error('Invalid value for collision_field_model')
        return
     end select
  endif

  if (collision_model == 4) then

     select case (collision_ene_diffusion)
     case(0)
        ! Collision energy diffusion  : off
     case(1)
        ! Collision energy diffusion  : on
        ctag(2) = 'x'
     case default
        call cgyro_error('Invalid value for collision_ene_diffusion')
        return
     end select

     select case (collision_ene_restore)
     case(0)
        ! Collision energy restoring  : off
     case(1)
        ! Collision energy restoring  : on
        ctag(4) = 'x'
     case default
        call cgyro_error('Invalid value for collision_ene_restore')
        return
     end select

     select case (collision_kperp)
     case(0)
        ! Collision kperp corrections : off
     case(1)
        ! Collision kperp corrections : on
        ctag(5) = 'x'
     case default
        call cgyro_error('Invalid value for collision_kperp')
        return
     end select

  endif

  if (collision_model /= 5 .and. collision_model /= 1) then
     select case (collision_ion_model)
     case(0)
        ! Collisional ions: on
        ctag(6) = 'x'
     case(1)
        ! Collisional ions: off
     case default
        call cgyro_error('Invalid value for collision_ion_model')
        return
     end select
  endif

  call cgyro_info('Collision terms: L D Rm Re kp ions field')
  call cgyro_info('               '// &
       '  '//ctag(1)// &
       ' '//ctag(2)// &
       ' '//ctag(3)// &
       '  '//ctag(4)// &
       '  '//ctag(5)// &
       '   '//ctag(6)// &
       '     '//ctag(7))
  !------------------------------------------------------------------------


  if (rotation_model > 1 .and. collision_model == 5) then
     call cgyro_error('Simple collisions not available with rotation_model > 1')
     return
  endif

  !------------------------------------------------------------------------
  ! Check profile parameters
  !
  do is=1,n_species
     if (dens(is) <= 0.0) then
        call cgyro_error('Densities must be positive.')
        return
     endif
     if (temp(is) <= 0.0) then
        call cgyro_error('Temperatures must be positive.')
        return
     endif
     if (nu(is) < 0.0) then
        call cgyro_error('Collision frequencies must be non-negative.')
        return
     endif
     if (abs(z(is)) < epsilon(0.0)) then
        call cgyro_error('Charges must be non-zero.')
     endif
  enddo
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Dissipation checks
  !
  if (nup_theta < 1 .or. nup_theta > 3) then
     call cgyro_error('Invalid value for nup_theta')
     return
  endif
  if (nup_radial < 1 .or. nup_theta > 4) then
     call cgyro_error('Invalid value for nup_radial')
     return
  endif
  if (nup_alpha < 1 .or. nup_alpha > 4) then
     call cgyro_error('Invalid value for nup_alpha')
     return
  endif
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  ! Check for existence of restart file
  !
  if (restart_mode == 1) then

     inquire(file=trim(path)//runfile_restart,exist=lfe)
     if (lfe .eqv. .false.) then
        call cgyro_error('Missing restart file')
        return
     endif

  endif
  !------------------------------------------------------------------------

end subroutine cgyro_check
