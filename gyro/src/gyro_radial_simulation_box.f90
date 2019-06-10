!-----------------------------------------------------------
! gyro_radial_simulation_box.f90 [caller make_profiles]
!
! PURPOSE:
!  Determine box length from input parameters.   
!
! NOTES:
!  See documentation at:
!
!      https://fusion.gat.com/theory/BOX_MULTIPLIER
!------------------------------------------------------------------------

subroutine gyro_radial_simulation_box

  use gyro_globals
  use math_constants

  !--------------------------
  implicit none
  !--------------------------

  !------------------------------------------------------
  ! Information and sanity check for automated selection 
  ! of TOROIDAL_REFERENCE
  !
  if (n_ref <= 0) then
     call send_message('INFO: (GYRO) You are using TOROIDAL_REFERENCE=-1')
     call send_message('INFO: (GYRO) TOROIDAL_REFERENCE reset to lowest nonzero mode.')
     n_ref = n0
     if (n_ref == 0) n_ref = d_n
  endif
  !------------------------------------------------------

  !---------------------------------------------------------
  ! If box_multiplier < 0, recompute d_n and box_multiplier
  ! based on values set by l_x and l_y.
  !
  if (box_multiplier < 0.0 .and. radial_profile_method /= 3) then

     if (nonlinear_flag == 1) then

        ! NONLINEAR: use l_x,l_y

        rho_star = abs(2*pi*r0/(n_ref*q0*l_y))*(-btccw)
        box_multiplier = abs(2*pi*s0*l_x/l_y)

     else

        ! LINEAR: use k_theta rho_s -> l_y

        rho_star = abs(l_y*r0/(n_ref*q0))*(-btccw)
        box_multiplier = abs(box_multiplier)

     endif

  endif
  !---------------------------------------------------------

  !------------------------------------------------------
  ! Set box length (might be reset later)
  !
  x_length = abs(box_multiplier*r0/(n_ref*q0*s0))
  !------------------------------------------------------

  !------------------------------------------------------
  ! Select radial grid-spacing and gridpoints
  !
  select case (boundary_method)

  case (1) 

     ! PERIODIC

     ! Box ends at periodic point
     !
     ! x----x----x----x----p
     !           |
     !          r=r0

     d_x = x_length/n_x

  case (2)

     ! NONPERIODIC

     ! Box ends at zero-function boundaries; 
     ! all points interior.
     !
     ! x-----x-----x-----x
     !          |
     !         r=r0

     d_x = x_length/(n_x-1.0)

  end select

  ! General max-width "buffers"
  i1_buffer = 1-max(m_gyro,m_dx)
  i2_buffer = n_x+max(m_gyro-i_gyro,m_dx-i_dx)

  ! Derivative "buffers"
  i1_dx = 1-m_dx
  i2_dx = n_x+m_dx-i_dx
  !
  do i=1,n_x
     r(i) = r0-0.5*x_length+d_x*(i-1.0)
  enddo
  !
  ! Check for legal box size
  !
  if (r(n_x) > 1.0) then
     if (radial_profile_method == 3 .and. boundary_method == 1) then
        x_length = x_length/10.0
        d_x = x_length/n_x
        do i=1,n_x
           r(i) = r0-0.5*x_length+d_x*(i-1.0)
        enddo
        call send_message('INFO: (GYRO) Shrinking domain by factor 10 avoid domain-size error.')
     else
        call catch_error('ERROR: (GYRO) Radial domain too large (r/a > 1).')
     endif
  endif
  !------------------------------------------------------

  !----------------------------------------------
  ! Initialize r to the uniform grid.  This will 
  ! change if nonuniform_grid_flag > 0
  !
  r_e(:) = r(:)
  !----------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_simulation_box done]'
  endif

end subroutine gyro_radial_simulation_box
