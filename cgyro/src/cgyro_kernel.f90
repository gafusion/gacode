!-----------------------------------------------------------------
! cgyro_kernel.f90
!
! PURPOSE:
!  Subroutinized main cgyro program.  
!
! NOTES:
!  This can be called directly using the driver routine cgyro 
!  (in which case input data will read from input.dat) or called 
!  as a subroutine using cgyro_sub.
!-----------------------------------------------------------------

subroutine cgyro_kernel

  use timer_lib
  use mpi

  use cgyro_globals
  use cgyro_io

  implicit none

  character(len=30) :: final_msg

  i_time = 0

  ! Need to initialize the info runfile very early
  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//runfile_info,status='replace')
  endif

  ! 1. MPI setup
  call cgyro_mpi_grid
  if (error_status > 0) goto 100

  ! 2. Profile setup
  call cgyro_make_profiles
  if (error_status > 0) goto 100

  ! 3. Parameter consistency checks
  call cgyro_check
  if (error_status > 0) goto 100

  ! 4. Array initialization and construction
  !    NOTE: On exit, field_old = field 

  call cgyro_init_manager
  if (test_flag == 1) return

  !---------------------------------------------------------------------------
  !
  ! Time-stepping
  n_time = nint(max_time/delta_t)

  if (restart_flag == 0) then
     io_control = 1*(1-silent_flag)
  else
     io_control = 3*(1-silent_flag)
  endif

  call timer_lib_in('io_init')
  call cgyro_write_timedata
  call timer_lib_out('io_init')
  call write_timers(trim(path)//runfile_timers)

  io_control = 2*(1-silent_flag)

  do i_time=1,n_time

     call timer_lib_in('TOTAL')

     !------------------------------------------------------------
     ! Time advance
     !
     t_current = t_current+delta_t

     ! Collisionless step: returns new h_x, cap_h_x, fields 
     if (integration_error(2) > adapt_tol .and. nonlinear_flag == 1) then
        ! Trigger adaptive step
        delta_t = delta_t/4
        call cgyro_step_gk
        call cgyro_step_gk
        call cgyro_step_gk
        call cgyro_step_gk
        delta_t = 4*delta_t
     else
        ! Normal timestep
        call cgyro_step_gk
     endif

     ! Collision step: returns new h_x, cap_h_x, fields
     if (collision_model == 5) then
        call cgyro_step_collision_simple
     else
        call cgyro_step_collision
     endif

     ! Hammett method for ExB shear
     if (shear_method == 1) then
        ! Discrete shift (Hammett) 
        call timer_lib_in('shear')
        call cgyro_shear_hammett
        call timer_lib_out('shear')
     endif
     !------------------------------------------------------------

     !------------------------------------------------------------
     ! Diagnostics
     !
     ! NOTE: Fluxes are calculated in cgyro_write_timedata

     ! Error estimate
     call cgyro_error_estimate
     !------------------------------------------------------------

     !---------------------------------------
     ! IO
     !
     call timer_lib_in('io')

     ! Write simulation data
     call cgyro_write_timedata

     ! Write restart data
     call cgyro_write_restart

     call timer_lib_out('io')
     !---------------------------------------

     call timer_lib_out('TOTAL')

     ! Don't wrap timer output in a timer
     if (mod(i_time,print_step) == 0) call write_timers(trim(path)//runfile_timers)

     ! Trap error or convergence
     if (error_status > 0 .or. signal == 1) exit

  enddo
  !---------------------------------------------------------------------------

100 continue

  ! Manage exit message

  if (error_status == 0) then
     if (nonlinear_flag == 1) then
        final_msg = 'Normal'
     else
        if (signal == 1) then
           final_msg = 'Linear converged'
        else
           final_msg = 'Linear terminated at max time'
        endif
     endif
     if (silent_flag == 0 .and. i_proc == 0) then
        open(unit=io,file=trim(path)//runfile_info,status='old',position='append')
        write(io,'(a)') 'EXIT: (CGYRO) '//trim(final_msg)
        close(io)
     endif
  endif
 
  if(allocated(theta))          deallocate(theta)
  if(allocated(thetab))         deallocate(thetab)
  if(allocated(w_theta))        deallocate(w_theta)
  if(allocated(g_theta))        deallocate(g_theta)
  if(allocated(g_theta_geo))    deallocate(g_theta_geo)
  if(allocated(bmag))           deallocate(bmag)
  if(allocated(btor))           deallocate(btor)
  if(allocated(bpol))           deallocate(bpol)
  if(allocated(k_perp))         deallocate(k_perp)
  if(allocated(k_x))            deallocate(k_x)
  if(allocated(bigr))           deallocate(bigr)
  if(allocated(bigr_r))         deallocate(bigr_r)
  if(allocated(omega_stream))   then
!$acc exit data delete(omega_stream)
     deallocate(omega_stream)
  endif
  if(allocated(omega_trap))     deallocate(omega_trap)
  if(allocated(omega_rdrift))   deallocate(omega_rdrift)
  if(allocated(omega_adrift))   deallocate(omega_adrift)
  if(allocated(omega_aprdrift)) deallocate(omega_aprdrift)
  if(allocated(omega_cdrift))   deallocate(omega_cdrift)
  if(allocated(omega_cdrift_r)) deallocate(omega_cdrift_r)
  if(allocated(omega_gammap))   deallocate(omega_gammap)

  if(allocated(lambda_rot))          deallocate(lambda_rot)
  if(allocated(dlambda_rot))         deallocate(dlambda_rot)
  if(allocated(dens_rot))            deallocate(dens_rot)
  if(allocated(dens_ele_rot))        deallocate(dens_ele_rot)
  if(allocated(omega_rot_trap))      deallocate(omega_rot_trap)
  if(allocated(omega_rot_u))         deallocate(omega_rot_u)
  if(allocated(omega_rot_drift))     deallocate(omega_rot_drift)
  if(allocated(omega_rot_drift_r))   deallocate(omega_rot_drift_r)
  if(allocated(omega_rot_edrift))    deallocate(omega_rot_edrift)
  if(allocated(omega_rot_edrift_r))  deallocate(omega_rot_edrift_r)
  if(allocated(omega_rot_star))      deallocate(omega_rot_star)

  if(allocated(px))            deallocate(px)
  if(allocated(energy))        then
!$acc exit data delete(energy)
     deallocate(energy)
  endif
  if(allocated(w_e))           deallocate(w_e)
  if(allocated(e_deriv1_mat))  deallocate(e_deriv1_mat)
  if(allocated(xi))            then
!$acc exit data delete(xi)
     deallocate(xi)
  endif
  if(allocated(w_xi))          deallocate(w_xi)
  if(allocated(xi_lor_mat))    deallocate(xi_lor_mat)
  if(allocated(xi_deriv_mat))  deallocate(xi_deriv_mat)
  if(allocated(h_x))           deallocate(h_x)
  if(allocated(cap_h_c))       deallocate(cap_h_c)
  if(allocated(cap_h_v))       deallocate(cap_h_v)
  if(allocated(field))         deallocate(field)
  if(allocated(field_loc))     deallocate(field_loc)
  if(allocated(field_old))     deallocate(field_old)
  if(allocated(hzf))           deallocate(hzf)
  if(allocated(xzf))           deallocate(xzf)

  if (allocated(cmat)) deallocate(cmat)

end subroutine cgyro_kernel

