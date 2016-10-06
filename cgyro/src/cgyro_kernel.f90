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

  !---------------------------------------------------------------------------
  !
  ! Time-stepping
  n_time = nint(max_time/delta_t)
  
  if (restart_flag == 0) then
     io_control = 1*(1-silent_flag)
  else
     io_control = 3*(1-silent_flag)
  endif
  call cgyro_write_timedata
  io_control = 2*(1-silent_flag)

  do i_time=1,n_time

     call timer_lib_in('TOTAL')
     
     !------------------------------------------------------------
     ! Time advance
     !
     t_current = t_current+delta_t

     ! Collisionless step: returns new h_x, cap_h_x, fields 
     call cgyro_step_gk
        
     ! Collisionless implicit streaming term step
     ! : returns new h_x, cap_h_x, fields 
     call cgyro_step_implicit_gk

     ! Collision step: returns new h_x, cap_h_x, fields
     call cgyro_step_collision
     !------------------------------------------------------------

     !------------------------------------------------------------
     ! Diagnostics
     !
     ! NOTE: Fluxes are calculated in cgyro_write_timedata

     ! Error estimate
     call cgyro_error_estimate
     !------------------------------------------------------------

     !------------------------------------------------------------
     ! Spectral ExB shear
     call cgyro_shear
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

     if (abs(signal) == 1 .or. error_status > 0) exit

  enddo
  !---------------------------------------------------------------------------

100 continue

  if(allocated(theta))          deallocate(theta)
  if(allocated(thetab))         deallocate(thetab)
  if(allocated(w_theta))        deallocate(w_theta)
  if(allocated(bmag))           deallocate(bmag)
  if(allocated(k_perp))         deallocate(k_perp)
  if(allocated(omega_stream))   then
!$acc exit data delete(omega_stream)
      deallocate(omega_stream)
  endif
  if(allocated(omega_trap))     deallocate(omega_trap)
  if(allocated(omega_rdrift))   deallocate(omega_rdrift)
  if(allocated(omega_adrift))   deallocate(omega_adrift)
  if(allocated(omega_aprdrift)) deallocate(omega_aprdrift)
  if(allocated(omega_cdrift))   deallocate(omega_cdrift)
  if(allocated(omega_gammap))   deallocate(omega_gammap)

  if(allocated(indx_xi))       deallocate(indx_xi)
  if(allocated(px))            deallocate(px)
  if(allocated(energy))        then
!$acc exit data delete(energy)
    deallocate(energy)
  endif
  if(allocated(w_e))           deallocate(w_e)
  if(allocated(e_deriv1_mat))  deallocate(e_deriv1_mat)
  if(allocated(e_deriv2_mat))  deallocate(e_deriv2_mat)
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
  if(allocated(f_balloon))     deallocate(f_balloon)
  if(allocated(hzf))           deallocate(hzf)
  if(allocated(xzf))           deallocate(xzf)
  if(allocated(pvec_outr))     deallocate(pvec_outr)
  if(allocated(pvec_outi))     deallocate(pvec_outi)

  if(allocated(cmat))       then
!$acc exit data delete(cmat)
    deallocate(cmat)
  endif

  call GEO_alloc(0)

  call cgyro_clean_implicit_gk

end subroutine cgyro_kernel

