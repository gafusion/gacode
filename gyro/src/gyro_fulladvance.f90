!--------------------------------------------------------
! gyro_fulladvance.f90
!
! PURPOSE:
!  Control routine for a single timestep of length dt.
!--------------------------------------------------------

subroutine gyro_fulladvance

  use gyro_globals
  use gyro_pointers

  !----------------------------------------
  implicit none
  !----------------------------------------

  h_old(:,:,:,:) = h(:,:,:,:)

  h_cap_old2(:,:,:,:) = h_cap_old(:,:,:,:)
  h_cap_old(:,:,:,:)  = h_cap(:,:,:,:)

  field_blend_old2(:,:,:) = field_blend_old(:,:,:)
  field_blend_old(:,:,:)  = field_blend(:,:,:)

  field_tau_old2(:,:,:,:) = field_tau_old(:,:,:,:)
  field_tau_old(:,:,:,:)  = field_tau(:,:,:,:)

  gyro_uv_old2(:,:,:,:,:) = gyro_uv_old(:,:,:,:,:)
  gyro_uv_old(:,:,:,:,:)  = gyro_uv(:,:,:,:,:)

  !----------------------
  ! Start time step cycle
  !----------------------

  !------------------------------------------------------
  ! MANAGE pitch-angle scattering:
  !
  if (collision_flag == 1) call gyro_collision_main
  !------------------------------------------------------

  !------------------------------------------------------
  ! MANAGE time-integrator
  !
  select case (electron_method) 

  case (1,3,4)

     ! All species explicit

     call gyro_timestep_explicit

  case (2)

     ! Drift-kinetic electrons

     call gyro_timestep_implicit

  end select
  !------------------------------------------------------

  call gyro_timer_in('Diagnos.-allstep')

  call gyro_timestep_error

  !------------------------------------------------------
  ! MANAGE diagnostics
  !
  ! * Note that (2,3,4) are called once after making 
  !   initial conditions (gyro_initial_condition).
  !
  ! 1. Compute flux-surface average of (phi,a)
  ! 
  call gyro_field_fluxave
  !
  ! 2. Compute (phi,a) for plotting; we want to call 
  !    gyro_field_plot every timestep because it is 
  !    averaged (with a weight) over every timestep.
  !
  !    To compute E_parallel
  !
  call gyro_field_time_derivative
  !
  !    then 
  !
  call gyro_field_plot
  !
  ! 3. Compute (n,T) moments for plotting if user has 
  !    selected OUTPUT_METHOD > 1.
  !
  call gyro_moments_plot 
  !------------------------------------------------------

  !------------------------------------------------------
  ! MANAGE time:
  !
  ! Increment simulation time; this is the actual time 
  ! in the gyrokinetic equation:
  !  
  t_current = t_current+dt
  !------------------------------------------------------

  !-------------------------------------------------------------------
  ! MANAGE data output: 
  !
  call gyro_timer_out('Diagnos.-allstep')
  call gyro_timer_in('Diagnos.-datastep')


  if (modulo(step,time_skip) == 0) then

     ! Counter for number of data output events.
     if (modulo(step,time_skip) == 0) then
        data_step = data_step+1
     endif

     ! Main data I/O handler

     io_control = 2*output_flag

     if (modulo(step,time_skip) == 0) then
        call gyro_write_timedata
     endif

     !--------------------------------------------------
     ! Update diffusivity and flux time-record for TGYRO 
     !--------------------------------------------------

     ! Restart test:

     if (modulo(data_step,restart_data_skip) == 0 &
          .and. restart_method >= 0) then

        call gyro_write_restart

     endif

  endif

  call gyro_timer_out('Diagnos.-datastep')
  !-------------------------------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print '(t2,2(a,i5,3x))','-> step =',step,'data_step =',data_step
     print *,'**[gyro_fulladvance done]'
  endif

end subroutine gyro_fulladvance
