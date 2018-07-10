!-----------------------------------------------------------------
! gyro_do.f90
!
! PURPOSE:
!  Subroutinized main gyro program. 
!-----------------------------------------------------------------

subroutine gyro_do

  use mpi
  use gyro_globals
  use gyro_pointers
  use gyro_interface
  use math_constants
  use geo

  !--------------------------------------
  implicit none
  !
  logical :: rfe
  !--------------------------------------

  ! Begin with clean exit status
  !
  ! gyro_exit_status=0 : success
  ! gyro_exit_status=1 : failure
  ! gyro_exit_status=2 : warning
  !
  gyro_exit_status  = 0
  gyro_exit_message = 'unset'

  ! Prepend path:
  runfile  = trim(path)//trim(baserunfile)
  precfile = trim(path)//trim(baseprecfile)

  if (baserunfile == 'out.gyro.run')  then
     if (i_proc==0 .and. output_flag==1) THEN
        inquire(file=trim(runfile),exist=rfe)
        if (.not.rfe) then
           open(unit=1,file=trim(runfile),status='unknown')
           close(1)
        endif
     endif
  endif

  ! Dump interface variables to record exact state of GYRO.
  call map_global2interface
  call interfacelocaldump

  !--------------------------------------------------------------
  ! Early initializations:
  !
  total_memory  = 0.0
  alltime_index = 0
  cpu_maxindx   = 0
  cpu           = -1.0
  cpu_in        = 0.0
  !
  ! TIMER NOTES: 
  ! - print order follows init. order below,
  ! - strings will be split on '-' character.
  !
  call gyro_timer_init('Field-interp.a')
  call gyro_timer_init('Field-interp.b')
  call gyro_timer_init('Velocity-sum')
  call gyro_timer_init('Field-explicit')
  call gyro_timer_init('Field-implicit')
  call gyro_timer_init('Gyroave-h')
  call gyro_timer_init('Implicit-he')
  call gyro_timer_init('RHS-total')
  call gyro_timer_init('Coll.-step')
  call gyro_timer_init('Coll.-comm')
  call gyro_timer_init('Nonlinear-step')
  call gyro_timer_init('Nonlinear-comm')
  call gyro_timer_init('Diagnos.-allstep')
  call gyro_timer_init('Diagnos.-datastep')
  call gyro_timer_init('Full-step')
  !--------------------------------------------------------------

  !----------------------------------------------------------------
  ! If running in eigensolve mode.
  !
  if (linsolve_method == 2) then
     if (nonlinear_flag == 1) then
        if (i_proc==0 .and. gkeigen_j_set==0) then
           print '(a)','INFO: (gyro) Eigensolver unavailable in nonlinear mode.'
        endif
        stop 
     else
        if (i_proc==0 .and. gkeigen_j_set==0) then
           print '(a)', 'INFO: (gyro) GYRO is running in eigensolve mode.'
        endif
        eigensolve_restart_flag = restart_method
        restart_method = 0
        if (electron_method /= 1) then
           electron_method = 4
        endif
     endif
  endif
  !----------------------------------------------------------------

  !------------------------------------------------------------
  ! The order of these routines is critical:
  !
  ! Startup timer:
  startup_time = MPI_Wtime()
  !
  ! Sort through and check all combinations of operational modes
  !
  call gyro_select_methods
  !
  ! Set parameters connected with timestepping.
  !
  call gyro_initialize_timestep
  !
  ! Generate poloidal (theta) grid dimensions (no operators yet).
  !
  call gyro_theta_grid
  !
  ! Parallel setup 
  !
  if (gyrotest_flag == 0) then
     call gyro_mpi_grid
  else
     n_n_1     = 1
     i_group_1 = 0
  endif
  !
  ! Read, generate or otherwise construct equilibrium profiles.  If experimental 
  ! profiles are used, GEO will be allocate/deallocated with all settings 
  ! determined in EXPRO.
  !
  call gyro_alloc_profile_sim(1)
  call gyro_profile_init
  !
  ! Construct pitch-angle (lambda) and energy grids
  ! and integration weights using Gauss-Legendre rules:
  !
  call gyro_alloc_velocity(1)
  call gyro_alloc_orbit(1)
  !
  ! Set geometry (GEO) library control variables
  !
  GEO_nfourier_in = n_fourier_geo 
  GEO_model_in    = geometry_method
  GEO_signb_in    = -btccw
  !
  ! Lambda (pitch-angle) weights (GEO needed again, so just reallocate)
  call gyro_lambda_grid

  !
  ! Energy weights
  !
  call energy_integral(n_energy,energy_max,energy,w_energy)
  !
  ! Compute myriad arrays (for example, map between poloidal
  ! angle and orbit time) for poloidal discretization.
  !
  call gyro_banana_operators
  !
  ! Total velocity-space-tau weights.
  !
  call gyro_set_phase_space(trim(path)//'out.gyro.phase_space',1)
  !
  do i=1,n_x
     call gyro_to_geo(i)
  enddo
  !
  ! Generate geometry-dependent factors using model or Miller equilibrium:
  !
  call gyro_geometry_arrays
  !
  if (gyrotest_flag == 0) then

     call gyro_set_pointer_dim
     call gyro_alloc_distrib(1)

     ! Make pointers for use with parallelization scheme.
     call gyro_set_pointers

     ! Compute drift and diamagnetic frequency coefficients
     ! for GKE solution.
     call gyro_omegas

     ! Generate discrete Fourier and Fourier-sine transform 
     ! matrices for spectral Poisson and adaptive source methods, 
     ! compute radial derivative and gyroaverage operators,
     ! and then print them.

     call gyro_radial_operators
     !
     ! Precomputation of arrays which depend on blending 
     ! coefficients.  These are used in the Maxwell solves.
     call gyro_set_blend_arrays

     if (electron_method == 2) then

        ! Make advection operators for electrons
        call gyro_make_implicit_advect(0)

     endif
     !
     ! Collision setup (if required; note PNL hook)
     !
     if (collision_flag == 1) then
        call gyro_collision_setup
     endif

     !------------------------------------------------------
     ! FIELD MATRICES:
     !
     ! Manage building of field-solve matrices.  Must have 
     ! built the gyro-operators beforehand:
     !
     ! Explicit (Poisson,Ampere)
     !
     select case (n_field)
     case (1)
        call gyro_make_poisson_matrix
     case (2)
        call gyro_make_poisson_matrix
        call gyro_make_ampere_matrix
     case (3)
        call gyro_make_poissonaperp_matrix
        call gyro_make_ampere_matrix
     end select

     ! Implicit (Maxwell)
     !
     if (electron_method == 2) then
        call gyro_make_maxwell_matrix
     endif
     !------------------------------------------------------

     !------------------------------------------------------
     ! Make a few trivial quantities for the nonlinear step:
     !
     if (nonlinear_flag == 1) then
        call gyro_alloc_nl(1)
        call gyro_nl_setup
     endif
     !
     ! Write information about radial stencils:
     !
     call gyro_write_radial_op(trim(path)//'out.gyro.radial_op',1)
     !
     ! Allocate most large arrays:
     !
     call gyro_alloc_big(1)
     call gyro_initialize_arrays
     !
     ! Make attempt to estimate storage:
     !
     call gyro_memory_usage(trim(path)//'out.gyro.memory',30)
     !
     !------------------------------------------------------------

     !------------------------------------------------------------
     ! Open files and write values for t=0:
     !
     call gyro_read_restart
     !------------------------------------------------------------

     !---------------------------------------
     ! Write precision data:
     !
     call gyro_write_precision(10,sum(abs(gbflux)))
     !---------------------------------------

  endif

  startup_time = MPI_Wtime()-startup_time

  !---------------------------------------
  ! Dump input parameters runfile
  !
  call gyro_write_input
  !---------------------------------------

  !---------------------------------------------------------------
  ! I/O control for time-independent initial data
  !
  call gyro_write_initdata(&
       trim(path)//'out.gyro.profile',&
       trim(path)//'out.gyro.units',&
       trim(path)//'out.gyro.geometry_arrays',1)
  !
  !---------------------------------------------------------------

  !------------------------------------------------------------
  ! By this point, we are finished with running the
  ! 'gyro -t' (test) mode.
  !
  if (gyrotest_flag == 1) then
     call gyro_write_efficiency(trim(path)//'out.gyro.efficiency',1)
     call gyro_set_exit_status('test complete',0)
     return
  endif
  !------------------------------------------------------------

  if (restart_method == 0 .or. restart_method == 3) then
     ! Open (new or init)
     io_control = output_flag*1
  else
     ! Rewind (cont)
     io_control = output_flag*3
  endif
  if (gkeigen_j_set == 0) then
     call gyro_write_timedata
  endif

  !-------------------------------------------------
  ! NEW SIMULATION ONLY: write *initial conditions*
  !

  if (restart_method /= 1) then
     ! Write to output files.
     io_control = output_flag*2
     if (gkeigen_j_set == 0) then
        call gyro_write_timedata
     endif
  endif
  !--------------------------------------------

  select case (linsolve_method)

  case (1)

     !------------------------------
     ! INITIAL VALUE SOLVER
     !------------------------------

     !-----------------------------------------------
     ! start time-stepping

     do step=1,nstep

        call gyro_timer_in('Full-step')
        call gyro_fulladvance
        call gyro_timer_out('Full-step')

        !----------------------------------------------------
        ! Check for premature exit conditions
        ! 
        call catch_blowup
        call catch_halt_signal
        if (gyro_exit_status > 0) return
        !
        ! Convergence check for single-n simulation:
        ! freq_err calculated in gyro_write_freq
        !  
        if (freq_err < freq_tol) then
           call gyro_set_exit_status('converged',0)
           return
        endif
        !----------------------------------------------------

     enddo
     ! end time-stepping   
     !-----------------------------------------------

  case (2)

     !------------------------------
     ! GK EIGENVALUE SOLVER
     !------------------------------

     ! Determine true range of p_nek_loc for each processor
     n_nek_loc_true = n_nek_1/n_proc_1
     if (i_proc_1 < mod(n_nek_1,n_proc_1) ) then
        n_nek_loc_true = n_nek_loc_true + 1
     endif
     ! Length of total state vector
     h_length = n_kinetic * n_nek_1 * n_x * n_stack
     h_length_loc = n_kinetic * n_nek_loc_true * n_x * n_stack
     h_width_loc = h_length / gkeigen_proc_mult
     h_length_block = h_length_loc / gkeigen_proc_mult
     h_length_block_t = h_width_loc / n_proc
     seq_length = h_length_block * h_width_loc
     seq_length_t = h_length_block_t * h_length_loc

     call GKEIGEN_do                                       

  case (3)

     !------------------------------
     ! FIELD EIGENVALUE SOLVER
     !------------------------------

     call gyro_fieldeigen

  end select

end subroutine gyro_do
