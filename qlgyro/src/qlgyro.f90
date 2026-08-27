!-----------------------------------------------------------------
! qlgyro.f90
!
! PURPOSE:
!  Main program wrapper for *standalone* QLGYRO usage.
!-----------------------------------------------------------------

program qlgyro

  use mpi
  use qlgyro_globals
  use qlgyro_cgyro_interface

  !-----------------------------------------------------------------
  implicit none
  !
  integer :: supported, n_omp
  integer, external :: omp_get_max_threads
  !-----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Query OpenMP for threads
  !
  n_omp = omp_get_max_threads()
  !-----------------------------------------------------------------
  
  !-----------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator, including support for 
  ! funneled threading (needed if OpenMP is enabled).
  !
  if (n_omp > 1) then
     call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,supported,ierr)
     if (supported < MPI_THREAD_FUNNELED) then
        call qlgyro_catch_error('ERROR: (GYRO) Multi-threaded MPI not supported.')
     endif
  else 
     call MPI_INIT_THREAD(MPI_THREAD_SINGLE,supported,ierr)
  endif
  !-----------------------------------------------------------------
  
  !-----------------------------------------------------------------
  ! Path is cwd:
  !
  path= './'
  !-----------------------------------------------------------------
  
  !-----------------------------------------------------------------
  ! Query MPI for dimensions
  !
  QLGYRO_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(QLGYRO_COMM_WORLD,i_proc_global,ierr)
  call MPI_COMM_SIZE(QLGYRO_COMM_WORLD,n_proc_global,ierr)
  n_inst = n_proc_global
  !-----------------------------------------------------------------
  
  !-----------------------------------------------------------------
  ! Standard standalone operation
  transport_method = 0
  !

  call qlgyro_allocate_globals
  
  call qlgyro_read_input

  call qlgyro_comm_setup

  if (code .eq. 0) then
     ! Set GYRO communicator to local communicator (per instance).
     !GYRO_COMM_WORLD = qlgyro_comm
  
     !call gyro_init(path, GYRO_COMM_WORLD)

     !call qlgyro_run_gyro
     print*, "GYRO not supported in QLGYRO"
     call exit

  else

     CGYRO_COMM_WORLD = qlgyro_comm
     call cgyro_init(path, CGYRO_COMM_WORLD)

     ! Map CGYRO variables to their interface values
     call map_global2interface
     call allocate_cgyro_interface

     if (cgyro_n_toroidal_in .gt. 1) then
        call qlgyro_run_cgyro_fluxtube
     else
        call qlgyro_run_cgyro_balloon
     end if

  end if
  ! Synchronise all cores
  call qlgyro_comm_sync
  
  ! Apply saturation rule
  if (sat_rule .gt. 0) then
     call qlgyro_sat1
  else if (sat_rule .eq. -1) then
     call qlgyro_sat_mg
  else
     open(unit=1,file=trim(runfile),position='append')
     write(1,*) '-----------------------------------------'
     write(1,*) 'Sat rule ', sat_rule, 'not a valid option'
     write(1,*) '-----------------------------------------'
     stop
  end if
  
  ! Sum flux over ky spectrum
  call qlgyro_sum_fluxes

  if (i_proc_global .eq. 0) then
     ! write ky spectrum to file out.tglf.ky_spectrum
     call write_qlgyro_ky_spectrum
     call write_qlgyro_QL_weight_spectrum
     call write_qlgyro_eigenvalue_spectrum
     call write_qlgyro_flux_spectrum
     call write_qlgyro_field_spectrum
     call write_qlgyro_gbflux
     call write_qlgyro_units
     call write_qlgyro_sat_geo_spectrum
     call write_qlgyro_kxrms_spectrum
  
     print*, 'Completed QLGYRO run'
     
  end if
  
  call MPI_FINALIZE(ierr)
end program qlgyro
