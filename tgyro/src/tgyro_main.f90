program tgyro_main

  use mpi
  use tgyro_globals

  !-----------------------------------------------------------------
  implicit none
  !
  integer :: supported
  integer :: i_omp,n_omp
  integer, external :: omp_get_max_threads, omp_get_thread_num

  !-----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Query OpenMP for dimensions
  !
  i_omp = omp_get_thread_num()
  n_omp = omp_get_max_threads()
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Initialize MPI_COMM_WORLD communicator, including support for 
  ! funneled threading (needed if OpenMP is enabled).
  !
  if (n_omp > 1) then
     call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,supported,ierr)
     if (supported < MPI_THREAD_FUNNELED) then
        call tgyro_catch_error('ERROR: (TGYRO) Multi-threaded MPI not supported.')
     endif
  else 
     call MPI_INIT_THREAD(MPI_THREAD_SINGLE,supported,ierr)
  endif
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  call MPI_COMM_RANK(MPI_COMM_WORLD,i_proc_global,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc_global,ierr)
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Get input paths and broadcast them
  !
  call tgyro_read_input   
  !
  ! At this stage, paths/procs information is only on process 0.
  !
  call tgyro_comm_setup
  !-----------------------------------------------------------------

  select case (tgyro_mode)

  case (1)

     ! Local transport model

     call tgyro_iteration_driver 

     if (i_proc_global == 0) then
        open(unit=1,file=trim(runfile),position='append')
        write(1,*) error_msg
        close(1)
     endif

  case (3)

     ! Multi-job utility

     call tgyro_multi_driver

  case (4)

     call tgyro_cgyro_iterate
     
  case default

     call tgyro_catch_error('ERROR: (TGYRO) Bad value for tgyro_mode')

  end select

  call MPI_FINALIZE(ierr)

end program tgyro_main
