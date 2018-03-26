program cgyro

  use mpi
  use cgyro_globals, only : path, CGYRO_COMM_WORLD, i_proc, n_proc, i_err, n_omp
  use cgyro_io
  use timer_lib

  implicit none
  integer :: supported
  integer, external :: omp_get_max_threads

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
     call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,supported,i_err)
     if (supported < MPI_THREAD_FUNNELED) then
        write (*,*) "ERROR: Multi-threaded MPI not supported." 
        call cgyro_error('Multi-threaded MPI not supported.')
     endif
  else 
     call MPI_INIT_THREAD(MPI_THREAD_SINGLE,supported,i_err)
  endif
  !-----------------------------------------------------------------

  ! Path is cwd:
  path= './'

  !-----------------------------------------------------------------
  ! Set the world MPI communicator
  !
  CGYRO_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(CGYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(CGYRO_COMM_WORLD,n_proc,i_err)
  !-----------------------------------------------------------------

  call timer_lib_init('input')
  call timer_lib_in('input')
  call cgyro_read_input
  call timer_lib_out('input')
  
  call cgyro_kernel

  call MPI_FINALIZE(i_err)

end program cgyro
