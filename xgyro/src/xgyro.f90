program cgyro

  use mpi
  use xgyro_globals
  use cgyro_globals
  use xgyro_io
  use timer_lib

  implicit none
  
  integer :: supported
  integer, external :: omp_get_max_threads
  character(len=32) :: arg

  !----------------------------------------------------------------
  ! Find value of test flag
  call get_command_argument(1,arg)
  if (trim(arg) == '0') then
     test_flag = 0
  else
     test_flag = 1
  endif
  !----------------------------------------------------------------
  
  ! Base path is cwd:
  path= './'
  ! create xgyro info file ASAP, so we can report error
  open(unit=io,file=trim(path)//xgyro_runfile_info,status='replace')

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
        call xgyro_error('Multi-threaded MPI not supported.')
     endif
  else 
     call MPI_INIT_THREAD(MPI_THREAD_SINGLE,supported,i_err)
  endif
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Set the world MPI communicator
  !
  XGYRO_COMM_WORLD = MPI_COMM_WORLD
  !
  ! Query rank and size
  !
  call MPI_COMM_RANK(XGYRO_COMM_WORLD,xgyro_i_proc,i_err)
  call MPI_COMM_SIZE(XGYRO_COMM_WORLD,xgyro_n_proc,i_err)
  !write(*,*) "MPI size", xgyro_i_proc,xgyro_n_proc,i_err
  !-----------------------------------------------------------------

  ! read the xgyro input
  call timer_lib_init('input')
  call timer_lib_in('input')
  call xgyro_read_input
  call timer_lib_out('input')
  if (error_status /= 0) then
    write(*,*) "ERROR while reading input: ", error_message
    call MPI_ABORT(XGYRO_COMM_WORLD,1,i_err)
    call MPI_FINALIZE(i_err)
  endif

  ! 
  call xgyro_mpi_setup
  if (error_status /= 0) then
    write(*,*) "ERROR in initial MPI setup: ", error_message
    call MPI_ABORT(XGYRO_COMM_WORLD,1,i_err)
    call MPI_FINALIZE(i_err)
  endif

  ! my CGYRO path is from read params
  path = trim(xgyro_dir_name(xgyro_i_dir))//'/'
  !write(*,*) xgyro_i_proc,i_proc,n_proc,xgyro_i_dir,path

  ! now read the cgyro input
  !call timer_lib_init('input')
  !call timer_lib_in('input')
  !call cgyro_read_input
  !call timer_lib_out('input')

!
!  call cgyro_init_kernel
!  if (error_status == 0) then
!        call cgyro_kernel
!        call cgyro_final_kernel
!  endif
  call MPI_FINALIZE(i_err)

end program cgyro
