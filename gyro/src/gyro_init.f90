!--------------------------------------------------------------
! gyro_init.f90
!
! PURPOSE:
!  Initialize external GYRO interface.
!---------------------------------------------------------------

subroutine gyro_init(path_in, mpi_comm_in)

  use mpi
  use gyro_globals
  use gyro_interface
  use ompdata

  implicit none

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in
  integer, intent(in) :: mpi_comm_in

  ! Local variable
  logical :: inputdat_flag

  integer, external :: omp_get_max_threads

  ! Set appropriate global variables
  path = path_in
  GYRO_COMM_WORLD = mpi_comm_in

  !----------------------------------------------------------------
  ! Query OpenMP for dimensions
  !
  n_omp = omp_get_max_threads()
  !-----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Query MPI for dimensions
  !
  call MPI_COMM_RANK(GYRO_COMM_WORLD, i_proc, i_err)
  call MPI_COMM_SIZE(GYRO_COMM_WORLD, n_proc, i_err)
  !-----------------------------------------------------------------

  ! Check if input.gyro.gen file exists and set:
  !
  !   inputdat_flag=TRUE  if input.gyro.gen does exist, 
  !                 FALSE if input.gyro.gen does NOT exist
  !
  if (i_proc == 0) then
     inquire(file=trim(path)//'input.gyro.gen',exist=inputdat_flag)
  endif

  call MPI_BCAST(inputdat_flag, 1, MPI_LOGICAL, 0, GYRO_COMM_WORLD, i_err)

  if (i_proc == 0) then
     open(unit=1,file=trim(path)//trim(baserunfile),position='append')
     if (inputdat_flag .eqv. .true.) then
        write(1,'(a,a,a)') 'INFO: (GYRO) gyro_init reading ',&
             trim(path),'input.gyro.gen'
     else
        write(1,'(a,a,a)') 'INFO: (GYRO) gyro_init NOT reading ',&
             trim(path),'input.gyro.gen'
     endif
     close(1)
  endif

  if (inputdat_flag .eqv. .true.) then

     ! Only call read_input subroutine if input.gyro.gen file exists

     call gyro_read_input

     ! initialize OpenMP parameters including block indices
     call gyro_init_ompdata()

     ! Map GLOBAL variables -> INTERFACE parameters
     call map_global2interface

  else

     ! If input.gyro.gen file does NOT exist then initialize arrays
     ! NOTE: Vector initialization is normally done by read_input subroutine
     mu_vec(0:10)           = 0.0
     dlnndr_vec(0:10)       = 0.0
     dlntdr_vec(0:10)       = 0.0
     sdlnndr_vec(0:10)      = 0.0
     sdlntdr_vec(0:10)      = 0.0
     n_vec(0:10)            = 0.0
     t_vec(0:10)            = 0.0
     eps_dlnndr_vec(0:10)   = 0.0
     eps_dlntdr_vec(0:10)   = 0.0
     z_vec(0:10)            = 0.0
     orbit_upwind_vec(0:10) = 0.0

     ! Default parameter values explicitly specified by read_input subroutine
     ! and NOT definable by user
     mu_vec(1)      =  1.0
     n_vec(0)       =  1.0
     t_vec(0)       =  1.0
     z_vec(0)       = -1.0

  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *, '[gyro_init done]'
  endif

end subroutine gyro_init
