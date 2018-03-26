!--------------------------------------------------------------
! cgyro_init.f90
!
! PURPOSE:
!  Initialize external CGYRO interface.
!---------------------------------------------------------------

subroutine cgyro_init(path_in,mpi_comm_in)

  use mpi
  use cgyro_globals
  use timer_lib
 
  implicit none

  ! Input parameters (IN) - REQUIRED
  character(len=*), intent(in) :: path_in
  integer, intent(in) :: mpi_comm_in

  ! Local variable
  logical :: inputdat_flag

  integer, external :: omp_get_max_threads

  ! Set appropriate global variables
  path = path_in
  CGYRO_COMM_WORLD = mpi_comm_in

  !----------------------------------------------------------------
  ! Query OpenMP for dimensions
  !
  n_omp = omp_get_max_threads()
  !-----------------------------------------------------------------

  !----------------------------------------------------------------
  ! Query MPI for dimensions
  !
  call MPI_COMM_RANK(CGYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(CGYRO_COMM_WORLD,n_proc,i_err)
  !-----------------------------------------------------------------

  ! Check if input.cgyro.gen file exists and set:
  !
  !   inputdat_flag=TRUE  if input.cgyro.gen does exist, 
  !                 FALSE if input.cgyro.gen does NOT exist
  !
  if (i_proc == 0) then
     inquire(file=trim(path)//'input.cgyro.gen',exist=inputdat_flag)
  endif

  call MPI_BCAST(inputdat_flag,1,MPI_LOGICAL,0,CGYRO_COMM_WORLD,i_err)

  if (inputdat_flag .eqv. .false.) then
     call MPI_FINALIZE(i_err)
     stop
  endif

  call timer_lib_init('input')
  call timer_lib_in('input')
  call cgyro_read_input
  call timer_lib_out('input')

end subroutine cgyro_init
