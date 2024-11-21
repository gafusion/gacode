!-----------------------------------------------------------------
! xgyro_globals.f90
!
! PURPOSE:
!  XGYRO global variables.  The idea is to have a primary, large
!  module containing all essential XGYRO arrays and scalars.
!-----------------------------------------------------------------

module xgyro_globals

  !---------------------------------------------------------------
  ! I/O and error management variables
  !
  character(len=80) :: xgyro_path
  character(len=14) :: xgyro_runfile_info    = 'out.xgyro.info'

  !---------------------------------------------------------------
  ! Input parameters:
  !
  integer :: xgyro_mpi_rank_order
  ! all the arrays are xgyro_n_dirs in size
  integer :: xgyro_n_dirs
  integer, dimension(:), allocatable :: xgyro_n_mpi
  ! note: dir_name has same length as path in cgyro_globals
  character(len=80), dimension(:), allocatable :: xgyro_dir_name

  !---------------------------------------------------------------
  ! MPI variables and pointers
  ! 
  integer :: xgyro_i_proc
  integer :: xgyro_n_proc
  integer :: xgyro_i_dir
  integer :: XGYRO_COMM_WORLD

end module xgyro_globals
