!-----------------------------------------------------------------
! xgyro_mpi_setup.f90
!
! PURPOSE:
!  Split main MPI communicator in per-cgyro ones  
!-----------------------------------------------------------------

subroutine xgyro_mpi_setup

  use timer_lib
  use mpi

  use xgyro_globals
  use cgyro_globals, only : i_err, i_proc, n_proc, CGYRO_COMM_WORLD
  use xgyro_io

  implicit none

  integer :: i_group, splitkey
  character(len=192) :: msg
  integer :: i, tmpi

  ! Local group indices:

  if (xgyro_mpi_rank_order == 1) then
     i_group = 1
     tmpi = 0
     do i=1,xgyro_n_dirs
       tmpi = tmpi + xgyro_n_mpi(i)
       if (tmpi <= xgyro_i_proc) i_group = i_group + 1 ! we could break, but keep it simple
     enddo
     call xgyro_info('MPI rank alignment 1')
  else
     call xgyro_error('Unsupported MPI rank alignment, not 1')
     return
  endif
  xgyro_i_dir = i_group;

  !------------------------------------------------

  !-----------------------------------------------------------
  ! Split up GYRO_COMM_WORLD into groups and adjoint:
  !
  !             NEW_COMM_1  and  NEW_COMM_2
  !
  splitkey = xgyro_i_proc
  call MPI_COMM_SPLIT(XGYRO_COMM_WORLD,&
       i_group,& 
       splitkey,&
       CGYRO_COMM_WORLD, &
       i_err)
  if (i_err /= 0) then
     call xgyro_error('CGYRO_COMM_WORLD not created')
     return
  endif

  !
  ! Query cgyro rank and size
  !
  call MPI_COMM_RANK(CGYRO_COMM_WORLD,i_proc,i_err)
  call MPI_COMM_SIZE(CGYRO_COMM_WORLD,n_proc,i_err)

end subroutine xgyro_mpi_setup

