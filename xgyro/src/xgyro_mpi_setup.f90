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
  integer :: i, j, i_mpi
  integer, dimension(:), allocatable :: groups

  ! Local group indices:

  if (xgyro_mpi_rank_order == 1) then
     ! convert from incremental xgyro_n_mpi(i) to abs mapping in groups
     ! e.g.
     !  1 -> 1
     !  2 -> 2
     !       3
     !  1 -> 4
     !
     allocate(groups(xgyro_n_proc))
     i_mpi = 0
     do i=1,xgyro_n_dirs
       do j=1,xgyro_n_mpi(i)
         i_mpi = i_mpi + 1 ! keep it 1-based
         groups(i_mpi) = i
       enddo
     enddo

     ! then use the i_mpi -> group mapping
     i_group = groups(xgyro_i_proc+1)  !xgyro_i_proc is 0-based

     deallocate(groups)
     call xgyro_info('MPI rank alignment 1')
  else
     ! find the 
     ! e.g.
     !   1 -> 1 -
     !   2 -> 2 4
     !   1 -> 3 -
     allocate(groups(xgyro_n_dirs))
     do i=1,xgyro_n_dirs
        groups(i) = xgyro_n_mpi(i) ! make a copy of xgyro_n_mpi, so we save the orig values
     enddo
     ! grab the first i_mpi that was not claimed by lower i
     i_group = 0
     do i=0,xgyro_i_proc
       ! try to pick the next group
       i_group = i_group+1
       if (i_group>xgyro_n_dirs) i_group = 1 ! wraparound
       ! if the current group is not available, find the next avaialble one
       do while ( groups(i_group) == 0 )
               i_group = i_group+1
               if (i_group>xgyro_n_dirs) i_group = 1 ! wraparound
       enddo
       ! remember the choice
       groups(i_group) = groups(i_group) - 1
     enddo
     deallocate(groups)
     call xgyro_info('MPI rank alignment 2')
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

