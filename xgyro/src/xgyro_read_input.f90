subroutine xgyro_read_input

  use mpi
  use xgyro_globals
  use cgyro_globals

  implicit none

  character(len=40) :: msg
  integer :: i
  integer :: total_mpi

  if (xgyro_i_proc == 0) open(unit=1,file=trim(path)//'input.xgyro.gen',status='old')

  call xgyro_readbc_int(xgyro_mpi_rank_order,'MPI_RANK_ORDER')
  call xgyro_readbc_int(xgyro_n_dirs,'N_DIRS')
  allocate(xgyro_n_mpi(xgyro_n_dirs))
  allocate(xgyro_dir_name(xgyro_n_dirs))
  total_mpi = 0
  do i=1,xgyro_n_dirs
    write (msg, "(A,I0)") 'N_MPI_', i
    call xgyro_readbc_int(xgyro_n_mpi(i),msg)
    total_mpi = total_mpi + xgyro_n_mpi(i)
    write (msg, "(A,I0)") 'DIR_', i
    call xgyro_readbc_string(xgyro_dir_name(i),msg)
  enddo
  if (total_mpi /= xgyro_n_proc) then
          if (xgyro_i_proc == 0) write(*,*) 'ERROR: Invalid MPI rank size', total_mpi, ' != ',xgyro_n_proc
          call MPI_ABORT(XGYRO_COMM_WORLD,1,i_err)
          STOP 'Invalid MPI rank size'
  endif


  if (xgyro_i_proc == 0) close(1)

end subroutine xgyro_read_input

!------------------------------------------------------------
! Service routines: 
!
! (1) read and broadcast an integer:
!
subroutine xgyro_readbc_int(p,label)

  use mpi
  use xgyro_globals, only : xgyro_i_proc,XGYRO_COMM_WORLD
  use cgyro_globals, only : i_err

  implicit none
  integer, intent(inout) :: p
  character (len=*), intent(in) :: label
  character (len=40) :: actual_label

  if (xgyro_i_proc == 0) then
       read(1,*) p,actual_label
       if (label /= actual_label) then
          write(*,*) 'ERROR: Invalid label found in input file', actual_label, ' != ', label
          call MPI_ABORT(XGYRO_COMM_WORLD,1,i_err)
          STOP 'Invalid label found in input file'
       endif
  endif

  call MPI_BCAST(p,1,MPI_INTEGER,0,XGYRO_COMM_WORLD,i_err)

end subroutine xgyro_readbc_int
!
! (2) read and broadcast a real:
!
subroutine xgyro_readbc_real(x)
  
  use mpi
  use xgyro_globals, only : xgyro_i_proc,XGYRO_COMM_WORLD
  use cgyro_globals, only : i_err

  implicit none
  real, intent(inout) :: x

  if (xgyro_i_proc == 0) read(1,*) x

  call MPI_BCAST(x,1,MPI_DOUBLE_PRECISION,0,XGYRO_COMM_WORLD,i_err)

end subroutine xgyro_readbc_real

! (3) read and broadcast a string:
!
subroutine xgyro_readbc_string(p,label)

  use mpi
  use xgyro_globals, only : xgyro_i_proc,XGYRO_COMM_WORLD
  use cgyro_globals, only : i_err

  implicit none
  character (len=80), intent(inout) :: p
  character (len=*), intent(in) :: label
  character (len=40) :: actual_label

  if (xgyro_i_proc == 0) then
       read(1,*) p,actual_label
       if (label /= actual_label) then
          write(*,*) 'ERROR: Invalid label found in input file', actual_label, ' != ', label
          call MPI_ABORT(XGYRO_COMM_WORLD,1,i_err)
          STOP 'Invalid label found in input file'
       endif
  endif

  call MPI_BCAST(p,80,MPI_CHARACTER,0,XGYRO_COMM_WORLD,i_err)

end subroutine xgyro_readbc_string

!------------------------------------------------------------
