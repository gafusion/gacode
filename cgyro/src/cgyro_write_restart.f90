!------------------------------------------------
! cgyro_write_restart.f90
!
! PURPOSE:
!  This is the master file controlling output of
!  restart data using MPI-IO.
!------------------------------------------------

subroutine cgyro_write_restart

  use mpi
  use cgyro_globals
  use cgyro_io

  !----------------------------------------------
  implicit none
  !
  ! Required for MPI-IO: 
  !
  integer :: filemode
  integer :: finfo
  integer :: fhv
  integer :: fstatus(MPI_STATUS_SIZE)
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer(kind=MPI_OFFSET_KIND) :: offset1
  !
  ! File chunking variables
  !
  integer :: i1,i2,j,n_loc,n_remain
  complex, dimension(:,:), allocatable :: h_chunk
  !----------------------------------------------

  ! Print this data on restart steps only; otherwise exit now
  if (mod(i_time,restart_step*print_step) /= 0) return

  !-----------------------------------------------
  ! Dump h and blending coefficients:
  !
  filemode = IOR(MPI_MODE_RDWR,MPI_MODE_CREATE)
  disp     = 0

  n_loc = nc/n_chunk
  n_remain = nc-n_loc*n_chunk
  allocate(h_chunk(n_loc+n_remain,nv_loc))

  offset1 = size(h_chunk)*i_proc

  do j=1,n_chunk

     call MPI_INFO_CREATE(finfo,i_err)

     call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart//rtag(j),&
          filemode,&
          finfo,&
          fhv,&
          i_err)

     call MPI_FILE_SET_VIEW(fhv,&
          disp,&
          MPI_COMPLEX16,&
          MPI_COMPLEX16,&
          'native',&
          finfo,&
          i_err)

     i1 = 1+(j-1)*n_loc
     i2 = j*n_loc
     if (j == n_chunk) then
        i2 = nc
        n_loc = n_loc+n_remain
     endif

     h_chunk(1:n_loc,:) = h_x(i1:i2,:)

     call MPI_FILE_WRITE_AT(fhv,&
          offset1,&
          h_chunk,&
          size(h_chunk),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)

     call MPI_FILE_SYNC(fhv,i_err)
     call MPI_FILE_CLOSE(fhv,i_err)

  enddo

  deallocate(h_chunk)

  !---------------------------------------------------------
  ! Dump restart parameters to ASCII file.
  !
  if (i_proc == 0) then

     open(unit=io,file=trim(path)//runfile_restart_tag,status='replace')
     write(io,*) i_current
     write(io,fmtstr) t_current
     close(io)

  endif
  !---------------------------------------------------------

end subroutine cgyro_write_restart
