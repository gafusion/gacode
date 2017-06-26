!------------------------------------------------------
! cgyro_read_restart.f90
!
! PURPOSE:
!  This is the master file controlling the restart
!  via MPI-IO.
!------------------------------------------------------

subroutine cgyro_read_restart

  use mpi
  use cgyro_globals
  use cgyro_io

  !---------------------------------------------------
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
  !---------------------------------------------------

  !---------------------------------------------------------
  ! Read restart parameters from ASCII file.
  !
  if (i_proc == 0) then

     open(unit=io,&
          file=trim(path)//runfile_restart_tag,&
          status='old')

     read(io,*) i_current
     read(io,fmtstr) t_current
     close(io)

  endif

  ! Broadcast to all cores.

  call MPI_BCAST(i_current,&
       1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)

  call MPI_BCAST(t_current,&
       1,MPI_DOUBLE_PRECISION,0,CGYRO_COMM_WORLD,i_err)
  !---------------------------------------------------------

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

     call MPI_FILE_READ_AT(fhv,&
          offset1,&
          h_chunk,&
          size(h_chunk),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)

     i1 = 1+(j-1)*n_loc
     i2 = j*n_loc
     if (j == n_chunk) then
        i2 = nc
        n_loc = n_loc+n_remain
     endif

     h_x(i1:i2,:) = h_chunk(1:n_loc,:)

     call MPI_FILE_CLOSE(fhv,i_err)

  enddo

  deallocate(h_chunk)

end subroutine cgyro_read_restart
