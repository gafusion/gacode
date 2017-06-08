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
  integer :: i1,i2,j
  character (len=1) :: i_tag
  integer(kind=MPI_OFFSET_KIND) :: disp
  integer(kind=MPI_OFFSET_KIND) :: offset1
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

  allocate(h_chunk(nc/n_chunk,nv_loc))

  do j=1,n_chunk

     if (n_chunk > 1) then
        i_tag = achar(j+iachar("0"))
     else
        i_tag = ''
     endif

     call MPI_INFO_CREATE(finfo,i_err)

     call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart//i_tag,&
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

     offset1 = size(h_chunk)*i_proc

     call MPI_FILE_READ_AT(fhv,&
          offset1,&
          h_chunk,&
          size(h_chunk),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)

     i1 = 1+(j-1)*(nc/n_chunk)
     i2 = j*(nc/n_chunk)
     if (j == n_chunk) i2 = nc

     h_x(i1:i2,:) = h_chunk

     call MPI_FILE_CLOSE(fhv,i_err)

  enddo

  deallocate(h_chunk)

end subroutine cgyro_read_restart
