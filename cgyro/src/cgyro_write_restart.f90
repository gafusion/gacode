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

  ! Print this data on restart steps only; otherwise exit now
  if (mod(i_time,restart_step*print_step) /= 0) return

  ! Write to single or multiple files
  if (mpiio_num_files == 1) then
     call cgyro_write_restart_one
  else
     call cgyro_write_restart_many
  endif

  ! Write restart tag
  if (i_proc == 0) then
     open(unit=io,file=trim(path)//runfile_restart_tag,status='replace')
     write(io,*) i_current
     write(io,fmtstr) t_current
     close(io)
  endif

end subroutine cgyro_write_restart

subroutine cgyro_write_restart_one

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
  !----------------------------------------------

  !-----------------------------------------------
  ! Dump h and blending coefficients:
  !
  filemode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
  disp     = 0

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*i_proc
  if (offset1<0) then
     call cgyro_error('ERROR: (CGYRO) overflow detected in cgyro_write_restart_one')
     return
  endif

  ! TODO Error handling
  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_INFO_SET(finfo,"striping_factor",mpiio_stripe_str,i_err)

  call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart,&
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

  call MPI_FILE_WRITE_AT(fhv,&
          offset1,&
          h_x,&
          size(h_x),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)

  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_WRITE_AT in cgyro_write_restart_one failed')
     return
  endif

  call MPI_FILE_SYNC(fhv,i_err)
  call MPI_FILE_CLOSE(fhv,i_err)
  call MPI_INFO_FREE(finfo,i_err)

end subroutine cgyro_write_restart_one

subroutine cgyro_write_restart_many

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
  !----------------------------------------------

  !-----------------------------------------------
  ! Dump h and blending coefficients:
  !
  filemode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
  disp     = 0

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*i_proc_restart_io
  if (offset1 < 0) then
     call cgyro_error('ERROR: (CGYRO) overflow detected in cgyro_write_restart_v2')
     return
  endif

  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_INFO_SET(finfo,"striping_factor",mpiio_stripe_str,i_err)

  call MPI_FILE_OPEN(NEW_COMM_RESTART_IO,&
          trim(path)//runfile_restart//rtag(i_group_restart_io),&
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

  call MPI_FILE_WRITE_AT(fhv,&
          offset1,&
          h_x,&
          size(h_x),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)

  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_WRITE_AT in cgyro_write_restart_v2 failed')
     return
  endif

  call MPI_FILE_SYNC(fhv,i_err)
  call MPI_FILE_CLOSE(fhv,i_err)
  call MPI_INFO_FREE(finfo,i_err)

  ! Wait for all groups to finish
  call MPI_Barrier(CGYRO_COMM_WORLD,i_err)

end subroutine cgyro_write_restart_many
