!------------------------------------------------
! gyro_write_restart.mpiio.f90
!
! PURPOSE:
!  This is the master file controlling output of
!  restart data.
!------------------------------------------------

subroutine gyro_write_restart

  use mpi
  use gyro_globals
  use gyro_pointers

  !----------------------------------------------
  implicit none
  !
  real :: t_current_old
  !
  integer :: io
  integer :: data_step_old
  integer :: n_proc_old
  integer :: i_restart_old
  character (len=1) :: i_tag
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

  !---------------------------------------------
  ! Restarting without writing new restart files
  !
  if (restart_method == 2) return
  !---------------------------------------------

  io = io_restart

  !------------------------------
  ! Determine which file to use:
  !
  i_restart = 1-i_restart 
  !------------------------------

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[Saving to ',trim(path)//file_restart(i_restart),']'
  endif

  !-----------------------------------------------
  ! Dump h and blending coefficients:
  !
  filemode = IOR(MPI_MODE_RDWR,MPI_MODE_CREATE)
  disp     = 0

  if (huge_restart == 0) then

     call MPI_INFO_CREATE(finfo,i_err)

     call MPI_FILE_OPEN(GYRO_COMM_WORLD,&
          trim(path)//file_restart(i_restart),&
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

     offset1 = size(h)*i_proc

     call MPI_FILE_WRITE_AT(fhv,&
          offset1,&
          h,&
          size(h),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)

     call MPI_FILE_SYNC(fhv,i_err)
     call MPI_FILE_CLOSE(fhv,i_err)

  else

     do is=1,n_kinetic
        i_tag = achar(is-1+iachar("0"))

        call MPI_INFO_CREATE(finfo,i_err)

        call MPI_FILE_OPEN(GYRO_COMM_WORLD,&
             trim(path)//file_restart(i_restart)//i_tag,&
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

        offset1 = size(h(:,:,:,is))*i_proc

        call MPI_FILE_WRITE_AT(fhv,&
             offset1,&
             h(:,:,:,is),&
             size(h(:,:,:,is)),&
             MPI_COMPLEX16,&
             fstatus,&
             i_err)

        call MPI_FILE_SYNC(fhv,i_err)
        call MPI_FILE_CLOSE(fhv,i_err)

     enddo

  endif

  !---------------------------------------------------------
  ! Dump restart parameters
  !
  if (i_proc == 0) then

     if (t_current > restart_data_skip*time_skip*dt*1.1) then

        open(unit=io,file=trim(path)//file_tag_restart(1),status='old')
        read(io,*) data_step_old
        read(io,fmtstr) t_current_old
        read(io,*) n_proc_old
        read(io,*) i_restart_old
        close(io)

        open(unit=io,file=trim(path)//file_tag_restart(0),status='replace')
        write(io,*) data_step_old
        write(io,fmtstr) t_current_old
        write(io,*) n_proc_old
        write(io,*) i_restart_old
        close(io)

     endif

     open(unit=io,file=trim(path)//file_tag_restart(1),status='replace')
     write(io,*) data_step
     write(io,fmtstr) t_current
     write(io,*) n_proc
     write(io,*) i_restart
     close(io)

  endif


end subroutine gyro_write_restart
