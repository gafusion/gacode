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

  implicit none

  !---------------------------------------------------------
  ! Read restart parameters from ASCII file.
  !
  if (restart_flag == 1) then
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

  endif

  call cgyro_read_restart_one

end subroutine cgyro_read_restart

subroutine cgyro_read_restart_one

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
  !---------------------------------------------------

  character(8)  :: sdate
  character(10) :: stime
  character(5)  :: szone
  integer(KIND=8) :: start_time,cp_time
  integer(KIND=8) :: count_rate, count_max
  real :: cp_dt
  integer :: statusfd

  ! use system_clock to be consistent with cgyro_kernel
  call system_clock(start_time,count_rate,count_max)

  filemode = MPI_MODE_RDONLY
  disp     = 0

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*i_proc
  if (offset1 < 0) then
     call cgyro_error('ERROR: (CGYRO) overflow detected in cgyro_read_restart_one')
     return
  endif

  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart,&
          filemode,&
          finfo,&
          fhv,&
          i_err)

  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_OPEN in cgyro_read_restart_one failed')
     return
  endif

  call MPI_FILE_SET_VIEW(fhv,&
          disp,&
          MPI_COMPLEX16,&
          MPI_COMPLEX16,&
          'native',&
          finfo,&
          i_err)

  call MPI_FILE_READ_AT(fhv,&
          offset1,&
          h_x,&
          size(h_x),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)

  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_READ_AT in cgyro_read_restart_one failed')
     return
  endif

  call MPI_FILE_CLOSE(fhv,i_err)
  call MPI_INFO_FREE(finfo,i_err)

  call system_clock(cp_time,count_rate,count_max)
  if (cp_time.gt.start_time) then
    cp_dt = (cp_time-start_time)/real(count_rate)
  else
    cp_dt = (cp_time-start_time+count_max)/real(count_rate)
  endif

  if (i_proc == 0) then
    call date_and_time(sdate,stime,szone);
    open(NEWUNIT=statusfd,FILE=trim(path)//runfile_startups,action="write",access="append",status="unknown")
    write(statusfd, '(a,a,a,a,a,a,a,a,a,a,a,a,a,1pe10.3)') sdate(1:4),"/",sdate(5:6),"/",sdate(7:8)," ", &
                    stime(1:2),":",stime(3:4),":",stime(5:10), szone, ' [READ CHECKPOINT] Restart checkpoint read time: ', cp_dt
    close(statusfd)
  endif

end subroutine cgyro_read_restart_one

