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

  integer :: restart_fmt

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

  if (i_proc == 0) then
     open(unit=io,&
          file=trim(path)//runfile_restart,&
          status='old',iostat=i_err)
     close(io)
     if (i_err == 0) then
        restart_fmt = 1
     else
        ! just assume it is the old format
        restart_fmt = 2
     endif
  endif

  call MPI_BCAST(restart_fmt         ,&
       1,MPI_INTEGER,0,CGYRO_COMM_WORLD,i_err)

  if (restart_fmt == 1) then
     call cgyro_read_restart_one
  else
     call cgyro_read_restart_one_old
  endif

end subroutine cgyro_read_restart

subroutine cgyro_read_restart_one_old

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
     call cgyro_error('ERROR: (CGYRO) overflow detected in cgyro_read_restart_one_old')
     return
  endif

  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart_old,&
          filemode,&
          finfo,&
          fhv,&
          i_err)

  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_OPEN in cgyro_read_restart_one_old failed')
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
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_READ_AT in cgyro_read_restart_one_old failed')
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
    open(NEWUNIT=statusfd,FILE=trim(path)//runfile_startups,action="write",status="unknown",position='append')
    write(statusfd, '(a,a,a,a,a,a,a,a,a,a,a,a,a,1pe10.3)') sdate(1:4),"/",sdate(5:6),"/",sdate(7:8)," ", &
                    stime(1:2),":",stime(3:4),":",stime(5:10), szone, ' [READ CHECKPOINT] Restart checkpoint read time: ', cp_dt
    close(statusfd)
  endif

  ! re-create the restart file using the new format
  if (i_proc == 0) then
     RENAME(trim(path)//runfile_restart_old, trim(path)//runfile_restart_old//".old")
  endif

  call MPI_BARRIER(CGYRO_COMM_WORLD,i_err)
  call cgyro_write_restart

end subroutine cgyro_read_restart_one_old

subroutine cgyro_read_restart_verify
  use cgyro_globals
  use cgyro_io

  !---------------------------------------------------
  implicit none

  integer :: magic, version
  integer :: t_n_theta,t_n_radial,t_n_species,t_n_xi,t_n_energy,t_n_toroidal

  open(unit=io,&
       file=trim(path)//runfile_restart,&
       status='old')

  read(io,"(i10)") magic
  if ( magic /= restart_magic) then
     close(io)
     call cgyro_error('ERROR: (CGYRO) Wrong magic number in restart header')
     return
  endif
   
  read(io,"(i4)") version
  if ( version /= 2) then
     close(io)
     call cgyro_error('ERROR: (CGYRO) Wrong version in restart header, only v2 supported')
     return
  endif

  read(io,"(6(i8,1x)") t_n_theta,t_n_radial,t_n_species,t_n_xi,t_n_energy,t_n_toroidal
  if ( (t_n_theta/=n_theta) .or. (t_n_radial/=n_radial) .or. &
       (t_n_species/=n_species) .or. (t_n_xi/=n_xi) .or. &
       (t_n_energy/=n_energy) .or. (t_n_toroidal/=n_toroidal) ) then
     close(io)
     call cgyro_error('ERROR: (CGYRO) Wrong geometry in restart header')
     return
  endif

  ! follow MPI params... will ignore them for v2, as they do not change the file format

  close(io)
end subroutine cgyro_read_restart_verify


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

  ! First read the header, and verify that it is compatible with current setup
  
  if (i_proc == 0) then
     call cgyro_read_restart_verify
     if (error_status /=0 ) return
  endif

  filemode = MPI_MODE_RDONLY
  disp     = 0

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*(i_proc_1+i_proc_2*n_proc_1) + restart_header_size
  if (offset1 < restart_header_size) then
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
    open(NEWUNIT=statusfd,FILE=trim(path)//runfile_startups,action="write",status="unknown",position='append')
    write(statusfd, '(a,a,a,a,a,a,a,a,a,a,a,a,a,1pe10.3)') sdate(1:4),"/",sdate(5:6),"/",sdate(7:8)," ", &
                    stime(1:2),":",stime(3:4),":",stime(5:10), szone, ' [READ CHECKPOINT] Restart checkpoint read time: ', cp_dt
    close(statusfd)
  endif

end subroutine cgyro_read_restart_one


