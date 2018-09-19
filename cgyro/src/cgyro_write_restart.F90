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

  implicit none
  
  ! Print this data on restart steps only; otherwise exit now
  if (mod(i_time,restart_step*print_step) /= 0) return

  call cgyro_write_restart_one

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
#ifdef __INTEL_COMPILER
  ! ifort defined rename in the ifport module
  use ifport
#endif

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

#ifndef __INTEL_COMPILER
  integer :: rename
#endif

  character(8)  :: sdate
  character(10) :: stime
  character(5)  :: szone
  integer(KIND=8) :: start_time,cp_time
  integer(KIND=8) :: count_rate, count_max
  real :: cp_dt
  integer :: statusfd

  ! use system_clock to be consistent with cgyro_kernel
  call system_clock(start_time,count_rate,count_max)

  !-----------------------------------------------
  ! Dump h and blending coefficients:
  !
  filemode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
  disp     = 0

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*(i_proc_1+i_proc_2*n_proc_1) + restart_header_size
  if (offset1 < restart_header_size) then
     call cgyro_error('ERROR: (CGYRO) overflow in cgyro_write_restart')
     return
  endif

  ! TODO Error handling
  call MPI_INFO_CREATE(finfo,i_err)

  ! write to a temp file name first, so we don't end up with partially written files
  call MPI_INFO_SET(finfo,"striping_factor",mpiio_stripe_str,i_err)

  call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart//".part",&
          filemode,&
          finfo,&
          fhv,&
          i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_OPEN in cgyro_write_restart failed')
     return
  endif

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
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_WRITE_AT in cgyro_write_restart failed')
     return
  endif

  call MPI_FILE_SYNC(fhv,i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_SYNC in cgyro_write_restart failed')
     return
  endif

  call MPI_FILE_CLOSE(fhv,i_err)
  if (i_err /= 0) then
     call cgyro_error('ERROR: (CGYRO) MPI_FILE_CLOSE in cgyro_write_restart failed')
     return
  endif

  call MPI_INFO_FREE(finfo,i_err)

  ! now update the header
  call MPI_BARRIER(CGYRO_COMM_WORLD,i_err)
  if (i_proc == 0) then 
     call cgyro_write_restart_header_part
     if (error_status /=0 ) return
  endif

  ! now that we know things worked well, move the file in its final location
  if (i_proc == 0) then 
    ! but first try to save any existing file
    i_err = RENAME(trim(path)//runfile_restart, trim(path)//runfile_restart//".old")
    ! NOTE: We will not check if it succeeded... not important, may not even exist (yet)

    i_err = RENAME(trim(path)//runfile_restart//".part", trim(path)//runfile_restart)
    if (i_err /= 0) then
       call cgyro_error('ERROR: (CGYRO) Final rename in cgyro_write_restart failed')
       return
    endif
  endif

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
                    stime(1:2),":",stime(3:4),":",stime(5:10), szone, ' [CHECKPOINTED] Checkpoint time: ', cp_dt
    close(statusfd)
  endif

end subroutine cgyro_write_restart_one

subroutine cgyro_write_restart_header_part
  use cgyro_globals
  use cgyro_io

  !---------------------------------------------------
  implicit none

  integer :: magic, version

  open(unit=io,&
       file=trim(path)//runfile_restart,&
       status='old')

  write(io,"(i10)") restart_magic
  write(io,"(i4)") 2
  write(io,"(6(i8,1x)") n_theta,n_radial,n_species,n_xi,n_energy,n_toroidal
  write(io,"(i2)") mpi_rank_order
  write(io,"(i10)") n_proc
  write(io,"(i2)") 0 ! just to have a clean end

  ! follow MPI params... will ignore them for v2, as they do not change the file format

  close(io)
end subroutine cgyro_write_restart_header_part
