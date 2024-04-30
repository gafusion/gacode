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
  character(len=64) :: platform
  integer(KIND=8) :: start_time,cp_time
  integer(KIND=8) :: count_rate, count_max
  real :: cp_dt
  integer :: j,ic0,statusfd

  ! use system_clock to be consistent with cgyro_kernel
  call system_clock(start_time,count_rate,count_max)

  !-----------------------------------------------
  ! Write h_x [filling (0,0) with source]
  !
  filemode = IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE)
  disp     = 0

  ! Pack source into h(0,0)
  if (source_flag == 1 .and. nt1 == 0) then
     ic0 = (n_radial/2)*n_theta
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd
#elif defined(_OPENACC)
!$acc parallel loop present(h_x,source)
#endif
     do j=1,n_theta
        h_x(ic0+j,:,0) = source(j,:,0)
     enddo
  endif

#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target update from(h_x)
#elif defined(_OPENACC)
!$acc update host(h_x) async(2)
#endif

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*(i_proc_1+i_proc_2*n_proc_1) + restart_header_size
  if (offset1 < restart_header_size) then
     call cgyro_error('Overflow in cgyro_write_restart')
     return
  endif

  if ((i_proc == 0) .and. (restart_preservation_mode<4)) then 
    ! Anything but restart_preservation_mode == 4
    ! User does not want high guarntees for the old file
    ! So, remove .old file, if it exists
    i_err = UNLINK(trim(path)//runfile_restart//".old")
    ! NOTE: We will not check if it succeeded... not important, may not even exist (yet)
  endif

  if ((i_proc == 0) .and. (restart_preservation_mode<2)) then 
    ! restart_preservation_mode == 1
    ! User wants to save disk space
    ! So, remove existing restart file, if it exists
    i_err = UNLINK(trim(path)//runfile_restart)
    ! NOTE: We will not check if it succeeded... not important, may not even exist (yet)
  endif

  ! TODO Error handling
  call MPI_INFO_CREATE(finfo,i_err)

  if (mpiio_stripe_factor > 0) then
    ! user asked us to explicitly set the MPI IO striping factor
    call MPI_INFO_SET(finfo,"striping_factor",mpiio_stripe_str,i_err)
  endif

  ! write to a temp file name first, so we don't end up with partially written files
  call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart//".part",&
          filemode,&
          finfo,&
          fhv,&
          i_err)
  if (i_err /= 0) then
     call cgyro_error('MPI_FILE_OPEN in cgyro_write_restart failed')
     return
  endif

  call MPI_FILE_SET_VIEW(fhv,&
          disp,&
          MPI_COMPLEX16,&
          MPI_COMPLEX16,&
          'native',&
          finfo,&
          i_err)

  ! need h_x here
#if defined(OMPGPU)
  ! no async for OMPGPU for now
#elif defined(_OPENACC)
!$acc wait(2)
#endif

  call MPI_FILE_WRITE_AT(fhv,&
          offset1,&
          h_x,&
          size(h_x),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)

  if (i_err /= 0) then
     call cgyro_error('MPI_FILE_WRITE_AT in cgyro_write_restart failed')
     return
  endif

  call MPI_FILE_SYNC(fhv,i_err)
  if (i_err /= 0) then
     call cgyro_error('MPI_FILE_SYNC in cgyro_write_restart failed')
     return
  endif

  call MPI_FILE_CLOSE(fhv,i_err)
  if (i_err /= 0) then
     call cgyro_error('MPI_FILE_CLOSE in cgyro_write_restart failed')
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
    if (restart_preservation_mode>2) then 
       ! restart_preservation_mode == 3 or 4
       ! First try to save any existing restart file as old
       i_err = RENAME(trim(path)//runfile_restart, trim(path)//runfile_restart//".old")
       ! NOTE: We will not check if it succeeded... not important, may not even exist (yet)
    endif

    ! Rename part into the final expected file name
    i_err = RENAME(trim(path)//runfile_restart//".part", trim(path)//runfile_restart)
    if (i_err /= 0) then
       call cgyro_error('Final rename in cgyro_write_restart failed')
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
    call date_and_time(sdate,stime);
    call get_environment_variable('GACODE_PLATFORM',platform)
    open(NEWUNIT=statusfd,FILE=trim(path)//runfile_startups,action='write',status='unknown',position='append')
    write(statusfd,'(14(a),f7.3)') &
         sdate(1:4),'/',sdate(5:6),'/',sdate(7:8),' ', &
         stime(1:2),':',stime(3:4),':',stime(5:6),' ', &
         trim(platform),' [CHECKPOINT WRITE] Time =',cp_dt
    close(statusfd)
  endif

  ! Re-set h(0,0)=0
  if (source_flag == 1 .and. nt1 == 0) then
     ic0 = (n_radial/2)*n_theta
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd
#elif defined(_OPENACC)
!$acc parallel loop present(h_x)
#endif
     do j=1,n_theta
        h_x(ic0+j,:,0) = 0.0
     enddo
  endif
  
end subroutine cgyro_write_restart_one

subroutine cgyro_write_restart_header_part
  use cgyro_globals
  use cgyro_io

  !---------------------------------------------------
  implicit none

  integer, parameter :: version = 2
  integer :: recid

  ! Different compilers have a different semantics for RECL... must use inquire
  integer :: reclen
  integer, dimension(1) :: recltest

  inquire(iolength=reclen) recltest

  open(unit=io,&
       file=trim(path)//runfile_restart//".part",&
       status='old',access='DIRECT',RECL=reclen)

  recid = 1
  write(io,REC=recid) restart_magic
  recid = recid + 1
  write(io,REC=recid) version
  recid = recid + 1
  write(io,REC=recid) n_theta
  recid = recid + 1
  write(io,REC=recid) n_radial
  recid = recid + 1
  write(io,REC=recid) n_species
  recid = recid + 1
  write(io,REC=recid) n_xi
  recid = recid + 1
  write(io,REC=recid) n_energy
  recid = recid + 1
  write(io,REC=recid) n_toroidal
  recid = recid + 1
  write(io,REC=recid) mpi_rank_order
  recid = recid + 1
  write(io,REC=recid) n_proc
  recid = recid + 1
  write(io,REC=recid) restart_magic ! just to have a clean end

  close(io)
end subroutine cgyro_write_restart_header_part
