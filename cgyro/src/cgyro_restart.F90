!------------------------------------------------
! cgyro_restart.f90
!
! PURPOSE:
!  This is the master file controlling input and output of
!  restart data using MPI-IO.
!------------------------------------------------

module cgyro_restart

  implicit none

  integer, parameter :: restart_header_size = 1024
  integer, parameter :: restart_magic = 140974129
  integer, parameter :: restart_version = 3

  integer, private :: t_velocity_order
  integer, private :: t_nt_loc
  integer, private :: t_nv_loc

contains

subroutine cgyro_write_restart

  use mpi
  use cgyro_globals
  use cgyro_io

  implicit none
  
  integer :: j,ic0

  ! Print this data on restart steps only; otherwise exit now
  if (mod(i_time,restart_step*print_step) /= 0) return

  call cgyro_write_restart_one

  ! Unpack h(0,0) into source 
  if (source_flag == 1 .and. nt1 == 0) then
     ic0 = (n_radial/2)*n_theta
     do j=1,n_theta
        source(j,:,0) = h_x(ic0+j,:,0)
        h_x(ic0+j,:,0) = 0.0
     enddo
     sa = 0.0
     do j=1,nint(t_current/delta_t)
        sa = 1.0+exp(-delta_t/tau_ave)*sa
     enddo
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
  integer :: UNLINK
#endif

  character(8)  :: sdate
  character(10) :: stime
  character(len=64) :: platform
  integer(KIND=8) :: start_time,cp_time
  integer(KIND=8) :: count_rate, count_max
  real :: cp_dt
  integer :: j,ic0,statusfd
  integer :: ierr

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
    ierr  = UNLINK(trim(path)//runfile_restart//".old")
    ! NOTE: We will not check if it succeeded... not important, may not even exist (yet)
  endif

  if ((i_proc == 0) .and. (restart_preservation_mode<2)) then 
    ! restart_preservation_mode == 1
    ! User wants to save disk space
    ! So, remove existing restart file, if it exists
    ierr = UNLINK(trim(path)//runfile_restart)
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
     if (error_status > 0) return
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
  write(io,REC=recid) restart_version
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
  write(io,REC=recid) velocity_order
  recid = recid + 1
  write(io,REC=recid) nt_loc
  recid = recid + 1
  write(io,REC=recid) nv_loc
  recid = recid + 1
  write(io,REC=recid) nc_loc
  recid = recid + 1
  write(io,REC=recid) n_toroidal_procs
  recid = recid + 1
  write(io,REC=recid) n_proc
  recid = recid + 1
  write(io,REC=recid) mpi_rank_order
  recid = recid + 1
  write(io,REC=recid) delta_t_method
  recid = recid + 1
  write(io,REC=recid) nonlinear_flag
  recid = recid + 1
  write(io,REC=recid) n_jtheta
  recid = recid + 1
  write(io,REC=recid) nsplit
  recid = recid + 1
  write(io,REC=recid) nsplitA
  recid = recid + 1
  write(io,REC=recid) nsplitB
  recid = recid + 1
  write(io,REC=recid) nup_theta
  recid = recid + 1
  write(io,REC=recid) nv
  recid = recid + 1
  write(io,REC=recid) nc
  recid = recid + 1
  write(io,REC=recid) restart_magic ! just to have a clean end

  close(io)
end subroutine cgyro_write_restart_header_part
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

        open(unit=io,file=trim(path)//runfile_restart_tag,status='old')

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
     open(unit=io,file=trim(path)//runfile_restart,status='old',iostat=i_err)
     close(io)
  endif

  call cgyro_read_restart_one

end subroutine cgyro_read_restart


subroutine cgyro_read_restart_verify

  use cgyro_globals
  use cgyro_io

  !---------------------------------------------------
  implicit none

  integer :: magic, version, recid
  integer :: t_n_theta,t_n_radial,t_n_species,t_n_xi,t_n_energy,t_n_toroidal
  integer :: restart_magic_v2

  ! Different compilers have a different semantics for RECL... must use inquire
  integer :: reclen
  integer, dimension(1) :: recltest

  inquire(iolength=reclen) recltest

  open(unit=io,&
       file=trim(path)//runfile_restart,&
       status='old',access='DIRECT',RECL=reclen,ACTION='READ')

  recid = 1

  read(io,REC=recid) magic
  recid = recid + 1
  if ( magic == restart_magic) then
     read(io,REC=recid) version
     recid = recid + 1
     if ( version /= restart_version) then
        close(io)
        call cgyro_error('Wrong version in restart header, expected 3')
        return
     endif
  else
     ! older versions had different magic
     if (velocity_order == 1) then
       ! traditional ordering
       restart_magic_v2 = 140906808
     else
       ! alternative ordering, needed different magic
       restart_magic_v2 = 140916753
     endif
     if ( magic /= restart_magic_v2) then
        ! we don't support anything else, abort
        close(io)
        call cgyro_error('Wrong magic number in restart header')
        return
     endif
     read(io,REC=recid) version
     recid = recid + 1
     if ( version /= 2) then
        close(io)
        call cgyro_error('Wrong version in restart header, expected 2')
        return
     endif
  endif

  read(io,REC=recid) t_n_theta
  recid = recid + 1
  read(io,REC=recid) t_n_radial
  recid = recid + 1
  read(io,REC=recid) t_n_species
  recid = recid + 1
  read(io,REC=recid) t_n_xi
  recid = recid + 1
  read(io,REC=recid) t_n_energy
  recid = recid + 1
  read(io,REC=recid) t_n_toroidal
  recid = recid + 1
  if ( (t_n_theta/=n_theta) .or. (t_n_radial/=n_radial) .or. &
       (t_n_species/=n_species) .or. (t_n_xi/=n_xi) .or. &
       (t_n_energy/=n_energy) .or. (t_n_toroidal/=n_toroidal) ) then
     close(io)
     call cgyro_error('Wrong geometry in restart header')
     return
  endif

  if ( magic == restart_magic) then
     read(io,REC=recid) t_velocity_order
     recid = recid + 1
     read(io,REC=recid) t_nt_loc
     recid = recid + 1
     read(io,REC=recid) t_nv_loc
     recid = recid + 1
  else
     ! not recorded, so assume they are the same as current run
     t_velocity_order = velocity_order
     t_nt_loc = nt_loc
     t_nv_loc = nv_loc
  endif

  ! follow other params... will ignore them, as they do not change the file format

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
  character(len=64) :: platform
  integer(KIND=8) :: start_time,cp_time
  integer(KIND=8) :: count_rate, count_max
  real :: cp_dt
  integer, dimension(3) :: mpibuf
  integer :: ic0,j,statusfd

  ! use system_clock to be consistent with cgyro_kernel
  call system_clock(start_time,count_rate,count_max)

  ! First read the header, and verify that it is compatible with current setup
  
  if (i_proc == 0) then
     call cgyro_read_restart_verify
     if (error_status /=0 ) return
     mpibuf(1) = t_nt_loc
     mpibuf(2) = t_nv_loc
     mpibuf(3) = t_velocity_order
  endif

  call MPI_Bcast(mpibuf, 3, MPI_INTEGER, 0, CGYRO_COMM_WORLD, i_err)

  if (i_err /= 0) then
     call cgyro_error('MPI in restart failed')
     return
  endif

  if (i_proc /= 0) then
     t_nt_loc = mpibuf(1)
     t_nv_loc = mpibuf(2)
     t_velocity_order = mpibuf(3)
  endif

  if ( (t_nt_loc /= nt_loc) .or. (t_nv_loc /= nv_loc) .or. (t_velocity_order /= velocity_order) ) then
     ! the layout on disk does not align to current process distribution
     ! use the more complicated (and slower) logic
     call cgyro_read_restart_slow
     return
  endif

  filemode = MPI_MODE_RDONLY
  disp     = 0

  offset1 = size(h_x,kind=MPI_OFFSET_KIND)*(i_proc_1+i_proc_2*n_proc_1) + restart_header_size
  if (offset1 < restart_header_size) then
     call cgyro_error('Overflow detected in cgyro_read_restart_one')
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
     call cgyro_error('MPI_FILE_OPEN in cgyro_read_restart_one failed')
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
     call cgyro_error('MPI_FILE_READ_AT in cgyro_read_restart_one failed')
     return
  endif

  call MPI_FILE_CLOSE(fhv,i_err)
  call MPI_INFO_FREE(finfo,i_err)

  call system_clock(cp_time,count_rate,count_max)
  if (cp_time > start_time) then
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
         trim(platform),' [ CHECKPOINT READ] Time =',cp_dt
    close(statusfd)
  endif

end subroutine cgyro_read_restart_one

subroutine cgyro_read_restart_slow

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
  integer(kind=MPI_OFFSET_KIND) :: offset1,off_block,off_row,off_col
  !---------------------------------------------------

  integer, dimension(:,:,:), allocatable :: t_iv_v
  character(8)  :: sdate
  character(10) :: stime
  character(len=64) :: platform
  integer(KIND=8) :: start_time,cp_time
  integer(KIND=8) :: count_rate, count_max
  real :: cp_dt
  integer :: ic0,j,statusfd
  integer :: ie,ix,is,itor
  integer :: t_iv

  ! we already checked what the format is
  allocate(t_iv_v(n_energy,n_xi,n_species))

  ! Velocity pointers for the restart velocity_order
  t_iv = 0
  if (t_velocity_order==1) then
    !original
    do ie=1,n_energy
      do ix=1,n_xi
        do is=1,n_species
           t_iv = t_iv+1
           t_iv_v(ie,ix,is) = t_iv
        enddo
      enddo
    enddo
  else if (t_velocity_order==2) then
    ! optimized for minimizing species
    do is=1,n_species
      do ie=1,n_energy
        do ix=1,n_xi
           t_iv = t_iv+1
           t_iv_v(ie,ix,is) = t_iv
        enddo
      enddo
    enddo
  else
     call cgyro_error('Unknown restart VELOCITY_ORDER.')
     return
  endif


  ! use system_clock to be consistent with cgyro_kernel
  call system_clock(start_time,count_rate,count_max)

  filemode = MPI_MODE_RDONLY
  disp     = 0

  call MPI_INFO_CREATE(finfo,i_err)

  call MPI_FILE_OPEN(CGYRO_COMM_WORLD,&
          trim(path)//runfile_restart,&
          filemode,&
          finfo,&
          fhv,&
          i_err)

  if (i_err /= 0) then
     call cgyro_error('MPI_FILE_OPEN in cgyro_read_restart_one failed')
     return
  endif

  call MPI_FILE_SET_VIEW(fhv,&
          disp,&
          MPI_COMPLEX16,&
          MPI_COMPLEX16,&
          'native',&
          finfo,&
          i_err)

  do itor=nt1,nt2  ! itor is 0-based
     do iv=nv1,nv2
        t_iv = t_iv_v(ie_v(iv),ix_v(iv),is_v(iv)) 
        iv_loc = iv-nv1+1
        off_block = modulo(t_iv-1,t_nv_loc) + modulo(itor,t_nt_loc)*t_nv_loc
        off_row = ((t_iv-1)/t_nv_loc)*(t_nv_loc*t_nt_loc)
        off_col = (itor/t_nt_loc)*(nv*t_nt_loc)
        offset1 = size(h_x(:,iv_loc,itor),kind=MPI_OFFSET_KIND)*(off_block+off_row+off_col)
        offset1 = offset1 + restart_header_size
        if (offset1 < restart_header_size) then
           call MPI_FILE_CLOSE(fhv,i_err)
           call MPI_INFO_FREE(finfo,i_err)
           call cgyro_error('Overflow detected in cgyro_read_restart_one')
           return
        endif

        call MPI_FILE_READ_AT(fhv,&
          offset1,&
          h_x(:,iv_loc,itor),&
          size(h_x(:,iv_loc,itor)),&
          MPI_COMPLEX16,&
          fstatus,&
          i_err)

        if (i_err /= 0) then
           call MPI_FILE_CLOSE(fhv,i_err)
           call MPI_INFO_FREE(finfo,i_err)
           call cgyro_error('MPI_FILE_READ_AT in cgyro_read_restart_slow failed')
           return
        endif
     enddo
  enddo

  call MPI_FILE_CLOSE(fhv,i_err)
  call MPI_INFO_FREE(finfo,i_err)

  call system_clock(cp_time,count_rate,count_max)
  if (cp_time > start_time) then
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
         trim(platform),' [ CHECKPOINT READ (slow)] Time =',cp_dt
    close(statusfd)
  endif

  deallocate(t_iv_v)

end subroutine cgyro_read_restart_slow

end module cgyro_restart

