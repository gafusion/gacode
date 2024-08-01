!-----------------------------------------------------------------
! cgyro_init_kernel.f90
!  
! PURPOSE:  
!  Initializations/allocation for cgyro_kernel
!
! NOTES:
!  This can be called directly using the driver routine cgyro 
!  (in which case input data will read from input.dat) or called 
!  as a subroutine using cgyro_sub.
!-----------------------------------------------------------------

subroutine cgyro_init_kernel

  use timer_lib
  use mpi
  use cgyro_globals
  use cgyro_step
  use cgyro_io

  implicit none

  character(8)  :: sdate
  character(10) :: stime
  character(len=64)  :: platform
  integer(KIND=8) :: aftermpi_time, beforetotal_time
  real :: mpi_dt,init_dt
  integer :: statusfd

  ! Need to initialize the info runfile very early
  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//runfile_info,status='replace')
     close(io)
  endif
  
  ! initialize tiny float
  small = tiny(0.0)
  
  ! the time_lib relies on MPI being initalized, so need to use lower level functions for this
  call system_clock(kernel_start_time,kernel_count_rate,kernel_count_max)
  
  i_time = 0

  ! Number of fluxes to output 
  nflux = 3+exch_flag

  ! 1. MPI setup
  call cgyro_mpi_grid
  if (error_status > 0) then
     call cgyro_final_kernel
     return
  endif
  
  call system_clock(aftermpi_time,kernel_count_rate,kernel_count_max)
  if (aftermpi_time > kernel_start_time) then
     mpi_dt = (aftermpi_time-kernel_start_time)/real(kernel_count_rate)
  else
     mpi_dt = (aftermpi_time-kernel_start_time+kernel_count_max)/real(kernel_count_rate)
  endif

  ! 2. Profile setup
  call cgyro_make_profiles

  ! 3. Parameter consistency checks
  call cgyro_check
  if (error_status > 0) then
     call cgyro_final_kernel
     return
  endif

  ! 4. Array initialization and construction
  !    NOTE: On exit, field_old = field 
  call cgyro_init_manager
  if (error_status > 0) then
     ! something went terribly wrong, hard abort, as things may be
     ! in weird state
     call cgyro_final_kernel
     return
  endif

  if (test_flag == 1) return
  
  if (restart_flag == 1) then
     io_control = 3*(1-silent_flag)
  else
     io_control = 1*(1-silent_flag)
  endif

  call timer_lib_in('io_init')
  call cgyro_write_timedata
  call timer_lib_out('io_init')
  call write_timers(trim(path)//runfile_timers)

  io_control = 2*(1-silent_flag)

  call system_clock(beforetotal_time,kernel_count_rate,kernel_count_max)
  if (beforetotal_time > kernel_start_time) then
     init_dt = (beforetotal_time-kernel_start_time)/real(kernel_count_rate)
  else
     init_dt = (beforetotal_time-kernel_start_time+kernel_count_max)/real(kernel_count_rate)
  endif

  if (i_proc == 0) then
     call date_and_time(sdate,stime);
     call get_environment_variable('GACODE_PLATFORM',platform)
     open(NEWUNIT=statusfd,FILE=trim(path)//runfile_startups,action='write',status='unknown',position='append')
     write(statusfd,'(14(a),f7.3,a,f7.3,a)') &
          sdate(1:4),'/',sdate(5:6),'/',sdate(7:8),' ', &
          stime(1:2),':',stime(3:4),':',stime(5:6),' ', &
          trim(platform),' [  INITIALIZATION] Time =',init_dt,' (mpi init: ', mpi_dt, ')'
     close(statusfd)
  endif

  ! GPU versions of step_gk and coll work on the following in the GPU memory
  call timer_lib_in('str_mem')
#if defined(OMPGPU)
!$omp target update to(field,cap_h_c,h_x,source)
#elif defined(_OPENACC)
!$acc update device(field,cap_h_c,h_x,source)
#endif
  call timer_lib_out('str_mem')

  ! Initialize adaptive time-stepping parameter
  delta_t_gk = delta_t
  
end subroutine cgyro_init_kernel
