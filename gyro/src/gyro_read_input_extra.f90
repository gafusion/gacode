!---------------------------------------------------------------
! gyro_read_input_extra.f90
!
! PURPOSE:
!  Complete read of extra (restart and test) parameters.
!---------------------------------------------------------------

subroutine gyro_read_input_extra

  use mpi
  use gyro_globals
  use expro

  !-------------------------
  implicit none
  !
  logical :: lfe
  !-------------------------

  !-------------------------------------------------------------------
  ! Parse restart.dat.
  ! 
  ! -1=no save/restart
  !  0=begin new simulation
  !  1=restart, continue writing restart files
  !  2=restart, but do not write restart files
  !  3=restart, set t=0
  !
  if (i_proc == 0) then
     inquire(file=trim(path)//'restart.dat',exist=lfe)
     if (lfe .eqv. .true.) then
        open(unit=1,file=trim(path)//'restart.dat',status='old')
        read(1,*) restart_method
     else
        restart_method = 0
     endif
  endif
  call MPI_BCAST(restart_method,1,MPI_INTEGER,0,GYRO_COMM_WORLD,i_err)
  if (i_proc == 0 .and. lfe .eqv. .true.) close(1)
  !------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Check for gyrotest flag by reading a single integer.
  !
  ! 0: not a test
  ! 1: run in single-processor test mode
  !
  if (i_proc == 0) then
     inquire(file=trim(path)//'gyrotest_flag',exist=lfe)
     if (lfe .eqv. .true.) then
        open(unit=1,file=trim(path)//'gyrotest_flag',status='old')
        read(1,*) gyrotest_flag
     else
        gyrotest_flag = 0
     endif
  endif
  call MPI_BCAST(gyrotest_flag,1,MPI_INTEGER,0,GYRO_COMM_WORLD,i_err)
  if (i_proc == 0 .and. lfe .eqv. .true.) close(1)
  !------------------------------------------------------------------

  !-------------------------------------------------------- 
  ! Read geometry Fourier coefficients if they exist.
  ! 
  if (i_proc == 0) then
     inquire(file=trim(path)//'input.geo',exist=lfe)
     if (lfe .eqv. .true.) then
        open(unit=1,file=trim(path)//'input.geo',status='old')
        call expro_skip_header(1)
        read(1,*) n_fourier_geo
        read(1,*) a_fourier_geo(:,0:n_fourier_geo)
        close(1)
     endif
  endif
  call MPI_BCAST(lfe,1,MPI_LOGICAL,0,GYRO_COMM_WORLD,i_err)
  if (lfe .eqv. .true.) then
     call MPI_BCAST(n_fourier_geo,1,MPI_INTEGER,0,GYRO_COMM_WORLD,i_err)
     call MPI_BCAST(a_fourier_geo,size(a_fourier_geo),MPI_DOUBLE_PRECISION,0,GYRO_COMM_WORLD,i_err)
  else
     n_fourier_geo = 0
     a_fourier_geo(:,:) = 0.0
  endif
  !--------------------------------------------------------

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[read_input_extra done]'
  endif

 
end subroutine gyro_read_input_extra
