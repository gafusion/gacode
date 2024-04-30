!--------------------------------------------------------------
! gyro_read_input.f90
!
! PURPOSE:
!  Complete read of input parameters.
!
! NOTES:
!  Input parameters are automatically broadcast 
!  to entire GYRO_COMM_WORLD. 
!
!   =============================================
!   HOW TO ADD INPUT VARIABLES 
!
!   Adding a new input variable is now simple.
!   Just append a line analogous to one of those 
!   below.  The variable will be read from the 
!   file 'input.dat' and broadcast to all members 
!   of GYRO_COMM_WORLD.
!   =============================================
!---------------------------------------------------------------

subroutine qlgyro_read_input

  use qlgyro_globals

  !-------------------------
  implicit none
  !
  integer :: i_kypx0
  !-------------------------

  if (i_proc_global == 0) open(unit=1,file=trim(path)//'input.qlgyro.gen',status='old')

  !--------------------------------------------------------
  ! BEGIN reading data:
  !
  call ql_readbc_int(n_parallel)
  call ql_readbc_int(n_runs)
  call ql_readbc_real(gamma_exb)
  call ql_readbc_int(code)
  call ql_readbc_int(kygrid_model)
  call ql_readbc_real(ky_in)
  call ql_readbc_int(n_ky)
  call ql_readbc_int(xnu_model)
  call ql_readbc_int(auto_box_size)
  call ql_readbc_real(kx_max_box)
  call ql_readbc_int(sat_rule)
  call ql_readbc_int(n_px0)
  call ql_readbc_int(px0grid_model)
  call ql_readbc_int(restart_mode)

  ! DONE reading data.
  !--------------------------------------------------------
  if (mod(n_inst, n_parallel) .ne. 0) then
     open(unit=1,file=trim(runfile),position='append')
     write(1,*) '------------------------------------------------------------'
     write(1,*) 'No. of processors cant be divided in number of parallel jobs'
     write(1,*) '------------------------------------------------------------'
     stop
  end if

  if (sat_rule .eq. 1 .and. n_px0 .ne. 1) then
     open(unit=1,file=trim(runfile),position='append')
     write(1,*) '---------------------------------------------------------'
     write(1,*) 'Saturation rule 1 can only be run with a single PX0 value'
     write(1,*) '---------------------------------------------------------'
     stop
  end if

  if (gamma_exb .ne. 0.0) then
     qlgyro_rotation_flag = 1
  endif

  procs = n_inst/n_parallel
  
  if (i_proc_global == 0) close(1)

  if (debug_flag == 1 .and. i_proc_global == 0) then
     print *,'[qlgyro_read_input done]'
  endif

end subroutine qlgyro_read_input

!------------------------------------------------------------
! Service routines: 
!
! (1) read and broadcast an integer:
!
subroutine ql_readbc_int(p)

  use mpi
  use qlgyro_globals

  implicit none
  integer, intent(inout) :: p

  if (i_proc_global == 0) read(1,*) p

  call MPI_BCAST(p,1,MPI_INTEGER,0,QLGYRO_COMM_WORLD,ierr)

end subroutine ql_readbc_int
!
! (2) read and broadcast a real:
!
subroutine ql_readbc_real(x)
  
  use mpi
  use qlgyro_globals

  implicit none
  real, intent(inout) :: x

  if (i_proc_global == 0) read(1,*) x

  call MPI_BCAST(x,1,MPI_DOUBLE,0,QLGYRO_COMM_WORLD,ierr)

end subroutine ql_readbc_real
!------------------------------------------------------------
