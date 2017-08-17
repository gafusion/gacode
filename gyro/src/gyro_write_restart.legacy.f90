!------------------------------------------------
! gyro_write_restart.f90 [caller: gyro_fulladvance]
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
  integer :: io
  !
  integer :: data_step_old
  integer :: n_proc_old
  integer :: i_restart_old
  !
  real :: t_current_old
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

  !-----------------------------------------------
  ! Dump h:
  !
  if (i_proc == 0) then
     open(unit=io,file=trim(path)//file_restart(i_restart),status='replace')
     write(io,fmtstr2) h
  endif

  do i_proc_w=1,n_proc-1

     if (i_proc == 0) then

        call MPI_RECV(h_0,&
             size(h_0),&
             MPI_DOUBLE_COMPLEX,&
             i_proc_w,&
             i_proc_w,&
             GYRO_COMM_WORLD,&
             recv_status,&
             i_err)

        write(io,fmtstr2) h_0

        ! Determine correct integer format for printing
        ! brackets:

        if (debug_flag == 1) then
           if (i_proc_w < 10) then
              write(6,'(a,i1,a)',advance='no') '[',i_proc_w,']'
           else if (i_proc_w < 100) then
              write(6,'(a,i2,a)',advance='no') '[',i_proc_w,']'
           else
              write(6,'(a,i3,a)',advance='no') '[',i_proc_w,']'
           endif
           if (modulo(i_proc_w,8) == 0) print *
        endif

     else if (i_proc == i_proc_w) then

        call MPI_SEND(h,&
             size(h),&
             MPI_DOUBLE_COMPLEX,&
             0,&
             i_proc_w,&
             GYRO_COMM_WORLD,&
             i_err)

     endif

  enddo

  if (i_proc == 0) close(io)
  !------------------------------------------------

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
  !---------------------------------------------------------

  ! ** Keep this consistent with gyro_read_restart.f90

end subroutine gyro_write_restart
