!-----------------------------------------------------------------
! qlgyro_trap_component_error.f90
!
! PURPOSE:
!  Manage return status from components (NEO, TGLF, GYRO).
! 
!  success:
!   status=0
!   message='null'
!
!  failure
!   status=1
!   message=<error description>
!
!  warning
!   status=2
!   message=<warning description>
!-------------------------------------------------------------------

subroutine qlgyro_trap_component_error(status,message)

  use mpi
  use qlgyro_globals

  implicit none

  integer :: i
  ! INPUT
  integer, intent(in) :: status
  character (len=80), intent(in) :: message
  ! INTERNAL
  integer, dimension(n_inst) :: status_vec
  character (len=80), dimension(n_inst) :: message_vec


  ! Collect distributed error messages into global vectors
  call MPI_ALLGATHER(status,1,MPI_INTEGER,status_vec,1,MPI_INTEGER,qlgyro_adj,ierr)
  call MPI_ALLGATHER(message,80,MPI_CHARACTER,message_vec,80,MPI_CHARACTER,qlgyro_adj,ierr)

  ! Broadcast vectors to all members of a given instance 
  ! (Required for GYRO runs, not NEO or TGLF):
  call MPI_BCAST(status_vec,n_inst,MPI_INTEGER,0,qlgyro_comm,ierr)
  call MPI_BCAST(message_vec,80*n_inst,MPI_CHARACTER,0,qlgyro_comm,ierr)
  
 ! Serial logic from here onward

  
  if (sum(status_vec) > 0) then
     if (i_proc_global == 0) then
        open(unit=1,file=trim(runfile),position='append')
        do i=1,n_inst
           write(1,*) trim(path),": ",trim(message_vec(i))
        enddo
        close(1)
     endif
  endif

end subroutine qlgyro_trap_component_error
