!-----------------------------------------------------------
! qlgyro_catch_error.f90
!
! PURPOSE:
!  Routine to print error message, finalize MPI, and
!  stop program execution gracefully. 
!-----------------------------------------------------------

subroutine qlgyro_catch_error(message)

  use qlgyro_globals

  implicit none

  character (len=*), intent(in) :: message

  if (i_proc_global == 0) then
     open(unit=1,file=trim(runfile),status='old',position='append')
     write(1,10) message
     close(1)
     print '(a)','** TGYRO HALTED ON ERROR ** (see out.tgyro.run)'
  endif

  call MPI_finalize(ierr)
  stop

10 format(a)

end subroutine qlgyro_catch_error
