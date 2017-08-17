!---------------------------------
! GKEIGEN_do.stub.f90
!
! PURPOSE:
!  Stub file for systems without
!  PETSc and SLEPc availability.
!
!---------------------------------

subroutine GKEIGEN_do

  use gyro_globals


  If (i_proc == 0) Then

    open(unit=10,file=trim(precfile),status='replace')
    write(10,*) 'GKEIGEN not supported.'
    close(10)

    print *, "Run terminated: Eigensolver unavailable on this platform"

    open(unit=1,file=trim(runfile),status='old',position='append')

    write (1,*) "Eigensolver unavailable on this platform"

    close (unit=1)

  EndIf

end subroutine GKEIGEN_do
