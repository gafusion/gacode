!--------------------------------------------------------------
! prgen_read_null.f90
!
! PURPOSE:
!  Minimal setup for gfile reading.
!--------------------------------------------------------------

subroutine prgen_read_null

  use prgen_globals

  implicit none

  nx = 40

  call allocate_internals
 
end subroutine prgen_read_null
