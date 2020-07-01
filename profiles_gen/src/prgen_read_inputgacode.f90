!--------------------------------------------------------------
! prgen_read_inputgacode.f90
!
! PURPOSE:
!  Read input.gacode
!--------------------------------------------------------------

subroutine prgen_read_inputgacode

  use prgen_globals
  use expro

  implicit none

  expro_ctrl_quasineutral_flag = 0
  expro_ctrl_numeq_flag = 0

  call expro_read(file_state)

  nx = expro_n_exp

  call prgen_allocate

  ! Needed for diagnostic printing
  rmin(:) = expro_rmin(:)
  rmaj(:) = expro_rmaj(:)

end subroutine prgen_read_inputgacode
