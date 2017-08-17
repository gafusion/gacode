!--------------------------------------------------------------
! prgen_read_inputprofiles.f90
!
! PURPOSE:
!  Read input.profiles.gen
!--------------------------------------------------------------

subroutine prgen_read_inputprofiles

  use prgen_globals
  use EXPRO_interface

  implicit none

  EXPRO_ctrl_quasineutral_flag = 0
  EXPRO_ctrl_z(1:3) = 1.0 
  EXPRO_ctrl_numeq_flag = 0 
  EXPRO_ctrl_n_ion = 10
  
  call EXPRO_alloc('./',1) 
  call EXPRO_read

  nx    = EXPRO_n_exp
  n_ion = EXPRO_n_ion

  call allocate_internals

  rmin(:) = EXPRO_rmin(:)
  rmaj(:) = EXPRO_rmaj(:)

end subroutine prgen_read_inputprofiles
