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

  integer :: i,ierr

  EXPRO_ctrl_quasineutral_flag = 0
  EXPRO_ctrl_z(1:3) = 1.0 
  EXPRO_ctrl_numeq_flag = 0 
  EXPRO_ctrl_n_ion = 10
  
  call EXPRO_alloc('./',1) 
  call EXPRO_read
  
  nx    = EXPRO_n_exp
  n_ion = EXPRO_n_ion

  call allocate_internals

  ! Needed for disagnostic printing
  rmin(:) = EXPRO_rmin(:)
  rmaj(:) = EXPRO_rmaj(:)

  ! Need to close then reopen as usual in map
  call EXPRO_alloc('./',0) 

  open(unit=1,file='profile_header',status='old')
  do i=1,n_ion
     read(1,*) ion_z(i),ion_mass(i),ion_type(i)
  enddo
  close(1)
     
end subroutine prgen_read_inputprofiles
