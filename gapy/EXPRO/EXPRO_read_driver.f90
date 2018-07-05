!--------------------------------------------------------------
! EXPRO_read_driver.f90
!
! PURPOSE:
!  Low-level read control for EXPRO.
!--------------------------------------------------------------

subroutine EXPRO_read_driver

  use EXPRO_globals
  use EXPRO_interface

  implicit none

  integer, parameter :: io=1
  integer :: ierr

  real, dimension(:), allocatable :: dummy

  !--------------------------------------------------------------
  ! READ generated (stripped) version of input.profiles:
  !
  open(unit=io,&
       file=trim(path)//'input.profiles.gen',&
       status='old',&
       iostat=ierr)

  read(io,*) EXPRO_shot
  read(io,*) EXPRO_n_ion
  read(io,*) EXPRO_n_exp
  read(io,*) EXPRO_b_ref
  read(io,*) EXPRO_ip_exp
  read(io,*) EXPRO_rvbv
  read(io,*) EXPRO_arho
  
  allocate(dummy(EXPRO_n_exp))

  ! 1-5
  read(io,*) EXPRO_rho(:)
  read(io,*) EXPRO_rmin(:)
  read(io,*) EXPRO_polflux(:)
  read(io,*) EXPRO_q(:)    
  read(io,*) EXPRO_w0(:) 

  ! 6-10
  read(io,*) EXPRO_rmaj(:)
  read(io,*) EXPRO_zmag(:)
  read(io,*) EXPRO_kappa(:)
  read(io,*) EXPRO_delta(:)
  read(io,*) EXPRO_zeta(:)

  ! 11-15
  read(io,*) EXPRO_ne(:)
  read(io,*) EXPRO_te(:)
  read(io,*) EXPRO_ptot(:)
  read(io,*) EXPRO_z_eff(:)
  read(io,*) dummy(:)

  ! 16-20
  read(io,*) EXPRO_ni(1,:)
  read(io,*) EXPRO_ni(2,:)
  read(io,*) EXPRO_ni(3,:)
  read(io,*) EXPRO_ni(4,:)
  read(io,*) EXPRO_ni(5,:)

  ! 21-25
  read(io,*) EXPRO_ni(6,:)
  read(io,*) EXPRO_ni(7,:)
  read(io,*) EXPRO_ni(8,:)
  read(io,*) EXPRO_ni(9,:)
  read(io,*) EXPRO_ni(10,:)

  ! 26-30
  read(io,*) EXPRO_ti(1,:)
  read(io,*) EXPRO_ti(2,:)
  read(io,*) EXPRO_ti(3,:)
  read(io,*) EXPRO_ti(4,:)
  read(io,*) EXPRO_ti(5,:)

  ! 31-35
  read(io,*) EXPRO_ti(6,:)
  read(io,*) EXPRO_ti(7,:)
  read(io,*) EXPRO_ti(8,:)
  read(io,*) EXPRO_ti(9,:)
  read(io,*) EXPRO_ti(10,:)

  ! 36-40
  read(io,*) EXPRO_vtor(1,:)
  read(io,*) EXPRO_vtor(2,:)
  read(io,*) EXPRO_vtor(3,:)
  read(io,*) EXPRO_vtor(4,:)
  read(io,*) EXPRO_vtor(5,:)

  ! 41-45
  read(io,*) EXPRO_vtor(6,:)
  read(io,*) EXPRO_vtor(7,:)
  read(io,*) EXPRO_vtor(8,:)
  read(io,*) EXPRO_vtor(9,:)
  read(io,*) EXPRO_vtor(10,:)

  ! 46-50
  read(io,*) EXPRO_vpol(1,:)
  read(io,*) EXPRO_vpol(2,:)
  read(io,*) EXPRO_vpol(3,:)
  read(io,*) EXPRO_vpol(4,:)
  read(io,*) EXPRO_vpol(5,:)

  ! 51-55
  read(io,*) EXPRO_vpol(6,:)
  read(io,*) EXPRO_vpol(7,:)
  read(io,*) EXPRO_vpol(8,:)
  read(io,*) EXPRO_vpol(9,:)
  read(io,*) EXPRO_vpol(10,:)

  ! 56-60
  read(io,*) EXPRO_flow_beam(:)
  read(io,*) EXPRO_flow_wall(:)
  read(io,*) EXPRO_flow_mom(:)
  read(io,*) dummy(:)
  read(io,*) dummy(:)

  ! 61-65
  read(io,*) EXPRO_pow_e(:)
  read(io,*) EXPRO_pow_i(:)
  read(io,*) EXPRO_pow_ei(:)
  read(io,*) EXPRO_pow_e_aux(:)
  read(io,*) EXPRO_pow_i_aux(:)

  ! 66-70
  read(io,*) EXPRO_pow_e_fus(:)
  read(io,*) EXPRO_pow_i_fus(:)
  read(io,*) EXPRO_pow_e_sync(:)
  read(io,*) EXPRO_pow_e_brem(:)
  read(io,*) EXPRO_pow_e_line(:)

  ! 71-75
  read(io,*,iostat=ierr) EXPRO_sbeame(:)
  if (ierr == 0) then
     read(io,*) EXPRO_sbcx(:)
     read(io,*) EXPRO_sscxl(:)
     read(io,*) dummy(:)
     read(io,*) dummy(:)
  else
     print('(a)'),'INFO: (EXPRO_read_driver) This is an OLD format input.profiles.  Please regenerate.'
  endif

  close(io)
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! READ general shape coefficients if they exist: 
  !
  if (EXPRO_nfourier > 0) then
     call EXPRO_read_geo
  endif
  !--------------------------------------------------------------

  deallocate(dummy)

end subroutine EXPRO_read_driver
 
 
