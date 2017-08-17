!------------------------------------------------------------
! create.f90
!
! PURPOSE:
!------------------------------------------------------------

program create

  use create_globals
  use EXPRO_interface

  implicit none

  integer :: ierr, i

  call create_read_input 
  call create_check

  EXPRO_ctrl_quasineutral_flag = 0
  EXPRO_ctrl_z(1:3) = exm_z(1:3)
  EXPRO_ctrl_numeq_flag = 0 
  EXPRO_ctrl_n_ion = 3

  ! We're going to see if this file exists
  open(unit=1,file='input.profiles.gen',status='old',iostat=ierr)

  if (ierr > 0) then

     ! Starting from scratch (no input.profiles)

     ! Setting EXPRO_n_exp > 0 skips reading input.profiles

     EXPRO_n_exp  = nx

     call EXPRO_alloc('./',1) 

     set_exm_b_ref  =1
     set_exm_arho   =1
     set_exm_rho    =1
     set_exm_rmin   =1
     set_exm_rmaj   =1
     set_exm_q      =1
     set_exm_kappa  =1
     set_exm_delta  =1
     set_exm_te     =1
     set_exm_ne     =1
     set_exm_z_eff  =1
     set_exm_w0     =1
     set_exm_zeta   =1
     set_exm_zmag   =1
     do i=1,nions_max
        set_exm_ni(i)=1
        set_exm_ti(i)=1
     enddo
     set_exm_flow_mom  =1
     set_exm_pow_e     =1
     set_exm_pow_i     =1
     set_exm_pow_ei    =1
     set_exm_flow_beam =1
     set_exm_flow_wall =1
     set_exm_ptot      =1
     set_exm_polflux   =1
     set_exm_vpol(:)   =1
     set_exm_vtor(:)   =1

  else

     ! Starting from pre-existing input.profiles

     close(1)
     call EXPRO_alloc('./',1) 
     call EXPRO_read

  endif

  call create_set
  call create_write
  call EXPRO_alloc('./',0) 

end program create
