!-----------------------------------------------------------
! vgen_harvest.f90
!
! PURPOSE:
!  subroutine to send vgen_compute_neo input and output parameters
!  to harvest
!---------------------------------------------------------

module vgen_harvest

  implicit none
  
  integer :: ierr
  character(len=65507) :: harvest_sendline
  character(len=255) :: harvest_tag
  character NUL
  parameter(NUL = char(0))
  character(len=2) :: NUM

  include 'harvest_lib.inc'

contains
    
  subroutine vgen_harvest_inputandoutput

    use mpi
    use vgen_globals
    use neo_interface
    use EXPRO_interface

    implicit none 

    ierr=init_harvest('Vgen_jbs'//NUL,harvest_sendline,len(harvest_sendline))
    ierr=set_harvest_protocol('UDP'//NUL)
    ierr=set_harvest_port(41000)
    ierr=set_harvest_payload_int(harvest_sendline,'vel_method'//NUL,vel_method)
    ierr=set_harvest_payload_int(harvest_sendline,'erspecies_indx'//NUL,erspecies_indx)
    ierr=set_harvest_payload_int(harvest_sendline,'nth_min'//NUL,nth_min)
    ierr=set_harvest_payload_int(harvest_sendline,'nth_max'//NUL,nth_max)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_RHO'//NUL,EXPRO_rho,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'PFLUX_SUM'//NUL,pflux_sum,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'JBS_NEO'//NUL,jbs_neo,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'JBS_SAUTER'//NUL,jbs_sauter,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'JBS_NCLASS'//NUL,jbs_nclass,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'JBS_KOH'//NUL,jbs_koh,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'KY_SPECTRUM'//NUL,jbs_neo,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_rmin'//NUL,EXPRO_rmin,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_rmaj'//NUL,EXPRO_rmaj,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_q'//NUL,EXPRO_q,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_kappa'//NUL,EXPRO_kappa,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_delta'//NUL,EXPRO_delta,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_z_eff'//NUL,EXPRO_z_eff,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_w0'//NUL,EXPRO_w0,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_flow_mom'//NUL,EXPRO_rmin,EXPRO_flow_mom)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_sbcx'//NUL,EXPRO_sbcx,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_sbeame'//NUL,EXPRO_sbeame,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_sscxl'//NUL,EXPRO_sscxl,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_pow_e'//NUL,EXPRO_pow_e,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_pow_i'//NUL,EXPRO_pow_i,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_pow_ei'//NUL,EXPRO_pow_ei,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_zeta'//NUL,EXPRO_zeta,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_flow_beam'//NUL,EXPRO_flow_beam,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_flow_wall'//NUL,EXPRO_flow_wall,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_zmag'//NUL,EXPRO_zmag,EXPRO_n_exp)                                    
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_ptot'//NUL,EXPRO_ptot,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_polflux'//NUL,EXPRO_polflux,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_pow_e_fus'//NUL,EXPRO_pow_e_fus,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_pow_i_fus'//NUL,EXPRO_pow_i_fus,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_pow_e_sync'//NUL,EXPRO_pow_e_sync,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_pow_e_brem'//NUL,EXPRO_pow_e_brem,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_pow_e_line'//NUL,EXPRO_pow_e_line,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_pow_e_aux'//NUL,EXPRO_pow_e_aux,EXPRO_n_exp)
    ierr=set_harvest_payload_dbl_array(harvest_sendline,'EXPRO_pow_i_aux'//NUL,EXPRO_pow_i_aux,EXPRO_n_exp)

    ierr=harvest_send(harvest_sendline)

  end subroutine vgen_harvest_inputandoutput

end module vgen_harvest
