!-----------------------------------------------------------
! neo_run.f90
!
! PURPOSE:
!  Manage call to local NEO simulation.
!---------------------------------------------------------

subroutine neo_run()

  use neo_globals
  use neo_interface
  use neo_harvest

  implicit none

  integer :: is

  ! Map INTERFACE parameters -> GLOBAL variables
  call map_interface2global
  call neo_harvest_input
  ! Can exit if we are in test mode
  if (neo_test_flag_in == 1) return
  
  !----------------------------------------------------------------------
  ! Initialize error status
  error_status = 0
  error_message = '(NEO) completed successfully'

  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io_neoout,file=trim(path)//runfile_neoout,status='replace')
     close(io_neoout)
  endif
  !----------------------------------------------------------------------

  neo_dke_out=0.0
  neo_dke_1d_out= 0.0
  neo_gv_out=0.0
  neo_th_out=0.0
  neo_nclass_out=0.0
  neo_thHS_out=0.0

  if (neo_nn_flag_in .ne. 1) then
    ! Run NEO
    call neo_do  
  else
    call neo_jbs_nn
  endif

  ! ----------------------------------------------------------------------
  ! OUTPUT NORMALIZATION
  !
  ! NEO particle flux            : Gamma_norm = n0_norm vt_norm
  ! NEO energy flux              : Q_norm     = n0_norm vt_norm T0_norm
  ! NEO momentum flux (toroidal) : Pi_norm    = n0_norm T0_norm a_norm
  ! NEO boostrap current         : j_norm     = e n0_norm vt_norm B_unit
  ! where vt_norm = sqrt(T_norm / mass_norm)

  !!!!!!!!!!!!!!!!!!
  ! theory results
  !!!!!!!!!!!!!!!!!!

  ! Hinton-Hazeltine and Chang-Hinton fluxes
  neo_pflux_thHH_out   = neo_th_out(1)               ! Gamma_i/Gamma_norm 
  ! (HH - ambipolar)
  neo_eflux_thHHi_out  = neo_th_out(2)               ! Q_i/Q_norm (HH)
  neo_eflux_thHHe_out  = neo_th_out(3)               ! Q_e/Q_norm (HH)
  neo_eflux_thCHi_out  = neo_th_out(4)               ! Q_i/Q_norm (CH)

  ! Hirshman-Sigmar (multi-species) fluxes
  do is=1,11
     neo_pflux_thHS_out(is) = neo_thHS_out(is,1)     ! Gamma_s/Gamma_norm 
     neo_eflux_thHS_out(is) = neo_thHS_out(is,2)     ! Q_s/Q_norm
  enddo

  ! Hinton-Hazeltine and Sauter bootstrap current
  neo_jpar_thS_out    = neo_th_out(5)                ! <j B> / j_norm (Sauter)
  neo_jpar_thK_out    = neo_th_out(6)                ! <j B> / j_norm (Koh)
  neo_jpar_thN_out    = neo_th_out(7)                ! <j B> / j_norm (NCLASS)

  !!!!!!!!!!!!!!!!!!
  ! NCLASS results
  !!!!!!!!!!!!!!!!!!
  do is=1,11
     neo_pflux_nclass_out(is)    = neo_nclass_out(is,1)    ! Gamma_is/Gamma_norm
     neo_efluxtot_nclass_out(is) = neo_nclass_out(is,2)    ! Q_is/Q_norm total
     neo_vpol_nclass_out(is)     = neo_nclass_out(is,3)    ! v_theta/vt_norm (at theta=0)
     neo_vtor_nclass_out(is)     = neo_nclass_out(is,4)    ! v_phi/vt_norm (at theta=0)
     neo_nclassvis_out(is)       = neo_nclass_out(is,5)    ! ymu_s(1,1) / (m_norm*dens_norm*vth_norm/a_norm)
  enddo
  neo_jpar_nclass_out = neo_nclass_1d_out                  ! <j B> / j_norm
 
  !!!!!!!!!!!!!!!!!!
  ! dke results
  !!!!!!!!!!!!!!!!!!
  do is=1,11
     neo_pflux_dke_out(is)    = neo_dke_out(is,1)    ! Gamma_is/Gamma_norm
     neo_efluxtot_dke_out(is) = neo_dke_out(is,2)    ! Q_is/Q_norm total
     neo_mflux_dke_out(is)    = neo_dke_out(is,3)    ! Pi_is/Pi_norm
     neo_efluxncv_dke_out(is) = neo_dke_out(is,4)    ! Q_is/Q_norm non-convective
     neo_vpol_dke_out(is)     = neo_dke_out(is,5)    ! v_theta/vt_norm (at theta=0)
     neo_vtor_dke_out(is)     = neo_dke_out(is,6)    ! v_phi/vt_norm (at theta=0)
  enddo
  neo_jpar_dke_out = neo_dke_1d_out                  ! <j B> / j_norm

  ! gyro-viscosity results
  do is=1,11
     neo_pflux_gv_out(is)     = neo_gv_out(is,1)     ! Gamma_is/Gamma_norm
     neo_efluxtot_gv_out(is)  = neo_gv_out(is,2)     ! Q_is/Q_norm total
     neo_mflux_gv_out(is)     = neo_gv_out(is,3)     ! Pi_is/Pi_norm
     neo_efluxncv_gv_out(is)  = neo_gv_out(is,4)     ! Q_is/Q_norm non-convective
  enddo

  neo_error_status_out  = error_status
  neo_error_message_out = error_message
  call neo_harvest_output
end subroutine neo_run
