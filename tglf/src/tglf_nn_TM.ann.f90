!-----------------------------------------------------------------
!
     SUBROUTINE tglf_nn_TM
!
!  Execution of the NN
!
     use tglf_pkg
     use tglf_interface
     use tglf_dimensions
     use tglf_global
     use tglf_mpi
!
     IMPLICIT NONE
!
     character(len=1000) :: tglfnn_model

     integer :: j,is,ierr
     real :: v_bar0, phi_bar0
     real :: pflux0(nsm,3),eflux0(nsm,3)
     real :: stress_par0(nsm,3),stress_tor0(nsm,3)
     real :: exch0(nsm,3)
     real :: nsum0(nsm),tsum0(nsm)
!
     real, DIMENSION(10) :: tmp
     integer, DIMENSION(10) :: ions_order
!
     real :: OUT_ENERGY_FLUX_1_RNG, OUT_ENERGY_FLUX_i_RNG
     real :: OUT_PARTICLE_FLUX_1_RNG, OUT_STRESS_TOR_i_RNG
     real :: INPUT_PARAMETERS_2IONS(24)
     real :: OUTPUT_PARAMETERS_2IONS(6)
     real :: INPUT_PARAMETERS_3IONS(26)
     real :: OUTPUT_PARAMETERS_3IONS(7)

!     real :: start, finish
     integer :: n, i

     CHARACTER NUL
     PARAMETER(NUL = CHAR(0))

     include 'brainfusetf_lib.inc'

!     call cpu_time(start)
!
! initialize fluxes
!
     do is=1,ns
       do j=1,3 
          particle_flux_out(is,j) = 0.0
          energy_flux_out(is,j) = 0.0
          stress_par_out(is,j) = 0.0
          stress_tor_out(is,j) = 0.0
          exchange_out(is,j) = 0.0
          pflux0(is,j) = 0.0
          eflux0(is,j) = 0.0
          stress_par0(is,j) = 0.0
          stress_tor0(is,j) = 0.0
          exch0(is,j) = 0.0
       enddo
       q_low_out(is) = 0.0
     enddo

!
!sort ions by A, Z, Te/Ti
!electrons always first (as required by TGLF)
!main ion always second (affects GB normalization)
!
     tmp=0.0
     tmp(1)=1E10
     tmp(2)=1E9
     DO i = 3,tglf_ns_in
         tmp(i)=tglf_mass_in(i)*100000+tglf_zs_in(i)*1000+1./(1+tglf_rlts_in(i))
     ENDDO
     DO i = 1,tglf_ns_in
         j=MAXLOC(tmp, DIM=1)
         ions_order(i)=j
         tmp(j)=0
     ENDDO

!
!fill in input parameters array
!
    if      ( (tglf_ns_in.eq.3 ) .and. &
              (abs(tglf_mass_in(2)-1)   .lt. 1E-3) .and. &
              (abs(tglf_zs_in(2)-1)     .lt. 1E-3) .and. &
              (abs(tglf_mass_in(3)-6)   .lt. 1E-3) .and. &
              (abs(tglf_zs_in(3)-6)     .lt. 1E-3) )then

       INPUT_PARAMETERS_2IONS( 1)=tglf_as_in(2)         ! AS_2
       INPUT_PARAMETERS_2IONS( 2)=tglf_as_in(3)         ! AS_3
       INPUT_PARAMETERS_2IONS( 3)=tglf_betae_in         ! BETAE
       INPUT_PARAMETERS_2IONS( 4)=tglf_debye_in         ! DEBYE
       INPUT_PARAMETERS_2IONS( 5)=tglf_delta_loc_in     ! DELTA_LOC
       INPUT_PARAMETERS_2IONS( 6)=tglf_drmajdx_loc_in   ! DRMAJDX_LOC
       INPUT_PARAMETERS_2IONS( 7)=tglf_kappa_loc_in     ! KAPPA_LOC
       INPUT_PARAMETERS_2IONS( 8)=tglf_p_prime_loc_in   ! P_PRIME_LOC
       INPUT_PARAMETERS_2IONS( 9)=tglf_q_loc_in         ! Q_LOC
       INPUT_PARAMETERS_2IONS(10)=tglf_q_prime_loc_in   ! Q_PRIME_LOC
       INPUT_PARAMETERS_2IONS(11)=tglf_rlns_in(1)       ! RLNS_1
       INPUT_PARAMETERS_2IONS(12)=tglf_rlns_in(2)       ! RLNS_2
       INPUT_PARAMETERS_2IONS(13)=tglf_rlns_in(3)       ! RLNS_3
       INPUT_PARAMETERS_2IONS(14)=tglf_rlts_in(1)       ! RLTS_1
       INPUT_PARAMETERS_2IONS(15)=tglf_rlts_in(2)       ! RLTS_2
       INPUT_PARAMETERS_2IONS(16)=tglf_rmaj_loc_in      ! RMAJ_LOC
       INPUT_PARAMETERS_2IONS(17)=tglf_rmin_loc_in      ! RMIN_LOC
       INPUT_PARAMETERS_2IONS(18)=tglf_s_kappa_loc_in   ! S_KAPPA_LOC
       INPUT_PARAMETERS_2IONS(19)=tglf_taus_in(2)       ! TAUS_2
       INPUT_PARAMETERS_2IONS(20)=tglf_vexb_shear_in    ! VEXB_SHEAR
       INPUT_PARAMETERS_2IONS(21)=tglf_vpar_in(1)       ! VPAR_1
       INPUT_PARAMETERS_2IONS(22)=tglf_vpar_shear_in(1) ! VPAR_SHEAR_1
       INPUT_PARAMETERS_2IONS(23)=tglf_xnue_in          ! XNUE
       INPUT_PARAMETERS_2IONS(24)=tglf_zeff_in          ! ZEFF

       call get_environment_variable('TGLFNN_MODEL_2IONS',tglfnn_model)
       ierr=btf_run(TRIM(tglfnn_model)//NUL, &
            INPUT_PARAMETERS_2IONS, SIZE(INPUT_PARAMETERS_2IONS), &
            OUTPUT_PARAMETERS_2IONS, SIZE(OUTPUT_PARAMETERS_2IONS))

       energy_flux_out(1,1)   = OUTPUT_PARAMETERS_2IONS(1)
       energy_flux_out(2,1)   = OUTPUT_PARAMETERS_2IONS(2)/2.
       energy_flux_out(3,1)   = OUTPUT_PARAMETERS_2IONS(2)/2.
       particle_flux_out(1,1) = OUTPUT_PARAMETERS_2IONS(3)
       particle_flux_out(2,1) = OUTPUT_PARAMETERS_2IONS(4)
       particle_flux_out(3,1) = OUTPUT_PARAMETERS_2IONS(5)
       stress_tor_out(3,1)    = OUTPUT_PARAMETERS_2IONS(6)

    else if ( (tglf_ns_in.eq.4 ) .and. &
              (abs(tglf_mass_in(2)-1)   .lt. 1E-3) .and. &
              (abs(tglf_zs_in(2)-1)     .lt. 1E-3) .and. &
              (abs(tglf_mass_in(3)-1.6) .lt. 1E-3) .and. &
              (abs(tglf_zs_in(3)-2)     .lt. 1E-3) .and. &
              (abs(tglf_mass_in(4)-8)   .lt. 1E-3) .and. &
              (abs(tglf_zs_in(4)-10)    .lt. 1E-3) )then

       INPUT_PARAMETERS_3IONS( 1)=tglf_as_in(2)         ! AS_2
       INPUT_PARAMETERS_3IONS( 2)=tglf_as_in(4)         ! AS_3 # 3<->4
       INPUT_PARAMETERS_3IONS( 3)=tglf_as_in(3)         ! AS_4 # 3<->4
       INPUT_PARAMETERS_3IONS( 4)=tglf_betae_in         ! BETAE
       INPUT_PARAMETERS_3IONS( 5)=tglf_debye_in         ! DEBYE
       INPUT_PARAMETERS_3IONS( 6)=tglf_delta_loc_in     ! DELTA_LOC
       INPUT_PARAMETERS_3IONS( 7)=tglf_drmajdx_loc_in   ! DRMAJDX_LOC
       INPUT_PARAMETERS_3IONS( 8)=tglf_kappa_loc_in     ! KAPPA_LOC
       INPUT_PARAMETERS_3IONS( 9)=tglf_p_prime_loc_in   ! P_PRIME_LOC
       INPUT_PARAMETERS_3IONS(10)=tglf_q_loc_in         ! Q_LOC
       INPUT_PARAMETERS_3IONS(11)=tglf_q_prime_loc_in   ! Q_PRIME_LOC
       INPUT_PARAMETERS_3IONS(12)=tglf_rlns_in(1)       ! RLNS_1
       INPUT_PARAMETERS_3IONS(13)=tglf_rlns_in(2)       ! RLNS_2
       INPUT_PARAMETERS_3IONS(14)=tglf_rlns_in(4)       ! RLNS_3 # 3<->4
       INPUT_PARAMETERS_3IONS(15)=tglf_rlns_in(3)       ! RLNS_4 # 3<->4
       INPUT_PARAMETERS_3IONS(16)=tglf_rlts_in(1)       ! RLTS_1
       INPUT_PARAMETERS_3IONS(17)=tglf_rlts_in(2)       ! RLTS_2
       INPUT_PARAMETERS_3IONS(18)=tglf_rmaj_loc_in      ! RMAJ_LOC
       INPUT_PARAMETERS_3IONS(19)=tglf_rmin_loc_in      ! RMIN_LOC
       INPUT_PARAMETERS_3IONS(20)=tglf_s_kappa_loc_in   ! S_KAPPA_LOC
       INPUT_PARAMETERS_3IONS(21)=tglf_taus_in(2)       ! TAUS_2
       INPUT_PARAMETERS_3IONS(22)=tglf_vexb_shear_in    ! VEXB_SHEAR
       INPUT_PARAMETERS_3IONS(23)=tglf_vpar_in(1)       ! VPAR_1
       INPUT_PARAMETERS_3IONS(24)=tglf_vpar_shear_in(1) ! VPAR_SHEAR_1
       INPUT_PARAMETERS_3IONS(25)=tglf_xnue_in          ! XNUE
       INPUT_PARAMETERS_3IONS(26)=tglf_zeff_in          ! ZEFF

       call get_environment_variable('TGLFNN_MODEL_3IONS',tglfnn_model)
       ierr=btf_run(TRIM(tglfnn_model)//NUL, &
            INPUT_PARAMETERS_3IONS, SIZE(INPUT_PARAMETERS_3IONS), &
            OUTPUT_PARAMETERS_3IONS, SIZE(OUTPUT_PARAMETERS_3IONS))

       energy_flux_out(1,1)   = OUTPUT_PARAMETERS_3IONS(1)
       energy_flux_out(2,1)   = OUTPUT_PARAMETERS_3IONS(2)/3.
       energy_flux_out(3,1)   = OUTPUT_PARAMETERS_3IONS(2)/3.
       energy_flux_out(4,1)   = OUTPUT_PARAMETERS_3IONS(2)/3.
       particle_flux_out(1,1) = OUTPUT_PARAMETERS_3IONS(3)
       particle_flux_out(2,1) = OUTPUT_PARAMETERS_3IONS(4)
       particle_flux_out(3,1) = OUTPUT_PARAMETERS_3IONS(6) !# 3<->4
       particle_flux_out(4,1) = OUTPUT_PARAMETERS_3IONS(5) !# 3<->4
       stress_tor_out(2,1)    = OUTPUT_PARAMETERS_3IONS(7)
    else

        write(*,*)'A',tglf_mass_in(2:tglf_ns_in)
        write(*,*)'Z',tglf_zs_in(2:tglf_ns_in)
        write (*,*)'NN trained only with D,C or DT,He4,Ne ions'
        stop
    endif

!
!switch between TGLF and the NN depending on 'nn_max_error' values (disabled)
     OUT_ENERGY_FLUX_i_RNG=0
     OUT_ENERGY_FLUX_1_RNG=0
     OUT_PARTICLE_FLUX_1_RNG=0
     OUT_STRESS_TOR_i_RNG=0

     if ((OUT_ENERGY_FLUX_1_RNG   < tglf_nn_max_error_in) .and. &
         (OUT_ENERGY_FLUX_i_RNG   < tglf_nn_max_error_in) .and. &
         (OUT_PARTICLE_FLUX_1_RNG < tglf_nn_max_error_in) .and. &
         (OUT_STRESS_TOR_i_RNG    < tglf_nn_max_error_in)) then
        valid_nn = .TRUE.
        !write(*,*) '============>    NN    RUNNING    <=============='
     else
        !write(*,*) '------------>   TGLF    RUNNING   <--------------'
     endif
!
!     call cpu_time(finish)
!     print '("Time = ",f12.6," seconds.")',finish-start

      END SUBROUTINE tglf_nn_TM
