!-----------------------------------------------------------------
!
     SUBROUTINE tglf_nn_TM
!
!  Execution of the NN
!
     use tglf_pkg
     use tglf_interface
     USE tglf_dimensions
     USE tglf_global
!
     IMPLICIT NONE
!
     character(len=1000) :: nn_executable
     character(len=1000) :: nn_files

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
     real :: INPUT_PARAMETERS(23)
     integer :: n, i

     CHARACTER NUL
     PARAMETER(NUL = CHAR(0))

     include 'brainfuse_lib.inc'
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
!fill in input parameters array
!
    INPUT_PARAMETERS( 1)=tglf_as_in(2)         ! AS_2
    INPUT_PARAMETERS( 2)=tglf_as_in(3)         ! AS_3
    INPUT_PARAMETERS( 3)=tglf_betae_in         ! BETAE
    INPUT_PARAMETERS( 4)=tglf_delta_loc_in     ! DELTA_LOC
    INPUT_PARAMETERS( 5)=tglf_drmajdx_loc_in   ! DRMAJDX_LOC
    INPUT_PARAMETERS( 6)=tglf_kappa_loc_in     ! KAPPA_LOC
    INPUT_PARAMETERS( 7)=tglf_p_prime_loc_in   ! P_PRIME_LOC
    INPUT_PARAMETERS( 8)=tglf_q_loc_in         ! Q_LOC
    INPUT_PARAMETERS( 9)=tglf_q_prime_loc_in   ! Q_PRIME_LOC
    INPUT_PARAMETERS(10)=tglf_rlns_in(1)       ! RLNS_1
    INPUT_PARAMETERS(11)=tglf_rlns_in(2)       ! RLNS_2
    INPUT_PARAMETERS(12)=tglf_rlns_in(3)       ! RLNS_3
    INPUT_PARAMETERS(13)=tglf_rlts_in(1)       ! RLTS_1
    INPUT_PARAMETERS(14)=tglf_rlts_in(2)       ! RLTS_2
    INPUT_PARAMETERS(15)=tglf_rmaj_loc_in      ! RMAJ_LOC
    INPUT_PARAMETERS(16)=tglf_rmin_loc_in      ! RMIN_LOC
    INPUT_PARAMETERS(17)=tglf_s_kappa_loc_in   ! S_KAPPA_LOC
    INPUT_PARAMETERS(18)=tglf_taus_in(2)       ! TAUS_2
    INPUT_PARAMETERS(19)=tglf_vexb_shear_in    ! VEXB_SHEAR
    INPUT_PARAMETERS(20)=tglf_vpar_in(1)       ! VPAR_1
    INPUT_PARAMETERS(21)=tglf_vpar_shear_in(1) ! VPAR_SHEAR_1
    INPUT_PARAMETERS(22)=tglf_xnue_in          ! XNUE
    INPUT_PARAMETERS(23)=tglf_zeff_in          ! ZEFF

!    WRITE(*,*)INPUT_PARAMETERS

    call get_environment_variable('TGLFNN_MODEL',nn_files)
    ierr=load_anns(TRIM(nn_files)//NUL,'brainfuse'//NUL)
    ierr=run_anns()

    energy_flux_out(1,1)   = get_anns_avg(0)
    energy_flux_out(3,1)   = get_anns_avg(1)
    particle_flux_out(1,1) = get_anns_avg(2)
    stress_tor_out(3,1)    = get_anns_avg(3)

!    write(*,*)    energy_flux_out(1,1),   &
!                  energy_flux_out(3,1),   &
!                  particle_flux_out(1,1), &
!                  stress_tor_out(3,1)

!
! '#---------------------------------------------------'
! '# Sort ions by A, Z, Te/Ti'
! '#---------------------------------------------------'
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

    if ( (tglf_mass_in(2) .ne. 1) .or. &
         (tglf_zs_in(2)   .ne. 1) .or. &
         (tglf_mass_in(3) .ne. 6) .or. &
         (tglf_zs_in(3)   .ne. 6) ) then
        write(*,*)'A',tglf_mass_in(:tglf_ns_in)
        write(*,*)'Z',tglf_zs_in(:tglf_ns_in)
        write (*,*)'NN trained only with D C ions'
        stop
    endif

    OUT_ENERGY_FLUX_i_RNG=0
    OUT_ENERGY_FLUX_1_RNG=0
    OUT_PARTICLE_FLUX_1_RNG=0
    OUT_STRESS_TOR_i_RNG=0

!
! Switch between TGLF and the NN depending on 'nn_max_error' values
!
     if ((OUT_ENERGY_FLUX_1_RNG   < tglf_nn_max_error_in) .and. &
         (OUT_ENERGY_FLUX_i_RNG   < tglf_nn_max_error_in) .and. &
         (OUT_PARTICLE_FLUX_1_RNG < tglf_nn_max_error_in) .and. &
         (OUT_STRESS_TOR_i_RNG    < tglf_nn_max_error_in)) then
        valid_nn = .TRUE.
        write(*,*) '============>    NN    RUNNING    <=============='
     else
        write(*,*) '------------>   TGLF    RUNNING   <--------------'
     endif
!
      END SUBROUTINE tglf_nn_TM
