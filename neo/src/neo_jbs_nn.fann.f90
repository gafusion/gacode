!---------------------------------------------------------------------------------------------------------------------------
! neo_jbs_nn 
!
! PURPOSE: 
!  Driver for the  (bootstrap current generator ) capability of NEO-NN.
!  This at long therm   will write a new input.profiles with NEO-NN-computed bootstrap current and sauter current
!  The generator for the moment call and run the artificial neural network and calculates the coefficients for the Currents
! 
!----------------------------------------------------------------------------------------------------------------------------


     SUBROUTINE neo_jbs_nn
!
!  Execution of the NN
!
     use mpi
     use neo_globals
     use neo_interface
     use EXPRO_interface
     use jbsnn_globals
     
!
     IMPLICIT NONE
!
     character(len=1000) :: jbsnn_model

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
     real :: OUT_CNEO_CTD,OUT_CNEO_CTe,OUT_CNEO_CTi,OUT_CNEO_CnD,OUT_CNEO_Cne,OUT_CNEO_Cni
     real :: OUT_CSAU_CTD,OUT_CSAU_CTe,OUT_CSAU_CTi,OUT_CSAU_CnD,OUT_CSAU_Cne,OUT_SAU_Cni

     real(4) :: INPUT_PARAMETERS(5)
     real(4) :: OUTPUT_PARAMETERS(12)
!    real :: start, finish
     integer :: n, i

     CHARACTER NUL
     PARAMETER(NUL = CHAR(0))

     include 'brainfuse_lib.inc'




!fill in input parameters array


  ! Read jbsnn control parameters

    
  
     open(unit=1,file='jbsnn.dat',status='old')
     
     read(1,*) impurities_spec
     read(1,*) plasma_spec

     close(1)

     select case(impurities_spec)
     case(1)
     
        print '(a)','INFO: (JBSNN) Computing Bootstrap Current  (JbsNN) with Carbon impurity'
    
        INPUT_PARAMETERS( 1)=neo_rmin_over_a_in    ! rmin
        INPUT_PARAMETERS( 2)=neo_q_in              ! q
        INPUT_PARAMETERS( 3)=neo_nu_1_in           ! nuee
        INPUT_PARAMETERS( 4)=neo_dens_in(1)        ! ni1/ne
        INPUT_PARAMETERS( 5)=neo_temp_in(1)        ! ti1/te
    
       ! INPUT_PARAMETERS( 6)=neo_mass_in(2)       ! mi2/mD
       ! INPUT_PARAMETERS( 7)=neo_z_in(2)          ! zi2
    

        WRITE(*,*)INPUT_PARAMETERS



        !run NN

        call get_environment_variable('JBSNN_MODEL_DIR',jbsnn_model)
    
        ierr=load_anns(0, TRIM(jbsnn_model)//NUL,'brainfuse'//NUL)
        ierr=load_anns_inputs(INPUT_PARAMETERS)
        ierr=run_anns()

    
        ierr=get_anns_avg_array(OUTPUT_PARAMETERS)
    
        OUT_CNEO_Cne= OUTPUT_PARAMETERS(1)
        OUT_CNEO_CnD= OUTPUT_PARAMETERS(2)
        OUT_CNEO_Cni= OUTPUT_PARAMETERS(3)
        OUT_CNEO_CTe= OUTPUT_PARAMETERS(4)
        OUT_CNEO_CTD= OUTPUT_PARAMETERS(5)
        OUT_CNEO_CTi= OUTPUT_PARAMETERS(6)
        OUT_CSAU_Cne= OUTPUT_PARAMETERS(7)
        OUT_CSAU_CnD= OUTPUT_PARAMETERS(8)
        OUT_CSAU_Cni= OUTPUT_PARAMETERS(9)
        OUT_CSAU_CTe= OUTPUT_PARAMETERS(10)
        OUT_CSAU_CTD= OUTPUT_PARAMETERS(11)
        OUT_CSAU_CTi= OUTPUT_PARAMETERS(12)

        WRITE(*,*)OUTPUT_PARAMETERS

       !open(unit=1,file='input.profiles.jbs',status='replace')
       ! input.profiles.jbs is to be filled by computed NN computed jbs. 
       !write(1,'(a)') '#'
       ! write(1,'(a)') '# expro_rho
       !write(1,'(a)') '# jbs_neo_nn    (MA/m^2)'
       !write(1,'(a)') '# jbs_sauter (MA/m^2)'
   
       !write(1,'(a)') '#'


    ! do i=1,EXPRO_n_exp
    !    write(1,'(6(1pe14.7,2x))') EXPRO_rho(i), pflux_sum(i), &
    !         jbs_neo_nn(i), jbs_sauter_nn(i),
    ! enddo

 
     close(1)
        

     END SUBROUTINE neo_jbs_nn
