!---------------------------------------------------------------------------------------------------------------------------
! neo_jbs_nn 
!
! PURPOSE: 
!  Driver for the  (bootstrap current generator ) capability of NEO-NN.
!  This at long therm   will write a new input.profiles with NEO-NN-computed bootstrap and sauter  current and neoclassical fuxes
!  The generator for the moment call and run the artificial neural network and calculates the bootstrap Current, from coefficients
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
     use vgen_globals
!    use jbsnn_globals, To use Later when NN will compute fluxes, rotation,...
     
!
     IMPLICIT NONE
!
     character(len=1000) :: jbsnn_model

     integer :: ierr
     real :: K_0, K_e, K_i1, Ki2
     
     
!
     real :: OUT_CNEO_CTD,OUT_CNEO_CTe,OUT_CNEO_CTi,OUT_CNEO_CnD,OUT_CNEO_Cne,OUT_CNEO_Cni
     real :: OUT_CSAU_CTD,OUT_CSAU_CTe,OUT_CSAU_CTi,OUT_CSAU_CnD,OUT_CSAU_Cne,OUT_SAU_Cni

     real(4) :: INPUT_PARAMETERS(5)
     real(4) :: OUTPUT_PARAMETERS(12)
     real :: start, finish
    

     CHARACTER NUL
     PARAMETER(NUL = CHAR(0))

     include 'brainfuse_lib.inc'

!fill in input parameters array

     call cpu_time(start)
     print '(a)','INFO: (JBSNN) Computing Bootstrap Current with NEO-NN'
    
     INPUT_PARAMETERS( 1)=neo_rmin_over_a_in/ neo_rmaj_over_a_in     ! rmin
     INPUT_PARAMETERS( 2)=neo_q_in                                   ! q
     INPUT_PARAMETERS( 3)=neo_nu_1_in                                ! nuee nor
     INPUT_PARAMETERS( 4)=neo_dens_in(1)                             ! ni1/ne
     INPUT_PARAMETERS( 5)=neo_temp_in(1)                             ! ti1/te
    
     ! INPUT_PARAMETERS( 6)=neo_mass_in(2)                           ! mi2/mD
     ! INPUT_PARAMETERS( 7)=neo_z_in(2)                              ! zi2
     ! parameters mi2/mD are to be added when using model with random Z impurity...
    

     WRITE(*,*)INPUT_PARAMETERS



     !run NN

     call get_environment_variable('JBSNN_MODEL_DIR',jbsnn_model)
    
     ierr=load_anns(0, TRIM(jbsnn_model)//NUL,'brainfuse'//NUL)
     ierr=load_anns_inputs(INPUT_PARAMETERS)
     ierr=run_anns()

    
     ierr=get_anns_avg_array(OUTPUT_PARAMETERS)
    
     OUT_CNEO_Cne= OUTPUT_PARAMETERS(1)
     OUT_CNEO_CnD= OUTPUT_PARAMETERS(2)
     OUT_CNEO_CnC= OUTPUT_PARAMETERS(3)
     OUT_CNEO_CTe= OUTPUT_PARAMETERS(4)
     OUT_CNEO_CTD= OUTPUT_PARAMETERS(5)
     OUT_CNEO_CTC= OUTPUT_PARAMETERS(6)
     OUT_CSAU_Cne= OUTPUT_PARAMETERS(7)
     OUT_CSAU_CnD= OUTPUT_PARAMETERS(8)
     OUT_CSAU_CnC= OUTPUT_PARAMETERS(9)
     OUT_CSAU_CTe= OUTPUT_PARAMETERS(10)
     OUT_CSAU_CTD= OUTPUT_PARAMETERS(11)
     OUT_CSAU_CTC= OUTPUT_PARAMETERS(12)


     WRITE(*,*)OUTPUT_PARAMETERS

     !!!!!!!!! Bootsrap Current Calculation !!!!!!!!!!!!!
     
     ! K_0 in A.m^2
     K_0=charge_norm_fac*nn_vnorm*nn_anorm*nn_I_over_phi_prime*nn_rho_star
     !K_e in 1/m^2
     K_e=abs(zfac(3))*dens_norm*((OUT_CNEO_CTe*nn_1_over_Lte)+(OUT_CNEO_Cne*nn_1_over_Lne))
     ! likewise for K_i1, Ki2
     K_i1=zfac(1)*nn_ni1_dens*((OUT_CNEO_CTD*nn_1_over_LtD)+(OUT_CNEO_CnD*nn_1_over_LnD))
     K_i2=zfac(2)*nn_ni2_dens*((OUT_CNEO_CTC*nn_1_over_LtC)+(OUT_CNEO_CnC*nn_1_over_LnC))
     
     
     
        
        !!!!!OUTPUT!!!!
        
    !neo_th_out(5) = sauter ! Need to fill in sum over terms

     neo_dke_1d_out_nn_normalized = K_0*(K_e+K_i1+K_i2)    !jbs_neo  in (A)  Need to fill in sum over terms

     neo_dke_1d_out = neo_dke_1d_out_nn_normalized/jbs_norm 

     
        !!!!!!!!!!!!!!!
        
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

 
     !close(1)
        
     call cpu_time(finish)
     print '("Time = ",f12.6," seconds.")',finish-start

     
   END SUBROUTINE neo_jbs_nn
