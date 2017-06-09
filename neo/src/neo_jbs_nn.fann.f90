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
     

     IMPLICIT NONE

     character(len=1000) :: jbsnn_model

     integer :: ierr
     real :: K_0, K_e, K_i1, Ki2
     
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
    
     INPUT_PARAMETERS(4)=nn_rmin_in                                ! rmin
     INPUT_PARAMETERS(3)=nn_q_in                                   ! q
     INPUT_PARAMETERS(2)=log10(nn_nuee_in)                         ! nuee nor
     INPUT_PARAMETERS(1)=nn_ni1_ne_in                              ! ni1/ne
     INPUT_PARAMETERS(5)=nn_ti1_te_in                              ! ti1/te
    
     ! INPUT_PARAMETERS( 6)=neo_mass_in(2)                          ! mi2/mD
     ! INPUT_PARAMETERS( 7)=neo_z_in(2)                             ! zi2
     ! parameters mi2/mD are to be added when using model with random Z impurity...
    

     WRITE(*,*)INPUT_PARAMETERS



     !run NN

     call get_environment_variable('JBSNN_MODEL_DIR',jbsnn_model)
    
     ierr=load_anns(0, TRIM(jbsnn_model)//NUL,'brainfuse'//NUL)
     ierr=load_anns_inputs(INPUT_PARAMETERS)
     ierr=run_anns()

    
     ierr=get_anns_avg_array(OUTPUT_PARAMETERS)
    
     OUT_CNEO_CnC= OUTPUT_PARAMETERS(1)
     OUT_CNEO_CnD= OUTPUT_PARAMETERS(2)
     OUT_CNEO_Cne= OUTPUT_PARAMETERS(3)
     OUT_CNEO_CTC= OUTPUT_PARAMETERS(4)
     OUT_CNEO_CTD= OUTPUT_PARAMETERS(5)
     OUT_CNEO_CTe= OUTPUT_PARAMETERS(6)
     OUT_CSAU_CnC= OUTPUT_PARAMETERS(7)
     OUT_CSAU_CnD= OUTPUT_PARAMETERS(8)
     OUT_CSAU_Cne= OUTPUT_PARAMETERS(9)
     OUT_CSAU_CTC= OUTPUT_PARAMETERS(10)
     OUT_CSAU_CTD= OUTPUT_PARAMETERS(11)
     OUT_CSAU_CTe= OUTPUT_PARAMETERS(12)


     WRITE(*,*)OUTPUT_PARAMETERS

     !!!!!!!!! Bootsrap Current Calculation !!!!!!!!!!!!!
     
     ! K_0 in MA.m^2
     K_0=charge_norm_fac*nn_vnorm*nn_anorm*nn_I_over_phi_prime*nn_rho_star/1e6
     
     !K_e in 1/m^4
     K_e_NEO=abs(zfac(3))*dens_norm*((OUT_CNEO_CTe*nn_1_over_Lte)+(OUT_CNEO_Cne*nn_1_over_Lne))
     K_e_SAU=abs(zfac(3))*dens_norm*((OUT_CSAU_CTe*nn_1_over_Lte)+(OUT_CSAU_Cne*nn_1_over_Lne))
     
     ! likewise for K_i1, Ki2
     K_i1_NEO=zfac(1)*nn_ni1_dens*((OUT_CNEO_CTD*nn_1_over_LtD)+(OUT_CNEO_CnD*nn_1_over_LnD))
     K_i1_SAU=zfac(1)*nn_ni1_dens*((OUT_CSAU_CTD*nn_1_over_LtD)+(OUT_CSAU_CnD*nn_1_over_LnD))
     
     K_i2_NEO=zfac(2)*nn_ni2_dens*((OUT_CNEO_CTC*nn_1_over_LtC)+(OUT_CNEO_CnC*nn_1_over_LnC))
     K_i2_SAU=zfac(2)*nn_ni2_dens*((OUT_CSAU_CTC*nn_1_over_LtC)+(OUT_CSAU_CnC*nn_1_over_LnC))
     
     
     
        
        !!!!!OUTPUT!!!!
        
    !neo_th_out(5) = sauter ! Need to fill in sum over terms

     nn_NEO_jbs_in_A_m2 = K_0*(K_e_NEO+K_i1_NEO+K_i2_NEO)             !jbs_neo  in (A/m^2)  Need to fill in sum over terms
     nn_SAU_jbs_in_A_m2 = K_0*(K_e_SAU+K_i1_SAU+K_i2_SAU)             !jbs_sau  in (A/m^2)  Need to fill in sum over terms


     neo_dke_1d_out(1) = nn_NEO_jbs_in_A_m2/jbs_norm                  ! dimensionless again to be consistent with full NEO.
     neo_th_out(5) = nn_SAU_jbs_in_A_m2/jbs_norm  
     
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
