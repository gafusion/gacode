    MODULE nfreya_run

      USE nrtype,                                   ONLY : DP,I4B

      USE Nfreya_routines,                          ONLY : freyas

      USE ions_gcnmp,                               ONLY : atw

      USE neutral_beams,                            ONLY : npart_all_beamlines,   &
                                                           no_injectors,iexcit,   &
                                                           beam_sim_time_start,   &
                                                           beam_sim_time_end,     &
                                                           time_now,              &
                                                           no_physical_injectors, &
                                                           pseudo_injectors,      &
                                                           npart_pseudo,          &
                                                           no_birth_points,       &
                                                           nfreya_vb,             &
                                                           master_injector_no,    &
                                                           inj_list,              &
                                                           injctr_wk_list,        &
                                                           phys_injctr_wk_list

      USE tport_mhd_grid_data,                      ONLY : psi_g_cm2,rmhdgrid,    &
                                                           zmhdgrid,ymagn1,       &
                                                           zplasmin,zplasmax,     &
                                                           rplasmax

      USE common_constants,                         ONLY : izero,zeroc,cm2m

      USE error_handler,                            ONLY : terminate,iomaxerr,lerrno

      USE io_gcnmp,                                 ONLY : ncrt,nlog



#if defined USEMPI
      USE mpi

#endif

      IMPLICIT NONE

      CHARACTER *36 codeid

      CONTAINS

      SUBROUTINE simulation_time_advance
!------------------------------------------------------------------------------------
!
!----------------------------------------------------------------HSJ-----------------
       USE MPI_data,                   ONLY : numprocs
        IMPLICIT NONE 

        time_now = beam_sim_time_start

        IF(numprocs  == 1)THEN 

           CALL nfreya_single_proc 

        ELSE                  ! user requested more than 1 cpu,run  beamlines in parallel 
#ifdef USEMPI                              ! using master/slave paradigm
           CALL nfreya_master_slave
#endif
        ENDIF 

      END SUBROUTINE simulation_time_advance

      SUBROUTINE nfreya_single_proc
!-----------------------------------------------------------------
! -- Use this routine to run nfreya in single processor mode
! -- But note that compilation/linking is always done with MPI
! -- We pass all injectors to freyas at onece
! -- and let freyas loop over the injectors internally
! -- the number of pseudo neutrals to be followed into the 
! -- plasma is given by npart_pseudo(ie,ib) in all cases
! -- (npart_pseudo was preloaded in sub beam_prop)
! -- npart_all_beamlines is used in freyas only to check bounds
! -- on output arrays
!-----------------------------------------------------------------
     USE nf_param,                                       ONLY : kbe,ke,kb

      IMPLICIT NONE

      INTEGER(I4B) injs_to_do,injctr_id

           injs_to_do = no_injectors
           injctr_id  = izero                  ! not relevent in this case

           CALL freyas(injs_to_do, injctr_id, codeid,            &
                   psi_g_cm2, rmhdgrid,rmhdgrid(1), rplasmax,    &
                   zmhdgrid,ymagn1, zplasmin,zplasmax, iexcit,   &
                   time_now,beam_sim_time_start)

          RETURN
      END SUBROUTINE nfreya_single_proc




      SUBROUTINE nfreya_update_master
!----------------------------------------------------------------------------------
! -- THIS ROUTINE IS NO LONGER USED ! KEEP it here just in case Iwant it back
! --   Put data generated on master into proper slots for output
! --   This is done by routine nfreya_store_packed for the slaves
! --   NOTE: ROUTINES THAT MUST BE KEPT IN SYNC:
! --       nfreya_pack_data
! --       nfreya_update_master
! --       nfreya_store_packed 
! --       allocate_mpi_buffer
!-----------------------------------------------------------------------------------

        USE nrtype,                                   ONLY :  DP,I4B

        USE neutral_beams,                            ONLY :                        &
                                                             master_injector_no,    &
                                                             inj_list,              &
                                                             no_physical_injectors, &
                                                             no_injectors,          &
                                                             pseudo_injectors,      &
                                                             nmbrz,vx_izpt,vy_izpt, &
                                                             vz_izpt,x_izpt,y_izpt, &
                                                             z_izpt,r_izpt,pitch_a, &
                                                             bion,bneut,pbeam,fap,  &
                                                             fwall,forb,fb11,fb00,  &
                                                             fb01,fb10,wb11,wb10,   &
                                                             wb01,wb00,ftrapfi,     &
                                                             hibrz,hdepz,zetaz,     &
                                                             angmpz,fber,hicmz


        IMPLICIT NONE

        INTEGER(I4B) index,k,l,j,ir,ic,i3

             IF(master_injector_no .GT. no_physical_injectors)THEN
                index = pseudo_injectors(master_injector_no - no_physical_injectors)
             ELSE
                index = master_injector_no
             ENDIF
    
             nmbrz(:,index)          =  nmbrz(:,index)  +  nmbrz(:,master_injector_no)


             vx_izpt(:,:,index)      = vx_izpt(:,:,master_injector_no) + vx_izpt(:,:,index)
             vy_izpt(:,:,index)      = vy_izpt(:,:,master_injector_no) + vy_izpt(:,:,index)
             vz_izpt(:,:,index)      = vz_izpt(:,:,master_injector_no) + vz_izpt(:,:,index)

             x_izpt(:,:,index)       = x_izpt(:,:,master_injector_no)  + x_izpt(:,:,index)
             y_izpt(:,:,index)       = y_izpt(:,:,master_injector_no)  + y_izpt(:,:,index)
             z_izpt(:,:,index)       = z_izpt(:,:,master_injector_no)  + z_izpt(:,:,index)
             r_izpt(:,:,index)       = r_izpt(:,:,master_injector_no)  + r_izpt(:,:,index)

             pitch_a(:,:,index)      = pitch_a(:,:,master_injector_no) + pitch_a(:,:,index)


             bion(:,index)           = bion(:,master_injector_no)  !+ bion(:,index) check ??
             bneut(:,index)          = bneut(:,master_injector_no) !+ bneut(:,index)
             pbeam(:,index)          = pbeam(:,master_injector_no) !+ pbeam(:,index)

             fap(:,index)            = fap(:,master_injector_no)   + fap(:,index)   
             fwall(:,index)          = fwall(:,master_injector_no) + fwall(:,index)  
             forb(:,index)           = forb(:,master_injector_no)  + forb(:,index)

             fb11(:,index)           =  fb11(:,master_injector_no) + fb11(:,index) 
             fb10(:,index)           =  fb10(:,master_injector_no) + fb10(:,index) 
             fb01(:,index)           =  fb01(:,master_injector_no) + fb01(:,index) 
             fb00(:,index)           =  fb00(:,master_injector_no) + fb00(:,index) 
             fber(:,index)           =  fber(:,master_injector_no) + fber(:,index)

             wb11(:,index)           =  wb11(:,master_injector_no) + wb11(:,index)
             wb10(:,index)           =  wb10(:,master_injector_no) + wb10(:,index)
             wb01(:,index)           =  wb01(:,master_injector_no) + wb01(:,index)
             wb00(:,index)           =  wb00(:,master_injector_no) + wb00(:,index)

             ftrapfi(:,:,index)      =  ftrapfi(:,:,master_injector_no)  + ftrapfi(:,:,index)
             hibrz(:,:,index)        =  hibrz(:,:,master_injector_no)    + hibrz(:,:,index)
             hdepz(:,:,index)        =  hdepz(:,:,master_injector_no)    + hdepz(:,:,index)
             zetaz(:,:,index)        =  zetaz(:,:,master_injector_no)    + zetaz(:,:,index)
             angmpz(:,:,index)       =  angmpz(:,:,master_injector_no)   + angmpz(:,:,index)*1.E-7_dp  ! CONVERT TO KG M^2/SEC 
             ir = SIZE(hicmz,1) ; ic = SIZE(hicmz,2) ; i3 = SIZE(hicmz,4) 
             DO j=1,i3              
                DO k=1,ic
                   DO l =1,ir
                        hicmz(l,k,index,j) = hicmz(l,k,master_injector_no,j) + hicmz(l,k,index,j)
                   ENDDO
                ENDDO
            ENDDO
 

        RETURN

      END SUBROUTINE nfreya_update_master




      SUBROUTINE load_plasma_prop
!-----------------------------------------------------------------------------
! -- phys_injctr is real physical injector number
! -- the quantities indexed by k below are all stored on the master at this point
! -- (by sub nfreya_store_packed. The index of, for example fap(:,j) 
! -- is for pseudo injector number j. Here we collect this info into 
! -- physical injector number p_injct
! -- **_stdev arrays were allocated and zeroed (even for single processor case)
! -- But in single processor case we do not calculate these arrays. Instead
! -- they are left at the intialized values and printed to the statefile
! -- as such.
!-----------------------------------------------------------------------------

        USE nrtype,                                   ONLY : DP,I4B

        USE nf_param,                                 ONLY : ke,kcm

        USE MPI_data,                                 ONLY : numprocs

        USE neutral_beams,                            ONLY : nmbrz,                 &
                                                             fap,fwall,forb,bneut,  &
                                                             bion,fb00,fb01,fb10,   &
                                                             fb11,wb00,wb01,wb10,   &
                                                             wb11,pbeam,ebeam,vbeam,&
                                                             fber,ftrapfi,hibrz,    &
                                                             hdepz,zetaz,angmpz,    &
                                                             hicmz,calc_variance,   &
                                                             nsample_izpt

        USE Plasma_properties,                        ONLY : neut_beam

        USE zonal_data,                               ONLY :  mf,mfm1

        IMPLICIT NONE

        INTEGER (I4B) dup_ct,i,j,ir,ic,k,p_injct



 nprocs:       IF(numprocs .GT.1)THEN 


    phys_injct:     DO p_injct = 1, no_physical_injectors
            dup_ct = izero
            Do k =1, no_injectors

 
               IF(phys_injctr_wk_list(k) == p_injct)THEN ! pick off pseudo injectors that match physical injector p_injct

                  dup_ct = dup_ct + 1

                  neut_beam%pbeam(:,p_injct) = pbeam(:,k) ! power to apperatures
 
                  neut_beam%ebeam(:,p_injct) = ebeam(:,k)

  
                  neut_beam%nmbrz(:,p_injct)        = neut_beam%nmbrz(:,p_injct)     + nmbrz(:,k)

                  neut_beam%nsample_izpt(:,p_injct) = neut_beam%nsample_izpt(:,p_injct) + nsample_izpt(:,k) 


                  neut_beam%fap  (:,p_injct)        = neut_beam%fap(:,p_injct)       + fap  (:,k)


                  neut_beam%fwall(:,p_injct)        = neut_beam%fwall(:,p_injct)     + fwall(:,k)


                  neut_beam%forb (:,p_injct)        = neut_beam%forb(:,p_injct)      + forb (:,k) 

                  neut_beam%bneut(:,p_injct)         = neut_beam%bneut(:,p_injct)     + bneut(:,k)
                  neut_beam%bion (:,p_injct)         = neut_beam%bion(:,p_injct)      + bion (:,k)

                  neut_beam%vbeam(:,p_injct)         = neut_beam%vbeam(:,p_injct)     + vbeam (:,k)


                  neut_beam%fb00 (:,p_injct)         = neut_beam%fb00(:,p_injct)      + fb00 (:,k)
                  neut_beam%fb01 (:,p_injct)         = neut_beam%fb01(:,p_injct)      + fb01 (:,k)
                  neut_beam%fb10 (:,p_injct)         = neut_beam%fb10(:,p_injct)      + fb10 (:,k) 
                  neut_beam%fb11 (:,p_injct)         = neut_beam%fb11(:,p_injct)      + fb11 (:,k)
                  neut_beam%fber (:,p_injct)         = neut_beam%fber(:,p_injct)      + fber (:,k)

                  ! set orbit widths to meters for output to netcdf file:
                  neut_beam%wb00 (:,p_injct)         = neut_beam%wb00(:,p_injct)      + cm2m*wb00 (:,k)
                  neut_beam%wb01 (:,p_injct)         = neut_beam%wb01(:,p_injct)      + cm2m*wb01 (:,k)
                  neut_beam%wb10 (:,p_injct)         = neut_beam%wb10(:,p_injct)      + cm2m*wb10 (:,k)
                  neut_beam%wb11 (:,p_injct)         = neut_beam%wb11(:,p_injct)      + cm2m*wb11 (:,k)

                  neut_beam%ftrapfi(:,:,p_injct)     = neut_beam%ftrapfi(:,:,p_injct) + ftrapfi (:,:,k)
                  neut_beam%hibrz (:,:,p_injct)      = neut_beam%hibrz(:,:,p_injct)   + hibrz (:,:,k)

                  neut_beam%hdepz (:,:,p_injct)      = neut_beam%hdepz(:,:,p_injct)   + hdepz (:,:,k)
                  neut_beam%zetaz (:,:,p_injct)      = neut_beam%zetaz(:,:,p_injct)   + zetaz (:,:,k)
                  neut_beam%angmpz(:,:,p_injct)      = neut_beam%angmpz(:,:,p_injct)  + angmpz (:,:,k)*1.E-7_dp  ! CONVERT TO KG M^2/SEC 


                  neut_beam%hicmz (:,:,p_injct,:)    = neut_beam%hicmz(:,:,p_injct,:) + hicmz (:,:,k,:)


                  IF(calc_variance)THEN
                     ir = SIZE(nmbrz,1) 
                     CALL INT_sumsq(nmbrz(:,k),ir,neut_beam%stdev_nmbrz(:,p_injct))

                     ir = SIZE(fap,1) 
                     CALL sumsq(bion(:,k),ir,neut_beam%stdev_bion(:,p_injct))
                     CALL sumsq(bneut(:,k),ir,neut_beam%stdev_bneut(:,p_injct))
                     CALL sumsq(fap(:,k),ir,neut_beam%stdev_fap(:,p_injct))
                     CALL sumsq(fwall(:,k),ir,neut_beam%stdev_fwall(:,p_injct))
                     CALL sumsq(forb(:,k),ir,neut_beam%stdev_forb(:,p_injct))
                     CALL sumsq(fber(:,k),ir,neut_beam%stdev_fber(:,p_injct))
                     CALL sumsq(fb11(:,k),ir,neut_beam%stdev_fb11(:,p_injct))
                     CALL sumsq(fb00(:,k),ir,neut_beam%stdev_fb00(:,p_injct))
                     CALL sumsq(fb01(:,k),ir,neut_beam%stdev_fb01(:,p_injct))
                     CALL sumsq(fb10(:,k),ir,neut_beam%stdev_fb10(:,p_injct))
                     CALL sumsq(wb11(:,k),ir,neut_beam%stdev_wb11(:,p_injct))
                     CALL sumsq(wb00(:,k),ir,neut_beam%stdev_wb00(:,p_injct))
                     CALL sumsq(wb01(:,k),ir,neut_beam%stdev_wb01(:,p_injct))
                     CALL sumsq(wb10(:,k),ir,neut_beam%stdev_wb10(:,p_injct))

                     ir = mfm1 
                     DO j = 1,ke
                        CALL sumsq(ftrapfi(:,j,k),ir,neut_beam%stdev_ftrapfi(:,j,p_injct))
                        CALL sumsq(hibrz(:,j,k),ir,neut_beam%stdev_hibrz(:,j,p_injct))
                        CALL sumsq(hdepz(:,j,k),ir,neut_beam%stdev_hdepz(:,j,p_injct))
                        CALL sumsq(zetaz(:,j,k),ir,neut_beam%stdev_zetaz(:,j,p_injct))
                        CALL sumsq(angmpz(:,j,k),ir,neut_beam%stdev_angmpz(:,j,p_injct))
                     ENDDO

                     ir = mfm1 ;
                     DO j = 1,kcm
                        DO i = 1,ke
                           CALL sumsq(hicmz(:,i,k,j),ir,neut_beam%stdev_hicmz(:,i,p_injct,j))
                        ENDDO
                     ENDDO
                  ENDIF

               ENDIF
            ENDDO
            IF(dup_ct .LT. 1_I4B)THEN
               lerrno = iomaxerr + 233_I4B
               CALL terminate(lerrno,nlog)
            ELSEIF(dup_ct .GT. 1_I4B)THEN

               neut_beam%nmbrz(:,p_injct)          = neut_beam%nmbrz(:,p_injct)/dup_ct   
               neut_beam%fap  (:,p_injct)          = neut_beam%fap(:,p_injct)/dup_ct      
               neut_beam%fwall(:,p_injct)          = neut_beam%fwall(:,p_injct)/dup_ct    
               neut_beam%forb (:,p_injct)          = neut_beam%forb(:,p_injct) /dup_ct    

               neut_beam%bneut(:,p_injct)          = neut_beam%bneut(:,p_injct)/dup_ct   
               neut_beam%bion (:,p_injct)          = neut_beam%bion(:,p_injct)/dup_ct    

               neut_beam%vbeam(:,p_injct)          = neut_beam%vbeam(:,p_injct)/dup_ct   


               neut_beam%fb00 (:,p_injct)          = neut_beam%fb00(:,p_injct)/dup_ct    
               neut_beam%fb01 (:,p_injct)          = neut_beam%fb01(:,p_injct)/dup_ct    
               neut_beam%fb10 (:,p_injct)          = neut_beam%fb10(:,p_injct)/dup_ct    
               neut_beam%fb11 (:,p_injct)          = neut_beam%fb11(:,p_injct)/dup_ct    
               neut_beam%fber (:,p_injct)          = neut_beam%fber(:,p_injct)/dup_ct    


               neut_beam%wb00 (:,p_injct)          = neut_beam%wb00(:,p_injct)/dup_ct    
               neut_beam%wb01 (:,p_injct)          = neut_beam%wb01(:,p_injct)/dup_ct    
               neut_beam%wb10 (:,p_injct)          = neut_beam%wb10(:,p_injct)/dup_ct    
               neut_beam%wb11 (:,p_injct)          = neut_beam%wb11(:,p_injct)/dup_ct    

               neut_beam%ftrapfi(:,:,p_injct)      = neut_beam%ftrapfi(:,:,p_injct)/dup_ct 
               neut_beam%hibrz  (:,:,p_injct)      = neut_beam%hibrz(:,:,p_injct)/dup_ct 
               neut_beam%hdepz  (:,:,p_injct)      = neut_beam%hdepz(:,:,p_injct)/dup_ct 
               neut_beam%zetaz  (:,:,p_injct)      = neut_beam%zetaz(:,:,p_injct)/dup_ct 
               neut_beam%angmpz (:,:,p_injct)      = neut_beam%angmpz(:,:,p_injct)/dup_ct 
               neut_beam%hicmz  (:,:,p_injct,:)    = neut_beam%hicmz(:,:,p_injct,:)/dup_ct

               IF(calc_variance)THEN
                  neut_beam%stdev_nmbrz(:,p_injct) = SQRT( neut_beam%stdev_nmbrz(:,p_injct)/(dup_ct-1._DP) &
                       -(dup_ct/(dup_ct-1._DP))* neut_beam%nmbrz(:,p_injct)**2)

                  neut_beam%stdev_bion(:,p_injct)   = SQRT( neut_beam%stdev_bion(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))* neut_beam%bion (:,p_injct)**2)

                  neut_beam%stdev_bneut(:,p_injct)   = SQRT( neut_beam%stdev_bneut(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%bneut (:,p_injct)**2)

                  neut_beam%stdev_fap(:,p_injct)   = SQRT( neut_beam%stdev_fap(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%fap(:,p_injct)**2)

                  neut_beam%stdev_fwall(:,p_injct)   = SQRT( neut_beam%stdev_fwall(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%fwall(:,p_injct)**2)

                  neut_beam%stdev_forb(:,p_injct)   = SQRT( neut_beam%stdev_forb(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%forb(:,p_injct)**2)

                  neut_beam%stdev_fber(:,p_injct)   = SQRT( neut_beam%stdev_fber(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%fber(:,p_injct)**2)

                  neut_beam%stdev_fb11(:,p_injct)   = SQRT( neut_beam%stdev_fb11(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%fb11(:,p_injct)**2)

                  neut_beam%stdev_fb00(:,p_injct)   = SQRT( neut_beam%stdev_fb00(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%fb00(:,p_injct)**2)

                  neut_beam%stdev_fb10(:,p_injct)   = SQRT( neut_beam%stdev_fb10(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%fb10(:,p_injct)**2)

                  neut_beam%stdev_fb01(:,p_injct)   = SQRT( neut_beam%stdev_fb01(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%fb01(:,p_injct)**2)


                  neut_beam%stdev_wb01(:,p_injct)   = SQRT( neut_beam%stdev_wb01(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%wb01(:,p_injct)**2)

                  neut_beam%stdev_wb00(:,p_injct)   = SQRT( neut_beam%stdev_wb00(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%wb00(:,p_injct)**2)

                  neut_beam%stdev_wb10(:,p_injct)   = SQRT( neut_beam%stdev_wb10(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%wb10(:,p_injct)**2)

                  neut_beam%stdev_wb11(:,p_injct)   = SQRT( neut_beam%stdev_wb11(:,p_injct)/(dup_ct-1._DP) &
                       - (dup_ct/(dup_ct-1._DP))*neut_beam%wb11(:,p_injct)**2)


                  ir = mfm1 
                  DO j = 1,ke
                     neut_beam%stdev_ftrapfi(:,j,p_injct) = SQRT( neut_beam%stdev_ftrapfi(:,j,p_injct)/(dup_ct-1._DP) &
                          - (dup_ct/(dup_ct-1._DP))*neut_beam%ftrapfi(:,j,p_injct)**2)

                     neut_beam%stdev_hibrz(:,j,p_injct) = SQRT( neut_beam%stdev_hibrz(:,j,p_injct)/(dup_ct-1._DP) &
                          - (dup_ct/(dup_ct-1._DP))*neut_beam%hibrz(:,j,p_injct)**2)

                     neut_beam%stdev_hdepz(:,j,p_injct) = SQRT( neut_beam%stdev_hdepz(:,j,p_injct)/(dup_ct-1._DP) &
                          - (dup_ct/(dup_ct-1._DP))*neut_beam%hdepz(:,j,p_injct)**2)

                     neut_beam%stdev_zetaz(:,j,p_injct) = SQRT( neut_beam%stdev_zetaz(:,j,p_injct)/(dup_ct-1._DP) &
                          - (dup_ct/(dup_ct-1._DP))*neut_beam%zetaz(:,j,p_injct)**2)

                     neut_beam%stdev_angmpz(:,j,p_injct) = SQRT( neut_beam%stdev_angmpz(:,j,p_injct)/(dup_ct-1._DP) &
                          - (dup_ct/(dup_ct-1._DP))*neut_beam%angmpz(:,j,p_injct)**2)

                  ENDDO

                  ir = mfm1 ;
                  DO j = 1,kcm
                     DO i = 1,ke
                        neut_beam%stdev_hicmz(:,i,p_injct,j) = SQRT( neut_beam%stdev_hicmz(:,i,p_injct,j)/(dup_ct-1._DP) &
                             - (dup_ct/(dup_ct-1._DP))*neut_beam%hicmz(:,i,p_injct,j)**2)
                     ENDDO
                  ENDDO
               ENDIF

            ENDIF

         ENDDO phys_injct
      ELSE  nprocs
            ! single processor case load  output arrays:
 
                  neut_beam%pbeam(:,:)         = pbeam(:,:) ! power to apperatures
                  neut_beam%ebeam(:,:)         = ebeam(:,:)
                  neut_beam%nmbrz(:,:)         = nmbrz(:,:)
                  neut_beam%nsample_izpt(:,:)  = nsample_izpt(:,:) 
                  neut_beam%fap  (:,:)         = fap(:,:)
                  neut_beam%fwall(:,:)         = fwall(:,:)
                  neut_beam%forb (:,:)         = forb(:,:)
                  neut_beam%bneut(:,:)         = bneut(:,:)
                  neut_beam%bion (:,:)         = bion(:,:)
                  neut_beam%vbeam(:,:)         = vbeam(:,:)
                  neut_beam%fb00 (:,:)         = fb00(:,:)
                  neut_beam%fb01 (:,:)         = fb01(:,:)
                  neut_beam%fb10 (:,:)         = fb10(:,:) 
                  neut_beam%fb11 (:,:)         = fb11(:,:)
                  neut_beam%fber (:,:)         = fber(:,:)
                  neut_beam%wb00 (:,:)         = cm2m*wb00(:,:)
                  neut_beam%wb01 (:,:)         = cm2m*wb01(:,:)
                  neut_beam%wb10 (:,:)         = cm2m*wb10(:,:)
                  neut_beam%wb11 (:,:)         = cm2m*wb11(:,:)
                  neut_beam%ftrapfi(:,:,:)     = ftrapfi(:,:,:)
                  neut_beam%hibrz (:,:,:)      = hibrz(:,:,:)
                  neut_beam%hdepz (:,:,:)      = hdepz(:,:,:)
                  neut_beam%zetaz (:,:,:)      = zetaz(:,:,:)
                  neut_beam%angmpz(:,:,:)      = angmpz(:,:,:)  *1.E-7_dp  ! CONVERT TO KG M^2/SEC 
                  neut_beam%hicmz (:,:,:,:)    = hicmz(:,:,:,:)

         
      ENDIF nprocs



      RETURN

    END SUBROUTINE load_plasma_prop



    SUBROUTINE int_sumsq(int_vec_in,ir,real_sq_out)
!--------------------------------------------------------------
! --
!--------------------------------------------------------------
    USE nrtype,                                     ONLY :DP, I4B
    IMPLICIT NONE
    INTEGER(I4B)ir,j
    INTEGER(I4B) int_vec_in(ir)
    REAL(DP) real_sq_out(ir)
       DO j=1,ir
             real_sq_out(j) = FLOAT(int_vec_in(j)**2) + real_sq_out(j)
       ENDDO
       RETURN
    END SUBROUTINE int_sumsq



    SUBROUTINE sumsq(vec_in,ir,sq_out)
!--------------------------------------------------------------
! --
!--------------------------------------------------------------
    USE nrtype,                                     ONLY :DP, I4B
    IMPLICIT NONE
    INTEGER(I4B)ir,i
    REAL(DP)  vec_in(ir),sq_out(ir)
          DO i =1,ir
             sq_out(i) = vec_in(i)**2 + sq_out(i)
          ENDDO
       RETURN
    END SUBROUTINE sumsq
#ifdef USEMPI
      SUBROUTINE nfreya_master_slave
!-----------------------------------------------------------------
! -- Use this routine to run nfreya in multiple  processor mode
! -- Note that compilation/linking is always done with mpi libs assuming multiple 
! -- processors. The parallelism is done over injectors (not over
! -- number of pseudo particles injected). At this time we use only
! -- a single communicator group. Passing multiple cpus to each 
! -- injector may be implemented at a later time.
!-------------------------------------------------------------HSJ----
        USE nf_param,                  ONLY : kb
        USE MPI_data,                  ONLY : mpi_start_time,mpi_end_time,         & 
                                              numprocs,master,myid,mpiierr,        &
                                              mpi_status,mpi_buffer,buff_size
         IMPLICIT NONE
         INTEGER(I4B)injctr,j,k,slave_no,my_mpi_tag,injctr_done,injrcv,                &
                     nfreya_max_pack,nfreya_receive_packed,nfreya_send_packed,     &
                     count,slave_id,last_slave,excess_proc,nstore,                 &
                     duplicate_injector_no,orig_injector_number,injctr_none,       &
                     n_m_recv,p_injct

         LOGICAL master_done,all_sent,no_more_work

         REAL(DP) t_start,t_end,t_send,t_recv,t_master_com, &
                  t_slave_compt,t_slave_end,t_slave_start,  &
                  t_store_packed_start,t_store_packed_end
         REAL(DP), ALLOCATABLE,DIMENSION(:) :: t_m_recv,t_store_packed 
 
         master_injector_no = -1
         nfreya_send_packed = 1 ; nfreya_receive_packed = 1  ; nfreya_max_pack = 1
         n_m_recv = izero

         IF(myid == master)THEN        
                !----------------------------------------------------------------
                ! master passes out injectors to as many processors as possible:
                !
                ! if no processors <  no injectors then at the end of this loop
                ! we still have some injectors to hand out. These will be
                ! done below by master and freed up slaves
                !
                ! if no processors = no injectors then at the end of this loop
                ! all slaves are  assigned out and the master does the last injector
                !
                ! if no processors > no injectors then the master process will send
                ! injector number  =  no_injectors +1  to  the remaining slaves 
                ! This is used as a flag to tell the extraneous processes that they 
                ! to sit at the MPI_Barrier until the other processes finish.
                ! This is wasteful and one could split injectors into multiple
                ! injectors, each with an appropriate fraction of the pseudo
                ! particles. But this is not done  at present.
                !----------------------------------------------------------------
                IF(nfreya_vb)PRINT *,'No physical injectors =',no_physical_injectors
                IF(nfreya_vb)PRINT *,'No computational injectors =',no_injectors
                IF(nfreya_vb)PRINT *,'No of cpus  available =',numprocs
                ALLOCATE(injctr_wk_list(no_injectors))
                ALLOCATE(phys_injctr_wk_list(no_injectors))
                ALLOCATE(t_m_recv(no_injectors)) ; t_m_recv(:) = zeroc
                ALLOCATE(t_store_packed(no_injectors)) ; t_store_packed(:) = zeroc
                t_start = MPI_WTIME() ; t_master_com = zeroc ; t_store_packed = zeroc
                injctr = izero ; injctr_done = izero 
                master_done = .FALSE. ; all_sent = .FALSE. ; nstore = izero
                DO slave_no = 0,MIN(no_injectors,numprocs-1) ! note numprocs can be less than no_injectors
                   IF(slave_no .NE. master)THEN
                      injctr = injctr + 1          ! get an injctr to send to slave
                      my_mpi_tag = injctr
                      t_send = MPI_WTIME()
                      CALL MPI_SEND(injctr,1,MPI_INTEGER,slave_no,my_mpi_tag,&
                           MPI_COMM_WORLD,mpiierr)
                      t_recv = MPI_WTIME()
                      t_master_com = t_master_com + t_recv - t_send
                      last_slave = slave_no
                      IF(nfreya_vb)PRINT  *,'mstra sent injctr to slave ',injctr,slave_no
                      injctr_wk_list(injctr) = slave_no
                   ENDIF
                 ENDDO
                 IF(injctr == no_injectors )all_sent = .TRUE.
                 IF(all_sent .AND. numprocs .GT. no_injectors)THEN
                    ! arrive here if no processors > no injectors
                    ! this also accounts for the situation where the number
                    ! of injetors was increased over the input number so that
                    ! a maximum of kb processors is involved
                    ! this is wasteful of cpu resources but is allowed here.
                    ! send message to excess slaves so that they
                    ! can enter the receive loop  below and act accordingly
                    injctr_none = no_injectors + 1
                    my_mpi_tag = injctr_none
                    DO slave_no = last_slave +1, numprocs -1
                       CALL MPI_SEND(MPI_BOTTOM,izero,MPI_INTEGER,slave_no,my_mpi_tag, &
                                                                MPI_COMM_WORLD,mpiierr ) 
                    ENDDO
                 ENDIF
              
                 DO WHILE (injctr_done .LT. no_injectors)
                    IF( .NOT. all_sent .AND. .NOT. master_done)THEN
                       injctr = injctr + 1   ! this injctr will be done by master
                       master_injector_no = injctr
                       IF(nfreya_vb)PRINT *,'master doing injector ',injctr
                       injctr_wk_list(injctr) = myid
                       CALL launch_nfreya(injctr)
                       master_done  = .TRUE. ! master did one injector
                       inj_list = 'single_injector'
                       injctr_done  = injctr_done +1 ! at this point master does not know
                       ! if slaves have completed so this is
                       ! correct
                    ENDIF

                    !----------------------------------------------------------------
                    ! now master receives results from any slaves
                    !----------------------------------------------------------------
                    t_send = MPI_WTIME()
                    CALL MPI_RECV(mpi_buffer,buff_size,MPI_PACKED,      &
                                  MPI_ANY_SOURCE,MPI_ANY_TAG,             &
                                  MPI_COMM_WORLD,mpi_status,mpiierr)
                    t_recv = MPI_WTIME()
                    n_m_recv = n_m_recv + 1_I4B
                    t_m_recv(n_m_recv) = t_recv - t_send
                    t_master_com = t_master_com + t_recv - t_send
                    my_mpi_tag     = mpi_status(MPI_TAG)
                    !print *,'mstrc  recevd reply from slave ',mpi_status(MPI_SOURCE)
                    IF(my_mpi_tag > 0)THEN 
                       injctr_done = injctr_done + 1     
                       slave_id  = mpi_status(MPI_SOURCE)
                       t_store_packed_start = MPI_WTIME()
                       CALL nfreya_store_packed(slave_id)  ! put data received from slaves
                       nstore = nstore +1
                       t_store_packed_end = MPI_WTIME()
                       t_store_packed(nstore) =   t_store_packed_end - t_store_packed_start 
                    ENDIF                                
                    slave_id  = mpi_status(MPI_SOURCE)
                    !send next injctr to idle processor:
                    IF(injctr .LT. no_injectors)THEN
                       injctr = injctr + 1 
                       my_mpi_tag = injctr
                       injctr_wk_list(injctr) = slave_id
                       t_send = MPI_WTIME()
                       CALL MPI_SEND(injctr,1,MPI_INTEGER,slave_id,my_mpi_tag,        &
                                                                MPI_COMM_WORLD,mpiierr)
                       t_recv = MPI_WTIME()
                       t_master_com = t_master_com + t_recv - t_send
                    ELSE   
                      ! send no work to be done message by setting tag 
                      ! to no_injectors +1 with zero length message 
                      ! (tags must be  >=0, <= MPI_TAG_UB)                   
                       my_mpi_tag = no_injectors +1        
                       t_send = MPI_WTIME()
                       CALL MPI_SEND(MPI_BOTTOM,izero,MPI_INTEGER,slave_id,my_mpi_tag, &
                                                                MPI_COMM_WORLD,mpiierr ) 
                        t_recv = MPI_WTIME()
                       t_master_com = t_master_com + t_recv - t_send
                    ENDIF
                 ENDDO

         ELSE IF(myid .NE. master)THEN    
                 ! only   slaves  execute this section
                 no_more_work = .FALSE. 
                 t_slave_compt = zeroc ; injrcv = izero
                 DO WHILE(no_more_work .eqv. .FALSE.)
                    !slave waits for injector number  message from master
                    CALL MPI_RECV(injctr,1,MPI_INTEGER,master,MPI_ANY_TAG,             &
                                                      MPI_COMM_WORLD,mpi_status,mpiierr)

                    IF(MPI_STATUS(MPI_TAG) < no_injectors + 1)THEN 
                       ! received valid injector request  from master
                       injrcv = injrcv + 1
                       t_slave_start = MPI_WTIME()
                       CALL launch_nfreya(injctr)
                       t_slave_end   = MPI_WTIME()
                       t_slave_compt = t_slave_compt + t_slave_end - t_slave_start

                       !pack the data and send results of calculations to master:
                       CALL  nfreya_pack_data(injctr,t_slave_compt,mpi_tag) 
                    ELSE
                       ! master sent invalid injector number because there are no
                       ! more injectors to be processed:
                       no_more_work = .TRUE.
                    ENDIF
                 ENDDO
         ENDIF



         IF(myid == master)THEN
            t_end = MPI_WTIME()
            IF(nfreya_vb)PRINT *,"Elapsed time in master/slave module =",t_end - t_start
            IF(nfreya_vb)PRINT *,'master communication time =',t_master_com
            !IF(nfreya_vb)PRINT *,'t_store_packed =',t_store_packed
            IF(nfreya_vb)PRINT *,'t_m_recv =',t_m_recv
            IF(nfreya_vb)PRINT *,'avg receive time =',SUM(t_m_recv)/n_m_recv
            IF(nfreya_vb)PRINT *,'injctr_wk_list =', injctr_wk_list
               DO j=1,SIZE(injctr_wk_list)
                   IF(j .GT. no_physical_injectors)THEN
                     p_injct = pseudo_injectors(j - no_physical_injectors)
                  ELSE 
                     p_injct = j
                  ENDIF
                  phys_injctr_wk_list(j) = p_injct
               ENDDO
               IF(ALLOCATED(pseudo_injectors)) PRINT *,'pseudo_injectors =',pseudo_injectors
               PRINT *,'phys_injctr_wk_list =',phys_injctr_wk_list
         ELSE
            IF(nfreya_vb)PRINT *,'slave compute  time =',t_slave_compt
         ENDIF




           CALL MPI_BARRIER(MPI_COMM_WORLD,mpiierr) ! wait until all processors arrive here


         !---------------------------------------------------------------------------------
         ! NOTE:
         ! data collected on master is NOT rebroadcast to slaves.
         ! Hence only master process should be used for further processing of data beyond this point
         !---------------------------------------------------------------------------------

      END SUBROUTINE nfreya_master_slave




      SUBROUTINE nfreya_pack_data(injctr,t_slave_compt,mpi_tag) 
!---------------------------------------------------------------------------------
! -- pack the data into a single list for shipping to master process
! -- NOTE: ROUTINES THAT MUST BE KEPT IN SYNC:
! --       nfreya_pack_data
! --       nfreya_update_master
! --       nfreya_store_packed 
! --       allocate_mpi_buffer
!
!  the calculated  quantities are:
!
!     bion(ie,ib)     intensity of ion beam (particles/s)
!     bneut(ie,ib)    intensity of neutral beam (particles/s)
!     pbeam(ie,ib)     beam power to aperture (w)
!     fap(ie,ib)      fraction of beam stopped by aperture
!     fwall(ie,ib)    fraction of beam incident on wall (shinethrough)
!     forb(ie,ib)     fraction of beam lost on orbits
!     ftrapfi(i,ie,ib)fraction of trapped fast ions in each zone
!     fb11(ie,ib)     fraction of ions passing and axis-encircling
!     fb10(ie,ib)     fraction of ions passing and not encircling
!     fb01(ie,ib)     fraction of ions trapped and axis-encircling
!     fb00(ie,ib)     fraction of ions trapped and not encircling
!     fber(ie,ib)     fraction of ions trapped for which error was detected
!     hibrz(i,ie,ib)  normalized hot ion birth rate
!     hdepz(i,ie,ib)  normalized hot ion deposition rate
!     hicmz(i,ie,ib,imd)
!                     hot ion creation mode (i.e. ionization or cx)
!     n_izpt            number of birth points to be plotted
!     x_izpt(ii,ie,ib)        x coordinate of birth point
!     y_izpt(ii,ie,ib)        y coordinate of birth point
!     z_izpt(ii,ie,ib)        z coordinate of birth point
!     r_izpt(ii,ie,ib) 
!     vx_izpt(ii,ie,ib)          x component of birth velocity
!     vy_izpt(ii,ie,ib)          y component of birth velocity
!     vz_izpt(ii,ie,ib)          z component of birth velocity
!     pitch_a(ii,ie,ib)          pitch angle of ion         
!     nmbrz(1:mfmi,ib)   counts ions born in zone i
!     wb11(ie,ib)     orbit width of fb11 ions (cm)
!     wb10(ie,ib)     orbit width of fb10 ions (cm)
!     wb01(ie,ib)     orbit width of fb01 ions (cm)
!     wb00(ie,ib)     orbit width of fb00 ions (cm)
!     zetaz(i,ie,ib)  average pitch angle cosine of deposited hot ions
!     angmpz(i,ie,ib) average toroidal angular momentum of a single ion
!                     born in zone i.  no orbit smearing is done.  for a
!                     pencil beam, angmpz will be constant across all
!                     zones.  the orbit smearing of the toroidal
!                     momentum is accomplished (in postnub) by
!                     multiplying angmpf by the smeared deposition rate
!                     of fast ions.
!---------------------------------------------------------------------------------

      USE neutral_beams,                               ONLY :                              &
                                                               vbeam,cangv ,cangh ,sangv,  &
                                                               sangh ,thetp , npskip,      &
                                                               thetpp ,costp ,sintp ,      &
                                                               costpp ,sintpp,nsourc,      &
                                                               iatype,kt,npulse,pbeamOn,   &
                                                               pbeamOff,beam_on,            &
                                                               beam_end,source2_phase,     &
                                                               nap,nashape,bvofset,        &
                                                               bhofset,bion,bneut,pbeam,   &
                                                               fwall,fap,forb,             &
                                                               fb11,fb01,fb10,fb00,fber,   &
                                                               wb11,wb00,wb01,wb10,bleni,  &
                                                               ftrapfi,hibrz,hdepz,angmpz, &
                                                               zetaz,hicmz,olossc,nbshape, &
                                                               nashape,bheigh,bwidth,bhfoc,&
                                                               bvfoc,bhdiv,bvdiv,sfrac1,   &
                                                               naptr,aheigh,awidth,alen,   &
                                                               blenp,rpivot,zpivot,n_izpt, &
                                                               vx_izpt,vy_izpt,vz_izpt,    &
                                                               x_izpt,y_izpt,z_izpt,       &
                                                               r_izpt,pitch_a,nmbrz,       &
                                                               nsample_izpt
                                                            

       USE MPI_data,                                    ONLY : myid,mpiierr,mpi_start_time,&
                                                               mpi_end_time,numprocs,master, &
                                                               mpi_buffer,buff_size

        INTEGER(I4B) j,k,l,injctr,ir,ic,i3,count,int_size,    &
                     position,double_size,testing,one,        &
                     mpi_tag,no_0di_l,no_3dr_p_l ,            &
                     no_2di_l,no_2dr_l,no_3dr_l,no_4dr_l,     &
                     no_2di_l2
        REAL(DP) t_slave_compt

        REAL(DP),ALLOCATABLE,DIMENSION(:,:,:)  :: work
       ! NOTE mpi_buffer size is set in sub allocate_mpi_buffer,file nfreya_init.f90

             position = izero ; one = 1_I4B
             !--------------------------------------------------------------------------
             ! scalars (also see end of routine for more)
             !--------------------------------------------------------------------------
             no_0di_l = izero
             CALL MPI_pack(injctr,one,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)
             no_0di_l = no_0di_l + one
 
             CALL MPI_pack(n_izpt,one,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)
             no_0di_l = no_0di_l + one


             !--------------------------------------------------------------------------
             ! 2d integer arrays (each group of arrays must have common size)
             !--------------------------------------------------------------------------
              no_2di_l = izero
              ir = SIZE(nmbrz,1) ! ir = # zones
              CALL MPI_pack(nmbrz(1,injctr),ir,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)
              no_2di_l = no_2di_l + 1

              no_2di_l2 = izero
              ir = SIZE(nsample_izpt,1) ! ir = ke
              CALL MPI_pack(nsample_izpt(1,injctr),ir,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)
              no_2di_l2 = no_2di_l2 + 1

             !--------------------------------------------------------------------------
             ! 3d real  arrays of size (no_birth_points,ke,no_injectors)
             ! for a given injector
             !--------------------------------------------------------------------------
 
             no_3dr_p_l = izero
             ir = size(vx_izpt,1) ; ic = SIZE(vx_izpt,2)
             count = ir*ic
             CALL MPI_pack(vx_izpt(:,:,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
              no_3dr_p_l = no_3dr_p_l + 1

             CALL MPI_pack(vy_izpt(:,:,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
              no_3dr_p_l = no_3dr_p_l + 1

             CALL MPI_pack(vz_izpt(:,:,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
              no_3dr_p_l = no_3dr_p_l + 1

             CALL MPI_pack(x_izpt(:,:,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
              no_3dr_p_l = no_3dr_p_l + 1

             CALL MPI_pack(y_izpt(:,:,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
              no_3dr_p_l = no_3dr_p_l + 1

             CALL MPI_pack(z_izpt(:,:,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
              no_3dr_p_l = no_3dr_p_l + 1

             CALL MPI_pack(r_izpt(:,:,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
              no_3dr_p_l = no_3dr_p_l + 1

             CALL MPI_pack(pitch_a(:,:,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
              no_3dr_p_l = no_3dr_p_l + 1

             !--------------------------------------------------------------------------
             ! only column injctr  of 2d arrays  is filled by this slave so 
             ! send only this column ==>fb11(1,injctr) for example
             !--------------------------------------------------------------------------
             ir = SIZE(fb11,1) ! all 2d  arrays to be packed have this size
             count = ir
             no_2dr_l = izero

             CALL MPI_pack(bion(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(bneut(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(pbeam(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1


             CALL MPI_pack(fap(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(fwall(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(forb(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1


             CALL MPI_pack(fb11(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(fb10(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(fb01(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(fb00(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(fber(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1



             CALL MPI_pack(wb11(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(wb10(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(wb01(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1
             CALL MPI_pack(wb00(1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_2dr_l = no_2dr_l +1



             !--------------------------------------------------------------------------
             ! only 2d array labelled by injector   of 3d arrays  is filled by this slave so 
             ! send only this face ==>hibrz(1:,1:,injctr) for example
             !--------------------------------------------------------------------------
             ir = SIZE(hibrz,1) 
             ic = SIZE(hibrz,2)
             count = ir*ic
             no_3dr_l = izero

             CALL MPI_pack(ftrapfi(1,1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_3dr_l = no_3dr_l +1

             CALL MPI_pack(hibrz(1,1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_3dr_l = no_3dr_l +1

             CALL MPI_pack(hdepz(1,1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_3dr_l = no_3dr_l +1

             CALL MPI_pack(zetaz(1,1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_3dr_l = no_3dr_l +1

             CALL MPI_pack(angmpz(1,1,injctr),count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
             no_3dr_l = no_3dr_l +1


             !--------------------------------------------------------------------------
             ! only 3d array labelled by injector  of 4d arrays  is filled by this slave so 
             ! send only this 3d array  ==> hicmz(:,:,injctr,:), imd =1,2,3
             !--------------------------------------------------------------------------
             ir = SIZE(hicmz,1) ; ic = SIZE(hicmz,2) ; i3 = SIZE(hicmz,4) 
             IF(ALLOCATED(work))DEALLOCATE(work)
             ALLOCATE(work(ir,ic,i3))
               DO j=1,i3               ! I dont like mpi stride routines
                  DO k=1,ic
                     DO l =1,ir
                        work(l,k,j) = hicmz(l,k,injctr,j)
                     ENDDO
                  ENDDO
               ENDDO
               count = ir*ic*i3
               no_4dr_l = izero
               CALL MPI_pack(work,count,MPI_DOUBLE_PRECISION,mpi_buffer,buff_size,position,&
                           MPI_COMM_WORLD,mpiierr)
               no_4dr_l = no_4dr_l + 1
               DEALLOCATE(work)

             !-----------------------------------------------------------------------------
             ! send some error checking info:
             ! {6 = no_0_di_l,no_1d_l,no_2dr_l,no_3dr_l,no_4dr_l}
             !-----------------------------------------------------------------------------
              no_0di_l = no_0di_l + 7  ! eg 7 entries below
   
              CALL MPI_pack(no_0di_l,1,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)

              CALL MPI_pack(no_2di_l,1,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)

             CALL MPI_pack(no_2di_l2,1,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)

              CALL MPI_pack(no_3dr_p_l,1,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)

              CALL MPI_pack(no_2dr_l,1,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)

              CALL MPI_pack(no_3dr_l,1,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)

              CALL MPI_pack(no_4dr_l,1,MPI_INTEGER,mpi_buffer,buff_size,position, &
                           MPI_COMM_WORLD,mpiierr)


             !--------------------------------------------------------------------------
             ! done packing data, send the buffer to maser process:
             !--------------------------------------------------------------------------
             CALL MPI_send(mpi_buffer,buff_size,MPI_PACKED,master,mpi_tag,MPI_COMM_WORLD,mpiierr)

        RETURN

      END SUBROUTINE nfreya_pack_data
#endif


     SUBROUTINE launch_nfreya(injctr)
!---------------------------------------------------------------------------------
! -- preprocessing for nfreya call
! -- Here we call freyas with one injector only. The injector to be
! -- processed is indified in freyas list injs_to_do.
! -- This routine is called from master_slave routine.
! -- Not used if only a single cpu is selected by the user. 
! -- "INPUT"
!      injctr_id      injector number, real or pseudo injector
!---------------------------------------------------------------------------------
     USE nf_param,                                           ONLY : kb

        IMPLICIT NONE

        INTEGER(I4B) injctr,j,k,real_injctr_no,           &
                     injs_to_do,injctr_id

        injs_to_do = 1
        injctr_id  = injctr
        CALL freyas(injs_to_do, injctr_id,codeid,                 &
                   psi_g_cm2, rmhdgrid,rmhdgrid(1), rplasmax,     &
                   zmhdgrid,ymagn1, zplasmin,zplasmax, iexcit,    &
                   time_now,beam_sim_time_start)
         
        

        RETURN

     END SUBROUTINE launch_nfreya

#ifdef USEMPI
     SUBROUTINE nfreya_store_packed(slave_id)
!---------------------------------------------------------------------------------
! -- master receives packed data from slaves and puts it into proper slots
! -- NOTE: ROUTINES THAT MUST BE KEPT IN SYNC:
! --       nfreya_pack_data
! --       nfreya_update_master
! --       nfreya_store_packed 
! --       allocate_mpi_buffer
!---------------------------------------------------------------------------------
      USE MPI_data,                                     ONLY : myid,mpiierr,mpi_start_time,&
                                                               mpi_end_time,numprocs,master, &
                                                               mpi_buffer,buff_size

      USE neutral_beams,                                ONLY :                              &
                                                               vbeam,cangv ,cangh ,sangv,  &
                                                               sangh ,thetp , npskip,      &
                                                               thetpp ,costp ,sintp ,      &
                                                               costpp ,sintpp,nsourc,      &
                                                               iatype,kt,npulse,pbeamOn,   &
                                                               pbeamOff,beam_on,           &
                                                               beam_end,source2_phase,     &
                                                               nap,nashape,bvofset,        &
                                                               bhofset,bion,bneut,pbeam,   &
                                                               fap,fwall,forb,             &
                                                               fb11,fb01,fb10,fb00,fber,   &
                                                               wb11,wb00,wb01,wb10,bleni,  &
                                                               ftrapfi,hibrz,hdepz,angmpz, &
                                                               zetaz,hicmz,olossc,nbshape, &
                                                               nashape,bheigh,bwidth,bhfoc,&
                                                               bvfoc,bhdiv,bvdiv,sfrac1,   &
                                                               naptr,aheigh,awidth,alen,   &
                                                               blenp,rpivot,zpivot,        &
                                                               pseudo_injectors,           &
                                                               no_physical_injectors,      &
                                                               n_izpt,                     &
                                                               vx_izpt,vy_izpt,vz_izpt,    &
                                                               x_izpt,y_izpt,z_izpt,       &
                                                               r_izpt,pitch_a,nmbrz,       &
                                                               no_0di,no_3dr_p,no_2di,     &
                                                               no_2dr,no_3dr,no_4dr,       &
                                                               nsample_izpt,no_2di2


        IMPLICIT NONE

        INTEGER(I4B) slave_id
        INTEGER(I4B) injctr,ir,ic,i3,j,k,l,count,int_size,index,        &
                     position,elmt_cnt,testing,one,tag,        &
                     n_izpt_l,no_0di_l,no_2di_l,no_3dr_p_l,      &
                     no_2dr_l,no_3dr_l,no_4dr_l,               &
                     no_0di_s,no_2di_s,no_3dr_p_s,no_2dr_s,      &
                     no_3dr_s, no_4dr_s,no_2di_l2,no_2di_s2
        ! no_**_s*  are the number of arrays (of the type 0d,1d,2d,etc,
        ! with i ==> integer,r ==> real8)  that are sent by the slaves
        ! no_**_l* are the corresponding number picked up here.
        ! they must match otherwise we have gotten out of sync.
        ! this is checked for below.
        ! Also the initial buffer allocation (sub allocate_mpi_buffer)
        ! sets up the corresponding no_** arrays (eg no_2di,etc) that
        ! nust match as well 



        REAL(DP),     ALLOCATABLE,DIMENSION(:)          :: work

        REAL(DP),     ALLOCATABLE,DIMENSION(:,:)        :: work2d

        REAL(DP),     ALLOCATABLE,DIMENSION(:,:,:)      :: work3d



        REAL(DP) t_strt,t_end

        !---------------------------------------------------------------------
        ! -- the unpacking sequence must be correlated with the packing
        ! -- sequence ( see nfreya_pack_data)
        !---------------------------------------------------------------------
        position = izero ; one = 1_I4B

        t_strt = MPI_WTIME()


        !--------------------------------------------------------------------------
        ! scalars (also see end of routine for more)
        !--------------------------------------------------------------------------
        no_0di_l = izero
        CALL MPI_unpack(mpi_buffer,buff_size,position,injctr,one,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr)
        no_0di_l = no_0di_l +1

        CALL MPI_unpack(mpi_buffer,buff_size,position,n_izpt_l,one,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr)            
        no_0di_l = no_0di_l +1

        index  = injctr

       !--------------------------------------------------------------------------
       ! 2d integer arrays
       !--------------------------------------------------------------------------
        no_2di_l = izero
        count = SIZE(nmbrz,1)


        CALL MPI_unpack(mpi_buffer,buff_size,position,nmbrz(1,injctr),count,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr) 
        no_2di_l = no_2di_l +1
 
        no_2di_l2 = izero
        count = SIZE(nsample_izpt,1)
        CALL MPI_unpack(mpi_buffer,buff_size,position,nsample_izpt(1,injctr),count,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr) 
        no_2di_l2 = no_2di_l2 +1



       !--------------------------------------------------------------------------
       ! 3d real  arrays of size (no_birth_points,ke,no_injectors)
       ! for a given injector
       !--------------------------------------------------------------------------
        no_3dr_p_l = izero
        ir = size(vx_izpt,1) ; ic = SIZE(vx_izpt,2)
        !IF(ALLOCATED(work2d))DEALLOCATE(work2d)
        !ALLOCATE(WORK2d(ir,ic))
        count = ir*ic
 
         
        CALL MPI_unpack(mpi_buffer,buff_size,position,vx_izpt(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! vx_izpt
        no_3dr_p_l = no_3dr_p_l + 1

        CALL MPI_unpack(mpi_buffer,buff_size,position,vy_izpt(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! vy_izpt
        no_3dr_p_l = no_3dr_p_l + 1

        CALL MPI_unpack(mpi_buffer,buff_size,position,vz_izpt(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! vz_izpt
        no_3dr_p_l = no_3dr_p_l + 1




        CALL MPI_unpack(mpi_buffer,buff_size,position,x_izpt(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! x_izpt
        no_3dr_p_l = no_3dr_p_l + 1
 
        CALL MPI_unpack(mpi_buffer,buff_size,position,y_izpt(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! y_izpt
        no_3dr_p_l = no_3dr_p_l + 1
 
        CALL MPI_unpack(mpi_buffer,buff_size,position,z_izpt(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! z_izpt
        no_3dr_p_l = no_3dr_p_l + 1

        CALL MPI_unpack(mpi_buffer,buff_size,position,r_izpt(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! r_izpt
        no_3dr_p_l = no_3dr_p_l + 1




        CALL MPI_unpack(mpi_buffer,buff_size,position,pitch_a(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! pitch_a
        no_3dr_p_l = no_3dr_p_l + 1
 


       !--------------------------------------------------------------------------
       ! only column injctr  of 2d arrays  is filled by  slaves  
       !--------------------------------------------------------------------------
        no_2dr_l = izero
        elmt_cnt = SIZE(fb11,1)
       ! IF(ALLOCATED(work))DEALLOCATE(work)
       ! ALLOCATE(work(elmt_cnt))


        CALL MPI_unpack(mpi_buffer,buff_size,position,bion(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! bion
        no_2dr_l = no_2dr_l +1


        CALL MPI_unpack(mpi_buffer,buff_size,position,bneut(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! bneut
        no_2dr_l = no_2dr_l +1
 

        CALL MPI_unpack(mpi_buffer,buff_size,position,pbeam(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! pbeam
        no_2dr_l = no_2dr_l +1




        CALL MPI_unpack(mpi_buffer,buff_size,position,fap(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! fap
        no_2dr_l = no_2dr_l +1
 
 
        CALL MPI_unpack(mpi_buffer,buff_size,position,fwall(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! fwall
        no_2dr_l = no_2dr_l +1

        CALL MPI_unpack(mpi_buffer,buff_size,position,forb(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! forb
        no_2dr_l = no_2dr_l +1
 



        CALL MPI_unpack(mpi_buffer,buff_size,position,fb11(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! fb11
        no_2dr_l = no_2dr_l +1

        CALL MPI_unpack(mpi_buffer,buff_size,position,fb10(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! fb10
        no_2dr_l = no_2dr_l +1
 
        CALL MPI_unpack(mpi_buffer,buff_size,position,fb01(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! fb01
        no_2dr_l = no_2dr_l +1
        
        CALL MPI_unpack(mpi_buffer,buff_size,position,fb00(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! fb00
        no_2dr_l = no_2dr_l +1

        CALL MPI_unpack(mpi_buffer,buff_size,position,fber(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! fber
        no_2dr_l = no_2dr_l +1





        CALL MPI_unpack(mpi_buffer,buff_size,position,wb11(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! wb11
        no_2dr_l = no_2dr_l +1
 
        CALL MPI_unpack(mpi_buffer,buff_size,position,wb10(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! wb10
        no_2dr_l = no_2dr_l +1
 
        CALL MPI_unpack(mpi_buffer,buff_size,position,wb01(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! wb01
        no_2dr_l = no_2dr_l +1

        CALL MPI_unpack(mpi_buffer,buff_size,position,wb00(1,index),elmt_cnt,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! wb00
        no_2dr_l = no_2dr_l +1
  


       !--------------------------------------------------------------------------
       ! only face injector  of 3d arrays  is filled by this slave so 
       ! send only this face ==>hibrz(1:,1:,injctr) for example
       !--------------------------------------------------------------------------
        ir = SIZE(hibrz,1) 
        ic = SIZE(hibrz,2)
        count = ir*ic
        no_3dr_l = izero
        IF(ALLOCATED(work2d))DEALLOCATE(work2d)
        ALLOCATE(work2d(ir,ic))
        CALL MPI_unpack(mpi_buffer,buff_size,position,ftrapfi(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! ftrapfi
        no_3dr_l = no_3dr_l +1


        CALL MPI_unpack(mpi_buffer,buff_size,position,hibrz(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! hibrz
        no_3dr_l = no_3dr_l +1

        CALL MPI_unpack(mpi_buffer,buff_size,position,hdepz(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! hdepz
        no_3dr_l = no_3dr_l +1


        CALL MPI_unpack(mpi_buffer,buff_size,position,zetaz(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! zetaz
        no_3dr_l = no_3dr_l +1

        CALL MPI_unpack(mpi_buffer,buff_size,position,angmpz(1,1,index),count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! angmpz
        no_3dr_l = no_3dr_l +1



       !--------------------------------------------------------------------------
       ! only 3d array labelled by injector  of 4d arrays  is filled by this slave so 
       ! send only this 3d array  ==> hicmz(:,:,injctr,:), imd =1,2,3
       !--------------------------------------------------------------------------
       ir = SIZE(hicmz,1) ; ic = SIZE(hicmz,2) ; i3 = SIZE(hicmz,4) 
       count = ir*ic*i3
       IF(ALLOCATED(work3d))DEALLOCATE(work3d)
       ALLOCATE(work3d(ir,ic,i3))
       no_4dr_l = izero

       CALL MPI_unpack(mpi_buffer,buff_size,position,work3d,count,MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD,mpiierr) ! hicmz

       no_4dr_l = no_4dr_l +1

       DO j=1,i3               ! I dont like mpi stride routines
         DO k=1,ic
           DO l =1,ir
             hicmz(l,k,index,j) =  work3d(l,k,j) 
           ENDDO
         ENDDO
       ENDDO


       CALL MPI_unpack(mpi_buffer,buff_size,position,no_0di_s,one,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr)
       no_0di_l = no_0di_l +1

       CALL MPI_unpack(mpi_buffer,buff_size,position,no_2di_s,one,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr)
       no_0di_l = no_0di_l +1

       CALL MPI_unpack(mpi_buffer,buff_size,position,no_2di_s2,one,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr)
       no_0di_l = no_0di_l +1


       CALL MPI_unpack(mpi_buffer,buff_size,position,no_3dr_p_s,one,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr)
       no_0di_l = no_0di_l +1

       CALL MPI_unpack(mpi_buffer,buff_size,position,no_2dr_s,one,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr)
       no_0di_l = no_0di_l +1

       CALL MPI_unpack(mpi_buffer,buff_size,position,no_3dr_s,one,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr)
       no_0di_l = no_0di_l +1

       CALL MPI_unpack(mpi_buffer,buff_size,position,no_4dr_s,one,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpiierr)
       no_0di_l = no_0di_l +1

       ! check for errors in number of items to be transmitted
       IF((no_0di_l .NE. no_0di_s) .OR. (no_0di_l .NE. no_0di) &
                                    .OR. (no_0di_s .NE. no_0di)) THEN
         WRITE(ncrt,FMT='("error in unpacking 0di data,sub nfreya_store_packed")')
         WRITE(ncrt,10)no_0di,no_0di_l,no_0di_s
10       FORMAT(2x,"no_0di,no_0di_l,no_0di_s =",3(2x,i5))
         lerrno = iomaxerr + 234_I4B
         CALL terminate(lerrno,nlog)
       ENDIF

       IF((no_2di_l .NE. no_2di_s) .OR. (no_2di_l .NE. no_2di) &
                                    .OR. (no_2di_s .NE. no_2di)) THEN
         WRITE(ncrt,FMT='("error in unpacking 1di data,sub nfreya_store_packed")')
         WRITE(ncrt,15)no_2di,no_2di_l,no_2di_s
15       FORMAT(2x,"no_2di,no_2di_l,no_2di_s =",3(2x,i5))
         lerrno = iomaxerr + 235_I4B
         CALL terminate(lerrno,nlog)
       ENDIF

       IF((no_2di_l2 .NE. no_2di_s2) .OR. (no_2di_l2 .NE. no_2di2) &
                                    .OR. (no_2di_s2 .NE. no_2di2)) THEN
         WRITE(ncrt,FMT='("error in unpacking 2di data,sub nfreya_store_packed")')
         WRITE(ncrt,17)no_2di,no_2di_l2,no_2di_s2
17       FORMAT(2x,"no_2di,no_2di_l2,no_2di_s2 =",3(2x,i5))
         lerrno = iomaxerr + 235_I4B
         CALL terminate(lerrno,nlog)
       ENDIF

       IF((no_3dr_p_l .NE. no_3dr_p_s) .OR. (no_3dr_p_l .NE. no_3dr_p) &
                                    .OR. (no_3dr_p_s .NE. no_3dr_p)) THEN
         WRITE(ncrt,FMT='("error in unpacking 1dr data,sub nfreya_store_packed")')
         WRITE(ncrt,20)no_3dr_p,no_3dr_p_l,no_3dr_p_s
20       FORMAT(2x,"no_3dr_p,no_3dr_p_l,no_3dr_p_s =",3(2x,i5))
         lerrno = iomaxerr + 236_I4B
         CALL terminate(lerrno,nlog)
       ENDIF

       IF((no_2dr_l .NE. no_2dr_s) .OR. (no_2dr_l .NE. no_2dr) &
                                    .OR. (no_2dr_s .NE. no_2dr)) THEN
         WRITE(ncrt,FMT='("error in unpacking 2dr data,sub nfreya_store_packed")')
         WRITE(ncrt,25)no_2dr,no_2dr_l,no_2dr_s
25       FORMAT(2x,"no_2dr,no_2dr_l,no_2dr_s =",3(2x,i5))
         lerrno = iomaxerr + 237_I4B
         CALL terminate(lerrno,nlog)
       ENDIF

      IF((no_3dr_l .NE. no_3dr_s) .OR. (no_3dr_l .NE. no_3dr) &
                                    .OR. (no_3dr_s .NE. no_3dr)) THEN
         WRITE(ncrt,FMT='("error in unpacking 3dr data,sub nfreya_store_packed")')
         WRITE(ncrt,30)no_3dr,no_3dr_l,no_3dr_s
30       FORMAT(2x,"no_3dr,no_3dr_l,no_3dr_s =",3(2x,i5))
         lerrno = iomaxerr + 238_I4B
         CALL terminate(lerrno,nlog)
      ENDIF

      IF((no_4dr_l .NE. no_4dr_s) .OR. (no_4dr_l .NE. no_4dr) &
                                    .OR. (no_4dr_s .NE. no_4dr)) THEN
         WRITE(ncrt,FMT='("error in unpacking 4dr data,sub nfreya_store_packed")')
         WRITE(ncrt,35)no_4dr,no_4dr_l,no_4dr_s
35       FORMAT(2x,"no_4dr,no_4dr_l,no_4dr_s =",3(2x,i5))
         lerrno = iomaxerr + 239_I4B
         CALL terminate(lerrno,nlog)
      ENDIF

       RETURN
     END SUBROUTINE nfreya_store_packed
#endif



    END MODULE nfreya_run
