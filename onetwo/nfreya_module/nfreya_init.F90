  SUBROUTINE nfreya_init
!----------------------------------------------------------------------------
! -- initialize  beams and zonal quantities and cross sections
! -- define allocatable arrays. Note that namelist data arrays 
! -- are allocateable however.
!-------------------------------------------------------------HSJ 1/14/2011--
    USE nrtype,                                     ONLY : DP,I4B

    USE common_constants,                           ONLY : izero,Rad_per_Deg

    USE error_handler

    USE prep_zones,                                 ONLY : prenub

    USE neutral_beams,                              ONLY : iterate_beam,fast_ion_target,   &
                                                           mi,mj,npat,drpat,dzpat,nrpat,   &
                                                           nzpat, inubpat,no_injectors,    &
                                                           split_injectors,time_dep_beam,  &
                                                           no_physical_injectors 

    USE xsct,                                       ONLY : setup_xsct 

    USE tport_mhd_grid_data,                        ONLY : rmhdgrid, zmhdgrid

    USE Nfreya_routines,                            ONLY : beam_prop

    USE nfreya_run,                                 ONLY : codeid

    USE MPI_data,                                   ONLY : myid,numprocs

    USE io_gcnmp,                                   ONLY : nlog

    IMPLICIT NONE
    INTEGER(I4B) ib,icall_prenub,task,error
  
    DATA icall_prenub /1/

      fast_ion_target =0
      IF(iterate_beam)fast_ion_target = 1


      time_dep_beam = -1
      IF(numprocs == 1)THEN 
         codeid = '2d_single_proc'
      ELSE
         codeid = '2d_multiple_proc'
      ENDIF

!      task = izero  ; error = izero
!      CALL setup_plot_file(task,error)              ! creates P_NF_bpltfil 
                                                     ! moved to  nfreya_load_single_source
                                                     ! because  inpt namelist may change
                                                     ! due to redefinitions casued by nubeam
                                                     ! namelist read
!      IF(error .NE. izero)THEN
!         lerrno = 4_I4B
!         CALL terminate(lerrno,nlog)
!      ENDIF
 





      CALL prenub(icall_prenub)                     ! prep_zonal_input.f90
      icall_prenub = icall_prenub+1
      mi = SIZE(rmhdgrid) ; mj = SIZE(zmhdgrid)     ! mi,mj can be subset of nw,nh

      IF (inubpat .EQ. 1)THEN
         npat(1) = mi ; npat(2) = mj
         CALL setrz (npat,rmhdgrid(1),rmhdgrid(mi),zmhdgrid(1),zmhdgrid(mj),&
                     drpat,dzpat,nrpat,nzpat)
      ENDIF
!

      no_physical_injectors = no_injectors
      CALL sort_injectors      ! by beam times,powers,etc.
      IF(numprocs .gt. no_injectors .AND. split_injectors ) call reset_no_injectors
 

!     allocate beam arrays (except vx_izpt,vy_izpt,vz_izpt,x_izpt,y_izpt,
!     z_izpt,r_izpt,pitch_izpt which are done in sub birth_point_allocation) 
!     These arrays are allocated after we get the number
!     of psuedo particles in sub beam_prop

      CALL beam_allocation

!     calculate beam power; account for neutralizer efficiency
!     allocate some arrays needed in setup_xsct,etc.
      CALL beam_prop

      CALL birth_point_allocation

      IF(numprocs .GT. 1)THEN                 ! must be called after beam_allocation
         !     setup send/receive buffer,mpi_buffer
#ifdef USEMPI
               CALL allocate_mpi_buffer
#endif
      ENDIF


! determine the cross section sets needed

      CALL setup_xsct !  called by all but only master does calcs (xsct_class.F90)
                      !  distribute_xsct passes out data (xsct_class.F90)


      RETURN
  END SUBROUTINE nfreya_init


  SUBROUTINE  beam_allocation
! ------------------------------------------------------------------------------------------
! -- setup allocatable arrays
! ------------------------------------------------------------------------------------------
   USE nrtype,                                     ONLY :   DP,I4B

   USE nf_param,                                   ONLY :   ke,kcm

   USE neutral_beams,                              ONLY:    cangv,cangh,sangv,sangh,        &
                                                            anglev,angleh,no_injectors,     &
                                                            ebeam,vbeam,pbeam,bion,bneut,   &
                                                            forb, fb11,fb10, fb01, fb00,    &
                                                            fber,no_physical_injectors,fap, &
                                                            bencap, fbe, fbi,fbth,          &
                                                            bke, bki,enb, enbsav,hicmz,     &
                                                            enbav,hdep,hdepz, pb0,          &
                                                            ftrapfi,angmpz,vx_izpt,vy_izpt, &
                                                            vz_izpt,x_izpt,y_izpt,z_izpt,   &
                                                            r_izpt,pitch_a,nmbrz,           &
                                                            ppb, ppbsav, ppbav,qbsav,       &
                                                            sbsav,spbsav, taupb,            &
                                                            tauppb, taueb, wb, wbsav,       &
                                                            wbav,zetaz,hibr,hibrz,          &
                                                            wb11,wb10,wb01,wb00,fwall,      &
                                                            no_0di,                         &
                                                            no_2di,no_2dr,no_3dr,no_4dr,    &
                                                            npart_all_beamlines,no_2di2,    &
                                                            calc_variance,nsample_izpt

   USE zonal_data,                                 ONLY :  mf,mfm1

   USE common_constants,                           ONLY :  izero,zeroc,Rad_per_Deg

   USE grid_class,                                 ONLY :  nj

   USE Plasma_properties,                          ONLY :  neut_beam

   IMPLICIT NONE         
   INTEGER(I4B) ib

      IF(.NOT. ALLOCATED(cangv))ALLOCATE(cangv(no_injectors))
      IF(.NOT. ALLOCATED(cangh))ALLOCATE(cangh(no_injectors))
      IF(.NOT. ALLOCATED(sangv))ALLOCATE(sangv(no_injectors))
      IF(.NOT. ALLOCATED(sangh))ALLOCATE(sangh(no_injectors))
      DO ib = 1,no_injectors
        cangv(ib) = COS (anglev(ib)*Rad_per_Deg)
        cangh(ib) = COS (angleh(ib)*Rad_per_Deg)
        sangv(ib) = SIN (anglev(ib)*Rad_per_Deg)
        sangh(ib) = SIN (angleh(ib)*Rad_per_Deg)
      END DO

      !no_0di is set in allocate_mpi_buffer
      no_2di = izero                   ! counts arrays(which are the same size)
                                       ! to be put in buffer for message passing

      no_2di2 = izero                  ! same for oter size integer arrays (eg nsample_izpt)

      IF(ALLOCATED(nmbrz))DEALLOCATE(nmbrz)
      ALLOCATE(nmbrz(mfm1,no_injectors))
      no_2di = no_2di +1
      IF(.NOT. ASSOCIATED(neut_beam%nmbrz))ALLOCATE(neut_beam%nmbrz(mfm1,no_physical_injectors))
      nmbrz(:,:)           = izero
      neut_beam%nmbrz(:,:) = izero
      
     ! NOTE:  for 1dr arrays see sub birth_point_allocation 

      no_2dr = izero  ! counts arrays to be put in buffer for message passing

      IF(ALLOCATED(nsample_izpt))DEALLOCATE(nsample_izpt)
      ALLOCATE(nsample_izpt(ke,no_injectors))
      no_2di2 = no_2di2 +1
      IF(.NOT. ASSOCIATED(neut_beam%nsample_izpt))ALLOCATE(neut_beam%nsample_izpt(ke,no_physical_injectors))
      nsample_izpt(:,:)           = izero
      neut_beam%nsample_izpt(:,:) = izero
      
     ! NOTE:  for 1dr arrays see sub birth_point_allocation 

      no_2dr = izero  ! counts arrays to be put in buffer for message passing




      IF(ALLOCATED(bion)) DEALLOCATE(bion)
      ALLOCATE(bion(ke,no_injectors))
      bion(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%bion))ALLOCATE(neut_beam%bion(ke,no_physical_injectors))
      neut_beam%bion(:,:) = zeroc


      IF(ALLOCATED(bneut)) DEALLOCATE(bneut)
      ALLOCATE(bneut(ke,no_injectors))
      bneut(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%bneut))ALLOCATE(neut_beam%bneut(ke,no_physical_injectors))
      neut_beam%bneut(:,:) = zeroc


      IF(.NOT. ASSOCIATED(neut_beam%fbcur))ALLOCATE(neut_beam%fbcur(ke,no_physical_injectors))
      neut_beam%fbcur(:,:) = zeroc ! 88888899999
 
      IF(ALLOCATED(pbeam)) DEALLOCATE(pbeam)
      ALLOCATE(pbeam(ke,no_injectors))
      pbeam(:,:) = zeroc
      no_2dr = no_2dr + 1

      IF(.NOT. ASSOCIATED(neut_beam%pbeam))ALLOCATE(neut_beam%pbeam(ke,no_physical_injectors))
      neut_beam%pbeam(:,:) = zeroc
 
      IF(ALLOCATED(fap)) DEALLOCATE(fap)
      ALLOCATE(fap(ke,no_injectors))
      fap(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%fap))ALLOCATE(neut_beam%fap(ke,no_physical_injectors))
      neut_beam%fap(:,:) = zeroc


      IF(ALLOCATED(fwall)) DEALLOCATE(fwall)
      ALLOCATE(fwall(ke,no_injectors))
      fwall(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%fwall))ALLOCATE(neut_beam%fwall(ke,no_physical_injectors))
      neut_beam%fwall(:,:) = zeroc


      IF(ALLOCATED(forb)) DEALLOCATE(forb)
      ALLOCATE(forb(ke,no_injectors))
      forb(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%forb))ALLOCATE(neut_beam%forb(ke,no_physical_injectors))
      neut_beam%forb(:,:) = zeroc






      IF(ALLOCATED(fb11)) DEALLOCATE(fb11)
      ALLOCATE(fb11(ke,no_injectors))
      fb11(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%fb11))ALLOCATE(neut_beam%fb11(ke,no_physical_injectors))
      neut_beam%fb11(:,:)= zeroc

      IF(ALLOCATED(fb10)) DEALLOCATE(fb10)
      ALLOCATE(fb10(ke,no_injectors))
      fb10(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%fb10))ALLOCATE(neut_beam%fb10(ke,no_physical_injectors))
      neut_beam%fb10(:,:)= zeroc

      IF(ALLOCATED(fb01)) DEALLOCATE(fb01)
      ALLOCATE(fb01(ke,no_injectors))
      fb01(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%fb01))ALLOCATE(neut_beam%fb01(ke,no_physical_injectors))
      neut_beam%fb01(:,:)= zeroc

      IF(ALLOCATED(fb00)) DEALLOCATE(fb00)
      ALLOCATE(fb00(ke,no_injectors))
      fb00(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%fb00))ALLOCATE(neut_beam%fb00(ke,no_physical_injectors))
      neut_beam%fb00(:,:)= zeroc

      IF(ALLOCATED(fber)) DEALLOCATE(fber)
      ALLOCATE(fber(ke,no_injectors))
      fber(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%fber))ALLOCATE(neut_beam%fber(ke,no_physical_injectors))
      neut_beam%fber(:,:)= zeroc





      IF(ALLOCATED(wb11)) DEALLOCATE(wb11)
      ALLOCATE(wb11(ke,no_injectors))
      wb11(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%wb11))ALLOCATE(neut_beam%wb11(ke,no_physical_injectors))
      neut_beam%wb11(:,:) = zeroc

      IF(ALLOCATED(wb10)) DEALLOCATE(wb10)
      ALLOCATE(wb10(ke,no_injectors))
      wb10(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%wb10))ALLOCATE(neut_beam%wb10(ke,no_physical_injectors))
      neut_beam%wb10(:,:) = zeroc

      IF(ALLOCATED(wb01)) DEALLOCATE(wb01)
      ALLOCATE(wb01(ke,no_injectors))
      wb01(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%wb01))ALLOCATE(neut_beam%wb01(ke,no_physical_injectors))
      neut_beam%wb01(:,:) = zeroc

      IF(ALLOCATED(wb00)) DEALLOCATE(wb00)
      ALLOCATE(wb00(ke,no_injectors))
      wb00(:,:) = zeroc
      no_2dr = no_2dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%wb00))ALLOCATE(neut_beam%wb00(ke,no_physical_injectors))
      neut_beam%wb00(:,:) = zeroc

! end of list of 2d arrays to be bufferd (in mpi_buffer)


      IF(ALLOCATED(ebeam)) DEALLOCATE(ebeam)
      ALLOCATE(ebeam(ke,no_injectors))
      ebeam(:,:) = zeroc
      IF(.NOT. ASSOCIATED(neut_beam%ebeam))ALLOCATE(neut_beam%ebeam(ke,no_physical_injectors))

      IF(ALLOCATED(vbeam)) DEALLOCATE(vbeam)
      ALLOCATE(vbeam(ke,no_injectors))
      vbeam(:,:) = zeroc
      IF(.NOT. ASSOCIATED(neut_beam%vbeam))ALLOCATE(neut_beam%vbeam(ke,no_physical_injectors))
      neut_beam%vbeam(:,:) = zeroc






      !----------------------------------------------------------------------
      !no_3d  counts arrays to be put in buffer for message passing
      !----------------------------------------------------------------------
      no_3dr = izero
      IF(ALLOCATED(ftrapfi)) DEALLOCATE(ftrapfi)
      ALLOCATE(ftrapfi(mfm1,ke,no_injectors))
      ftrapfi(:,:,:) = zeroc
      no_3dr = no_3dr +1
      IF(.NOT. ASSOCIATED(neut_beam%ftrapfi))ALLOCATE(neut_beam%ftrapfi(mfm1,ke,no_injectors))
      neut_beam%ftrapfi(:,:,:) = zeroc


      IF(ALLOCATED(hibrz)) DEALLOCATE(hibrz)
      ALLOCATE(hibrz(mfm1,ke,no_injectors))
      hibrz(:,:,:) = zeroc
      no_3dr = no_3dr +1
      IF(.NOT. ASSOCIATED(neut_beam%hibrz))ALLOCATE(neut_beam%hibrz(mfm1,ke,no_injectors))
      neut_beam%hibrz(:,:,:) = zeroc

      IF(ALLOCATED(hdepz)) DEALLOCATE(hdepz)
      ALLOCATE(hdepz(mfm1,ke,no_injectors))
      hdepz(:,:,:) = zeroc
      no_3dr = no_3dr +1
      IF(.NOT. ASSOCIATED(neut_beam%hdepz))ALLOCATE(neut_beam%hdepz(mfm1,ke,no_injectors))
      neut_beam%hdepz(:,:,:) = zeroc

      IF(ALLOCATED(zetaz)) DEALLOCATE(zetaz)
      ALLOCATE(zetaz(mfm1,ke,no_injectors))
      zetaz(:,:,:) = zeroc
      no_3dr = no_3dr +1
      IF(.NOT. ASSOCIATED(neut_beam%zetaz))ALLOCATE(neut_beam%zetaz(mfm1,ke,no_injectors))
      neut_beam%zetaz(:,:,:) = zeroc

      IF(ALLOCATED(angmpz)) DEALLOCATE(angmpz)
      ALLOCATE(angmpz(mfm1,ke,no_injectors))
      angmpz(:,:,:) = zeroc
      no_3dr = no_3dr +1
      IF(.NOT. ASSOCIATED(neut_beam%angmpz))ALLOCATE(neut_beam%angmpz(mfm1,ke,no_injectors))
      neut_beam%angmpz(:,:,:) = zeroc


! end of list of 3dq arrays to be bufferd (in mpi_buffer)




! nj arrays exist on transport grid, not nfreya zonal grid:

      IF(ALLOCATED(bencap)) DEALLOCATE(bencap)
      ALLOCATE(bencap(nj,ke,no_injectors))
      bencap(:,:,:) = zeroc

      IF(ALLOCATED(fbe)) DEALLOCATE(fbe)
      ALLOCATE(fbe(nj,ke,no_injectors))
      fbe(:,:,:) = zeroc

      IF(ALLOCATED(fbi)) DEALLOCATE(fbi)
      ALLOCATE(fbi(nj,ke,no_injectors))
      fbi(:,:,:) = zeroc

      IF(ALLOCATED(fbth)) DEALLOCATE(fbth)
      ALLOCATE(fbth(nj,ke,no_injectors))
      fbth(:,:,:) = zeroc

      IF(ALLOCATED(bke)) DEALLOCATE(bke)
      ALLOCATE(bke(nj,ke,no_injectors))
      bke(:,:,:) = zeroc

      IF(ALLOCATED(bki)) DEALLOCATE(bki)
      ALLOCATE(bki(nj,ke,no_injectors))
      bki(:,:,:) = zeroc

      IF(ALLOCATED(enb)) DEALLOCATE(enb)
      ALLOCATE(enb(nj,ke,no_injectors))
      enb(:,:,:) = zeroc

      IF(ALLOCATED(enbsav)) DEALLOCATE(enbsav)
      ALLOCATE(enbsav(nj,ke,no_injectors))
      enbsav(:,:,:) = zeroc

      IF(ALLOCATED(enbav)) DEALLOCATE(enbav)
      ALLOCATE(enbav(nj,ke,no_injectors))
      enbav(:,:,:) = zeroc

      IF(ALLOCATED(hdep)) DEALLOCATE(hdep)
      ALLOCATE(hdep(nj,ke,no_injectors))
      hdep(:,:,:) = zeroc


      IF(ALLOCATED(hibr)) DEALLOCATE(hibr)
      ALLOCATE(hibr(nj,ke,no_injectors))
      hibr(:,:,:) = zeroc


      IF(ALLOCATED(pb0)) DEALLOCATE(pb0)
      ALLOCATE(pb0(nj,ke,no_injectors))
      pb0(:,:,:) = zeroc


      IF(ALLOCATED(ppb)) DEALLOCATE(ppb)
      ALLOCATE(ppb(nj,ke,no_injectors))
      ppb(:,:,:) = zeroc

      IF(ALLOCATED(ppbsav)) DEALLOCATE(ppbsav)
      ALLOCATE(ppbsav(nj,ke,no_injectors))
      ppbsav(:,:,:) = zeroc

      IF(ALLOCATED(ppbav)) DEALLOCATE(ppbav)
      ALLOCATE(ppbav(nj,ke,no_injectors))
      ppbav(:,:,:) = zeroc

      IF(ALLOCATED(qbsav)) DEALLOCATE(qbsav)
      ALLOCATE(qbsav(nj,ke,no_injectors))
      qbsav(:,:,:) = zeroc

      IF(ALLOCATED(sbsav)) DEALLOCATE(sbsav)
      ALLOCATE(sbsav(nj,ke,no_injectors))
      sbsav(:,:,:) = zeroc

      IF(ALLOCATED(spbsav)) DEALLOCATE(spbsav)
      ALLOCATE(spbsav(nj,ke,no_injectors))
      spbsav(:,:,:) = zeroc

      IF(ALLOCATED(taupb)) DEALLOCATE(taupb)
      ALLOCATE(taupb(nj,ke,no_injectors))
      taupb(:,:,:) = zeroc

      IF(ALLOCATED(tauppb)) DEALLOCATE(tauppb)
      ALLOCATE(tauppb(nj,ke,no_injectors))
      tauppb(:,:,:) = zeroc

      IF(ALLOCATED(taueb)) DEALLOCATE(taueb)
      ALLOCATE(taueb(nj,ke,no_injectors))
      taueb(:,:,:) = zeroc

      IF(ALLOCATED(wb)) DEALLOCATE(wb)
      ALLOCATE(wb(nj,ke,no_injectors))
      wb(:,:,:) = zeroc

      IF(ALLOCATED(wbsav)) DEALLOCATE(wbsav)
      ALLOCATE(wbsav(nj,ke,no_injectors))
      wbsav(:,:,:) = zeroc

      IF(ALLOCATED(wbav)) DEALLOCATE(wbav)
      ALLOCATE(wbav(nj,ke,no_injectors))
      wbav(:,:,:) = zeroc


   
      !----------------------------------------------------------------------
      !no_4d  counts arrays to be put in buffer for message passing
      !----------------------------------------------------------------------
      no_4dr = izero
      IF(ALLOCATED(hicmz)) DEALLOCATE(hicmz)
      ALLOCATE(hicmz(mfm1,ke,no_injectors,kcm))
      no_4dr = no_4dr + 1
      IF(.NOT. ASSOCIATED(neut_beam%hicmz))ALLOCATE(neut_beam%hicmz(mfm1,ke,no_injectors,kcm))
      neut_beam%hicmz(:,:,:,:) = zeroc


! end of list of 4d arrays to be bufferd (in mpi_buffer)

  
!      IF(calc_variance)THEN  ! always allocate and set to zero
 
         ALLOCATE(neut_beam%stdev_nmbrz(mfm1,no_physical_injectors))
         neut_beam%stdev_nmbrz(:,:) = izero

         ALLOCATE(neut_beam%stdev_bion(ke,no_physical_injectors))
         neut_beam%stdev_bion(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_bneut(ke,no_physical_injectors))
         neut_beam%stdev_bneut(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_fap(ke,no_physical_injectors))
         neut_beam%stdev_fap(:,:)  = zeroc
 
         ALLOCATE(neut_beam%stdev_fwall(ke,no_physical_injectors))
         neut_beam%stdev_fwall(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_forb(ke,no_physical_injectors))
         neut_beam%stdev_forb(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_fber(ke,no_physical_injectors))
         neut_beam%stdev_fber(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_fb11(ke,no_physical_injectors))
         neut_beam%stdev_fb11(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_fb00(ke,no_physical_injectors))
         neut_beam%stdev_fb00(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_fb01(ke,no_physical_injectors))
         neut_beam%stdev_fb01(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_fb10(ke,no_physical_injectors))
         neut_beam%stdev_fb10(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_wb11(ke,no_physical_injectors))
         neut_beam%stdev_wb11(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_wb00(ke,no_physical_injectors))
         neut_beam%stdev_wb00(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_wb01(ke,no_physical_injectors))
         neut_beam%stdev_wb01(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_wb10(ke,no_physical_injectors))
         neut_beam%stdev_wb10(:,:)  = zeroc

         ALLOCATE(neut_beam%stdev_ftrapfi(mfm1,ke,no_physical_injectors))
         neut_beam%stdev_ftrapfi(:,:,:) = zeroc

         ALLOCATE(neut_beam%stdev_hibrz(mfm1,ke,no_physical_injectors))
         neut_beam%stdev_hibrz(:,:,:) = zeroc

         ALLOCATE(neut_beam%stdev_hdepz(mfm1,ke,no_physical_injectors))
         neut_beam%stdev_hdepz(:,:,:) = zeroc

         ALLOCATE(neut_beam%stdev_zetaz(mfm1,ke,no_physical_injectors))
         neut_beam%stdev_zetaz(:,:,:) = zeroc

         ALLOCATE(neut_beam%stdev_angmpz(mfm1,ke,no_physical_injectors))
         neut_beam%stdev_angmpz(:,:,:) = zeroc

         ALLOCATE(neut_beam%stdev_hicmz(mfm1,ke,no_injectors,kcm))
         neut_beam%stdev_hicmz(:,:,:,:) = zeroc

!      ENDIF

 
      RETURN

  END SUBROUTINE beam_allocation

  SUBROUTINE birth_point_allocation
! ------------------------------------------------------------------------------------------
! --these arrays drive the mpi packing/unpacking time s we want them
! -- to be as small as posible. To avoid intricate codeing we make them
! -- the same size on all injectors
! ------------------------------------------------------------------------------------------  
   USE nrtype,                                     ONLY :   DP,I4B

   USE nf_param,                                   ONLY :   ke

   USE MPI_data,                                   ONLY :   numprocs

   USE neutral_beams,                              ONLY:    vx_izpt,vy_izpt,     &
                                                            vz_izpt,x_izpt,      &
                                                            y_izpt,z_izpt,       &
                                                            r_izpt,pitch_a,      &
                                                            no_3dr_p,            &
                                                            nsample_izpt,        &
                                                            npart_pseudo,        &
                                                            no_birth_points,     &
                                                            npart_all_beamlines, &
                                                            no_injectors

   USE common_constants,                           ONLY :   izero,zeroc

      IMPLICIT NONE

      INTEGER(I4B) maxbirth

      !-------------------------------------------------------------------------
      ! -- in the single cpu case there is no message passing.
      ! -- hence only the **_izpt arrays on master are required.
      ! -- these can be large and must account for all beamlines,energies
      ! -- In the multiple cpu  case the **_izpt arrays  must be passed
      ! -- back to the master and contain only the current beam line.
      ! -- The master then collects the data
      ! -- Note that no_injectors = no_physical injectors 
      ! --                          if numprocs <= no_physical injectors
      ! -- If if numprocs >  no_physical injectors then no_injectors = numprocs
      !-------------------------------------------------------------------------

      maxbirth = MAXVAL(npart_pseudo)  ! ke energies per beamline
      maxbirth = maxbirth/1000.        ! keep one in every 1000
      maxbirth = 200                   ! make this an input ?
      no_birth_points = maxbirth

      no_3dr_p = izero           ! counts arrays(which are the same size)
      IF(ALLOCATED(vx_izpt))DEALLOCATE(vx_izpt)
      ALLOCATE(vx_izpt(no_birth_points,ke,no_injectors))
      no_3dr_p = no_3dr_p +1

      IF(ALLOCATED(vy_izpt))DEALLOCATE(vy_izpt)
      ALLOCATE(vy_izpt(no_birth_points,ke,no_injectors))
      no_3dr_p = no_3dr_p +1

      IF(ALLOCATED(vz_izpt))DEALLOCATE(vz_izpt)
      ALLOCATE(vz_izpt(no_birth_points,ke,no_injectors))
      no_3dr_p = no_3dr_p +1

      IF(ALLOCATED(x_izpt))DEALLOCATE(x_izpt)
      ALLOCATE(x_izpt(no_birth_points,ke,no_injectors))
      no_3dr_p = no_3dr_p +1

      IF(ALLOCATED(y_izpt))DEALLOCATE(y_izpt)
      ALLOCATE(y_izpt(no_birth_points,ke,no_injectors))
      no_3dr_p = no_3dr_p +1

      IF(ALLOCATED(z_izpt))DEALLOCATE(z_izpt)
      ALLOCATE(z_izpt(no_birth_points,ke,no_injectors))
      no_3dr_p = no_3dr_p +1

      IF(ALLOCATED(r_izpt))DEALLOCATE(r_izpt)
      ALLOCATE(r_izpt(no_birth_points,ke,no_injectors))
      no_3dr_p = no_3dr_p +1

      IF(ALLOCATED(pitch_a))DEALLOCATE(pitch_a)
      ALLOCATE(pitch_a(no_birth_points,ke,no_injectors))
      no_3dr_p = no_3dr_p +1

      IF(ALLOCATED(nsample_izpt))DEALLOCATE(nsample_izpt)
      ALLOCATE(nsample_izpt(ke,no_injectors))


      RETURN

  END SUBROUTINE birth_point_allocation




#ifdef USEMPI
       SUBROUTINE  allocate_mpi_buffer
! ------------------------------------------------------------------------------------------
! -- create allocatable send/receive buffer , mpi_buffer(buff_size)
! -- NOTE: ROUTINES THAT MUST BE KEPT IN SYNC:
! --       nfreya_pack_data
! --       nfreya_update_master
! --       nfreya_store_packed 
! --       allocate_mpi_buffer
! ------------------------------------------------------------------------------------------

         USE nrtype,                                   ONLY :  DP,I4B

         USE MPI_data,                                 ONLY :  mpi_buffer,buff_size,mpiierr

         USE MPI

         USE common_constants,                         ONLY : izero

         USE neutral_beams,                            ONLY :                              &
                                                               vbeam,cangv ,cangh ,sangv,  &
                                                               sangh ,thetp , npskip,      &
                                                               thetpp ,costp ,sintp ,      &
                                                               costpp ,sintpp,nsourc,      &
                                                               iatype,kt,npulse,pbeamon,   &
                                                               pbeamOff,beam_on,nmbrz,      &
                                                               beam_end,source2_phase,     &
                                                               nap,nashape,bvofset,        &
                                                               bhofset, vx_izpt,           &
                                                               fb11,fb01,fb10,fb00,fber,   &
                                                               wb11,wb00,wb01,wb10,bleni,  &
                                                               ftrapfi,hibrz,hdepz,angmpz, &
                                                               zetaz,hicmz,olossc,nbshape, &
                                                               nashape,bheigh,bwidth,bhfoc,&
                                                               bvfoc,bhdiv,bvdiv,sfrac1,   &
                                                               naptr,aheigh,awidth,alen,   &
                                                               blenp,rpivot,zpivot,no_0di, &
                                                               no_3dr_p,no_2di,no_2dr,     &
                                                               no_4dr,npart_all_beamlines, &
                                                               no_2di2,nsample_izpt,no_3dr

         IMPLICIT NONE

         INTEGER(I4B) ir,ic,i3,count,int_size,                &
                      double_size,testing,one,tag



             buff_size = izero ; one = 1_I4B

             !--------------------------------------------------------------------------
             ! get buff size of integers (eq 0d arrays)  to be sent:
             !--------------------------------------------------------------------------
             no_0di = izero
             ir = one
             CALL MPI_PACK_SIZE(ir,MPI_INTEGER,MPI_COMM_WORLD,int_size,mpiierr)
             buff_size = buff_size + int_size                    ! for injctr
             no_0di = no_0di + 1
             ir = one
             CALL MPI_PACK_SIZE(ir,MPI_INTEGER,MPI_COMM_WORLD,int_size,mpiierr)
             buff_size = buff_size + int_size                    ! for n_izpt
             no_0di = no_0di + 1

             buff_size = buff_size + 7*int_size  !7 = {no_0_di,no_3dr_p,no_2di,no_2dr,no_3dr,no_4dr,no_2di2} 
             no_0di = no_0di + 7
 
             !--------------------------------------------------------------------------
             ! get buff size of integer  1d arrays(of size mfm1)  to be sent:
             !--------------------------------------------------------------------------
             ir = size(nmbrz)
             ! note no_2di was set in sub beam_allocation
             CALL MPI_PACK_SIZE(ir,MPI_INTEGER,MPI_COMM_WORLD,int_size,mpiierr)
             buff_size = buff_size + no_2di*int_size

             ir = size(nsample_izpt,1)
             CALL MPI_PACK_SIZE(ir,MPI_INTEGER,MPI_COMM_WORLD,int_size,mpiierr)
             buff_size = buff_size + no_2di2*int_size

             !--------------------------------------------------------------------------
             ! get buff size for  real  3d arrays(of common size (no_birth_points,ke,no_injectors))  to be sent:
             !--------------------------------------------------------------------------
             ir = size(vx_izpt,1) ; ic = SIZE(vx_izpt,2)
             count = ir*ic
             CALL MPI_PACK_SIZE(count,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD, &
                                double_size,mpiierr)
             buff_size = buff_size + no_3dr_p * double_size




             !--------------------------------------------------------------------------
             ! only column injctr  of 2d arrays  is filled by slaves  so 
             ! send only this column ==>fb11(1,injctr) for example
             ! see sub nfreya_pack data for list of arrays to be sent
             ! there are no_2d arrays of this type to be sent
             ! no_2d was set in beam_allocation
             !--------------------------------------------------------------------------
             ir = SIZE(fb11,1)
             count = ir
             CALL MPI_PACK_Size(count,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,double_size,mpiierr)
             ! there are no_2d arrays of this type to be sent
             buff_size = buff_size + no_2dr*double_size  




            !--------------------------------------------------------------------------
            ! -- 3d arrays
            !--------------------------------------------------------------------------
             ir = SIZE(ftrapfi,1)
             ic = SIZE(ftrapfi,2)
             count = ir*ic
             CALL MPI_PACK_Size(count,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,double_size,mpiierr)
             ! there are no_3d arrays of this type to be sent
             buff_size = buff_size + no_3dr*double_size  


            !--------------------------------------------------------------------------
            ! -- 4d arrays
            !--------------------------------------------------------------------------
             ir = SIZE(hicmz,1) ; ic = SIZE(hicmz,2) ; i3 = SIZE(hicmz,4) 
             count = ir*ic*i3
             CALL MPI_PACK_Size(count,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,double_size,mpiierr)
             buff_size = buff_size + no_4dr*double_size  



         IF(ALLOCATED(mpi_buffer))DEALLOCATE(mpi_buffer)
         ALLOCATE(mpi_buffer(buff_size))
      
         RETURN
         
       END SUBROUTINE allocate_mpi_buffer
#endif




      SUBROUTINE setrz (ndim, rmin, rmax, zmin, zmax, dr, dz, nr, nz)
! ----------------------------------------------------------------------
! this subroutine sets up (r,z) grid parameters for FREYA when the
! optional two-dimensional deposition calculation is used.
! default parameters correspond to original eqdsk grid.
! ----------------------------------------------------------------------
!

      USE nrtype,                              ONLY : Dp,I4B
      IMPLICIT NONE

      REAL(DP) dr,dz,rmin,rmax,zmin,zmax
      INTEGER(I4B)  ndim(2),nr,nz
!
      nr =  ndim(1)
      nz =  ndim(2)
      dr = (rmax-rmin)/(nr-1)
      dz = (zmax-zmin)/(nz-1)
      RETURN
!
      END SUBROUTINE setrz



      SUBROUTINE      setup_plot_file(task,error)
! ----------------------------------------------------------------------
! -- intialize file to store data for NF_nubplt
! -- file is handled only by master process
! ----------------------------------------------------------------------
      USE nrtype,                            ONLY : DP,I4B,I2B

      USE neutral_beams,                     ONLY : nfreya_plot_unit,          &
                                                    nfreya_plot_file_name,     &
                                                    no_physical_injectors,iborb,npitch,       &
                                                    atw_beam,npitch

      USE Nfreya_routines,                   ONLY : NF_nubplt

      USE Plasma_properties ,                ONLY : mhd_dat,dischg                                                    

      USE MPI_data,                          ONLY : myid,master

      USE io_gcnmp,                          ONLY : namelist_filename,nlog,ncrt

      USE common_constants,                  ONLY : izero,m2cm

      USE nfreya_run,                        ONLY : codeid

      USE tport_mhd_grid_data,               ONLY : rmhdgrid, zmhdgrid,ymagn1,   &
                                                    rcontour,zcontour,ncontour

      USE zonal_data,                        ONLY : mfm1,nw,nh

      USE grid_class,                        ONLY : nj

      USE nfreya_version              

      USE Nfreya_namelist

      USE Vector_class,                      ONLY : vector_min_value,vector_max_value

      USE shot_info,                         ONLY : shot_id  

      IMPLICIT NONE

      REAL(DP) xdim,ydim,rmajor,rminor,kappa,rin,rout
      INTEGER(I4B) i,task,error,ix,iy
      INTEGER(I2B) nml_unit 

      error = izero


      IF(myid .ne. master)RETURN           

          IF(task == izero)THEN              ! echo input namelist onto nfreya_plot_file_name

             nfreya_plot_unit = get_next_io_unit ()
           !  nml_unit         = get_next_io_unit ()
             OPEN (unit = nfreya_plot_unit, file = nfreya_plot_file_name, status = 'UNKNOWN')
           !  OPEN (unit = nml_unit, file = namelist_filename, status = 'OLD',err = 1)
           !  CALL copy_file(nml_unit,nfreya_plot_unit) ! simply copy namelist to plot file
           !  WRITE(nfreya_plot_unit,FMT='("stop_nml")')
           !  CLOSE(UNIT  = nml_unit)
           !  process namelist properly here with redifined values from nubeam included

                WRITE(nfreya_plot_unit,10)
10              FORMAT("dump of nubeam namelist modified nfreya input namelist",/, &
                "-----------------------------------------------")
                ! Note: ---- above is scanned for in p_Nf_nubplt
                WRITE(nfreya_plot_unit,NML = run_data)
                WRITE(nfreya_plot_unit,FMT='("stop_nml")')

             RETURN
1                error = 1_I4B
             RETURN


          ELSEIF(task == 1)THEN
             xdim = rmhdgrid(nw) - rmhdgrid(1) 
             ydim = zmhdgrid(nw) - zmhdgrid(1)
             rmajor = dischg%rmajor*m2cm ; rminor = (dischg%rplasmin-dischg%zplasmin)*0.5_DP*m2cm
             kappa  = dischg%kappa
             rin    = Vector_min_value(dischg%rlimiter)  ;  rout = Vector_max_value(dischg%rlimiter)
             rin    = rin*m2cm                           ;  rout = rout*m2cm

             ! rin and rout are inside and outside major radius of vaccum vessel, used for
             ! plotting only.
             WRITE (nfreya_plot_unit, 8001)  ADJUSTL(codeid(1:LEN_TRIM(codeid))),p_nfreya_ver

             WRITE (nfreya_plot_unit, 8010)  no_physical_injectors,mfm1,nj,nw,nh,iborb
             WRITE (nfreya_plot_unit, 8010)  npitch,shot_id%shot_nmbr
             WRITE (nfreya_plot_unit, 8020)  xdim,ydim,ymagn1       !  in m ??
             WRITE (nfreya_plot_unit, 8020) (rmhdgrid(ix),ix=1,nw)  !  in cm 
             WRITE (nfreya_plot_unit, 8020) (zmhdgrid(iy),iy=1,nh)  !  in cm
             WRITE (nfreya_plot_unit, 8020)  rmajor,rminor,rout,rin,atw_beam,kappa
             WRITE (nfreya_plot_unit, 8020)  mhd_dat%tot_cur,mhd_dat%btor
             RETURN
8001         format(a,2x,'version =',a)
8010         format (6i10)
8011         format (2x, 1e12.4, 2x, l5, 2x, i5, 2x, i5, 2x, i5)
8012         format (2x, 1e12.4, 2x, i5)
8020         format (6(1pe14.6))  !this format must match other formats for ntrplt


          ELSEIF(task == 2)THEN 
             !WRITE (nfreya_plot_unit, '(i6)')izero
             WRITE (nfreya_plot_unit, '(" the end**")')


          ELSEIF(task == 3)THEN
              call NF_nubplt

          ELSEIF(task == -1)THEN
              CLOSE (unit = nfreya_plot_unit)

          ENDIF

          RETURN



      END SUBROUTINE      setup_plot_file

      SUBROUTINE sort_injectors
!----------------------------------------------------------------------
! -- create  a sorted list of physical injectors based on various injector
! -- attributes.
! -- This list is used only if number of processors assigned exceeds the
! -- number of physical beams lines the array sorted_physical_injectors
! -- is zero based so the MOD function can be applied directly, see sub
! -- reset_no_injectors
!----------------------------------------------------------------------
      USE nrtype,                            ONLY : DP,I4B,I2B
  
      USE neutral_beams,                     ONLY : sorted_physical_injectors, &
                                                    no_physical_injectors

      IMPLICIT NONE

      INTEGER(I4b) j

      IF(.NOT. ALLOCATED(sorted_physical_injectors))      &
           ALLOCATE(sorted_physical_injectors(0 : no_physical_injectors -1))

      ! for now just use input sequence
      DO j=1,no_physical_injectors
        sorted_physical_injectors(j-1) = j
      ENDDO

      END SUBROUTINE sort_injectors 



      SUBROUTINE copy_file(ioin,ioout)
! ----------------------------------------------------------------------
! -- copy file on unit  ioin to file on unit ioout
! -- files must exist and be open no error checking
! -- NOTE: P_NF_nubplt reads the P_NF_bpltfil file and checks certain items
! -- which means we need to maintain the ADJUSTL(line(1:LEN_TRIM(line)))
! -- below so as nto avoid parsing problems when the file is read
! ----------------------------------------------------------------------

      USE nrtype,                            ONLY : DP,I4B,I2B

      USE neutral_beams,                     ONLY : maxchr

        IMPLICIT NONE
        INTEGER(I2B) ioin,ioout
        CHARACTER(LEN = maxchr) line

        DO WHILE(.TRUE.)
           READ(ioin,'(a)',END = 10)line 
           WRITE(ioout,'(a)')ADJUSTL(line(1:LEN_TRIM(line)))
        ENDDO

10      RETURN

      END SUBROUTINE copy_file



      SUBROUTINE reset_no_injectors
!-----------------------------------------------------------------------------
! -- if we have more cpus than injectors then provided split_injectors is true 
! -- we will use the  additional cpus to do copies of the existing injectors
! -- this is set up here. 
!-----------------------------------------------------------------------------
   USE nrtype,                                      ONLY : DP,I4B

   USE nf_param,                                    ONLY : kb

   USE neutral_beams,                               ONLY : no_injectors,split_injectors,    &
                                                           no_physical_injectors,           &
                                                           pseudo_injectors,                &
                                                           sorted_physical_injectors

   USE MPI_data,                                    ONLY : master,myid,numprocs

   USE common_constants,                            ONLY : izero,zeroc

   IMPLICIT NONE

   INTEGER(I4B)      excess_proc, duplicate_injector_no,orig_injector_number, &
                     extra_injctr,k


         excess_proc = numprocs - no_injectors ! more processes than injectors?
         no_physical_injectors = no_injectors 


         IF(excess_proc .GT. izero .AND. split_injectors)THEN
            ALLOCATE(pseudo_injectors(excess_proc))
             pseudo_injectors(1:excess_proc) = -1 
             extra_injctr = 0
             !create pseudo injectors so number injectors = number processes
             k = -1 
             DO WHILE(no_injectors .LT. numprocs .AND. no_injectors .LT. kb)
                no_injectors  = no_injectors + 1
                extra_injctr  = extra_injctr + 1
                k = k+1
                duplicate_injector_no  = extra_injctr + no_physical_injectors
                orig_injector_number = MOD(k,no_physical_injectors)
                orig_injector_number  = sorted_physical_injectors(orig_injector_number)
                pseudo_injectors(duplicate_injector_no - no_physical_injectors) = orig_injector_number
                CALL copy_injector_attributes(orig_injector_number,duplicate_injector_no)

             ENDDO

             
 
         ENDIF

      END SUBROUTINE reset_no_injectors




      SUBROUTINE set_beam_id
!--------------------------------------------------------------------------------
! -- set nlco and beam_id for graphics output for cases that dont use the ufile
! -- type input.
! -- Note the convention:
! --      Positive current (from statefile) means that plasma current is
! --      in counter clockwise direction when viewed  from above the torus.
! --      For a standard rhs system (R,phi,Z) (z vertivcally upward)
! --      this means that psi is minimum 
! --      on axis (grad psi points outward from amg exis) .
! --      In Nfreya co injection (injection in direction of current)
! --      is set by angleh which is the beam centerline angle relative 
! --      to a vertical plane through the magnetic axis. Positive angleh means that the
! --      beam centerlines projection onto the phi direction is positive.
! --      Thus we have (Jphi is toroidal current density)
! --           Jphi            angleh           injections
! --          positive       positive            co
! --          positive       negative            ctr
! --          negative       positive            ctr
! --          negative       negative            co
! --      In other words if the sign of Jphi and angleh match
! --      then we have co injection. If the signs are opposite
! --      we have ctr injection.
! --      In the calcultiosn below negative tangency radii imply that
! --      the beam source is aimed in such away that its projection
! --      on the postive phi direction is negative. (eg the source
! --      center beamline points in the counter clockwise dircetion
! --      when the torus is viewd from the top.
!------------------------------------------------------------------HSJ--------------
   USE nrtype,                                      ONLY : DP,I4B

   USE neutral_beams,                               ONLY : nbeams,angleh,        &
                                                           rpivot,blenp,bhofset, &
                                                           rtan_ls,rtan_rs

   USE P_nfreya_interface,                          ONLY : beam_data

   USE Plasma_properties ,                          ONLY : mhd_dat

   USE common_constants,                            ONLY : izero,zeroc,Rad_Per_Deg

   IMPLICIT NONE

   INTEGER(I4B) jb
   REAL(DP) ah,x,y 
   CHARACTER*7 outform_ls,outform_rs


  RETURN  ! no longer used 8/15/11


      IF(ASSOCIATED(beam_data%nlco))DEALLOCATE(beam_data%nlco)
      ALLOCATE(beam_data%nlco(nbeams))
      IF(ASSOCIATED(beam_data%beam_id))DEALLOCATE(beam_data%beam_id)
      ALLOCATE(beam_data%beam_id(nbeams))
      IF(ALLOCATED(rtan_ls))DEALLOCATE(rtan_ls)
      ALLOCATE(rtan_ls(nbeams))
      IF(ALLOCATED(rtan_rs))DEALLOCATE(rtan_rs)
      ALLOCATE(rtan_rs(nbeams))




      DO jb = 1,nbeams
         go to 25 ! use nubeam naelist values for nlco
         IF(mhd_dat%tot_cur .GT. zeroc)THEN
            IF(angleh(jb) .GT. zeroc)THEN
               beam_data%nlco(jb) = .TRUE.
            ELSE
               beam_data%nlco(jb) = .FALSE.
            ENDIF

         ELSE
            IF(angleh(jb) .GT. zeroc)THEN
               beam_data%nlco(jb) = .FALSE.
            ELSE
               beam_data%nlco(jb) = .TRUE.
            ENDIF
         ENDIF


25       CONTINUE


         ! get tangency radius of left and right sources:
         ah = angleh(jb)*Rad_Per_Deg
         y = blenp(jb)*SIN(ah) + bhofset(jb)
         x = blenp(jb)*COS(ah)

         rtan_ls(jb) = rpivot(jb)* SIN(ATAN2(y,x))
         y = y  - 2._DP*bhofset(jb)

         rtan_rs(jb) = rpivot(jb)* SIN(ATAN2(y,x))
         Write(outform_ls,FMT ='(f7.2)')rtan_ls(jb)
         Write(outform_rs,FMT ='(f7.2)')rtan_rs(jb)

         beam_data%beam_id(jb) = 'rtan ls '//outform_ls//",rs "//outform_rs
      ENDDO

      RETURN

      END SUBROUTINE set_beam_id
