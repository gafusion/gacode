     SUBROUTINE distribute_statefile
!----------------------------------------------------------------
! -- this data should be packed into a single array and then 
! -- bcast. However only one call per run is made to this routine
! -- so it is not worthwhile to make the change. 
! ------------------------------------------------------HSJ-------

#if defined USEMPI
      USE mpi
#endif


      USE bc_values_gcnmp,           ONLY : totcur_bc,time_bc,                   &
                                              fix_edge_te_bc,fix_edge_ti_bc,     &
                                              fix_edge_rot_bc,fix_edge_ni_bc,    &
                                              zeff_bc,ene_bc,te_bc,ti_bc,        &
                                              angrot_bc,en_bc,vloop_bc,flux_bc


      USE curden_terms,              ONLY : curden,curohm,currf,curbeam,         &
                                            curboot,ibcur,irfc,curpar

      USE dep_var,                   ONLY : dp4                                 
                                         

      USE fast_ion_data_gcnmp,       ONLY : walp,enalp,w_alpha

      USE neutral_beams,             ONLY : nbion,nameb,enbeam,storqueb,wbeam,   &
                                            enbeam_tot,fd_beam,bptor

      USE fdyimp,                    ONLY : dfdt,dgdt,dhdt
    
      USE grid_class,                ONLY : nj,psir_grid,rho_grid,rho_gridn,     &
                                            eps,xhm2,xi11,xi33,xips,xhm20,xi110, &
                                            xi330,xips0,rcap,r2capi,r2cap

      USE iterdbmd_gcnmp,            ONLY : iterdsc,irwflag,iterdb_file_name,    &
                                            iterdb_outpt

      USE ions_gcnmp,                ONLY : nion,nprim,nprimp1,nneu,nimp,namep,  &
                                            namei,namen,zeff,zsq,z,              &
                                            name_size,fi_index,fd_thermal,atw,   &
                                            atomno

      USE neutral_data,              ONLY : enn,ennw,volsn,ennv


      USE nrtype,                    ONLY : DP,I4B,I2B

      USE Plasma_properties,         ONLY : dischg,profile,mhd_dat,get_values,   &
                                            diffuse,pwrden, prtcl_src,           &
                                            pellet,wpdot,fus_prod,neut_beam,     &
                                            plasma_frequencies

      USE tglfin,                    ONLY : tglf_p_output,tglf_e_output,         &
                                            tglf_m_output,                       &
                                            sign_It_tg_state,sign_Bt_tg_state

      USE glf23_gcnmp,               ONLY : glf_etg_output,glf_p_output,         &
                                            glf_e_output,glf_m_output,           &
                                            glf_eff_chi_output



      USE shot_info,                 ONLY : shot_id


      USE solcon_gcnmp,              ONLY : create_plot_file,time,tGCNMf,eqtime, &
                                            tGCNMs

      USE MPI_data,                  ONLY : mpi_start_time,mpi_end_time,         & 
                                            numprocs,master,myid,mpiierr,        &
                                            jac_comm, proc_map ,proc_time,       &
                                            comm_time

      USE source_terms_gcnmp,        ONLY : stsource,scx,sion,srecom,sbcx,sbeame, &
                                            dudtsv,sbeam,spellet

      USE vector_class,              ONLY : new_Vector,get_element,vector,       &
                                            real_mult_Vector,list,size_Vector,   &
                                            load_vector,length_Vector,           &
                                            delete_Vector_nf



      USE common_constants,          ONLY : izero,zeroc

    IMPLICIT NONE 
    INTEGER(I4B) j,jj,njnprim,ntot,njnneu,npsi,njnbion,njntot,njnion,  &
                 nwh,nsend,ke,ngW20,oknf
    REAL(DP), ALLOCATABLE,DIMENSION(:) :: work,work_npsi
    CHARACTER(len = name_size) dumy
    INTERFACE
       SUBROUTINE  bcast_vector(vector_in,nj,workt)
            USE nrtype,          ONLY : DP,I4B
            USE vector_class       
            TYPE(vector),INTENT(INOUT) :: vector_in
            INTEGER(I4B) nj
            REAL(DP),DIMENSION(:)      :: workt
       END  SUBROUTINE  bcast_vector
    END INTERFACE 


    
    jac_comm = 0
    IF(.NOT. ALLOCATED(proc_map))ALLOCATE(proc_map(0:numprocs-1))
    IF(.NOT. ALLOCATED(proc_time))ALLOCATE(proc_time(0:numprocs-1))
    IF(.NOT. ALLOCATED(comm_time))ALLOCATE(comm_time(0:numprocs-1))


    proc_map(:) = izero  ; proc_time(:) = zeroc ;comm_time(:) = zeroc
    eqtime = time ; tGCNMs = time 
    ibcur = 1_I4B ; irfc = 1_I4B
    IF(ABS(mhd_dat%totbeam_cur) .LT. 1.e-5)ibcur = 0
    IF(ABS(mhd_dat%totrf_cur)   .LT. 1.e-5)irfc  = 0


    IF(numprocs == 1) THEN
       !do this for 1 process here because nbion isnt known
       !on other processes until we pass it out below:
       IF(.NOT. ASSOCIATED(fi_index))ALLOCATE (fi_index(nbion))
       RETURN
    ENDIF

!  -----------------------------------------------------------------------------
!  following only gets executed if numprocs > 1 (otherwise returns above)
!  -----------------------------------------------------------------------------

#if defined USEMPI
    !at some point combine these scalars into a single bcast:

    CALL MPI_BCAST(iterdsc,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(shot_id%shot_nmbr,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(nj,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(fd_beam,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(fd_thermal,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    IF(myid == master) npsi  = mhd_dat%npsi
    CALL MPI_BCAST(npsi,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    IF(myid .NE. master) mhd_dat%npsi =npsi
    CALL MPI_BCAST(mhd_dat%mhd_info_avail,1,MPI_LOGICAL,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(nion,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(nprim,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(nprimp1,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(nimp,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(nneu,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(nbion,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(neut_beam%nbeams,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    IF(myid == master)ke = SIZE(neut_beam%fbcur,1)
    CALL MPI_BCAST(ke,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)

       IF( .NOT. ASSOCIATED(fi_index))ALLOCATE (fi_index(nbion))
       IF( .NOT. ASSOCIATED(namep))ALLOCATE ( namep(nprim))
       IF( .NOT. ASSOCIATED(namen))ALLOCATE ( namen(nneu) )
       IF( .NOT. ASSOCIATED(namei))ALLOCATE ( namei(nimp) )
       IF( .NOT. ASSOCIATED(nameb))ALLOCATE ( nameb(nbion) )
       IF( .NOT. ASSOCIATED(profile%en))ALLOCATE (profile%en(nion)) !2d array
       IF( .NOT. ASSOCIATED(wpdot%dpidt))ALLOCATE(wpdot%dpidt(nion))
       IF( .NOT. ASSOCIATED(prtcl_src%srecom))ALLOCATE(prtcl_src%srecom(nion))
       IF( .NOT. ASSOCIATED(pwrden%brems_nions))ALLOCATE(pwrden%brems_nions(nion))
       IF( .NOT. ASSOCIATED(diffuse%vpinch_nclass))ALLOCATE(diffuse%vpinch_nclass(nion+1))


       DO j=1,nimp 
          dumy(:) = namei(j)(:) !meaningful on process 0
          CALL MPI_BCAST(dumy,name_size,MPI_CHARACTER,master,MPI_COMM_WORLD,mpiierr)
          namei(j) = dumy(:) !now on all processes
       ENDDO

       DO j=1,nprim
          dumy(:) = namep(j)(:) !meaningful on process 0
          CALL MPI_BCAST(dumy,name_size,MPI_CHARACTER,master,MPI_COMM_WORLD,mpiierr)
          namep(j) = dumy(:) !now on all processes
       ENDDO

       DO j=1,nneu
          dumy(:) = namen(j)(:) !meaningful on process 0
          CALL MPI_BCAST(dumy,name_size,MPI_CHARACTER,master,MPI_COMM_WORLD,mpiierr)
          namen(j) = dumy(:) !now on all processes
       ENDDO

       DO j=1,nbion
          dumy(:) = nameb(j)(:) !meaningful on process 0
          CALL MPI_BCAST(dumy,name_size,MPI_CHARACTER,master,MPI_COMM_WORLD,mpiierr)
          nameb(j) = dumy(:) !now on all processes
       ENDDO

        IF(myid == master) dumy(:) = pellet%name(:)
        CALL MPI_BCAST(dumy,name_size,MPI_CHARACTER,master,MPI_COMM_WORLD,mpiierr)
        IF(myid .NE. master) pellet%name = dumy(:) !now on all processes


    CALL MPI_BCAST(fi_index,nbion,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(nameb,8,MPI_CHARACTER,master,MPI_COMM_WORLD,mpiierr)

    CALL MPI_BCAST(time,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    shot_id%shot_time = time 

    CALL MPI_BCAST(dischg%nr_mhd,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%nz_mhd,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)

    CALL MPI_BCAST(dischg%rplasmin,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%rplasmax,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    CALL MPI_BCAST(dischg%zplasmin,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%zplasmax,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%rgeom,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   CALL MPI_BCAST(dischg%btgeom,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%rmag,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   CALL MPI_BCAST(dischg%rma,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   CALL MPI_BCAST(dischg%zma,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%rmajor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%kappa,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%deltao,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%pindento,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%volo,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%areao,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%rminor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%rsep,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%zsep,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dischg%circum,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    nwh = dischg%nr_mhd * dischg%nz_mhd



    IF(.NOT. ASSOCIATED(mhd_dat%psi))ALLOCATE(mhd_dat%PSI(dischg%nr_mhd,dischg%nz_mhd))
    CALL MPI_BCAST(mhd_dat%psi,nwh,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%btor,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%psibdry,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%psiaxis,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%tot_cur,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%totohm_cur,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%totboot_cur,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%totbeam_cur,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%totrf_cur,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST( mhd_dat%betap,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%beta,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%ali,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mhd_dat%R0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(profile%te0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(profile%ti0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)



   IF( .NOT. ALLOCATED(work))ALLOCATE(work(nj))
   IF( .NOT. ALLOCATED(work_npsi))ALLOCATE(work_npsi(npsi))

! MPI 1 at least does not work with passing derived types that
! contain pointers. Using a subroutine  to do the following sequence
! of steps also failed. BCAST must be called by all processes according 
! to  the MPI 1 standard. Hence the awkward but workable:


   IF(myid .EQ. master)work(:) = get_values(pwrden%qexch)
!  nj is known on all processes (see above)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qexch = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(psir_grid)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)psir_grid = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(rho_grid)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)rho_grid = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(rho_gridn)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)rho_gridn = new_Vector(nj,work)
    

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%fcap)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%fcap = new_Vector(nj,work)
    

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%hcap)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%hcap = new_Vector(nj,work)
    
   IF(myid .EQ. master)work(:) = get_values(mhd_dat%gcap)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%gcap = new_Vector(nj,work)

  IF(myid .EQ. master)work(:) = get_values(mhd_dat%btotrmaj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%btotrmaj = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%curden)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%curden = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%curpar)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%curpar  = new_Vector(nj,work)

   IF(myid .EQ. master)work_npsi(:) = get_values(mhd_dat%qpsinpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%qpsinpsi  = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(mhd_dat%pressnpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%pressnpsi  = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(mhd_dat%ffprimnpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%ffprimnpsi  = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(mhd_dat%pprimnpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%pprimnpsi  = new_Vector(npsi,work_npsi)


   IF(myid .EQ. master)work(:) = get_values(mhd_dat%curohm)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%curohm = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%curboot)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%curboot = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%q_value)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%q_value = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%rbp)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%rbp = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%bprmaj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%bprmaj = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%ffprim)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%ffprim = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%pprim)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%pprim = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%rbp)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%rbp = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%r2cap)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%r2cap = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%r2capi)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%r2capi = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%rcap)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%rcap = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%rcapi)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%rcapi = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(mhd_dat%betan)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%betan = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(wpdot%dnedt)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)wpdot%dnedt = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(wpdot%dpedt)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)wpdot%dpedt = new_Vector(nj,work)


   DO j = 1,nion
      IF(myid == master)work(:) =  wpdot%dpidt(j)%data(:)
      CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      IF(myid .NE. master)wpdot%dpidt(j) = new_Vector(nj,work)
      IF(myid == master)work(:) =  prtcl_src%srecom(j)%data(:)
      CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      IF(myid .NE. master)prtcl_src%srecom(j) = new_Vector(nj,work)
   ENDDO


   IF(myid .EQ. master)work(:) = get_values(profile%te)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%te = new_Vector(nj,work)
   
   IF(myid .EQ. master)work(:) = get_values(profile%ti)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%ti = new_Vector(nj,work)
    
   IF(myid .EQ. master)work(:) = get_values(profile%ene)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%ene = new_Vector(nj,work)
    
   IF(myid .EQ. master)work(:) = get_values(profile%zeff)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%zeff = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(profile%vpol)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%vpol = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(profile%vpol_nclass)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%vpol_nclass = new_Vector(nj,work)


   IF(myid .EQ. master)work(:) = get_values(profile%vpar)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%vpar = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(profile%vpar_nclass)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%vpar_nclass = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(profile%er_tot_nclass)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%er_tot_nclass = new_Vector(nj,work)


   IF(myid .EQ. master)work(:) = get_values(profile%press)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%press = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(profile%etor)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%etor = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(profile%angrot)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%angrot = new_Vector(nj,work)

   DO j=1,nion
       IF(myid .EQ. master)work(:) = get_values(profile%en(j))
       CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
       IF(myid .NE. master)profile%en(j) = new_Vector(nj,work)
   END DO

   ntot = nion+dp4
   IF( .NOT. ASSOCIATED(profile%flux))ALLOCATE (profile%flux(ntot)) !2d array
   IF( .NOT. ASSOCIATED(profile%flux_conv))ALLOCATE (profile%flux_conv(ntot)) !2d array
   DO j=1,ntot
       IF(myid .EQ. master)work(:) = get_values(profile%flux(j))
       CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
       IF(myid .NE. master)profile%flux(j) = new_Vector(nj,work)
       IF(myid .EQ. master) work(:) = get_values(profile%flux_conv(j))
       CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
       IF(myid .NE. master)profile%flux_conv(j) = new_Vector(nj,work)
   END DO
   IF(myid .EQ. master)work(:) = get_values(profile%fluxe)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)profile%fluxe = new_Vector(nj,work)




    IF(.NOT. ALLOCATED(flux_bc))ALLOCATE(flux_bc(nj,ntot))
  

    IF( .NOT. ALLOCATED(stsource))THEN
!              ALLOCATE(stsource(nion,nj),sion(nj,nion),srecom(nj,nion),       &
!                       sbcx(nj,nion),scx(nj,nion),dudtsv(ntot,nj))
              ALLOCATE(stsource(nion,nj),sion(nj,nion),                        &
                       sbcx(nj,nion),scx(nj,nion),dudtsv(ntot,nj))
    ENDIF

    njnprim = nj*nprim
    njntot  = nj*ntot
    njnion  = nj*nion

    CALL MPI_BCAST(sion,njnion,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
!    CALL MPI_BCAST(srecom,njnion,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(scx,njnion,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(sbcx,njnion,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(stsource,njnion,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dudtsv,njntot,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(flux_bc,njntot,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF( .NOT. ALLOCATED(enbeam)) ALLOCATE(enbeam(nj,nbion))
    njnbion = nj * nbion
    CALL MPI_BCAST(enbeam,njnbion,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    IF( .NOT. ALLOCATED(enbeam_tot))THEN 
       ALLOCATE(enbeam_tot(nj))
       enbeam_tot(:) = zeroc
    ENDIF
    CALL MPI_BCAST(enbeam_tot,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

  ! bptor ==> beam power into torus 
    CALL MPI_BCAST(bptor,neut_beam%nbeams,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)



    IF(.NOT. ASSOCIATED(neut_beam%fbcur))ALLOCATE(neut_beam%fbcur(ke,neut_beam%nbeams))
    CALL MPI_BCAST(neut_beam%fbcur,ke*neut_beam%nbeams,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF(.NOT. ASSOCIATED(neut_beam%prompt_pwr_in_plasma))    &
         ALLOCATE(neut_beam%prompt_pwr_in_plasma(ke,neut_beam%nbeams))
    CALL MPI_BCAST(neut_beam%prompt_pwr_in_plasma,ke*neut_beam%nbeams,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
 
    IF(.NOT. ASSOCIATED(neut_beam%pbeam))    &
         ALLOCATE(neut_beam%pbeam(ke,neut_beam%nbeams))
    CALL MPI_BCAST(neut_beam%pbeam,ke*neut_beam%nbeams,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF(.NOT. ASSOCIATED(neut_beam%ebeam))    &
         ALLOCATE(neut_beam%ebeam(ke,neut_beam%nbeams))
    CALL MPI_BCAST(neut_beam%ebeam,ke*neut_beam%nbeams,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF(.NOT. ALLOCATED(enn))ALLOCATE(enn(nj,nneu))
    njnneu = nj*nneu
    CALL MPI_BCAST(enn,njnneu,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF(.NOT. ALLOCATED(ennw))ALLOCATE(ennw(nj,nneu))
    CALL MPI_BCAST(ennw,njnneu,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF(.NOT. ALLOCATED(ennv))ALLOCATE(ennv(nj,nneu))
    CALL MPI_BCAST(ennv,njnneu,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF(.NOT. ALLOCATED(volsn))ALLOCATE(volsn(nj,nneu))
    CALL MPI_BCAST(volsn,njnneu,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF(.NOT. ALLOCATED(sbeame))ALLOCATE(sbeame(nj))
    CALL MPI_BCAST(sbeame,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF(nbion >  0)THEN
       IF(.NOT. ALLOCATED(sbeam))ALLOCATE(sbeam(nj,nbion))
       njnbion = nj*nbion
       CALL MPI_BCAST(sbeam,njnbion,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    ENDIF
        
    IF(.NOT. ALLOCATED(curbeam))ALLOCATE(curbeam(nj))
    CALL MPI_BCAST(curbeam,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

        
    IF(.NOT. ALLOCATED(currf))ALLOCATE(currf(nj))
    CALL MPI_BCAST(currf,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

!    IF(.NOT. ALLOCATED(rbp))ALLOCATE(rbp(nj))
    IF(.NOT. ASSOCIATED(zeff))ALLOCATE(zeff(nj))
    IF(myid == master) work(:) = zeff(:)
    CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    zeff(:) = work(:)

   IF(myid .EQ. master)work(:) = get_values(diffuse%chieinv)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)diffuse%chieinv = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(diffuse%chiinv)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)diffuse%chiinv = new_Vector(nj,work)

 

   IF(myid .EQ. master)work(:) = get_values(diffuse%xkineo)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)diffuse%xkineo = new_Vector(nj,work)


   IF(myid .EQ. master)work(:) = get_values(diffuse%eta)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)diffuse%eta = new_Vector(nj,work)


   IF(myid .EQ. master)work(:) = get_values(diffuse%ftrap)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)diffuse%ftrap = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(diffuse%xnuse)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)diffuse%xnuse = new_Vector(nj,work)

   IF(myid .NE. master) ALLOCATE(diffuse%xnus(nion))
   DO j=1,nion
      IF(myid .EQ. master)work(:) = get_values(diffuse%xnus(j))
      CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      IF(myid .NE. master)diffuse%xnus(j) = new_Vector(nj,work)
   ENDDO
 
    !---------------------------------------------------------------------------
    ! plasma frequencies
    ! master had these allocated in statefile read.
    !---------------------------------------------------------------------------
    IF(.NOT. ASSOCIATED(plasma_frequencies%omega_pi)) &
          ALLOCATE(plasma_frequencies%omega_pi(nprimp1))
    IF(.NOT. ASSOCIATED(plasma_frequencies%omega_ci)) &
          ALLOCATE(plasma_frequencies%omega_ci(nprimp1))
    IF(.NOT. ASSOCIATED(plasma_frequencies%omega_lh)) &
          ALLOCATE(plasma_frequencies%omega_lh(nprimp1))
    IF(.NOT. ASSOCIATED(plasma_frequencies%omega_uh)) &
          ALLOCATE(plasma_frequencies%omega_uh(nprimp1))
   DO j=1,nprimp1
      IF(myid .EQ. master)work(:) = get_values(plasma_frequencies%omega_pi(j))
      CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      IF(myid .NE. master)plasma_frequencies%omega_pi(j) = new_Vector(nj,work)
      IF(myid .EQ. master)work(:) = get_values(plasma_frequencies%omega_ci(j))
      CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      IF(myid .NE. master)plasma_frequencies%omega_ci(j) = new_Vector(nj,work)
      IF(myid .EQ. master)work(:) = get_values(plasma_frequencies%omega_lh(j))
      CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      IF(myid .NE. master)plasma_frequencies%omega_lh(j) = new_Vector(nj,work)
      IF(myid .EQ. master)work(:) = get_values(plasma_frequencies%omega_uh(j))
      CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      IF(myid .NE. master)plasma_frequencies%omega_uh(j) = new_Vector(nj,work)
   ENDDO

      IF(myid == master)work(:) = get_values(plasma_frequencies%omega_ce)
      CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      IF(myid .NE. master)plasma_frequencies%omega_ce = new_Vector(nj,work)
      IF(myid == master)work(:) = get_values(plasma_frequencies%omega_pe)
      CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      IF(myid .NE. master)plasma_frequencies%omega_pe = new_Vector(nj,work)

 !------------------------------------------------------------------------------
! NOTE: glf,tglf,mmm *_output arrays are read in and bcast becasue we may
! use some of the information in the current run to construc fixed chis,etc.
! -- glf23 related
!-------------------------------------------------------------------------------
     IF(myid == master) THEN
       IF(.NOT. ALLOCATED(glf_eff_chi_output))THEN
          ALLOCATE(glf_eff_chi_output(ntot,nj))
          glf_eff_chi_output(:,:) = zeroc
       ENDIF
     ENDIF
     IF(myid .NE. master)THEN
        IF(ALLOCATED(glf_p_output))DEALLOCATE(glf_p_output)
        ALLOCATE(glf_p_output(nj,3))
        IF(ALLOCATED(glf_e_output))DEALLOCATE(glf_e_output)
        ALLOCATE(glf_e_output(nj,3))
        IF(ALLOCATED(glf_etg_output))DEALLOCATE(glf_etg_output)
        ALLOCATE(glf_etg_output(nj))
        IF(ALLOCATED(glf_eff_chi_output))DEALLOCATE(glf_eff_chi_output)
        ALLOCATE(glf_eff_chi_output(ntot,nj))
     ENDIF


      DO j=1,3
         IF(myid == master)work(:) = glf_p_output(:,j)
         CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
         IF(myid .NE. master)glf_p_output(:,j) = work(:)

         IF(myid == master)work(:) = glf_e_output(:,j)
         CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
         IF(myid .NE. master)glf_e_output(:,j) = work(:)

         IF(myid == master)work(:) = glf_p_output(:,j)
         CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      ENDDO
      IF(myid == master)work = glf_etg_output(:)
      CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      IF(myid .NE. master)glf_etg_output(:) = work(:)

      CALL MPI_BCAST(glf_eff_chi_output,njntot,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)


   !----------------------------------------------------------------------------
   ! -- multimode related
   !----------------------------------------------------------------------------

      ngW20 = 4
      IF(myid .NE. master)THEN
         IF(ASSOCIATED(diffuse%mmm_gammaW20))DEALLOCATE(diffuse%mmm_gammaW20)
         ALLOCATE(diffuse%mmm_gammaW20(ngW20))     
         IF(ASSOCIATED(diffuse%mmm_omegaW20))DEALLOCATE(diffuse%mmm_omegaW20)
         ALLOCATE(diffuse%mmm_omegaW20(ngW20))    
         IF(ASSOCIATED(diffuse%mmm_vflux))DEALLOCATE(diffuse%mmm_vflux)
         ALLOCATE(diffuse%mmm_vflux(ngW20))    
         IF(ASSOCIATED(diffuse%mmm_vconv))DEALLOCATE(diffuse%mmm_vconv)
         ALLOCATE(diffuse%mmm_vconv(ngW20+2))    
      ENDIF
       DO jj =1,ngW20+2
           IF(jj .LE. ngW20)THEN
              IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_gammaW20(jj))
              CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
              IF(myid .NE. master)diffuse%mmm_gammaW20(jj) = new_Vector(nj,work)

              IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_omegaW20(jj))
              CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
              IF(myid .NE. master)diffuse%mmm_omegaW20(jj) = new_Vector(nj,work)

              IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_vflux(jj))
              CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
              IF(myid .NE. master)diffuse%mmm_vflux(jj) = new_Vector(nj,work)
           ENDIF
           IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_vconv(jj))
           CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
           IF(myid .NE. master)diffuse%mmm_vconv(jj) = new_Vector(nj,work)
           IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_vconv(jj))
           CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
           IF(myid .NE. master)diffuse%mmm_vconv(jj) = new_Vector(nj,work)
        ENDDO

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_gammaDBM)

   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_gammaDBM)
      diffuse%mmm_gammaDBM = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_omegaDBM)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_omegaDBM)
      diffuse%mmm_omegaDBM = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xti)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xtiDBM)
      diffuse%mmm_xti = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xdi)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xdi)
      diffuse%mmm_xdi = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xte)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xte)
      diffuse%mmm_xte = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xdz)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xdz)
      diffuse%mmm_xdz = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xvt)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xvt)
      diffuse%mmm_xvt = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xvp)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xvp)
      diffuse%mmm_xvp = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xtiW20)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xtiW20)
      diffuse%mmm_xtiW20 = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xdiW20)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xdiW20)
      diffuse%mmm_xdiW20 = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xteW20)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xteW20)
      diffuse%mmm_xteW20 = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xtiDBM)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xtiDBM)
      diffuse%mmm_xtiDBM = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xdiDBM)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xdiDBM)
      diffuse%mmm_xdiDBM = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xteDBM)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xteDBM)
      diffuse%mmm_xteDBM = new_Vector(nj,work)
   ENDIF

   IF(myid .EQ. master)work(:) = get_values(diffuse%mmm_xteETG)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)THEN
      oknf = delete_Vector_nf(diffuse%mmm_xteETG)
      diffuse%mmm_xteETG = new_Vector(nj,work)
   ENDIF

   ! end multimode -------------------------------------------------------------

   ! tgl realted ---------------------------------------------------------------
      IF( myid .NE. master)THEN
         IF(ALLOCATED(tglf_p_output))DEALLOCATE(tglf_p_output)
         ALLOCATE(tglf_p_output(nj,3))
         IF(ALLOCATED(tglf_e_output))DEALLOCATE(tglf_e_output)
         ALLOCATE(tglf_e_output(nj,3))
         IF(ALLOCATED(tglf_m_output))DEALLOCATE(tglf_m_output)
         ALLOCATE(tglf_m_output(nj,3))
      ENDIF
      DO j=1,3
         IF(myid .EQ. master)work(:) = tglf_p_output(:,j)
         CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
         IF(myid .NE. master)tglf_p_output(:,j) = work(:)

         IF(myid .EQ. master)work(:) = tglf_e_output(:,j)
         CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
         IF(myid .NE. master)tglf_e_output(:,j) = work(:)

         IF(myid .EQ. master)work(:) = tglf_m_output(:,j)
         CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
         IF(myid .NE. master)tglf_m_output(:,j) = work(:)
      ENDDO
   ! end tgl realted -----------------------------------------------------------



   IF(myid .EQ. master)work(:) = get_values(fus_prod%neutr_ddn_th)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)fus_prod%neutr_ddn_th = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(fus_prod%neutr_ddn_beam_beam)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)fus_prod%neutr_ddn_beam_beam = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(fus_prod%neutr_ddn_beam_thermal)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)fus_prod%neutr_ddn_beam_thermal = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(fus_prod%neutr_ddn_knock)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)fus_prod%neutr_ddn_knock = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(fus_prod%neutr_ddn_tot)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)fus_prod%neutr_ddn_tot = new_Vector(nj,work)

  
   CALL MPI_BCAST(fus_prod%total_neutr_ddn_th,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   CALL MPI_BCAST(fus_prod%total_neutr_ddn_beam_beam,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   CALL MPI_BCAST(fus_prod%total_neutr_ddn_beam_thermal,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   CALL MPI_BCAST(fus_prod%total_neutr_ddn,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qconde)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qconde= new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qcondi)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qcondi = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qconve)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qconve = new_Vector(nj,work)


   IF(myid .EQ. master)work(:) = get_values(pwrden%qconvi)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qconvi = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qbeame)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qbeame = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qbeami)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qbeami = new_Vector(nj,work)


   IF(myid .EQ. master)work(:) = get_values(pwrden%qdelt)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qdelt = new_Vector(nj,work)


   IF(myid .EQ. master)work(:) = get_values(pwrden%qrfe)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qrfe = new_Vector(nj,work)


   IF(myid .EQ. master)work(:) = get_values(pwrden%qrfi)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qrfi = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qione)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qione = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qioni)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qioni = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qcx)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qcx = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qe2d)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qe2d = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qi2d)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qi2d = new_Vector(nj,work)


   IF(myid .EQ. master)work(:) = get_values(pwrden%qfuse)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qfuse = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qfusi)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qfusi = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qbfuse)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qbfuse = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qbfusi)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qbfusi = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qmag)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qmag = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qsawe)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qsawe = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qsawi)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qsawi = new_Vector(nj,work)
 
   IF(myid .EQ. master)work(:) = get_values(pwrden%qrad)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qrad = new_Vector(nj,work)
   DO j=1,nion
       IF(myid .EQ. master)work(:) = get_values(pwrden%brems_nions(j))
       CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
       IF(myid .NE. master)pwrden%brems_nions(j) = new_Vector(nj,work)
   END DO
   DO j=1,nion+1
       IF(myid .EQ. master)work(:) = get_values(diffuse%vpinch_nclass(j))
       CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
       IF(myid .NE. master)diffuse%vpinch_nclass(j) = new_Vector(nj,work)
   END DO

   IF(myid .EQ. master)work(:) = get_values(pwrden%qohm)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qohm = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%omegale)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%omegale = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qomegapi)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qomegapi = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%qangce)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%qangce = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%sprcxre)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%sprcxre = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%sprcxree)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%sprcxree = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(pwrden%spreimpe)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)pwrden%spreimpe = new_Vector(nj,work)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%rmajavnpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%rmajavnpsi = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%rminavnpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%rminavnpsi = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(mhd_dat%psivalnpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)mhd_dat%psivalnpsi = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%psivolpnpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%psivolpnpsi = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%elongxnpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%elongxnpsi = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%triangnpsi_u)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%triangnpsi_u = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%triangnpsi_l)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%triangnpsi_l = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%sfareanpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%sfareanpsi = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%cxareanpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%cxareanpsi = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%pindentnpsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%pindentnpsi = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%grho1npsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%grho1npsi = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work_npsi(:) = get_values(dischg%grho2npsi)
   CALL MPI_BCAST(work_npsi,npsi,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%grho2npsi = new_Vector(npsi,work_npsi)

   IF(myid .EQ. master)work(:) = get_values(prtcl_src%stfuse)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)prtcl_src%stfuse = new_Vector(nj ,work)

   IF(myid .EQ. master)work(:) = get_values(prtcl_src%sbfuse)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)prtcl_src%sbfuse = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(prtcl_src%spellet)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)prtcl_src%spellet = new_Vector(nj,work)

   DEALLOCATE(work,work_npsi)
   ALLOCATE(work(nj))
   IF(myid .EQ. master)work(:) = get_values(dischg%rmajavnj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%rmajavnj = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(dischg%rminavnj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%rminavnj = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(dischg%psivolpnj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%psivolpnj = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(dischg%elongxnj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%elongxnj = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(dischg%triangnj_u)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%triangnj_u = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(dischg%triangnj_l)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%triangnj_l = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(dischg%pindentnj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%pindentnj = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(dischg%sfareanj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%sfareanj = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(dischg%cxareanj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%cxareanj = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(dischg%grho1nj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%grho1nj = new_Vector(nj,work)

   IF(myid .EQ. master)work(:) = get_values(dischg%grho2nj)
   CALL MPI_BCAST(work,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%grho2nj = new_Vector(nj,work)



   CALL MPI_BCAST(dischg%nplasbdry,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)

   DEALLOCATE(work)
   ALLOCATE(work(dischg%nplasbdry))

   IF(myid .EQ. master)work(:) = get_values(dischg%rplasbdry)
   CALL MPI_BCAST(work,dischg%nplasbdry,MPI_DOUBLE_PRECISION,master, &
        MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%rplasbdry = new_Vector(dischg%nplasbdry,work)


   IF(myid .EQ. master)work(:) = get_values(dischg%zplasbdry)
   CALL MPI_BCAST(work,dischg%nplasbdry,MPI_DOUBLE_PRECISION,master, &
        MPI_COMM_WORLD,mpiierr)
   IF(myid .NE. master)dischg%zplasbdry = new_Vector(dischg%nplasbdry,work)


    IF(.NOT. ALLOCATED(storqueb))ALLOCATE(storqueb(nj))
    CALL MPI_BCAST(storqueb,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)


    IF(.NOT. ASSOCIATED(atw))ALLOCATE(atw(nion))  ! atw,atomno were set in set_ion_prop2
    IF(.NOT. ASSOCIATED(atomno))ALLOCATE(atomno(nion))
    CALL MPI_BCAST(atw,nion,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(atomno,nion,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)


      IF(.NOT. ALLOCATED(eps))ALLOCATE(eps(nj))
      IF(.NOT. ALLOCATED(te_bc))ALLOCATE(te_bc(nj))
      IF(.NOT. ALLOCATED(ti_bc))ALLOCATE(ti_bc(nj))
      IF(.NOT. ALLOCATED(ene_bc))ALLOCATE(ene_bc(nj))
      IF(.NOT. ALLOCATED(zeff_bc))ALLOCATE(zeff_bc(nj))
      IF(.NOT. ALLOCATED(angrot_bc))ALLOCATE(angrot_bc(nj))
      IF(.NOT. ALLOCATED(en_bc))ALLOCATE(en_bc(nj,nion))
      IF(.NOT. ALLOCATED(wbeam))ALLOCATE(wbeam(nj))
      IF(.NOT. ALLOCATED(walp))ALLOCATE(walp(nj))
      IF(.NOT. ALLOCATED(w_alpha))ALLOCATE(w_alpha(nj))
      IF(.NOT. ALLOCATED(enalp))ALLOCATE(enalp(nj))
      IF(.NOT. ASSOCIATED(zsq))ALLOCATE(zsq(nj,nion))
      IF(.NOT. ASSOCIATED(z))ALLOCATE(z(nj,nion))
      IF(.NOT. ALLOCATED(dfdt))ALLOCATE(dfdt(nj))
      IF(.NOT. ALLOCATED(dgdt))ALLOCATE(dgdt(nj))
      IF(.NOT. ALLOCATED(dhdt))ALLOCATE(dhdt(nj))
      IF(.NOT. ALLOCATED(xips))ALLOCATE(xips(nj))
      IF(.NOT. ALLOCATED(xips0))ALLOCATE(xips0(nj))
      IF(.NOT. ALLOCATED(xhm2))ALLOCATE(xhm2(nj))
      IF(.NOT. ALLOCATED(xhm20))ALLOCATE(xhm20(nj))
      IF(.NOT. ALLOCATED(xi11))ALLOCATE(xi11(nj))
      IF(.NOT. ALLOCATED(xi110))ALLOCATE(xi110(nj))
      IF(.NOT. ALLOCATED(xi33))ALLOCATE(xi33(nj))
      IF(.NOT. ALLOCATED(xi330))ALLOCATE(xi330(nj))
      IF(.NOT. ALLOCATED(fix_edge_ni_bc))ALLOCATE(fix_edge_ni_bc(nion))

      CALL MPI_BCAST(tGCNMf,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      CALL MPI_BCAST(time_bc,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      CALL MPI_BCAST(totcur_bc,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      CALL MPI_BCAST(vloop_bc,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      CALL MPI_BCAST(fix_edge_te_bc,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      CALL MPI_BCAST(fix_edge_ti_bc,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      CALL MPI_BCAST(fix_edge_rot_bc,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
      CALL MPI_BCAST(fix_edge_ni_bc,nion,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)


    CALL MPI_BCAST(dfdt,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dgdt,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(dhdt,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(te_bc,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ti_bc,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ene_bc,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(zeff_bc,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(angrot_bc,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    DO j=1,nion
       CALL MPI_BCAST(en_bc(1,j),nj,MPI_DOUBLE_PRECISION,master,&
                MPI_COMM_WORLD,mpiierr)
       CALL MPI_BCAST(z(1,j),nj,MPI_DOUBLE_PRECISION,master,&
                MPI_COMM_WORLD,mpiierr)

       CALL MPI_BCAST(zsq(1,j),nj,MPI_DOUBLE_PRECISION,master,&
                MPI_COMM_WORLD,mpiierr)
    END DO
    CALL MPI_BCAST(wbeam,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF(myid == master)walp(:) = get_values(profile%walp)
    CALL MPI_BCAST(walp,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    IF(myid .NE. master)profile%walp   = new_Vector(nj,walp)  

    CALL MPI_BCAST(enalp,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(eps,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    CALL MPI_BCAST(xhm2,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(xi11,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(xi33,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(xips,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    IF(myid .NE. master)THEN
       xhm20(:) = xhm2(:)
       xi110(:) = xi11(:)
       xi330(:) = xi33(:)
       xips0(:) = xips(:)
    ENDIF
   

    
    IF(SIZE(work) .LT. MAX(dischg%nr_mhd,dischg%nz_mhd))THEN
       DEALLOCATE(work)
       ALLOCATE(work( MAX(dischg%nr_mhd,dischg%nz_mhd)))
    ENDIF
    IF(myid .EQ. master)work(1:dischg%nr_mhd) = get_Values(dischg%rmhdgrid)
    CALL MPI_BCAST(work,dischg%nr_mhd,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    IF(myid .NE. MASTER)dischg%rmhdgrid = new_Vector(dischg%nz_mhd,work)
    IF(myid .EQ. master)work(1:dischg%nz_mhd) = get_Values(dischg%zmhdgrid)
    CALL MPI_BCAST(work,dischg%nz_mhd,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    IF(myid .NE. MASTER)dischg%zmhdgrid = new_Vector(dischg%nz_mhd,work)


    
    CALL MPI_BCAST(dischg%nlimiter,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    IF(SIZE(work) .LT. dischg%nlimiter)THEN
       DEALLOCATE(work)
       ALLOCATE(work(dischg%nlimiter))
    ENDIF
    IF(myid .EQ. master)work(1:dischg%nlimiter) = get_Values(dischg%rlimiter)
    CALL MPI_BCAST(work,dischg%nlimiter,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    IF(myid .NE. MASTER)dischg%rlimiter = new_Vector(dischg%nlimiter,work)
    IF(myid .EQ. master)work(1:dischg%nlimiter) = get_Values(dischg%zlimiter)
    CALL MPI_BCAST(work,dischg%nlimiter,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    IF(myid .NE. master)dischg%zlimiter = new_Vector(dischg%nlimiter,work)


    nsend = izero
    IF(myid == master)THEN
       nsend = length_Vector(mhd_dat%fpsinpsi)
    ENDIF
    CALL MPI_BCAST(nsend,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    IF(SIZE(work) .LT. nsend)THEN
       DEALLOCATE(work)
       ALLOCATE(work(nsend))
    ENDIF
    IF(myid .EQ. master)work(1:nsend) = get_Values(mhd_dat%fpsinpsi)
    CALL MPI_BCAST(work,nsend,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    IF(myid .NE. master)mhd_dat%fpsinpsi = new_Vector(nsend,work)
    DEALLOCATE(WORK)



    IF(myid == master)THEN
       sign_It_tg_state = SIGN(1._DP,mhd_dat%tot_cur)
       sign_Bt_tg_state = SIGN(1._DP,mhd_dat%btor)
    ENDIF
    CALL MPI_BCAST(sign_It_tg_state,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(sign_Bt_tg_state,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)


    RETURN
    END SUBROUTINE distribute_statefile
   

    SUBROUTINE  bcast_vector(vector_in,nj,workt)
! ---------------------------------------------------------------------------------
! the derived type vector uses pointers (because allocatable arrays are not
! allowed in derived types) But MPI (version 1) does not allow for typemaps
! that will properly broadcast this type of variable. (MPI_STRUCT does not
! work !!). This subroutine  is left in for my reference but is not called.
! ------------------------------------------------------------------HSJ--01/19/06---
   USE nrtype,                         ONLY : DP,I4B,I2B

   USE vector_class,                   ONLY : vector,new_Vector,get_values,      &
                                                 length_Vector,load_Vector

   USE MPI_data,                       ONLY : myid,master

   USE mpi

   IMPLICIT NONE
   TYPE(vector),INTENT(INOUT) :: vector_in
   INTEGER(I4B) nj,njj,mpiierr
   REAL(DP),DIMENSION(:) :: workt

   njj = length_Vector(vector_in)
   IF(nj .NE. njj) PRINT *,'myid njj prob',myid,njj,nj
   workt(:) = get_values(vector_in)
   CALL MPI_BCAST(nj,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
   CALL MPI_BCAST(workt,nj,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
#endif


   RETURN
   END
