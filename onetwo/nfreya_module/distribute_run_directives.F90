    SUBROUTINE distribute_run_directives(read_ufile,read_nubeam_namelist)
!  --------------------------------------------------------------------------
!  bring processes other than the master up to date with the namelist  input
!  (which was read only by the master process)
!  We could  pack the data using mpi type def but since these
!  messages (and similar ones in distribute_iterdb) are accessed
!  only onece at the beginning of a run it is probably not worthwhile.
!  -------------------------------------------------------------HSJ-----------

#if defined (USEMPI)
    USE nf_param,              ONLY : kb,ke,kimp,nap

    USE nrtype,                ONLY : DP,I4B,I2B

    USE Plasma_properties,     ONLY : mhd_dat,pwrden

 
    USE MPI_data,              ONLY : myid,master,mpiierr,numprocs



    USE P_nfreya_interface,    ONLY : d_fast_ion

    USE io_gcnmp,              ONLY : namelist_filename,nlog,ncrt

    USE ions_gcnmp,            ONLY : nprim,nimp,name_size,ni_sc,atw,atomno

    USE nub,                   ONLY : bfr_neutrlz

    USE MPI

    USE neutral_beams,         ONLY  : timbplt, beam_on, beam_time,           &
                                       nbeams,no_injectors,namelist_nameb,    &
                                       anglev, angleh,relnub,nashape,bcur,    &
                                       bptor,aheigh,awidth,blenp,bleni,bheigh,&
                                       bwidth,bhfoc,bvfoc,bhdiv,bvdiv,nbshape,&
                                       ebkev,fbcur,alen,bvofset,bhofset,      &
                                       sfrac1,fe_tk,rpivot,zpivot,naptr,      &
                                       nsourc,npart_all_beamlines,npart,      &
                                       npskip,randomize_seed,fionx,fdbeam,    &
                                       mstate,nbion,iterate_beam,             &
                                       izstrp,iexcit,ilorent,npskip,iborb,    &
                                       ncont,ngl,ngh,kdeni,kdenz,ksvi,ksvz,   &
                                       ksve,krad,kdene,ibion,atw_beam,        &
                                       beam_sim_time_start,beam_sim_time_end, &
                                       split_injectors,nfreya_vb,use_ufile,   &
                                       ne_tk
                                       

    USE zonal_data,            ONLY :  mf,mfm1

    USE xsct,                  ONLY:   adas_xsct_path


    IMPLICIT NONE

    INTEGER(I4B) itrmax,i,j,nsend,sz

    CHARACTER(len = name_size) dumy

    LOGICAL ldum,read_ufile,read_nubeam_namelist

    REAL(DP) rdum

!----------------------------------------------------------------------------------------------
! fast ion diffusion using parameters in module fast_ion_diffusion are not part of
! the rn directives at this time. But the off condition of the fas tion diffusion model must
! be signalled to the remaining processors. So we do it here in case nubeam namelist
! is not used in this run . Note that these a rebroadcast in distribute nubeam namelist
! if nubeam file is read.
!----------------------------------------------------------------------------HSJ----------------
         sz = SIZE(d_fast_ion%adiff_a)        ! the adiff params are hard wired to dimension of 3
                                              ! sz =3 is known on all processors but send it anyway
                                              ! for uniformity in codeing.

         CALL MPI_BCAST(d_fast_ion%adiff_a,sz,MPI_DOUBLE_PRECISION,master, &
                        MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(d_fast_ion%adiff_0,sz,MPI_DOUBLE_PRECISION,master, &
                        MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(d_fast_ion%adiff_xpin,sz,MPI_DOUBLE_PRECISION,master, &
                        MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(d_fast_ion%adiff_xpout,sz,MPI_DOUBLE_PRECISION,master, &
                        MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(d_fast_ion%fidif_on,1,MPI_INTEGER,master, &
                        MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(d_fast_ion%nkdifb,1,MPI_INTEGER,master, &
                        MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(d_fast_ion%ndifbe,1,MPI_INTEGER,master, &
                        MPI_COMM_WORLD,mpiierr)





    CALL MPI_BCAST(nbeams,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mstate,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(iexcit,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ilorent,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(naptr,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(nsourc,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(mf,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    mfm1 = mf -1
    CALL MPI_BCAST(ne_tk,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(npart,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
                   npart_all_beamlines = npart
    CALL MPI_BCAST(npskip,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(iborb,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ibion,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(nbion,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ncont,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ngl,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ngh,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(kdeni,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(kdene,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(kdenz,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ksvi,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ksvz,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ksve,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(krad,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)

    CALL MPI_BCAST(izstrp,kimp,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
    nsend = SIZE(ni_sc)
    CALL MPI_BCAST(ni_sc,nsend,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)



    CALL MPI_BCAST(timbplt,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(atw_beam,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(relnub,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(fe_tk,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(fionx,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(fdbeam,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(beam_sim_time_start,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(beam_sim_time_end,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)


    CALL MPI_BCAST(beam_on,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(beam_time,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(anglev,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(angleh,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bcur,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bptor,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(blenp,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bleni,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bheigh,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bwidth,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bhfoc,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bvfoc,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bhdiv,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bvdiv,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(ebkev,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bvofset,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bhofset,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(sfrac1,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(rpivot,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(zpivot,kb,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)



    nsend = ke*kb
    CALL MPI_BCAST(fbcur,nsend,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    nsend = nap*kb
    CALL MPI_BCAST(aheigh,nsend,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(awidth,nsend,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(alen,nsend,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

    IF(myid .ne. master)no_injectors = nbeams

    CALL MPI_BCAST(randomize_seed,1,MPI_LOGICAL,master,MPI_COMM_WORLD,mpiierr)

    CALL MPI_BCAST(iterate_beam,1,MPI_LOGICAL,master,MPI_COMM_WORLD,mpiierr)

    CALL MPI_BCAST(split_injectors,1,MPI_LOGICAL,master,MPI_COMM_WORLD,mpiierr)

    CALL MPI_BCAST(nfreya_vb,1,MPI_LOGICAL,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(use_ufile,1,MPI_LOGICAL,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(read_nubeam_namelist,1,MPI_LOGICAL,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(read_ufile,1,MPI_LOGICAL,master,MPI_COMM_WORLD,mpiierr)
    CALL MPI_BCAST(bfr_neutrlz,1,MPI_LOGICAL,master,MPI_COMM_WORLD,mpiierr)

    DO j=1,SIZE(namelist_nameb)
       dumy(:) = namelist_nameb(j)(:) ! meaningful on process 0
       CALL MPI_BCAST(dumy,name_size,MPI_CHARACTER,master,MPI_COMM_WORLD,mpiierr)
       namelist_nameb(j) = dumy(:) !now on all processes
    ENDDO

    DO j=1,kb
       dumy(:) = nbshape(j)(:) ! meaningful on process 0
       CALL MPI_BCAST(dumy,name_size,MPI_CHARACTER,master,MPI_COMM_WORLD,mpiierr)
       nbshape(j) = dumy(:)        !now on all processes
    ENDDO

    
    DO i=1,SIZE(nashape,1)
     DO j=1,SIZE(nashape,2)
        dumy(:) = nashape(i,j)(:)
        CALL MPI_BCAST(dumy,name_size,MPI_CHARACTER,master,MPI_COMM_WORLD,mpiierr)
        nashape(i,j) = dumy(:)
     ENDDO
    ENDDO

    nsend =LEN(adas_xsct_path)
    CALL MPI_BCAST(adas_xsct_path,nsend,MPI_CHARACTER,master,MPI_COMM_WORLD,mpiierr)

    !nfreya_plot_file_name used only on master

    ! distribute fast ion diffusion data here since the run may not use
    ! the nubeam namelist
    


#endif
              


    RETURN
    END SUBROUTINE distribute_run_directives 
