    SUBROUTINE distribute_ufile
!--------------------------------------------------------------------------
! -- pass data read from ufile by master process to other processes
! -- associated with this run
!---------------------------------------------------------------------------
#ifdef USEMPI

    USE nrtype,                ONLY : DP,I4B,I2B

    USE MPI

    USE MPI_data,              ONLY : myid,master,mpiierr,numprocs

    USE neutral_beams,         ONLY : bptor,ebkev,fbcur,sfrac1

    USE P_nfreya_interface,    ONLY : name_beam => beam_data

    USE beam_structure,        ONLY : beam_structure_allocate

    USE ions_gcnmp,            ONLY : nprim,nion,nimp

    IMPLICIT NONE
     
    INTEGER(I4B)  ngmax,nrhix,nbeam,nfast,ucdim,lcdim,utdim,ltdim,sz




         IF(myid == master)THEN
            utdim  = UBOUND(name_beam%beam_times,1)
            ltdim  = LBOUND(name_beam%beam_times,1)
            ucdim  = UBOUND(name_beam%beam_chan,1)
            lcdim  = LBOUND(name_beam%beam_chan,1)
         ENDIF
         CALL MPI_BCAST(utdim,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(ltdim,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(ucdim,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(lcdim,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(name_beam%nbeam,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(name_beam%nrhix,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(name_beam%ngmax,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)

         nbeam = name_beam%nbeam
         IF(.NOT. name_beam%data_allocated) THEN
               nrhix =  nimp
               ngmax =  nion
               nfast =  2
               CALL beam_structure_allocate(name_beam,nbeam,nrhix,ngmax,nfast) ! beam_structure.f90

         ENDIF



         IF(myid .NE. master)THEN ! non master processes may have been initialized 
                                  ! with default sizes. Correct this here:
            IF( ASSOCIATED(name_beam%beam_times))   DEALLOCATE(name_beam%beam_times)
                ALLOCATE(name_beam%beam_times(ltdim:utdim))

            IF(  ASSOCIATED(name_beam%beam_inject)) DEALLOCATE(name_beam%beam_inject)
                ALLOCATE(name_beam%beam_inject(ltdim:utdim,lcdim:ucdim))

            IF( ASSOCIATED( name_beam%beam_chan))  DEALLOCATE(name_beam%beam_chan)
                ALLOCATE(name_beam%beam_chan(lcdim:ucdim))



         ENDIF

         CALL MPI_BCAST(name_beam%nbeam,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(name_beam%beam_inject,SIZE(name_beam%beam_inject), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
         CALL MPI_BCAST(name_beam%beam_times,SIZE(name_beam%beam_times), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

         CALL MPI_BCAST(name_beam%beam_chan,SIZE(name_beam%beam_chan), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)



         
        CALL MPI_BCAST(name_beam%einja,SIZE(name_beam%einja), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
        CALL MPI_BCAST(name_beam%pinja,SIZE(name_beam%pinja), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

        CALL MPI_BCAST(name_beam%ffulla,SIZE(name_beam%ffulla), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
        CALL MPI_BCAST(name_beam%fhalfa,SIZE(name_beam%fhalfa), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

        CALL MPI_BCAST(bptor,SIZE(bptor), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
        CALL MPI_BCAST(ebkev,SIZE(ebkev), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
        CALL MPI_BCAST(fbcur,SIZE(fbcur), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
        CALL MPI_BCAST(sfrac1,SIZE(sfrac1), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

        CALL MPI_BCAST(name_beam%tbona,SIZE(name_beam%tbona), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
        CALL MPI_BCAST(name_beam%tboffa,SIZE(name_beam%tboffa), &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

  If (myid == master)sz = SIZE(name_beam%xzimpx)
     CALL MPI_BCAST(sz,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
     If (myid .NE. master)THEN
        IF(ASSOCIATED(name_beam%xzimpx))DEALLOCATE(name_beam%xzimpx)
        ALLOCATE(name_beam%xzimpx(sz))
     ENDIF
     CALL MPI_BCAST(name_beam%xzimpx,sz, &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

  If (myid == master)sz = SIZE(name_beam%aimpx)
     CALL MPI_BCAST(sz,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
     If (myid .NE. master)THEN
        IF(ASSOCIATED(name_beam%aimpx))DEALLOCATE(name_beam%aimpx)
        ALLOCATE(name_beam%aimpx(sz))
     ENDIF
     CALL MPI_BCAST(name_beam%aimpx,sz, &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

  If (myid == master)sz = SIZE(name_beam%backz)
     CALL MPI_BCAST(sz,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
     If (myid .NE. master)THEN
        IF(ASSOCIATED(name_beam%backz))DEALLOCATE(name_beam%backz)
        ALLOCATE(name_beam%backz(sz))
     ENDIF
     CALL MPI_BCAST(name_beam%backz,sz, &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

  If (myid == master)sz = SIZE(name_beam%aplasm)
     CALL MPI_BCAST(sz,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
     If (myid .NE. master)THEN
        IF(ASSOCIATED(name_beam%aplasm))DEALLOCATE(name_beam%aplasm)
        ALLOCATE(name_beam%aplasm(sz))
     ENDIF
     CALL MPI_BCAST(name_beam%aplasm,sz, &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

  If (myid == master)sz = SIZE(name_beam%beam_times)
     CALL MPI_BCAST(sz,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
     If (myid .NE. master)THEN
        IF(ASSOCIATED(name_beam%beam_times))DEALLOCATE(name_beam%beam_times)
        ALLOCATE(name_beam%beam_times(sz))
     ENDIF
     CALL MPI_BCAST(name_beam%beam_times,sz, &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
  sz =0
  If (myid == master .AND. ASSOCIATED(name_beam%beam_switch_times)) &
       sz = SIZE(name_beam%beam_switch_times)
     CALL MPI_BCAST(sz,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpiierr)
     If (myid .NE. master .AND. sz .NE.0)THEN
        IF(ASSOCIATED(name_beam%beam_switch_times))DEALLOCATE(name_beam%beam_switch_times)
        ALLOCATE(name_beam%beam_switch_times(sz))
     ENDIF
     IF(sz .NE. 0)CALL MPI_BCAST(name_beam%beam_switch_times,sz, &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

  sz =1
  CALL MPI_BCAST(name_beam%beam_power_rise_time,sz, &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
  CALL MPI_BCAST(name_beam%beam_power_decay_time,sz, &
                       MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

#endif
        RETURN
    END SUBROUTINE distribute_ufile
