  PROGRAM NFREYA_PARALLEL
!---------------------------------------------------------------------------
! -- Stand alone version of nfreya with independent beamline timing
! -- Each cpu is assigned a beamline on a per needed basis using
! -- the master/slave programming paradigm.
! -- Uses statefile input from Onetwo  to get kinetic data
! -- Uses namelist input to get necessary beam data
! -- Uses ufile input on option
!------------------------------------------------------------HSJ 1/6/2011---


    USE nrtype,                     ONLY : I4B,DP

    USE common_constants,           ONLY : izero

    USE error_handler,              ONLY : load_errno,lerrno,iomaxerr,terminate

 
    USE io_gcnmp,                   ONLY : ncrt,nlog

    USE nfreya_run,                 ONLY : simulation_time_advance,        &
                                           load_plasma_prop

    USE ran_num_gen,                ONLY : set_random_number_seed

    USE neutral_beams,              ONLY : write_performance_data

    USE prep_tport_outpt,           ONLY : set_profiles

    USE Nfreya_routines,            ONLY : process_nubplt 

    USE P_nfreya_interface,         ONLY : beam_data


#if defined (USEMPI)
    USE mpi
    USE MPI_data,                   ONLY : myid,mpiierr,mpi_start_time,    &
                                           mpi_end_time,numprocs,master,   &
                                           MPI_profile,local_proc_name,    &
                                           mpiversion,mpisubversion,       &
                                           mstr_window,mstr_active
#endif


#ifndef USEMPI 
   USE MPI_data,                   ONLY : myid,mpiierr,mpi_start_time,     &


                                          mpi_end_time,numprocs,master,    &
                                          mstr_active
#endif





    IMPLICIT NONE
 
    INTEGER(I4B) ptyp,rtyp,strsize,task,error,include_tglf
    REAL(DP)  start_time,end_time,run_time



 
      master = izero                                 ! process 0 is master
      myid   = master                                ! if compiled without MPI
      numprocs = 1                                   ! interface
      mstr_active = .FALSE.                          ! set to true below if rma window is used                                 
#if defined (USEMPI)

      CALL MPI_Init(mpiierr)
!     initialized_mpi in module MPI  is now .TRUE. query its value 
!     with routine MPI_initialized

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,mpiierr) !get processor id
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, mpiierr )
      CALL MPI_GET_PROCESSOR_NAME(local_proc_name,strsize,mpiierr)
      CALL MPI_GET_VERSION(mpiversion,mpisubversion,mpiierr)
      IF(mpiierr .NE. 0)THEN
         WRITE(ncrt,FMT='("ERROR initiaizing mpi in PROGRAM GCNM_PARALLEL")')
         lerrno =0
         CALL MPI_ABORT(MPI_COMM_WORLD,lerrno,mpiierr)
      ENDIF
      IF(myid == master) start_time = MPI_WTIME()
#endif

!   nullify pointers:
         CALL null_pointers                  ! utils.F90
 
!   load the error handler :
         CALL load_errno                     ! error_handler.F90

! --  get input info from  command line arguments.
! --  input consists  of run directives and beam info  from a namelists  and  
! --  kinetic information from statefile.

        CALL load_nfreya_input               ! load_nfreya_input.f90

        CALL set_random_number_seed          ! ran_num_gen.F90



        CALL nfreya_init                     ! nfreya_init.f90: initialize beam geometry,
                                             ! set zonal values,
                                             ! determine cross sections


  
        CALL simulation_time_advance         ! advance the beams from beam_sim_time_start 
                                             ! to beam_sim_time_end (nfreya_run.F90)

 

        IF(myid == master)THEN
    
           !check of beamlets is doen here because npart_pseudo was not known
           ! until afer call to freyas
           CALL Check_P_Nfreya_beamlets(beam_data%nbeam)


           CALL load_plasma_prop             ! nfreya_run.F90 
                                             ! adjust units of some beam related parameters
                                             ! after this call do not use pseudo injectors anymore
 
           CALL set_profiles                 ! prep_transport_ouput.f90
                                             ! prepare output profiles,
                                             ! may interpolate onto different rho grid
                                             ! Adjust units to mks
 
           CALL write_statefile(rtyp)        ! nfreya_write_statefile.f90
                                             ! Write statefile output on grid nj_out

           CALL process_nubplt               ! 


        ENDIF


#if defined (USEMPI)
      include_tglf = 0
      CALL MPI_profile(ncrt,include_tglf)
      CALL MPI_profile(nlog,include_tglf)
      IF(myid == master)THEN
         end_time = MPI_WTIME()
         run_time = end_time - start_time 
         WRITE(ncrt,10)numprocs,run_time
         WRITE(nlog,10)numprocs,run_time
         IF(write_performance_data .NE. 'NONE') &
              call write_perf_data(numprocs,run_time)
!         WRITE(ncrt,20)icalln
!         WRITE(nlog,20)icalln
         CALL timestamp(ncrt)
         CALL timestamp(nlog)
      ENDIF
      IF(mstr_active)CALL MPI_WIN_FREE(mstr_window, mpiierr) 
      CALL MPI_FINALIZE(mpiierr)
#endif


#ifndef USEMPI
!      WRITE(ncrt,20)icalln
!      WRITE(nlog,20)icalln
 
      CALL timestamp(ncrt)
      CALL timestamp(nlog)
#endif

 10   FORMAT(2X,'num procs, elapsed time =',2X,i5,1pe14.2)
 20   FORMAT(2x,'num model evaluations     :',2x,i5)
 30   FORMAT(2x,'num Jacobian evaluations  :',2x,i5)


  END PROGRAM NFREYA_PARALLEL


    SUBROUTINE debug_stop(message) 
!-----------------------------------------------------------------------------
! -- used to stop code in various places for debugging
! -- (A poor mans debugger for multiple processor codes)
! -- see  sub terminate in module error_handler for permanent error stops
!-----------------------------------------------------------------------------
#ifdef USEMPI
    USE mpi
#endif
    USE MPI_data,                   ONLY : myid,mpiierr,master
    USE io_gcnmp,                   ONLY : ncrt,nlog
    IMPLICIT NONE
    CHARACTER *(*) message
    IF(myid == master) THEN
       WRITE(ncrt,1)message
       WRITE(nlog,1)message
    ENDIF
1   FORMAT(2x,'code called debug stop with message:  ',a)

#ifdef USEMPI
      mpiierr =-1
      CALL MPI_FINALIZE(mpiierr)
#endif
      CALL EXIT


    END     SUBROUTINE debug_stop





    SUBROUTINE write_perf_data(numprocs,run_time)
!------------------------------------------------------------------
! -- for keeping track of mpi performance issues this
! -- routine tracks the timming
!------------------------------------------------------------------
      USE nrtype,                      ONLY : DP,I4B,I2B
   
      USE nf_param,                    ONLY : kb

      USE neutral_beams,               ONLY : no_physical_injectors,  &
                                              write_performance_data, &
                                              npart_all_beamlines,    &
                                              split_injectors

      USE error_handler,               ONLY : terminate,iomaxerr,lerrno

      USE io_gcnmp,                    ONLY : ncrt,nlog

      USE nfreya_version,              ONLY : p_nfreya_ver



      IMPLICIT NONE

      INTEGER(I4B) numprocs
      INTEGER(I2B) iounit, get_next_io_unit 
      REAL(DP) run_time
      LOGICAL exists

      iounit = get_next_io_unit()

      ! does write_performance_data file exist ?
      INQUIRE (FILE = write_performance_data, EXIST  = exists)
      IF(exists)THEN            ! open for appending
         WRITE(ncrt,FMT='(" appending to  file " ,a)')TRIM(write_performance_data)
         OPEN(unit = iounit, file = write_performance_data ,status = 'OLD', &
                                   POSITION = 'APPEND',err = 1)
      ELSE                  ! create file
         WRITE(ncrt,FMT='(" Creating   file ",a )')TRIM(write_performance_data)
         OPEN(unit = iounit, file = write_performance_data , status = 'UNKNOWN',err = 2)
         !WRITE(iounit,15)
      ENDIF
      CALL timestamp (iounit )
      WRITE(iounit,10)no_physical_injectors,numprocs,kb,npart_all_beamlines,run_time,split_injectors 

      RETURN

1     lerrno = iomaxerr + 240_I4B
      CALL terminate(lerrno,nlog)
 
2     lerrno = iomaxerr + 241_I4B
      CALL terminate(lerrno,nlog)

10    FORMAT(2x,"No injectors = ",i3," No cpus used/allowed =",i3,1x,i3," tl npart =",i8," run_time,sec =",1pe12.2," si = ",l1)
15    FORMAT(2x,"date         P_Nfreya version  no injectors  no cpus  run time")
    END     SUBROUTINE write_perf_data
