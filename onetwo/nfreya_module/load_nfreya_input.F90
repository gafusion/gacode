
  SUBROUTINE load_nfreya_input
!-------------------------------------------------------------------------
! read the initial conditions from an statefile type file
! The name of the state file is  taken from the first argument of the command line.
! The second argument gives  the name of a namelist input file that contains
! run directives.
! The third argument is the name of the logfile for errors
!-----------------------------------------------------------------HSJ 1/13/11----
   USE nrtype,    ONLY : I2B , DP,I4B
   
   USE nfreya_version

   USE error_handler,          ONLY : lerrno,terminate,iomaxerr

   USE iterdbmd_gcnmp,         ONLY : iterdb_file_name,irwflag,idb_sufix

   USE io_gcnmp,               ONLY : runlog_filename,namelist_filename,nlog

   USE ions_gcnmp,             ONLY : nprim,nimp,nion,fi_index,ni_sc
  
   USE neutral_beams,          ONLY : nbion,use_ufile,ibion

   USE Plasma_properties,      ONLY : mhd_dat,neut_beam

   USE solcon_gcnmp,           ONLY : time0, time
                        
   USE MPI_data,               ONLY : myid,numprocs,master,mpiierr

   USE grid_class,             ONLY : meshgen,mhd2tport

   USE bc_values_gcnmp,        ONLY : set_bc,bc_conditions_init,set_bc_type,  &
                                      gammat_bc
   
   USE nfreya_namelist,         ONLY : read_run_directives_namelist

   USE P_nfreya_interface,      ONLY : beam_data,get_nf_beam_data

   IMPLICIT NONE
   INTEGER(I2B),PARAMETER :: arglist = 9
   INTEGER(I2b) get_next_io_unit
   INTEGER(I2B) numarg,i,iounit,j,m
   CHARACTER(len=256) carg(0:arglist) 
   LOGICAL set_dmassden,read_ufile,read_nubeam_namelist

     INTERFACE
        SUBROUTINE my_getarg(carg,numarg,unit)
          USE nrtype,  ONLY : DP,I4B,I2B
          INTEGER(I2B),INTENT(OUT) :: numarg
          INTEGER(I4B),INTENT(IN),OPTIONAL :: unit
          CHARACTER ( len = *),DIMENSION(0:),INTENT(OUT) ::  carg
        END  SUBROUTINE my_getarg
        SUBROUTINE to_upper_case(string)
          USE nrtype,            ONLY : I4B
          IMPLICIT NONE
          INTEGER(I4B) l
          CHARACTER*(*), INTENT (INOUT) :: string
        END SUBROUTINE to_upper_case
     END INTERFACE


    nlog = 6    !initial output is only to screen
                !reset below once logfile name is known

!-------------------------------------------------------------
! set global namelist defaults for run directives namelist
! (all processors call this)
!-------------------------------------------------------------
    CALL nfreya_namelist_defaults               ! nfreya_namelist_defaults.f90


!-------------------------------------------------------------
!   get command line arguments and process the input files:
!-------------------------------------------------------------
    carg(0:arglist)(:) = ' '
    IF(myid == master) THEN
           CALL my_getarg(carg,numarg)

           IF(numarg .LT. 3) THEN
              lerrno = 200_I4B + iomaxerr
              CALL terminate(lerrno,nlog)
           ENDIF
  
  
           carg(1) =ADJUSTL(carg(1))

           iterdb_file_name = carg(1)(1:LEN_TRIM(carg(1)))
           carg(2) = ADJUSTL(carg(2))
           namelist_filename = carg(2)(1:LEN_TRIM(carg(2)))

           carg(3) =ADJUSTL(carg(3))
           runlog_filename = carg(3)(1:LEN_TRIM(carg(3)))


           !create the output  file:
           nlog = get_next_io_unit ()
           OPEN (unit = nlog , file = runlog_filename, status = 'UNKNOWN',err = 1)
           WRITE(nlog,10)carg(0)(1:LEN_TRIM(carg(0))),p_nfreya_ver
           WRITE(nlog,11)carg(1)(1:LEN_TRIM(carg(1)))
           WRITE(nlog,12)carg(2)(1:LEN_TRIM(carg(2)))
           WRITE(nlog,13)carg(3)(1:LEN_TRIM(carg(3)))

 10        FORMAT(2X,'executable name    : ',a,X,a)
 11        FORMAT(2X,'iterdb file name   : ',a)
 12        FORMAT(2X,'input  file name   : ',a)
 13        FORMAT(2X,'log file name      : ',a)

           ! read namelist that tells P_Nfreya how to run
           CALL read_run_directives_namelist     ! nfreya_namelist.F90


!-----------------------------------------------------------------------------
!          read the statefile file. Note that nprim,nimp,nneu are obtained
!          here and will be needed to read the namelist data in sub
!          read_namelist. Hence the order of reading the files is 
!          critical
!-----------------------------------------------------------------------------

           irwflag = 1 
           j = LEN_TRIM(iterdb_file_name)
           idb_sufix = iterdb_file_name(j-2:j+1)
           CALL to_upper_case(idb_sufix)
           IF ( idb_sufix == '.NC')THEN
              CALL iter_dbase_nc                  ! rw_iterdb_netcdf.F90
           ELSE
              CALL iter_dbase_txt                 ! rw_iterdb.F90
           ENDIF

           ! transfer mhd grid npsi metrics to nj psir grid:
           CALL mhd2tport                         ! grid_class.F90

           CALL set_ion_prop2      ! set_ion_prop2.f90
                                   ! primary and beam ion name processing (set_ion_prop2.f90)



           !----------------------------------------------------------------------
           ! check on usage of nlbdat:
           ! if use_ufile is set in run directives namelist
           ! but nlbdat is not consistent then
           ! use_ufile takes precendence
           !----------------------------------------------------------------------
           beam_data%nlbdat_set    = .FALSE.
           IF( use_ufile)THEN
              beam_data%nlbdat     = .TRUE.
              beam_data%nlbdat_set = .TRUE.
           ENDIF



           !-----------------------------------------------------------
           ! read beam geometry  namelist in file beam_data_namelist
           ! and possibly a ufile named  beam_data_ufile.
           ! The names, beam_data_namelist,beam_data_ufile,and use_ufile,
           ! were set in run directives file read above.
           !------------------------------------------------------------
           beam_data%data_allocated = .FALSE.
           read_nubeam_namelist     = .FALSE.
           CALL get_nf_beam_data(beam_data,use_ufile,     &    ! P_nfeya_interface.F90 
                      read_ufile,read_nubeam_namelist)         ! loads beam_data structure 
                                                               ! defined as TYPE(neutral_beam)
                                                               ! uses atomno,atw from set_ion_prop2
 
          IF(neut_beam%nbeams .LE. 0)neut_beam%nbeams = beam_data%nbeam
    ENDIF !(myid == master)


!-------------------------------------------------------------
! --   if this is a multiple processor run pass out the input to
! --   other processors. Allocate some arrays for all processors.
! --   (Must be called by all processes) 
! --    NOTE that distribute_xsct is also used to distribute xsct data 
! --    later on.
!-------------------------------------------------------------
    IF (numprocs > 1)THEN
#ifdef USEMPI
       CALL distribute_run_directives(read_ufile,read_nubeam_namelist)     ! distribute_run_directives.F90

       CALL distribute_statefile          ! distribute_statefile.F90

       IF(read_nubeam_namelist) CALL distribute_nubeam_namelist

       IF(read_ufile) CALL distribute_ufile

 
       ! CALL debug_stop('debug stop line 178 load_nfreya') 
!       IF(use_ufile)THEN
!          CALL distribute_nubeam_namelist ! distribute_nubeam_namelist.f90
!          CALL debug_stop('debug stop after dstrib load_nfreya') 
!          CALL distribute_ufile           ! distribute_ufile.f90
!       ENDIF
#endif
    ENDIF



!   setup required meshes
    CALL meshgen         ! grid_class.F90


!   set dependent variables species  to be run in analysis/simulation:
!    nion = nprim + nimp
    CALL allocate_itran_species(nprim,nimp)      ! load_nfreya_input.f90


!   Scale the dependent variable profiles
!   1) densities (profile%ene,profile%ni
!       CALL scale_den                    ! scale_profiles.f90
!   2) te
!       CALL scale_te                     ! scale_profiles.f90
!   3) ti
!       CALL scale_ti                     ! scale_profiles.f90
!   4) toroidal rotation
!       CALL scale_trot                   ! scale_profiles.f90

    set_dmassden = .TRUE. 
!    CALL set_ion_prop(set_dmassden)

    RETURN
1   lerrno = 5
    CALL terminate(lerrno,nlog)
    RETURN
  END SUBROUTINE load_nfreya_input


  SUBROUTINE allocate_itran_species(nprim,nimp)
    USE nrtype,                      ONLY : DP,I4B
    USE solcon_gcnmp,                ONLY : itenp,iteni,ntot
    USE bc_values_gcnmp,             ONLY : bc_type,mult_den,mult_flux,gammat_bc
    USE ions_gcnmp,                  ONLY : ni_sc
    USE dep_var,                     ONLY : dp4
    IMPLICIT NONE

    INTEGER(I4B),INTENT(IN) :: nprim,nimp

    IF(.NOT. ALLOCATED(bc_type))ALLOCATE(bc_type(nprim+nimp+dp4))
    IF(.NOT. ALLOCATED(gammat_bc))ALLOCATE(gammat_bc(nprim+nimp))
    IF(.NOT. ALLOCATED(mult_den))ALLOCATE(mult_den(nprim+nimp))
    IF(.NOT. ALLOCATED(mult_flux))ALLOCATE(mult_flux(nprim+nimp))
    IF(.NOT. ALLOCATED(iteni)) ALLOCATE(iteni(nimp))
    IF(.NOT. ALLOCATED(itenp)) ALLOCATE(itenp(nprim))
    IF(.NOT. ASSOCIATED(ni_sc)) ALLOCATE(ni_sc(4))
    RETURN

  END  SUBROUTINE allocate_itran_species



    SUBROUTINE strip_path(path,prefix,extension,filename,task,ishot,limit)
!--------------------------------------------------------------------------
! return path, shot number and filename (with path removed)
! task   input values
!        if task   =0 then try to get shot number,
!        if task   =1 return after determining path and file names.
!------------------------------------------------------------HSJ-11/18/03--

   USE error_handler,                          ONLY : lerrno,terminate,iomaxerr
  
   USE io_gcnmp,                               ONLY : ncrt,nlog

   IMPLICIT NONE
   INTEGER j,k,l,ll,lll,task,lp
   INTEGER, INTENT (out) :: ishot
   CHARACTER*(*), INTENT(inout):: filename
   CHARACTER*(*), INTENT(out)  :: path,prefix,extension !JMP
   INTEGER,  INTENT(in):: limit      !limit search to this many characters
                                     !ufiles has 64 character limit
   CHARACTER(len=6):: shot,buf

   l = LEN_TRIM(filename)

   lp = LEN(path)
   IF(l .LE. limit)THEN
     DO j=l-1,2,-1
        k = j
        IF(filename(j:j) .EQ. '/')EXIT
     ENDDO
     IF(k .EQ. 2)THEN
        !no path specification found set to local path:
        IF(lp .GE. 2)THEN
           path = "./"
        ELSE
           WRITE(ncrt,3)path
3          FORMAT(2x,'path character string too short',a)
           lerrno = 246 + iomaxerr
           CALL  terminate(lerrno,nlog)
        ENDIF
     ELSE  ! found '/' in filename, assume it is path specifier
        IF(lp .GE. k)THEN
           path = filename(1:k)
        ELSE
           WRITE(ncrt,3)path
           lerrno = 246 + iomaxerr
           CALL  terminate(lerrno,nlog)
        ENDIF
        !rotate filename until it is at the begining of the string:
        filename(1:k) =' '
        filename = ADJUSTL(filename) 
        l = LEN_TRIM(filename)
     ENDIF 

     If(task == 1) RETURN ! return with just path and filename set
     ll =SCAN(filename(1:l),'0123456789')
     prefix = filename(1:ll-1)
     lll = SCAN(filename(ll:l),'.')
     extension = filename(lll+2:l) !JMP
     IF (lll-ll .LE. 6)THEN
            shot = filename(ll:lll)
            !convert to integer
            WRITE(buf,FMT='(a)')shot
            READ(buf,FMT='(i6)')ishot

     ELSE
        WRITE(ncrt,2)filename(1:l)
 2      FORMAT(2x,'could not determine shot no. for',a)
        lerrno = 246 + iomaxerr
        CALL  terminate(lerrno,nlog)
     ENDIF

   ELSE

     WRITE(ncrt,1)filename(1:LEN_TRIM(filename))
 1   FORMAT(2x,'Error, ufile routines require that ',a,/, &
            2x,'be 64 characters or less in length')
        lerrno = 247 + iomaxerr
        CALL  terminate(lerrno,nlog)
   ENDIF

   RETURN
   END SUBROUTINE strip_path
