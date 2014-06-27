
#ifdef ONETWO
  SUBROUTINE load_input
   ! dumy link for Onetwo
  END SUBROUTINE load_input
#endif
#ifdef GCNMP
  SUBROUTINE load_input
!-------------------------------------------------------------------------
! read the initial conditions from an statefile type file
! The name of the state file is  taken from the first argument of the command line.
! The second argument gives  the name of a namelist input file that contains
! run directives.
! The third argument is the name of the logfile for errors
!-----------------------------------------------------------------HSJ-----
   USE nrtype,    ONLY : I2B , DP,I4B
   
   USE gcnmp_version

   USE error_handler

   USE iterdbmd_gcnmp,         ONLY : iterdb_file_name,irwflag,idb_sufix

   USE io_gcnmp,               ONLY : runlog_filename,namelist_filename,nlog

   USE ions_gcnmp,             ONLY : nprim,nimp,nion,fi_index,ni_sc
  
   USE neutral_beams,          ONLY : nbion

   USE Plasma_properties,      ONLY : mhd_dat

   USE solcon_gcnmp,           ONLY : time0, time
                        
   USE MPI_data,               ONLY : myid,numprocs,master,mpiierr

   USE grid_class,             ONLY : meshgen,mhd2tport

   USE bc_values_gcnmp,        ONLY : set_bc,bc_conditions_init,set_bc_type,  &
                                      gammat_bc
   
   USE gcnmp_namelist,         ONLY : read_namelist

   USE  common_constants,      ONLY : zeroc 

   USE plasma_shape,           ONLY : process_shape


   IMPLICIT NONE
   INTEGER(I2B),PARAMETER :: arglist = 9
   INTEGER(I2b) get_next_io_unit
   INTEGER(I2B) numarg,i,iounit,j,m
   CHARACTER(len=256) carg(0:arglist) 
   LOGICAL set_dmassden

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
    carg(0:arglist)(:) = ' '
!   get command line arguments:
    if(myid == master) THEN
       CALL my_getarg(carg,numarg)
       IF(numarg .LT. 3) THEN
          lerrno =1_I4B
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
       WRITE(nlog,10)carg(0)(1:LEN_TRIM(carg(0))),gcnmp_ver
       WRITE(nlog,11)carg(1)(1:LEN_TRIM(carg(1)))
       WRITE(nlog,12)carg(2)(1:LEN_TRIM(carg(2)))
       WRITE(nlog,13)carg(3)(1:LEN_TRIM(carg(3)))

 10     FORMAT(2X,'executable name    : ',a,X,a)
 11     FORMAT(2X,'iterdb file name   : ',a)
 12     FORMAT(2X,'input  file name   : ',a)
 13     FORMAT(2X,'log file name      : ',a)


       !----------------------------------------------------------------
       !   read the state file. Note that nprim,nimp,nneu are obtained
       !   here and will be needed to read the namelist data in sub
       !   read_namelist. Hence the order of reading the files is 
       !   critical
       !-----------------------------------------------------------------
       irwflag = 1 
       j = LEN_TRIM(iterdb_file_name)
       idb_sufix = iterdb_file_name(j-2:j+1)
       CALL to_upper_case(idb_sufix)
       if ( idb_sufix == '.NC')THEN
          CALL iter_dbase_nc
       ELSE
          CALL iter_dbase_txt
       ENDIF

       !transfer mhd grid npsi metrics to nj psir grid:
       CALL mhd2tport

    ENDIF !(myid == master)

!   if this is a multiple processor run pass out the input to
!   other procesors. Allocate some arrays for all processors.
!   (Must be called by all processes) 

    CALL distribute_statefile  ! distribute_statefile.F90

!   setup required meshes
    CALL meshgen
!   set dependent variables species  to be run in analysis/simulation:
!    nion = nprim + nimp
    CALL allocate_itran_species(nprim,nimp)

!   read the run directives from the namelist input file
    IF(myid == master) THEN
       CALL read_namelist ! has dependency on nj through 
                          ! currf and curbeam scale factors

       CALL set_bc_type   ! sets bc_type (Robin or Dirichlet) for densities only

    ENDIF


!   note that some quantites are used only
!   by the master process and hence are not spread to other 
!   processes
    if(numprocs > 1) CALL distribute_namelist

! some quantitites are not in statefile (eq rminor) calculate it here:
    CALL process_shape



!   set dependent variables to be run in analysis/simulation:
    CALL set_itran        !assumes  distribute_namelist was called
    
    set_dmassden = .TRUE. 
    CALL set_ion_prop(set_dmassden)

!   setup  boundary condition items that dont change during this call
!   to GCNMP 
    set_bc = .FALSE.
    CALL bc_conditions_init


!   Scale the dependent variable profiles
!   1) densities: profile%ene,profile%ni
!      ni_sc(1) = nictr, ni_sc(2) = niavg,ni_sc(3) = nistd,ni_sc(4) = nirho

       If(   ni_sc(1) .GT. zeroc .AND.  ni_sc(2) .GT. zeroc   &
             .AND. ni_sc(3) .GT. zeroc .AND. ni_sc(4) .GT. zeroc)THEN
           CALL scale_den
       ELSEIF(  ni_sc(1) .GT. zeroc .AND. mhd_dat%ni_spline .GT. 0 )THEN

           CALL scale_betan_ped

           ! densities at location of pedestal are now consistent
           ! with given input betan_ped
           ! now scale density to get necessary fusion power
           ! and recalculate ene,zeff
           CALL calc_densities

       ENDIF

!   2) te
       CALL scale_te ! not active

!   3) ti
       CALL scale_ti ! not active


!   4) toroidal rotation
       CALL scale_trot ! not active


    RETURN


1   lerrno = 5
    CALL terminate(lerrno,nlog)


    RETURN
  END SUBROUTINE load_input


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
#endif
