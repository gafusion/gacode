
  SUBROUTINE set_ion_prop2
!----------------------------------------------------------------------------------------
!     This is an updated version of the gcnmp set_ion_prop routine. It should be
!     back ported into gcnmp eventually.
!     The important update is that the names of the ions do not have to
!     appear anywhere except in the original definition. So future changes will
!     be much less troublesome.
!     NEUTRALS done  elesewhere ?
!  1) check all ions against respective master lists. From index in maser list we can
!     get mass and charge numbers of ions used in an actual  run.
!  2) check beam species name as given  in namelist, This may or may not match beam species
!     read from statefile !!
!     Assumes  master list of primary,impurity and posible beams species  in lower case
!  INPUTS (from statefile,and namelist)
!    namep(1:nion)
!    namei(1:nimp)
!    nameb(1:nbion)
!  INPUTS from data statments in modules (ions_gcnmp,neutral_beams)
!    namep_ml(1:nprim_ml)
!    namei_ml(1:nimp_ml)
!    nameb_ml(1:nbeam_ml)
!    atwi_ml(1:nimp_ml)
!    atwp_ml(1:nprim_ml)
!    atomnoi_ml(1:nimp_ml)
!    atomnoi_ml(1:npriim_ml)
!
!  OUTPUTS:
!    atw(1:nion)
!    atmno(1:nion)
!    nameb_index
!    nimp_index
!    ibion
! -- namep_ml(1:nprim_ml),nprim_index,namep are input names of primary ions
! -- namei_ml(1:nimp_ml) , nimp_index,namei are input names of impurity ions
! -- nameb_ml(1:nbeam_ml), nameb_index,nameb is input species  name of beam
! all ions must be elements in the respective master list
! this is checked for in this routine and appropriate indecies,atomc wiehgts
! and atomic numbers are assigned.
!-----------------------------------------------------------------------------------------
       USE nrtype,                    ONLY : DP,I4B

       USE ions_gcnmp,                ONLY : nprim,namep_ml,nprim_ml,namep,name_size,   &
                                             nprim_index,h_index ,t_index,he_index,     &
                                             dt_index,d_index,fi_index,atw,atomno,nion, &
                                             nimp,namei,fd_thermal,namei_ml,nimp_ml,    &
                                             atwp_ml,atomnop_ml,nimp_index,atwi_ml,     &
                                             atomnoi_ml

       USE neutral_beams,             ONLY : nameb,ibion,nameb_ml,nameb_index,nbeam_ml, &
                                             namelist_nameb,nbion,atw_beam
                                            ! NOTE: fd_beam is defined in statfile read
                                            ! fdbeam is namelist input (which is used here)
                                            ! nameb(1:nbion) defined in statfile and 
                                            ! in namelist input

       USE error_handler

       USE io_gcnmp,                   ONLY : nlog,ncrt

       USE MPI_data,                   ONLY : myid,master
       USE common_constants,           ONLY : izero

       IMPLICIT NONE

       INTEGER(I4B) i,j,k,l,ml_beam_index
       CHARACTER(len = name_size) chdum,chdum2




       IF(.NOT. ASSOCIATED(atw))ALLOCATE(atw(nion))
       IF(.NOT. ASSOCIATED(atomno))ALLOCATE(atomno(nion))

!-------------------------------------------------------------------------------------------------
! -- process beam species by discovering index in primary ion master list
! -- Assumes single beam species input in namelist.
!-------------------------------------------------------------------------------------------------
         ml_beam_index = -1 
         nbion = 1
         chdum2 = ADJUSTL(namelist_nameb(nbion))       ! assumes lower case in namelist_nameb,
                                                       ! index of 1 because we allow only 
                                                       ! single  species beam

         ! beam species thermalizes into primary ion species so must be in namep_ml:
 	 DO j=1,nprim_ml
           chdum = ADJUSTL(namep_ml(j))                ! assumes lower case in namep_ml
           IF(TRIM(chdum)  == TRIM(chdum2)) ml_beam_index = j  ! namelist_nameb is in master 
                                                               ! primary ion list
         ENDDO
         IF(ml_beam_index .LT. 1)THEN
            IF(myid == master)     & 
                WRITE(ncrt,FMT='("SUBROUTINE:  check_names, name of beam species ",/, &
                "  not in primary ion master list")')
            lerrno = 214  + iomaxerr
            CALL terminate(lerrno,nlog)
         ENDIF

!-------------------------------------------------------------------------------------------------
! -- process primary ion species by discovering index in primary ion master list
! -- and assigning appropriate atomic and mass numbers
!-------------------------------------------------------------------------------------------------
 

         IF( ASSOCIATED(nprim_index))DEALLOCATE(nprim_index)
         ALLOCATE(nprim_index(nprim))
         IF( ASSOCIATED(fi_index))DEALLOCATE(fi_index)
         ALLOCATE(fi_index(nbion))
         l = izero
         DO k = 1,nprim_ml
            chdum2 = ADJUSTL(namep_ml(k))
            DO j = 1,nprim                 ! find all primary ions in primary ion master list
               chdum = ADJUSTL(namep(j))
               IF(TRIM(chdum) == TRIM(chdum2))THEN ! thermal species index in master list
                  l = l + 1
                  nprim_index(j) = k
               ENDIF
            ENDDO
         ENDDO

         IF(l .ne. nprim)THEN
            IF(myid == master)     & 
                WRITE(ncrt,FMT='("SUBROUTINE:  check_names, name of primary ion species ",/, &
                "  not in primary ion master list")')
            lerrno = 215  + iomaxerr
            CALL terminate(lerrno,nlog)
         ENDIF



         h_index  = izero  ; t_index  = izero ; he_index = izero 
         dt_index = izero  ; d_index  = izero ; fi_index = izero
         l = izero
         DO i=1,nprim
            k = nprim_index(i)           ! k is index in namep_ml
            IF(k == ml_beam_index)THEN   ! both primary ion and beam species point to same element
                                         ! in master primary ion list.
               l = l+1
               fi_index(l) = i
               atw_beam    = atwp_ml(k)
            ENDIF
 
           atomno(i) = atomnop_ml(k)
           atw(i)    = atwp_ml(k)
         ENDDO

         ibion = izero
         ibion = fi_index(1)                        ! ibion could be 'dt' 

!-------------------------------------------------------------------------------------------------
! -- process impurity  ion species in same way as above:
!-------------------------------------------------------------------------------------------------
        IF (nimp .NE. 0)THEN
           IF( ASSOCIATED(nimp_index))DEALLOCATE(nimp_index)
           ALLOCATE(nimp_index(nimp))
           l = izero 
           DO k = 1,nimp_ml
              chdum2 = ADJUSTL(namei_ml(k))
              DO j = 1,nimp                          ! find all impurity  ions in impurity ion master list
                 chdum = ADJUSTL(namei(j))
                 IF(TRIM(chdum) == TRIM(chdum2))THEN ! impurity species index found in master list
                    l = l + 1
                    nimp_index(l) = k
                 ENDIF
              ENDDO
           ENDDO
      
         IF(l .ne. nimp)THEN
            IF(myid == master)     & 
                WRITE(ncrt,FMT='("SUBROUTINE:  check_names, name of impurity ion species ",/, &
                "  not in impurity  ion master list")')
            lerrno = 225  + iomaxerr
            CALL terminate(lerrno,nlog)
         ENDIF


           DO i=1,nimp
              k      = nprim + i
              atw(k) = atwi_ml(nimp_index(i))
              atomno(k) =  atomnoi_ml(nimp_index(i))
              IF (  atw(k) .eq.  0.0)THEN
                 k = nprim -k 
                 IF(myid == master)WRITE(nlog,                                             &
                    FMT='(2x,"Error, no info for impurity,namei =",i5,2x,a)')nimp,namei(k)
                 IF(myid == master)WRITE(ncrt,                                             &
                    FMT='(2x,"Error, no info for impurity,namei =",i5,2x,a)')nimp,namei(k)
                 lerrno  = iomaxerr + 31
                 CALL terminate(lerrno,nlog)
              ENDIF
           ENDDO
        ENDIF


         RETURN
         END  SUBROUTINE set_ion_prop2
