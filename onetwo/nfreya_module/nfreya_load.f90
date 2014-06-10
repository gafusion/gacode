
 SUBROUTINE nfreya_load_single_source(name_beam)
! ------------------------------------------------------------------------
!  setup nfreya type input for beams
!  this data comes from  the nubeam namelist input.
!
!  Construct beams suitable for use in NFreya
!  by analyzing the beam input from  nubeam namelist.
!  Beams have many attributes.
!  info for fast species(injected and fusion products) is here:
!  name_beam%label_thions(:)
!  name_beam%label_fastions(:)
!  name_beam%label_impions(:)
!  name_beam%nfusion_species 
!  name_beam%nbeam_species
!  number of fast species = name_beam%nbeam_species 
!                                   + name_beam%nfusion_specie
!  info for beam geometry:
!  In DIII-D the left source of a beamline is the "tangential"
!  source. It has a tangency radius of 114.6cm
!  The right source of a beamline is the "perpendicular" source
!  it has atangency radius  of 74.94 cm
!  Sfrac1(ib) means the perpendicular source of beam line ib (ib =1,2..kb)
!  The tangential source supplies  1.-sfrac1(ib) fraction of the ions.
!  For P_Nfreya we have up to kb beams each with one  source
!  only 1 beam species is allowed. ( DT mixture is
!  treated as one effective species in Freya.)
!  the sources determine the tangency radius.
! --------------------------------------------------HSJ-4/4/11----------

   USE nrtype,                             ONLY : DP,I4B,SP

   USE nf_param,                           ONLY : kb

   USE P_nfreya_interface,                 ONLY : sort_unique,nptcls

   USE neutral_beams,                      ONLY : nameb,fdbeam,fbcur,ebkev,bptor,      &
                                                  sfrac1,nbeams,beam_on,beam_time,     &
                                                  beam_end,iexcit,izstrp,no_injectors, &
                                                  ebkev,npart,npart_all_beamlines

   USE error_handler,                      ONLY : lerrno,terminate,iomaxerr

   USE io_gcnmp,                           ONLY : nlog

   USE beam_structure,                     ONLY : neutral_beam

   USE common_constants,                   ONLY : izero,zeroc

   IMPLICIT NONE 
 
   TYPE(neutral_beam)name_beam

   INTEGER j,l,k,ks,ufx,ufy,lfx,lfy,uxx,lxx,uyx,             &
           lyx,match,nfbeam,kbspec,ksort,task,error

   LOGICAL  set_diiid

   REAL (DP) rtangcy(name_beam%nbeam)            !temporary array
 



  INTERFACE 
 
     SUBROUTINE set_Pnfreya_beam(nfb,set_diiid,match,name_beam)
       USE nrtype,                   ONLY : DP,I4B,SP
       USE beam_structure,           ONLY : neutral_beam
       IMPLICIT NONE
       TYPE(neutral_beam),INTENT(INOUT)   :: name_beam
       INTEGER,INTENT(in) :: nfb,match
       LOGICAL set_diiid
     END SUBROUTINE set_Pnfreya_beam
   END INTERFACE

   sfrac1(:) = zeroc

 
   IF(npart .GT. nptcls)THEN
     nptcls = npart
   ELSE
     npart = nptcls
   ENDIF
   npart_all_beamlines = npart 

      !decide on unique beam species
       IF(name_beam%nbeam_species .GT. 2)THEN
          PRINT *,'More than 2 injected beam species not allowed' 
          lerrno =298 + iomaxerr
          CALL terminate(lerrno,nlog)
       ELSE IF(name_beam%nbeam_species .GT. 1)THEN
          !max of 2 allowed, must be 'd' , 't'
          kbspec=2
          DO j=1,name_beam%nbeam_species
             match=0
!             DO l=1,kbspec
!                IF(name_beam%label_injected_ions(j) == 'd'      &
!                 .or.  name_beam%label_injected_ions(j) == 't')THEN
!                   match=1
!                   EXIT
!                ENDIF
!              ENDDO
              IF(match == 0) THEN
                   PRINT *,'If two beam species are given'
                   PRINT *,'they must be "d" and "t" '
                   lerrno = 299 + iomaxerr
                   CALL terminate(lerrno,nlog)
              ENDIF
          ENDDO
       ELSE IF(name_beam%nbeam_species == 1)THEN !only single species beam
          kbspec =1 
       ELSE
          PRINT *,'beam species number not set correctly'
          lerrno =298 + iomaxerr
          CALL terminate(lerrno,nlog)
       ENDIF


      !decide on sources based on tangecy radii:
      rtangcy(1:name_beam%nbeam) = name_beam%rtcena(1:name_beam%nbeam)



       CALL sort_unique(rtangcy,ks)               ! There are ks unique tangency radii
                                                 ! but could be co/counter difference
                                                 ! which isnt checked for here?

       PRINT *,' Found ',ks,' unique beamline tangency radii'
 
!      IF(ks .GT. 2)THEN
!         PRINT *,'Not set up to handle more than 2 tangency radii'
!         PRINT *,'Number of sources specified =',ks
!         PRINT *,'with tangency radii :', name_beam%rtcena(:)
!         lerrno = 301 + iomaxerr
!         CALL terminate(lerrno,nlog)
!      ENDIF
 



       ! check the times we allow only one beam on time and one beam off time
       ! for nfreya. OD Fokker Planck time dependent beam model in Onetwo
       ! is not accounted for here at this time.

       beam_on(1:name_beam%nbeam) = name_beam%tbona(:)  !if there is only one beam on time
                                    !then the first element of name_beam%tbona
                                    !will have it

       beam_time(1:name_beam%nbeam) = name_beam%tboffa(:) - name_beam%tbona(1)

       beam_end(1:name_beam%nbeam)  = name_beam%tboffa(:)

       ! OK, we have ks unique sources and kbspec(<=2)
       ! unique beam species, and kbe unique injection
       ! full energies and neutral beam current fractions
       ! and only one beam on time and one beam off time.
       ! info is sufficient to construct Freya input:




      IF(name_beam%nbeam_species  .GT. 1)THEN
         lerrno = 302 + iomaxerr
         CALL terminate(lerrno,nlog)
      ELSE !one beam species.
           !get number of distinct sources 
           !(eg perpendicular and/or parallel sources)
         nbeams =0
         DO j=1,name_beam%nbeam

         IF(name_beam%abeams(1) .LT. 2.)THEN
            nameb='h'
            fdbeam        = 0.150e-3_DP    ! isotopic content of d in h  
         ELSE IF(name_beam%abeams(1) .LT. 3.)THEN
            nameb='d'
            fdbeam = 1.0
         ELSE IF(name_beam%abeams(1) .LT. 4.)THEN
            nameb='t'
            fdbeam =0.0
         ELSE 
            PRINT *,'beam species not recognized'
            PRINT *,'A,Z =',name_beam%abeams(1),name_beam%xzbeams(1)
            lerrno = 303 + iomaxerr
            CALL terminate(lerrno,nlog)
         ENDIF
         ENDDO
      ENDIF


      
      nfbeam =0
      set_diiid = .FALSE.
      no_injectors = name_beam%nbeam 
      DO j=1,name_beam%nbeam                             ! loop over no. of beams 
               match  = izero                            ! assume each beam is unique
               nfbeam = nfbeam+1                         ! requires new freya beam
               IF(nfbeam .LE. kb)THEN

                 CALL set_Pnfreya_beam(nfbeam,set_diiid,match,name_beam)

               ELSE
                 PRINT *,' not enough beams available for P_NfREYA'
                 PRINT *,' max beams allowed =',kb
                 PRINT *,' requested =',nfbeam
                 lerrno = 304 + iomaxerr
                 CALL terminate(lerrno,nlog)
               ENDIF
      ENDDO
      nbeams = nfbeam        !nbeam is final # of Freya beams

!      CALL Check_P_Nfreya_beamlets(name_beam%nbeam)
     

      task = izero
      CALL setup_plot_file(task,error)  ! write name;ist to plot file
      IF(error .NE. izero)THEN
         lerrno = 4_I4B
         CALL terminate(lerrno,nlog)
      ENDIF

       !CALL  set_beam_id               ! nfreya_init.f90 ,no longer used
                                        ! load nlco,beam_id for non ufile
                                        ! type runs
                                        ! moved to nfreya_load_single source

      RETURN



 END SUBROUTINE nfreya_load_single_source




   SUBROUTINE set_Pnfreya_beam(nfb,set_diiid,match,name_beam)
! --------------------------------------------------------------------
! set nfreya beam nfb to attributes of nubeam beam nfb
! if match > 0 then append nfb attributes to nfb.
! if match = 0 then create new nfb entry
! if match = -1 then only differenc is in tangency radius of source
! account for this by adjusting sfrac and appending to nfb
! ---------------------------------------------------------------HSJ-03/24/04
   USE nrtype,                             ONLY : DP,I4B,SP

   USE nf_param,                           ONLY : kb


   USE neutral_beams,                      ONLY : nameb,fdbeam,fbcur,ebkev,bptor, &
                                                  sfrac1,nbeams,beam_on,          &
                                                  beam_time,beam_end,             &
                                                  iexcit,izstrp,anglev,angleh,    &
                                                  nsourc,npart,npskip,relnub,     &
                                                  iterate_beam,naptr,bvofset,     &
                                                  bhofset,bleni,bcur,nbshape,     &
                                                  bheigh,bwidth,bhdiv,bhfoc,bvfoc,&
                                                  nashape,aheigh,awidth,alen,     &
                                                  bvdiv,blenp,rpivot,zpivot,      &
                                                  beam_sim_time_start,            &
                                                  beam_sim_time_end,no_injectors, &
                                                  no_physical_injectors

   USE error_handler,                      ONLY : lerrno,terminate,iomaxerr

   USE io_gcnmp,                           ONLY : nlog

   USE zonal_data,                         ONLY : mf,mfm1

   USE beam_structure,                     ONLY : neutral_beam

   USE common_constants,                   ONLY : izero,zeroc

   USE Plasma_properties ,                 ONLY : neut_beam

   IMPLICIT NONE

   TYPE(neutral_beam),INTENT(INOUT) :: name_beam
   INTEGER,INTENT(in) :: nfb,match
   INTEGER i,ke,ks,npts,j,kt
   REAL(DP) pwravg,einjavg,fullavg,halfavg,stdpwravg,stdeinjavg, &
            stdfullavg,stdhalfavg,dt_sum,dt_incp,dt_incm,dt_inc, &
            dummy
   LOGICAL set_diiid
   


      !---------------------------------------------------------------------
      ! set data from nubeam namelist input:
      !--------------------------------------------------------------------

             mf = name_beam%nzone_nb
             mfm1 = mf-1

             beam_on(1:name_beam%nbeam)   = name_beam%tbona(:)
             beam_time(1:name_beam%nbeam) = name_beam%tboffa(:) - name_beam%tbona(:)
             beam_end(1:name_beam%nbeam)  = name_beam%tboffa(:)



      !---------------------------------------------------------------------
      !  Set some default data:
      !---------------------------------------------------------------------
          nsourc= 2 


      !---------------------------------------------------------------------
      !  Remainder of namelist input set from P_nfreya run directives file
      !----------------------------------------------------------------------
!       npart = 50000 get this from  P_nfreya run directives namelist input
!       npskip=5 
!       iexcit = 5
        izstrp(1)  = 1 ; izstrp(2) =1
        relnub = 0.01
        iterate_beam=.TRUE.      !dont do this if beam is initially off
  

      !-----------------------------------------------------------------------
      ! Set power,energy,full,half,third current fractions.
      !  name_beam%beam_inject,name_beam%beam_times,name_beam%beam_chan
      ! were set in set_beam_data if ufile was read. Otherwise these
      ! items were set in ufile_bypass.
      !-----------------------------------------------------------------------
! print *,'nfreya_load line 301:'  ! 8888899999
! print *,'name_beam%beam_inject =',name_beam%beam_inject
! print *,'name_beam%beam_times =',name_beam%beam_times
! print *,'name_beam%beam_chan =',name_beam%beam_chan
! print *,'name_beam%tbona =',name_beam%tbona
! print *,'name_beam%tbona =',name_beam%tboffa
! print *,'name_beam%beam_times=',name_beam%beam_times

      ! 1) select a time index range  based on current time interval
           ! (beam_sim_time_start ,beam_sim_time_end )
           !-----------------------------------------

           kt = SIZE(name_beam%beam_times)
           IF(kt .LT. 2)THEN
              lerrno = 325 + iomaxerr
              CALL terminate(lerrno,nlog)
           ENDIF
           ke = izero  ; ks = izero
           DO j=1,kt
              IF(name_beam%beam_times(j) .LE. beam_sim_time_start)ks =j
              IF(name_beam%beam_times(j) .LE. beam_sim_time_end)  ke =j
           ENDDO
           ks  = MAX(ks,2) ;   ! ke = MIN(ke,SIZE(name_beam%beam_times))
           ke = Max(ke,ks)
           npts = ke-ks+1                         ! at least one point required
           IF(npts .LT. 1 )THEN
                 PRINT *,"ERROR beam simualtion times must be set "
                 PRINT *,"  so that at least one valid data point can be found"
                 PRINT *,"  beam_sim_time_end = ",beam_sim_time_end
                 PRINT *,"  beam_sim_time_start = ",beam_sim_time_start
                 lerrno = iomaxerr + 314_I4B
                 CALL terminate(lerrno,nlog)
           ENDIF

              pwravg       = zeroc
              einjavg      = zeroc 
              fullavg      = zeroc
              halfavg      = zeroc 
              stdpwravg    = zeroc
              stdeinjavg   = zeroc
              stdfullavg   = zeroc
              stdhalfavg   = zeroc
              ! here we assume that dt may be nonuniform
              dt_sum       = zeroc


              DO j=ks,ke                 ! makes at least one pass (if ks =ke)

                 dt_incp  = name_beam%beam_times(j+1) - name_beam%beam_times(j)
                 dt_incm  = name_beam%beam_times(j) - name_beam%beam_times(j-1)
                 dt_inc = 0.5_DP*(dt_incp+dt_incm)
                 dt_sum  = dt_sum  + dt_inc
                 pwravg  = pwravg  + name_beam%beam_inject(j,nfb)*dt_inc
                 stdpwravg  = stdpwravg  + (name_beam%beam_inject(j,nfb)**2)*dt_inc

                 einjavg = einjavg + name_beam%beam_inject(j,no_injectors + nfb)*dt_inc
                 stdeinjavg = stdeinjavg + (name_beam%beam_inject(j,no_injectors + nfb)**2)*dt_inc

                 fullavg = fullavg + name_beam%beam_inject(j,2*no_injectors + nfb)*dt_inc
                 stdfullavg = stdfullavg + (name_beam%beam_inject(j,2*no_injectors + nfb)**2)*dt_inc

                 halfavg = halfavg + name_beam%beam_inject(j,3*no_injectors + nfb)*dt_inc
                 stdhalfavg = stdhalfavg + (name_beam%beam_inject(j,3*no_injectors + nfb)**2)*dt_inc

              ENDDO

                 pwravg     = pwravg/dt_sum
                 einjavg    = einjavg/dt_sum
                 fullavg    = fullavg/dt_sum
                 halfavg    = halfavg/dt_sum                
                 IF(npts .gt.1)THEN
                    dummy = MAX(stdpwravg-dt_sum*pwravg**2,zeroc)
                    stdpwravg     = SQRT(dummy/dt_sum)
                    dummy = MAX(stdeinjavg-dt_sum*einjavg**2,zeroc)
                    stdeinjavg    = SQRT(dummy/dt_sum)
                    dummy = MAX(stdfullavg-dt_sum*stdfullavg**2,zeroc)
                    stdfullavg    = SQRT(dummy/dt_sum)
                    dummy = MAX(stdhalfavg-dt_sum*stdhalfavg**2,zeroc)
                    stdhalfavg    = SQRT(dummy/dt_sum)
                 ELSE

                    stdpwravg     = zeroc
                    stdeinjavg    = zeroc
                    stdfullavg    = zeroc
                    stdhalfavg    = zeroc   
                 ENDIF


      ! 2) store the data at this time:
              name_beam%pinja(nfb)   = pwravg  ! select average power over time interval
              name_beam%einja(nfb)   = einjavg
              name_beam%ffulla(nfb)  = fullavg
              name_beam%fhalfa(nfb)  = halfavg
         
      bptor(nfb)     =  name_beam%pinja(nfb)
      ebkev(nfb)     =  name_beam%einja(nfb)/1000.
      fbcur(1,nfb)   =  name_beam%ffulla(nfb)
      fbcur(2,nfb)   =  name_beam%fhalfa(nfb)
      fbcur(3,nfb)   =  1._DP  -fbcur(1,nfb) -fbcur(2,nfb)
      fbcur(3,nfb)   =  MAX(0.0D0,fbcur(3,nfb)) !jmp.ibm

! print *,'ASSOCIATED neut_beam%fbcur =',ASSOCIATED(neut_beam%fbcur) ! 88889999
! print *,'no_physical_injectors =',no_physical_injectors,no_injectors
      IF(.NOT. ASSOCIATED(neut_beam%fbcur))ALLOCATE(neut_beam%fbcur(3,no_injectors))
      neut_beam%fbcur(1,nfb) = fbcur(1,nfb)  ! name_beam%fbcur for output to statefile
      neut_beam%fbcur(2,nfb) = fbcur(2,nfb)
      neut_beam%fbcur(3,nfb) = fbcur(3,nfb)


      CALL set_P_Nfreya_geometry(nfb,name_beam)        

    RETURN
   END    SUBROUTINE set_Pnfreya_beam



    SUBROUTINE check_beam_attrib(k,j,match,name_beam)
! --------------------------------------------------------
! check attributes of beam j against those of beam k
! set match =1 IF perfect match
! set match =-1   IF ONLY difference is in
! tangency radius of source for the two beams 
! otherwise RETURN match = 0
! ------------------------------------------------------------
    USE nrtype,                              ONLY : DP,I4B,SP

    USE beam_structure,                       ONLY : neutral_beam
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k,j
    REAL(DP) ::  dtmatch = 0.001


    INTEGER match
    TYPE(neutral_beam)name_beam


    match = 1                                         !assume perfect match

    IF(ABS(name_beam%rtcena(j) -  name_beam%rtcena(k)) .GT. 0.01) match = -1


    IF(ABS(name_beam%tbona(j) - name_beam%tbona(k)) .GT. dtmatch) match =0
    IF(ABS(name_beam%tboffa(j) - name_beam%tboffa(k)) .GT. dtmatch) match =0


    IF(ABS(name_beam%abeama(j) -  name_beam%abeama(k)) .GT. 0.01)match =0
    IF(ABS(name_beam%xzbeama(j) -  name_beam%xzbeama(k)) .GT. 0.01)match =0



!!  IF(ABS(name_beam%aplasm(j) -  name_beam%aplasm(k)) .GT. 0.01)match =0
!!  IF(ABS(name_beam%backz(j) -  name_beam%backz(k)) .GT. 0.01)match =0
!!  IF(ABS(name_beam%xzimpx(j) -  name_beam%xzimpx(k)) .GT. 0.01)match = 0
!!  IF(ABS(name_beam%aimpx(j) -  name_beam%aimpx(k)) .GT. 0.01)match = 0


    IF(ABS(name_beam%bmwidra(j) -  name_beam%bmwidra(k)) .GT. 0.01)match =0
    IF(ABS(name_beam%bmwidza(j) -  name_beam%bmwidza(k)) .GT. 0.01)match =0

    IF(ABS(name_beam%xlbtna(j) -  name_beam%xlbtna(k)) .GT. 0.01)match = 0
    IF(ABS(name_beam%xybsca(j) -  name_beam%xybsca(k)) .GT. 0.01)match = 0

    IF(ABS(name_beam%divza(j) -  name_beam%divza(k)) .GT. 0.01)match = 0
    IF(ABS(name_beam%divra(j) -  name_beam%divra(k)) .GT. 0.01)match = 0


    IF(ABS(name_beam%foclza(j) -  name_beam%foclza(k)) .GT. 0.01)match = 0
    IF(ABS(name_beam%foclra(j) -  name_beam%foclra(k)) .GT. 0.01)match = 0


    IF(ABS(name_beam%rapedga(j) -  name_beam%rapedga(k)) .GT. 0.01)match = 0
    IF(ABS(name_beam%xzpedga(j) -  name_beam%xzpedga(k)) .GT. 0.01)match = 0


    IF(ABS(name_beam%xlbapa2(j) -  name_beam%xlbapa2(k)) .GT. 0.01)match = 0
    IF(ABS(name_beam%xlbapa(j) -  name_beam%xlbapa(k)) .GT. 0.01)match = 0


    IF(ABS(name_beam%xybapa(j) -  name_beam%xybapa(k)) .GT. 0.01)match = 0
    IF(ABS(name_beam%rapedg2(j) -  name_beam%rapedg2(k)) .GT. 0.01)match = 0

    IF(ABS(name_beam%xzpedg2(j) -  name_beam%xzpedg2(k)) .GT. 0.01)match = 0



    IF(ABS(name_beam%xbzeta(j) -  name_beam%xbzeta(k)) .GT. 0.01)match = 0

    IF(ABS(name_beam%pinja(j) -  name_beam%pinja(k)) .GT. 0.01)match = 0
    IF(ABS(name_beam%einja(j) -  name_beam%einja(k)) .GT. 0.01)match = 0

    IF(ABS(name_beam%ffulla(j) -  name_beam%ffulla(k)) .GT. 0.01)match = 0
    IF(ABS(name_beam%fhalfa(j) -  name_beam%fhalfa(k)) .GT. 0.01)match = 0


    IF(ABS(name_beam%nbshapa(j) -  name_beam%nbshapa(k)) .GT. 0.01)match = 0
    IF(ABS(name_beam%nbapsha(j) -  name_beam%nbapsha(k)) .GT. 0.01)match = 0


!    IF(name_beam%nlco(j) .NE.  name_beam%nlco(k))match = 0
    IF(name_beam%nlco(j))THEN
       IF(.NOT. name_beam%nlco(k))match =0
    ELSEIF( .NOT. name_beam%nlco(j))THEN
       IF(name_beam%nlco(k))match =0
    ENDIF
    IF(ABS(name_beam%nbapsha(j) -  name_beam%nbapsha(k)) .GT. 0.01)match = 0

    END SUBROUTINE check_beam_attrib




   SUBROUTINE nbdrive_nubeam_checknaml(ierr,name_beam)
!--------------------------------------------------------------------------
!  perform additional namelist checks.  Call a NUBEAM routine to
!  pre-fetch the fast ion distribution function grids.
!  The following checks are performed:
!
!  a) no duplicates in thermal species lists
!  b) for each beam injected specie a matching thermal specie must
!     be present
!  c) for each fusion product specie a matching thermal specie is
!     not required, but a warning will be issued.
!
!  call ancillary NUBEAM routines to check ion species lists; return 
!  indices to specific isotopes, total no. of fast ion species, error 
!  flags if duplicate species are given, etc.
!  the NUBEAM calls in this routine do not affect the NUBEAM module's
!  internal state-- no NUBEAM f90 modules or COMMONs are used-- 
!  and they can be called at any time.  They provide information
!  which is available both to NUBEAM and to NUBEAM driver codes, 
!  and which depends only on the routines' input arguments.
!--------------------------------------------------------------------------

  USE nrtype,                         ONLY : DP,I4B,SP

  USE P_nfreya_interface,             ONLY : nlfst,nlfhe3,nlfhe4, &
                                              nlfsp,nlusf3,sort_unique

  USE io_gcnmp,                       ONLY : lun_msgs => ncrt

   USE error_handler,                 ONLY : lerrno,terminate,iomaxerr

   USE io_gcnmp,                      ONLY : nlog

   USE beam_structure,ONLY : neutral_beam

  IMPLICIT NONE

  TYPE(Neutral_beam)name_beam

  REAL(DP), PARAMETER :: CZERO = 0.0d0

  INTEGER, INTENT(out) :: ierr

 
  INTEGER  iwarn, ierr_p, i, j, iA, iZ, iA2, iZ2, imatch, ilin,ngmax,nsfast, &
           nuqa,nuqm,nuqz
  REAL(DP) zelin,zeinjmax

  REAL(DP) List_abeam(name_beam%nbeam)
  REAL(DP) List_xzbeam(name_beam%nbeam)

 


  ierr = 0
  ngmax = name_beam%ngmax

  !-------------------------------------------------------------
  !  check thermal ion species list...
  !  input: species list; output: indices, error/warning flags
  !--------------------------------------------------------------
  DO i=1,ngmax-1
     iA =  name_beam%aplasm(i) + 0.1d0
     iZ =  name_beam%backz(i) + 0.1d0
     DO j=i+1,ngmax
        iA2 = name_beam%aplasm(j) + 0.1d0
        iZ2 = name_beam%backz(j) + 0.1d0
        IF((iA.EQ.iA2).AND.(iZ.EQ.iZ2)) ierr=ierr+1
     ENDDO
  ENDDO

  IF(ierr.NE.0) THEN
     WRITE(lun_msgs,*) '  duplicate species in thermal species list:'
     WRITE(lun_msgs,*) '  name_beam%aplasm(1:ngmax) = ',name_beam%aplasm(1:ngmax)
     WRITE(lun_msgs,*) '   name_beam%backZ(1:ngmax) = ', name_beam%backZ(1:ngmax)
     lerrno = 306 + iomaxerr
     CALL terminate(lerrno,nlog)

  ENDIF

  !  do same for impurities

  DO i=1,name_beam%nrhix-1
     iA = name_beam%aimpx(i)
     iZ = name_beam%xzimpx(i)
     DO j=i+1,name_beam%nrhix
        iA2 = name_beam%aimpx(j)
        iZ2 = name_beam%xzimpx(j)
        IF((iA.EQ.iA2).AND.(iZ.EQ.iZ2)) ierr=ierr+1
     ENDDO
  ENDDO

  IF(ierr.NE.0) THEN
     WRITE(lun_msgs,*) ' ?duplicate species in impurity species list:'
     WRITE(lun_msgs,*) '  Aimpx(1:name_beam%nrhix) = ', name_beam%Aimpx(1:name_beam%nrhix)
     WRITE(lun_msgs,*) '  xZimpx(1:name_beam%nrhix) = ',name_beam%xZimpx(1:name_beam%nrhix)
     lerrno = 307 + iomaxerr
     CALL terminate(lerrno,nlog)
  ENDIF


  !------------------------------------------------------------------
  !  build list of fast species from beam species...
  !  NOTE name_beam%abeams and name_beam%xzbeams are not passed to nubeam
  !  Instead name_beam%abeama and name_beam%xzbeama are passed
  !  This is due to the fact that nubeam constructs abeams and xzbeams
  !  Internally. We keep the local results in name_beam so that we may compare 
  !  with the output  of nubeam in sub load_12.
  !-------------------------------------------------------------------
   List_abeam(1:name_beam%nbeam)  = name_beam%abeama(1:name_beam%nbeam)
   List_xzbeam(1:name_beam%nbeam) = name_beam%xzbeama(1:name_beam%nbeam)

   nuqa = name_beam%nbeam         ; nuqz = name_beam%nbeam
   CALL sort_unique( List_abeam,nuqa)
   CALL sort_unique( List_xzbeam,nuqz)
   nuqm = MAX(nuqa,nuqz)
 


  nsfast = nuqm
  IF(nsfast == 1)THEN
     name_beam%abeams(:)   = LIST_abeam(1)
     name_beam%xzbeams(:)  = LIST_xzbeam(1)
  ELSEIF(nsfast == 2)THEN

  ELSE

  ENDIF

!  GO TO 100
  nsfast=0


!  DO i=1,name_beam%nbeam              ! each beam may be a different species
!     iA=name_beam%abeama(i)           ! abeama,xzbeama are read from nubeam namelist file
!     iZ=name_beam%xzbeama(i)          ! in sub nblist_read
!     imatch=0

     !DO j=1,SIZE(name_beam%xzbeams)   ! Z of beam species,xzbeams,abeams are from ?

!        iA2=INT(name_beam%abeams(j))
!        iZ2=INT(name_beam%xzbeams(j))

!        IF(iZ2.GT.2) THEN
!           WRITE(lun_msgs,*) ' ?nbdrive: beams must be isotopes of H or He'
!           ierr = ierr+1
!           imatch = 99
!           EXIT
!        ELSE IF((iA.EQ.iA2).AND.(iZ.EQ.iZ2)) THEN
!           imatch=j
!           EXIT
!        ENDIF
!     ENDDO

!100 CONTINUE



    ! IF(imatch.EQ.0) THEN

        nsfast = nsfast + 1
        name_beam%abeams(nsfast)  = name_beam%abeama(i)
        name_beam%xzbeams(nsfast) = name_beam%xzbeama(i)

        !  require that corresponding thermal specie exist

        imatch = 0
        DO j=1,ngmax
           iA2 = name_beam%aplasm(j) + 0.1d0
           iZ2 = name_beam%backz(j) + 0.1d0
           IF((iA.EQ.iA2).AND.(iZ.EQ.iZ2)) THEN
              imatch=j
              EXIT
           ENDIF
        ENDDO

        IF(imatch.EQ.0) THEN
           WRITE(lun_msgs,*) ' ?nbdrive_nubeam_checknaml: beam no. ',i
           WRITE(lun_msgs,*) '  injects specie A=',name_beam%abeama(i),' Z=',name_beam%xzbeama(i)
           WRITE(lun_msgs,*) '  but there is no such thermal specie:'
           WRITE(lun_msgs,*) '  Aplasm(1:ngmax) = ',name_beam%aplasm(1:ngmax)
           WRITE(lun_msgs,*) '   backZ(1:ngmax) = ', name_beam%backZ(1:ngmax)
           ierr = ierr+1
        ENDIF
     !ENDIF
!  ENDDO







  IF(ierr.GT.0) THEN
     WRITE(lun_msgs,*) &
          ' NUBEAM nbi_datchkb detected problem with beam species.'
      lerrno = 308 + iomaxerr
      CALL terminate(lerrno,nlog)
  ENDIF




   name_beam%nbeam_species = nsfast     ! number of unique  injected beam species
                                       ! regardless of beam number

!------------------------------------------------------------------------
! --  check fusion product species
! --  also compute array of maximum energies for distribution functions.
!-------------------------------------------------------------------------
  zeinjmax = MAXVAL(name_beam%einja(1:MAX(1,name_beam%nbeam)))
  name_beam%ebdmaxa = CZERO
  IF(name_beam%ebdmax .LE. 0.0)name_beam%ebdmax = 1.1*zeinjmax
  name_beam%ebdmaxa(1:nsfast) = name_beam%ebdmax   ! max energy for BEAM species (eV) from namelist
   name_beam%nznbmea = 0
   name_beam%nznbmea(1:nsfast) = name_beam%nznbme   ! no. of energy zones for BEAM species dist fcns

    name_beam%nlfprod(:) = 0            ! fusion product species flags

  !  input: flags for specific fusion product species
  !         information arrays for beam species
  ! output: information arrays extended to cover fusion product species
  !         indices for specific fusion product isotopes


  !  note setting of maximum energy for fusion product distribution
  !  functions (eV) using NUBEAM utility call "r8_nbfusn_emax(iZ,iA,emax)"

  IF(nlfst) THEN
     nsfast = nsfast + 1
     name_beam%xzbeams(nsfast)=1
     name_beam%abeams(nsfast)=3
     CALL r8_nbfusn_emax(1,3,name_beam%ebdmaxa(nsfast))
                        ! safely above the 1.01 MeV fusion T birth energy
      name_beam%nznbmea(nsfast)=name_beam%nznbme
      name_beam%nlfprod(nsfast)= 1
  ENDIF

  IF(nlfhe3) THEN
     nsfast = nsfast + 1
     name_beam%xzbeams(nsfast)=2
     name_beam%abeams(nsfast)=3
     CALL r8_nbfusn_emax(2,3,name_beam%ebdmaxa(nsfast))
                        ! safely above the 0.82 MeV fusion He3 birth energy
      name_beam%nznbmea(nsfast)=name_beam%nznbme
      name_beam%nlfprod(nsfast)= 1
  ENDIF

  IF(nlfhe4) THEN
     nsfast = nsfast + 1
     name_beam%xzbeams(nsfast)=2
     name_beam%abeams(nsfast)=4
     CALL r8_nbfusn_emax(2,4,name_beam%ebdmaxa(nsfast))
                        ! safely above the alpha (He4) birth energy
                        ! for any fusion reaction
     name_beam%nznbmea(nsfast)=name_beam%nznbme
     name_beam%nlfprod(nsfast)= 1
  ENDIF

   name_beam%nfusion_species = nsfast - name_beam%nbeam_species




  !-------------------------------
  !  there is also an r8_nbi_datckrf call for RF minorities -- not
  !    currently used in nbdrive.
  !-------------------------------
  !  allocate energy grids for the fast ion distributions; call a NUBEAM
  !  routine to fill these in.  Although not necessary to make NUBEAM work,
  !  this allows the NUBEAM caller to see the distribution function grids
  !  ahead of time.
  !    note that the above calls to "r8_nbfusn_emax(iZ,iA,emax)" to set
  !    the max energy.

  !  in nbdrive, all species distribution functions have the same number
  !  of energy grid points, but different maximum energies.  In general,
  !  it is allowed for the number of energy grid points to vary depending
  !  on fast specie, however.

  ALLOCATE(name_beam%efbm(name_beam%nznbme,nsfast),name_beam%efbmb(name_beam%nznbme+1,nsfast))
  name_beam%efbm=CZERO; name_beam%efbmb=CZERO
  DO i=1,nsfast
     CALL r8_nbi_efbm(name_beam%nlfprod(i), name_beam%nznbmea(i),name_beam%ebdmaxa(i), &   ! ** NUBEAM CALL **
          ilin,zelin,name_beam%efbmb(1: name_beam%nznbmea(i)+1,i), &
          name_beam%efbm(1: name_beam%nznbmea(i),i))
     ! (r8_)nbi_efbm outputs:
     ! efbm -- grid zone centers
     ! efbmb -- grid boundaries
     ! ilen,zelin -- index and energy to which grid is linearly (rather
     ! than logarithmically) spaced.
  ENDDO



  !-------------------------------
  !  the following is for the benefit of nbdrive output routines.
  !  it is placed here because the fast ion species list has now
  !  been defined.
  !-------------------------------
  !  now gen some labels for species
  
  IF(ierr.EQ.0) THEN
     ALLOCATE(name_beam%label_thions(name_beam%ngmax))
     DO i=1,name_beam%ngmax
        CALL zlabel( name_beam%backz(i), name_beam%aplasm(i),'thermal plasma ion',name_beam%label_thions(i))
     ENDDO

     ALLOCATE(name_beam%label_impions(name_beam%nrhix))
     DO i=1,name_beam%nrhix
        CALL zlabel( name_beam%xzimpx(i), name_beam%aimpx(i),'impurity plasma ion',name_beam%label_impions(i))
     ENDDO

     ALLOCATE(name_beam%label_fastions(nsfast))
     DO i=1,nsfast
        IF(name_beam%nlfprod(i) .gt. 0 ) THEN
           CALL zlabel( name_beam%xzbeams(i), name_beam%abeams(i),'fusion product ion', &
                name_beam%label_fastions(i))
        ELSE
           CALL zlabel( name_beam%xzbeams(i), name_beam%abeams(i),'beam ion', &
                name_beam%label_fastions(i))
        ENDIF
     ENDDO

     ALLOCATE(name_beam%label_beamions(name_beam%nbeam_species))
     i=0
     DO j=1,nsfast
         IF(name_beam%nlfprod(j) .le. 0 )THEN
           i=i+1
           IF(i .gt. name_beam%nbeam_species) THEN     
              lerrno = 309 + iomaxerr
              CALL terminate(lerrno,nlog)
           ENDIF
           IF(name_beam%abeams(j) == 1.0 )THEN
              name_beam%label_beamions(i) = 'h' !note lower case for onetwo
           ELSE IF(name_beam%abeams(j) == 2.0)THEN
              name_beam%label_beamions(i) = 'd'
           ELSE
              lerrno = 310 + iomaxerr
              CALL terminate(lerrno,nlog)
           ENDIF
        ENDIF
     ENDDO
  ENDIF

  IF(ierr.EQ.0) THEN
     WRITE(lun_msgs,*) ' ---> nbdrive_nubeam_checknaml: success.'
  ENDIF

  CONTAINS
    SUBROUTINE zlabel(zz,aa,suffix,label)

      USE periodic_table_mod

      !  generate labels for ion species

      REAL(DP), INTENT(in) :: zz,aa   ! Z & A of species
      CHARACTER*(*), INTENT(in) :: suffix   ! label suffix
      CHARACTER*(*), INTENT(out) :: label   ! generated label.

      REAL(SP) aaa
      INTEGER iz,ia

      CHARACTER*12 prefix

      !----------------

      aaa = SNGL(aa)
      iz = zz + 0.1
      ia = aa + 0.1

      prefix='?'
      IF(iz.EQ.1) THEN
         IF(ia.EQ.1) THEN
            prefix='H'
         ELSE IF(ia.EQ.2) THEN
            prefix='D'
         ELSE IF(ia.EQ.3) THEN
            prefix='T'
         ENDIF
      ELSE IF(iz.EQ.2) THEN
         IF(ia.EQ.3) THEN
            prefix='He3'
         ELSE IF(ia.EQ.4) THEN
            prefix='He4'
         ENDIF
      ELSE

         prefix = to_periodic_table(iz,aaa,-1,0)

      ENDIF

      label = TRIM(prefix)//' '//TRIM(suffix)

    END SUBROUTINE zlabel

END SUBROUTINE nbdrive_nubeam_checknaml

 SUBROUTINE nfrpar_from_nubpar(xybsca,             &
                               xybapa,             &
                               xlbapa,             &
                               xlbtna,             &
                               rtcena,             &
                               rpivot,             &
                               anglev,             &
                               angleh,             &
                               zpivot,             &
                               blenp               )
!---------------------------------------------------------------------
! -- This subroutine (from R.Prater) maps nubeam parameters to
! -- Nfreya beam parameters
! -- INPUTS
!        xybsca             
!        xybapa             
!        xlbapa            
!        xlbtna            
!        rtcena             
!        rpivot  <==NOTE
!
! -- OUTPUTS
!        anglev             
!        angleh             
!        zpivot             
!        blenp               
!---------------------------------------------------------------------
   USE nrtype,                             ONLY : DP,I4B,SP

   USE common_constants,                   ONLY : Deg_Per_Rad 
  
   IMPLICIT NONE

   REAL(DP) xybsca,xybapa,xlbapa,xlbtna,rtcena

   REAL(DP) angleh,anglev,rpivot,zpivot,blenp
 
   REAL(DP) ysource,ypivot,a,b,c

   go to 10
   anglev  = -ASIN((xybsca-xybapa)/xlbapa)  ! ADDED - SIGN hsj 
   angleh  = ASIN(rtcena/rpivot)
   ysource = xlbtna*COS(anglev) 
   ypivot  = rpivot*COS(angleh) 
   blenp   = (ysource -ypivot)/COS(anglev)
!   zpivot  = xybsca -blenp*SIN(anglev)
   zpivot  = xybsca +blenp*SIN(anglev)     !hsj 3/7/2012

!-----------------new (Prater)--------3/30/20012---------
10 CONTINUE
   ANGLEV=-ASIN((XYBSCA-XYBAPA)/XLBAPA)
   ANGLEH=ASIN(RTCENA/RPIVOT)
   ! parameterize the beam line. Parameter t=0 is the source location.
   ! Find t where x2+y2=RPIVOT2
   ! x=rtcena
   ! y=-xlbtna*cos(anglev)+cos(anglev)*t
   ! z=xybsca-sin(anglev)*t
   C= (RTCENA**2 - RPIVOT**2)/cos(anglev)**2  + XLBTNA**2
   B=-2._DP*XLBTNA
   A=1.0_DP
   BLENP=(-B-SQRT(B**2 -4._DP*A*C))/(2._DP*A)	! take nearer solution
   ZPIVOT=XYBSCA+BLENP*SIN(ANGLEV)

   anglev  = Deg_Per_Rad*anglev
   angleh  = Deg_Per_Rad*angleh


   RETURN

 END SUBROUTINE nfrpar_from_nubpar


 SUBROUTINE set_P_Nfreya_geometry(nfb,name_beam)
 !-------------------------------------------------------------------------
 ! Setup  for running P_nfreya  with single source per beamline:
 ! INPUT:
 !    nfb         beamlet no
 !                We assume that the beamlets are in sequence
 !            1   30LT
 !            2   30RT
 !            3   150LT UP
 !            4   150LT MU
 !            5   150LT ML
 !            6   150LT LO
 !            7   150RT UP
 !            8   150RT MU
 !            9   150RT ML
 !           10   150RT LO
 !           11    21LT
 !           12    21RT
 !           13    33LT
 !           14    33RT
 !-------------------------------------------------------------------------

   USE nrtype,                             ONLY : DP,I4B,SP

   USE error_handler,                      ONLY : lerrno,terminate,iomaxerr

   USE io_gcnmp,                           ONLY : nlog

   USE beam_structure,                     ONLY : neutral_beam

   USE common_constants,                   ONLY : izero,zeroc,Deg_Per_Rad
 
   USE neutral_beams,                      ONLY : anglev, angleh,nashape,aheigh,  &
                                                  awidth, bcur, bptor, blenp,     &
                                                  nbshape,bleni,bheigh, bwidth,   &
                                                  bhfoc, bvfoc, bhdiv, bvdiv,     &
                                                  ebkev, fbcur, nbeams,naptr,     &
                                                  alen, bvofset, bhofset, nsourc, &
                                                  sfrac1, rpivot,zpivot,nap


   IMPLICIT NONE


   INTEGER(I4B) nfb,iap


   TYPE(neutral_beam)name_beam

 


!------------------------------------------------------------------------------! 
! FOLLOWING GEOMETRY ATTRIBUTES HAVE TO BE SET FROM nubeam input.
!      DEFAULT VALUES for Nfreya input were set in nfreya_namelist_defaults.f90
!      and are stored in module neutral_beams:
!      DIII-D beam input.  DEFAULT IS LONG-PULSE-SOURCE SPECIFICATIONS
!      naptr = 4                                                          
!     DO i=1,kb                      
!        anglev(i)    =   0.0                                                  
!        angleh(i)    =  19.5
!        bvofset(i)   =   0.0                                     
!        bhofset(i)   =  42.074                                                
!        bleni(i)     = 556.808 
!        bcur(i)      = 110.0 
!        bptor(i)     =  0.0e6                                                 
!        nbshape(i)   = 'rect-lps' 
!        bheigh(i)    =  48.0                                                  
!        bwidth(i)    =  12.0                                                  
!        bhdiv(i)     =   0.50       !degrees                                  
!        bvdiv(i)     =   1.3        !degrees 
!        fbcur(1,i)   =   0.7                                                  
!        fbcur(2,i)   =   0.2                                                  
!        fbcur(3,i)   =   0.1                                                  
!        bhfoc(i)     =   1.0d100  
!        bvfoc(i)     =   1.0d3                                                
!        ebkev(i)     =  75.0 
!        sfrac1(i)    =   0.5                                                  
!        nashape(1,i) = 's-rect' 
!        nashape(2,i) = 's-rect' 
!        nashape(3,i) = 'b-d3d' 
!        nashape(4,i) = 'b-circ' 
!        aheigh(1,i)  =  47.8
!        awidth(1,i)  =  13.8
!        alen(1,i)    = 186.1 
!        aheigh(2,i)  =  48.0 
!        awidth(2,i)  =  17.7
!        alen(2,i)    = 346.0 
!        alen(3,i)    = 449.0 
!        awidth(4,i)  =  50.9  
!        alen(4,i)    = 500.0 
!        blenp(i)     = 539.0                                                  
!        rpivot(i)    = 286.6 
!        zpivot(i)    =   0.0                                                  
!      END DO  
!----------------------------------------------------------------------------                                                  

! P_Nfreya always uses 2 sources one of which is turned off to
! get the single source per beamline. (Using nsourc =1 would require
! change logic in Nfreya to account for beamline ctr versus
! offset source locations HSJ
           nsourc = 2
           SELECT CASE(nfb)               ! fixed input sequence assumed
              CASE (1)                    ! 30 LT
                 sfrac1(nfb) = zeroc
                 name_beam%beam_id(nfb)  =  "30LT"
              CASE(2)                     ! 30RT
                 sfrac1(nfb) = 1._DP
                 name_beam%beam_id(nfb)  =  "30RT"
              CASE(3)                     ! 150LT UP  ! need to verify these
                 sfrac1(nfb) = 0._DP      
                 name_beam%beam_id(nfb)  =  "150LT UP"
              CASE(4)                     ! 150LT MU
                 sfrac1(nfb) = 0._DP
                 name_beam%beam_id(nfb)  =  "150LT MU"
              CASE(5)                     ! 150LT ML
                 sfrac1(nfb) = 0._DP
                 name_beam%beam_id(nfb)  =  "150LT ML"
              CASE(6)                     ! 150LT LO
                 sfrac1(nfb) = 0._DP
                 name_beam%beam_id(nfb)  =  "150LT LO"
              CASE(7)                     ! 150RT UP
                 sfrac1(nfb) = 1._DP
                 name_beam%beam_id(nfb)  =  "150RT UP"  
              CASE(8)                     ! 150RT MU
                 sfrac1(nfb) = 1._DP
                 name_beam%beam_id(nfb)  =  "150RT MU" 
              CASE(9)                     ! 150RT ML
                 sfrac1(nfb) = 1._DP
                 name_beam%beam_id(nfb)  =  "150RT ML"  
              CASE(10)                    ! 150RT LO
                 sfrac1(nfb) = 1._DP        
                 name_beam%beam_id(nfb)  =  "150RT LO"      
              CASE(11)                    ! 21LT
                sfrac1(nfb) = zeroc     
                name_beam%beam_id(nfb)  =  "21LT" 
              CASE(12)                    ! 21RT
                sfrac1(nfb) = 1._DP      
                name_beam%beam_id(nfb)  =  "21RT"
              CASE(13)                    ! 33LT
                sfrac1(nfb) = zeroc      
                name_beam%beam_id(nfb)  =  "33LT"
              CASE(14)                    ! 33RT
                sfrac1(nfb) = 1._DP 
                name_beam%beam_id(nfb)  =  "33RT"
              CASE DEFAULT
                PRINT *,'ERROR, beam number ',nfb
                PRINT *,'Is not one of the 14 beamlets recognized'
                lerrno = 359  + iomaxerr
                CALL terminate(lerrno,nlog)
           END SELECT


         
      CALL nfrpar_from_nubpar(name_beam%xybsca(nfb),  &
                              name_beam%xybapa(nfb),  &
                              name_beam%xlbapa(nfb),  &
                              name_beam%xlbtna(nfb),  &
                              name_beam%rtcena(nfb),  &
                              rpivot(nfb),            &   ! note input need to fix
                              anglev(nfb),            &   ! outputs ! 888889999
                              angleh(nfb),            &
                              zpivot(nfb),            &
                              blenp(nfb)              ) 
      IF( .NOT. name_beam%nlco(nfb) ) angleh(nfb) = -angleh(nfb)
      bheigh(nfb) = 2._DP*name_beam%bmwidza(nfb)
      bwidth(nfb) = 2._DP*name_beam%bmwidra(nfb)
      bhdiv(nfb)  = Deg_Per_Rad*name_beam%divra(nfb)
      bvdiv(nfb)  = Deg_Per_Rad*name_beam%divza(nfb)
      bhfoc(nfb)  = name_beam%foclra(nfb)
      bvfoc(nfb)  = name_beam%foclza(nfb)

       DO  iap =1,nap
         ! for table printout indicae that thes values are not used
         IF(nashape(iap,nfb) .EQ. 'b-d3d')awidth(iap,nfb) =   HUGE(1._DP)
         IF(nashape(iap,nfb) .EQ. 'b-circ' .OR.  nashape(iap,nfb) == 'b-d3d') aheigh(iap,nfb) =  HUGE(1._DP)
       ENDDO


 ! try the following
      !bvofset(nfb) = name_beam%xybsca(nfb)
      ! this is a guess
      !bhofset(nfb) = 0.0 ! ? related to xbzeta ??
      !if(nfb == 3) bvofset(nfb) =18.
      !if(nfb == 4) bvofset(nfb) =6.
      !if(nfb == 5) bvofset(nfb) =-6.
      !if(nfb == 6) bvofset(nfb) =-18.
      !if(nfb == 7) bvofset(nfb) =18.
      !if(nfb == 8) bvofset(nfb) =6.
      !if(nfb == 9) bvofset(nfb) =-6.
      !if(nfb == 10) bvofset(nfb) =-18.

      RETURN


  END SUBROUTINE set_P_Nfreya_geometry


  SUBROUTINE Check_P_Nfreya_beamlets(no_beams)
!--------------------------------------------------------------------------------------
!-- print out Nfreya beamlet specification. 
!--------------------------------------------------------------------------------------

  USE nrtype,                              ONLY :  DP,I4B,I2B

  USE MPI_data,                            ONLY :  myid,master

  USE nf_param,                            ONLY :  kb,ke

  USE neutral_beams,                       ONLY :  anglev, angleh,nashape, aheigh, &
                                                   awidth, bcur, bptor, blenp,    &
                                                   nbshape,bleni,bheigh, bwidth,  &
                                                   bhfoc, bvfoc, bhdiv, bvdiv,    &
                                                   ebkev, fbcur, nbeams,naptr,    &
                                                   alen, bvofset, bhofset, nsourc,&
                                                   sfrac1, rpivot,zpivot,npart_pseudo

  USE P_nfreya_interface,                  ONLY :  beam_data,d_fast_ion

  USE error_handler,                       ONLY :  lerrno,terminate,iomaxerr

  USE io_gcnmp,                            ONLY :  nlog

  USE shot_info,                           ONLY :  shot_id

  USE nub,                                 ONLY :  hdepsmth,fidiff_on

  USE fast_ion_diffusion,                  ONLY :  get_d_fast_ion

  IMPLICIT NONE

  INTEGER(I4B)    no_beams,nprint,iap,l,m,k,j,nob,ie,jb
  INTEGER(I2b) :: iout,get_next_io_unit
  LOGICAL cont


  If(nbeams .NE. no_beams)THEN
     lerrno = 375
  ENDIF


  IF(myid == master) THEN
     iout   = get_next_io_unit ()

     OPEN(unit=iout,file='Nfreya_beamlet_table',status='unknown')


     nprint = 7 ! <= 10 due to format below
     k      = 1
     cont   = .TRUE.

     DO WHILE(cont)
     WRITE(iout,6)shot_id%shot_nmbr
6    FORMAT(30x,"shot: ",i6)
     WRITE(iout,5)  
5    FORMAT(15x,'-----------------------------------Nfreya beamlet no-------------------------------------------')
        nob    = MIN(k+nprint-1,nbeams)
        WRITE(iout,10)(j,j=k,nob)
10      FORMAT(15x,10(2x,i3,10x))

        WRITE(iout,40)
40      FORMAT(2x,/,'Beamline related settings:')
        WRITE(iout,60)
60      FORMAT(2x,  '--------------------------')
        WRITE(iout,48)(beam_data%beam_id(l),l=k,nob)
48      FORMAT(2x,'beam id  =',t20,10(2x,a12))
        WRITE(iout,45)(anglev(l),l=k,nob)
45      FORMAT(2x,'anglev = ',t15,10(2x,1pe12.4))
        WRITE(iout,46)(angleh(l),l=k,nob)
46      FORMAT(2x,'angleh = ',t15,10(2x,1pe12.4))
        WRITE(iout,41)(nbshape(l),l=k,nob)
41      FORMAT(2x,'nbshape = ',t15,10(2x,a12))
        WRITE(iout,42)(rpivot(l),l=k,nob)
42      FORMAT(2x,'rpivot = ',t15,10(2x,1pe12.4))
        WRITE(iout,43)(zpivot(l),l=k,nob)
43      FORMAT(2x,'zpivot = ',t15,10(2x,1pe12.4))
        WRITE(iout,44)(blenp(l),l=k,nob)
44      FORMAT(2x,'blenp,cm = ',t15,10(2x,1pe12.4))
        WRITE(iout,29)(bleni(l),l=k,nob)
29      FORMAT(2x,'bleni,cm =',t15,10(2x,1pe12.4))
        WRITE(iout,50)(beam_data%rtcena(l),l=k,nob)
50      FORMAT(2x,'r tan,cm',t15,10(2x,1pe12.4))
        WRITE(iout,49)(beam_data%nlco(l),l=k,nob)
49      FORMAT(2x,'co inject  =',t15,10(2x,l12))



        WRITE(iout,30)
30      format(2X,/,'Beam input related settings:') 
        WRITE(iout,61)
61      FORMAT(2x,'----------------------------')
        IF(beam_data%NLBDAT)THEN
           WRITE(iout,12)
12         FORMAT(2x,'used ufile data input')
           WRITE(iout,65)
65         FORMAT(2x,'---------------------')
        ELSE
           WRITE(iout,13)
13         FORMAT(2x,'used namelist data input instead of ufile')
           WRITE(iout,66)
66         FORMAT(2x,'-----------------------------------------')
        ENDIF
        WRITE(iout,31)(bptor(l),l=k,nob)
31      FORMAT(2x,'bptor,watts =',t15,10(2x,1pe12.4))
        WRITE(iout,32)(ebkev(l),l=k,nob)
32      FORMAT(2x,'ebkev,KEV =',t15,10(2x,1pe12.4))
        DO m =1,ke
           WRITE(iout,25)(fbcur(m,l),l=k,nob)
           WRITE(iout,70)(npart_pseudo(m,l),l=k,nob) 
        ENDDO

25      FORMAT(2x,'fbcur =',t15,10(2x,1pe12.4))
70      FORMAT(2x,'nptcls =',t15,10(2x,i12))
        WRITE(iout,11)naptr
11      FORMAT(2x,/,'Aperture related settings,using (',i2,') apertures:')
        WRITE(iout,62)
62      FORMAT(2x,  '-----------------------------------------------------')
        WRITE(iout,80)
80      FORMAT(2x,"apertures must be in sequence,first two source ctr last two beam ctr !!")
        WRITE(iout,81)
81      FORMAT(2x,'-----------------------------------------------------------------------')
        WRITE(iout,82)
82      FORMAT(2x,'for b-circ aheigh is not used,for b-d3d aheigh,awidth are hard coded internally')
        WRITE(iout,83)
83      FORMAT(2x,'-------------------------------------------------------------------------------')
        DO iap = 1,naptr
           WRITE(iout,15)(alen(iap,l),l=k,nob)
           WRITE(iout,16)(aheigh(iap,l),l=k,nob)
           WRITE(iout,17)(awidth(iap,l),l=k,nob)
           WRITE(iout,18)(nashape(iap,l),l=k,nob)
15         FORMAT(2x,'alen =',t15,10(2x,1pe12.4))
16         FORMAT(2x,'aheigh =',t15,10(2x,1pe12.4))
17         FORMAT(2x,'awidth =',t15,10(2x,1pe12.4))
18         FORMAT(2x,'nashape =',t15,10(2x,a12))
           WRITE(iout,19)
19         FORMAT(/)
        ENDDO
  
        WRITE(iout,20)
20      format(2X,/,'Source related settings:')
        WRITE(iout,63)
63      FORMAT(2x,  '------------------------')
        WRITE(iout,21)(bheigh(l),l=k,nob)
21      FORMAT(2x,'bheigh =',t15,10(2x,1pe12.4))
        WRITE(iout,22)(bwidth(l),l=k,nob)
22      FORMAT(2x,'bwidth =',t15,10(2x,1pe12.4))
        WRITE(iout,23)(bvdiv(l),l=k,nob)
23      FORMAT(2x,'bvdiv =',t15,10(2x,1pe12.4))
        WRITE(iout,53)(bhdiv(l),l=k,nob)
53      FORMAT(2x,'bhdiv =',t15,10(2x,1pe12.4))
        WRITE(iout,24)(bvfoc(l),l=k,nob)
24      FORMAT(2x,'bvfoc =',t15,10(2x,1pe12.4))
        WRITE(iout,54)(bhfoc(l),l=k,nob)
54      FORMAT(2x,'bhfoc =',t15,10(2x,1pe12.4))
        WRITE(iout,26)(sfrac1(l),l=k,nob)
26      FORMAT(2x,'sfrac1 =',t15,10(2x,1pe12.4))
        WRITE(iout,27)(bvofset(l),l=k,nob)
27      FORMAT(2x,'bvofset =',t15,10(2x,1pe12.4))
        WRITE(iout,28)(bhofset(l),l=k,nob)
28      FORMAT(2x,'bhofset =',t15,10(2x,1pe12.4))


        k = k + nprint
        k=MIN(k,nbeams)
        IF(nob .GE. nbeams)cont = .FALSE.

     END DO ! do while cont = true



        IF(hdepsmth .GE. 0.0 .AND. fidiff_on )THEN

           !Load d_fast_ion diffusion parameters so they can
           ! be printed out. 
           DO jb =1,nbeams
              DO ie =1,ke
                 CALL get_d_fast_ion(d_fast_ion,ie,jb)
              ENDDO
           ENDDO

              WRITE(iout,715)d_fast_ion%anal_time
              WRITE(iout,716)
              WRITE(iout,717)
              WRITE(iout,718)
              WRITE(iout,711)(d_fast_ion%adiff_al,m=1,ke)  ! currently  no energy dependance
              WRITE(iout,712)(d_fast_ion%adiff_0l,m=1,ke) 
              WRITE(iout,713)(d_fast_ion%adiff_xpinl,m=1,ke)
              WRITE(iout,714)(d_fast_ion%adiff_xpoutl,m=1,ke)
        ENDIF
711      FORMAT(2x,'adiff_a=',t15,10(2x,1pe12.4),2x,"M^2/sec")
712      FORMAT(2x,'adiff_0=',t15,10(2x,1pe12.4),2x,"M^2/sec")
713      FORMAT(2x,'xpin=',t15,10(2x,1pe12.4))
714      FORMAT(2x,'xpout=',t15,10(2x,1pe12.4))
715      FORMAT(//,2x,'Profile smoothing info: evaluated  at time =',1pe14.6,',applied to all beams')

716      FORMAT(2x,   '------------------------------------------------------------------------------')
717      FORMAT(2x,'energy             full          half         third     )')
718      FORMAT(2x,'                   ----          ----         -----')
     CLOSE(iout)

  ENDIF ! myid = master

  RETURN

  END SUBROUTINE Check_P_Nfreya_beamlets
