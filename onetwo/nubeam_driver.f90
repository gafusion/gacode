

      SUBROUTINE  NTCC_driver
! -----------------------------------------------------------------
!  (a) create input file that will run with stand alone nubeam program
!  (b) spawn nubeam
!  (c) read results back into Onetwo
! -----------------------------------------------------HSJ 12/03------

        USE transp
!JMP    USE nbi_types
!JMP    USE nbi_dimensions
        USE solcon,ONLY : time,dt
        USE nub,ONLY  : pbeam,ebeam
        USE nub2, ONLY: bion,bneut
        USE nub4,ONLY : vbeam
        USE io,ONLY : ncrt,nout
        USE ext_prog_info, ONLY : get_nubeam,nubeam_path
        IMPLICIT NONE
        INTEGER ISHELL,io2,len_str,istat,i,j
        REAL*8 x,xxx
        LOGICAL first_time
        CHARACTER *256 nubeam_path_out,command, &
        nubeamout_spawn,nubeam_setup_out
        CHARACTER(len=24) nubeam_namelist
        CHARACTER (len=12) atime,intfl
        CHARACTER(len= 256) nubeam_output ! must be > nubeam_root
        DATA first_time /.TRUE./


        len_str = LEN_TRIM(nubeam_path)
        IF (nubeam_restart .EQ. 1)THEN
                   first_time = .FALSE.
           CALL get_nubeam(ncrt,nout,nubeam_path_out, &
                           nubeam_setup_out,len_str)
        ENDIF
        IF(first_time)THEN
          print *,'len-str =',len_str
          print *,'nubeam_path=',nubeam_path
         !get fully qualified name of nubeam to run:
           CALL get_nubeam(ncrt,nout,nubeam_path_out, &
                           nubeam_setup_out,len_str)

           IF(len_str .LE. 0)THEN
              PRINT *,'sub get_nubeam did not find required path'
              PRINT *,'nubeam_path_out = ',nubeam_path_out
              CALL STOP('nubeam path problem',1)
           ENDIF


           !set some environment variables required by nubeam module:
           !because nubeam is run as a separate process I have not found a  way to
           !make this information available in nubeam_driver.
           !command = ' '//ADJUSTL(nubeam_setup_out(1:LEN_TRIM( nubeam_setup_out)))
           !WRITE(6, FMT ='(" running : ",a)')command
           !IF (ISHELL (command) .LT. 0)         &
           !    CALL STOP ('subroutine nubeam_driver env set', 1)

           nubeam_dt = beam_data%beam_power_rise_time




!          nubeam_calls may be allocated in sub beam_prof_init
!          (in module nbi_restart called at begining of tport)
!           IF(.NOT. ALLOCATED(nubeam_calls))
           IF(.NOT. ASSOCIATED(nubeam_calls))                           &
           ALLOCATE (nubeam_calls(0:10),STAT = istat) !start with 10
                                                      !and increas size 
                                                      !as necessary
           IF(istat .NE. 0)                             &
             CALL allocate_error("nubeam_calls,sub NTCC_drivr ",0,istat)
          nubeam_calls(:) =0.0
          first_time = .FALSE.
     

       ELSE
           !ALSO ENTER HERE IF THIS IS THE FIRT CALL TO 
           !NUBEAM FOR A RESTART CASE
           !get nubeam_dt here
           xxx = 1.e30
           DO j=1,SIZE(beam_data%beam_times)
              x= beam_data%beam_times(j)-nubeam_calls(nubeam_steps-1)
              IF(x .GT.0.0)THEN
                 xxx = MIN(x,xxx)
              ENDIF
           ENDDO
           nubeam_dt = MIN(nubeam_dt,xxx) ! nubeam_dt is set to nubeam_dt0
                                          ! at start of each time step
           nubeam_dt = MAX(nubeam_dt,                                    &
                     beam_data%beam_power_rise_time*0.5)
        ENDIF

        IF (ifix_nubeam_dt .EQ. 1) THEN !JMP BLOCK START 
       	  nubeam_dt = nubeam0_dt
        END IF !JMP BLOCK END

!       write input file for Nubeam:
        nubeam_namelist = nubeam_root(1:LEN_TRIM(nubeam_root))//'_namelist.dat'
        CALL wrt_nubeam_in(nubeam_namelist,time)

        if (time .ge. nubeam_fix_t) then !JMP BLOCK START
          print *,'NUBEAM CALL SKIPPED'
          return
        end if  !JMP BLOCK END

        command = ADJUSTL(nubeam_path_out(1:LEN_TRIM( nubeam_path_out)))
        command = command(1:LEN_TRIM(command))                        &
           //' '//nubeam_root(1:LEN_TRIM(nubeam_root))
!         execute nubeam  code  and then connect to file 
!

      WRITE  (6, '(/ '' ---- NUBEAM  started'')')
      WRITE(6, FMT ='(" running : ",a)')command
!
         i=0
         i = ISHELL (command) 
         IF(i .LT. 0)THEN      
            PRINT *,'return code from Nubeam is ',i
            CALL STOP ('subroutine nubeam_driver spawned NUBEAM', 67)
         ENDIF
!
      WRITE  (6, '(/ '' ---- NUBEAM finished'')')
      IF(save_nubeam_input .EQ. 1)THEN
          WRITE(intfl,FMT='(f12.6)')time
          READ(intfl,FMT='(a)')atime
          atime = ADJUSTL(atime)
          command = 'mv '//nubeam_namelist(1:LEN_TRIM(nubeam_namelist))  &
            //' '//nubeam_namelist(1:LEN_TRIM(nubeam_namelist))//'_'//   &
                 atime(1:LEN_TRIM(atime))
          IF (ISHELL (command) .LT. 0)                           &
          CALL STOP ('sub NTCC_driver: failure of spawned mv command', 67)
      ENDIF



!
!     Read data produced by nubeam code,
!     and put into form for output variables in this subroutine
!
       nubeam_output = nubeam_root(1:LEN_TRIM(nubeam_root))//'_details.dat'
       CALL getioun(io2,42)
       OPEN (unit = io2, file = nubeam_output, status = 'OLD')

       CALL load_12(io2)       !put info from nubeam into Onetwo slots
 
!       CALL load_restart_12


       CLOSE (unit = io2)
       CALL giveupus(io2)


      RETURN
      END





      SUBROUTINE wrt_nubeam_in(nubeam_namelist,time_eqdsk)
! -------------------------------------------------------------------------
!    subroutine writes namelist for nubeam code
!
! ------------------------------------------------------------------------
       USE param, ONLY : kj,kprim,kimp
       USE ename ,ONLY : eqdskfilename,eqfile,eqdsk_tdem
       USE transp,ONLY : beam_data,nubeam_nclass,nlfhe3,nlfhe4,   &
                         nlfst,nlfsp, nlusf3,xdatsfa, xdatsf3,    &
                         xdatsft,xdatsfp,nlusfa,nlusfp,nlusft,    &
                         plfhe3,plfhe4,plfst,plfsp,dxbsmoo,wghta, &
                         nptcls,nptclf,ndep0,goocon,dtn_orbit,    &
                         nsigexc,nlminsv,nlbbcx,nlebei,dn0out,    &
                         nmsigx,xcfanbi,xdfanbi,xefanbi,xcfafus,  &
                         xdfafus,xefafus,nlbcde,nlbcoh,nlbcpa,    &
                         nlorbo,nlbflr,nlfbmflr,nlbgflr,gflr_min, &
                         nkdifb,fdifbe,edifbe,plhgt,nlcprb,       &
                         nmcurb,xdepmod,xp_nt1,inzri,             &
                         ndifbep,nubeam_restart,         &
                         nubeam_dt,nubeam_calls,nubeam_steps,     &
                         d_fast_ion
       USE solcon,    ONLY   : time,dtt,dt,time_tol,time0 !JMP
       USE numbrs,    ONLY   : nj,nprim,nimpc=>nimp 
       USE ions ,     ONLY   : zeffc => zeff
       USE soln ,     ONLY   : etor,ene,te,ti,en      !en(kj,kion)
       USE tordlrot,  ONLY   : iangrot, angrot
       USE machin,    ONLY   : rmajor,btor
       USE constnts,  ONLY   : pi
       USE neut ,     ONLY   : ennw,ennv
       USE nbnamelist
       IMPLICIT NONE
       INTEGER lenge,io12,err_state,err_xplasma,j,klc,ngmaxlc,     &
               nrhixlc

       CHARACTER(len = 64) runid_nubeam
       CHARACTER(*) nubeam_namelist
       CHARACTER (len = 64) eqdsk_name
       REAL*8 time_eqdsk,time_pwr

   INTERFACE local

    SUBROUTINE REALLOCATE1DD0(A,KM) 
      IMPLICIT NONE
      REAL *8,  DIMENSION(:),POINTER :: a
      INTEGER km,sizea
    END SUBROUTINE REALLOCATE1DD0

   END INTERFACE





       nstep_start = nubeam_steps



       !normally we want to take a nubeam time step of size
       !nubeam0_dt which is typically much greater than dt.
       !but we dont want to skip over events. The current value of
       !nubeam_dt is set in chekdt and nubeam0_dt is the max that
       !nubeam_dt can be.

  
       kj_check = kj
       time_start = time
       time_stop  = time +nubeam_dt



   IF( nubeam_restart .GT. 0)THEN 

         !user wishes to use existing restart files on the 
         !first call to nubeam. Make this possible by 
         !increasing nstep_start from 0 to 1.
         !but first check for valid files:
         
         CALL check_restart(err_state,err_xplasma)
              nstep_start = nstep_start+1
         IF(ABS(time_start-time_stop) .LE. 2.*time_tol)THEN
            time_stop= time_start + beam_data%beam_power_rise_time
            nubeam_dt = beam_data%beam_power_rise_time
         ENDIF
         nubeam_restart = 0
   ENDIF



       time_pwr = time + 0.5*nubeam_dt  ! get the beam power in the middle
                                        ! of the analysis interval
       tr_delta_t =0.0                  ! means no time averaging

       tr_sol_width = 5.0d0             ! 5cm  presumed scrape-off region width

       nsnccwi_tr =1                    !=1 if I  is ccw from top view
       nsnccwb_tr = INT(SIGN(1.0d0,btor)) !=1 if btor is ccw
       nRgrid_tr = 100 ;   nZgrid_tr = 100    ! (R,Z) grid sizes -- default: 100x100 OK.
 
       !load current beam powers into beam_data%pinja:
 
       CALL  beam_power_interp(time_pwr)


       nsteps      = 1             !take one step of duration dtt
                                   !(eg, time_stop-time_start)
                                   !actually not used in  onetwo version
                                   !nbdrive_main.f90


 
        nstep_stop  = 0             !not used 









! eqdsk input :-------------------------------------------------------
        runid_nubeam ='NUBEAM  input file created by ONETWO'
        eqfile = eqdskfilename
        IF(eqdsk_tdem .NE. 'tdem' ) THEN
             lenge = LEN_TRIM(eqdskfilename)
             eqdsk_name  = eqdskfilename(1:lenge)
        ELSE
             if((nubeam_restart.eq.-1).and.(time .lt. time0)) then !JMP BLOCK START
             	eqdsk_name='g0.restart'
             	print *,time,time0,'g0.restart'
             else
                CALL wrt_tdem_eqdsk(time_eqdsk,eqdsk_name)
                print *,time,time0,'neweq' 
             endif !JMP BLOCK END
             eqfile = eqdsk_name
             PRINT *,'eqdsk file created in sub write_nubeam_in'
        ENDIF
        PRINT *,'eqdsk_name =', eqdsk_name
        PRINT *,'eqfile =',eqfile
        PRINT *,'eqdsk_tdem =',eqdsk_tdem
        IF (eqfile .EQ. 'none')  eqdsk_name =  'eqdskin'
        lenge = LEN_TRIM(eqdsk_name)
        lenge = lenge + 2 ! +EFIT: + ' + '
        IF(lenge .GT. LEN( mhdpath))THEN
           CALL STOP('dim of mhdpath too small',1)
        ELSE
           mhdpath ='"'// 'EFIT:'//eqdsk_name(1:LEN_TRIM(eqdsk_name))//'"'
        ENDIF




 
! assign values from Onetwo input
        IF(goocon .EQ. 0.0) goocon          = 10.0d0 
        dtn             = dtn_orbit
        nseed           = beam_data%nseed
        nclass          = nubeam_nclass
        nznbma          = beam_data%nznbma 
        nznbme          = beam_data%nznbme
        nzones          = beam_data%nzone_nb
        ebdmax          = beam_data%ebdmax 
!        ngmax           = beam_data%ngmax
        nbeam           = beam_data%nbeam
        nbbcal          = beam_data%nbbcal
        nlhvion         = .FALSE.   !for heavy ion injection only
! anomalous diffusion of fast ions
        adiff_0 = d_fast_ion%adiff_0
        adiff_a = d_fast_ion%adiff_a
        adiff_xpin = d_fast_ion%adiff_xpin
        adiff_xpout = d_fast_ion%adiff_xpout
        IF(adiff_0*adiff_a .GT. 0.0)d_fast_ion%fidif_on = 1
        !       energy dependent p[art of fi diffusion:
        ndifbe   = d_fast_ion%ndifbe
        nkdifb   = d_fast_ion%nkdifb
        IF(ndifbe .GT. 0)THEN
           fdifbe(1:ndifbe) = d_fast_ion%fdifbe(1:ndifbe)
           edifbe(1:ndifbe) = d_fast_ion%edifbe(1:ndifbe)
        ELSE
           fdifbe(:) =0.0
           edifbe(:) =0.0
        ENDIF

! first check impurities. In Onetwo an impurity may have mass no of 2
! this is not allowed in nubeam. SO we have to do a little dance here:
        nrhixlc  = beam_data%nrhix
        ngmaxlc  = beam_data%ngmax
        klc = 0
        IF(beam_data%nrhix  .LE. nrhixm)THEN
           DO j =1,beam_data%nrhix
               IF(beam_data%xzimpx(j) .LT. 3.) THEN
                  nrhixlc = nrhixlc -1
                  ngmaxlc = ngmaxlc +1
                  IF(ngmaxlc .GT. kprim)              &
                     CALL STOP('kprim too small for ni_kj, in nubeam_driver',1)
                  !impurity species in Onetwo will be represented as
                  !a primary species in Nubeam:
                  aplasm(ngmaxlc) = beam_data%aimpx(j)
                  backz(ngmaxlc) =   beam_data%xzimpx(j)
                  ni_kj(:,ngmaxlc) = en(:,nprim+j)
               ELSE
                  klc =klc+1
                  xzimpx(klc)    = beam_data%xzimpx(j)
                  aimpx(klc)     = beam_data%aimpx(j)
                  nimp_kj(:,klc) = en(:,nprim+j)
               ENDIF
           ENDDO
        ELSE
           PRINT *,'beam_data%nrhix =',beam_data%nrhix
           PRINT *,'nrhixm in nubeam_driver =',nrhixm
           CALL STOP('ERROR: beam_data%nrhix  too large',1)
        ENDIF

        nrhix = klc        !nrhix will be passsed to nubeam
        ngmax = ngmaxlc    !ngmax will be pased to nubeam







!profile_input: --------------------------------------------------------
! parabolic inputparameters are not used by Onetwo created file
!so just set them arbitrarily (but nubeam checks  these so
!reasonable values need to be given)
        IF(xp_nt1 .EQ. 0)xp_nt1 = 65
        ntheta = xp_nt1
        nrho = nj
        nzone_fb = inzri
        zeff(1:kj) = zeffc(1:kj)
! temperatures (te_kj,ti_kj in eV for namelist input in nubeam_driver)
        te_kj(:) = 1000.*te(:) 
        ti_kj(:) = 1000.* ti(:) 

! electron density (ne_kj in cm**-3 for namelist input in nubeam_driver)
        ne_kj(:) = ene(:)

!primary ion density(position beyond nprim were loaded above):
        ni_kj(:,1:nprim) = en(:,1:nprim)



!impurity  density is taken care of above

!neutral density:

        n0w_kj(:,:) = ennw(:,:)
        n0v_kj(:,:) = ennv(:,:)
! toroidal rotation, rad/sec:
      IF(iangrot .GT. 0)THEN

            omega_kj(:) = angrot(:)
      ELSE
            omega_kj(:) = 0.0d0
      ENDIF

!       loop voltage, volts. We have V = 2 Pi R0 Etor
        vloop_kj(:) = 2.*pi*rmajor*etor(:)


!       electrostatic radial potential,volts. Onetwo does not
!       have this quantity so we use the parabolic form:
        epot_0=1000.    ; epot_a = 0.0 ; epot_xpin = 2.0 ; epot_xpout = 1.



!    following parabolic profiles are not used. They need to be set however
!    to bypass nubeam error checking:
        te_0 =1000.0 ; te_a =10.0   ; te_xpin =2.0    ; te_xpout =1.
        ti_0 =1000.0 ; ti_a =10.0   ; ti_xpin =2.0    ; ti_xpout =1.
        ne_0 =1.0e13 ; ne_a = 1.e12 ; ne_xpin =2.0    ; ne_xpout =1.
        n0w_0 =0.0   ; n0w_a =0.0   ; n0w_xpin =2.0   ; n0w_xpout =1.
        n0w_0 =0.0   ; n0w_a =0.0   ; n0w_xpin =2.0   ; n0w_xpout =1.
        omeg_0= 100. ; omeg_a= 10.  ; omeg_xpin=2.0   ; omeg_xpout=1.
        vloop_0=1.   ; vloop_a =1.  ; vloop_xpin= 2.0 ; vloop_xpout=1.













        IF(beam_data%ngmax .LE. ngmaxx)THEN
           backz(1:beam_data%ngmax) = beam_data%backz(1:beam_data%ngmax)
           aplasm(1:beam_data%ngmax) = beam_data%aplasm(1:beam_data%ngmax)
        ELSE
           CALL STOP('ERROR: beam_data%ngmax too large',1)
        ENDIF
        IF(beam_data%nbeam .LE. nbeamx)THEN
           xzbeama(1:beam_data%nbeam)  = beam_data%xzbeama(1:beam_data%nbeam)
           abeama(1:beam_data%nbeam)   = beam_data%abeama(1:beam_data%nbeam)
           nlco(1:beam_data%nbeam)     = beam_data%nlco(1:beam_data%nbeam)
           ntrace(1:beam_data%nbeam)   = beam_data%ntrace(1:beam_data%nbeam)
           rtcena(1:beam_data%nbeam)   = beam_data%rtcena(1:beam_data%nbeam)
           xlbtna(1:beam_data%nbeam)   = beam_data%xlbtna(1:beam_data%nbeam)
           xybsca(1:beam_data%nbeam)   = beam_data%xybsca(1:beam_data%nbeam)
           xbzeta(1:beam_data%nbeam)   = beam_data%xbzeta(1:beam_data%nbeam)
           nbshapa(1:beam_data%nbeam)  = beam_data%nbshapa(1:beam_data%nbeam)
           bmwidra(1:beam_data%nbeam)  = beam_data%bmwidra(1:beam_data%nbeam)
           bmwidza(1:beam_data%nbeam)  = beam_data%bmwidza(1:beam_data%nbeam)
           divra(1:beam_data%nbeam)    = beam_data%divra(1:beam_data%nbeam)
           divza(1:beam_data%nbeam)    = beam_data%divza(1:beam_data%nbeam)
           foclra(1:beam_data%nbeam)   = beam_data%foclra(1:beam_data%nbeam)
           foclza(1:beam_data%nbeam)   = beam_data%foclza(1:beam_data%nbeam)
           nbapsha(1:beam_data%nbeam)  = beam_data%nbapsha(1:beam_data%nbeam)
           rapedga(1:beam_data%nbeam)  = beam_data%rapedga(1:beam_data%nbeam)
           xzpedga(1:beam_data%nbeam)  = beam_data%xzpedga(1:beam_data%nbeam)
           xlbapa(1:beam_data%nbeam)   = beam_data%xlbapa(1:beam_data%nbeam)
           xybapa(1:beam_data%nbeam)   = beam_data%xybapa(1:beam_data%nbeam)
           nbapsh2(1:beam_data%nbeam)  = beam_data%nbapsh2(1:beam_data%nbeam)
           rapedg2(1:beam_data%nbeam)  = beam_data%rapedg2(1:beam_data%nbeam)
           xzpedg2(1:beam_data%nbeam)  = beam_data%xzpedg2(1:beam_data%nbeam)
           xlbapa2(1:beam_data%nbeam)  = beam_data%xlbapa2(1:beam_data%nbeam)

           einja(1:beam_data%nbeam)    = beam_data%einja(1:beam_data%nbeam)
           pinja(1:beam_data%nbeam)    = beam_data%pinja(1:beam_data%nbeam)
           ffulla(1:beam_data%nbeam)   = beam_data%ffulla(1:beam_data%nbeam)
           fhalfa(1:beam_data%nbeam)   = beam_data%fhalfa(1:beam_data%nbeam)
           tbona(1:beam_data%nbeam)   = beam_data%tbona(1:beam_data%nbeam)
           tboffa(1:beam_data%nbeam)  = beam_data%tboffa(1:beam_data%nbeam)
        ELSE
           CALL STOP('ERROR: beam_data%nbeam too large',1)
        ENDIF



       CALL getioun(io12,42)
!       OPEN  (unit = io12,                                           &
!            file = nubeam_namelist(1:LEN_TRIM(nubeam_namelist)),     &
!            status = 'UNKNOWN')
       OPEN  (unit = io12,                                           &
            file = nubeam_namelist(1:LEN_TRIM(nubeam_namelist)),     &
            status = 'UNKNOWN',RECL = 100 )   ! RECL required for lf95
       WRITE (unit = io12, fmt = '(3x, a)') runid_nubeam
       WRITE (unit = io12, nml = nbdrive_naml)

       CLOSE (unit = io12)
       CALL giveupus(io12)

       IF(nubeam_steps .GT. SIZE(nubeam_calls) - 1)                  &
            CALL reallocate1dd0(nubeam_calls,nubeam_steps+10)
       nubeam_calls(nubeam_steps) = time_stop
       IF(nubeam_steps .EQ. 1)                                       &
                    nubeam_calls(nubeam_steps-1)=time_start
       nubeam_steps = nubeam_steps + 1     
                                   !set for next call to this routine
                                   !(each new time interval starts with the
                                   !saved file of the previous time interval
                                   !for the first time this routine is called
                                   !there is no saved file, indicated by
                                   !nstep_start = 1)
       PRINT *,'wrt_nubeam_in run nubeam from ',time_start
       PRINT *,'to ',time_stop
       PRINT *,'that is ',nubeam_calls(nubeam_steps-2),              &
                          nubeam_calls(nubeam_steps-1)
       PRINT *,'nubeam_steps is now ',nubeam_steps
       !JMP PRINT *,'nubeam_calls =',nubeam_calls
       RETURN
       END SUBROUTINE wrt_nubeam_in










       SUBROUTINE  check_restart (err_state,err_xplasma)
! ---------------------------------------------------------------------------
! check for valid restart files for nubeam
! Besides the two files nubeam_state_path_nubeam_state.cdf and
! nubeam_xplasma_path_xplasma_state.cdf created by nubeam (and read
! by nubeam) we also require a Onetwo restart file that has all the profiles
! in it needed to initialize the fast ion related profiles in Onetwo.
! This latter file is called nubeam_state_path(1:k)_profile_file
! and is created by sub  write_restart_profs and  read by
! sub  read_restart_profs in Onetwo proper.
! -----------------------------------------------------------------HSJ-------
       USE transp, ONLY : nubeam_state_path, nubeam_xplasma_path, nubeam_root, &
                          nubeam_profile_path
       IMPLICIT NONE

       INTEGER,INTENT(OUT)    :: err_state,err_xplasma
       INTEGER*4  GETCWD,VERIFY,ier,k,l,ISHELL
       CHARACTER  *256 state_file
       CHARACTER  *256 xplasma_file
       CHARACTER  *256 profile_file
       CHARACTER  *256 cwd,filename,new_filename,command
       LOGICAL     ex,ex1

       err_state = 0 ; err_xplasma = 0 ;ier =0
       nubeam_state_path = ADJUSTL(nubeam_state_path)
       nubeam_xplasma_path = ADJUSTL(nubeam_xplasma_path)


       ! file names  used by nubeam, dont change  unless nubeam changes:
       state_file = ADJUSTL(nubeam_root(1:LEN_TRIM(nubeam_root)))//       &
                   "_nubeam_state.cdf"

       xplasma_file =  ADJUSTL(nubeam_root(1:LEN_TRIM(nubeam_root)))//   &
                   "_xplasma_state.cdf"   


       ! file name used by Onetwo, sub read_restart_profs and 
       ! write_restart_profs
        profile_file = ADJUSTL(nubeam_root(1:LEN_TRIM(nubeam_root)))//       &
                   "_restart_profs.txt"


       !get current working directory:
        ier = getcwd(cwd)
        IF(ier .NE. 0)CALL STOP('getcwd error in check_restart_profs',1)
        cwd =cwd(1:LEN_TRIM(cwd))//'/'
        cwd = ADJUSTL(cwd)

!        print *,'nubeam_state_path =',nubeam_state_path(1:LEN_TRIM(nubeam_state_path))
       !does nubeam_state_path directory exist, and is it the same as cwd ??
        l =LEN_TRIM(nubeam_state_path)
        IF(nubeam_state_path(1:2) == './')          &   ! swap ./ for full path
            nubeam_state_path = cwd(1:LEN_TRIM(CWD))//nubeam_state_path(3:l)
        ier = -1
        l =LEN_TRIM(nubeam_state_path)
        ier  = VERIFY(cwd,nubeam_state_path)
        k = INDEX(nubeam_state_path,'/',BACK = .TRUE.)

        nubeam_profile_path = nubeam_state_path(1:k)//profile_file(1:LEN_TRIM(profile_file))

        filename = ADJUSTL(nubeam_state_path(k+1:l))
        INQUIRE(FILE =nubeam_state_path,EXIST  = ex)
        IF( .NOT. ex)THEN
          PRINT *,'ERROR,filename =',filename
          PRINT *,'in directory ',nubeam_state_path(1:k)
          PRINT *,'does not exist'
          CALL STOP('check_restart: file doesnt exist',1)
        ENDIF
        IF(ier .EQ. 0)THEN      !cwd is the same as path part of nubeam_state_path
             INQUIRE(FILE =filename,EXIST  = ex)
             IF(ex)THEN          !file exists, copy it to working file name  if 
                                 !necessary .
                 IF(filename(1:LEN_TRIM(filename))  ==                     & 
                                 state_file(1:LEN_TRIM(state_file))) THEN
                      new_filename = filename(1:LEN_TRIM(filename))//"_orig"
                      command = "cp "//filename(1:LEN_TRIM(filename))//" "// &
                                      new_filename(1:LEN_TRIM(new_filename))
!                       print *,'new_filename 1 =',new_filename(1:LEN_TRIM(new_filename))
                 ELSE       !filename not the same as state_filenmae
                       new_filename = state_file(1:LEN_TRIM(state_file))
                       command ="cp "//filename(1:LEN_TRIM(filename))//" "// &
                                              state_file(1:LEN_TRIM(state_file))
!                       print *,'new_filename 2 =',new_filename(1:LEN_TRIM(new_filename))
                 ENDIF

             ENDIF
             
        ELSE                    !nubeam_state_path points to different directory/file ,
                                !copy file to here
            new_filename = cwd(1:LEN_TRIM(cwd))//state_file(1:LEN_TRIM(state_file))
            command ="cp " // nubeam_state_path(1:LEN_TRIM(nubeam_state_path)) &
                  //" "// new_filename(1:LEN_TRIM(new_filename))

        ENDIF 
!        print *,'command =',command(1:LEN_TRIM(command))
        ier = ISHELL (command)
        IF(ier .NE.0)CALL STOP('ISHELL cp command error',1)
        INQUIRE(FILE =state_file,EXIST  = ex)
        IF( .NOT. ex)THEN
          PRINT *,'ERROR, file ',state_file(1:LEN_TRIM(state_file))
          PRINT *,'This file must exist in cwd'
          CALL STOP('check_restart: file does not exist',1)
        ENDIF




!        print *,'nubeam_profile_path =',nubeam_profile_path(1:LEN_TRIM(nubeam_profile_path))
       !does nubeam_profile_path directory exist, and is it the same as cwd ??
        l =LEN_TRIM(nubeam_profile_path)
        IF(nubeam_profile_path(1:2) == './')          &   ! swap ./ for full path
            nubeam_profile_path = cwd(1:LEN_TRIM(CWD))//nubeam_profile_path(3:l)
        ier = -1
        l =LEN_TRIM(nubeam_profile_path)
        ier  = VERIFY(cwd,nubeam_profile_path)
        k = INDEX(nubeam_profile_path,'/',BACK = .TRUE.)
        filename = ADJUSTL(nubeam_profile_path(k+1:l))
        INQUIRE(FILE =nubeam_profile_path,EXIST  = ex)
        IF( .NOT. ex)THEN
          PRINT *,'ERROR,filename =',filename
          PRINT *,'in directory ',nubeam_profile_path(1:k)
          PRINT *,'does not exist'
          CALL STOP('check_restart:file doesnt exist',1)
        ENDIF
        IF(ier .EQ. 0)THEN      !cwd is the same as path part of nubeam_profile_path
             INQUIRE(FILE =filename,EXIST  = ex)
             IF(ex)THEN          !file exists, copy it to working file name  if 
                                 !necessary .
                 IF(filename(1:LEN_TRIM(filename))  ==                     & 
                                 profile_file(1:LEN_TRIM(profile_file))) THEN
                      new_filename = filename(1:LEN_TRIM(filename))//"_orig"
                      command = "cp "//filename(1:LEN_TRIM(filename))//" "// &
                                      new_filename(1:LEN_TRIM(new_filename))
!                       print *,'new_filename 1 =',new_filename(1:LEN_TRIM(new_filename))
                 ELSE       !filename not the same as profile_filenmae
                       new_filename = profile_file(1:LEN_TRIM(profile_file))
                       command ="cp "//filename(1:LEN_TRIM(filename))//" "// &
                                              profile_file(1:LEN_TRIM(profile_file))
!                       print *,'new_filename 2 =',new_filename(1:LEN_TRIM(new_filename))
                 ENDIF

             ENDIF
             
        ELSE                    !nubeam_profile_path points to different directory/file ,
                                !copy file to here
            new_filename = cwd(1:LEN_TRIM(cwd))//profile_file(1:LEN_TRIM(profile_file))
            command ="cp " // nubeam_profile_path(1:LEN_TRIM(nubeam_profile_path)) &
                  //" "// new_filename(1:LEN_TRIM(new_filename))

        ENDIF 
!        print *,'command =',command(1:LEN_TRIM(command))
        ier = ISHELL (command)
        IF(ier .NE.0)CALL STOP('ISHELL cp command error',1)
        INQUIRE(FILE =profile_file,EXIST  = ex)
        IF( .NOT. ex)THEN
          PRINT *,'ERROR, file ',profile_file(1:LEN_TRIM(profile_file))
          PRINT *,'This file must exist in cwd'
          CALL STOP('check_restart: file does not exist',1)
        ENDIF







!
!        print *,'nubeam_xplasma_path =',nubeam_xplasma_path(1:LEN_TRIM(nubeam_xplasma_path))
       !does nubeam_xplasma_path directory exist, and is it the same as cwd ??
        l =LEN_TRIM(nubeam_xplasma_path)
        IF(nubeam_xplasma_path(1:2) == './')          &   ! swap ./ for full path
            nubeam_xplasma_path = cwd(1:LEN_TRIM(CWD))//nubeam_xplasma_path(3:l)
        ier = -1
        l =LEN_TRIM(nubeam_xplasma_path)
        ier  = VERIFY(cwd,nubeam_xplasma_path)
        k = INDEX(nubeam_xplasma_path,'/',BACK = .TRUE.)
        filename = ADJUSTL(nubeam_xplasma_path(k+1:l))
        INQUIRE(FILE =nubeam_xplasma_path,EXIST  = ex)
        IF( .NOT. ex )THEN
             PRINT *,'filename =',filename
             PRINT *,'in directory ',nubeam_xplasma_path(1:k)
             PRINT *,'does not exist'
!           CALL STOP('check_restart: file doesnt exist',1) file may not exist if new mhd case
        ENDIF
        IF(ier .EQ. 0 .AND. ex)THEN      !cwd is the same as path part of nubeam_xplasma_path
             INQUIRE(FILE =filename,EXIST  = ex)
             IF(ex)THEN          !file exists, copy it to working file name  if 
                                 !necessary .
                 IF(filename(1:LEN_TRIM(filename))  ==                     & 
                                 xplasma_file(1:LEN_TRIM(xplasma_file))) THEN
                      new_filename = filename(1:LEN_TRIM(filename))//"_orig"
                      command = "cp "//filename(1:LEN_TRIM(filename))//" "// &
                                      new_filename(1:LEN_TRIM(new_filename))
!                     print *,'new_filename 1 =',new_filename(1:LEN_TRIM(new_filename))
                 ELSE       !filename not the same as xplasma_filenmae
                      new_filename = xplasma_file(1:LEN_TRIM(xplasma_file))
                      command ="cp "//filename(1:LEN_TRIM(filename))//" "// &
                                              xplasma_file(1:LEN_TRIM(xplasma_file))
!                     print *,'new_filename 2 =',new_filename(1:LEN_TRIM(new_filename))
       
                 ENDIF
                 ier = ISHELL (command)

             ENDIF
             
        ELSEIF (ex) THEN                  !nubeam_xplasma_path points to different directory/file ,
                                !copy file to here
             new_filename = cwd(1:LEN_TRIM(cwd))//xplasma_file(1:LEN_TRIM(xplasma_file))
             command ="cp " // nubeam_xplasma_path(1:LEN_TRIM(nubeam_xplasma_path)) &
                  //" "// new_filename(1:LEN_TRIM(new_filename))
             ier = ISHELL (command)
        ELSE   
           ! file does not exist this should happen only if this is
           ! a new mhd case
           ier =0
        ENDIF 

        IF(ier .NE.0)CALL STOP('ISHELL cp command error',2)
!        print *,'command =',command(1:LEN_TRIM(command))
        IF(ex) INQUIRE(FILE =xplasma_file,EXIST  = ex1)
        IF(ex .AND.  .NOT. ex1)THEN
          PRINT *,'ERROR, file ',xplasma_file(1:LEN_TRIM(xplasma_file))
          PRINT *,'This file must exist in cwd'
          CALL STOP('check_restart: file does not exist',2)
        ENDIF


       RETURN
       END SUBROUTINE  check_restart


       SUBROUTINE set_nubeam_init_profs
!  ----------------------------------------------------------------------
!   Set the profiles determined by nubeam
!   Both restart and NEW cases are handled 
!  ----------------------------------------------------------------------

      USE transp,ONLY :                                                  &
                       qbeame_nubp,qbeami_nubp,                          &
                       qbeame_intg_nubp,qbeami_intg_nubp,                &
                       qbth_nubp,qbth_intg_nubp,curb_nubp,               &
                       storqueb_nubp,storqueb_intg_nubp,                 &
                       sprbeame_nubp,storqueb_nub,                       &
                       sprbeami_nubp,wbeam_nubp,                         &
                       enbeam_nubp,sprbeami_intg_nubp,                   &
                       sprbeame_intg_nubp,wbeam_intg_nubp,               &
                       pwf_tot_intg_nub,pwf_tot_intg_nubp,               &
                       enbeam_intg_nubp,                                 &
                       beam_thermal_dtntot_nubp,                         &
                       beam_thermal_ddntot_nubp,                         &
                       beam_thermal_ddptot_nubp,                         &
                       beam_thermal_tt2ntot_nubp,                        &
                       beam_beam_dtntot_nubp,                            &
                       beam_beam_ddntot_nubp,                            &
                       beam_beam_ddptot_nubp,beam_beam_tt2ntot_nubp,     &
                       beam_thermaltth_df_nubp,beam_thermalddp_nubp,     &
                       beam_thermalddn_nubp,beam_thermaltt2n_nubp,       &
                       beam_thermaldth_tf_nubp,beam_beamddn_nubp,        &
                       beam_beamdtn_nubp,beam_beamddp_nubp,              &
                       beam_beamtt2n_nubp,curb_intg_nubp,                &
                       nubeam_restart,nubeam_state_path,nubeam_root,     &
                       nubeam_profile_path,enbeam_species,               &
                       enbeam_species_p, &
                       sorbn0_nub,sorbn0_nubp,sorbh_nub,sorbh_nubp !jmp.den
                       

       USE restore_12,      ONLY : storqueb_nubr   ! HSJ 07/07/05

      IMPLICIT NONE
      INTEGER err_state,err_xplasma



            IF((nubeam_restart .EQ. 0).OR.(nubeam_restart .EQ. -1).OR. & !JMP
               (nubeam_restart .EQ.-2).OR.(nubeam_restart .EQ.-99))THEN !JMP
              qbeame_nubp(:)                = 0.0D0
              qbeami_nubp(:)                = 0.0D0
              qbth_nubp(:)                  = 0.0D0
              curb_nubp(:)                  = 0.0D0
              storqueb_nubp(:)              = 0.0D0
              sprbeame_nubp(:)              = 0.0D0
              sprbeami_nubp(:)              = 0.0D0
              wbeam_nubp(:)                 = 0.0D0
              enbeam_nubp(:)                = 0.0D0
              enbeam_species_p(:,:)         = 0.0D0
              beam_beamddn_nubp(:)          = 0.0D0
              beam_beamdtn_nubp(:)          = 0.0D0
              beam_beamddp_nubp(:)          = 0.0D0
              beam_beamtt2n_nubp(:)         = 0.0D0 
              beam_thermaltth_df_nubp(:)    = 0.0D0 
              beam_thermalddp_nubp(:)       = 0.0D0 
              beam_thermalddn_nubp(:)       = 0.0D0 
              beam_thermaltt2n_nubp(:)      = 0.0D0 
              beam_thermaldth_tf_nubp(:)    = 0.0D0 

              sorbn0_nubp(:)                = 0.0D0 !jmp.den
              sorbh_nubp(:)                 = 0.0D0 !jmp.den

              qbeame_intg_nubp              = 0.0D0
              qbeami_intg_nubp              = 0.0D0
              qbth_intg_nubp                = 0.0D0
              curb_intg_nubp                = 0.0D0
              storqueb_intg_nubp            = 0.0D0
              sprbeame_intg_nubp            = 0.0D0
              sprbeami_intg_nubp            = 0.0D0
              wbeam_intg_nubp               = 0.0D0
              enbeam_intg_nubp              = 0.0D0
              pwf_tot_intg_nubp             = 0.0D0
              beam_thermal_dtntot_nubp      = 0.0D0
              beam_thermal_ddntot_nubp      = 0.0D0
              beam_thermal_ddptot_nubp      = 0.0D0
              beam_thermal_tt2ntot_nubp     = 0.0D0
              beam_beam_dtntot_nubp         = 0.0D0
              beam_beam_ddntot_nubp         = 0.0D0
              beam_beam_ddptot_nubp         = 0.0D0
              beam_beam_tt2ntot_nubp        = 0.0D0
          ELSE
              !start nubeam from a previously generated
              !restart file. In this case we have to load the
              !**nubp quantities with values that were active at the
              !restart time.
              CALL check_restart(err_state,err_xplasma)
              IF(err_state + err_xplasma .GT. 0)THEN
                 PRINT *,'check_restart reports :'
                 PRINT *,'err_state,err_xplasma + ', err_state,err_xplasma
                 CALL STOP('set_nubeam_init_profile',1)
              ENDIF
              CALL read_restart_profs 
              storqueb_nub(:) = storqueb_nubr
          ENDIF
          
      RETURN

      END SUBROUTINE set_nubeam_init_profs



      SUBROUTINE  read_restart_profs
!  -----------------------------------------------------------------------
!   read in  the profile information required in Onetwo to restart nubeam
!   using the nubeam generated restart files (eg. *.cdf files ,which are read by
!   the nubeam code) and the profile file, which is read here. 
!           
!  -----------------------------------------------------------HSJ---------
      USE transp,ONLY :                                                  &
                       qbeame_nubp,qbeami_nubp,                          &
                       qbeame_intg_nubp,qbeami_intg_nubp,                &
                       qbth_nubp,qbth_intg_nubp,curb_nubp,               &
                       storqueb_nubp,storqueb_intg_nubp,                 &
                       sprbeame_nubp,                                    &
                       sprbeami_nubp,wbeam_nubp,                         &
                       enbeam_nubp,sprbeami_intg_nubp,                   &
                       sprbeame_intg_nubp,wbeam_intg_nubp,               &
                       pwf_tot_intg_nub,pwf_tot_intg_nubp,               &
                       enbeam_intg_nubp,                                 &
                       beam_thermal_dtntot_nubp,                         &
                       beam_thermal_ddntot_nubp,                         &
                       beam_thermal_ddptot_nubp,                         &
                       beam_thermal_tt2ntot_nubp,                        &
                       beam_beam_dtntot_nubp,                            &
                       beam_beam_ddntot_nubp,                            &
                       beam_beam_ddptot_nubp,beam_beam_tt2ntot_nubp,     &
                       beam_thermaltth_df_nubp,beam_thermalddp_nubp,     &
                       beam_thermalddn_nubp,beam_thermaltt2n_nubp,       &
                       beam_thermaldth_tf_nubp,beam_beamddn_nubp,        &
                       beam_beamdtn_nubp,beam_beamddp_nubp,              &
                       beam_beamtt2n_nubp,curb_intg_nubp,                &
                       nubeam_restart,nubeam_state_path,nubeam_root,     &
                       nubeam_profile_path,nubeam_restart_time,          &
                       beam_data,nubeam_steps,nubeam_calls,              &
                       qbeame_nub,qbeami_nub,                            &
                       qbeame_intg_nub,qbeami_intg_nub,                  &
                       qbth_nub,qbth_intg_nub,curb_nub,                  &
                       storqueb_nub,storqueb_intg_nub,                   &
                       sprbeame_nub,                                     &
                       sprbeami_nub,wbeam_nub,                           &
                       enbeam_nub,sprbeami_intg_nub,                     &
                       sprbeame_intg_nub,wbeam_intg_nub,                 &
                       pwf_tot_intg_nub,pwf_tot_intg_nub,                &
                       enbeam_intg_nub,                                  &
                       beam_thermal_dtntot_nub,                          &
                       beam_thermal_ddntot_nub,                          &
                       beam_thermal_ddptot_nub,                          &
                       beam_thermal_tt2ntot_nub,                         &
                       beam_beam_dtntot_nub,                             &
                       beam_beam_ddntot_nub,                             &
                       beam_beam_ddptot_nub,beam_beam_tt2ntot_nub,       &
                       beam_thermaltth_df_nub,beam_thermalddp_nub,       &
                       beam_thermalddn_nub,beam_thermaltt2n_nub,         &
                       beam_thermaldth_tf_nub,beam_beamddn_nub,          &
                       beam_beamdtn_nub,beam_beamddp_nub,                &
                       beam_beamtt2n_nub,curb_intg_nub, &
                       sorbn0_nub,sorbn0_nubp,sorbh_nub,sorbh_nubp !jmp.den

      USE io, ONLY : ncrt,nout,nitre
      USE solcon, ONLY : time0

      USE restore_12
      IMPLICIT NONE

      REAL *8 transpt_time,nubeam_end_pen_time,nubeam_end_time
      REAL *8 adif,adifm
      INTEGER    io12,j,nbcs,nbcs2,nbcs3
      INTEGER js
      CHARACTER(len=256)  profile_file
      LOGICAL     ex
      CHARACTER(len=256) runid

       CALL getioun(io12,42)
       OPEN  (unit = io12,                                     &
            file = nubeam_profile_path,status = 'OLD',ERR = 3)




       READ (unit = io12, fmt = '(3x, a)') runid

       READ (unit = io12,fmt=1)transpt_time,nubeam_restart_time !jmp.ibm and subsequent similar changes
       READ(unit = io12,fmt=2)beam_data%nbeam,nubeam_steps,nbcs,nbcs2,nbcs3
       DO j=1,beam_data%nbeam
          READ (unit = io12,fmt=1)beam_data%pinja(j)
          READ (unit = io12,fmt=1)beam_data%einja(j)
          READ (unit = io12,fmt=1)beam_data%ffulla(j)
          READ (unit = io12,fmt=1)beam_data%fhalfa(j)
       ENDDO

!       IF(ALLOCATED(nubeam_calls))DEALLOCATE(nubeam_calls)
       IF(ASSOCIATED(nubeam_calls))DEALLOCATE(nubeam_calls)
       ALLOCATE(nubeam_calls(0:nbcs-1))
       READ (unit = io12,fmt=1)nubeam_calls


       !nubeam_calls(0:nubeam_steps-1) contains the
       !intervals over which nubeam was run. For example
       !nubeam_calls(1) - nubeam_calls(0)  is the first interval
       !of time that nubeam was run. Note that nubeam_calls(0) is
       !the very first  time that any of the beams is turned on.
       !It is always required that the first Onetwo/Nubeam run 
       !has time0 <= to this time. Thereafter, for subsequqnt runs
       !using the restart option, we build upon this array as is shown
       !below. It is required that in a restart run time0 is in the
       !last interval that nubeam was called. That is, 
       !nubeam_calls(nubeam_steps-2) .lt. time0 
       !and   time0 .le. nubeam_calls(nubeam_steps -1)
       nubeam_end_time = nubeam_calls(nubeam_steps -1)
       nubeam_end_pen_time = nubeam_calls(nubeam_steps -2)


!       IF(ALLOCATED(beam_data%beam_times_prev))                        &
       IF(ASSOCIATED(beam_data%beam_times_prev))                        &
                              DEALLOCATE(beam_data%beam_times_prev)
       ALLOCATE(beam_data%beam_times_prev(nbcs2))
       READ(unit = io12,fmt=1)beam_data%beam_times_prev




!       IF(ALLOCATED(beam_data%beam_switch_times_prev))                  &
       IF(ASSOCIATED(beam_data%beam_switch_times_prev))                  &
                        DEALLOCATE(beam_data%beam_switch_times_prev)
       ALLOCATE(beam_data%beam_switch_times_prev(nbcs3))
       READ(unit = io12,fmt=1)beam_data%beam_switch_times_prev

       !beam_data%beam_switch_times_prev contains all of the beam
       !turnon and turnoff times that were active in the run that 
       !created the profile restart file that we are about to read.
       ! Note that if the previous run
       !ended at the same time that  a new beam would have been turned on
       !then that beam tunon time time is not included in this list.
       !Similarly if a beam is turned off exactly at that time then it is
       !also not included in the list. These effects must be allowed for
       !in the restart file. 
       
       !check if  nubeam_end_time is abeam switching time (beam on or off)
      ! adifm = 100.*timmax
       DO j = 1,SIZE(beam_data%beam_switch_times)
          adif = ABS(beam_data%beam_switch_times(j) -  nubeam_end_time)
          adifm = MIN(adif,adifm)
          IF(adif .EQ. adifm) js =j
       ENDDO




       !if adifm is sufficiently small we assume that the run that
       !created the restart file was run up to the beginning or end
       !of abeam turn on or beam turn off time (or both). If adifm is
       !is larger than the threshold value then no beam switching 
       !was involved.


       
       



        PRINT *,'nubeam_steps =',nubeam_steps
!        IF(ALLOCATED(nubeam_calls))THEN
        IF(ASSOCIATED(nubeam_calls))THEN
           !JMP PRINT *,'nubeam_calls',nubeam_calls
        ELSE
           PRINT *,'nubeam_calls not allocated'
        ENDIF
        !JMP PRINT *,'beam_data%beam_times_prev',beam_data%beam_times_prev
        !JMP PRINT *,'beam_data%beam_times',beam_data%beam_times
        !JMP PRINT *,'beam_data%beam_switch_times_prev',beam_data%beam_switch_times_prev
        !JMP PRINT *,'beam_data%beam_switch_times',beam_data%beam_switch_times

        

       IF(time0 .LE. nubeam_end_pen_time .OR. time0 .GT. nubeam_end_time)THEN
          WRITE(nitre,100)nubeam_end_pen_time,nubeam_end_time,time0  
          WRITE(nout,100)nubeam_end_pen_time,nubeam_end_time,time0  
          WRITE(ncrt,100)nubeam_end_pen_time,nubeam_end_time,time0  
100       FORMAT(2x,'ERROR, start time  of this run must be greater than:',&
               1pe14.8,/, 2x,'and less than or equal to:',1pe14.8,/,       &
               2x,'But Input value (time0) is ',1pe14.8,/,                 &
               2x,'The fast ion profiles are interpolated in time in this',/,&
               2x,'interval.',/,           &
               2x,'But note that the initial condition for the thermal profiles',/,&
               2x, 'will be set to the values at the end time of the ',/,   &
               2x,'previous run even if time0 does not match that end time !')
          CALL STOP('read_restart_file',1)
       ENDIF




       READ (unit = io12,fmt=1)qbeame_nubp
       READ (unit = io12,fmt=1)qbeame_nub

       READ (unit = io12,fmt=1)qbeami_nubp 
       READ (unit = io12,fmt=1)qbeami_nub 
                           
       READ (unit = io12,fmt=1)qbth_nubp
       READ (unit = io12,fmt=1)qbth_nub
                  
       READ (unit = io12,fmt=1)curb_nubp 
       READ (unit = io12,fmt=1)curb_nub

       READ (unit = io12,fmt=1)storqueb_nubp                                     
       READ (unit = io12,fmt=1)storqueb_nub
        
       READ (unit = io12,fmt=1)sprbeame_nubp
       READ (unit = io12,fmt=1)sprbeame_nub

       READ (unit = io12,fmt=1)sprbeami_nubp
       READ (unit = io12,fmt=1)sprbeami_nubp

       READ (unit = io12,fmt=1)wbeam_nubp  
       READ (unit = io12,fmt=1)wbeam_nub 
            
       READ (unit = io12,fmt=1)enbeam_nubp  
       READ (unit = io12,fmt=1)enbeam_nub

       READ (unit = io12,fmt=1)beam_beamddn_nubp  
       READ (unit = io12,fmt=1)beam_beamddn_nub

     
       READ (unit = io12,fmt=1)beam_beamdtn_nubp    
       READ (unit = io12,fmt=1)beam_beamdtn_nub

 
       READ (unit = io12,fmt=1)beam_beamddp_nubp  
       READ (unit = io12,fmt=1)beam_beamddp_nub
  
       READ (unit = io12,fmt=1)beam_beamtt2n_nubp
       READ (unit = io12,fmt=1)beam_beamtt2n_nub
     
       READ (unit = io12,fmt=1)beam_thermaltth_df_nubp
       READ (unit = io12,fmt=1)beam_thermaltth_df_nub

       READ (unit = io12,fmt=1)beam_thermalddp_nubp 
       READ (unit = io12,fmt=1)beam_thermalddp_nub
  
       READ (unit = io12,fmt=1)beam_thermalddn_nubp  
       READ (unit = io12,fmt=1)beam_thermalddn_nub

       READ (unit = io12,fmt=1)beam_thermaltt2n_nubp 
       READ (unit = io12,fmt=1)beam_thermaltt2n_nub

       READ (unit = io12,fmt=1)beam_thermaldth_tf_nubp 
       READ (unit = io12,fmt=1)beam_thermaldth_tf_nub


       READ (unit = io12,fmt=1)qbeame_intg_nubp,qbeami_intg_nubp
       READ (unit = io12,fmt=1)qbeame_intg_nub, qbeami_intg_nub

       READ (unit = io12,fmt=1)qbth_intg_nubp,curb_intg_nubp
       READ (unit = io12,fmt=1)qbth_intg_nub, curb_intg_nub

       READ (unit = io12,fmt=1)storqueb_intg_nubp,pwf_tot_intg_nubp
       READ (unit = io12,fmt=1)storqueb_intg_nub, pwf_tot_intg_nub

       READ (unit = io12,fmt=1)sprbeame_intg_nubp,sprbeami_intg_nubp
       READ (unit = io12,fmt=1)sprbeame_intg_nub, sprbeami_intg_nub

       READ (unit = io12,fmt=1)wbeam_intg_nubp,enbeam_intg_nubp
       READ (unit = io12,fmt=1)wbeam_intg_nub, enbeam_intg_nub

       READ (unit = io12,fmt=1)beam_thermal_dtntot_nubp,beam_thermal_ddntot_nubp
       READ (unit = io12,fmt=1)beam_thermal_dtntot_nub,beam_thermal_ddntot_nub

       READ (unit = io12,fmt=1)beam_thermal_ddptot_nubp ,beam_thermal_tt2ntot_nubp
       READ (unit = io12,fmt=1)beam_thermal_ddptot_nub ,beam_thermal_tt2ntot_nub

       READ (unit = io12,fmt=1)beam_beam_dtntot_nubp,beam_beam_ddntot_nubp
       READ (unit = io12,fmt=1)beam_beam_dtntot_nub, beam_beam_ddntot_nub

       READ (unit = io12,fmt=1)beam_beam_ddptot_nubp,beam_beam_tt2ntot_nubp
       READ (unit = io12,fmt=1)beam_beam_ddptot_nub, beam_beam_tt2ntot_nub

       READ (unit = io12,fmt=1)storque_nub
       READ (unit = io12,fmt=1)storqueb_nubr
       READ (unit = io12,fmt=1)angrot_nub
       READ (unit = io12,fmt=1)te_nub
       READ (unit = io12,fmt=1)ti_nub
       READ (unit = io12,fmt=1)rbp_nub
       READ (unit = io12,fmt=1)ene_nub
       READ (unit = io12,fmt=1)en_nub
       READ (unit = io12,fmt=1)curden_nub
       READ (unit = io12,fmt=1)etor_nub
       READ (unit = io12,fmt=1)curdri_nub
       READ (unit = io12,fmt=1)curohm_nub
       READ (unit = io12,fmt=1)currf_nub 
       READ (unit = io12,fmt=1)curboot_nub
       READ (unit = io12,fmt=1)hcap_nub
       READ (unit = io12,fmt=1)gcap_nub
       READ (unit = io12,fmt=1)fcap_nub


       CLOSE (unit = io12)
       CALL giveupus(io12)

      

1      FORMAT(2x,1PE16.8,2x,1PE16.8)
2      FORMAT(6(2x,i5)) 
       RETURN



3       PRINT *, 'need this file for restart :',nubeam_profile_path
      CALL STOP('read_restart_profs: file not found',1) 
      END SUBROUTINE read_restart_profs








      SUBROUTINE  write_restart_profs
!  ----------------------------------------------------------------
!   write out the information required in Onetwo to restart nubeam
!   at a later time . This file is required in addition to the two
!   restart file created by nubeam .  The fully qualified file name is
!   ./nubeam_root//"_restart_profs.txt" 
! 
!  --------------------------------------------------------HSJ-----


      USE solcon,ONLY : time
      
      USE transp,ONLY :                                                  &
                       qbeame_nub,qbeami_nub,use_nubeam,                 &
                       qbeame_intg_nub,qbeami_intg_nub,                  &
                       qbth_nub,qbth_intg_nub,curb_nub,                  &
                       storqueb_nub,storqueb_intg_nub,                   &
                       sprbeame_nub,                                     &
                       sprbeami_nub,wbeam_nub,                           &
                       enbeam_nub,sprbeami_intg_nub,                     &
                       sprbeame_intg_nub,wbeam_intg_nub,                 &
                       pwf_tot_intg_nub,pwf_tot_intg_nub,                &
                       enbeam_intg_nub,                                  &
                       beam_thermal_dtntot_nub,                          &
                       beam_thermal_ddntot_nub,                          &
                       beam_thermal_ddptot_nub,                          &
                       beam_thermal_tt2ntot_nub,                         &
                       beam_beam_dtntot_nub,                             &
                       beam_beam_ddntot_nub,                             &
                       beam_beam_ddptot_nub,beam_beam_tt2ntot_nub,       &
                       beam_thermaltth_df_nub,beam_thermalddp_nub,       &
                       beam_thermalddn_nub,beam_thermaltt2n_nub,         &
                       beam_thermaldth_tf_nub,beam_beamddn_nub,          &
                       beam_beamdtn_nub,beam_beamddp_nub,                &
                       beam_beamtt2n_nub,curb_intg_nub,                  &
                       nubeam_restart,nubeam_state_path,nubeam_root,     &
                       nubeam_profile_path,wrt_restart_file_time,        &
                       nubeam_xplasma_path,profile_save_file,            &
                       state_file,xplasma_file,xplasma_save_file,        &
                       state_save_file,profile_file,nubeam_restart_time, &
                       beam_data, nubeam_steps,nubeam_calls,             &
                       qbeame_nubp,qbeami_nubp,                          &
                       qbeame_intg_nubp,qbeami_intg_nubp,                &
                       qbth_nubp,qbth_intg_nubp,curb_nubp,               &
                       storqueb_nubp,storqueb_intg_nubp,                 &
                       sprbeame_nubp,                                    &
                       sprbeami_nubp,wbeam_nubp,                         &
                       enbeam_nubp,sprbeami_intg_nubp,                   &
                       sprbeame_intg_nubp,wbeam_intg_nubp,               &
                       pwf_tot_intg_nubp,pwf_tot_intg_nubp,              &
                       enbeam_intg_nubp,                                 &
                       beam_thermal_dtntot_nubp,                         &
                       beam_thermal_ddntot_nubp,                         &
                       beam_thermal_ddptot_nubp,                         &
                       beam_thermal_tt2ntot_nubp,                        &
                       beam_beam_dtntot_nubp,                            &
                       beam_beam_ddntot_nubp,                            &
                       beam_beam_ddptot_nubp,beam_beam_tt2ntot_nubp,     &
                       beam_thermaltth_df_nubp,beam_thermalddp_nubp,     &
                       beam_thermalddn_nubp,beam_thermaltt2n_nubp,       &
                       beam_thermaldth_tf_nubp,beam_beamddn_nubp,        &
                       beam_beamdtn_nubp,beam_beamddp_nubp,              &
                       beam_beamtt2n_nubp,curb_intg_nubp, &
                       sorbn0_nub,sorbn0_nubp,sorbh_nub,sorbh_nubp !jmp.den

       USE restore_12
       USE tordlrot,    ONLY : storque,storqueb,angrot
       USE soln,        ONLY : te,ti,rbp,ene,en,curden,etor,curden
       USE sourc,       ONLY : curdri,curohm,currf,curboot
       USE geom,        ONLY : fcap,gcap,hcap

      IMPLICIT NONE
      INTEGER    io12,nbcs,nbcs2,nbcs3,j,ishell
      CHARACTER(len = 512) command
      CHARACTER(len=12) atime,intfl
      LOGICAL     ex


      print *,'use_nubeam =',use_nubeam
      IF(.NOT. use_nubeam)RETURN


      nbcs =0 ; nbcs2 =0 ; nbcs3 =0
!      if(ALLOCATED(nubeam_calls))nbcs = SIZE(nubeam_calls)
      IF(ASSOCIATED(nubeam_calls))nbcs = SIZE(nubeam_calls)
!      IF(ALLOCATED(beam_data%beam_times))nbcs2 = SIZE(beam_data%beam_times)
      IF(ASSOCIATED(beam_data%beam_times))nbcs2 = SIZE(beam_data%beam_times)
!      IF(ALLOCATED(beam_data%beam_switch_times))nbcs3 = SIZE(beam_data%beam_switch_times)
      IF(ASSOCIATED(beam_data%beam_switch_times))nbcs3 = SIZE(beam_data%beam_switch_times)

!      IF(nbcs*nbcs2*nbcs3 == 0) RETURN ! cant save until nubeam is called
     ! beam_data%beam_switch_times disabled (by JM ???)
     IF(nbcs*nbcs2 == 0) RETURN ! cant save until nubeam is called
      profile_file = './'//ADJUSTL(nubeam_root(1:LEN_TRIM(nubeam_root)))//       &
                   "_restart_profs.txt"
     

      INQUIRE(FILE = profile_file, EXIST  = ex)
      IF( ex ) CALL DESTROY (profile_file)
      CALL getioun(io12,42)
      OPEN  (unit = io12,                                           &
            file = profile_file(1:LEN_TRIM(profile_file)),     &
            status = 'NEW')
 
      Print *,'writing restart file :',profile_file(1:LEN_TRIM(profile_file))
       WRITE (unit = io12, fmt = '(3x, a)')'profile data from nubeam created by Onetwo'

       WRITE (unit = io12,fmt=1)time,nubeam_restart_time
       WRITE (unit = io12,fmt=2)beam_data%nbeam,nubeam_steps,nbcs,nbcs2,nbcs3
       DO j=1,beam_data%nbeam
          WRITE (unit = io12,fmt=1)beam_data%pinja(j)
          WRITE (unit = io12,fmt=1)beam_data%einja(j)
          WRITE (unit = io12,fmt=1)beam_data%ffulla(j)
          WRITE (unit = io12,fmt=1)beam_data%fhalfa(j)
       ENDDO
       WRITE (unit = io12,fmt=1)nubeam_calls
       WRITE (unit = io12,fmt=1)beam_data%beam_times

       ! added 12/8/2011 HSJ to avoid undefined 
       ! beam_data%beam_switch_times. Apparently these are not used in
       ! JMS's scheme snce he eliminated the creation of 
       ! beam_data%beam_switch_times!!??
       IF( .NOT. ASSOCIATED(beam_data%beam_switch_times))THEN
          ALLOCATE(beam_data%beam_switch_times(SIZE(beam_data%beam_times)))
          beam_data%beam_switch_times(:) = 0.0_dp
       ENDIF
       WRITE (unit = io12,fmt=1)beam_data%beam_switch_times

       WRITE (unit = io12,fmt=1)qbeame_nubp
       WRITE (unit = io12,fmt=1)qbeame_nub

       WRITE (unit = io12,fmt=1)qbeami_nubp
       WRITE (unit = io12,fmt=1)qbeami_nub

       WRITE (unit = io12,fmt=1)qbth_nubp
       WRITE (unit = io12,fmt=1)qbth_nub   

       WRITE (unit = io12,fmt=1)curb_nubp               
       WRITE (unit = io12,fmt=1)curb_nub 

       WRITE (unit = io12,fmt=1)storqueb_nubp
       WRITE (unit = io12,fmt=1)storqueb_nub                                        
       WRITE (unit = io12,fmt=1)sprbeame_nubp
       WRITE (unit = io12,fmt=1)sprbeame_nub

       WRITE (unit = io12,fmt=1)sprbeami_nubp
       WRITE (unit = io12,fmt=1)sprbeami_nub

       WRITE (unit = io12,fmt=1)wbeam_nubp
       WRITE (unit = io12,fmt=1)wbeam_nub

       WRITE (unit = io12,fmt=1)enbeam_nubp 
       WRITE (unit = io12,fmt=1)enbeam_nub
          
       WRITE (unit = io12,fmt=1)beam_beamddn_nubp                       
       WRITE (unit = io12,fmt=1)beam_beamddn_nub

       WRITE (unit = io12,fmt=1)beam_beamdtn_nubp
       WRITE (unit = io12,fmt=1)beam_beamdtn_nub  

       WRITE (unit = io12,fmt=1)beam_beamddp_nubp 
       WRITE (unit = io12,fmt=1)beam_beamddp_nub 

       WRITE (unit = io12,fmt=1)beam_beamtt2n_nubp  
       WRITE (unit = io12,fmt=1)beam_beamtt2n_nub

       WRITE (unit = io12,fmt=1)beam_thermaltth_df_nubp
       WRITE (unit = io12,fmt=1)beam_thermaltth_df_nub

       WRITE (unit = io12,fmt=1)beam_thermalddp_nubp   
       WRITE (unit = io12,fmt=1)beam_thermalddp_nub 

       WRITE (unit = io12,fmt=1)beam_thermalddn_nubp
       WRITE (unit = io12,fmt=1)beam_thermalddn_nub 

       WRITE (unit = io12,fmt=1)beam_thermaltt2n_nubp
       WRITE (unit = io12,fmt=1)beam_thermaltt2n_nub 


       WRITE (unit = io12,fmt=1)beam_thermaldth_tf_nubp
       WRITE (unit = io12,fmt=1)beam_thermaldth_tf_nub 


       WRITE (unit = io12,fmt=1)qbeame_intg_nubp,qbeami_intg_nubp
       WRITE (unit = io12,fmt=1)qbeame_intg_nub,qbeami_intg_nub

       WRITE (unit = io12,fmt=1)qbth_intg_nubp,curb_intg_nubp
       WRITE (unit = io12,fmt=1)qbth_intg_nub,curb_intg_nub

       WRITE (unit = io12,fmt=1)storqueb_intg_nubp,pwf_tot_intg_nubp
       WRITE (unit = io12,fmt=1)storqueb_intg_nub,pwf_tot_intg_nub

       WRITE (unit = io12,fmt=1)sprbeame_intg_nubp,sprbeami_intg_nubp
       WRITE (unit = io12,fmt=1)sprbeame_intg_nub,sprbeami_intg_nub

       WRITE (unit = io12,fmt=1)wbeam_intg_nubp,enbeam_intg_nubp
       WRITE (unit = io12,fmt=1)wbeam_intg_nub,enbeam_intg_nub

       WRITE (unit = io12,fmt=1)beam_thermal_dtntot_nubp,beam_thermal_ddntot_nubp
       WRITE (unit = io12,fmt=1)beam_thermal_dtntot_nub, beam_thermal_ddntot_nub

       WRITE (unit = io12,fmt=1)beam_thermal_ddptot_nubp, beam_thermal_tt2ntot_nubp
       WRITE (unit = io12,fmt=1)beam_thermal_ddptot_nub,beam_thermal_tt2ntot_nub

       WRITE (unit = io12,fmt=1)beam_beam_dtntot_nubp,beam_beam_ddntot_nubp
       WRITE (unit = io12,fmt=1)beam_beam_dtntot_nub, beam_beam_ddntot_nub

       WRITE (unit = io12,fmt=1)beam_beam_ddptot_nubp,beam_beam_tt2ntot_nubp
       WRITE (unit = io12,fmt=1)beam_beam_ddptot_nub, beam_beam_tt2ntot_nub


       WRITE (unit = io12,fmt=1)storque 
       WRITE (unit = io12,fmt=1)storqueb 
       WRITE (unit = io12,fmt=1)angrot 
       WRITE (unit = io12,fmt=1)te 
       WRITE (unit = io12,fmt=1)ti 
       WRITE (unit = io12,fmt=1)rbp 
       WRITE (unit = io12,fmt=1)ene 
       WRITE (unit = io12,fmt=1)en 
       WRITE (unit = io12,fmt=1)curden 
       WRITE (unit = io12,fmt=1)etor 
       WRITE (unit = io12,fmt=1)curdri 
       WRITE (unit = io12,fmt=1)curohm 
       WRITE (unit = io12,fmt=1)currf 
       WRITE (unit = io12,fmt=1)curboot 
       WRITE (unit = io12,fmt=1)hcap 
       WRITE (unit = io12,fmt=1)gcap 
       WRITE (unit = io12,fmt=1)fcap 

       CLOSE (unit = io12)
       CALL giveupus(io12)


       IF(time .LE. wrt_restart_file_time)THEN 
          !special section save files at time which is closest to
          !but not greater than wrt_restart_file_time. SPecial uses
          !ubclude possibly saving these files for mhd restart,etc.
          !first destroy old file (if it exists)
          INQUIRE(FILE = profile_save_file, EXIST  = ex)
          IF( ex ) CALL DESTROY (profile_save_file)
          !now create new file by copying current restart file :         
          WRITE(intfl,FMT='(1pe12.6)')time
          READ(intfl,FMT='(a)')atime
          profile_save_file = profile_file(1:LEN_TRIM(profile_file))//'_'//atime
          command = 'cp '//profile_file(1:LEN_TRIM(profile_file))//'  ' &
                        //profile_save_file(1:LEN_TRIM(profile_save_file))  
          IF (ISHELL (command) .LT. 0)THEN
           PRINT *,'Error, wrt_restart_file_time =',wrt_restart_file_time
           PRINT *,'did not create profile  file ',profile_file
          ENDIF
          
          !the two restart files (nubeam_root_nubeam_state.cdf,
          !nubeam_root_xplasma_state.cdf)  created by nubeam also 
          !have to be saved here:

          INQUIRE(FILE = state_save_file, EXIST  = ex)
          IF( ex ) CALL DESTROY (state_save_file)
          state_file = './'//ADJUSTL(nubeam_root(1:LEN_TRIM(nubeam_root)))//       &
                   "_nubeam_state.cdf"
          state_save_file =state_file(1:LEN_TRIM(state_file))//'_'//atime
          command = 'cp '//state_file(1:LEN_TRIM(state_file))//'  ' &
                        //state_save_file(1:LEN_TRIM(state_save_file)) 
          IF (ISHELL (command) .LT. 0)THEN
           PRINT *,'Error, wrt_restart_file_time =',wrt_restart_file_time
           PRINT *,'did not create state  file ' ,state_file
          ENDIF

 

          INQUIRE(FILE = xplasma_save_file, EXIST  = ex)
          IF( ex ) CALL DESTROY (xplasma_save_file)
           xplasma_file = './'//ADJUSTL(nubeam_root(1:LEN_TRIM(nubeam_root)))//       &
                   "_xplasma_state.cdf"
          xplasma_save_file =xplasma_file(1:LEN_TRIM(xplasma_file))//'_'//atime
          command = 'cp '//xplasma_file(1:LEN_TRIM(xplasma_file))//'  ' &
                        //xplasma_save_file(1:LEN_TRIM(xplasma_save_file))  
          IF (ISHELL (command) .LT. 0)THEN
           PRINT *,'Error, wrt_restart_file_time =',wrt_restart_file_time
           PRINT *,'did not create xplasma  file '
          ENDIF
       ENDIF


1      FORMAT(2x,1PE16.8,2x,1PE16.8)
2      FORMAT(6(2x,i5))
       RETURN
      END SUBROUTINE write_restart_profs



       SUBROUTINE update_nubeam_profs!
!  -------------------------------------------------------------------------


!  ---------------------------------------------------------------HSJ-------


      USE numbrs,ONLY :  nj,nion
      USE transp,ONLY :                                                  &
                       qbeame_nub,qbeami_nub,                            &
                       qbeame_intg_nub,qbeami_intg_nub,                  &
                       qbth_nub,qbth_intg_nub,curb_nub,                  &
                       storqueb_nub,storqueb_intg_nub,                   &
                       sprbeame_nub,                                     &
                       sprbeami_nub,wbeam_nub,                           &
                       enbeam_nub,sprbeami_intg_nub,                     &
                       sprbeame_intg_nub,wbeam_intg_nub,                 &
                       pwf_tot_intg_nub,pwf_tot_intg_nub,                &
                       enbeam_intg_nub,                                  &
                       beam_thermal_dtntot_nub,                          &
                       beam_thermal_ddntot_nub,                          &
                       beam_thermal_ddptot_nub,                          &
                       beam_thermal_tt2ntot_nub,                         &
                       beam_beam_dtntot_nub,                             &
                       beam_beam_ddntot_nub,                             &
                       beam_beam_ddptot_nub,beam_beam_tt2ntot_nub,       &
                       beam_thermaltth_df_nub,beam_thermalddp_nub,       &
                       beam_thermalddn_nub,beam_thermaltt2n_nub,         &
                       beam_thermaldth_tf_nub,beam_beamddn_nub,          &
                       beam_beamdtn_nub,beam_beamddp_nub,                &
                       beam_beamtt2n_nub,curb_intg_nub,                  &
                       qbeame_nubp,qbeami_nubp,                          &
                       qbeame_intg_nubp,qbeami_intg_nubp,                &
                       qbth_nubp,qbth_intg_nubp,curb_nubp,               &
                       storqueb_nubp,storqueb_intg_nubp,                 &
                       sprbeame_nubp,                                    &
                       sprbeami_nubp,wbeam_nubp,                         &
                       enbeam_nubp,sprbeami_intg_nubp,                   &
                       sprbeame_intg_nubp,wbeam_intg_nubp,               &
                       pwf_tot_intg_nub,pwf_tot_intg_nubp,               &
                       enbeam_intg_nubp,                                 &
                       beam_thermal_dtntot_nubp,                         &
                       beam_thermal_ddntot_nubp,                         &
                       beam_thermal_ddptot_nubp,                         &
                       beam_thermal_tt2ntot_nubp,                        &
                       beam_beam_dtntot_nubp,                            &
                       beam_beam_ddntot_nubp,                            &
                       beam_beam_ddptot_nubp,beam_beam_tt2ntot_nubp,     &
                       beam_thermaltth_df_nubp,beam_thermalddp_nubp,     &
                       beam_thermalddn_nubp,beam_thermaltt2n_nubp,       &
                       beam_thermaldth_tf_nubp,beam_beamddn_nubp,        &
                       beam_beamdtn_nubp,beam_beamddp_nubp,              &
                       beam_beamtt2n_nubp,curb_intg_nubp,                &
                       nubeam_restart,nubeam_state_path,nubeam_root,     &
                       nubeam_profile_path,enbeam_species,enbeam_species_p, &
                       nubeam_back_average, & !JMP
                       sorbn0_nub,sorbn0_nubp,sorbh_nub,sorbh_nubp !jmp.den

       USE tordlrot,    ONLY : storque,storqueb,angrot
       USE soln,        ONLY : te,ti,rbp,ene,en,curden,etor,curden
       USE sourc,       ONLY : curdri,curohm,currf,curboot
       USE param,       ONLY : kj !JMP
       USE solcon,      ONLY : time !JMP

       IMPLICIT NONE

!JMP START

       TYPE :: nubeam_data

         REAL*8 :: time
         REAL*8 :: curb(kj),wbeam(kj),enbeam(kj),qbeame(kj),qbeami(kj)
         TYPE(nubeam_data), POINTER :: next

       END TYPE

       TYPE(nubeam_data), pointer :: head
       TYPE(nubeam_data), pointer :: tail
       TYPE(nubeam_data), pointer :: p

       INTEGER :: count,istat

       !nubeam_back_average = 0.1

       IF (nubeam_back_average .gt. 0.0 ) THEN
  
         IF (.not. ASSOCIATED(head)) THEN

            ALLOCATE(head,STAT=istat)
            tail => head
            NULLIFY(tail%next)

            head%time = time
            head%curb(1:nj) = curb_nub(1:nj)
            head%wbeam(1:nj) = wbeam_nub(1:nj)
            head%enbeam(1:nj) = enbeam_nub(1:nj)
            head%qbeame(1:nj) = qbeame_nub(1:nj)
            head%qbeami(1:nj) = qbeami_nub(1:nj)
    
         ELSE
    
            ALLOCATE(tail%next,STAT=istat)
            tail => tail%next
            NULLIFY(tail%next)
    
            tail%time = time
            tail%curb(1:nj) = curb_nub(1:nj)
            tail%wbeam(1:nj) = wbeam_nub(1:nj)
            tail%enbeam(1:nj) = enbeam_nub(1:nj)
            tail%qbeame(1:nj) = qbeame_nub(1:nj)
            tail%qbeami(1:nj) = qbeami_nub(1:nj)
  
         ENDIF
  
         DO
            IF ( (tail%time-head%time) .gt. nubeam_back_average ) THEN
               p => head
               head => p%next
               DEALLOCATE(p,STAT=istat)
            ELSE
               EXIT
            ENDIF

         END DO  

         count = 0
         curb_nub(1:nj) = 0.0
         wbeam_nub(1:nj) = 0.0 
         enbeam_nub(1:nj) = 0.0
         qbeame_nub(1:nj) = 0.0
         qbeami_nub(1:nj) = 0.0
         p => head

         DO
           count = count + 1
           
           curb_nub(1:nj)   = curb_nub(1:nj)  + tail%curb(1:nj)
           wbeam_nub(1:nj)  = wbeam_nub(1:nj) + tail%wbeam(1:nj)
           enbeam_nub(1:nj) = enbeam_nub(1:nj)+ tail%enbeam(1:nj)
           qbeame_nub(1:nj) = qbeame_nub(1:nj)+ tail%qbeame(1:nj)
           qbeami_nub(1:nj) = qbeami_nub(1:nj)+ tail%qbeami(1:nj)
           
           !if ( p .eq. tail ) EXIT
           p => p%next
           IF(.NOT. ASSOCIATED(p)) EXIT
        
         END DO  

         curb_nub(1:nj) = curb_nub(1:nj) / count     
         wbeam_nub(1:nj) =  wbeam_nub(1:nj) / count   
         enbeam_nub(1:nj) = enbeam_nub(1:nj) / count      
         qbeame_nub(1:nj) = qbeame_nub(1:nj) / count 
         qbeami_nub(1:nj) = qbeami_nub(1:nj) / count
  
       ENDIF

!JMP END

                qbeame_nubp(1:nj)          = qbeame_nub(1:nj)
                qbeami_nubp(1:nj)          = qbeami_nub(1:nj)
                qbth_nubp(1:nj)            = qbth_nub(1:nj)
                curb_nubp(1:nj)            = curb_nub(1:nj)
                storqueb_nubp(1:nj)        = storqueb_nub(1:nj)
                sprbeame_nubp(1:nj)        = sprbeame_nub(1:nj)
                sprbeami_nubp(1:nj)        = sprbeami_nub(1:nj)
                wbeam_nubp(1:nj)           = wbeam_nub(1:nj)
                enbeam_nubp(1:nj)          = enbeam_nub(1:nj)
                enbeam_species_p(1:nj,:)   = enbeam_species(1:nj,:)
                beam_beamddn_nubp(:)       = beam_beamddn_nub(:)  
                beam_beamdtn_nubp(:)       = beam_beamdtn_nub(:)
                beam_beamddp_nubp(:)       = beam_beamddp_nub(:) 
                beam_beamtt2n_nubp(:)      = beam_beamtt2n_nub(:) 
                beam_thermaltth_df_nubp(:) = beam_thermaltth_df_nub(:)
                beam_thermalddp_nubp(:)    = beam_thermalddp_nub(:) 
                beam_thermalddn_nubp(:)    = beam_thermalddn_nub(:)
                beam_thermaltt2n_nubp(:)   = beam_thermaltt2n_nub(:)  
                beam_thermaldth_tf_nubp(:) = beam_thermaldth_tf_nub(:) 
 
                sorbn0_nubp(:) = sorbn0_nub(:)
                sorbh_nubp(:) = sorbh_nub(:)

! scalar beam quantities:
                qbeame_intg_nubp          = qbeame_intg_nub
                qbeami_intg_nubp          = qbeami_intg_nub
                qbth_intg_nubp            = qbth_intg_nub
                curb_intg_nubp            = curb_intg_nub
                storqueb_intg_nubp        = storqueb_intg_nub
                pwf_tot_intg_nubp         = pwf_tot_intg_nub
                sprbeame_intg_nubp        = sprbeame_intg_nub
                sprbeami_intg_nubp        = sprbeami_intg_nub
                wbeam_intg_nubp           = wbeam_intg_nub
                enbeam_intg_nubp          = enbeam_intg_nub
                beam_thermal_dtntot_nubp  = beam_thermal_dtntot_nub
                beam_thermal_ddntot_nubp  = beam_thermal_ddntot_nub
                beam_thermal_ddptot_nubp  = beam_thermal_ddptot_nub
                beam_thermal_tt2ntot_nubp = beam_thermal_tt2ntot_nub
                beam_beam_dtntot_nubp     = beam_beam_dtntot_nub
                beam_beam_ddntot_nubp     = beam_beam_ddntot_nub
                beam_beam_ddptot_nubp     = beam_beam_ddptot_nub
                beam_beam_tt2ntot_nubp    = beam_beam_tt2ntot_nub
 

      RETURN

      END SUBROUTINE update_nubeam_profs
