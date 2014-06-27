
   SUBROUTINE nblist_read(lun,nbname,ndifbe,ierr)
!  -------------------------------------------------------------------
!  THIS is not included in P_nfreya_inerface module because it 
!  would lead to a circular reference nbnamelist ==> uses P_nfreya_interface
!  read beam geometry data from beam_data_namelist
!  the data consist of up to 32 sources (each source is a beamline
!  in nubeam) The ordering is not unique apparently. Instead each
!  input "block"  completely specfies the beamline. beamlines that are
!  off are not output in the namelist file and the on  beams
!  are numbered contiguously from 1,2,3..  
!   we can identify
!  left sources by (tangency radii) rtcena(j) = 114.6cm and 
!  right sources by rtcena(j) = 76.2
!  There appears to be no information that allows us to distinguish
!  between beamlines however
!  NBLIST.FOR !   Note: 8 sources in order( HSJ what order?):
!  NBLIST.FOR !   30LEFT, 30RIGHT, 150LEFT, 150RIGHT,
!  NBLIST.FOR !   210LEFT, 210RIGHT, 330LEFT, 330RIGHT
!  ----------------------------------------------------HSJ-2-24-12----
   USE nrtype,                               ONLY : DP,i4B,I2B

   USE beam_structure,                       ONLY : neutral_beam

   USE error_handler,                        ONLY : lerrno,terminate,iomaxerr
        
   USE io_gcnmp,                             ONLY : ncrt,nlog

   USE common_constants,                     ONLY : izero,zeroc,im22icm2 
        

   ! originaly p_nfreya was designed with P_Nfreya_interface
   ! New overlay of nubeam 2011-11 requires this approach now: HSJ 
   USE P_nfreya_interface,                   ONLY :                                &
                                                    beam_data,                     &
                                                    dxbsmoo0       => dxbsmoo,     &
                                                    dtn_orbit,                     &
                                                    dn0out0        => dn0out,      &
                                                    d_fast_ion,                    &
                                                    fdifbe,                        &
                                                    edifbe,                        &
                                                    gflr_min0      => gflr_min,    &
                                                    goocon0        => goocon,      &
                                                    inzri,                         &
                                                    nsigexc0       => nsigexc,     &
                                                    nlminsv0       => nlminsv,     &
                                                    nlbbcx0        => nlbbcx,      &
                                                    nlebei0        => nlebei,      &
                                                    nptcls0        => nptcls,      &
                                                    nptclf0        => nptclf,      &
                                                    ndep00         => ndep0,       &
                                                    nmsigx,                        &
                                                    ndifbep,                       &
                                                    nubeam_restart,                &
                                                    nubeam_dt0     => nubeam_dt,   &
                                                    nubeam0_dt,                    &
                                                    nlbcde0        => nlbcde,      &
                                                    nlbcoh0        => nlbcoh,      &
                                                    nlbcpa0        => nlbcpa,      &
                                                    nlorbo0        => nlorbo,      &
                                                    nlbflr0        => nlbflr,      &
                                                    nlfbmflr0      => nlfbmflr,    &
                                                    nlbgflr0       => nlbgflr,     &
                                                    nlcprb0        => nlcprb,      &
                                                    nmcurb0        => nmcurb,      &
                                                    nkdifb0        => nkdifb,      &
                                                    nubeam_nclass,                 &
                                                    nlfhe3,                        &
                                                    nlfhe4,                        &
                                                    nlfst,                         &
                                                    nlfsp0         => nlusfp,      & 
                                                    nlusf30        => nlusf3,      &
                                                    nlusfa0        => nlusfa,      &
                                                    nlusfp0        => nlusfp,      &
                                                    nlusft0        => nlusft,      &
                                                    plfhe30        => plfhe3,      &
                                                    plfhe40        => plfhe4,      &
                                                    plfst0         => plfst,       &
                                                    plfsp0         => plfsp,       &
                                                    plhgt,                         &
                                                    wghta0         => wghta,       &
                                                    xdepmod0       => xdepmod,     &
                                                    xp_nt1,                        &
                                                    xcfanbi0       => xcfanbi,     &
                                                    xdfanbi0       => xdfanbi,     &
                                                    xefanbi0       => xefanbi,     &
                                                    xcfafus0       => xcfafus,     &  
                                                    xdfafus0       => xdfafus,     &
                                                    xefafus0       => xefafus,     &
                                                    xdatsfa0       => xdatsfa,     & 
                                                    xdatsf30       => xdatsf3,     &
                                                    xdatsft0       => xdatsft,     &
                                                    xdatsfp0       => xdatsfp


   USE fast_ion_diffusion,                   ONLY : check_d_fast_ion_inpt

   USE dnubeam_mod,                          ONLY : ntimemax,                                &
                                                    abeama,adiff_0,adiff_a,                  &
                                                    adiff_ntime,adiff_time,adiff_xpin,       &
                                                    adiff_xpout,aimp,asrd ,                  &
                                                    blk_mpi_bcast_int,blk_mpi_bcast_r8,      &
                                                    bmwidra,bmwidza,bsrd,cxpcon,             &
                                                    cxsplt,divra,divza,dn0out,dtn,dt_acc,    &
                                                    dxbsmoo,dxbsm_nc,                        &
                                                    edbfac,einja,erngfi,                     &
                                                    fbemax,fbemin,fbltim,                    &
                                                    fbtrap_depth,fdtnxy,ffulla,              &
                                                    fhalfa,foclra,foclza,                    &
                                                    fporcelli,fppcon,frac_depmax,            &
                                                    frac_depmin,frac_dep_lim,frac_orbrr,     &
                                                    fshper,fshwid,fvpvmn,                    &
                                                    fvpvmx,gflr_ll,gflr_min,                 &
                                                    gflr_op,gflr_rl,gflr_xv,                 &
                                                    goocon,levmod_halo,lev_nbidep,           &
                                                    lmidtube,lunnbx,lunres,                  &
                                                    mrstrt,nbapsh2,nbapsha,                  &
                                                    nbbcal,nbbcx_avg,nbbcx_bb,               &
                                                    nbbox,nbeam,nbebox,                      &
                                                    nbetube,nbfallgr,nbsbox,                 &
                                                    nbshapa,nbstube,nbtube,                  &
                                                    nchdvp,nclass,ncx0,                      &
                                                    ndep0,ndep0_max,ndepbox,                 &
                                                    ndeptube,ndep_set_beam,ndep_set_ien,     &
                                                    ndtorb,nerngfi,nfbon_species,            &
                                                    nfbon_vpvopt,ngoocon_ebin,ngoocon_rbin,  &
                                                    ngoocon_vpvbin,ngradd_opt,ngyro,         &
                                                    nkdifb,nlbbcx,nlbcde,                    &
                                                    nlbcoh,nlbcpa,                           &
                                                    nlbeamcx,nlbflr,nlbfpp,                  &
                                                    nlbgflr,nlbout,nlbox,                    &
                                                    nlco,nlcprb,nldep0_gather,               &
                                                    nlebei,nlfatom,nlfbmflr,                 &
                                                    nlfbon,nlfdep,nlhvion,                   &
                                                    nlminsv,nlorbo,nlpsirz,                  &
                                                    nlsawb,nlsawf,nlsym2b,                   &
                                                    nltest_output,nltrk_dep0,nlusf3,         &
                                                    nlusfa,nlusfp,nlusft,                    &
                                                    nmcurb,nmimp,nmix_kdsaw,                 &
                                                    nonlin,nper_cx,nptclf,                   &
                                                    nptclh,nptcls,nptcl_max,                 & 
                                                    nrip,nsdbgb,nseed,                       &
                                                    nsegtube,nsigexc,nth0,                   &
                                                    ntrace,nxbox,nybox,                      &
                                                    nznbma,nznbme,nzones,                    &
                                                    nzone_fb,only_io,orbrzv_option,          &
                                                    orbrzv_rf_option,orbrzv_toric_frnm,      &
                                                    orbrzv_toric_prftot0,                    &
                                                    orbrzv_rf_option,orbrzv_toric_frnm,      &
                                                    orbrzv_toric_prftot0,                    &
                                                    orbrzv_zzerr_con,phitube,pinja,          &
                                                    plfhe3,plfhe4,plfsp,                     &
                                                    plfst,quasi_check,rapedg2,               &
                                                    rapedga,ref_namelist,rhotube,            &
                                                    rtcena,rtube,sawflag,                    &
                                                    taurip,                                  &
                                                    tfshof,tfshon,thetatube,                 &
                                                    wghta,xboxhw,xbzeta,                     &
                                                    xcfafus,xcfanbi,xdatsf3,                 &
                                                    xdatsfa,xdatsfp,xdatsft,                 &
                                                    xdepmod,xdfafus,xdfanbi,                 &
                                                    xefafus,xefanbi,xfishmax,                &
                                                    xfishmin,xl1tube,xl2tube,                &
                                                    xlbapa,xlbapa2,xlbox1,                   &
                                                    xlbox2,xlbtna,xpand_nptcl,               &
                                                    xplasma_in_memory,xrapoff2,xrapoffa,     &
                                                    xswfrac_allfast,xswfrac_beam,            &
                                                    xswfrac_fusn,                            &
                                                    xybapa,xybsca,xzapoff2,                  &
                                                    xzapoffa,xzbeama,xzetatube,              &
                                                    xzimp,xzpedg2,xzpedga,                   &
                                                    yboxhw,ytube 
                                                    
                                                    
                                                    

                                                    
                                                    



   USE param,                                ONLY : kb,ke

   IMPLICIT NONE
   !TYPE(neutral_beam),INTENT(inout)  :: nbname
   TYPE(neutral_beam) nbname
   !INTEGER , INTENT (in )  :: lun
   INTEGER(I2b) lun
   !INTEGER , INTENT (out ) :: ierr
   INTEGER(I4B)  ierr
   INTEGER ll,j,k,kk
   INTEGER nimp,inta,intmax
   CHARACTER(Len=256) filename
   INTEGER iostat,ndifbe,ie
   LOGICAL fl_open,exists

   INTEGER,      PARAMETER         :: nbeamx = kb 
   REAL(DP),     DIMENSION(nbeamx) :: tbona, tboffa
!   INTEGER(I4B), DIMENSION(nbeamx) :: 
   REAL(DP)   :: nubeam_dt
   REAL(DP)   :: ebdmax
   LOGICAL    :: nlbdat
   INTEGER(i4B) nbnsz
 
   INCLUDE '../shared_modules/nbnamelist_post_201112.inc' ! takes place  of  USE nbnamelist


     ierr   = 0
     nbeam  = 0
     nimp   = 0
     intmax = 8  
     nlbdat =.FALSE.


!----------------------------------------------------------------------------------
! THEse FOUR  PARAMETERS DETERMINE IF UFILE MUST BE READ BY THE FACT THAT
! THEY ARE ZERO. (SEE USB UFILE_BYPASS)
! HERE WE DEFAULT THEM TO ZERO AND READ THE NAMELIST WHICH WILL
! GIVE PROPER VALUES IF THE NAMELIST FILE HAS THEM.
! OTHERWISE THE ZERO VALUES REMAIN UNCAHNGED AND WE WILL
! BE FORCED TO READ THE UFILE TO GET THIS INFORMATION
!----------------------------------------------------------------------------------
     pinja(:)  = zeroc
     einja(:)  = zeroc
     ffulla(:) = zeroc
     fhalfa(:) = zeroc

     xrapoffa(:)     = zeroc        
     xzapoffa(:)     = zeroc        
     xrapoff2(:)     = zeroc        
     xzapoff2(:)     = zeroc                      
     ntrace(:)       = zeroc
     adiff_time(:)   = zeroc
     adiff_ntime     = izero
     adiff_0(:)      = zeroc
     adiff_a(:)      = zeroc
     adiff_xpin(:)   = zeroc
     adiff_xpout(:)  = zeroc

     ! read namelist,typically in file NUBEAM_shot.DAT
     READ(lun,nml = nbdrive_naml)  ! ,END =15,ERR  =10) 
     !more informative to just let it fail



     !--------------------------------------------------------------
     ! set parameters in P_nfreya_interface to those in dnubeam_mod
     ! (some parameters are not defined in dnubeam_mod however)  
     !--------------------------------------------------------------
             dxbsmoo0       =  dxbsmoo
!            dtn_orbit
             dn0out0        =  dn0out
!            d_fast_ion
!            fdifbe
!            edifbe
             gflr_min0      =  gflr_min
             goocon0        =  goocon
!            inzri
             nsigexc0       =  nsigexc 
             nlminsv0       =  nlminsv 
             nlbbcx0        =  nlbbcx
             nlebei0        =  nlebei
             nptcls0        =  nptcls
             nptclf0        =  nptclf
             ndep00         =  ndep0
!            nmsigx
!            ndifbep
!            nubeam_restart
             nubeam_dt0     =  nubeam_dt
!            nubeam0_dt
             nlbcde0        =  nlbcde
             nlbcoh0        =  nlbcoh
             nlbcpa0        =  nlbcpa
             nlorbo0        =  nlorbo
             nlbflr0        =  nlbflr
             nlfbmflr0      =  nlfbmflr
             nlbgflr0       =  nlbgflr 
             nlcprb0        =  nlcprb
             nmcurb0        =  nmcurb
             nkdifb0        =  nkdifb
!            nubeam_nclass
!            nlfhe3
!            nlfhe4
!            nlfst
             nlfsp0         =  nlusfp 
             nlusf30        =  nlusf3
             nlusfa0        =  nlusfa
             nlusfp0        =  nlusfp  
             nlusft0        =  nlusft     
             plfhe30        =  plfhe3      
             plfhe40        =  plfhe4      
             plfst0         =  plfst       
             plfsp0         =  plfsp       
!            plhgt,                         
             wghta0         =  wghta       
             xdepmod0       =  xdepmod 
!            xp_nt1,                        
             xcfanbi0       =  xcfanbi 
             xdfanbi0       =  xdfanbi 
             xefanbi0       =  xefanbi 
             xcfafus0       =  xcfafus   
             xdfafus0       =  xdfafus 
             xefafus0       =  xefafus 
             xdatsfa0       =  xdatsfa  
             xdatsf30       =  xdatsf3 
             xdatsft0       =  xdatsft 
             xdatsfp0       =  xdatsfp

     !--------------------------------------------------------

     ! check input for expected consistency:
     !--------------------------------------------------------
!     CALL nbdrive_namel_check 


!    set those quantities read from the namelist that are used
!    by Onetwo and are stored in the nbname data structure:


     nbname%pinja(:)   = pinja(1:nbeam)
     nbname%einja(:)   = einja(1:nbeam)
     nbname%ffulla(:)  = ffulla(1:nbeam)
     nbname%fhalfa(:)  = fhalfa(1:nbeam) 
     nbname%abeama(:)  = abeama(1:nbeam)
     nbname%xzbeama(:) = xzbeama(1:nbeam)
     nbname%bmwidra(:) = bmwidra(1:nbeam)
     nbname%bmwidza(:) = bmwidza(1:nbeam) 
     nbname%rtcena(:)  = rtcena(1:nbeam)
     nbname%xlbtna(:)  = xlbtna(1:nbeam)
     nbname%xybsca(:)  = xybsca(1:nbeam)
     nbname%divza(:)   = divza(1:nbeam)
     nbname%divra(:)   = divra(1:nbeam) 
     nbname%foclza(:)  = foclza(1:nbeam)
     nbname%foclra(:)  = foclra(1:nbeam)
     nbname%rapedga(:) = rapedga(1:nbeam)
     nbname%xzpedga(:) = xzpedga(1:nbeam)
     nbname%xlbapa2(:) = xlbapa2(1:nbeam)
     nbname%xlbapa(:)  = xlbapa(1:nbeam)
     nbname%xybapa(:)  = xybapa(1:nbeam)
     nbname%rapedg2(:) = rapedg2(1:nbeam)
     nbname%xzpedg2(:) = xzpedg2(1:nbeam)
     nbname%xbzeta(:)  = xbzeta(1:nbeam) 

     nbname%ntrace(:)  = ntrace(1:nbeam)
     nbname%nbshapa(:) = nbshapa(1:nbeam)
     nbname%nbapsha(:) = nbapsha(1:nbeam)
     nbname%nlco(:)    = nlco(1:nbeam)
     nbname%nbeam      = nbeam
     nbname%nbbcal     = nbbcal
     nbname%nseed      = nseed 
     nubeam_nclass     = nclass
     nbname%nzone_nb   = nzones
     nbname%tbona(:)   = tbona(1:nbeam)
     nbname%tboffa(:)  = tboffa(1:nbeam)


     IF( .NOT. nbname%nlbdat_set)nbname%nlbdat      = nlbdat  ! use_ufile takes precedence

     nbnsz = SIZE(nbname%xrapoffa)
     ! added nbnsz 2/28/13 HSJ
     nbname%xrapoffa(:) = xrapoffa(1:nbnsz)       ! added 8/9/11 HSJ  
     nbname%xzapoffa(:) = xzapoffa(1:nbnsz)       ! added 8/9/11 HSJ  
     nbname%xrapoff2(:) = xrapoff2(1:nbnsz)       ! added 8/9/11 HSJ  
     nbname%xzapoff2(:) = xzapoff2(1:nbnsz)       ! added 8/9/11 HSJ  



     nubeam0_dt        = nubeam_dt 
     !----------------------------------------------------------------------------------
     ! in dnubeam_mod.f90 adiff_a(1:ntimemax),etc are fucntions of time only.
     ! In P_nfreya adiff_a(1:ntimemax,1:ke) are allowed to be fucntion of injection
     ! enrgy as well but at present there is no way to set the energy values
     ! differently
     !----------------------------------------------------------------------------------HSJ
     DO ie=1,ke
        d_fast_ion%adiff_a(:,ie)       = adiff_a(:) * im22icm2   ! namelist input is in cm2/sec
        d_fast_ion%adiff_0(:,ie)       = adiff_0(:) * im22icm2   ! we want m2/sec for use in rep_transport_output.f90
        d_fast_ion%adiff_xpin(:,ie)    = adiff_xpin
        d_fast_ion%adiff_xpout(:,ie)   = adiff_xpout
     ENDDO
     d_fast_ion%adiff_time(:)         = adiff_time(:)
     d_fast_ion%adiff_ntime           = adiff_ntime
     d_fast_ion%nkdifb                = nkdifb
     d_fast_ion%ndifbe                = ndifbe 

     CALL check_d_fast_ion_inpt(d_fast_ion,fdifbe,edifbe,ndifbep)


!     nptcls  these are also read from the namelist 
!             but they are stored in tranp.mod  rather than in nbname
!     nptclf
!     ndep0

!     finally there are some quantities that are accepted by the
!     namelist read for purposes of compatibility
!     only. These quantities are not sued in Onetwo at all.
!     Examples are frac_ions,frac)imp, etc.


     IF(nbeam .GT. nbeamx)THEN
        lerrno = 258 + iomaxerr
        CALL terminate(lerrno,nlog)
     ENDIF

   RETURN


10   ierr =1
     PRINT *,'"Error in reading file "',filename(1:LEN_TRIM(filename)),'"'
   RETURN
15   ierr = 1
     PRINT *,"Ecountered end  file without finding namelist ",filename(1:LEN_TRIM(filename)),'"'
   END SUBROUTINE nblist_read



   SUBROUTINE nblist_read_parse(lun,nbname,ierr)
!  -------------------------------------------------------------------
!  read beam geometry data from beam_data_namelist
!  the data consist of up to 14 sources (each source is a beamline
!  in nubeam) The ordering is not unique apparently. Instead each
!  input "block"  completely specfies the beamline. beamlines that are
!  off are not output in the namelist file and the ona beams
!  are numbered contiguously from 1,2,3..  
!   we can identify
!  left sources by (tangency radii) rtcena(j) = 114.6cm and 
!  right sources by rtcena(j) = 76.2
!  There appears to be no information that allows us to distinguish
!  between beamlines however
!  NBLIST.FOR !   Note: 8 sources in order:
!  NBLIST.FOR !   30LEFT, 30RIGHT, 150LEFT, 150RIGHT,
!  NBLIST.FOR !   210LEFT, 210RIGHT, 330LEFT, 330RIGHT

!  ----------------------------------------------------HSJ-11/18/03----

   USE nrtype,            ONLY : DP,i4B,I2B
   USE beam_structure,    ONLY : neutral_beam

   USE string_util,       ONLY : strip_path,to_upper_case1,to_lower_case1,       &
                                 set_r8_1darray_value,set_i4_1darray_value,      &
                                 set_L_1darray_value,set_r4_value,               &
                                 set_integer_value, set_r4_1darray_value,        &
                                 set_r8_value,set_L_value

   USE error_handler,     ONLY : lerrno,terminate,iomaxerr

   USE io_gcnmp,          ONLY : ncrt,nlog

   USE P_nfreya_interface, ONLY : nptcls,nptclf,ndep0,nubeam_nclass,nubeam_dt,   &
                                  d_fast_ion,nubeam0_dt,ndifbep,fdifbe,edifbe

   USE fast_ion_diffusion, ONLY : check_d_fast_ion_inpt

   IMPLICIT NONE
   !TYPE(neutral_beam),INTENT(out)  :: nbname
   TYPE(neutral_beam)  nbname
   !INTEGER , INTENT (in ) :: lun
   INTEGER(I2B)  lun
   !INTEGER , INTENT (out ) :: ierr
   INTEGER(I4B)  ierr
   INTEGER ll,j,k,kk,ndifbe
   INTEGER nbeam,nimp,inta,intmax
   CHARACTER(len=256) :: line,line_prev
   LOGICAL no_cont,last_line
   REAL(DP) adiff0_l,adiffa_l,xpin_l,xpout_l

   




   !the file actually does not define a formal namelist so we have to 
   !punt to fetch out the data:
   nbeam = 0   ; d_fast_ion%ndifbe =0
   nimp = 0
   intmax = 8  
   ierr = 0
   line_prev(1:1) = '!'
   last_line = .FALSE.

   DO WHILE(1 .GT.  0)
     no_cont = .TRUE.
     line = line_prev

     IF( .NOT. last_line)  THEN                            
        READ(lun,FMT='(A)',ERR = 15,END=10)line_prev
        go to 11
 10     last_line = .TRUE.
        !check if line just read is a continuation of previous line:
 11     line_prev  = ADJUSTL(line_prev)
        !line_prev = to_upper_case(line_prev)
        CALL to_upper_case(line_prev)
        ll = LEN_TRIM(line_prev) 
        DO j=1,ll
           IF(line_prev(j:j) == ' ')CYCLE
           IF(line_prev(j:j) == '!')EXIT
           IF(line_prev(j:j) == ';')EXIT
           !first character on line must be '=','+','-' or digit
           !otherwise this is not a continuation line
           IF(line_prev(j:j) == '=' .OR.            &
              line_prev(j:j) == '+' .OR.            &
              line_prev(j:j) == '-' .OR.            &
              line_prev(j:j) == '0' .OR.            &
              line_prev(j:j) == '1' .OR.            &
              line_prev(j:j) == '2' .OR.            &
              line_prev(j:j) == '3' .OR.            &
              line_prev(j:j) == '4' .OR.            &
              line_prev(j:j) == '5' .OR.            &
              line_prev(j:j) == '6' .OR.            &
              line_prev(j:j) == '7' .OR.            &
              line_prev(j:j) == '8' .OR.            &
              line_prev(j:j) == '9')THEN
                  no_cont = .FALSE.
                  PRINT *,'lie_prev(j:j) =',line_prev(j:j)
                  PRINT *,'NOT SET UP TO READ CONTINUATION LINES'
                  PRINT *,'Line encountered is:'
                  PRINT*,line_prev(1:ll)
                  lerrno = 259 + iomaxerr
                  CALL terminate(lerrno,nlog)
            ELSE
               EXIT
            ENDIF
        ENDDO
     ENDIF


     line  = ADJUSTL(line)
     !line = to_upper_case(line)
     CALL to_upper_case(line)
     ll = LEN_TRIM(line)                    
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
!     IF(k == 0)k=ll                                 !string does not have '!'
     IF(k ==0 ) THEN
         k =ll
     ELSE
         k = k-1
     ENDIF


     kk= INDEX(line(1:k),'ABEAMA')
     IF(kk .NE. 0)                                            &        
        CALL set_r8_1darray_value(k,kk,line,nbname%abeama)

     kk= INDEX(line(1:k),'XZBEAMA')
     IF(kk .NE. 0)                                             &
        CALL set_r8_1darray_value(k,kk,line,nbname%xzbeama)

     kk= INDEX(line(1:k),'APLASM')
     IF(kk .NE. 0)                                            &   
        CALL set_r8_1darray_value(k,kk,line,nbname%aplasm)

     kk= INDEX(line(1:k),'BACKZ')
     IF(kk .NE. 0)                                            &   
        CALL set_r8_1darray_value(k,kk,line,nbname%backz)


     kk= INDEX(line(1:k),'BMWIDRA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%bmwidra)



     kk= INDEX(line(1:k),'BMWIDZA')
     IF(kk .NE. 0)                                            &
          CALL set_r8_1darray_value(k,kk,line,nbname%bmwidza)



     kk= INDEX(line(1:k),'RTCENA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%rtcena)


     kk= INDEX(line(1:k),'XLBTNA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xlbtna)


     kk= INDEX(line(1:k),'XYBSCA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xybsca)


     kk= INDEX(line(1:k),'DIVZA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%divza)

     kk= INDEX(line(1:k),'DIVRA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%divra)


     kk= INDEX(line(1:k),'FOCLZA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%foclza)


     kk= INDEX(line(1:k),'FOCLRA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%foclra)



     kk= INDEX(line(1:k),'RAPEDGA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%rapedga)

     kk= INDEX(line(1:k),'XZPEDGA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xzpedga)

    kk= INDEX(line(1:k),'XRAPOFFA') ! HSJ 8/9/11
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xrapoffa)

    kk= INDEX(line(1:k),'XZAPOFFA') ! HSJ 8/9/11
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xzapoffa)


    kk= INDEX(line(1:k),'XRAPOFF2') ! HSJ 8/9/11
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xrapoff2)

    kk= INDEX(line(1:k),'XZAPOFF2') ! HSJ 8/9/11
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xzapoff2)



     kk= INDEX(line(1:k),'XLBAPA2')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xlbapa2)



     kk= INDEX(line(1:k),'XLBAPA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xlbapa)



     kk= INDEX(line(1:k),'XYBAPA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xybapa)



     kk= INDEX(line(1:k),'RAPEDG2')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%rapedg2)



     kk= INDEX(line(1:k),'XZPEDG2')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xzpedg2)


!     kk= INDEX(line(1:k),'XZIMPX')
!     IF(kk .NE. 0)                                            &
!         CALL set_r8_1darray_value(k,kk,line,nbname%xzimpx)



!     kk= INDEX(line(1:k),'AIMPX')
!     IF(kk .NE. 0)                                            &
!         CALL set_r8_1darray_value(k,kk,line,nbname%aimpx)


     kk= INDEX(line(1:k),'XBZETA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%xbzeta)

     kk= INDEX(line(1:k),'PINJA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%pinja)


     kk= INDEX(line(1:k),'EINJA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%einja)


     kk= INDEX(line(1:k),'FFULLA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%ffulla)

     kk= INDEX(line(1:k),'FHALFA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%fhalfa)


     kk= INDEX(line(1:k),'NTRACE')
     IF(kk .NE. 0)                                            &
         CALL set_i4_1darray_value(k,kk,line,nbname%ntrace)



     kk= INDEX(line(1:k),'NPTCLS')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nptcls)



     kk= INDEX(line(1:k),'NPTCLF')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nptclf)


     kk= INDEX(line(1:k),'NDEP0')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,ndep0)


     kk= INDEX(line(1:k),'NBSHAPA')
     IF(kk .NE. 0)                                            &
         CALL set_i4_1darray_value(k,kk,line,nbname%nbshapa)

    kk= INDEX(line(1:k),'NBAPSHA')
     IF(kk .NE. 0)                                            &
         CALL set_i4_1darray_value(k,kk,line,nbname%nbapsha)


    kk= INDEX(line(1:k),'NLCO')
     IF(kk .NE. 0)                                            &
         CALL set_L_1darray_value(k,kk,line,nbname%nlco)



    kk= INDEX(line(1:k),'NLBDAT')
     IF(kk .NE. 0)    THEN     
                                    
          CALL set_L_value(k,kk,line,nbname%nlbdat)

     ENDIF


    kk= INDEX(line(1:k),'NBEAM')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nbname%nbeam)

!    kk= INDEX(line(1:k),'NGMAX')
!     IF(kk .NE. 0)                                            &
!         CALL set_integer_value(k,kk,line,nbname%ngmax)

!    kk= INDEX(line(1:k),'NRHIX')
!     IF(kk .NE. 0)                                            &
!         CALL set_integer_value(k,kk,line,nbname%nrhix)


    kk= INDEX(line(1:k),'NBBCAL')
     IF(kk .NE. 0)                                            & 
         CALL set_integer_value(k,kk,line,nbname%nbbcal)
 
    kk= INDEX(line(1:k),'NSEED')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nbname%nseed)
 
    kk= INDEX(line(1:k),'NCLASS')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nubeam_nclass)


 
    kk= INDEX(line(1:k),'NZONES')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,nbname%nzone_nb)

    kk= INDEX(line(1:k),'TBONA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%tbona)
 


    kk= INDEX(line(1:k),'TBOFFA')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,nbname%tboffa)

 

    kk= INDEX(line(1:k),'NUBEAM_DT')
     IF(kk .NE. 0)                                            &
         CALL set_r8_value(k,kk,line,nubeam0_dt)



    kk= INDEX(line(1:k),'ADIFF_0')
     IF(kk .NE. 0)THEN                                            
         !CALL set_r8_value(k,kk,line,d_fast_ion%adiff_0)
         CALL set_r8_value(k,kk,line,adiff0_l)
         d_fast_ion%adiff_0(:,:) = adiff0_l
     ENDIF

    kk= INDEX(line(1:k),'ADIFF_A')
     IF(kk .NE. 0)THEN                                            
         !CALL set_r8_value(k,kk,line,d_fast_ion%adiff_a)
         CALL set_r8_value(k,kk,line,adiffa_l)
         d_fast_ion%adiff_a(:,:) = adiffa_l
     ENDIF



    kk= INDEX(line(1:k),'ADIFF_XPIN')
     IF(kk .NE. 0) THEN                                           
         !CALL set_r8_value(k,kk,line,d_fast_ion%adiff_xpin)
         CALL set_r8_value(k,kk,line,xpin_l)
         d_fast_ion%adiff_xpin(:,:) = xpin_l
     ENDIF

    kk= INDEX(line(1:k),'ADIFF_XPOUT')
     IF(kk .NE. 0)THEN                                            
         !CALL set_r8_value(k,kk,line,d_fast_ion%adiff_xpout)
         CALL set_r8_value(k,kk,line,xpout_l)
         d_fast_ion%adiff_xpout(:,:) = xpout_l
     ENDIF

    kk= INDEX(line(1:k),'NDIFBE')
     IF(kk .NE. 0)                                            & 
         CALL set_integer_value(k,kk,line,d_fast_ion%ndifbe)

    kk= INDEX(line(1:k),'NKDIFB')
     IF(kk .NE. 0)                                            &
         CALL set_integer_value(k,kk,line,d_fast_ion%nkdifb)

     kk= INDEX(line(1:k),'FDIFBE')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,fdifbe)

     kk= INDEX(line(1:k),'EDIFBE')
     IF(kk .NE. 0)                                            &
         CALL set_r8_1darray_value(k,kk,line,edifbe)

    IF(last_line)EXIT
   END DO
   
   CALL check_d_fast_ion_inpt(d_fast_ion,fdifbe,edifbe,ndifbep)


   RETURN

15 ierr = 1
   RETURN

   END SUBROUTINE nblist_read_parse

