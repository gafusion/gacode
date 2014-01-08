!----------------------------------------------------------------------
!   d3d nubeam driver/interface for onetwo
!   last modified, jul 2011
!   JM
!   Feb 2012 HSJ

SUBROUTINE dnubeam_driver

  USE dnubeam_mod
  USE numbrs, ONLY : nj,nion
  USE solcon, ONLY : time,dt
  USE transp, ONLY : beam_data,nubeam_dt,use_ufile
  USE nub,    ONLY : pbeam,ebeam
  USE nub2,   ONLY : bion,bneut
  USE nub4,   ONLY : vbeam
  USE io,     ONLY : ncrt,nout
  USE ext_prog_info, ONLY : get_nubeam,nubeam_path

  IMPLICIT NONE

  REAL*8  :: nubeam_dt_step

  INTEGER :: ierr,ishell,len_str,k,k0
  REAL*8  :: x,xxx

  CHARACTER(len=256) :: nubeam_path_out,nubeam_setup_out
  CHARACTER(len=256) :: command
  CHARACTER(len=24)  :: nubeam_namelist

  LOGICAL :: first_time
  DATA first_time /.TRUE./

  ! --------------------------------------------------------------------
  ! entry

  WRITE(*, fmt = '("\n\n==========================================")' )
  WRITE(*, fmt = '("+ entering call_d3d_nubeam_driver")' )

  ! --------------------------------------------------------------------
  ! get nubeam driver path

  len_str = LEN_TRIM(nubeam_path)
  IF (first_time) THEN

     WRITE(*, fmt = '("+ nubeam_path: ",a)' ),TRIM(nubeam_path)
     CALL get_nubeam(ncrt,nout,nubeam_path_out,nubeam_setup_out,len_str)
     IF(len_str .LE. 0)THEN
         PRINT *,'sub get_nubeam did not find required path'
         PRINT *,'nubeam_path_out = ',nubeam_path_out
         CALL STOP('nubeam path problem',1)
     ENDIF

  ENDIF

  ! ------------------------------------------------------------------
  ! initialize

  IF (first_time) THEN

     t0_nubeam = time
     t1_nubeam = time

     CALL init_nubeam_data()

  ENDIF

  ! ------------------------------------------------------------------
  ! load nubeam output from previous nubeam call if t <= t1_nubeam

  IF ( (time .LE. t1_nubeam) .AND. first_time .EQ. .FALSE.) THEN
      PRINT *, 'load_nubeam from the previous nubeam run'
      CALL load_nubeam(time)
      first_time = .FALSE.
      RETURN
  ENDIF

  ! ------------------------------------------------------------------
  ! nubeam_dt_step

  IF (use_ufile) THEN
    xxx = 1.e30
    DO k=1,SIZE(beam_data%beam_times)
      x= beam_data%beam_times(k)-t1_nubeam
      IF(x .GT.0.0)THEN
         k0 = k
         EXIT 
      ENDIF
    ENDDO
    nubeam_dt_step = beam_data%beam_times(k0)-beam_data%beam_times(k0-1)
  ELSE
    nubeam_dt_step = nubeam_dt
  ENDIF

  t0_nubeam = t1_nubeam
  t1_nubeam = t0_nubeam+nubeam_dt_step

  WRITE(*, fmt = '("+ nubeam_dt:",f8.4)' ) nubeam_dt
  WRITE(*, fmt = '("+ nubeam_dt_step:",f8.4)' ) nubeam_dt_step
  WRITE(*, fmt = '("+ time,t0,t1:",3f8.4,"\n")' ) time, t0_nubeam, t1_nubeam

  ! ------------------------------------------------------------------
  ! write nubeam input

  CALL write_nubeam_input(TRIM(nubeam_namelist),time)

  ! ------------------------------------------------------------------
  ! call ubeam driver

  IF (first_time) THEN
    command = TRIM(nubeam_path_out)//' '//'init'
  ELSE
    command = TRIM(nubeam_path_out)//' '//'step 1'
  ENDIF

  WRITE(*, fmt = '("\n\n==========================================")' )
  WRITE(*, fmt = '("+ nubeam started")' )
  WRITE(*, fmt = '("+ running : ",a,"\n\n")'   ) TRIM(command)
  ierr = ishell(command) 
  IF(ierr .LT. 0)THEN      
    WRITE(6, fmt = '("error in nubeam",i)' ) ierr
    CALL STOP ('subroutine nubeam_driver ', 501)
  ENDIF

  ! ------------------------------------------------------------------
  ! read output

  CALL read_nubeam_output()

  ! ------------------------------------------------------------------
  ! load nubeam output to onetwo

  PRINT *,'load nubeam' 
  CALL load_nubeam(time)

  ! ------------------------------------------------------------------
  ! return

  first_time = .FALSE.

  PRINT *, 'leaving dnubeam'

  RETURN

END

SUBROUTINE write_nubeam_input(nubeam_namelist)

  USE dnubeam_mod

  USE constnts,   ONLY : pi

  USE param,      ONLY : kj,kprim,kimp
  USE numbrs,     ONLY : nj,nprim,nimp 

  USE transp,     ONLY : beam_data

  USE ename ,     ONLY : eqdskfilename,eqfile,eqdsk_tdem

  USE soln ,      ONLY : etor,ene,te,ti,en !en(kj,kion)
  USE ions ,      ONLY : zeff
  USE tordlrot,   ONLY : iangrot, angrot
  USE machin,     ONLY : rmajor,btor
  USE neut ,      ONLY : ennw,ennv,tn,fluxn,in
  USE nub,        ONLY : nbeams,sbcx,ibion
  USE geom,       ONLY : sfarea

  USE solcon,      ONLY : time0

  IMPLICIT NONE

  CHARACTER*(*), INTENT(in) :: nubeam_namelist

  CHARACTER(len=64)  :: eqdsk_name
  CHARACTER(len=256) :: command
  REAL*8 time_eqdsk,time_pwr
  INTEGER :: ishell,ierr,j,len_str

  ! ------------------------------------------------------------------
  ! beam power, energy, species mix
 
  time_pwr = 0.5*(t0_nubeam+t1_nubeam) !time + 0.5*nubeam_dt

 
  CALL beam_power_interp(time_pwr)

  beam_data%pwf_tot_intg = 0.0
  DO j=1,beam_data%nbeam 
      pinja (j) = beam_data%pinja(j)
      einja (j) = beam_data%einja (j)
      ffulla(j) = beam_data%ffulla(j)
      fhalfa(j) = beam_data%fhalfa(j)
      beam_data%pwf_tot_intg = beam_data%pwf_tot_intg + beam_data%pinja(j)
  END DO

  WRITE(*, fmt = '("\n\n==========================================")' )
  WRITE(*, fmt = '("+ time_pwr",f8.4)' ),time_pwr
  PRINT *,t0_nubeam,t1_nubeam
  PRINT *,beam_data%pwf_tot_intg

  ! ------------------------------------------------------------------
  ! time dependent anomalous diffusion

  IF ( time_pwr .LE. adiff_time(1) ) THEN
    difb_0   = adiff_0(1)
    difb_0   = adiff_0(1)
    difb_in  = adiff_xpin(1)
    difb_out = adiff_xpout(1)
  ELSE IF ( time_pwr .GE. adiff_time(adiff_ntime) ) THEN
    difb_0   = adiff_0(adiff_ntime)
    difb_a   = adiff_a(adiff_ntime)
    difb_in  = adiff_xpin(adiff_ntime)
    difb_out = adiff_xpout(adiff_ntime)
  ELSE 
    DO j=1,adiff_ntime-1
      IF ( time_pwr .GE. adiff_time(j) .AND. time_pwr .LE. adiff_time(j+1) ) THEN
        CALL interp_nubeam_0d(difb_0  ,adiff_0    (j:j+1),time_pwr,adiff_time(j),adiff_time(j+1))
        CALL interp_nubeam_0d(difb_a  ,adiff_a    (j:j+1),time_pwr,adiff_time(j),adiff_time(j+1))
        CALL interp_nubeam_0d(difb_in ,adiff_xpin (j:j+1),time_pwr,adiff_time(j),adiff_time(j+1))
        CALL interp_nubeam_0d(difb_out,adiff_xpout(j:j+1),time_pwr,adiff_time(j),adiff_time(j+1))
        EXIT
      ENDIF
    ENDDO
  ENDIF

  PRINT *,'diff:',difb_0,difb_a,difb_in,difb_out

  ! ------------------------------------------------------------------
  ! eqdsk

  eqfile = eqdskfilename
  IF(eqdsk_tdem .NE. 'tdem' ) THEN
     len_str = LEN_TRIM(eqdskfilename)
     eqdsk_name  = eqdskfilename(1:len_str)
  ELSE
     IF (t0_nubeam .LT. time0) THEN
        CALL wrt_tdem_eqdsk(time0,eqdsk_name)
     ELSE
        CALL wrt_tdem_eqdsk(t0_nubeam,eqdsk_name)
     ENDIF
     eqfile = eqdsk_name
     WRITE(*, fmt = '("+ eqdsk file created from tdem input")' )
  ENDIF
  WRITE(*, fmt = '("+ eqdsk_tdem: ",a)' ),eqdsk_tdem
  WRITE(*, fmt = '("+ eqdsk_name: ",a)' ),eqdsk_name
  WRITE(*, fmt = '("+ eqfile    : ",a)' ),eqfile

  ! ------------------------------------------------------------------
  ! dnubeam_naml

  IF ( nseed .EQ. 0 ) nseed    = 1094088621

  ! ------------------------------------------------------------------
  ! dnubeam_prof

   nrho_in = nj
   DO j = 1,nj
      rho_in(j)  = (j-1.0)/(nj-1.0)
   ENDDO
   zeff_in(1:nj) = zeff(1:nj)
   te_in(1:nj) = te(1:nj) 
   ti_in(1:nj) = ti(1:nj) 
   ne_in(1:nj) = ene(1:nj)*1.0e6
  !ni_in(1:nj) = en(1:nj,1) !en(:,1:nprim)
  !n0_in(1:nj) = (ennw(1:nj,1)+ennv(1:nj,1))*1.0e6
   n0w_in(1:nj) = ennw(1:nj,1)*1.0e6 ! m^-3
   n0v_in(1:nj) = ennv(1:nj,1)*1.0e6 ! m^-3
   t0_in (1:nj) = tn(1:nj,1) ! keV
   f0w_in = fluxn(1)*sfarea  ! #/sec

   IF (iangrot .GT. 0) THEN
       omega_in(1:nj) = angrot(1:nj)
   ELSE
       omega_in(1:nj) = 0.0d0
   ENDIF
   vloop_in(1:nj) = 2.*pi*rmajor*etor(1:nj)

  ! ------------------------------------------------------------------
  ! write files

   command ="cp " // TRIM(eqfile) //" dnubeam_eqdsk.dat"
   ierr = ishell(command) 

   CALL wrt_dnubeam_nbi()
   CALL wrt_dnubeam_nml()
   CALL wrt_dnubeam_prof()

END

SUBROUTINE read_nubeam_output()

  USE dnubeam_mod 
  USE io,                          ONLY : lun_nubeam
  USE Plasma_properties,           ONLY : neut_beam
  IMPLICIT NONE
  INTEGER ierr,k

  NAMELIST/dnubeam_out/                     &
    nrho_out,rho_out,                       &
    pbe_out,pbi_out,pbth_out,               &
    curb_out,ucurb_out,curbdotb_out,        &
    bdep_out,bdenss_out,                    &
    tqbsum_out,tqbi_out,tqbe_out,           &
    udenspl_out,udenspp_out,                &
    pbe_tot,pbi_tot,pbth_tot,               &
    curb_tot,ucurb_tot,                     &
    bdep_tot,                               &
    bdenss_tot,                             &
    tqbsum_tot,tqbi_tot,tqbe_tot,           &
    udenspl_tot,udenspp_tot,                &
    bbntot_out,btneut_out
   
  OPEN(unit=lun_nubeam,file='dnubeam_out.dat',status='old',iostat=ierr)
  IF(ierr.NE.0) THEN
     WRITE(*,*) ' namelist file open failure.'
     CALL STOP ('subroutine read_nubeam_output ', 501)
     RETURN
  ENDIF
  READ(lun_nubeam,dnubeam_out)
  CLOSE(lun_nubeam)

  rho_save     (1:nrho_out,1) = rho_save     (1:nrho_out,2)
  pbe_save     (1:nrho_out,1) = pbe_save     (1:nrho_out,2)    
  pbi_save     (1:nrho_out,1) = pbi_save     (1:nrho_out,2)
  curb_save    (1:nrho_out,1) = curb_save    (1:nrho_out,2)
  curbdotb_save(1:nrho_out,1) = curbdotb_save(1:nrho_out,2) 
  ucurb_save   (1:nrho_out,1) = ucurb_save   (1:nrho_out,2)
  bdep_save    (1:nrho_out,1) = bdep_save    (1:nrho_out,2)
  bdenss_save  (1:nrho_out,1) = bdenss_save  (1:nrho_out,2)
  tqbsum_save  (1:nrho_out,1) = tqbsum_save  (1:nrho_out,2)
  tqbe_save    (1:nrho_out,1) = tqbe_save    (1:nrho_out,2)
  tqbi_save    (1:nrho_out,1) = tqbi_save    (1:nrho_out,2)
  udenspl_save (1:nrho_out,1) = udenspl_save (1:nrho_out,2)
  udenspp_save (1:nrho_out,1) = udenspp_save (1:nrho_out,2)

  rho_save     (1:nrho_out,2) = rho_out     (1:nrho_out)
  pbe_save     (1:nrho_out,2) = pbe_out     (1:nrho_out)
  pbi_save     (1:nrho_out,2) = pbi_out     (1:nrho_out)
  curb_save    (1:nrho_out,2) = curb_out    (1:nrho_out)
  curbdotb_save(1:nrho_out,2) = curbdotb_out(1:nrho_out)
  ucurb_save   (1:nrho_out,2) = ucurb_out   (1:nrho_out)
  bdep_save    (1:nrho_out,2) = bdep_out    (1:nrho_out)
  bdenss_save  (1:nrho_out,2) = bdenss_out  (1:nrho_out) ! no species index ??
  tqbsum_save  (1:nrho_out,2) = tqbsum_out  (1:nrho_out)
  tqbe_save    (1:nrho_out,2) = tqbe_out    (1:nrho_out)
  tqbi_save    (1:nrho_out,2) = tqbi_out    (1:nrho_out)
  udenspl_save (1:nrho_out,2) = udenspl_out (1:nrho_out)
  udenspp_save (1:nrho_out,2) = udenspp_out (1:nrho_out)

  pbe_tot_save    (1) = pbe_tot_save    (2)
  pbi_tot_save    (1) = pbi_tot_save    (2)
  pbth_tot_save   (1) = pbth_tot_save   (2)    
  curb_tot_save   (1) = curb_tot_save   (2)
  ucurb_tot_save  (1) = ucurb_tot_save  (2)
  bdep_tot_save   (1) = bdep_tot_save   (2)
  bdenss_tot_save (1) = bdenss_tot_save (2)
  tqbsum_tot_save (1) = tqbsum_tot_save (2)
  tqbi_tot_save   (1) = tqbi_tot_save   (2)
  tqbe_tot_save   (1) = tqbe_tot_save   (2)
  udenspl_tot_save(1) = udenspl_tot_save(2)
  udenspp_tot_save(1) = udenspp_tot_save(2)
  bbntot_save     (1) = bbntot_save     (2)
  btneut_save     (1) = btneut_save     (2)

  pbe_tot_save    (2) = pbe_tot       
  pbi_tot_save    (2) = pbi_tot       
  pbth_tot_save   (2) = pbth_tot          
  curb_tot_save   (2) = curb_tot   
  ucurb_tot_save  (2) = ucurb_tot  
  bdep_tot_save   (2) = bdep_tot   
  bdenss_tot_save (2) = bdenss_tot 
  tqbsum_tot_save (2) = tqbsum_tot 
  tqbi_tot_save   (2) = tqbi_tot   
  tqbe_tot_save   (2) = tqbe_tot   
  udenspl_tot_save(2) = udenspl_tot
  udenspp_tot_save(2) = udenspp_tot
  bbntot_save     (2) = bbntot_out
  btneut_save     (2) = btneut_out
 
  IF(ASSOCIATED(neut_beam%rhog_beam))DEALLOCATE(neut_beam%rhog_beam)
  ALLOCATE(neut_beam%rhog_beam(nrho_out))
  neut_beam%rhog_beam(1:nrho_out) = rho_out(1:nrho_out)
  neut_beam%nj_beam   = nrho_out  

  RETURN

END   SUBROUTINE read_nubeam_output

SUBROUTINE load_nubeam(time)

   USE dnubeam_mod

   USE param,    ONLY : kj,kprim,kimp
   USE numbrs,   ONLY : nj,nprim,nimp,nion 
   USE sourc,    ONLY : qbeame,qbeami,qbeame_intg,qbeami_intg,    &
                        qbth,qbth_intg,                           &
                        curb,curb_intg,wbeam_intg,                &
                        enbeam_intg,sbeam,sbeam_intg
   USE fusion,   ONLY : beam_thermal_dtntot, beam_thermal_ddntot, &
                        beam_thermal_ddptot, beam_thermal_tt2ntot,&
                        beam_beam_dtntot, beam_beam_ddntot,       & 
                        beam_beam_tt2ntot,beam_beam_ddptot,       & 
                        beam_beamddn,beam_beamdtn,beam_beamddp,   &
                        beam_beamtt2n,beam_thermaltth_df,         &
                        beam_thermalddp, beam_thermalddn,         &
                        beam_thermaltt2n,beam_thermaldth_tf 
   USE nub,      ONLY : nbeams,sbcx,ibion
   USE nub2,     ONLY : enbeam,wbeam
   USE tordlrot, ONLY : storqueb,storqueb_intg,                   &
                        sprbeame,sprbeame_intg,                   &
                        sprbeami,sprbeami_intg,                   &
                        storque, angrot
   USE extra,    ONLY : rgeom,btgeom
   USE machin,   ONLY : rmajor,btor

   IMPLICIT NONE

   REAL*8, INTENT(in)  :: time

   REAL*8, PARAMETER :: conv_power   = 6.2415097e+9  ! watts/m^3 to kev/(cm^3sec)
   REAL*8, PARAMETER :: conv_current = 1.0e-4        ! A/m^2 to A/(cm^2)
   REAL*8, PARAMETER :: conv_density = 1.0e-6        ! #/m^3 to #/cm^3
   REAL*8, PARAMETER :: conv_torque  = 10.           ! NT-M/M**3 to dyne-cm/cm**3
   REAL*8, PARAMETER :: conv_energy  = 0.62415064e10 ! 0.62415064e22 ! J/M**3 to kev/cm**3

   REAL*8 :: wbeam_pl(kj),wbeam_pp(kj)
   REAL*8 :: wbeam_pl_intg, wbeam_pp_intg

   WRITE(*,*) 'load_nubeam: ',time,t0_nubeam,t1_nubeam
   WRITE(*,*) 'rmajor,btor,abs(btor): ',rmajor,btor,ABS(btor)

   CALL interp_nubeam(enbeam  ,bdenss_save  ,time,t0_nubeam,t1_nubeam,nj,kj,nrhomax)
   CALL interp_nubeam(sbeam  ,bdep_save  ,time,t0_nubeam,t1_nubeam,nj,kj,nrhomax)
   CALL interp_nubeam(qbeame  ,pbe_save     ,time,t0_nubeam,t1_nubeam,nj,kj,nrhomax)
   CALL interp_nubeam(qbeami  ,pbi_save     ,time,t0_nubeam,t1_nubeam,nj,kj,nrhomax)
   CALL interp_nubeam(curb    ,curbdotb_save,time,t0_nubeam,t1_nubeam,nj,kj,nrhomax)
   CALL interp_nubeam(storqueb,tqbsum_save  ,time,t0_nubeam,t1_nubeam,nj,kj,nrhomax)
   CALL interp_nubeam(sprbeame,tqbe_save    ,time,t0_nubeam,t1_nubeam,nj,kj,nrhomax)
   CALL interp_nubeam(sprbeami,tqbi_save    ,time,t0_nubeam,t1_nubeam,nj,kj,nrhomax)
   CALL interp_nubeam(wbeam_pl,udenspl_save ,time,t0_nubeam,t1_nubeam,nj,kj,nrhomax)
   CALL interp_nubeam(wbeam_pp,udenspp_save ,time,t0_nubeam,t1_nubeam,nj,kj,nrhomax)
   wbeam(1:nj) = wbeam_pl(1:nj)+wbeam_pp(1:nj) 

   CALL interp_nubeam_0d(enbeam_intg  ,bdenss_tot_save ,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(sbeam_intg  ,bdep_tot_save ,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(qbeame_intg  ,pbe_tot_save    ,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(qbeami_intg  ,pbi_tot_save    ,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(qbth_intg    ,pbth_tot_save   ,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(curb_intg    ,curb_tot_save   ,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(storqueb_intg,tqbsum_tot_save ,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(sprbeame_intg,tqbe_tot_save   ,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(sprbeami_intg,tqbi_tot_save   ,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(wbeam_pl_intg,udenspl_tot_save,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(wbeam_pp_intg,udenspp_tot_save,time,t0_nubeam,t1_nubeam)
   wbeam_intg =  wbeam_pl_intg + wbeam_pp_intg

   CALL interp_nubeam_0d(beam_thermal_ddntot,btneut_save,time,t0_nubeam,t1_nubeam)
   CALL interp_nubeam_0d(beam_beam_ddntot,bbntot_save,time,t0_nubeam,t1_nubeam)

   enbeam(1:nj)   = enbeam(1:nj)   * conv_density
   sbeam(1:nj)    = sbeam(1:nj)    * conv_density
   qbeame(1:nj)   = qbeame(1:nj)   * conv_power
   qbeami(1:nj)   = qbeami(1:nj)   * conv_power 
   curb(1:nj)     = curb(1:nj)     * conv_current/(ABS(btor)*1.e-4)
   storqueb(1:nj) = storqueb(1:nj) * conv_torque
   sprbeame(1:nj) = sprbeame(1:nj) * conv_torque
   sprbeami(1:nj) = sprbeami(1:nj) * conv_torque
   wbeam(1:nj)    = wbeam(1:nj)    * conv_energy

   enbeam_intg    = enbeam_intg    * conv_density
   sbeam_intg     = sbeam_intg     * conv_density
   qbeame_intg    = qbeame_intg    * conv_power
   qbeami_intg    = qbeami_intg    * conv_power
   qbth_intg      = qbth_intg      * conv_power
   curb_intg      = curb_intg      * conv_current
   storqueb_intg  = storqueb_intg  * conv_torque
   sprbeame_intg  = sprbeame_intg  * conv_torque
   sprbeami_intg  = sprbeami_intg  * conv_torque
   wbeam_intg     = wbeam_intg     * conv_energy

   RETURN

END SUBROUTINE load_nubeam

SUBROUTINE nubeam_beam_startup()

  USE dnubeam_mod

  USE solcon,   ONLY : time,time0
  USE param,    ONLY : kk,kj
  USE numbrs,   ONLY : nj,nk,nion
  USE sourc,    ONLY : curb,curbe,curbet,curbi,sbeam
  USE nub2,     ONLY : enbeam,enb,enbmin,enbsav,enbs,wbeam,wb  
  USE soln,     ONLY : u,usave,te,ti,ene,enesav,en,rbp
  USE tordlrot, ONLY : iangrot,angrot
  USE transp,   ONLY : nubeam_restart,nubeam_back_delt,nubeam_init, &
                       beam_data
  USE io_gcnmp, ONLY : nlog,ncrt

  IMPLICIT NONE
  REAL*8 :: time_save
  INTEGER :: j

  IF (nubeam_restart .NE. -1) RETURN 

  WRITE(*, fmt = '("\n\n==========================================")' )
  WRITE(*, fmt = '("+ calling nubeam for initial guess of beam density")' )
  IF(beam_data%nbeam_species .GT. 1)THEN
     WRITE(ncrt,10)beam_data%nbeam_species
     WRITE(nlog,10)beam_data%nbeam_species
10   FORMAT("ERROR: nubeam version after  201107 returns ",/, &
     "     only one species density input has",i5,/,          &
     "     species set - code must stop")
     CALL STOP('nubeam only 1 beam species allowed',1)
  ENDIF
  time_save = time0
  time = time0-nubeam_back_delt
  
  DO WHILE (time .LT. time_save)

    CALL dnubeam_driver()
    time = t1_nubeam+(t1_nubeam-t0_nubeam)*0.5

  ENDDO 

  time = time_save
  CALL load_nubeam(time)
  curbe(:) = 0.0 
  curbet(:) = 0.0 
  curbi(:) = curb(:)
  wb(:,1,1) = wbeam(:)/0.62415064e16
  enb(:,1,1) = enbeam(:)
  enbs(1:nj) = enbeam(1:nj)
  enbsav(:,1,1) = enbeam(:)
  DO j=1,SIZE(enb,dim=1)
    enb(j,1,1) = MAX(enb(j,1,1),enbmin)
  ENDDO

  CALL redate (u,en,te,ti,rbp,nk,nj,kj,kk,iangrot,angrot)
  CALL savesol (u, usave, nk, nj, kk) ! sets usave =u 
  CALL update (u, en, te, ti, rbp, nk, nj, kj, kk,iangrot,angrot)
  CALL copya (ene, enesav, nj)

  nubeam_init = .TRUE.

  WRITE(*, fmt = '("\n\n+ leaving nubeam_beam_startup")' )
  WRITE(*, fmt = '("==========================================\n\n")' )

  RETURN

END 

SUBROUTINE init_nubeam_data()

  USE dnubeam_mod

  IMPLICIT NONE

  rho_save(:,:)       = 0.0 
  pbe_save(:,:)       = 0.0
  pbi_save(:,:)       = 0.0
  curb_save(:,:)      = 0.0
  curbdotb_save(:,:)  = 0.0
  ucurb_save(:,:)     = 0.0
  bdep_save(:,:)      = 0.0
  bdenss_save(:,:)    = 0.0
  pbth_save(:,:)      = 0.0
  tqbsum_save(:,:)    = 0.0
  tqbi_save(:,:)      = 0.0
  tqbe_save(:,:)      = 0.0
  udenspl_save(:,:)   = 0.0
  udenspp_save(:,:)   = 0.0

  pbe_tot_save(:)     = 0.0
  pbi_tot_save(:)     = 0.0
  pbth_tot_save(:)    = 0.0
  curb_tot_save(:)    = 0.0
  ucurb_tot_save(:)   = 0.0
  bdep_tot_save(:)    = 0.0
  bdenss_tot_save(:)  = 0.0
  tqbsum_tot_save(:)  = 0.0
  tqbi_tot_save(:)    = 0.0
  tqbe_tot_save(:)    = 0.0
  udenspl_tot_save(:) = 0.0
  udenspp_tot_save(:) = 0.0
  bbntot_save(:)      = 0.0
  btneut_save(:)      = 0.0
  RETURN

END SUBROUTINE init_nubeam_data

SUBROUTINE interp_nubeam(fout,f,t,t0,t1,n,kj,kk)

    IMPLICIT NONE
    REAL*8, DIMENSION(kj), INTENT(out) :: fout
    REAL*8, DIMENSION(kk,2), INTENT(in) :: f
    REAL*8, INTENT(in) :: t,t0,t1
    INTEGER, INTENT(in) :: n,kj,kk
    INTEGER j

    DO j=1,n 
        fout(j) = f(j,1)+(f(j,2)-f(j,1))*(t-t0)/(t1-t0)
    ENDDO
    RETURN
END SUBROUTINE interp_nubeam

SUBROUTINE interp_nubeam_0d(fout,f,t,t0,t1)

    IMPLICIT NONE
    REAL*8 :: fout
    REAL*8, DIMENSION(2), INTENT(in) :: f
    REAL*8, INTENT(in) :: t,t0,t1
    fout = f(1)+(f(2)-f(1))*(t-t0)/(t1-t0)
    RETURN
END SUBROUTINE interp_nubeam_0d

SUBROUTINE read_nubeam_namelist(lun,nbname,ierr)
!---------------------------------------------------------------------
! -- orig created by JM PARK ?
! -- modifed by HSJ to use namelist module nbnamelist_201112.F90
! -- called in ufiles_12.f90 if nubeam version >= 201107
!-----------------------------------------------------HSJ-2-14-2012---
  USE dnubeam_mod
  USE beam_structure, ONLY : neutral_beam

  USE transp, ONLY : beam_data
  USE transp, ONLY : d_fast_ion
  USE transp, ONLY : nubeam_dt_old   => nubeam_dt
  USE transp, ONLY : nubeam0_dt_old  => nubeam0_dt
 
!  use nbnamelist, only : tbona, tboffa
  USE nbnamelist, ONLY : nbeam_old => nbeam
  IMPLICIT NONE

  TYPE(neutral_beam),INTENT(out)  :: nbname

  REAL*8, DIMENSION(nbmax) :: tbona, tboffa
  REAL*8  :: nubeam_dt
  REAL*8  :: ebdmax
  LOGICAL :: nlbdat

  ! this include file  contains all !***** info below HSJ
  ! note that P_Nfreya_nubeam_namelist also includes this file HSJ
  include './shared_modules/nbnamelist_post_201112.inc'

  !*****namelist/nbdrive_naml/                   &
  !*****    nbeam,                               &
  !*****    nlco,                                &
  !*****    abeama,xzbeama,                      &
  !*****    nbshapa,                             &
  !*****    rtcena,xlbtna,xybsca,                &
  !*****    bmwidra,bmwidza,                     &
  !*****    divra,divza,foclza,foclra,           &
  !*****    nbapsha,                             &
  !*****    rapedga,xzpedga,xlbapa,xybapa,       &
  !*****    xrapoffa,xzapoffa,                   &
  !*****    nbapsh2,                             &
  !*****    rapedg2,xzpedg2,xlbapa2,             &
  !*****    xrapoff2,xzapoff2,                   &
  !*****    xbzeta

  !*****namelist/nbdrive_naml/                 &
  !*****    pinja,einja,ffulla,fhalfa

  !*****namelist/nbdrive_naml/                 &
  !*****    tbona,tboffa

  !*****namelist/nbdrive_naml/                 &
  !*****    nubeam_dt

  !*****namelist/nbdrive_naml/                 &
  !*****    dn0out

  !*****namelist/nbdrive_naml/                 &
  !*****    nkdifb

  !*****namelist/nbdrive_naml/                 &
  !*****    adiff_a,adiff_0,adiff_xpin,adiff_xpout,adiff_time,adiff_ntime

  !*****namelist/nbdrive_naml/                 &
  !*****    ebdmax,nlbdat

  !*****namelist/nbdrive_naml/                                          &
  !*****     !mission, workpath,                                           &  
  !*****     orbrzv_toric_prftot0, orbrzv_toric_frnm,                     &
  !*****     xdatsf3, xdatsfa, xdatsft, xdatsfp,                          &
  !*****     plfhe3, plfhe4, plfst, plfsp,                                &
  !*****     cxsplt, dxbsmoo, dxbsm_nc, xpand_nptcl,                      &
  !*****     frac_depmax, frac_dep_lim, frac_depmin, frac_orbrr,          &
  !*****     aimp, xzimp, wghta, goocon,                                  &
  !*****     cxpcon, fppcon, fdtnxy, dtn,                                 &
  !*****     dt_acc, orbrzv_zzerr_con, xdepmod, xcfanbi,                  &
  !*****     xdfanbi, xefanbi, xcfafus, xdfafus,                          &
  !*****     xefafus, gflr_min, gflr_rl, gflr_ll,                         &
  !*****     gflr_xv, gflr_op, xswfrac_allfast, xswfrac_beam,             &
  !*****     xswfrac_fusn, fporcelli, taurip, asrd,                       &
  !*****    bsrd, erngfi, fbemin, fbemax,                                &
  !*****     fvpvmn, fvpvmx, fbltim, fshper,                              &
  !*****     fshwid, tfshon, tfshof, xfishmin,                            &
  !*****     xfishmax, fbtrap_depth, xboxhw, yboxhw,                      &
  !*****     xlbox1, xlbox2, rtube, ytube,                                &
  !*****     xzetatube, phitube, thetatube, xl1tube,                      &
  !*****     xl2tube, rhotube, edbfac, nonlin,                            &
  !*****     nseed, nlbout, lunnbx, lunres,                               &
  !*****     mrstrt, nclass, only_io, ref_namelist,                       &
  !*****     quasi_check, nltest_output, blk_mpi_bcast_int, blk_mpi_bcast_r8,  &
  !*****     nzones, nzone_fb, nlsym2b, nth0,                             &
  !*****     xplasma_in_memory, nznbma, nznbme, ngyro,                    &
  !*****     nlusf3, nlusfa, nlusft, nlusfp,                              &
  !*****     nlfatom, nlbfpp, nptcls, nptclf,                             &
  !*****     nptclh, nptcl_max, nltrk_dep0, nldep0_gather,                &
  !*****     ndep0, ndep0_max, ndep_set_beam, ndep_set_ien,               &
  !*****     nbbcal, ncx0, nper_cx, lev_nbidep,                           &
  !*****     levmod_halo, nmimp, ngoocon_vpvbin, ngoocon_ebin,            &
  !*****     ngoocon_rbin, ndtorb, nchdvp, orbrzv_option,                 &
  !*****     orbrzv_rf_option, nlminsv, nlebei, nsigexc,                  &
  !*****     nlbbcx, nbbcx_avg, nbbcx_bb, nbfallgr,                       &
  !*****     nlbeamcx, nlhvion, nlbcde, nlbcoh,                           &
  !*****     nlbcpa, nlorbo, nlbflr, nlfbmflr,                            &
  !*****     nlbgflr, sawflag, nlsawb, nlsawf,                            &
  !*****     nmix_kdsaw, ngradd_opt, nrip, nerngfi,                       &
  !*****     nsdbgb, nlfbon, nfbon_vpvopt, nfbon_species,                 &
  !*****     nbbox, nbsbox, nbebox, nxbox,                                &
  !*****     nybox, nlbox, ndepbox, nbtube,                               &
  !*****     nbstube, nbetube, lmidtube, nsegtube,                        &
  !*****     ndeptube, nlcprb, nlpsirz, nmcurb,                           &
  !*****     nlfdep

  INTEGER, INTENT (in) :: lun
  INTEGER, INTENT (out) :: ierr

  ierr = 0

  adiff_ntime = 1

  CALL set_nubeam_init_default()

  READ(lun,nml = nbdrive_naml)

  CALL logical2int()
 
  nbeam_old = nbeam
  nbname%xrapoffa (:) = xrapoffa(1:nbeam)
  nbname%xzapoffa (:) = xzapoffa(1:nbeam)
  nbname%xrapoff2 (:) = xrapoff2(1:nbeam)
  nbname%xzapoff2 (:) = xzapoff2(1:nbeam)
  nbname%nbapsh2  (:) = nbapsh2(1:nbeam)

  nbname%abeama   (:) = abeama  (1:nbeam)
  nbname%xzbeama  (:) = xzbeama (1:nbeam)
  nbname%bmwidra  (:) = bmwidra (1:nbeam)
  nbname%bmwidza  (:) = bmwidza (1:nbeam) 
  nbname%rtcena   (:) = rtcena  (1:nbeam)
  nbname%xlbtna   (:) = xlbtna  (1:nbeam)
  nbname%xybsca   (:) = xybsca  (1:nbeam)
  nbname%divza    (:) = divza   (1:nbeam)
  nbname%divra    (:) = divra   (1:nbeam) 
  nbname%foclza   (:) = foclza  (1:nbeam)
  nbname%foclra   (:) = foclra  (1:nbeam)
  nbname%rapedga  (:) = rapedga (1:nbeam)
  nbname%xzpedga  (:) = xzpedga (1:nbeam)
  nbname%xlbapa2  (:) = xlbapa2 (1:nbeam)
  nbname%xlbapa   (:) = xlbapa  (1:nbeam)
  nbname%xybapa   (:) = xybapa  (1:nbeam)
  nbname%rapedg2  (:) = rapedg2 (1:nbeam)
  nbname%xzpedg2  (:) = xzpedg2 (1:nbeam)
  nbname%xbzeta   (:) = xbzeta  (1:nbeam) 
  nbname%ntrace   (:) = ntrace  (1:nbeam)
  nbname%nbshapa  (:) = nbshapa (1:nbeam)
  nbname%nbapsha  (:) = nbapsha (1:nbeam)
  nbname%nlco     (:) = nlco    (1:nbeam)
  nbname%nbeam        = nbeam

  nbname%pinja    (:) = pinja   (1:nbeam)
  nbname%einja    (:) = einja   (1:nbeam)
  nbname%ffulla   (:) = ffulla  (1:nbeam)
  nbname%fhalfa   (:) = fhalfa  (1:nbeam) 
  nbname%tbona    (:) = tbona   (1:nbeam)
  nbname%tboffa   (:) = tboffa  (1:nbeam)

  nbname%nbbcal   = nbbcal
  nbname%nseed    = nseed 
  nbname%nzone_nb = nzones

  nubeam_dt_old   = nubeam_dt
  nubeam0_dt_old  = nubeam_dt

  d_fast_ion%adiff_a     = adiff_a(1)
  d_fast_ion%adiff_0     = adiff_0(1)
  d_fast_ion%adiff_xpin  = adiff_xpin(1)
  d_fast_ion%adiff_xpout = adiff_xpout(1)
  d_fast_ion%nkdifb      = nkdifb

  RETURN

END

SUBROUTINE set_nubeam_init_default()
 
     USE dnubeam_mod
     IMPLICIT NONE
 
     !mission              = "standard"
     !workpath             = " "
     orbrzv_toric_prftot0 = 0.0d0
     orbrzv_toric_frnm    = 1.0d0
     xdatsf3              = 0.0d0
     xdatsfa              = 0.0d0
     xdatsft              = 0.0d0
     xdatsfp              = 0.0d0
     plfhe3               = 100.0d0
     plfhe4               = 10000.0d0
     plfst                = 100.0d0
     plfsp                = 100.0d0
     cxsplt               = 2.0d0
     dxbsmoo              = 0.05d0
     dxbsm_nc             = 0.0d0
     xpand_nptcl          = 3.0d0
     frac_depmax          = 0.25d0
     frac_dep_lim         = 1.5d0
     frac_depmin          = 0.01d0
     frac_orbrr           = 0.10d0
     aimp                 = 16.d0
     xzimp                = 8.d0
     wghta                = 1.0d0
     goocon               = 10.0d0
     cxpcon               = 20.0d0
     fppcon               = 8.0d0
     fdtnxy               = 0.5d0
     dtn                  = 0.0d0
     dt_acc               = 1.0d-4
     orbrzv_zzerr_con     = 1.0d-4
     xdepmod              = 1.0d0
     xcfanbi              = 1.0d0
     xdfanbi              = 1.0d0
     xefanbi              = 1.0d0
     xcfafus              = 1.0d0
     xdfafus              = 1.0d0
     xefafus              = 1.0d0
     gflr_min             = 0.0d0
     gflr_rl              = 1.0d6
     gflr_ll              = 30.0d0
     gflr_xv              = (/1.5d0, 0.4d0/)
     gflr_op              = (/0.d0,1.d-4,1.d-5,0.d0,3.d0,10.d0,1.d0,0.2d0,0.05d0,1.d0,2.d0,0.d0,0.d0,0.d0,0.d0/)
     xswfrac_allfast      = 1.0d0
     xswfrac_beam         = 1.0d0
     xswfrac_fusn         = 1.0d0
     fporcelli            = 1.0d0
     taurip               = 0.0
     asrd                 = 0.0
     bsrd                 = 0.0
     erngfi               = 0.0
     fbemin               = 0.0
     fbemax               = 0.0
     fvpvmn               = 0.0d0
     fvpvmx               = 0.0d0
     fbltim               = 1.0d-7
     fshper               = 5.0d-3
     fshwid               = 1.0d-3
     tfshon               = -1000.0
     tfshof               = +1.0d34
     xfishmin             = 0.0d0
     xfishmax             = 1.0d0
     fbtrap_depth         = 0.0d0
     xboxhw               = 5.0d0
     yboxhw               = 5.0d0
     xlbox1               = 100.0d0
     xlbox2               = 500.0d0
     rtube                = 0.0d0
     ytube                = 0.0d0
     xzetatube            = 0.d0
     phitube              = 0.d0
     thetatube            = 0.d0
     xl1tube              = 0.d0
     xl2tube              = 0.d0
     rhotube              = 1.d0
     edbfac               = 10.0d0
  
     nonlin               = 6
     nseed                = 1
     nlbout               = 1 !.true.
     lunnbx               = 0
     lunres               = 0
     mrstrt               = 1
     nclass               = 0
     only_io              = 0 !.false.
     ref_namelist         = 0 !.false.
     quasi_check          = 1 !.true.
     nltest_output        = 1 !.true.
     blk_mpi_bcast_int    = 2000000
     blk_mpi_bcast_r8     = 30000000
     nzones               = 0
     nzone_fb             = 0
     nlsym2b              = 0 !.false.
     nth0                 = 2
     xplasma_in_memory    = 0
     nznbma               = 50
     nznbme               = 100
     ngyro                = 4
     nlusf3               = 0 !.false.
     nlusfa               = 0 !.false.
     nlusft               = 0 !.false.
     nlusfp               = 0 !.false.
     nlfatom              = 1 !.true.
     nlbfpp               = 0 !.false.
     nptcls               = 1000
     nptclf               = 1000
     nptclh               = 500
     nptcl_max            = 0
     nltrk_dep0           = 0 !.false.
     nldep0_gather        = 0 !.false.
     ndep0                = 500
     ndep0_max            = 100000
     ndep_set_beam        = 0
     ndep_set_ien         = 0
     nbbcal               = 1000000
     ncx0                 = 100
     nper_cx              = 4
     lev_nbidep           = 1
     levmod_halo          = 0
     nmimp                = 1
     ngoocon_vpvbin       = 0
     ngoocon_ebin         = 1
     ngoocon_rbin         = 1
     ndtorb               = 4
     nchdvp               = 3
     orbrzv_option        = 0
     orbrzv_rf_option     = 0
     nlminsv              = 1 !.true.
     nlebei               = 1 !.true.
     nsigexc              = 0
     nlbbcx               = 1 !.true.
     nbbcx_avg            = 0
     nbbcx_bb             = 0
     nbfallgr             = 0
     nlbeamcx             = 1 !.true.
     nlhvion              = 1 !.true.
     nlbcde               = 1 !.true.
     nlbcoh               = 1 !.true.
     nlbcpa               = 1 !.true.
     nlorbo               = 0 !.false.
     nlbflr               = 1 !.true.
     nlfbmflr             = 1 !.true.
     nlbgflr              = 0 !.false.
     sawflag              = 0 !.false.
     nlsawb               = 1 !.true.
     nlsawf               = 1 !.true.
     nmix_kdsaw           = 1
     ngradd_opt           = 0
     nrip                 = 0
     nerngfi              = 0
     nsdbgb               = 1
     nlfbon               = 0 !.false.
     nfbon_vpvopt         = 0
     nfbon_species        = 3
     nbbox                = 0
     nbsbox               = 0
     nbebox               = 1
     nxbox                = 1
     nybox                = 1
     nlbox                = 100
     ndepbox              = 500
     nbtube               = 0
     nbstube              = 0
     nbetube              = 1
     lmidtube             = 0
     nsegtube             = 100
     ndeptube             = 500
     nlcprb               = 1 !.true.
     nlpsirz              = 0 !.false.
     nmcurb               = 1
     nlfdep               = 0 !.false.
 
     RETURN
END

SUBROUTINE logical2int()

     USE dnubeam_mod
     IMPLICIT NONE

     nlbout               =  ABS(nlbout)
     only_io              =  ABS(only_io)
     ref_namelist         =  ABS(ref_namelist)
     quasi_check          =  ABS(quasi_check)
     nltest_output        =  ABS(nltest_output)
     nlsym2b              =  ABS(nlsym2b)
     nlusf3               =  ABS(nlusf3)
     nlusfa               =  ABS(nlusfa)
     nlusft               =  ABS(nlusft)
     nlusfp               =  ABS(nlusfp)
     nlfatom              =  ABS(nlfatom)
     nlbfpp               =  ABS(nlbfpp)
     nltrk_dep0           =  ABS(nltrk_dep0)
     nldep0_gather        =  ABS(nldep0_gather)
     nlminsv              =  ABS(nlminsv)
     nlebei               =  ABS(nlebei)
     nlbbcx               =  ABS(nlbbcx)
     nlbeamcx             =  ABS(nlbeamcx)
     nlhvion              =  ABS(nlhvion)
     nlbcde               =  ABS(nlbcde)
     nlbcoh               =  ABS(nlbcoh)
     nlbcpa               =  ABS(nlbcpa)
     nlorbo               =  ABS(nlorbo)
     nlbflr               =  ABS(nlbflr)
     nlfbmflr             =  ABS(nlfbmflr)
     nlbgflr              =  ABS(nlbgflr)
     sawflag              =  ABS(sawflag)
     nlsawb               =  ABS(nlsawb)
     nlsawf               =  ABS(nlsawf)
     nlcprb               =  ABS(nlcprb)
     nlpsirz              =  ABS(nlpsirz)
     nlfdep               =  ABS(nlfdep)

END 

SUBROUTINE wrt_dnubeam_nml()

  USE dnubeam_mod
  USE io,               ONLY : lun_nubeam
  IMPLICIT NONE

  NAMELIST/nbi_init/                                              &
     !mission, workpath,                                          &  
     orbrzv_toric_prftot0, orbrzv_toric_frnm,                     &
     xdatsf3, xdatsfa, xdatsft, xdatsfp,                          &
     plfhe3, plfhe4, plfst, plfsp,                                &
     cxsplt, dxbsmoo, dxbsm_nc, xpand_nptcl,                      &
     frac_depmax, frac_dep_lim, frac_depmin, frac_orbrr,          &
     aimp, xzimp, wghta, goocon,                                  &
     cxpcon, fppcon, fdtnxy, dtn,                                 &
     dt_acc, orbrzv_zzerr_con, xdepmod, xcfanbi,                  &
     xdfanbi, xefanbi, xcfafus, xdfafus,                          &
     xefafus, gflr_min, gflr_rl, gflr_ll,                         &
     gflr_xv, gflr_op, xswfrac_allfast, xswfrac_beam,             &
     xswfrac_fusn, fporcelli, taurip, asrd,                       &
     bsrd, erngfi, fbemin, fbemax,                                &
     fvpvmn, fvpvmx, fbltim, fshper,                              &
     fshwid, tfshon, tfshof, xfishmin,                            &
     xfishmax, fbtrap_depth, xboxhw, yboxhw,                      &
     xlbox1, xlbox2, rtube, ytube,                                &
     xzetatube, phitube, thetatube, xl1tube,                      &
     xl2tube, rhotube, edbfac, nonlin,                            &
     nseed, nlbout, lunnbx, lunres,                               &
     mrstrt, nclass, only_io, ref_namelist,                       &
     quasi_check, nltest_output, blk_mpi_bcast_int, blk_mpi_bcast_r8,  &
     nzones, nzone_fb, nlsym2b, nth0,                             &
     xplasma_in_memory, nznbma, nznbme, ngyro,                    &
     nlusf3, nlusfa, nlusft, nlusfp,                              &
     nlfatom, nlbfpp, nptcls, nptclf,                             &
     nptclh, nptcl_max, nltrk_dep0, nldep0_gather,                &
     ndep0, ndep0_max, ndep_set_beam, ndep_set_ien,               &
     nbbcal, ncx0, nper_cx, lev_nbidep,                           &
     levmod_halo, nmimp, ngoocon_vpvbin, ngoocon_ebin,            &
     ngoocon_rbin, ndtorb, nchdvp, orbrzv_option,                 &
     orbrzv_rf_option, nlminsv, nlebei, nsigexc,                  &
     nlbbcx, nbbcx_avg, nbbcx_bb, nbfallgr,                       &
     nlbeamcx, nlhvion, nlbcde, nlbcoh,                           &
     nlbcpa, nlorbo, nlbflr, nlfbmflr,                            &
     nlbgflr, sawflag, nlsawb, nlsawf,                            &
     nmix_kdsaw, ngradd_opt, nrip, nerngfi,                       &
     nsdbgb, nlfbon, nfbon_vpvopt, nfbon_species,                 &
     nbbox, nbsbox, nbebox, nxbox,                                &
     nybox, nlbox, ndepbox, nbtube,                               &
     nbstube, nbetube, lmidtube, nsegtube,                        &
     ndeptube, nlcprb, nlpsirz, nmcurb,                           &
     nlfdep

  NAMELIST/nbi_update/                                            &
    nltest_output


  OPEN (unit=lun_nubeam,file='dnubeam_nml.dat',status='unknown')
  WRITE(lun_nubeam,nml=nbi_init)
  WRITE(lun_nubeam,nml=nbi_update)
  CLOSE(lun_nubeam)

END

SUBROUTINE wrt_dnubeam_nbi()

  USE dnubeam_mod 
  USE io,                      ONLY : lun_nubeam
  IMPLICIT NONE

  NAMELIST/dnubeam_nbi/                     &
    nbeam,                                  &
    nlco,                                   &
    pinja,einja,ffulla,fhalfa,              &
    abeama,xzbeama,                         &
    nbshapa,nbapsha,                        &
    rtcena,xlbapa,xlbtna,xybapa,xybsca,     &
    xrapoffa,xzapoffa,                      &
    bmwidra,bmwidza,                        &
    rapedga,xzpedga,                        &
    divra,divza,foclza,foclra,              &
    nbapsh2,                                &
    rapedg2,xzpedg2,xlbapa2,                &
    xrapoff2, xzapoff2
   OPEN (unit=lun_nubeam,file='dnubeam_nbi.dat',status='unknown')
   WRITE(lun_nubeam,nml=dnubeam_nbi)
   CLOSE(lun_nubeam)

END

SUBROUTINE wrt_dnubeam_prof()

  USE dnubeam_mod 
  USe io,                        ONLY : lun_nubeam
  IMPLICIT NONE

  NAMELIST/dnubeam_prof/                    &
    t0_nubeam,t1_nubeam,                    &
    nrho_in,                                &
    rho_in,                                 &
    te_in,ti_in,ne_in,zeff_in,omega_in,     &
    vloop_in,                               &
    n0w_in,n0v_in,t0_in,f0w_in,             &
    dn0out

  NAMELIST/dnubeam_difb/                    &
    difb_0,difb_a,difb_in,difb_out,nkdifb

  OPEN (unit=lun_nubeam,file='dnubeam_prof.dat',status='unknown')
  WRITE(lun_nubeam,nml=dnubeam_prof)
  WRITE(lun_nubeam,nml=dnubeam_difb)
 !write(lun_nubeam,nml=dnubeam_etc)
  CLOSE(lun_nubeam)

END

