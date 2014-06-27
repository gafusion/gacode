MODULE kinetic_efit    
!          
  USE nrtype,                         ONLY     : DP,I4B
  USE param,                          ONLY     :  kj
 
    USE param   ,ONLY     :  kj
    USE solcon  ,ONLY     :  time,external_beam_cur   
    USE soln    ,ONLY     :  te,ti,en,ene,curden,rbp,curpar_soln
    USE yoka    ,ONLY     :  ishot,pbe,pbi  
    USE ions    ,ONLY     :  atw
    USE numbrs  ,ONLY     :  nj,nion
    USE nub     ,ONLY     :  ibion
    USE nub2    ,ONLY     :  enbeam,wbeam
    USE fusion  ,ONLY     :  walp,enalp
    USE machin  ,ONLY     :  rmajor,btor
    USE mhdcom  ,ONLY     :  mhdmethd
    USE geom    ,ONLY     :  rcapi,rcap,r2cap,r2capi,bsqncap,bsq_avg_cap, &
                             fcap,gcap,hcap,b_avg_cap
    USE bd_condtn,ONLY    :  totcur
    USE mesh     ,ONLY    :  r
    USE CONSTNTS ,ONLY    :  psimks ,twopi,u0
    USE psig     ,ONLY    :  cxareanpsi
    USE extra    ,ONLY    :  q
    USE rhog     ,ONLY    :  rmajorvec,bprmaj,btotrmaj,psir,pprim 
    USE sourc    ,ONLY    :  curboot,curbe,curbi,curbet,curdri,          &
                             curb_external,currf,curohm,totboot

  IMPLICIT NONE
  INTEGER,PARAMETER :: strlen = 64
  LOGICAL wrt_kinetic_efit,first_call        
  CHARACTER(strlen) runid_EFIT  
  CHARACTER(strlen) kin_efit_namelist
  CHARACTER(strlen) kine_message

  REAL(DP),DIMENSION(:),ALLOCATABLE  ::                                      &
                             palpha,r_loc,psir_loc,                          &
                             rcap_loc,rcapi_loc,curden_loc,r2capi_loc,       &
                             q_loc,dpsidrho_loc,I_loc,curb_loc,curboot_loc,  &
                             currf_loc, curohm_loc,rmajor_loc,               &
                             bprmaj_loc,btotrmaj_loc,pfscur,pfscur_avg,      &
                             jpar_rmaj,curpar_loc,dpdpsi_loc,dpdrho_loc,     &
                             pprim_loc,r2cap_loc

  ! namelist arrays not allocatable:
  REAL (DP) vzeroj(kj),sizeroj(kj),vzerobj(kj),rpress(kj),         &
              pressr(kj),sibeam(kj),pbeam(kj),dnbeam(kj),          &
              dmass(kj),sigpre(kj)
  REAL(DP) R0,Bt0
  DATA first_call /.TRUE./
  CONTAINS

  SUBROUTINE  wrt_kin_efit_naml(cgs,src_ident)
! -------------------------------------------------------------------------
!    subroutine writes namelist for EFIT  code
!    in file 12_Kshot.time

!    src_ident :    integer, identifies source of data to be written

!    call with cgs = 1 if rmajor,rcapi are  in cgs (cm) (as it is when this
!    routine is called from tport)
!    otherwise, cgs  =0 means input is in MKS
! -------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER i,j,io12,itime,npress,nbeam,nmass,kzeroj,cgs,src_ident,iend,istart
    REAL (DP)  pnbeam,convert,convert2,massden,masse,rmaj,rcapil, &
               jave,cxa,rmloc,btloc

    CHARACTER(len = 12) ishota,intfl,itimea
    CHARACTER(len=1),ALLOCATABLE,DIMENSION(:) :: temp_string
    CHARACTER(len=1) dummy
    NAMELIST /IN1/       ishot,itime,pnbeam,npress,nbeam,      &
                         rpress,pressr,sibeam,pbeam,dnbeam,    &
                         dmass,sigpre

    NAMELIST /INWANT/    kzeroj,vzeroj,sizeroj,vzerobj

    CALL kine_efit_alloc(kj)
!    IF(.NOT. first_call)src_ident = 0
    IF(src_ident == 1)THEN
       iend = LEN_TRIM(kine_message) ! iend is last non blank character
       istart = iend
       DO j=iend-1,1,-1
          IF(kine_message(j:j) .ne. ' ')istart = j
       ENDDO
          
       IF(ALLOCATED(temp_string))DEALLOCATE(temp_string)
       ALLOCATE(temp_string(iend-istart+1))
       i = 0
       DO j=istart,iend
         i = i +1
         temp_string(i) = kine_message(j:j)
       ENDDO
    ENDIF

    convert  = 1.6021765e-10       !kev/cm**3 to J/M**3
    convert2 = 1.6605387e-27       ! amu to kg
    masse   = 9.10938970e-31       ! electron mass   (kg)
    rmaj = rmajor
    cxa = cxareanpsi(1)
    curden_loc(:)  =  curden(:)     !curden is in amps/cm**2 if cgs =1, amps/m**2 if cgs=0
    IF(cgs  .EQ. 1 )THEN
       rmaj= rmajor*0.01               ! rmaj in meters
       cxa = cxareanpsi(1)*1.e-4       ! total cross section area m**2
       curden_loc(:) = curden(:)*1.e4   ! total current density in A/m**2
    ENDIF
    IF(.NOT. wrt_kinetic_efit) RETURN
    IF (mhdmethd .EQ. 'tdem') THEN
       PRINT *,'Error, kinetic efit output can not '
       PRINT *,'be written at this time because the'
       PRINT *,'tdem mode does not have calcualtions of <1/R>'
       CALL STOP('wrt_kin_efit_naml',1)
    ENDIF

    runid_EFIT ='EFIT input file created by ONETWO'
    itime  = 1000*time
    WRITE(intfl,FMT='(i7)')itime
    READ(intfl,FMT ='(a)')itimea
    itimea = ADJUSTL(itimea)
    WRITE(intfl,FMT='(i7)')ishot
    READ(intfl,FMT ='(a)')ishota
    ishota = ADJUSTL(ishota)
    IF(src_ident == 0)THEN
       kin_efit_namelist = '12_K'//ishota(1:LEN_TRIM(ishota))//'.' &
                         //itimea(1:LEN_TRIM(itimea))
    ELSE
       kin_efit_namelist = '12_K'//ishota(1:LEN_TRIM(ishota))//'.' &
                         //itimea(1:LEN_TRIM(itimea))//'_'
       DO i=1,iend
          dummy = temp_string(i)
          kin_efit_namelist = kin_efit_namelist(1:LEN_TRIM(kin_efit_namelist))//dummy
       ENDDO
    ENDIF
    pnbeam = 1.e6*(pbe + pbi) !beam power,watts, absorbed by electrons and ions
    npress =nj ; nbeam = nj ; nmass = nj ; kzeroj = nj
    DO j=1,nj
       rpress(j) = (-psir(j) + psir(1))/(psir(nj) -psir(1)) !negative [0, -1]
       pressr(j) =0.0
       massden =0.0
       DO i =1, nion
          pressr(j) = pressr(j) + ti(j)*en(j,i)
          massden   = massden   + atw(i)*en(j,i)
       ENDDO
       massden      = (massden   + atw(ibion)*enbeam(j))* convert2 !kg
       massden      = (massden   + masse*ene(j))*1.e6              !kg/m**3
       dmass(j)     = massden           !kg/m**3
       pressr(j)    = pressr(j) + ene(j)*te(j)
       pressr(j)    = pressr(j) * convert             !kev/cm**3 to J/M**3
       pbeam(j)     = 0.66*wbeam(j) * convert         !kev/cm**3 to J/M**3
       palpha(j)    = 0.66*walp(j)  * convert
       pressr(j)    = pressr(j) + pbeam(j)            !total pressure Pascals
       pressr(j)    = pressr(j) + palpha(j)
       sigpre(j)    = pressr(j) * 0.1                 !10 % uncertainty
       sibeam(j)    = ABS(rpress(j))
       dnbeam(j)    = enbeam(j) *1.e6                 !#/M**3
       sizeroj(j)   = ABS(rpress(j))
       rcapil       = rcapi(j)
       IF(cgs .EQ.1) rcapil = rcapil*100.             !rcapil in 1/m
       vzeroj(j)    = curden_loc(j)/(rmaj*rcapil)
       vzerobj(j)   = curboot(j)*1.e4/(rmaj*rcapil)
       IF(cxa .GT. 0.0)THEN
          jave =         totcur(1)/cxa
          vzeroj(j)    = vzeroj(j)/jave
          vzerobj(j)    = vzerobj(j)/jave
       ELSE
          vzeroj(j)    = 0.0
          vzerobj(j)    = 0.0
       ENDIF
    ENDDO
!    print *,'cgs = ',cgs
!    print *,'cxa = ',cxa
!    print *,'rmaj = ',rmaj
!    print *,'totcur(1) =', totcur(1)
!    print *,'rcapil(1), rcapi(nj) =',rcapil,rcapi(1)
!    print *,'dmass(1) =',dmass(1)

    CALL getioun(io12,42) 
    OPEN  (unit = io12,                                           &
            file = kin_efit_namelist(1:LEN_TRIM(kin_efit_namelist)),     &
            status = 'UNKNOWN')
    WRITE (unit = io12, fmt = '(3x, a)') runid_EFIT
    WRITE (unit = io12, nml = IN1)
    WRITE (unit = io12, nml = INWANT)

!   write individual outputs:  <U> means flux surface average of U
!     <R0**2/R**2>     ==> r2cap(1..nj)
!     <R**2>           ==> r2capi(1..nj)
!     <R>              ==> rcap(1..nj)
!     <1./R>           ==> rcapi(1..nj)
!     <Bt0**2/B**2>    ==> bsqncap
!     <B**2/Bt0**2>    ==> bsq_avg_cap
!     <B/Bt0>          ==> b_avg_cap
!     rmajorvec        ==> Rmajor on rho(ie r) grd at elevation of mag axis
!     bprmaj           ==> Bp at rmajorvec in KGAUSS if called from tport
!     btotrmaj         ==> total B at rmajorvec
!     pfscur           ==> Pfirsch - Schluter current on Rmajorgrid at outboard sde
!                           of mag netic axis at z = mag axis elevation
!     pfscur_avg      ==> <pfscur> flux surface average (Note that <B*pfscur> =0)
    rmloc = rmajor
    btloc = btor
    r_loc(:)       =  r(:)
    r_loc(1)       =  1.                                 !reset to 0.0 below
    psir_loc(:)    =  psir(:)
    r2capi_loc(:)  =  r2capi(:)
    rcap_loc(:)    =  rcap(:)
    rcapi_loc(:)   =  rcapi(:)
    I_loc(:)       =  twopi*rbp(:)/(u0*fcap(:))     !assumes rbp in tesla m
                                                   !I = 2pi*Gcap*hcap*rho*bp0/u0
                                                   !is the MKS relation
                                                   !note that rbp = FCAP*GCAP*HCAP*RHO*Bp0
    rmajor_loc(:)  = rmajorvec(:)
    bprmaj_loc(:)  = bprmaj(:)
    btotrmaj_loc(:) = btotrmaj(:)  !          btotl on rmajr at z = mag axis
    dpsidrho_loc(:) = rmloc*rbp(:)/         &
                       (fcap(:)*gcap(:)*hcap(:)*r_loc(:))      !gives dpsidrho_loc(1) =0.0
    curb_loc(:)    = curbe(:) + curbi(:) + curbet(:)
    curboot_loc(:)    = curboot(:)
    currf_loc(:)      = currf(:)
    curohm_loc(:)     = curohm(:)
    curpar_loc(:)     = curpar_soln(:)            !< J dot B/Bt0>
    q_loc(:)           =  ABS(q(:))
    pprim_loc(:)      = pprim(:)
    if(cgs .eq. 1)then
       rmloc         = rmajor*0.01                                    !cm to m
       btloc         = btor/1.e4                                      !gauss to tesla
       r_loc(:)      = r(:)*0.01                                    !cm to m
       r_loc(1)      = 1.                                           !reset to 0.0 below
       psir_loc(:)   = psir(:)*psimks                            !kgauss cm**2/rad to volt sec/rad
       r2cap_loc(:)  = r2cap(:)
       r2capi_loc(:) = r2capi(:)*1.e-4                         !cm**2 to m**2
       rcap_loc(:)   = rcap(:)*0.01                              !cm to m
       rcapi_loc(:)  = rcapi(:)*100.                            !1/cm to 1/m
                                                              !rbp is in gauss cm 
       I_loc(:)        =  5.*rbp(:)/fcap(:)                     !factor of 5 converts to amps
       dpsidrho_loc(:) = rmloc*rbp(:)*1.e-6/         &
                       (fcap(:)*gcap(:)*hcap(:)*r_loc(:)) 
       curb_loc(:)    = 1.e4*curb_loc(:)
       curboot_loc(:) = 1.e4*curboot_loc(:)
       currf_loc(:)   = 1.e4*currf_loc(:)
       curohm_loc(:)  = 1.e4*curohm_loc(:) 
       curpar_loc(:)  = 1.e4*curpar_loc(:)
       rmajor_loc(:)  = 1.e-2*rmajor_loc(:)
       bprmaj_loc(:)  = 1.e-1*bprmaj_loc(:)       !kgauss to tesla
       btotrmaj_loc(:) = 1.e-1*btotrmaj_loc(:)    !kgauss to tesla
       pprim_loc(:)    = 1.e7 * pprim_loc(:)      ! from gram/gauss*cm**3sec**2 to amps/m**3
    endif


    do  j=1,nj
        if(external_beam_cur .eq. 0)then
               curb_loc(j) = curbe(j) + curbi(j) + curbet(j)
        else
               curb_loc(j) = curb_external(j)
        endif
    enddo


    r_loc(1) =0.0
    dpdrho_loc(1) = 0.0
    dpdpsi_loc(1) = pprim_loc(1)  
    !use dp/dpsi = 2*c*q(1)/Bt0, where c is defined by P= P0 + c*r**2 +d*r**3 +...
    dpdpsi_loc(1) = q_loc(1)*2.*(pressr(2)-pressr(1))/r_loc(2)**2/btloc
    dpdrho_loc(nj) = (pressr(nj)-pressr(nj-1))/(r_loc(nj)-r_loc(nj-1)) ! backward difference
    dpdpsi_loc(nj) =  dpdrho_loc(nj)/dpsidrho_loc(nj)
    do  j=1,nj
        if(j .gt. 1 .and. j .lt. nj)then
           dpdrho_loc(j) = (pressr(j+1) - pressr(j-1))/(r_loc(j+1)-r_loc(j-1))
           dpdpsi_loc(j) = dpdrho_loc(j)/dpsidrho_loc(j)
        endif
           pfscur(j) = -(btloc*rmloc/fcap(j)/btotrmaj_loc(j))* (                     &
          1.- btotrmaj_loc(j)**2/(btloc*btloc*bsq_avg_cap(j)))
           pfscur(j) = pfscur(j)*dpdpsi_loc(j)
           pfscur_avg(j) = -(btloc*rmloc/fcap(j))*(1./(btloc*b_avg_cap(j))          &
                          -btloc*b_avg_cap(j)/(btloc*btloc*bsq_avg_cap(j)))*dpdpsi_loc(j)
           jpar_rmaj(j)  = curpar_loc(j)*btotrmaj_loc(j)/                             &
                         (btloc*btloc*bsq_avg_cap(j)) + pfscur(j)
    enddo

    WRITE (unit = io12,fmt='("*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')

    WRITE (unit = io12,fmt='("* nj,number of elements in following arrays ")')
    WRITE (unit = io12,fmt ='(2x,i6)')nj !jmp.den and similar changes in subsequenct routine

    WRITE (unit = io12,fmt='("* R0,m ")')
    WRITE (unit = io12,fmt=1)rmloc
  
    WRITE (unit = io12,fmt='("* Bt0, Tesla")')
    WRITE (unit = io12,fmt=1)btloc

    WRITE (unit = io12,fmt='("* psir, psi on rho(==r) grid  ,volt sec/rad")')
    WRITE (unit = io12,fmt=1)(psir_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* r (== rho) grid  ,m")')
    WRITE (unit = io12,fmt=1)(r_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* r2cap = <R0**2/R**2> ,dimensionless")')
    WRITE (unit = io12,fmt=1)(r2cap_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* r2capi = <R**2> ,m**2 ")')
    WRITE (unit = io12,fmt=1)(r2capi_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* rcap = <R> ,m")')
    WRITE (unit = io12,fmt=1)(rcap_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* rcapi = <1./R> ,m-1")')
    WRITE (unit = io12,fmt=1)(rcapi_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* bsqinvavg  = <Bt0**2/B**2> ,dimensionless")')
    WRITE (unit = io12,fmt=1)(bsqncap(j),j=1,nj)

    WRITE (unit = io12,fmt='("* bsq_avg  = <B**2/Bt0**2> ,dimensionless")')
    WRITE (unit = io12,fmt=1)(bsq_avg_cap(j),j=1,nj)

    WRITE (unit = io12,fmt='("* btot_avg  = <B/Bt0> ,dimensionless")')
    WRITE (unit = io12,fmt=1)(b_avg_cap(j),j=1,nj)
    if(src_ident == 0)THEN
       WRITE (unit = io12,fmt='("* curden = <Jphi R0/R> ==>(c/4piHr)*d/drho(G*H*rho*Bp0) " &
            &  " amps/m**2,from Faradys law, NOT from flux average")')
    ELSE
       WRITE(unit = io12,fmt='("* curden = <Jphi R0/R> : ",256(a))')temp_string
    ENDIF
    WRITE (unit = io12,fmt=1)(curden_loc(j),j=1,nj)
    WRITE (unit = io12,fmt='("* curpar = <J dot B/Bt0> ==>(c/4piFcap**2*H*r)*d/drho(F*G*H*rho*Bp0) " &
      &  " amps/m**2,from Faradys law, NOT from flux average")')
    WRITE (unit = io12,fmt=1)(curpar_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* curboot,amps/m**2==> <Jboot dot B/Bt0>  ")')
    WRITE (unit = io12,fmt=1)(curboot_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* curbeam, amps/m**2==> <Jbeam dot B/Bt0>  ")')
    WRITE (unit = io12,fmt=1)(curb_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* cur rf,amps/m**2 ==> <Jrf dot B/Bt0> ")')
    WRITE (unit = io12,fmt=1)(currf_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* residual ohmic current,amps/m**2 ==> <Johm dot B/Bt0>  ")')
    WRITE (unit = io12,fmt=1)(curohm_loc(j),j=1,nj)


    WRITE (unit = io12,fmt='("* Integrated current,INTG(2piHcap*rho*curden drho), amps")')
    WRITE (unit = io12,fmt=1)(I_loc(j),j=1,nj)


    WRITE (unit = io12,fmt='("* dpsi/drho ==> R0*Bp0, volt sec/m ==> tesla m ")')
    WRITE (unit = io12,fmt=1)(dpsidrho_loc(j),j=1,nj)
 

    WRITE (unit = io12,fmt='("* Pressure ==> sum (ni*Ti )+ne*Te + 0.66*Wbeam Nt/m**2 ")')
    WRITE (unit = io12,fmt=1)(pressr(j),j=1,nj)



    WRITE (unit = io12,fmt='("* dP/dpsi ==> (dP/drho)/(dpsi/drho) from TRANSPORT(NOT EFIT !!)")')
    WRITE (unit = io12,fmt=1)(dpdpsi_loc(j),j=1,nj)


    WRITE (unit = io12,fmt='("* dP/dpsi ==> pprim ")')
    WRITE (unit = io12,fmt=1)(pprim_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* q = ABS(rho*Bt0/(R*Bp0))")')
    WRITE (unit = io12,fmt=1)(q_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* hcap = (dV/drho)*(1./(4pi**2 * R0 *r))")')
    WRITE (unit = io12,fmt=1)(hcap(j),j=1,nj)

    WRITE (unit = io12,fmt='("* fcap = R0*Bt0/(R*BT),dimensionless ")')
    WRITE (unit = io12,fmt=1)(fcap(j),j=1,nj)

    WRITE (unit = io12,fmt='("* gcap = < (grad(r))**2 R0**2/R**2 >,dimensionless")')
    WRITE (unit = io12,fmt=1)(gcap(j),j=1,nj)


    WRITE (unit = io12,fmt='("* Rmajor,m at Z=zmag axis,over rho grid ")')
    WRITE (unit = io12,fmt=1)(rmajor_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* Bp,tesla  on rmajor grid at z = mag axis")')
    WRITE (unit = io12,fmt=1)(bprmaj_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("*B total,tesla on rmajor grid  at z = mag axis ")')
    WRITE (unit = io12,fmt=1)(btotrmaj_loc(j),j=1,nj)

    WRITE (unit = io12,fmt='("* Pfirsch Schluter parallel current on Rmajor, amps/m**2")')
    WRITE (unit = io12,fmt=1)(pfscur(j),j=1,nj)


    WRITE (unit = io12,fmt='("* Pfirsch Schluter flux surface average parallel  current on Rmajor,amps/m**2, <JPS> (NOT <JPS dot B > == 0)")')
    WRITE (unit = io12,fmt=1)(pfscur_avg(j),j=1,nj)


    WRITE (unit = io12,fmt='("* local J parallel  current on Rmajor, amps/m**2")')
    WRITE (unit = io12,fmt=1)(jpar_rmaj(j),j=1,nj)








    CLOSE (unit = io12)
    CALL giveupus(io12)

     first_call = .FALSE.    

    RETURN
1   format(5(2x,1pe16.8))
  END SUBROUTINE  wrt_kin_efit_naml



  SUBROUTINE  read_kin_efit(kin_file)
! -------------------------------------------------------------------------
!    subroutine processes kinetic efit files prodoced by wrt_kin_efit_naml
!    and creates eqdsk
!    INPUT : kin_file, character variable, give name of kinetic efit file 
!    to read
! -------------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER i,j,io12,itime,npress,nbeam,nmass,kzeroj,cgs
    REAL(DP) pnbeam,convert,convert2,massden,masse,rmaj,rcapil, &
            jave,cxa,btloc

    CHARACTER(len = 12) ishota,intfl,itimea
    CHARACTER(len = 256) dumy
    CHARACTER(*) kin_file

    NAMELIST /IN1/       ishot,itime,pnbeam,npress,nbeam,      &
                         rpress,pressr,sibeam,pbeam,dnbeam,    &
                         dmass,sigpre

    NAMELIST /INWANT/    kzeroj,vzeroj,sizeroj,vzerobj




    convert  = 1.6021765e-10       !kev/cm**3 to J/M**3
    convert2 = 1.6605387e-27       ! amu to kg
    masse    = 9.10938970e-31      ! electron mass   (kg)


 


    CALL getioun(io12,42)
    OPEN  (unit = io12, file = kin_file,status = 'UNKNOWN')


    READ (unit = io12, nml = IN1)
    READ (unit = io12, nml = INWANT)
    READ (unit = io12,fmt='(a)')dumy

!   READ (unit = io12,fmt='("* nj,number of elements in following arrays ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt ='(2x,i6)')nj

    CALL kine_efit_alloc(nj)

!   read individual inputs:  <U> means flux surface average of U
!     <R0**2/R**2>     ==> r2cap(1..nj)
!     <R**2>           ==> r2capi(1..nj)
!     <R>              ==> rcap(1..nj)
!     <1./R>           ==> rcapi(1..nj)
!     <Bt0**2/B**2>    ==> bsqncap
!     <B**2/Bt0**2>    ==> bsq_avg_cap
!     <B/Bt0>          ==> b_avg_cap
!     rmajorvec        ==> Rmajor on rho(ie r) grd at elevation of mag axis
!     bprmaj           ==> Bp at rmajorvec in KGAUSS if called from tport
!     btotrmaj         ==> total B at rmajorvec
!     pfscur           ==> Pfirsch - Schluter current on Rmajorgrid at outboard sde
!                           of mag netic axis at z = mag axis elevation
!     pfscur_avg       ==> <pfscur> flux surface average (Note that <B*pfscur> =0)

 ! remove the following 
 
  
!    READ (unit = io12,fmt='("* R0,m ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)R0
  
!    READ (unit = io12,fmt='("* Bt0, Tesla")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)BT0

!    READ (unit = io12,fmt='("* psir, psi on rho(==r) grid  ,volt sec/rad")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(psir_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* r (== rho) grid  ,m")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(r_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* r2cap = <R0**2/R**2> ,dimensionless")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(r2cap_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* r2capi = <R**2> ,m**2 ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(r2capi_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* rcap = <R> ,m")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(rcap_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* rcapi = <1./R> ,m-1")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(rcapi_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* bsqinvavg  = <Bt0**2/B**2> ,dimensionless")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(bsqncap(j),j=1,nj)

!    READ (unit = io12,fmt='("* bsq_avg  = <B**2/Bt0**2> ,dimensionless")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(bsq_avg_cap(j),j=1,nj)

!    READ (unit = io12,fmt='("* btot_avg  = <B/Bt0> ,dimensionless")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(b_avg_cap(j),j=1,nj)

!    READ (unit = io12,fmt='("* curden = <Jphi R0/R> ==>(c/4piHr)*d/drho(G*H*rho*Bp0) " &
!      &  " amps/m**2,from Faradys law, NOT from flux average")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(curden_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* curpar = <J dot B/Bt0> ==>(c/4piFcap**2*H*r)*d/drho(F*G*H*rho*Bp0) " &
!      &  " amps/m**2,from Faradys law, NOT from flux average")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(curpar_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* curboot,amps/m**2==> <Jboot dot B/Bt0>  ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(curboot_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* curbeam, amps/m**2==> <Jbeam dot B/Bt0>  ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(curb_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* cur rf,amps/m**2 ==> <Jrf dot B/Bt0> ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(currf_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* residual ohmic current,amps/m**2 ==> <Johm dot B/Bt0>  ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(curohm_loc(j),j=1,nj)


!    READ (unit = io12,fmt='("* Integrated current,INTG(2piHcap*rho*curden drho), amps")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(I_loc(j),j=1,nj)


!    READ (unit = io12,fmt='("* dpsi/drho ==> R0*Bp0, volt sec/m ==> tesla m ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(dpsidrho_loc(j),j=1,nj)
 

!    READ (unit = io12,fmt='("* Pressure ==> sum (ni*Ti )+ne*Te + 0.66*Wbeam Nt/m**2 ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(pressr(j),j=1,nj)



!    READ (unit = io12,fmt='("* dP/dpsi ==> (dP/drho)/(dpsi/drho) from TRANSPORT(NOT EFIT !!)")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(dpdpsi_loc(j),j=1,nj)


!    READ (unit = io12,fmt='("* dP/dpsi ==> pprim ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(pprim_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* q = ABS(rho*Bt0/(R*Bp0))")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(q_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* hcap = (dV/drho)*(1./(4pi**2 * R0 *r))")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(hcap(j),j=1,nj)

!    READ (unit = io12,fmt='("* fcap = R0*Bt0/(R*BT),dimensionless ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(fcap(j),j=1,nj)

!    READ (unit = io12,fmt='("* gcap = < (grad(r))**2 R0**2/R**2 >,dimensionless")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(gcap(j),j=1,nj)


!    READ (unit = io12,fmt='("* Rmajor,m at Z=zmag axis,over rho grid ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(rmajor_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* Bp,tesla  on rmajor grid at z = mag axis")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(bprmaj_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("*B total,tesla on rmajor grid  at z = mag axis ")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(btotrmaj_loc(j),j=1,nj)

!    READ (unit = io12,fmt='("* Pfirsch Schluter parallel current on Rmajor, amps/m**2")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(pfscur(j),j=1,nj)


!    READ (unit = io12,fmt='("* Pfirsch Schluter flux surface average parallel  current on Rmajor,amps/m**2, <JPS> (NOT <JPS dot B > == 0)")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(pfscur_avg(j),j=1,nj)


!    READ (unit = io12,fmt='("* local J parallel  current on Rmajor, amps/m**2")')
    READ (unit = io12,fmt='(a)')dumy
    READ (unit = io12,fmt=1)(jpar_rmaj(j),j=1,nj)

    CLOSE (unit = io12)
    CALL giveupus(io12)

    RETURN
1   format(5(2x,1pe16.8))
  END SUBROUTINE  read_kin_efit


  SUBROUTINE kine_efit_alloc(ndim)

    IMPLICIT NONE
    INTEGER(I4B) ndim,ns



    IF(ALLOCATED(palpha))THEN
       ns = SIZE(palpha)
       IF(ns .NE. nj)DEALLOCATE(palpha)
    ENDIF
    IF(.NOT. ALLOCATED(palpha))ALLOCATE(palpha(ndim))
 

    IF(ALLOCATED(r_loc))THEN
       ns = SIZE(r_loc)
       IF(ns .NE. nj)DEALLOCATE(r_loc)
    ENDIF
    IF(.NOT. ALLOCATED(r_loc))ALLOCATE(r_loc(ndim))

    IF(ALLOCATED(psir_loc))THEN
       ns = SIZE(psir_loc)
       IF(ns .NE. nj)DEALLOCATE(psir_loc)
    ENDIF
    IF(.NOT. ALLOCATED(psir_loc))ALLOCATE(psir_loc(ndim))

    IF(ALLOCATED(rcap_loc))THEN
       ns = SIZE(rcap_loc)
       IF(ns .NE. nj)DEALLOCATE(rcap_loc)
    ENDIF
    IF(.NOT. ALLOCATED(rcap_loc))ALLOCATE(rcap_loc(ndim))

    IF(ALLOCATED(r2cap_loc))THEN
       ns = SIZE(r2cap_loc)
       IF(ns .NE. nj)DEALLOCATE(r2cap_loc)
    ENDIF
    IF(.NOT. ALLOCATED(r2cap_loc))ALLOCATE(r2cap_loc(ndim))

    IF(.NOT. ALLOCATED(rcapi_loc))ALLOCATE(rcapi_loc(ndim))
    IF(ALLOCATED(rcapi_loc))THEN
       ns = SIZE(rcapi_loc)
       IF(ns .NE. nj)DEALLOCATE(rcapi_loc)
    ENDIF
    IF(.NOT. ALLOCATED(rcapi_loc))ALLOCATE(rcapi_loc(ndim))

    IF(ALLOCATED(curden_loc))THEN
       ns = SIZE(curden_loc)
       IF(ns .NE. nj)DEALLOCATE(curden_loc)
    ENDIF
    IF(.NOT. ALLOCATED(curden_loc))ALLOCATE(curden_loc(ndim))

    IF(ALLOCATED(r2capi_loc))THEN
       ns = SIZE(r2capi_loc)
       IF(ns .NE. nj)DEALLOCATE(r2capi_loc)
    ENDIF
    IF(.NOT. ALLOCATED(r2capi_loc))ALLOCATE(r2capi_loc(ndim))

    IF(ALLOCATED(q_loc))THEN
       ns = SIZE(q_loc)
       IF(ns .NE. nj)DEALLOCATE(q_loc)
    ENDIF
    IF(.NOT. ALLOCATED(q_loc))ALLOCATE(q_loc(ndim))


    IF(ALLOCATED(dpsidrho_loc))THEN
       ns = SIZE(dpsidrho_loc)
       IF(ns .NE. nj)DEALLOCATE(dpsidrho_loc)
    ENDIF
    IF(.NOT. ALLOCATED(dpsidrho_loc))ALLOCATE(dpsidrho_loc(ndim))


    IF(ALLOCATED(I_loc))THEN
       ns = SIZE(I_loc)
       IF(ns .NE. nj)DEALLOCATE(I_loc)
    ENDIF
    IF(.NOT. ALLOCATED(I_loc))ALLOCATE(I_loc(ndim))

    IF(ALLOCATED(curb_loc))THEN
       ns = SIZE(curb_loc)
       IF(ns .NE. nj)DEALLOCATE(curb_loc)
    ENDIF
    IF(.NOT. ALLOCATED(curb_loc))ALLOCATE(curb_loc(ndim))

    IF(ALLOCATED(curboot_loc))THEN
       ns = SIZE(curboot_loc)
       IF(ns .NE. nj)DEALLOCATE(curboot_loc)
    ENDIF
    IF(.NOT. ALLOCATED(curboot_loc))ALLOCATE(curboot_loc(ndim))

    IF(ALLOCATED(currf_loc))THEN
       ns = SIZE(currf_loc)
       IF(ns .NE. nj)DEALLOCATE(currf_loc)
    ENDIF
    IF(.NOT. ALLOCATED(currf_loc))ALLOCATE(currf_loc(ndim))

    IF(ALLOCATED(curohm_loc))THEN
       ns = SIZE(curohm_loc)
       IF(ns .NE. nj)DEALLOCATE(curohm_loc)
    ENDIF
    IF(.NOT. ALLOCATED(curohm_loc))ALLOCATE(curohm_loc(ndim))

    IF(ALLOCATED(rmajor_loc))THEN
       ns = SIZE(rmajor_loc)
       IF(ns .NE. nj)DEALLOCATE(rmajor_loc)
    ENDIF
    IF(.NOT. ALLOCATED(rmajor_loc))ALLOCATE(rmajor_loc(ndim))

    IF(ALLOCATED(bprmaj_loc))THEN
       ns = SIZE(bprmaj_loc)
       IF(ns .NE. nj)DEALLOCATE(bprmaj_loc)
    ENDIF
    IF(.NOT. ALLOCATED(bprmaj_loc))ALLOCATE(bprmaj_loc(ndim))

    IF(ALLOCATED(btotrmaj_loc))THEN
       ns = SIZE(btotrmaj_loc)
       IF(ns .NE. nj)DEALLOCATE(btotrmaj_loc)
    ENDIF
    IF(.NOT. ALLOCATED(btotrmaj_loc))ALLOCATE(btotrmaj_loc(ndim))

    IF(ALLOCATED(pfscur))THEN
       ns = SIZE(pfscur)
       IF(ns .NE. nj)DEALLOCATE(pfscur)
    ENDIF
    IF(.NOT. ALLOCATED(pfscur))ALLOCATE(pfscur(ndim))

    IF(ALLOCATED(pfscur_avg))THEN
       ns = SIZE(pfscur_avg)
       IF(ns .NE. nj)DEALLOCATE(pfscur_avg)
    ENDIF
    IF(.NOT. ALLOCATED(pfscur_avg))ALLOCATE(pfscur_avg(ndim))

    IF(ALLOCATED(jpar_rmaj))THEN
       ns = SIZE(jpar_rmaj)
       IF(ns .NE. nj)DEALLOCATE(jpar_rmaj)
    ENDIF
    IF(.NOT. ALLOCATED(jpar_rmaj))ALLOCATE(jpar_rmaj(ndim))

    IF(ALLOCATED(curpar_loc))THEN
       ns = SIZE(curpar_loc)
       IF(ns .NE. nj)DEALLOCATE(curpar_loc)
    ENDIF
    IF(.NOT. ALLOCATED(curpar_loc))ALLOCATE(curpar_loc(ndim))

    IF(ALLOCATED(dpdpsi_loc))THEN
       ns = SIZE(dpdpsi_loc)
       IF(ns .NE. nj)DEALLOCATE(dpdpsi_loc)
    ENDIF
    IF(.NOT. ALLOCATED(dpdpsi_loc))ALLOCATE(dpdpsi_loc(ndim))

    IF(ALLOCATED(dpdrho_loc))THEN
       ns = SIZE(dpdrho_loc)
       IF(ns .NE. nj)DEALLOCATE(dpdrho_loc)
    ENDIF
    IF(.NOT. ALLOCATED(dpdrho_loc))ALLOCATE(dpdrho_loc(ndim))

    IF(ALLOCATED(pprim_loc))THEN
       ns = SIZE(pprim_loc)
       IF(ns .NE. nj)DEALLOCATE(pprim_loc)
    ENDIF
    IF(.NOT. ALLOCATED(pprim_loc))ALLOCATE(pprim_loc(ndim))

    RETURN 
 END SUBROUTINE kine_efit_alloc

END MODULE kinetic_efit   
