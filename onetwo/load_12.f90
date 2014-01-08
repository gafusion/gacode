
SUBROUTINE load_12(io)
! ----------------------------------------------------------------------
!   read output from nubeam code on unit io
!   load respective Onetwo variables
! ---------------------------------------------------------------------
  USE transp
  USE param, ONLY :kj,kb
  USE numbrs,ONLY : nj
  USE mesh,ONLY :  roa,ra,r

  USE nub, ONLY : hibr,hibrs
  USE string_util,ONLY : to_upper_case1,readblock,read_scalar,                &
                         read_scalar1,read_scalar2


  !following are defined in fusion.mod but are not used otherwise in Onetwo::
  USE fusion, ONLY :  beam_thermalntot, beam_beamnstot, beam_beamntot  
  USE sourc ,  ONLY : curbe,curbi
  USE verbose, ONLY : fusionvb
  USE machin,  ONLY : rmajor,btor

  IMPLICIT NONE
  INTEGER io,j,k,ll,kk,kjm1,nsfast,ngmax, &
          nbeams,n_3,n_miefi,i,i1,i2,i3,i4,nterms,load_calls
  CHARACTER(len=256) :: line 
  CHARACTER(len=2) chrf
  CHARACTER(len =1) TYPE
  CHARACTER (len = 132)intfile,numbr
  CHARACTER *60 fplotname
  REAL *8 zbuf(kj),bbntnst(8),value_intg,cvrt,timei,timee,rbt
  REAL *8 bbntot,bbnta(14)            !14 reaction rates
  REAL *8 btneut,btnta(14)            !14 reaction rates
  REAL *8 sftot(14),sum
  CHARACTER*20, DIMENSION(:), ALLOCATABLE :: zflbls
  CHARACTER(len =4) zthlbls(5)
  LOGICAL skip_summary,swap_lbl
  INTEGER,SAVE :: torque_ct
  INTEGER fi_species
  DATA zthlbls / 'H','D','T','He3','He4' /
  DATa load_calls /0/
  kjm1 = kj-1
  skip_summary = .FALSE. ;swap_lbl = .FALSE.


  load_calls = load_calls +1
! position the file for the first read:
  DO WHILE(1 .GT. 0)
     READ(io,FMT='(A)',ERR = 14,END=14)line
     CALL to_upper_case1(line)
     line  = ADJUSTL(line)
     ll = LEN_TRIM(line)  
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
     IF(k == 0)k=ll                                 !string does not have '!'
     kk=0
     kk= INDEX(line(1:k),'NSFAST')
     IF(kk .EQ. 0)CYCLE
     chrf=line(10:11)
     WRITE(intfile,FMT='(a)')chrf
     READ(intfile,FMT='(i2)')nsfast
     IF(.NOT. ALLOCATED(zflbls))ALLOCATE(zflbls(nsfast))
     DO j=1,nsfast
       READ(io,FMT='(a)')zflbls(j)
       zflbls(j) = ADJUSTL(zflbls(j))
     ENDDO
     EXIT
14   CALL STOP('load_12 file read error 1',1)
  ENDDO


  DO WHILE(1 .GT. 0)
     READ(io,FMT='(A)',ERR = 15,END=15)line
     CALL to_upper_case1(line)
     line  = ADJUSTL(line)
     ll = LEN_TRIM(line)  
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
     IF(k == 0)k=ll                                 !string does not have '!'
     kk=0
     kk= INDEX(line(1:k),'NGMAX')
     IF(kk .EQ. 0)CYCLE
     chrf=line(10:11)
     WRITE(intfile,FMT='(a)')chrf
     READ(intfile,FMT='(i2)')ngmax
     EXIT
15   CALL STOP('load_12 file read error 1',2)
  ENDDO

  DO WHILE(1 .GT. 0)
     READ(io,FMT='(A)',ERR = 16,END=16)line
     CALL to_upper_case1(line)
     line  = ADJUSTL(line)
     ll = LEN_TRIM(line)  
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
     IF(k == 0)k=ll                                 !string does not have '!'
     kk=0
     kk= INDEX(line(1:k),'NBEAMS')
     IF(kk .EQ. 0)CYCLE
     chrf=line(10:11)
     WRITE(intfile,FMT='(a)')chrf
     READ(intfile,FMT='(i2)')nbeams
     EXIT
16   CALL STOP('load_12 file read error 1',3)
  ENDDO

  DO WHILE(1 .GT. 0)
     READ(io,FMT='(A)',ERR = 17,END=17)line
     CALL to_upper_case1(line)
     line  = ADJUSTL(line)
     ll = LEN_TRIM(line)  
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
     IF(k == 0)k=ll                                 !string does not have '!'
     kk=0
     kk= INDEX(line(1:k),'N_3')
     IF(kk .EQ. 0)CYCLE
     chrf=line(10:11)
     WRITE(intfile,FMT='(a)')chrf
     READ(intfile,FMT='(i2)')n_3
     EXIT
17   CALL STOP('load_12 file read error 1',4)
  ENDDO

  DO WHILE(1 .GT. 0)
     READ(io,FMT='(A)',ERR = 18,END=18)line
     CALL to_upper_case1(line)
     line  = ADJUSTL(line)
     ll = LEN_TRIM(line)  
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
     IF(k == 0)k=ll                                 !string does not have '!'
     kk=0
     kk= INDEX(line(1:k),'N_MIEFI')
     IF(kk .EQ. 0)CYCLE
     chrf=line(11:12)
     WRITE(intfile,FMT='(a)')chrf
     READ(intfile,FMT='(i2)')n_miefi
     EXIT
18   CALL STOP('load_12 file read error 1',5)
  ENDDO



  DO WHILE(1 .GT. 0)
     READ(io,FMT='(A)',ERR = 27,END=27)line
     CALL to_upper_case1(line)
     line  = ADJUSTL(line)
     ll = LEN_TRIM(line)  
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
     IF(k == 0)k=ll                                 !string does not have '!'
     kk=0
     kk= INDEX(line(1:k),'TIMEI')
     IF(kk .EQ. 0)CYCLE
     chrf=line(11:12)
     WRITE(intfile,FMT='(a)')line(9:22)
     READ(intfile,FMT='(1pe14.8)',ERR=271)timei
     EXIT
271  PRINT *,'load_12 file read error 1'
     PRINT *,'reading line :',line
     PRINT *,'chrf = :',chrf
     PRINT *,'IGNORED timei'
     EXIT
27   CALL STOP('load_12 file read error 1',5)
  ENDDO


  DO WHILE(1 .GT. 0)
     READ(io,FMT='(A)',ERR = 28,END=28)line
     CALL to_upper_case1(line)
     line  = ADJUSTL(line)
     ll = LEN_TRIM(line)  
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
     IF(k == 0)k=ll                                 !string does not have '!'
     kk=0
     kk= INDEX(line(1:k),'TIMEE')
     IF(kk .EQ. 0)CYCLE
     chrf=line(11:12)
     WRITE(intfile,FMT='(a)')line(9:22)
     READ(intfile,FMT='(1pe14.8)',ERR=281)timee
     EXIT
281  PRINT *,'load_12 file read error 1'
     PRINT *,'reading line :',line
     PRINT *,'chrf = :',chrf
     PRINT *,'IGNORED timee'
     EXIT
28   CALL STOP('load_12 file read error 1',5)
  ENDDO
!     PRINT *,'nsfast,n_3,nbeams,n_miefi,ngmax = ',nsfast,n_3,nbeams,n_miefi,ngmax
     PRINT *,'timei,timee =',timei,timee
      nubeam_restart_time = timee


  DO WHILE(1 .GT. 0)
     READ(io,FMT='(A)',ERR = 25,END=25)line
     CALL to_upper_case1(line)
     line  = ADJUSTL(line)
     ll = LEN_TRIM(line)  
     k = INDEX(line,'!')
     IF(k == 1)CYCLE                                !ignore comment lines
     IF(k == 0)k=ll                                 !string does not have '!'
     kk=0
     kk= INDEX(line(1:k),'VOLUME')
     IF(kk .NE. 0)EXIT
     CYCLE
25   CALL STOP('load_12 file read error 1',6)
  ENDDO


    TYPE ='V'
    cvrt = 1.e6
    skip_summary = .TRUE.
      CALL readblock('VOLUME',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                     i4,skip_summary,swap_lbl,zbuf,value_intg)
!    print *,'volume =',zbuf(1:nj)



    TYPE =' '
    cvrt = 1.e4
    skip_summary = .TRUE.
      CALL readblock('AREA',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                     i4,skip_summary,swap_lbl,zbuf,value_intg)
!    print *,'area =',zbuf(1:nj)


    TYPE =' '
    cvrt = 1.
    skip_summary = .TRUE.
      CALL readblock('RBT',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                     i4,skip_summary,swap_lbl,zbuf,value_intg)
!    print *,'RBT =',zbuf(1:nj)
     rbt = zbuf(1) 


    TYPE =' '
    cvrt = 1.
    skip_summary = .TRUE.
      CALL readblock('OMEGAG',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                     i4,skip_summary,swap_lbl,zbuf,value_intg)
!    print *,'omegag =',zbuf(1:nj)




    TYPE ='V '
    cvrt = 1.
    skip_summary = .TRUE.
      CALL readblock('ROAZ',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                     i4,skip_summary,swap_lbl,zbuf,value_intg)
    zbuf(1) =0.0     !make sure extrapolation hits exacty zero here
!   print *,'ROAZ =',zbuf(1:nj)

    TYPE =' '
    cvrt = 1.
    skip_summary = .TRUE.
      CALL readblock('RINV2AVG',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                     i4,skip_summary,swap_lbl,zbuf,value_intg)
!    print *,'R2INVAVG =',zbuf(1:nj)

    TYPE =' '
    cvrt = 1.
    skip_summary = .TRUE.
      CALL readblock('RINVAVG',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                     i4,skip_summary,swap_lbl,zbuf,value_intg)
!    print *,'RINVAVG =',zbuf(1:nj)

    TYPE =' '
    cvrt = 1.
    skip_summary = .TRUE.
      CALL readblock('B2AVG',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                     i4,skip_summary,swap_lbl,zbuf,value_intg)
!    print *,'B2AVG =',zbuf(1:nj)


!---------------------------------------------------------------
!  output block: fusion
!---------------------------------------------------------------
!itype1 = 1 beams, itype1 =2 FP, itype1 =3 rf
!gives the source of the first reagent.
!the second reagent is either beam ion 
!(output to "bb" arrays),or thermal target 
!(outupt to "bt" arrays)
!Fast ion - Fast ion reactions
! where the 2nd reagent is a fusion product
! or RF tail ion are currently not  computed
! ireact= 1  itype1= 1  reagent1,2= D   T  
! ireact= 2  itype1= 1  reagent1,2= D   HE3
! ireact= 3  itype1= 1  reagent1,2= D   D:P
! ireact= 4  itype1= 1  reagent1,2= D   D:N
! ireact= 5  itype1= 1  reagent1,2= T   T  
! ireact= 6  itype1= 1  reagent1,2= T   HE3
! ireact= 7  itype1= 1  reagent1,2= T   D  
! ireact= 8  itype1= 1  reagent1,2= HE3 D  
! ireact= 9  itype1= 2  reagent1,2= T   D  
! ireact= 10  itype1= 2  reagent1,2= HE3 D  
! ireact= 11  itype1= 3  reagent1,2= T   D  
! ireact= 12  itype1= 3  reagent1,2= HE3 D  
! ireact= 13  itype1= 3  reagent1,2= D   T  
! ireact= 14  itype1= 3  reagent1,2= D   HE3      !bbnta(14) comes from here
   fus_tablebb_12(1) = 'beam -beam   D T reaction'     !in module transp
   fus_tablebb_12(2) = 'beam - beam  D He3  reaction'
   fus_tablebb_12(3) = 'beam - beam  D(D,p)  reaction'
   fus_tablebb_12(4) = 'beam - beam  D(D,n)  reaction'
   fus_tablebb_12(5) = 'beam - beam  T T  reaction'
   fus_tablebb_12(6) = 'beam - beam  T He3  reaction'
   fus_tablebb_12(7) = fus_tablebb_12(1) 
   fus_tablebb_12(nfusbb) = fus_tablebb_12(2)
!  nfusbb =8 (transp.f90)




!  first load beam_beam reaction rates

   ! BEAM-BEAM NEUTRON EMISSION 1D PROFILE, by reaction n/sec (integrated profile)
   ! ALSO: INFO ON FUSION REACTIONS WHICH DO NOT MAKE NEUTRONS

    TYPE ='V'
    cvrt = 1.e-6
    skip_summary = .FALSE.
    DO WHILE(1 .GT. 0)
      CALL readblock('BBNTNS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                         i4,skip_summary,swap_lbl,zbuf,value_intg)
      IF(i2 .LE. 0)EXIT
      IF(i2 .GT. 5)CYCLE
      IF(i2 == 1)beam_beamdtn_nub(1:nj)  = zbuf(1:nj)
      IF(i2 == 3)beam_beamddp_nub(1:nj)  = zbuf(1:nj) 
      IF(i2 == 4)beam_beamddn_nub(1:nj)  = zbuf(1:nj)
      IF(i2 == 5)beam_beamtt2n_nub(1:nj) = zbuf(1:nj)
      WRITE(intfile,FMT='(i1)')i2
      READ(intfile,FMT='(A)')chrf
      chrf = ADJUSTL(chrf)
!      fplotname ='bbntns' //'_reaction_'//chrf
!      CALL plot_dumpd(r,nj,zbuf,fplotname)
    ENDDO




   cvrt = 1.e-6
   TYPE ='V'
!   bbntnt(mj)  BEAM-BEAM TOTAL NEUTRON EMISSION 1D PROFILE
!   only includes reactions that produce neutrons.
!   reactiions that produce protons are not included here
   CALL readblock('BBNTNT',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                       i4,skip_summary,swap_lbl,zbuf,value_intg)
   beam_beamntot(1:nj) = zbuf(1:nj)

   IF(fusionvb ==1)THEN
     PRINT *,'BEAM-BEAM TOTAL NEUTRON EMISSION 1D PROFILE'
     PRINT *,beam_beamntot(1:nj)
   ENDIF

! next load beam - thermal raction rates


    TYPE ='V'
    cvrt = 1.e-6
    beam_thermaltth_df_nub(1:nj) = 0.0D0
    beam_thermalddp_nub(1:nj)    = 0.0D0
    beam_thermalddn_nub(1:nj)    = 0.0D0
    beam_thermaltt2n_nub(1:nj)   = 0.0D0
    beam_thermaldth_tf_nub(1:nj) = 0.0D0
    DO WHILE(1 .GT. 0 )
      CALL readblock('BTNTNS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                     i4,skip_summary,swap_lbl,zbuf,value_intg)
      IF(i2 .LE. 0)EXIT
      IF(i2 .GT. 5)CYCLE
      IF(i2 == 1)beam_thermaltth_df_nub(1:nj) = zbuf(1:nj)
      IF(i2 == 3)beam_thermalddp_nub(1:nj)    = zbuf(1:nj)
      IF(i2 == 4)beam_thermalddn_nub(1:nj)    = zbuf(1:nj)
      IF(i2 == 5)beam_thermaltt2n_nub(1:nj)   = zbuf(1:nj)
      IF(i2 == 7)beam_thermaldth_tf_nub(1:nj) = zbuf(1:nj)
    ENDDO






   cvrt = 1.e-6
   TYPE ='V'
!   btntnt(mj)  BEAM-THERMAL  TOTAL NEUTRON EMISSION 1D PROFILE
   CALL readblock('BTNTNT',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                       i4,skip_summary,swap_lbl,zbuf,value_intg)
   beam_thermalntot(1:nj) = zbuf(1:nj)


! TOTAL BEAM-BEAM NEUTRON PRODUCTION, bbntot n/sec :
   CALL read_scalar('BBNTOT',line,io,beam_beamnstot)
!     print *,'BBNTOT =',beam_beamnstot


! read 14 beam-beam  scalar reaction rates:
      beam_beam_dtntot_nub  =  0.0     !neuts/sec,d(t,n)
      beam_beam_ddntot_nub  =  0.0    !neuts/sec,d(d,n)
      beam_beam_ddptot_nub  =  0.0    !neuts/sec,d(d,p)
      beam_beam_tt2ntot_nub =  0.0    !neuts/sec,t(t,2n)
  DO j=1,14
      CALL read_scalar1('BBNTA',line,io,i2,zflbls,bbnta(j))
      IF(i2 .EQ. 1) beam_beam_dtntot_nub  = bbnta(j)     !neuts/sec,d(t,n)
      IF(i2 .EQ. 4) beam_beam_ddntot_nub  = bbnta(j)     !neuts/sec,d(d,n)
      IF(i2 .EQ. 3) beam_beam_ddptot_nub  = bbnta(j)     !neuts/sec,d(d,p)
      IF(i2 .EQ. 5) beam_beam_tt2ntot_nub = bbnta(j)     !neuts/sec,t(t,2n)
  ENDDO

!  btneut  ! n/sec
! BEAM TARGET NEUTRONS PER SECOND (whatever that means)
   CALL read_scalar('BTNEUT',line,io,btneut)
!     print *,'BTNEUT =',btneut

! read 14 beam-thermal  scalar reaction rates:
      beam_thermal_dtntot_nub    = 0.0D0     !neuts/sec,d(t,n)
      beam_thermal_ddntot_nub    = 0.0D0     !neuts/sec,d(d,n)
      beam_thermal_ddptot_nub    = 0.0D0     !neuts/sec,d(d,p)
      beam_thermal_tt2ntot_nub   = 0.0D0     !neuts/sec,t(t,2n)
  DO j=1,14
      CALL read_scalar1('BTNTA',line,io,i2,zflbls,btnta(j))
!      print *,'j,btnta(j) =',j,btnta(j)
      IF(i2 .EQ. 1)beam_thermal_dtntot_nub = btnta(j)     !neuts/sec,d(t,n)
      IF(i2 .EQ. 4)beam_thermal_ddntot_nub = btnta(j)     !neuts/sec,d(d,n)
      IF(i2 .EQ. 3)beam_thermal_ddptot_nub = btnta(j)     !neuts/sec,d(d,p)
      IF(i2 .EQ. 5)beam_thermal_tt2ntot_nub = btnta(j)    !neuts/sec,t(t,2n)
  ENDDO

! read 14 total   scalar reaction rates:
  DO j=1,14
      CALL read_scalar1('SFTOT',line,io,i2,zflbls,sftot(j))
  ENDDO

!  read beam densities


    TYPE ='V'
    cvrt = 1.e-6
    enbeam_nub(:) =0.0d0
    enbeam_intg_nub = 0.0d0
    fi_species =0
    DO j = 1, nsfast
      CALL readblock('BDENSS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,  &
               cvrt,i1,i2,i3,i4,skip_summary,swap_lbl,zbuf,value_intg)
      IF(i2 .GT.0)THEN
         fi_species = fi_species + 1
         enbeam_nub(:) = enbeam_nub(:) + zbuf(:) 
         enbeam_species(:,fi_species) =  zbuf(:)
         enbeam_intg_nub = enbeam_intg_nub + value_intg
      ENDIF
      IF(i2 .LT. 0)CYCLE          !beam has decayed away
      WRITE(intfile,FMT='(i1)')i2
      READ(intfile,FMT='(A)')chrf
      chrf = ADJUSTL(chrf)
!      fplotname ='bdenss' //zflbls(i2)(1:LEN_TRIM(zflbls(i2)))//chrf
!      CALL plot_dumpd(r,nj,zbuf,fplotname)
    ENDDO

    DO j = 1, nsfast
      CALL readblock('BDENSSW',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                         i4,skip_summary,swap_lbl,zbuf,value_intg)
!      IF(i2 .LT. 0)CALL STOP('load_12 bdenssw',12)  
      IF(i2 .LT. 0)CYCLE                     !HSJ 04/08/2004 bdenssw is not 
                                             !always printed out by nubeam ??
      WRITE(intfile,FMT='(i1)')i2
      READ(intfile,FMT='(A)')chrf
      chrf = ADJUSTL(chrf)
!      fplotname ='bdensw' //zflbls(i2)(1:LEN_TRIM(zflbls(i2)))//chrf
!      CALL plot_dumpd(r,nj,zbuf,fplotname)
    ENDDO


    wbeam_nub(:) =0.0
    DO j = 1, nsfast
      CALL readblock('UDENSPL',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,    &
                      cvrt,i1,i2,i3,i4,skip_summary,swap_lbl,zbuf,value_intg) 
      IF(i2 .GT. 0)THEN
        DO i1 =1,nj
           wbeam_nub(i1) = wbeam_nub(i1) +zbuf(i1)* 0.62415064e16      !kev/cm**3
        ENDDO                                            !assumes cvrt = 1.e-6
      ENDIF
      IF(i2 .LT. 0)CYCLE            
      WRITE(intfile,FMT='(i1)')i2
      READ(intfile,FMT='(A)')chrf
      chrf = ADJUSTL(chrf)
!      fplotname ='udenspl' //zflbls(i2)(1:LEN_TRIM(zflbls(i2)))//chrf
!      CALL plot_dumpd(r,nj,zbuf,fplotname)
    ENDDO



    DO j = 1, nsfast
      CALL readblock('UDENSPP',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                           i4,skip_summary,swap_lbl,zbuf,value_intg)
      IF(i2 .GT. 0)THEN
        DO i1 =1,nj
           wbeam_nub(i1) = wbeam_nub(i1) +zbuf(i1)* 0.62415064e16      !kev/cm**3
        ENDDO                                            !assumes cvrt = 1.e-6
      ELSE
        CYCLE
      ENDIF           
      WRITE(intfile,FMT='(i1)')i2
      READ(intfile,FMT='(A)')chrf
      chrf = ADJUSTL(chrf)
!      fplotname ='udenspp' //zflbls(i2)(1:LEN_TRIM(zflbls(i2)))//chrf
!      CALL plot_dumpd(r,nj,zbuf,fplotname)
    ENDDO



    DO j = 1, nsfast
      CALL readblock('UDENSPLW',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                                       i4,skip_summary,swap_lbl,zbuf,value_intg)
!      IF(i2 .LT. 0)CALL STOP('load_12 UDENSPLW ',12)   
      IF(i2 .LT. 0)CYCLE            !UDENSPLW  may be missing sometimes      
      WRITE(intfile,FMT='(i1)')i2
      READ(intfile,FMT='(A)')chrf
      chrf = ADJUSTL(chrf)
!      fplotname ='udensplw' //zflbls(i2)(1:LEN_TRIM(zflbls(i2)))//chrf
!      CALL plot_dumpd(r,nj,zbuf,fplotname)
    ENDDO


    DO j = 1, nsfast
      CALL readblock('UDENSPPW',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                               i4,skip_summary,swap_lbl,zbuf,value_intg)
!      IF(i2 .LT. 0)CALL STOP('load_12 UDENSPPW',12) 
      IF(i2 .LT. 0)CYCLE            !UDENSPPW  may be missing sometimes              
      WRITE(intfile,FMT='(i1)')i2
      READ(intfile,FMT='(A)')chrf
      chrf = ADJUSTL(chrf)
!      fplotname ='udensppw' //zflbls(i2)(1:LEN_TRIM(zflbls(i2)))//chrf
!      CALL plot_dumpd(r,nj,zbuf,fplotname)
    ENDDO

    DO j = 1, 3   !PROBLEM HERE WITH NEG ION INJECTION
                  !BDEPE is not present if the beam is off.
      CALL readblock('BDEPE',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3,  &
                                        i4,skip_summary,swap_lbl,zbuf,value_intg)
!      IF(i2 .LT. 0)CALL STOP('load_12 bdepe',13)  
       IF(i2 .LT. 0)CYCLE           
      WRITE(intfile,FMT='(i1)')i2
      READ(intfile,FMT='(A)')chrf
      chrf = ADJUSTL(chrf)
!      fplotname ='bdepe' //'_E'//chrf
!      CALL plot_dumpd(r,nj,zbuf,fplotname)
    ENDDO

    DO j = 1, nsfast
      CALL readblock('BDEPS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3,  &
                                           i4,skip_summary,swap_lbl,zbuf,value_intg)
!      IF(i2 .LT. 0)CALL STOP('load_12 bdeps',14)     
       IF(i2 .LT. 0)CYCLE        
      WRITE(intfile,FMT='(i1)')i2
      READ(intfile,FMT='(A)')chrf
      chrf = ADJUSTL(chrf)
      fplotname ='bdeps' //zflbls(i2)(1:LEN_TRIM(zflbls(i2)))//chrf
      CALL plot_dumpd(r,nj,zbuf,fplotname)
    ENDDO

!   the following are not always present in the nubeam output file ??
!   sbbidpe(mj, 3, mibs)  ! n/sec (integrated profile)
        ! BEAM *DEPOSITION* BY BEAM-BEAM II, BY E.FRAC. & by species
!   sbbidps(mj, mibs)  ! n/sec (integrated profile)
        ! BEAM *DEPOSITION* BY BEAM-BEAM II, BY SPECIES
        ! summed over energy fractions
!
!   sbbxdpe(mj, 3, mibs)  ! n/sec (integrated profile)
        ! BEAM *DEPO* BY BEAM-BEAM CX, BY E.FRAC. & by species
 
!   sbbxdps(mj, mibs)  ! n/sec (integrated profile)
        ! BEAM *DEPOSITION* BY BEAM-BEAM CX, BY SPECIES
        ! summed over energy fractions
  
    TYPE ='V'
    cvrt = 1.e-6
    DO k=1,nsfast                    !loop over fast species
         DO i=1,3                    !beam energy componenets
            CALL readblock('SBBIDPE',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

         ENDDO 
    ENDDO

    DO k=1,nsfast                    !loop over fast species
            CALL readblock('SBBIDPS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO

    DO k=1,nsfast                    !loop over fast species
         DO i=1,3                    !beam energy componenets
            CALL readblock('SBBXDPE',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
!            IF(i2 .EQ. -1)THEN
!               PRINT *,'sbbxdpe not present'
!            ELSE
!               PRINT *,'sbbxdpe,i1 ,i2 =',i1,i2
!               PRINT *,'value_intg =',value_intg
!            ENDIF
         ENDDO 
    ENDDO

    DO k=1,nsfast                    !loop over fast species
            CALL readblock('SBBXDPS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO

! read   sbcxdpe(fast species,beam energy,thermal species ) 
      DO k=1,nsfast                  !loop over fast species
         DO i=1,3                    !beam energy componenets
           DO j = 1, ngmax           !loop over thermal species
             CALL readblock('SBCXDPE',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
         ENDDO
    ENDDO


! read   sbcxdps(fast species,thermal species ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, ngmax           !loop over thermal species
             CALL readblock('SBCXDPS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


! read   sbedepe(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('SBEDEPE',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


! read   sbedeps(fast species) 
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SBEDEPS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
    ENDDO


! read   sbiedepe(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('SBIEDEPE',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO

! read   sbiedeps(fast species) 
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SBIEDEPS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO

! read   sbildepe(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('SBILDEPE',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


! read   sbildeps(fast species) 
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SBILDEPS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
    ENDDO


! read   sbizdepe(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('SBIZDEPE',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


! read   sbizdeps(fast species) 
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SBIZDEPS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO

! bnb_ien,pbe_ien,pbi_ien,upb_ien,fpb_ien,edp_ien,ddp_ien may not be present
! these are special outputs for NCLASS
! read   bnb_ien(fast species,beam energy ) 
      DO k=1,nbeams                  !loop over # beams
           DO j = 1, 3               !loop over beam energies
             CALL readblock('BNB_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO

! read   pbe_ien(#beams,beam energy ) 

      DO k=1,nbeams                  !loop over # beams
           DO j = 1, 3               !loop over beam energies
             CALL readblock('PBE_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
           ENDDO
    ENDDO

! read   pbi_ien(#beams,beam energy ) 
      DO k=1,nbeams                  !loop over # beams
           DO j = 1, 3               !loop over beam energies
             CALL readblock('PBI_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO

! read   upb_ien(#beams,beam energy ) 
      DO k=1,nbeams                  !loop over # beams
           DO j = 1, 3               !loop over beam energies
             CALL readblock('UPB_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


! read   fpb_ien(#beams,beam energy ) 
      DO k=1,nbeams                  !loop over # beams
           DO j = 1, 3               !loop over beam energies
             CALL readblock('FPB_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


! read   edp_ien(#beams,beam energy ) 
      DO k=1,nbeams                  !loop over # beams
           DO j = 1, 3               !loop over beam energies
             CALL readblock('EDP_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
           ENDDO
    ENDDO


! read   ddp_ien(#beams ,beam energy ) 
      DO k=1,nbeams                  !loop over # beams
           DO j = 1, 3               !loop over beam energies
             CALL readblock('DDP_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
           ENDDO
    ENDDO


! read   scalars ein_ien 
    DO k=1,nbeams
           DO j = 1, 3               !loop over beam energies
             CALL read_scalar2('EIN_IEN',line,io,i1,i4,zbuf(1))
!             PRINT *,'EIN_EIN ,i1,i4 = ',i1,i4,zbuf(1)
           ENDDO
    ENDDO

! read   scalars dei_ien 
    DO k=1,nbeams
           DO j = 1, 3               !loop over beam energies
             CALL read_scalar2('DEI_IEN',line,io,i1,i4,zbuf(1))
!             PRINT *,'DEI_IEN ,i1,i4 = ',i1,i4,zbuf(1)
           ENDDO
    ENDDO

! read   bnbs_ien(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('BNBS_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO

! read   pbes_ien(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('PBES_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


! read   pbis_ien(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('PBIS_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


     skip_summary = .TRUE.
! read   upbs_ien(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('UPBS_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
           ENDDO
    ENDDO


! read   fpbs_ien(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('FPBS_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


! read   edps_ien(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('EDPS_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO
 

! read   ddps_ien(fast species,beam energy ) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('DDPS_IEN',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO
    skip_summary = .FALSE.





! read   bn00t(3,mj)  ! n/sec (integrated profile)
        ! TOTAL 1ST GEN BEAM NEUTRAL DENSITIES
        ! full/half/third energy, summed over all species

        DO j = 1, 3               !loop over beam energies
             CALL readblock('BN00T',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

        ENDDO

! read   bn00ts  (E1,D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, 3               !loop over beam energies
             CALL readblock('BN00TS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO
    

 
! read   bn00(# beams,beam energy ) (beam_001,E1)
      DO k=1,nbeams                  !loop over beams
           DO j = 1, 3               !loop over beam energies
             CALL readblock('BN00',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO

! read   bn0tot(fast species) (D_MCBEAM) 
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('BN0TOT',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO
    

! read   bn0cx(fast species) 
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('BN0CX',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO


! read   btrap (fast species) 
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('BTRAP',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,.TRUE.,swap_lbl,zbuf,value_intg)

    ENDDO


! read   btrap0 (fast species) 
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('BTRAP0',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,.TRUE.,swap_lbl,zbuf,value_intg)

    ENDDO

! read   scalars tbtrap
    DO k=1,nsfast
             CALL read_scalar1('TBTRAP',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,TBTRAP = ',i2,zbuf(1)
    ENDDO

! read   scalars tbtrap0
    DO k=1,nsfast
             CALL read_scalar1('TBTRAP0',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,TBTRAP0 = ',i2,zbuf(1)
    ENDDO


! btrape,desnfi not always present":
! read   btrape  fast species,N_miefi
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, n_miefi       !loop over ??
             CALL readblock('BTRAPE',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,.TRUE.,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


! read   densfi  fast species,N_miefi
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, n_miefi       !loop over ??
             CALL readblock('DENSFI',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,.TRUE.,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO

  skip_summary =  .TRUE.
! read ebapl2s(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('EBAPL2S',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
  
    ENDDO


! read ebapp2(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('EBAPP2S',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO

! read ebapl2sw(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('EBAPL2SW',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
    ENDDO

! read ebapp2sw(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('EBAPP2SW',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO

! read ebavps(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('EBAVPS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
    ENDDO


! read excit0(E1,D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
         DO j = 1, 3                 !loop over beam energies
             CALL readblock('EXCIT0',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
         ENDDO
    ENDDO

! excitr may not be present
! read excitr(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('EXCITR',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
    ENDDO

  skip_summary = .FALSE.
! read nmc(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('NMCFS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
!          print *,'nmcvs value intg =',value_intg
    ENDDO


! read wnmc(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('WNMC',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO

! read   scalars nmctot
    DO k=1,nsfast
             CALL read_scalar1('NMCTOT',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,nmctot = ',i2,zbuf(1)
    ENDDO

! read   scalars wnmctot
    DO k=1,nsfast
             CALL read_scalar1('WNMCTOT',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,WNMCTOT = ',i2,zbuf(1)
    ENDDO


! read   scalars nmcloss
    DO k=1,nsfast
             CALL read_scalar1('NMCLOSS',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,nmcloss = ',i2,zbuf(1)
    ENDDO

 

  skip_summary = .TRUE.
  ! read omegabs(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('OMEGABS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO

! read omegabsw(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('OMEGABSW',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO


  skip_summary = .FALSE.
! read pbcprs(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PBCPRS',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO
! read pbcprss(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PBCPRSS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO


! read pbdfbes(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PBDFBES',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO


! read pbdfbis(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PBDFBIS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

    ENDDO




! read   pbe(mj)  ! watts/m**3
! POWER TO ELECTRONS FROM BEAM HEATING -- for electron power balance:
    qbeame_intg_nub =0.0
    cvrt = 6.2415097e+9   ! watts/m**3 to kev/(cm**3sec)
    CALL readblock('PBE',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt, &
                 i1,i2,i3,i4,.FALSE.,swap_lbl,qbeame_nub,qbeame_intg_nub)
        IF(i2 .EQ. -1)qbeame_nub(:) = 0.0d0
        WHERE(qbeame_nub < 0.0)qbeame_nub =0.0


!pbes may ot be present
! read pbes(D_MCBEAM) watts/m**3
        ! POWER TO ELECTRONS FROM BEAM HEATING, BY FAST SPECIES
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PBES',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
      ENDDO




! read  pbi(mj)  ! watts/M**3, cvrt converts
        ! POWER TO IONS FROM BEAM SLOWING DOWN -- for ion power balance
       qbeami_intg_nub = 0.0
       cvrt = 6.2415097e+9   ! watts/m**3 to kev/(cm**3sec)
       CALL readblock('PBI',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                  i4,.FALSE.,swap_lbl,qbeami_nub,qbeami_intg_nub)
        IF(i2 .EQ. -1)qbeami_nub(:) = 0.0d0
        WHERE(qbeami_nub < 0.0)qbeami_nub =0.0




!pbis may ot be present
! read pbis(mj, mibs)  ! watts/m**3
! POWER TO IONS FROM BEAM HEATING, BY SPECIES
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PBIS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
!           print *,'pbis value intg =',value_intg
      ENDDO


! read pbohs(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PBOHS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! pbth ,pbths,pbthr2s,pfe,pfi may not be present
!note pbth is to be added to pbe to get total beam thermalization power:
! read pbth :
             qbth_nub(:) =0.0
             cvrt = 6.2415097e+9   ! watts/m**3 to kev/(cm**3sec)
             CALL readblock('PBTH',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,qbth_nub,qbth_intg_nub)

            WHERE(qbth_nub < 0.0)qbth_nub =0.0
            qbeami_nub(1:nj) = qbeami_nub(1:nj)  + qbth_nub(1:nj) 
                                            ! thermal energy of fast ions 
                                            ! that drop below thermal Ti 
                                            ! is added to qbeami

! read pbths(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PBTHS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO



! read pbthr2s(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PBTHR2S',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! read pfe:
             CALL readblock('PFE',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


! read pfi:
             CALL readblock('PFI',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)






! read curjbeam(mj)  ! beam ion driven current (shielded), amps/M**2
            TYPE ='A'
            cvrt = 1.e-4
             CALL readblock('CURJBEAM',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
! 
! read curjfusi:
            TYPE ='A'
             CALL readblock('CURJFUSI',.FALSE.,line,io,nj,'A',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


! read curjs(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('CURJS',.TRUE.,line,io,nj,'A',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! read ucurjs(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('UCURJS',.TRUE.,line,io,nj,'A',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! read ucurjco:
            TYPE ='A'
             CALL readblock('UCURJCO',.FALSE.,line,io,nj,'A',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


! read ucurjctr:
            TYPE ='A'
             CALL readblock('UCURJCTR',.FALSE.,line,io,nj,'A',zflbls, &
                  zthlbls,cvrt,i1,i2,i3,i4,skip_summary,swap_lbl,     &
                                                     zbuf,value_intg)


! read curb_nub   output by Nubeam as  < Jdot B > (amp tesla)/m**2:
            TYPE ='A'
               
             cvrt = 1.d0/SIGN(btor,rbt)  !from A/m**2 to A/cm**2 and tesla to gauss
                    !cvrt converts curb (=< Jdot B/Bt0>) to amps/cm**2
             CALL readblock('CURB',.FALSE.,line,io,nj,'A',zflbls,zthlbls, &
                    cvrt,i1,i2,i3,i4,.TRUE.,swap_lbl,curb_nub,curb_intg_nub)
!           curb_intg_nub is not put out by nubeam
            IF(i2 .EQ. -1)curb_nub(:) = 0.0
            curbe(:) = 0.0       !electron drag not calcualted separately
            curbi(:) =0.0        !ion
!            print *,'curb_intg , Amps = ',curb_intg_nub


!xjbfac   ! (dimensionless)
        ! --> fast ion CURRENT SHIELDING FACTOR (btw 0 and 1)
        ! (1-1/ZEFF FOR SPITZER, + NEOCLASSICAL TERMS as requested)
        ! cf NMCURB -- member of "miscellaneous" input control structure
            TYPE ='A'
            cvrt = 1.0
       CALL readblock('XJBFAC',.FALSE.,line,io,nj,'A',zflbls,zthlbls,     &
                       cvrt,i1,i2,i3,i4,.TRUE.,swap_lbl,zbuf,value_intg)
 







! read sbcx0s(D_MCBEAM)
      skip_summary =.FALSE.
      
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SBCX0S',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO



!
      TYPE ='V'
! read sorbn0(D_th)
      DO k=1,ngmax                 !loop thermal fast species
             CALL readblock('SORBN0',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
      
      if(i3.eq.2) sorbn0_nub(:) = 1.0e-6*zbuf(:) !jmp.den
      if(i3.eq.2) print *,sorbn0_nub
      !if(i3.eq.2) pause
      !print *,i3 !jmp.den
      !print *,zbuf !jmp.den
      !pause !jmp.den

      ENDDO


! read  sorbn0s  fast species,thermals species (H_th,D_MCBEAM) 
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, ngmax           !loop over thermal species
             CALL readblock('SORBN0S',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO


! read sbcxrbs(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SBCXRBS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
      ENDDO


! read sbcxres(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SBCXRES',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! read sbcxris(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SBCXRIS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
      ENDDO

! read sbcxrls(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SBCXRLS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
      ENDDO


! read  sbcxrxs  fast species,thermals species
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, ngmax           !loop over thermal species
             CALL readblock('SBCXRXS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO
   
! read sbcxrzs(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SBCXRZS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

    skip_summary = .TRUE.
! read sinb0i(D_th)
      DO k=1,ngmax                 !loop thermal species
             CALL readblock('SINB0I(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO


! read  sinb0is fast species,thermals species (D_th,D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, ngmax           !loop over thermal species
             CALL readblock('SINB0IS(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 
           ENDDO
    ENDDO

! read sinbn0(D_MCBEAM)
      DO k=1,ngmax                 !loop thermal species
             CALL readblock('SINBN0(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO


! read  sinbn0s fast species,thermals species D_th,D_MCBEAM
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, ngmax           !loop over thermal species
             CALL readblock('SINBN0S(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           ENDDO
    ENDDO

   skip_summary = .FALSE.

! read sorbe:
            TYPE ='V' !jmp.den
             CALL readblock('SORBE',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

           PRINT *,'sorbe,value intg =',value_intg
           !print *,zbuf !jmp.den
           !pause !jmp.den
! read sorbes(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SORBES(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! read sorbh(D_MCBEAM)
      TYPE ='V' !jmp.den
      DO k=1,ngmax                 !loop thermal species
             CALL readblock('SORBH(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
        PRINT *,'SORBH, value_intg =',value_intg

      if(i3.eq.2) sorbh_nub(:) = 1.0e-6*zbuf(:) !jmp.den
      !print *,i3 !jmp.den
      !print *,zbuf !jmp.den
      !pause !jmp.den

      ENDDO
! read  sorbhs fast species,thermals species
      DO k=1,nsfast                  !loop over fast species
           DO j = 1, ngmax           !loop over thermal species
             CALL readblock('SORBHS(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
             PRINT *,'SORBHS, value_intg ngmax,nsfast =',value_intg,ngmax,nsfast
           ENDDO
    ENDDO

! read sorbth:
            TYPE ='V'
             CALL readblock('SORBTH',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

       IF (i2 .EQ. -1)PRINT *, 'skipping reading of sorbth'
       !print *,zbuf !jmp.den
       !pause

! read sorbths(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('SORBTHS(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT *, 'skipping reading of sorbths'
      ENDDO

      skip_summary = .TRUE.
! read taubpa0(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBPA0(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! read taubpaco1(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBPACO1(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! read taubpaco2(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBPACO2(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! read taubpaco3(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBPACO3(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! read taubpactr1(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBPACTR1(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT *, 'skipping reading of taubpactr1'
      ENDDO

! read taubpactr2(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBPACTR2(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT *, 'skipping reading of taubpactr2'
      ENDDO

! read taubpactr3(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBPACTR3(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of taubpactr3'
      ENDDO

! read taubsl0(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBSL0(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of taubsl0'
      ENDDO

! read taubslco1(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBSLCO1',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of taubslco1'
      ENDDO

! read taubslco2(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBSLCO2(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of taubslco2'
      ENDDO

! read taubslco3(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBSLCO3(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of taubslco3'
      ENDDO

! read taubslctr1(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBSLCTR1',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of taubslctr1'
      ENDDO

! read taubslctr2(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBSLCTR2(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
        IF (i2 .EQ. -1)PRINT*, 'skipping reading of taubslctr2'
      ENDDO

! read taubslctr3(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TAUBSLCTR3(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT *, 'skipping reading of taubslctr3'
      ENDDO


! read paxfrc(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PAXFRC',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT *, 'skipping reading of paxfrc'
      ENDDO

! read pbcfrc(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('PBXFRC',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT *, 'skipping reading of pbxfrc'
      ENDDO

! read   tirabd          
             CALL readblock('TIRABD',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
        IF (i2 .EQ. -1)PRINT *, 'skipping reading of tirabd'
        PRINT *,'skip_summary =',skip_summary
! read tirabds(D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species
             CALL readblock('TIRABDS(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT *, 'skipping reading of tirabds'
        PRINT *,'skip_summary,in tirabds =',skip_summary
      ENDDO

! read   tirfshs         
             CALL readblock('TIRFSHS',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

       IF (i2 .EQ. -1)PRINT *, 'skipping reading of tirfshs'
! read tirorb(D_MCBEAM)

             CALL readblock('TIRORB',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of tirorb'
! read   tirorbs     (D_MCBEAM) 
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TIRORBS(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of tirorbs'
      ENDDO

  skip_summary =.FALSE.
! tirrip,tirorbs may ot be present
! read tirrip(D_MCBEAM)

             CALL readblock('TIRRIP',.FALSE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

       IF (i2 .EQ. -1)PRINT*, 'skipping reading of tirrip'
! read   tirrips     (D_MCBEAM) 
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TIRRIPS(',.TRUE.,line,io,nj,TYPE,zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

       IF (i2 .EQ. -1)PRINT*, 'skipping reading of tirrips'
      ENDDO

! read tqbsum A net total torque 
!   tqbsum = collisional torque + JxB torque + thermalization torque
!            + anom. diffusivity transfer - charge exchange loss
!          = tqbe + tqbi + tqjxb + tqbth + tqdifb - tqbcx
! for plasma momentum balance: the recommended momentum source
! is   tqbsum + c1*tqrpl - c2*tqohb   where c1 and c2 are
! adjustable parameters (e.g. TRANSP usually assumes c1=c2=0.0)
!      storqueb_intg_nub = 0.0
       cvrt =10. ! convert from NT-M/M**3 to dyne-cm/cm**3
       CALL readblock('TQBSUM',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt, &
                 i1,i2,i3,i4,.FALSE.,swap_lbl,storqueb_nub,storqueb_intg_nub)
            IF(i2 .EQ. -1)storqueb_nub(:) = 0.0
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of tqbsum'
        IF(i2 .GT. -1)PRINT *,'tqbsum intg_value =',storqueb_intg_nub

! average the results if called for:
        IF (avg_nubeam_torque .GT. 0)THEN
           IF(.NOT. ALLOCATED(storque_nub_avg)) THEN                           
                ALLOCATE(storque_nub_avg(avg_nubeam_torque,nj))
                storque_nub_avg(:,:) = 0.0
                torque_ct =0
           ENDIF

           torque_ct = torque_ct +1
           nterms = torque_ct
           if(load_calls .ge. avg_nubeam_torque) nterms = avg_nubeam_torque
           IF(torque_ct .GT. avg_nubeam_torque)torque_ct =1
           storque_nub_avg(torque_ct,1:nj) = storqueb_nub(1:nj)
           DO j =1,nj
               sum = 0.0
               DO i2 =1,nterms
                  sum =sum +storque_nub_avg(i2,j)
               ENDDO
               storqueb_nub(j) = sum/nterms
           ENDDO

        ENDIF



! read tqbcx

       CALL readblock('TQBCX',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt, &
                         i1,i2,i3,i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of tqbcx'

! read tqbcxs (D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TQBCXS',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
       IF (i2 .EQ. -1)PRINT*, 'skipping reading of tqbcxs'

      ENDDO

! read  tqbe(mj)  ! Nt-m/m**3 (integrated profile)
        ! TORQUE TO ELECTRONS FROM BEAM SLOWING DOWN
        sprbeame_intg_nub =0.0
        cvrt =0.1 ! convert from NT-M/M**3 to dyne-cm/cm**3
        CALL readblock('TQBE',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt, &
          i1,i2,i3,i4,skip_summary,swap_lbl,sprbeame_nub,sprbeame_intg_nub)
        IF(i2 .EQ. -1)sprbeame_nub(:) =0.0


! read tqbes (D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TQBES(',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


      ENDDO

! read   tqbi(mj)  ! Nt-m/M**3(integrated profile)
        ! TORQUE TO IONS FROM BEAM SLOWING DOWN
        sprbeami_intg_nub =0.0
        cvrt =0.1 ! convert from NT-M/M**3 to dyne-cm/cm**3
        CALL readblock('TQBI',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt, &
          i1,i2,i3, i4,skip_summary,swap_lbl,sprbeami_nub,sprbeami_intg_nub)
            IF(i2 .EQ. -1)sprbeami_nub(:) =0.0


! read tqbis (D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TQBIS(',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


      ENDDO

! read tqbth

             CALL readblock('TQBTH',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


! read tqbths (D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TQBTHS(',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


      ENDDO
   skip_summary =.TRUE.
! read tqbxf

             CALL readblock('TQBXF',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


! read tqbxfs (D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TQBXFS(',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


      ENDDO


  skip_summary =.FALSE.
! read tqdifb
             CALL readblock('TQDIFB',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
    

! read tqbdifbs (D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TQDIFBS(',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO


  skip_summary =.FALSE.
! read tqjxb
             CALL readblock('TQJXB',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)
 

! read tqbjxbs  (D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TQJXBS(',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


      ENDDO

  skip_summary =.FALSE.
! read tqohb
             CALL readblock('TQOHB',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


! read tqohbs  (D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TQOHBS(',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

         PRINT *,'TQOHBS value_integrated =',value_intg
      ENDDO


  skip_summary =.FALSE.
! read tqrpl tqrpls may ot be present
             CALL readblock('TQRPL',.FALSE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)


! read tqrpls  (D_MCBEAM)
      DO k=1,nsfast                  !loop over fast species    
             CALL readblock('TQRPLS(',.TRUE.,line,io,nj,'V',zflbls,zthlbls,cvrt,i1,i2,i3, &
                            i4,skip_summary,swap_lbl,zbuf,value_intg)

      ENDDO

! read   scalars bpcxers
    DO k=1,nsfast
             CALL read_scalar1('BPCXERS(',line,io,i2,zflbls,zbuf(1))
             PRINT *,'i2,BPCXERS = ',i2,zbuf(1)
    ENDDO

! read   scalars bpcxi0s
    DO k=1,nsfast
             CALL read_scalar1('BPCXI0S(',line,io,i2,zflbls,zbuf(1))
             PRINT *,'i2,BPCXI0S= ',i2,zbuf(1)
    ENDDO

! read   scalars bpcxins
    DO k=1,nsfast
             CALL read_scalar1('BPCXINS(',line,io,i2,zflbls,zbuf(1))
             PRINT *,'i2,BPCXINS= ',i2,zbuf(1)
    ENDDO

! read   scalars bpcxots
    DO k=1,nsfast
             CALL read_scalar1('BPCXOTS(',line,io,i2,zflbls,zbuf(1))
             PRINT *,'i2,BPCXOTSS= ',i2,zbuf(1)
    ENDDO

! read   scalars bpcxris
    DO k=1,nsfast
             CALL read_scalar1('BPCXRIS(',line,io,i2,zflbls,zbuf(1))
             PRINT *,'i2,BPCXRIS= ',i2,zbuf(1)
    ENDDO

! read   scalars bpcxrxs
    DO k=1,nsfast
             CALL read_scalar1('BPCXRXS(',line,io,i2,zflbls,zbuf(1))
             PRINT *,'i2,BPCXRXS= ',i2,zbuf(1)
    ENDDO

! read   scalars bpcxx0s
    DO k=1,nsfast
             CALL read_scalar1('BPCXX0S(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPCXX0S= ',i2,zbuf(1)
    ENDDO

! read bperr
   CALL read_scalar('BPERR',line,io,zbuf(1))
     PRINT *,'BPERR =',zbuf(1)


! read   scalars bperrs
    DO k=1,nsfast
             CALL read_scalar1('BPERRS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPERRS= ',i2,zbuf(1)
    ENDDO

! read   scalars bpfshs
    DO k=1,nsfast
             CALL read_scalar1('BPFSHS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPFSHS= ',i2,zbuf(1)
    ENDDO


! read   scalars bpls
    DO k=1,nsfast
             CALL read_scalar1('BPLS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPLS= ',i2,zbuf(1)
    ENDDO

! read   scalars bpohs
    DO k=1,nsfast
             CALL read_scalar1('BPOHS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPOHS= ',i2,zbuf(1)
    ENDDO

! read   scalars bprpls
    DO k=1,nsfast
             CALL read_scalar1('BPRPLS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPRPLS= ',i2,zbuf(1)
    ENDDO


! read   scalars bpshins
    DO k=1,nsfast
             CALL read_scalar1('BPSHINS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPSHINS= ',i2,zbuf(1)
    ENDDO


! read   scalars bpsts
    DO k=1,nsfast
             CALL read_scalar1('BPSTS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPSTS= ',i2,zbuf(1)
    ENDDO

! read   scalars bpths
    DO k=1,nsfast
             CALL read_scalar1('BPTHS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPTHS= ',i2,zbuf(1)
    ENDDO


! read pftota
   CALL read_scalar('PFTOTA',line,io,zbuf(1))
!     print *,'PFTOTA =',zbuf(1)

! read   scalars bpprma
    DO k=1,nbeams
             CALL read_scalar1('BPPRMA(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPPRMA= ',i2,zbuf(1)
    ENDDO

! read   scalars pinjs
    pwf_tot_intg_nub =0.0D0
! injected power by beams is in beam_data%pinja(1..beam_data%nbeam)
! If beams are different species then pinjs(beam species) gives the
! power injected by that species.
    DO k=1,nsfast
             CALL read_scalar1('PINJS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,PINJS= ',i2,zbuf(1)
        IF(i2 .GT. 0)THEN
            pwf_tot_intg_nub =  pwf_tot_intg_nub + ZBUF(1)
        ENDIF
    ENDDO



! read   scalars bph0was
    DO k=1,nsfast
             CALL read_scalar1('BPH0WAS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPH0WAS= ',i2,zbuf(1)
    ENDDO


! read   scalars bphcks
    DO k=1,nsfast
             CALL read_scalar1('BPHCKS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHCKS= ',i2,zbuf(1)
    ENDDO

! read   scalars bphcols
    DO k=1,nsfast
             CALL read_scalar1('BPHCOLS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHCOLS= ',i2,zbuf(1)
    ENDDO


! read   scalars bphcons
    DO k=1,nsfast
             CALL read_scalar1('BPHCONS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHCONS= ',i2,zbuf(1)
    ENDDO



! read   scalars bphcxls
    DO k=1,nsfast
             CALL read_scalar1('BPHCXLS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHCXLS= ',i2,zbuf(1)
    ENDDO

! read   scalars bphdeps
    DO k=1,nsfast
             CALL read_scalar1('BPHDEPS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHDEPS= ',i2,zbuf(1)
    ENDDO


! read   scalars bphdfbs
    DO k=1,nsfast
             CALL read_scalar1('BPHDFBS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHDFBS= ',i2,zbuf(1)
    ENDDO

! read bpherr
   CALL read_scalar('BPHERR',line,io,zbuf(1))
!     print *,'BPHERR =',zbuf(1)


! read   scalars bpherrs
    DO k=1,nsfast
             CALL read_scalar1('BPHERRS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHERRS= ',i2,zbuf(1)
    ENDDO


! read   scalars bphfshs
    DO k=1,nsfast
             CALL read_scalar1('BPHFSHS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHFSHS= ',i2,zbuf(1)
    ENDDO



! read   scalars bphjxbs
    DO k=1,nsfast
             CALL read_scalar1('BPHJXBS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHJXBS= ',i2,zbuf(1)
    ENDDO

! read   scalars bphohbs
    DO k=1,nsfast
             CALL read_scalar1('BPHOHBS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHOHBS= ',i2,zbuf(1)
    ENDDO



! read   scalars bpholds
    DO k=1,nsfast
             CALL read_scalar1('BPHOLDS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHOLDS= ',i2,zbuf(1)
    ENDDO


! read   scalars bphorbs
    DO k=1,nsfast
             CALL read_scalar1('BPHORBS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHORBS= ',i2,zbuf(1)
    ENDDO



! read   scalars bphrcas
    DO k=1,nsfast
             CALL read_scalar1('BPHRCAS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHRCAS= ',i2,zbuf(1)
    ENDDO



! read   scalars bphrpls
    DO k=1,nsfast
             CALL read_scalar1('BPHRPLS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHRPLS= ',i2,zbuf(1)
    ENDDO



! read   scalars bphsts
    DO k=1,nsfast
             CALL read_scalar1('BPHSTS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHSTS= ',i2,zbuf(1)
    ENDDO



! read   scalars bphths
    DO k=1,nsfast
             CALL read_scalar1('BPHTHS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPHTHS= ',i2,zbuf(1)
    ENDDO



! read   scalars bptcols
    DO k=1,nsfast
             CALL read_scalar1('BPTCOLS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPTCOLS= ',i2,zbuf(1)
    ENDDO


! read   scalars bptdfbs
    DO k=1,nsfast
             CALL read_scalar1('BPTDFBS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPTDFBS= ',i2,zbuf(1)
    ENDDO



! read   scalars bpthros
    DO k=1,nsfast
             CALL read_scalar1('BPTHROS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPTHROS= ',i2,zbuf(1)
    ENDDO


! read   scalars bpthsfs
    DO k=1,nsfast
             CALL read_scalar1('BPTHSFS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPTHSFS= ',i2,zbuf(1)
    ENDDO



! read   scalars bpjxbs
    DO k=1,nsfast
             CALL read_scalar1('BPJXBS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPJXBS= ',i2,zbuf(1)
    ENDDO



! read   scalars bpthafs
    DO k=1,nsfast
             CALL read_scalar1('BPTHAFS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPTHAFS= ',i2,zbuf(1)
    ENDDO



! read   scalars bxflws
    DO k=1,nsfast
             CALL read_scalar1('BXFLWS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BXFLWS= ',i2,zbuf(1)
    ENDDO


! read   scalars bxfshs
    DO k=1,nsfast
             CALL read_scalar1('BXFSHS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BPXFSHS= ',i2,zbuf(1)
    ENDDO



! read   scalars bxrpls
    DO k=1,nsfast
             CALL read_scalar1('BXRPLS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BXEPLS= ',i2,zbuf(1)
    ENDDO


! read   scalars sbcxxs
    DO k=1,nsfast
             CALL read_scalar1('SBCXXS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SBCXXS= ',i2,zbuf(1)
    ENDDO



! read   scalars sfcx0mc
    DO k=1,nsfast
             CALL read_scalar1('SFCX0MC(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFCX0MC = ',i2,zbuf(1)
    ENDDO


! read   scalars sfcxbal
    DO k=1,nsfast
             CALL read_scalar1('SFCXBAL(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,BSFCXBAL= ',i2,zbuf(1)
    ENDDO


! read   scalars sfcxesc
    DO k=1,nsfast
             CALL read_scalar1('SFCXESC(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFCXESC = ',i2,zbuf(1)
    ENDDO



! read   scalars sfcxrmc
    DO k=1,nsfast
             CALL read_scalar1('SFCXRMC(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFCXRMC = ',i2,zbuf(1)
    ENDDO




! read   scalars sfxcrr
    DO k=1,nsfast
             CALL read_scalar1('SFCXRR(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFCXRR= ',i2,zbuf(1)
    ENDDO


! read   scalars sfdbbcx
    DO k=1,nsfast
             CALL read_scalar1('SFDBBCX(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFDBBCX= ',i2,zbuf(1)
    ENDDO



! read   scalars sfdBBiz
    DO k=1,nsfast
             CALL read_scalar1('SFDBBIZ(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFDBBIZ = ',i2,zbuf(1)
    ENDDO


! read   scalars sdepba
    DO k=1,nsfast
             CALL read_scalar1('SFDEPBA(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFDEPBA = ',i2,zbuf(1)
    ENDDO


! read   scalars sfddepcx
    DO k=1,nsfast
             CALL read_scalar1('SFDEPCX(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFDEPCX= ',i2,zbuf(1)
    ENDDO



! read   scalars sfdepiz
    DO k=1,nsfast
             CALL read_scalar1('SFDEPIZ(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFDEPIZ = ',i2,zbuf(1)
    ENDDO



! read   scalars sdepmc
    DO k=1,nsfast
             CALL read_scalar1('SFDEPMC(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFDEPMC = ',i2,zbuf(1)
    ENDDO


! read   scalars sfddeprr
    DO k=1,nsfast
             CALL read_scalar1('SFDEPRR(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFDEPRR = ',i2,zbuf(1)
    ENDDO



! read   scalars sfdepsc
    DO k=1,nsfast
             CALL read_scalar1('SFDEPSC(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFDEPSC ',i2,zbuf(1)
    ENDDO


! read   scalars sfdtbmc
    DO k=1,nsfast
             CALL read_scalar1('SFDTBMC(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFDTBMC = ',i2,zbuf(1)
    ENDDO



! read   scalars sforbal
    DO k=1,nsfast
             CALL read_scalar1('SFORBAL(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFORBAL ',i2,zbuf(1)
    ENDDO


! read   scalars sforbrr
    DO k=1,nsfast
             CALL read_scalar1('SFORBRR(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFORBRR ',i2,zbuf(1)
    ENDDO



! read   scalars sfrbbcx
    DO k=1,nsfast
             CALL read_scalar1('SFRBBCX(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFRBBCX',i2,zbuf(1)
    ENDDO




! read   scalars sfrbbiz
    DO k=1,nsfast
             CALL read_scalar1('SFRBBIZ(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFRBBIZ',i2,zbuf(1)
    ENDDO



! read   scalars sfshine
    DO k=1,nsfast
             CALL read_scalar1('SFSHINE(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFSHINE',i2,zbuf(1)
    ENDDO





! read   scalars sfxrccx
    DO k=1,nsfast
             CALL read_scalar1('SFXRCCX(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFXRCCX',i2,zbuf(1)
    ENDDO




! read   scalars sfxrciz
    DO k=1,nsfast
             CALL read_scalar1('SFXRCIZ(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFXRCIZ',i2,zbuf(1)
    ENDDO



! read   scalars sfxrcsc
    DO k=1,nsfast
             CALL read_scalar1('SFXRCSC(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFXRCSC',i2,zbuf(1)
    ENDDO




! read   scalars snbcxos
    DO k=1,nsfast
             CALL read_scalar1('SNBCXOS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFNBCXOS',i2,zbuf(1)
    ENDDO



! read   scalars snbcxs
    DO k=1,nsfast
             CALL read_scalar1('SNBCXS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SNBCXS',i2,zbuf(1)
    ENDDO



! read   scalars snbxbb0
    DO k=1,nsfast
             CALL read_scalar1('SNBXBB0(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SNBXBB0',i2,zbuf(1)
    ENDDO



! read   scalars snbxbb1
    DO k=1,nsfast
             CALL read_scalar1('SNBXBB1(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SNBXBB1',i2,zbuf(1)
    ENDDO




! read   scalars snbxtot
    DO k=1,nsfast
             CALL read_scalar1('SNBXTOT(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SNBXTOT',i2,zbuf(1)
    ENDDO


! read   scalars snbxv0
    DO k=1,nsfast
             CALL read_scalar1('SNBXV0(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SFNBXV0',i2,zbuf(1)
    ENDDO



! read   scalars snbxw0
    DO k=1,nsfast
             CALL read_scalar1('SNBXW0(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,SNBXW0',i2,zbuf(1)
    ENDDO


! read   scalars xysbcxs
    DO k=1,nsfast
             CALL read_scalar1('XYSBCXS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,XYSBCXS',i2,zbuf(1)
    ENDDO


! read   scalars xninjs
    DO k=1,nsfast
             CALL read_scalar1('XNINJS(',line,io,i2,zflbls,zbuf(1))
!             PRINT *,'i2,XNINJS',i2,zbuf(1)
    ENDDO


! read sftota
   CALL read_scalar('SFTOTA',line,io,zbuf(1))
!     print *,'SFTOTA =',zbuf(1)



    DEALLOCATE(zflbls)

  RETURN
  END
