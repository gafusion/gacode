   MODULE xsct
       USE nrtype,                                    ONLY : DP,I4B,I2B,SP

       USE nf_param,                                  ONLY : kcmp1,kbe,ksge,ke,kb,kprim

       USE ions_gcnmp,                                ONLY : nion,nprim,nimp,namep,namei,   &
                                                             atw,nimp_ml,namei_ml,          &
                                                             nimp_index,atomnoi_ml

       USE common_constants,                          ONLY : zeroc,izero,AMU_Value,Pi

       USE error_handler
       USE io_gcnmp,                                  ONLY : nlog,ncrt

       USE MPI_data,                                  ONLY : myid,master,numprocs,mpiierr

       USE tension_spline,                            ONLY : seval,spline_coef

       USE zonal_data,                                ONLY : nw,nh,mfm1


#if defined USEMPI
      USE mpi
#endif


       IMPLICIT NONE

       REAL(DP)      sgxnloc(kbe),                                &
                     sgxnmi(ke,kb),hxfrac(ke,kb)

       REAL(DP),    ALLOCATABLE,DIMENSION(:,:,:,:) ::  sgxn

       REAL(DP),    ALLOCATABLE,DIMENSION(:)       ::  zniim, znipm

       INTEGER(I4B),ALLOCATABLE,DIMENSION(:)       ::  izatom,iz,izstrp

       CHARACTER*256  adas_xsct_path

       CONTAINS

   SUBROUTINE setup_xsct  
!----------------------------------------------------------------------------
! -- use zonal quantities to get cross sections.
!-------------------------------------------------------------HSJ 1/20/2011--

    
     USE zonal_data,                ONLY :       wnoperm,zenbeam,zenbeamold,                 &
                                                 zone_volume,                                &
                                                 zone_area,potsid,psif,                      &
                                                 zone_area,potsid,rinsid,psivol,             &
                                                 rotsid,pinsid,b1ins,b2ins,                  &
                                                 b1ots,b2ots,fpsio,fpsii,zne,                &
                                                 zte,zti,zangrot,zni,zzi,mfm1

      USE  neutral_beams,          ONLY :        ebkev,fe_tk,ibion,no_injectors,             &
                                                 ne_tk,vbeam,de_tk,                          &
                                                 atw_beam,ncorin

                                                 
      IMPLICIT NONE
      INTEGER(I2B) get_next_io_unit

      IF(ALLOCATED(iz))DEALLOCATE(iz)          ; ALLOCATE(iz(nimp))
      IF(ALLOCATED(izatom))DEALLOCATE(izatom)  ; ALLOCATE(izatom(nion))
      IF(ALLOCATED(izstrp))DEALLOCATE(izstrp)  ; ALLOCATE(izstrp(nimp))
      IF(ALLOCATED(zniim)) DEALLOCATE(zniim)   ; ALLOCATE(zniim(nimp))
      IF(ALLOCATED(znipm)) DEALLOCATE(znipm)   ; ALLOCATE(znipm(nprim))
      IF(ALLOCATED(sgxn))  DEALLOCATE(sgxn)    ; ALLOCATE(sgxn(kcmp1,mfm1,kbe,ksge))


      IF(myid == master)THEN
              ncorin =  get_next_io_unit()
              CALL nbsgxn (ebkev,fe_tk,ibion,no_injectors,ne_tk,vbeam, &
                           de_tk,hxfrac,atw_beam)

      ENDIF

      IF(numprocs .GT. 1)CALL  distribute_xsct

      END SUBROUTINE setup_xsct  


       FUNCTION ceef (i, j, te)
!--------------------------------------------------------------------------
! -- e-impact excitation rate coefficient from Vriens & Smeets
!     used only for excitations to final state with n > 2 (i.e. j > 3).
!-------------------------------------------------------------------------
     USE hexnb_data,                             ONLY : ms,mc,         &  ! param
                                                        kdene,kdeni,   &
                                                        kdenz,ksvi,    &
                                                        ksvz,ksve,     &
                                                        krad,ngh,ngl,  &  ! b1
                                                        nouthx,istart, &
                                                        ihxbug,        &  ! b2
                                                        en,dg,ae,be,   &
                                                        de1,de2,ge1,   &
                                                        ge2               ! b4

      IMPLICIT NONE 
      REAL(DP) te,ryd,ceef,ceef1
      INTEGER(I4B) i,j,ni,nj,id

      DATA        ryd/13.6_DP/

  
      ceef1 (ni, nj, te) = 1.6e-7 * SQRT (te)                                                &
                  * EXP (-(ryd / te) * (1.0_DP / FLOAT (ni)**2 - 1.0_DP / FLOAT (nj)**2))    &
                  * (ae(ni,nj) * LOG (0.3_DP * te / ryd + de2(ni,nj)) + be(ni,nj))           &
                  / (te + ge2(ni,nj) * LOG (1.0_DP + (FLOAT (ni)**3) * te / ryd))
      IF (i .GE. j .OR. j .LE. 3) THEN
        ihxbug = 8
        IF (nouthx .GT. 0 .AND. myid == master) THEN
          WRITE (nouthx, 10)  i, j
   10     FORMAT (/ ' ERROR in CEEF: i = ', i3, '    j = ', i3 /)
        END IF
        ceef = 0.0
      ELSE
        ni   = nfhx(i)
        nj   = j - 1
        ceef = ceef1 (ni, nj, te)
        id   = 2
        IF (ni .EQ. 1)  id = 1
      END IF

      RETURN

      END FUNCTION ceef



      FUNCTION cief (i, te)
!---------------------------------------------------------------------
! --  e-impact ionization rate coefficient from Vriens & Smeets,
!     Phys. Rev. a 22, 940 (1980).
!     used only for 2p (i = 3) or for n > 2 (i=n+1).
!-----------------------------------------------------------------------
      USE hexnb_data,                            ONLY : ms,mc,         &  ! param
                                                        kdene,kdeni,   &
                                                        kdenz,ksvi,    &
                                                        ksvz,ksve,     &
                                                        krad,ngh,ngl,  &  ! b1
                                                        nouthx,istart, &
                                                        ihxbug,        &  ! b2
                                                        en,dg,ae,be,   &
                                                        de1,de2,ge1,   &
                                                        ge2               ! b4

      IMPLICIT NONE
      REAL(DP) cief,te,ent
      INTEGER(I4B) i

      ent(i) = en(i) / te ! ent is inline function

      IF (i .LE. 2) THEN
        ihxbug = 7
        IF (nouthx .GT. 0 .AND. myid == master) THEN
          WRITE  (nouthx, 10) i
   10     FORMAT (/ ' ERROR in CIEF: i = ', i3 /)
        END IF
        cief = zeroc
      ELSE
        cief = 9.56e-6_DP * EXP (-ent(i))                           &
                      / (te * SQRT (te) * (ent(i)**2.33          &
                      + 4.38_DP * ent(i)**1.72_DP + 1.32_DP * ent(i)))
      END IF

      RETURN

      END FUNCTION cief



      SUBROUTINE hxinit
!---------------------------------------------------------------------------------
! -- 
!---------------------------------------------------------------------------------
      USE cpub_dat
   
      USE hexnb_data,                                     ONLY : ms,mc,       & ! param
                                                                 kdene,kdeni, &
                                                                 kdenz,ksvi,  &
                                                                 ksvz,ksve,   &
                                                                 krad,ngh,    &
                                                                 ngl,         & ! b1
                                                                 f,ar,        & ! b3
                                                                 en,dg,ae,be, &
                                                                 de1,de2,ge1, &
                                                                 ge2,         & ! b4
                                                                 iexcit,      &
                                                                 ilorent,     &
                                                                 mstate,ncont   ! b8


      IMPLICIT NONE
!


      REAL(DP)    be1(mc),ryd,rrate,an1,en1,en12,s,ajsq, &
                  tot,frac2s,an2
      REAL(SP)    cpub,cpua
      INTEGER(I4B) i,j,n1,n2,nj,ni

      DATA          ryd/13.6_DP/
!
! --- total radiation rate from level n2 to level n1:
!
      rrate(n2,n1) = (8.0323e9_DP) * &
                   (((1.0+DP/FLOAT (n1))**2-(1.0_DP/FLOAT (n2))**2) * &
                         (FLOAT (n1)/FLOAT (n2)))**2*f(n1,n2)
!
      CALL SECOND (cpua)
!
! --- tabulate oscillator strengths
!
      CALL hxosc
!
      dg(1) = 1.0
      dg(2) = 1.0
      dg(3) = 3.0
      en(1) = ryd
      en(2) = ryd / 4.0
      en(3) = ryd / 4.0
!
      DO i=4,mc+1
        dg(i) = (FLOAT (i-1))**2
        en(i) =  ryd / dg(i)
      END DO
!
      DO n1=1,mc
        an1     = FLOAT (n1)
        be1(n1) = 1.4 * LOG (an1) / an1 &
                - 0.70 / an1 - 0.51 / (an1*an1) + 1.16 / (an1**3) &
                - 0.55 / (an1**4)
      END DO
!
      DO 9 n1=1,mstate
        an1  = FLOAT (n1)
        en1  = ryd / (an1 * an1)
        DO 9 n2=n1+1,mc
          an2        = FLOAT (n2)
          en12       = en1 - ryd / (an2*an2)
          ae(n1,n2)  = 2.0 * ryd * f(n1,n2) / en12
          be(n1,n2)  = 4.0 * ryd * ryd * (1.0 / (en12*en12) &
                     + 4.0 * en1 / (3.0 * en12**3) &
                     + be1(n1) * en1 * en1 / (en12**4)) / (an2**3)
          de1(n1,n2) = EXP (-be(n1,n2) / ae(n1,n2)) - 0.4 * en12 / ryd
          s          = an2 - an1
          de2(n1,n2) = EXP (-be(n1,n2) / ae(n1,n2)) &
                     + 0.06 * s * s / (an1 * an1 * an2)
          ge1(n1,n2) = ryd * (8.0 + 23.0 * (s/an1)**2) &
                     / (8.0 + 1.1 * an2 * s + 0.8 / (s * s) &
                     + 0.4 * (s - 1.0) * SQRT (an2 * an2 * an2 / s))
          ge2(n1,n2) = ryd * (3.0 + 11.0 * (s / an1)**2) &
                     / (6.0 + 1.6 * an2 * s + 0.3 / (s * s) &
                     + 0.8 * (s - 0.6) * SQRT (an2 * an2 * an2 / s))
    9 CONTINUE
!
! --- radiation rates
!
      DO   i=1,mstate+1
        DO j=1,mstate+1
          ar(i,j) = 0.0
        END DO
      END DO
!
      IF (krad .EQ. 0)  go to 1000
!
      DO 210 j=4,mstate+1
      nj = j-1
      ar(j,1) = rrate(nj,1)
      DO 210 i=4,j-1
      ni = i-1
      ar(j,i) = rrate(nj,ni)
  210 CONTINUE
!
! --- 2s to 1s:
!
      ar(2,1) = zeroc
!
! --- 2p to 1s:
!
      ar(3,1) = (4.0_DP/3.0_DP)*rrate(2,1)
!
! --- 2p to 2s:
!
      ar(3,2) = zeroc
!
      DO 220 j=4,mstate+1
      nj   = j - 1
      ajsq = (FLOAT (nj))**2
!
! --- total radiation to 2s+2p
!
      tot = rrate(nj,2)
!
! --- fraction to 2s:
!
      frac2s  = 12.0_DP*(ajsq-1.0)*(ajsq-4.0_DP)
      frac2s  = frac2s/(frac2s+ajsq*(ajsq-4.0_DP)+32.0_DP*ajsq*(ajsq-1.0+DP))
!
      ar(j,2) = frac2s*tot
      ar(j,3) = (1.0_DP-frac2s)*tot
  220 CONTINUE
!
 1000 CALL SECOND (cpub)
      cpu1 = cpub - cpua
!
      RETURN
!
      END SUBROUTINE hxinit



      SUBROUTINE hxosc
!----------------------------------------------------------------------------------------
! -- calculate oscillator strengths
!----------------------------------------------------------------------------------------

      USE hexnb_data,                                 ONLY : ms,mc,               &  ! param
                                                             f,ar,                &  ! b3
                                                             iexcit,ilorent,      &
                                                             mstate,ncont            ! b8
      IMPLICIT NONE

      REAL(DP) const,aj,y,ai
 
      INTEGER(I4B) i,j

      const = 32.0_DP / (SQRT (27.0_DP) * pi)
!
      DO 10 i=1,mstate
      ai = FLOAT (i)
      DO 10 j=i+1,mc
      aj = FLOAT (j)
      y  = 1.0 - (ai/aj)**2
      f(i,j) = (const*ai/((aj*y)**3))*hxgf(i,y)
   10 CONTINUE
      RETURN
!
      END SUBROUTINE hxosc

      FUNCTION hxgf (n, y)
!------------------------------------------------------------------------------
! --
!------------------------------------------------------------------------------

      IMPLICIT NONE

      REAL(DP) g0,g1,g2,r,rn,y,hxgf
      INTEGER(I4B) n
!
      g0(r) =   0.9935_DP +0.2328_DP/r-0.1296_DP/(r*r)
      g1(r) = -(0.6282_DP -0.5598_DP/r+0.5299_DP/(r*r))/r
      g2(r) =  (0.3887_DP -1.1810_DP/r+1.4700_DP/(r*r))/(r*r)
!
      rn = FLOAT (n)
      IF (n .EQ. 1)  go to 11
      IF (n .EQ. 2)  go to 12
      hxgf = g0(rn)+g1(rn)/y+g2(rn)/(y*y)
      RETURN
!
   11 hxgf = 1.1330_DP-0.4059_DP/y+.07014_DP/(y*y)
      RETURN
!
   12 hxgf = 1.0785_DP -0.2319_DP/y+.02947_DP/(y*y)
      RETURN
!
      END FUNCTION hxgf



      SUBROUTINE hradin (ncor, numz, nz, izstrp, amz)
!

!
! ----------------------------------------------------------------------
!     if numz ne 0, read the post radiation tables (file 'coronb').
!
!     arrangement of arad(j,k,i) for given i:
!
!    1    2    3    4    5    6    7    8    9    10   11   12   13   14
! 1  t1   t2   t3   t4   0.   t1   t2   t3   t4   0.   t1   t2   t3   t4
! 2  t2   t3   t4   t5   0.   t2   t3   t4   t5   0.   t2   t3   t4   t5
! 3 a(0) a(0) a(0) a(0)  0.  b(0) b(0) b(0) b(0)  0.  c(0) c(0) c(0) c(0
! 4  .    .    .    .    .    .    .    .    .    .    .    .    .    .
! 5  .    .    .    .    .    .    .    .    .    .    .    .    .    .
! 6  .    .    .    .    .    .    .    .    .    .    .    .    .    .
! 7  .    .    .    .    .    .    .    .    .    .    .    .    .    .
! 8 a(5) a(6) a(5) a(5)  0.  b(5) b(5) b(5) b(5)  0.  c(5) c(5) c(5) c(5
! -------------------------------------------------------------------------------------
!


      USE hexnb_data,                                ONLY : mz,                   & ! param
                                                            nouthx,istart,ihxbug, & ! b2
                                                            arad                  ! locrad
 
      IMPLICIT NONE

!      argument list:
      INTEGER(I4B)    numz
      INTEGER(I2B)    ncor
      REAL(DP)        amz(*)
      INTEGER(I4B)    nz(*), izstrp(*)
      INTEGER(I4B)    nchars

!     local storage:

      INTEGER(I2B)    get_next_io_unit
      REAL(DP)        dum(8)
      INTEGER(I4B)    izssum,i,j,k,imp,iz,ia,izsum,iostat
      CHARACTER*2     spec, cdum, cmod
      CHARACTER*4     loop, loop2
!

!
      IF (numz .EQ. 0)  RETURN
      izsum = 0
      DO i=1,numz
        izsum = izsum + izstrp(i)
      END DO
      IF (izsum .EQ. numz)  RETURN
      IF ( numz .GT. mz  )  go to 900
!
      DO imp=1,numz
        DO k=1,15
          DO j=1,8
            arad(j,k,imp) = 0.0
          END DO
        END DO
      END DO
!

      ncor = get_next_io_unit ()
      nchars  =LEN_TRIM(ADJUSTL(adas_xsct_path))
      OPEN (unit = ncor, status = 'OLD', err = 901, iostat = iostat, &
            file = adas_xsct_path(1:nchars) // '/data/coronb')
!
      DO imp=1,numz
!
!       search for desired element in data table
!
        loop = 'cont'
        DO WHILE (loop .EQ. 'cont')
!
!         at first record of data set
!
          READ (ncor, 1010, err=902, END=902)  spec, cmod, iz, ia
          IF (iz .EQ. nz(imp)) THEN
!
!           found correct data set
!
            loop = 'exit'
!
          ELSE
!
!           search for first record of next data set
!
            loop2 = 'cont'
            DO WHILE (loop2 .EQ. 'cont')
              READ (ncor, 1020, err=902, END=902) (dum(j), j=1,8)
              READ (ncor, 1015, err=902, END=902)  cdum, cmod
              IF (cdum .NE. spec)  loop2 = 'exit'
            END DO
            BACKSPACE (unit = ncor)
!
          END IF
        END DO
!
!       found desired element
!
        amz(imp) = ia
!
!       read data for this species.
!       three normal exit conditions:
!         (1)  read last species in file (read to end of file);
!         (2)  number of vectors in data set exceeds array dimension;
!         (3)  first record of next species encountered;
!
        k    = 1
        loop = 'cont'
        DO WHILE (loop .EQ. 'cont')
          READ (ncor, 1020, err=903, END=903) (arad(j,k,imp), j=1,8)
          READ (ncor, 1015, err=903, END=120)  cdum, cmod
          k = k + 1
          IF (cmod .EQ. 'zb')  k = MAX0 (k,  6)
          IF (cmod .EQ. 'zs')  k = MAX0 (k, 11)
          IF (k .GT. 15 .OR. cdum .NE. spec)  loop = 'exit'
        END DO
!
  120   REWIND (unit = ncor)
      END DO

      CLOSE (unit = ncor)
      RETURN
!
! ----------------------------------------------------------------------
! error exits
! ----------------------------------------------------------------------
!
!     number of species exceeds array dimension
!
  900 IF (nouthx .GT. 0 .AND. myid == master)  WRITE (nouthx, 9000) numz
      go to 999
!
!     file does not exist
!
  901 IF (nouthx .GT. 0 .AND. myid == master)  WRITE (nouthx, 9010)  iostat
 9010 FORMAT (' ERROR opening file coronb, iostatus = ', i8)
      go to 999
!
!     could not find desired element
!
  902 IF (nouthx .GT. 0 .AND. myid == master)  WRITE (nouthx, 9020)  nz(imp), ncor
 9020 FORMAT (' element #', i10, ' not found on unit ', i5)
      go to 999
!
!     error reading species data set
!
  903 IF (nouthx .GT. 0 .AND. myid == master)  WRITE (nouthx, 9030)  nz(imp), ncor
!
  999 ihxbug = 1

      CLOSE (unit = ncor)
!
! ----------------------------------------------------------------------
! format statements
! ----------------------------------------------------------------------
!
 1010 FORMAT (a2, 10x, a2, 7x, 2i4)
 1015 FORMAT (a2, 10x, a2)
 1020 FORMAT (2e10.3, 3e15.6 / 3e15.6)
 9000 FORMAT (' ERROR in subroutine HRADIN' / &
           7x, 'number of species = ', i8, 'exceeds array dimensions')
 9030 FORMAT (' ERROR reading element number', i10, ' on unit ', i5)
      RETURN
!
      END SUBROUTINE hradin



      SUBROUTINE hxradi (teev, z, zsq, numimp, iz, izstrp, iwatch)
!------------------------------------------------------------------------
! -- evaluate coronal z, zsq, and radiation.
!------------------------------------------------------------------------

      USE cpub_dat

      USE hexnb_data,                                     ONLY : mz,mi,       & ! param
                                                                 arad           ! locrad

      IMPLICIT NONE
!
! --- evaluate coronal z, zsq, and radiation.

!     argument list 
      REAL(DP) teev,z(mz),zsq(mz) 
      INTEGER(I4B)numimp, iz(mz), izstrp(mz), iwatch

!     local storage:

      REAL(DP) rad(mz),b(3),erad,factor,tl,t1,tekev,bb,cc,dd
      REAL(SP) cpub,cpua
      INTEGER(I4B) j,jj,jp,ihigh,imp,kk,k
 
      DATA erad /1.0_DP/
!
      CALL SECOND (cpua)
      iwatch = 0
      ihigh  = 0
      tekev  = teev*1.0e-3
      DO 60 imp=1,numimp
      IF (izstrp(imp) .EQ. 1)  go to 51
!
! --- find temperature region. test for less than 5 intervals.
!
      DO 10 j=1,5
        IF (arad(2,j,imp) .EQ. 0.0)  go to 10
        jj = j
        IF (tekev .GE. arad(1,j,imp)  .AND. &
            tekev .LT. arad(2,j,imp))  go to 30
   10 CONTINUE
!
! --- temperature out of range
!
      iwatch = 1
      IF (tekev .GE. arad(1,1,imp))  go to 20
!
! --- temperature too low. return zero power
!
!***  rad(imp) = 0.0
!***  z  (imp) = 0.0
!***  zsq(imp) = 0.0
!
! --- interpolation for low te (cdb):
!
      t1 = LOG10 (arad(1,1,imp))
      DO 15 j=1,3
        jp   = 5*j - 4
        b(j) = arad(3,jp,imp) + t1*(arad(4,jp,imp) + t1*(arad(5,jp,imp) &
             + t1*(arad(6,jp,imp) + t1*(arad(7,jp,imp) &
             + t1*arad(8,jp,imp)))))
   15 CONTINUE
      factor   = (tekev/arad(1,1,imp))**erad
      rad(imp) = factor*(10.0**b(1))
      z(imp)   = factor*b(2)
      zsq(imp) = factor*b(3)
      go to 50
!
! --- temperature too high. compute value for tmax, then extrapolate.
!
   20 tl    = LOG10 (arad(2,jj,imp))
      ihigh = 1
   30 IF (ihigh .EQ. 0) &
      tl    = LOG10 (tekev)
!
! --- polynomial fits  bb = sum(a(ik)*log(te)**k) etc.
!
      bb = zeroc
      cc = zeroc
      dd = zeroc
!
      DO kk=1,6
        k  = 7-kk
        bb = bb*tl + arad((k+2),jj,imp)
        cc = cc*tl + arad((k+2),(jj+5),imp)
        dd = dd*tl + arad((k+2),(jj+10),imp)
      END DO
!
      rad(imp) = 0.0
      z(imp)   = cc
      zsq(imp) = dd
      IF (bb .LT. -38.0)  go to 50
      rad(imp) = 10.0**bb
!
! --- for te > tmax extrapolate as for pure Bremsstrahlung
!
      IF (ihigh .EQ. 1)  rad(imp) = &
                         rad(imp) * SQRT (tekev / arad(2,4,imp))
      ihigh = 0
   50 CONTINUE
      go to 60
!
! --- fully stripped option [ izstrp(imp) = 1 ]
!
   51 z(imp)   = iz(imp)
      zsq(imp) = iz(imp)**2
   60 CONTINUE
!
      CALL SECOND (cpub)
      cpu2 = cpu2 + cpub - cpua
!
      RETURN
!
      END SUBROUTINE hxradi



      SUBROUTINE hxsve
!--------------------------------------------------------------------
! -- calculate the rate coefficients for collisions with electrons.
!!--------------------------------------------------------------------
      USE cpub_dat

      USE hexnb_data,                          ONLY : ms,mc,mz,mi,   & ! param
                                                      er0, v0, te,   &
                                                      ti, ami, deni, &
                                                      amz,denz,zcor, &
                                                      zsqcor, dene,  &
                                                      ns, nc, numi,  &
                                                      numz, iz,      &
                                                      izstrp,        &  ! b0
                                                      kdene,kdeni,   &
                                                      kdenz,ksvi,    &
                                                      ksvz,ksve,krad,&
                                                      ngh,ngl,       &  ! b1
                                                      nouthx,istart, &
                                                      ihxbug,        &  ! b2
                                                      en,dg,ae,be,   &
                                                      de1,de2,ge1,   &
                                                      ge2,           &  ! b4
                                                      cii,cei,ciz,   &
                                                      cez,cie,cee,   &
                                                      ccxi,ccxz         ! b6
      IMPLICIT NONE


      REAL(DP) et(mc+1)
      REAL(DP) xgl(64),wgl(64)
      REAL(DP) sinh1,arg,sqpi,facte,ve
      REAL cpua,cpub
      INTEGER(I4B) i,j,ifail,nsp1,kk
      CHARACTER*2 typent

      sinh1(arg) = SINH (arg)/arg
!
      CALL SECOND (cpua)
!
! --- initialize:
!
      IF (istart .EQ. 0)  go to 5

      sqpi = SQRT ( pi )
!
!     e-impact thresholds (temp.)
!
      et(1) = (13.723/13.6)*en(1)
      DO i=2,mc+1
        et(i) = en(i)
      END DO
!
      ifail  =  izero
      typent = 'gl'
      CALL d01bbf (typent, ngl, wgl, xgl, ifail)

      IF ( ifail .EQ. izero) go to 5
      IF (nouthx .GT. izero .AND. myid == master) THEN
        WRITE  (nouthx, 3939) ifail
 3939   FORMAT (' ERROR in HXSVE calling D01BBF: ifail= ', i4)
      END IF
      ihxbug = 3
      RETURN
!
    5 nsp1 = ns+1
!
      DO 205 i=1,ns
      cie(i) = zeroc
      DO 205 j=i+1,nsp1
      cee(i,j) = zeroc
  205 CONTINUE
!
      ve = 5.931e7_DP * SQRT (te)
      IF (kdene .EQ. 0)  go to 1000
      go to (220, 210)  ksve + 1
!
! --- Gauss-Laguerre integrations:
!
  210 facte = (2.0_DP/sqpi)*ve * EXP (-v0*v0/(ve*ve))
      DO 211 kk=1,ngl
      DO 212 i=1,ns
      cie(i) = cie(i) &
        +wgl(kk)*(xgl(kk)+et(i)/te)                    &
        *sinh1(2.0_DP*v0 * SQRT (xgl(kk)+et(i)/te)/ve) &
        *sief(i,et(i)+te*xgl(kk))                      &
        *facte * EXP (-et(i)/te)
      DO 213 j=i+1,ns
      cee(i,j) = cee(i,j) &
               + wgl(kk)*(xgl(kk)+(et(i)-et(j))/te)                    &
                *sinh1(2.0_DP*v0 * SQRT (xgl(kk)+(et(i)-et(j))/te)/ve) &
                *seef(i,j,et(i)-et(j)+te*xgl(kk))                      &
                *facte * EXP (-(et(i)-et(j))/te)
  213 CONTINUE
      DO 214 j=nsp1,nc
      cee(i,nsp1) = cee(i,nsp1) &
        +wgl(kk)*(xgl(kk)+(et(i)-et(j))/te)                    &
        *sinh1(2.0_DP*v0 * SQRT (xgl(kk)+(et(i)-et(j))/te)/ve) &
        *seef(i,j,et(i)-et(j)+te*xgl(kk))                      &
        *facte * EXP (-(et(i)-et(j))/te)
  214 CONTINUE
  212 CONTINUE
  211 CONTINUE
      go to 1000
!
! --- Gauss-Laguerre averages for ionization of 1s, ionization of 2s,
!     and excitations among 1s, 2s, and 2p:
!
  220 facte = (2.0_DP/sqpi)*ve * EXP (-v0*v0/(ve*ve))
      DO 221 kk=1,ngl
      DO 222 i=1,MIN0 (2,ns)
      cie(i) = cie(i) &
        +wgl(kk)*(xgl(kk)+et(i)/te) &
        *sinh1(2.0_DP*v0 * SQRT (xgl(kk)+et(i)/te)/ve) &
        *sief(i,et(i)+te*xgl(kk)) &
        *facte * EXP (-et(i)/te)
      DO 222 j=i+1,MIN0 (3,ns)
      cee(i,j) = cee(i,j) &
        +wgl(kk)*(xgl(kk)+(et(i)-et(j))/te) &
        *sinh1(2.0_DP*v0 * SQRT (xgl(kk)+(et(i)-et(j))/te)/ve) &
        *seef(i,j,et(i)-et(j)+te*xgl(kk)) &
        *facte * EXP (-(et(i)-et(j))/te)
  222 CONTINUE
  221 CONTINUE
!
! --- rates from Vriens & Smeets for other reactions:
!     ionizations--
!
      DO 223 i=3,ns
      cie(i) = cief(i,te)
  223 CONTINUE
!
!     excitations--
!
      DO 224 i=1,ns
        DO 225 j=i+1,ns
          IF (j .LE. 3)  go to 225
          cee(i,j) = ceef(i,j,te)
  225   CONTINUE
        DO 226 j=nsp1,nc
          cee(i,nsp1) = cee(i,nsp1)+ceef(i,j,te)
  226   CONTINUE
  224 CONTINUE
!
 1000 CALL SECOND (cpub)
      cpu5 = cpu5 + cpub - cpua
!
      RETURN
!
      END SUBROUTINE hxsve


      SUBROUTINE hxsvi
!-------------------------------------------------------------------------------
! -- calculate the rate coefficients for collisions with ions and impurities
!-------------------------------------------------------------------------------

      USE cpub_dat

      USE hexnb_data,                                   ONLY : ms,mc,mz,mi,    & ! param
                                                               er0, v0, te, ti,&
                                                               ami, deni, amz, &
                                                               denz, zcor,     &
                                                               zsqcor, dene ,  &
                                                               ns, nc,numi,    &
                                                               numz, iz,       &
                                                               izstrp,         & ! b0
                                                               kdene,kdeni,    &
                                                               kdenz,ksvi,ksvz,&
                                                               ksve,krad,ngh,  &
                                                               ngl,            & ! b1
                                                               nouthx,istart,  &
                                                               ihxbug,         & ! b2
                                                               cii,cei,ciz,    &
                                                               cez,cie,cee,    &
                                                               ccxi,ccxz         ! b6
      IMPLICIT NONE
!

!


      REAL(DP)    vi(mi),vz(mz)

      REAL(DP)     xgh(64),wgh(64),w1(32),ei

      REAL(DP)     zz,zz1,s,sqpi,xp,xm,er0old,facti, &
                   ez,factz
      REAL(SP)    cpua,cpub,cpua2,cpua1
      INTEGER(I4B) ister0,ngh1,nsp1,i,j,ifail,ki,kk,kz

      CHARACTER*2   typent

      DATA          ister0 /1/, zz1 /1.0_DP/,er0old /0.0_DP/
!
! --- initialize:
!
      CALL SECOND (cpua)
!
      IF (istart .EQ. 0)  go to 5
!      pi   = ACOS (-1.0)
      sqpi = SQRT (pi)
!
      ngh1 = (ngh+1)/2
      DO i=1,ngh1
        w1(i) = 2.0
      END DO
      IF (2*ngh1-ngh .EQ. 1)  w1(ngh1) = 1.0
      ifail  = 0
      typent = 'gh'
      CALL d01bbf (typent, ngh, wgh, xgh, ifail)
      IF (ifail  .EQ. 0) go to 5
      IF (nouthx .GT. 0 .AND. myid == master) THEN
        WRITE (nouthx, 3939) ifail
 3939   FORMAT (' ERROR in HXSVI calling D01BBF: ifail = ', i4)
      END IF
      ihxbug = 2
      RETURN
!
    5 CONTINUE
      nsp1 = ns+1
      IF (er0 .NE. er0old) ister0 = 1
!
! --- ion reactions-----------------------------------------------
!
      CALL SECOND (cpua1)
!
! --- first zero out arrays:
!
      IF ( (ister0 .EQ. 1 .OR. istart .EQ. 1) .OR. ksvi .EQ. 1) THEN
        DO 10 ki=1,numi
        DO 10 i=1,ns
        cii(i,ki) = 0.0
        ccxi(i,ki) = 0.0
        DO 10 j=i+1,nsp1
        cei(i,j,ki) = 0.0
   10   CONTINUE
      END IF
!
      DO ki=1,numi
        vi(ki) = 1.3841e6 * SQRT (ti/ami(ki))
      END DO
!
      IF (kdeni .EQ. 0)  go to 100
      IF ( ksvi .EQ. 0)  go to 40
!
! --- gauss-hermite integrations:
!
      DO 20 kk=1,ngh1
      DO 20 ki=1,numi
      xp = (xgh(kk)+v0/vi(ki))**2
      xm = (xgh(kk)-v0/vi(ki))**2
      s = 1.0_DP
      IF (xgh(kk)-v0/vi(ki) .LT. 0.0)  s = -1.0
      ei = ti/ami(ki)
      DO 20 i=1,ns
      cii(i,ki) = cii(i,ki)+wgh(kk)*w1(kk) &
        *(xp*siif(i,ei*xp)-s*xm*siif(i,ei*xm))
      ccxi(i,ki) = ccxi(i,ki)+wgh(kk)*w1(kk) &
        *(xp*scxif(i,ei*xp)-s*xm*scxif(i,ei*xm))
      DO j=i+1,ns
        cei(i,j,ki) = cei(i,j,ki)+wgh(kk)*w1(kk) &
                    *(xp*sezf(i,j,zz1,ei*xp)-s*xm*sezf(i,j,zz1,ei*xm))
      END DO
      DO j=nsp1,nc
        cei(i,nsp1,ki) = cei(i,nsp1,ki)+wgh(kk)*w1(kk) &
                      *(xp*sezf(i,j,zz1,ei*xp)-s*xm*sezf(i,j,zz1,ei*xm))
      END DO
   20 CONTINUE
!
      DO 30 ki=1,numi
      facti = vi(ki)*vi(ki)/(2.0*sqpi*v0)
      DO 30 i=1,ns
      cii(i,ki) = facti*cii(i,ki)
      ccxi(i,ki) = facti*ccxi(i,ki)
      DO 30 j=i+1,nsp1
      cei(i,j,ki) = facti*cei(i,j,ki)
   30 CONTINUE
      go to 100
!
! --- simple multiplication instead of maxwellian averaging:
! ---      cii and cei are independent of plasma parameters in this
! ---      approximation       so as long as er0 doesn't change, use
! ---      the old values.  
! ---      ister0 = 1 ==> restart everything
! ---             = 0 ==> cruise in no-update mode
!
! --- redo everything if er0 changes
!
   40 IF ((ister0 .EQ. 1 .OR. istart .EQ. 1) .OR. ksvi .EQ. 1) THEN
      DO 41 ki=1,numi
      DO 41 i=1,ns
      cii (i,ki) =  siif(i,er0)*v0
      ccxi(i,ki) =  scxif(i,er0)*v0
      DO 42 j=i+1,ns
      cei(i,j,ki) = sezf(i,j,zz1,er0)*v0
   42 CONTINUE
      DO j=nsp1,nc
        cei(i,nsp1,ki) = cei(i,nsp1,ki)+sezf(i,j,zz1,er0)*v0
      END DO
   41 CONTINUE
      END IF
  100 CONTINUE
      CALL SECOND (cpua2)
      cpuii = cpuii + cpua2 - cpua1
!
! --- impurity reactions------------------------------------------
!
! --- zero out arrays:
!
      DO 110 kz=1,numz
      vz(kz) = 1.3841e6 * SQRT (ti/amz(kz))
      DO 110 i=1,ns
      ciz(i,kz) = 0.0
      ccxz(i,kz) = 0.0
      DO 110 j=i+1,nsp1
      cez(i,j,kz) = 0.0
  110 CONTINUE
!
      IF (kdenz .EQ. 0)  go to 1000
      IF (ksvz  .EQ. 0)  go to  140
!
! --- gauss-hermite integrations:
!
      DO 120 kk=1,ngh1
      DO 120 kz=1,numz
      xp = (xgh(kk)+v0/vz(kz))**2
      xm = (xgh(kk)-v0/vz(kz))**2
      s = 1.0
      IF (xgh(kk)-v0/vz(kz) .LT. 0.0)  s = -1.0
      ez = ti/amz(kz)
      zz = zcor(kz)
      DO 120 i=1,ns
      ciz(i,kz) = ciz(i,kz)+wgh(kk)*w1(kk) &
        *(xp*sizf(i,zz,ez*xp)-s*xm*sizf(i,zz,ez*xm))
      ccxz(i,kz) = ccxz(i,kz)+wgh(kk)*w1(kk) &
        *(xp*scxzf(i,zz,ez*xp)-s*xm*scxzf(i,zz,ez*xm))
      DO j=i+1,ns
        cez(i,j,kz) = cez(i,j,kz)+wgh(kk)*w1(kk) &
                    *(xp*sezf(i,j,zz,ez*xp)-s*xm*sezf(i,j,zz,ez*xm))
      END DO
      DO j=nsp1,nc
        cez(i,nsp1,kz) = cez(i,nsp1,kz)+wgh(kk)*w1(kk) &
                       *(xp*sezf(i,j,zz,ez*xp)-s*xm*sezf(i,j,zz,ez*xm))
      END DO
  120 CONTINUE
!
      DO kz=1,numz
        factz = vz(kz)**2/(2.0*sqpi*v0)
        DO i=1,ns
          ciz(i,kz) = factz*ciz(i,kz)
          ccxz(i,kz) = factz*ccxz(i,kz)
          DO j=i+1,nsp1
            cez(i,j,kz) = factz*cez(i,j,kz)
          END DO
        END DO
      END DO
!
      go to 1000
!
! --- simple multiplications:
!
  140 DO kz=1,numz
        DO i=1,ns
          ciz (i,kz) = sizf(i,zcor(kz),er0)*v0
          ccxz(i,kz) = scxzf(i,zcor(kz),er0)*v0
          DO j=i+1,ns
            cez(i,j,kz) = sezf(i,j,zcor(kz),er0)*v0
          END DO
          DO j=nsp1,nc
            cez(i,nsp1,kz) = cez(i,nsp1,kz)+sezf(i,j,zcor(kz),er0)*v0
          END DO
        END DO
      END DO
!
 1000 CALL SECOND (cpub)
      cpu4   = cpu4  + cpub - cpua
      cpuiz  = cpuiz + cpub - cpua2
      ister0 = 0
      er0old = er0
      RETURN
!
      END SUBROUTINE hxsvi

      SUBROUTINE lorent (v0, bperp, nc, ns)
!-----------------------------------------------------------------------
! --
!----------------------------------------------------------------------
      USE cpub_dat
      USE hexnb_data,                            ONLY : ms,        &   ! param
                                                        al,        &   ! b5
                                                        iexcit,    &
                                                        ilorent,   &
                                                        mstate,    &
                                                        ncont          ! b8
      IMPLICIT NONE
!

      REAL(DP)    sl(30),fl,vo,bperp,one, eight,ep,v0,w0,almin,almax,&
                  expl,an3,ani3,el
      REAL(SP)    cpua,cpub
      INTEGER(I4B) nc,ns,n,i,il,ncrit,i1,ni
!
      DATA          w0 /4.1341e16_DP/, almin /1.0e-10_DP/, almax /1.0e15_DP/, &
                    expl /4.0_DP/
      DATA sl/ &
         1.00000000e+00_DP, 2.50000000e-01_DP, 2.77777778e-02_DP, 1.73611111e-03_DP, &
         6.94444444e-05_DP, 1.92901235e-06_DP, 3.93675989e-08_DP, 6.15118733e-10_DP, &
         7.59405843e-12_DP, 7.59405843e-14_DP, 6.27608135e-16_DP, 4.35838982e-18_DP, &
         2.57892889e-20_DP, 1.31578005e-22_DP, 5.84791131e-25_DP, 2.28434036e-27_DP, &
         7.90429189e-30_DP, 2.43959626e-32_DP, 6.75788439e-35_DP, 1.68947110e-37_DP, &
         10*0.0_DP/
!
      fl (n, an3, el) = &
         ((4.0_DP / (an3*el))**(2*n-1)) * EXP (-2.0_DP / (3.0_DP*an3*el)) / an3
!
!     input:
!       v0
!       bperp
!       nc
!     output:
!       ns
!       al(i),i = 1,ns
!
      CALL SECOND (cpua)
!
      one   = 1.0_DP
      eight = 8.0_DP
!
      IF (nc .GE. 2)  go to 1
      ns    = 1
      al(1) = 0.0
      go to 30
!
    1 IF (ilorent .EQ. 0) THEN
        ns = mstate + 1
        DO il=1,ns
          al(il) = 0.0
        END DO
        go to 30
      END IF
!
      IF (bperp .EQ. 0.0)  go to 3
!
! --- ep = ABS (v x b)/(electric field at first Bohr radius)
!
      ep    = v0*bperp / 5.1417e13
      ncrit = 0.5 / (ep**(1.0/expl))
      IF (ncrit .LE. 1)  ns = 1
      IF (ncrit .GE. 2)  ns = ncrit + 1
      IF (ns .LE. (mstate+1))  go to 4
!
    3 ns = mstate + 1
!
    4 DO 10 i=1,ns
        IF (i .GE. 4)  go to 13
        go to (11, 12, 12), i
   11   al(1) = fl (1, one  , ep)
        go to 10
   12   al(i) = fl (2, eight, ep) * 0.5
        go to 10
   13   ni = i-1
        ani3  = (FLOAT (ni))**3
        al(i) = fl (ni, ani3, ep) * sl(ni)
   10 CONTINUE
!
! --- normalize rates; zero out those which are too small
!
      DO i=1,ns
        al(i) = al(i)*w0
        IF (al(i) .LE. almin) &
        al(i) = 0.0
      END DO
!
! --- if some rates are too high, truncate system
!
      DO i=1,ns
        i1 = i
        IF (al(i) .GT. almax)  go to 21
      END DO
      go to 30
!
   21 ns = i1 - 1
!
   30 CALL SECOND (cpub)
      cpu3 = cpu3 + cpub - cpua
      RETURN
!
      END SUBROUTINE lorent




      SUBROUTINE matri
!-----------------------------------------------------------------------------------
! -- form the matrix of rate coefficients which occurs
!    in the system of rate equations.
!    evaluate deexcitation via detailed balance.
! -- vectorized
!-----------------------------------------------------------------------------------

      USE cpub_dat

      USE hexnb_data,                            ONLY :  ms,mc,mz,mi,                       & ! param
                                                         er0, v0, te, ti, ami, deni, amz,   &
                                                         denz, zcor, zsqcor, dene ,ns, nc,  &
                                                         numi, numz, iz, izstrp,            & ! b0
                                                         f,ar,                              & ! b3
                                                         en,dg,ae,be,de1,de2,ge1,ge2,       & ! b4
                                                         al,                                & ! b5
                                                         cii,cei,ciz, cez,cie,cee,          &
                                                         ccxi,ccxz,                         & ! b6
                                                         q                                    ! b7
  
      IMPLICIT NONE

      REAL(DP)     dge(ms+1),dgi(ms+1)
      REAL(SP)     cpua,cpub
      INTEGER(I4B) kz,i,j,jp,ki,nsp1
!
!c*** dump c arrays
!      write (60,3310) cii
!3310  format (/ ' cii= ' / (1x,1p5e11.3))
!      write (60,3315) ccxii
!3315  format (/ ' ccxii= ' / (1x,1p5e11.3))
!      write (60,3320) cei
!3320  format (/ ' cei= ' / (1x,1p5e11.3))
!      write (60,3330) ciz
!3330  format (/ ' ciz= ' / (1x,1p5e11.3))
!      write (60,3335) ccxz
!3335  format (/ ' ccxz= ' / (1x,1p5e11.3))
!      write (60,3340) cez
!3340  format (/ ' cez= ' / (1x,1p5e11.3))
!      write (60,3350) cie
!3350  format (/ ' cie= ' / (1x,1p5e11.3))
!      write (60,3360) cee
!3360  format (/ ' cee= ' / (1x,1p5e11.3))
!
      CALL SECOND (cpua)
!
      nsp1 = ns+1
      DO 10 i=1,ns
      DO 10 j=1,ns
        q(i,j) = 0.0
   10 CONTINUE
!
! --- detailed-balance factors (note that en's are positive):
!
      DO i=1,ns
        dge(i) = dg(i) * EXP (en(i)/te)
        dgi(i) = dg(i) * EXP (en(i)/ti)
      END DO
!
! --- rates due to collisions with electrons
!     radiation rates       lorentz rates:
!
      DO 20 j=1,ns
      q(j,j) = -cie(j)*dene -al(j)
      DO 21 jp=j+1,nsp1
      q(j,j) = q(j,j) -cee(j,jp)*dene
   21 CONTINUE
      DO 22 jp=1,j-1
      q(j,j) = q(j,j) &
        -(dge(jp)/dge(j))*cee(jp,j)*dene &
        -ar(j,jp)
   22 CONTINUE
      DO 23 i=j+1,ns
      q(i,j) = q(i,j) &
        +cee(j,i)*dene
   23 CONTINUE
      DO 24 i=1,j-1
      q(i,j) = q(i,j) &
        +(dge(i)/dge(j))*cee(i,j)*dene &
        +ar(j,i)
   24 CONTINUE
   20 CONTINUE
!
! --- add rates due to collisions with ions:
!
      DO 30 ki=1,numi
      DO 30 j=1,ns
      q(j,j) = q(j,j) &
        -cii(j,ki)*deni(ki)
      DO 31 jp=j+1,nsp1
      q(j,j) = q(j,j) - cei(j,jp,ki)*deni(ki)
   31 CONTINUE
      DO 32 jp=1,j-1
      q(j,j) = q(j,j) - (dgi(jp)/dgi(j))*cei(jp,j,ki)*deni(ki)
   32 CONTINUE
      DO 33 i=j+1,ns
      q(i,j) = q(i,j) + cei(j,i,ki)*deni(ki)
   33 CONTINUE
      DO 34 i=1,j-1
      q(i,j) = q(i,j) + (dgi(i)/dgi(j))*cei(i,j,ki)*deni(ki)
   34 CONTINUE
   30 CONTINUE
!
! --- add rates due to collisions with impurities:
!
      DO 40 kz=1,numz
      DO 40 j=1,ns
      q(j,j) = q(j,j) &
        -ciz(j,kz)*denz(kz)
      DO 41 jp=j+1,nsp1
      q(j,j) = q(j,j) &
        -cez(j,jp,kz)*denz(kz)
   41 CONTINUE
      DO 42 jp=1,j-1
      q(j,j) = q(j,j) &
        -(dgi(jp)/dgi(j))*cez(jp,j,kz)*denz(kz)
   42 CONTINUE
      DO 43 i=j+1,ns
      q(i,j) = q(i,j) &
        +cez(j,i,kz)*denz(kz)
   43 CONTINUE
      DO 44 i=1,j-1
      q(i,j) = q(i,j) &
        +(dgi(i)/dgi(j))*cez(i,j,kz)*denz(kz)
   44 CONTINUE
   40 CONTINUE
!
      CALL SECOND (cpub)
      cpu6 = cpu6 + cpub - cpua
      RETURN
!
      END SUBROUTINE matri



      SUBROUTINE nbsgxn (ebkev, ebfac, ibion, mb, nebin, vbeam, &
                         debin, hxfrac,atw_beam)
! ----------------------------------------------------------------------------------
!
!  this subroutine calculates the neutral beam attenuation array, sgxn.
!
!     input:
!
!          ebkev(mb)      - maximum energy of mb-th neutral beamline
!                           (keV)
!          ebfac          - factor defining upper bound on energy bin
!                           range,  > 1.0
!          atw_beam         mass no. of beam
!
!          mb             - number of beamlines modeled
!          mfm1           - number of flux zones minus one
!          nebin          - number of energy bins (rotation case only)
!          vbeam(ie,mb)   - speed of ie-th energy group of the mb-th
!                           beamline (cm/sec)
!          zne(mfm1)      - local electron density (cm**-3)
!          zni(mfm1,nion) - local density of nion-th ion species
!                           (cm**-3)
!          zzi(mfm1,nion) - local average charge state of nion-th ion
!                           species
!     input from common /io/:
!          ncrt,nouthx,ncorin
!     input from common /ions/:
!          namei, atw
!     input from common /nub3/:
!          iexcit,ilorent,mstate,ncont,kdene,kdeni,kdenz,ksvi,ksvz,ksve,
!          krad,ngh,ngl,znipm,atwpm,iz,zti,izstrp
!     input from common /numbrs/:
!          nprim,nimp,nion
!
!     output:
!          sgxn(i,j,k,l)
!               i - mode index
!                   = 1, fraction of reactions producing electrons;
!                   = 2, fraction of reactions producing species 1 ion;
!                   = 3, fraction of reactions producing species 2 ion;
!                   = 4, total inverse mean free path;
!               j - FREYA zone index
!               k - beam index
!                   k = 3*(ib-1) + ie, where ib is the beam index and ie
!                                      is the beam energy group index
!               l - index for relative energy.  bins are equispaced in
!                   delta(energy) between 0 and max(ebkev(ib))*ebfac,
!                   with bin width given by
!                             delta(energy) = ebmax*ebfac/nebin.
!                   ebfac .ge. 1.0 and nebin are user supplied namelist
!                   variables.
!          sgxnmi(ie,mb)
!                 - minimum inverse mean free path for ie-th energy
!                   group of mb-th beamline.  calculated for stationary
!                   target case only.  otherwise, calculated for each
!                   neutral trajectory in subroutine INJECT.
!          hxfrac(ie,mb)
!                 - fractional density of n = 3 excited state neutral.
!                   calculated in stationary target case only.
!          debin
!                 - width of energy bins (keV/amu)
! ---------------------------------------------------------------------------------
!
      USE nf_param,                                           ONLY : ke,kcmp1,kbe

      USE neutral_beams,                                   ONLY : ilorent,mstate,iexcit,ncont, &
                                                                  ngl, kdeni, kdenz, ksvi,     &
                                                                  ksvz, ksve, krad, ngh,ncont, &
                                                                  kdene,ncorin

      USE  zonal_data,                                     ONLY : mfm1,zne,zni,zte,zzi,zti

      IMPLICIT NONE


!     argument list:
      REAL(DP)     ebkev(*),vbeam(ke,*)

      REAL(DP)     hxfrac(ke,*)
      REAL(DP)     sgxne(kbe)
      REAL(DP)     atw_beam,debin,ebfac
      INTEGER(I4B) ibion,nebin,mb


!     local vars
      INTEGER(I4B) k,i,istart,j,iz0,ib,ie,ind,ihxerr
      REAL(DP) atwpm(nion),atwim(nion),eova,vbin,teev,tiev,ebmax,     &
               bperp,rerate,rmfp,hexfrac,ebin
      INTEGER(I2B) nouthx,  get_next_io_unit


!
! ----------------------------------------------------------------------
! Go right to ADAS routines and return if  iexcit = 5
! ----------------------------------------------------------------------
!

      IF (iexcit .EQ. 5) THEN

         CALL adassgxn (ebkev,ebfac,ibion,mb,nebin,vbeam, &
                        debin,sgxnmi,atw_beam)

        RETURN
      END IF

!
! ---------------------------------------------------- HSJ-2/5/98 ------
! --- Boley's parameterization including MSI effects:
! ----------------------------------------------------------------------
!
      IF (iexcit .EQ. 6) THEN

        CALL wrap_xboley (atw,ebkev,ibion,mb,nebin,ebfac, &
                          debin,nion,sgxnmi,atw_beam)
        RETURN
      END IF
!
! ----------------------------------------------------------------------
! HEXNB initialization
! ----------------------------------------------------------------------
!

 
      IF (iexcit .NE. 0) THEN
         nouthx = get_next_io_unit ()
        IF (nouthx .GT. 0) THEN
          OPEN (unit = nouthx, file = 'hexnbout', status = 'UNKNOWN')
        END IF
        istart = 1
!
!       separate primary and impurity ions
!
        DO j=1,nion
          IF (j .LE. nprim) THEN
            atwpm(j) = atw(j)
          ELSE
            k        = j - nprim
            atwim(k) = atw(j)
          END IF
        END DO
!
!       get atomic number of impurities
!
 
!        iz(:) = izero
!        DO j=1,nimp_ml
!           DO i=1,nimp
!              IF(namei(i) == namei_ml(j))THEN       ! check for inclusion in master list
!                 IF (namei(i) .EQ. 'he')  iz(i) =  2
!                 IF (namei(i) .EQ. 'c ')  iz(i) =  6
!                 IF (namei(i) .EQ. 'o ')  iz(i) =  8
!                 IF (namei(i) .EQ. 'si')  iz(i) = 14
!                 IF (namei(i) .EQ. 'ar')  iz(i) = 18
!                 IF (namei(i) .EQ. 'cr')  iz(i) = 24
!                 IF (namei(i) .EQ. 'fe')  iz(i) = 26
!                 IF (namei(i) .EQ. 'ni')  iz(i) = 28
!                 IF (namei(i) .EQ. 'kr')  iz(i) = 36
!                 IF (namei(i) .EQ. 'mo')  iz(i) = 42
!                 IF (namei(i) .EQ. 'w ')  iz(i) = 74
!              ENDIF
!           ENDDO
!        END DO


        iz(:) = izero
        DO i=1,nimp
           iz(i) = atomnoi_ml(nimp_index(i))
        ENDDO
 

        iz0 = MINVAL(iz)
        IF(iz0 == izero)THEN
           IF(myid == master)    &
                WRITE(ncrt,FMT='("sub nbsgxn : master impurity list error" )')
           lerrno = 219 + iomaxerr
           CALL terminate(lerrno,nlog)
        ENDIF





!
!       no Lorentz ionization limit
!
        bperp = zeroc
      END IF
!
! ----------------------------------------------------------------------
! stationary plasma case
! ----------------------------------------------------------------------
!
      IF (nebin .EQ. izero) THEN
        DO   ib=1,mb
          DO ie=1,ke
            sgxnmi(ie,ib) = 0.0
            ind  = ke*(ib-1) + ie
            eova = 1.0e3_DP*ebkev(ib)/(FLOAT (ie)*atw_beam)
            vbin = vbeam(ie,ib)
            DO i=1,mfm1
              teev = 1.0e3_DP*zte(i)
              IF (iexcit .NE. 0) THEN
                tiev = 1.0e3_DP*zti(i)
                DO j=1,nion
                  IF (j .LE. nprim) THEN
                    znipm(j) = zni(i,j)
                  ELSE
                    k        = j - nprim
                    zniim(k) = zni(i,j)
                  END IF
                END DO
              END IF
              sgxne(ind) = fsgxne(vbin,teev,zne(i))
              sgxn(1,i,ind,1) = 0.0
              DO k=1,nion
                sgxn(1,i,ind,1) = sgxn(1,i,ind,1) &
                                + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
                IF ((k .LE. nprim) .AND. (k .LE. 2.0)) &
                  sgxn(k+1,i,ind,1) = fsgxncx(atw(k),eova,zni(i,k))
              END DO
              sgxn(1,i,ind,1) = sgxn(1,i,ind,1) + sgxne(ind)
              sgxn(4,i,ind,1) = sgxn(1,i,ind,1) + sgxn(2,i,ind,1) &
                              + sgxn(3,i,ind,1)
              sgxn(1,i,ind,1) = sgxn(1,i,ind,1)/sgxn(4,i,ind,1)
              sgxn(2,i,ind,1) = sgxn(2,i,ind,1)/sgxn(4,i,ind,1)
              sgxn(3,i,ind,1) = sgxn(3,i,ind,1)/sgxn(4,i,ind,1)
              IF (iexcit .NE. 0) THEN
                CALL hexnb (istart, 1, ilorent, mstate, ncont, eova, &
                            teev, tiev, nprim, atwpm, znipm, nimp, iz, &
                            atwim, izstrp, zniim, bperp, kdene, kdeni, &
                            kdenz, ksvi, ksvz, ksve, krad, ngh, ngl, &
                            nouthx, ncorin, rerate, rmfp, hexfrac, &
                            ihxerr)
                IF (ihxerr .NE. izero  .AND. myid ==  master) WRITE (ncrt, 1000) ihxerr

                IF (ihxerr .EQ. 1) THEN
                  IF(myid == master) &
                           WRITE (ncrt, 1010)
                  lerrno = 220 + iomaxerr
                  CALL terminate(lerrno,nlog)
                END IF

                sgxn(4,i,ind,1) = 1.0/rmfp
                IF (i .EQ. 1) hxfrac(ie,ib) = hexfrac
              END IF
              sgxnmi(ie,ib) = MAX (sgxnmi(ie,ib),sgxn(4,i,ind,1))
            END DO
          END DO
        END DO
        DO   ib=1,mb
          DO ie=1,ke
            sgxnmi(ie,ib) = 1.0 / sgxnmi(ie,ib)
          END DO
        END DO
      ELSE
!
! ----------------------------------------------------------------------
! rotating plasma case
! ----------------------------------------------------------------------
!
        ebmax = zeroc
        DO ib=1,mb
          ebmax = MAX (ebmax,ebkev(ib))
        END DO
        ebmax = ebmax/atw_beam
        debin = ebmax*ebfac/FLOAT (nebin)
!
        DO 340 i=1,mfm1
        teev = 1.0e3 * zte(i)
        IF (iexcit .NE. 0) THEN
          tiev = 1.0e3 * zti(i)
          DO j=1,nion
            IF (j .LE. nprim) THEN
              znipm(j      ) = zni(i,j)
            ELSE
              zniim(j-nprim) = zni(i,j)
            END IF
          END DO
        END IF
!
!       electron impact ionization independent of rotation speed
!
        DO   ib=1,mb
          DO ie=1,ke
            ind        = ke*(ib-1) + ie
            sgxne(ind) = fsgxne(vbeam(ie,ib),teev,zne(i))
          END DO
        END DO
!
        DO 340 j=1,nebin
        ebin = FLOAT (j) * debin
        eova =     1.0e3 * ebin
        sgxn(1,i,1,j) = 0.0
        DO k=1,nion
          sgxn(1,i,1,j) = sgxn(1,i,1,j) &
                        + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
          IF ((k .LE. nprim) .AND. (k .LE. 2.0)) &
            sgxn(k+1,i,1,j) = fsgxncx(atw(k),eova,zni(i,k))
        END DO
        IF (iexcit .NE. 0) THEN
          CALL hexnb (istart, 1, ilorent, mstate, ncont, eova, &
                      teev, tiev, nprim, atwpm, znipm, nimp, iz, &
                      atwim, izstrp, zniim, bperp, kdene, kdeni, &
                      kdenz, ksvi, ksvz, ksve, krad, ngh, ngl, nouthx, &
                      ncorin, rerate, rmfp, hexfrac, ihxerr)
          IF (ihxerr .NE. 0 .AND. myid == master) WRITE (ncrt, 1000) ihxerr
          IF (ihxerr .EQ. 1) THEN
            IF(myid == master)WRITE (ncrt, 1010)
                  lerrno = 220 + iomaxerr
                  CALL terminate(lerrno,nlog)
          END IF
        END IF
        DO 340 k=ind,1,-1
          sgxn(1,i,k,j) = sgxn(1,i,1,j) + sgxne(k)
          sgxn(4,i,k,j) = sgxn(1,i,k,j) + sgxn(2,i,1,j) + sgxn(3,i,1,j)
          sgxn(1,i,k,j) = sgxn(1,i,k,j)/sgxn(4,i,k,j)
          sgxn(2,i,k,j) = sgxn(2,i,1,j)/sgxn(4,i,k,j)
          sgxn(3,i,k,j) = sgxn(3,i,1,j)/sgxn(4,i,k,j)
          IF (iexcit .NE. 0)  sgxn(4,i,k,j) = 1.0/rmfp
  340   CONTINUE
      END IF
!
 1000 FORMAT (' subroutine NBSGXN reports a HEXNB return code of ', i5)
 1010 FORMAT (' ERROR: execution terminated - file "coronb" not found')
!
 
      RETURN
!
      END SUBROUTINE nbsgxn

      FUNCTION nfhx (i)

      IMPLICIT  NONE
      INTEGER(I4B) nfhx,i

      nfhx = i - 1
      IF (i .LE. 2) nfhx = i
      RETURN

      END FUNCTION nfhx




      SUBROUTINE nbsgold (ebkev,ebfac,ibion,mb,nebin,vbeam,    &
                          debin,sgxnmi,atw_beam)
!
!
! ----------------------------------------------------------------------
!
!     This subroutine calculates the neutral beam attenuation array sgxn
!     using the old cross sections based on Freeman & Jones (1972).
!     No excited state effects are included. These cross sections have
!     been found to be outdated and are not to used except for purposes
!     of comparison to the newer ADAS data.
!
!     Created:     12-JUL-1995     Daniel Finkenthal
!
!     input:
!
!          ebkev(mb)      - maximum energy of mb-th neutral beamline
!                           (keV)
!          ebfac          - factor defining upper bound on energy bin
!                           range,  > 1.0
!          ibion          - index of beam ion species
!                           (e.g., atwb = atw(ibion))
!          mb             - number of beamlines modeled
!          mfm1           - number of flux zones minus one
!          nebin          - number of energy bins (rotation case only)
!          vbeam(ie,mb)   - speed of ie-th energy group of the mb-th
!                           beamline (cm/sec)
!          zne(mfm1)      - local electron density (cm**-3)
!          zni(mfm1,nion) - local density of nion-th ion species
!                           (cm**-3)
!          zzi(mfm1,nion) - local average charge state of nion-th ion
!                           species
!     input from common /io/:
!          ncrt,nout,nouthx,ncorin
!     input from common /ions/:
!          namei, atw
!     input from common /numbrs/:
!          nprim,nimp,nion
!
!     output:
!          sgxn(i,j,k,l)
!               i - mode index
!                   = 1, fraction of reactions producing electrons;
!                   = 2, fraction of reactions producing species 1 ion;
!                   = 3, fraction of reactions producing species 2 ion;
!                   = 4, total inverse mean free path;
!               j - FREYA zone index
!               k - beam index
!                   k = 3*(ib-1) + ie, where ib is the beam index and ie
!                                      is the beam energy group index
!               l - index for relative energy.  bins are equispaced in
!                   delta(energy) between 0 and max(ebkev(ib))*ebfac,
!                   with bin width given by
!                             delta(energy) = ebmax*ebfac/nebin.
!                   ebfac .ge. 1.0 and nebin are user supplied namelist
!                   variables.
!          sgxnmi(ie,mb)
!                 - minimum inverse mean free path for ie-th energy
!                   group of mb-th beamline.  calculated for stationary
!                   target case only.  otherwise, calculated for each
!                   neutral trajectory in subroutine INJECT.
!          debin
!                 - width of energy bins (keV/amu)
! ----------------------------------------------------------------------
!
       USE nf_param,                                        ONLY : kcmp1,kbe,ke

       USE  zonal_data,                                     ONLY : mfm1,zne,zni,zte,zzi

      IMPLICIT NONE

!     argument list:

      REAL(DP) ebkev(*),vbeam(ke,*)


      REAL(DP) sgxnmi(ke,*)
 
      REAL(DP) ebfac,debin,atw_beam

      REAL(DP) sgxne(kbe)
       
      INTEGER(I4B) ibion,mb,nebin


!     local vars:

      INTEGER(I4B) i,ib,ie,ind,j,k

      REAL(DP) eova,vbin,teev,ebmax,ebin

 
!
! ----------------------------------------------------------------------
! stationary plasma case
! ----------------------------------------------------------------------
!
      IF (nebin .EQ. 0) THEN
        DO 220 ib=1,mb
        DO 220 ie=1,ke
        sgxnmi(ie,ib) = 0.0
        ind  = ke*(ib-1) + ie
        eova = 1.0e3*ebkev(ib)/(FLOAT (ie)*atw_beam)
        vbin = vbeam(ie,ib)
        DO 220 i=1,mfm1
        teev = 1.0e3 * zte(i)
        sgxne(    ind  ) = fsgxne(vbin,teev,zne(i))
        sgxn (1,i,ind,1) = 0.0
        DO k=1,nion
          sgxn(1,i,ind,1) = sgxn(1,i,ind,1) &
                          + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
          IF ((k .LE. nprim) .AND. (k .LE. 2.0)) &
            sgxn(k+1,i,ind,1) = fsgxncx(atw(k),eova,zni(i,k))
        END DO
        sgxn(1,i,ind,1) = sgxn(1,i,ind,1) + sgxne(ind)
        sgxn(4,i,ind,1) = sgxn(1,i,ind,1) + sgxn(2,i,ind,1) &
                        + sgxn(3,i,ind,1)
        sgxn(1,i,ind,1) = sgxn(1,i,ind,1)/sgxn(4,i,ind,1)
        sgxn(2,i,ind,1) = sgxn(2,i,ind,1)/sgxn(4,i,ind,1)
        sgxn(3,i,ind,1) = sgxn(3,i,ind,1)/sgxn(4,i,ind,1)
        sgxnmi(ie,ib) = MAX (sgxnmi(ie,ib),sgxn(4,i,ind,1))
  220   CONTINUE
        DO   ib=1,mb
          DO ie=1,ke
            sgxnmi(ie,ib) = 1.0 / sgxnmi(ie,ib)
          END DO
        END DO
      ELSE
!
! ----------------------------------------------------------------------
! rotating plasma case
! ----------------------------------------------------------------------
!
        ebmax = zeroc
        DO ib=1,mb
           ebmax = MAX (ebmax,ebkev(ib))
        END DO
        ebmax = ebmax/atw_beam
        debin = ebmax*ebfac/FLOAT (nebin)
!
        DO 340 i=1,mfm1
        teev = 1.0e3*zte(i)
!       electron impact ionization independent of rotation speed
        DO 320 ib=1,mb
        DO 320 ie=1,ke
        ind = ke*(ib-1) + ie
        sgxne(ind) = fsgxne(vbeam(ie,ib),teev,zne(i))
  320   CONTINUE
        DO 340 j=1,nebin
        ebin = FLOAT (j)*debin
        eova = 1.0e3*ebin
        sgxn(1,i,1,j) = zeroc
        DO 330 k=1,nion
        sgxn(1,i,1,j) = sgxn(1,i,1,j) &
                      + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
        IF ((k .LE. nprim) .AND. (k .LE. 2 )) &
          sgxn(k+1,i,1,j) = fsgxncx(atw(k),eova,zni(i,k))
  330   CONTINUE
           
        DO 340 k=ind,1,-1       ! (ind=ke*(ib-1) +ie,beam energy index)
          sgxn(1,i,k,j) = sgxn(1,i,1,j) + sgxne(k)
          sgxn(4,i,k,j) = sgxn(1,i,k,j) + sgxn(2,i,1,j) + sgxn(3,i,1,j)
          sgxn(1,i,k,j) = sgxn(1,i,k,j)/sgxn(4,i,k,j)
          sgxn(2,i,k,j) = sgxn(2,i,1,j)/sgxn(4,i,k,j)
          sgxn(3,i,k,j) = sgxn(3,i,1,j)/sgxn(4,i,k,j)
  340   CONTINUE
      END IF

      RETURN
!
      END  SUBROUTINE nbsgold

      SUBROUTINE eigen (xeig)
!------------------------------------------------------------------------------------
! --
!------------------------------------------------------------------------------------
      USE cpub_dat

      USE hexnb_data,                                   ONLY :  er0, v0, te, ti, ami, &
                                                                deni, amz,denz, zcor, &
                                                                zsqcor, dene,ns, nc,  &
                                                                numi, numz, iz,       &
                                                                izstrp,nouthx,istart, &
                                                                ihxbug, q,eigvr,      &
                                                                eigvl,cj,xfrac,ms,mz, &
                                                                mi

      USE replace_imsl,                                 ONLY :  my_eigrf

!
!      IMPLICIT  INTEGER (i-n), REAL*8 (a-h, o-z)
!
      IMPLICIT NONE

 
      REAL(DP) xeig,xtest,eigmin,sum,cond
      REAL(SP) cpua,cpub
      REAL(DP)  qq(ms+1,ms+1), work2(3*(ms+1))

      REAL(DP) w(2*(ms+1)),z(2*(ms+1),ms+1)

      REAL(DP) cjmat(ms+1,ms+1), workdc(ms+1)

      REAL(DP)   infac(ms+1,ms+1), lexp(ms+1), in(ms+1)

      INTEGER(I4B)  i,j,ijob,ifail1, ipvt(ms+1)
!
      DATA      xtest /1.0_DP/

!     on return from eigrf, w contains the complex eigenvalues
!     as w = [wr(1),wi(1), ... ,wr(ns),wi(ns)]
!     on return from eigrf, z contains the complex eigenvectors
!     as z(col1) =  [zr(1),zi(1), ... ,zr(ns),zi(ns)], etc.,
!     each column representing an eigenvector
!
      CALL SECOND (cpua)
!
      DO 10 i=1,ns
      DO 10 j=1,ns
        qq(i,j) = q(i,j)
   10 CONTINUE
      IF (nouthx .GT. 0 .AND. myid == master) THEN
        WRITE  (nouthx, 11)
   11   FORMAT ( 'q _ matrix')
        DO i=1,ns
          WRITE  (nouthx, 12)  (q(i,j), j=1,ns)
   12     FORMAT (1x, 1p7e11.3)
        END DO
      END IF

      ijob = 1
!      CALL eigrf (qq, ns, ms+1, ijob, w, z, ms+1, work2, ifail1)
      CALL my_eigrf (qq, ns, ms+1, ijob, w, z, ms+1, work2, ifail1)
      IF (ifail1 .EQ. 0)  go to 15
      ihxbug = 9
      IF (nouthx .GT. 0 .AND. myid == master ) THEN
        WRITE  (nouthx, 3939) ifail1
 3939   FORMAT (' ERROR in EIGRF: ifail1 = ', i4)
      END IF
      RETURN
!
   15 DO j=1,ns
        DO i=1,ns
          eigvr(i,j) = z(2*i-1,j)
        END DO
        eigvl(j) = w(j*2-1)
      END DO
!
      eigmin = 1.0e30_DP
      DO i=1,2*ns,2
        eigmin = MIN (eigmin, ABS (w(i)))
      END DO
!
!***  do 20 i=1,ns
!***    eigmin = MIN (eigmin, ABS (eigr(i)))
!**20 continue
!
      xeig = zeroc
      IF (eigmin .NE. zeroc )  xeig = v0 / eigmin
!
!     determine fraction of 3rd excited state (approximate)
!
      sum = zeroc
      DO i=1,ns
        sum = sum + eigvr(i,1)
      END DO
!
      xfrac = eigvr(3,1) / sum
      IF (nouthx .EQ. 0)  RETURN
!
      IF(myid == master)THEN
         WRITE  (nouthx, 1000)  xeig
1000     FORMAT (' length (eigenvalue):  ', 1pe11.4)
         WRITE  (nouthx, 1010)  (i, eigvl(i), i=1,ns)
1010     FORMAT ('    i= ', i2, '         eigenvalue= ', 1pe11.4)
         DO j=1,ns
            WRITE  (nouthx, 1020)  j, (eigvr(i,j), i=1,ns)
1020        FORMAT (' eigenvector[', i2, ']= ' / (10x, 1p5e11.4))
         END DO
      ENDIF
!
!     apply initial condition constraint
!
      DO 1200 i=1,ns
      DO 1200 j=1,ns
        cjmat(i,j) = eigvr(i,j)
 1200 CONTINUE
!
      cj(1) = 1.0
      DO i=2,ns
        cj(i) = 0.0
      END DO
!
      CALL decomp (ms+1, ns, cjmat, cond, ipvt, workdc)
      CALL solveq (ms+1, ns, cjmat,   cj, ipvt)
!
      IF(myid == master)WRITE  (nouthx, 1320)  (cj(i), i=1,ns)
 1320 FORMAT (' bc constants: ' / (10x,1p5e11.4))
!
      DO j=1,ns
        DO i=1,ns
          lexp(i)    = eigvl(i) / v0
          infac(j,i) = cj(i) * eigvr(j,i)
        END DO
        IF(myid == master)WRITE  (nouthx, 1400)  j, (infac(j,i), lexp(i), i=1, ns)
 1400   FORMAT ( ' *** i(',i2,') =    ' / &
               (12x, 1pe15.7, ' * EXP (', 1pe15.7, ' * x)') )
      END DO
!
      DO j=1,ns
        in(j) = 0.0
        DO i=1,ns
          in(j) = in(j) + infac(j,i) * EXP (lexp(i)*xtest)
        END DO
      END DO
!
      WRITE  (nouthx, 1450)
 1450 FORMAT (/// ' qnj * ij / v0 ...........' /)
!
      DO j=1,ns
        sum = 0.0
        DO i=1,ns
          sum = sum + q(j,i) * in(i)
        END DO
        IF(myid == master)WRITE  (nouthx, 1500)  j, sum/v0
 1500   FORMAT (15x,i3,5x,1pe15.7)
      END DO
!
      CALL SECOND (cpub)
      cpu7 = cpu7 + cpub - cpua
!
      RETURN
!
      END  SUBROUTINE eigen


      SUBROUTINE d01bbf (TYPE, inum, weight, abscis, ier)
!----------------------------------------------------------------------
! --
!----------------------------------------------------------------------
      USE d01bff_dat
      USE hexnb_data,                       ONLY : nouthx, istart, ihxbug ! b2

      IMPLICIT NONE

      INTEGER(I4B)  inum, num(4), iorg(4),i,ier,irel,index

      REAL(DP) weight(inum), abscis(inum)
!
      CHARACTER*2 TYPE



      DATA        num  /10, 16, 20, 24/
      DATA        iorg / 1, 11, 27, 47/



      DO i=1,4
        index = i
        IF (inum .EQ. num(i))  go to 60
      END DO
      ier = 1
      RETURN
!
   60 ier = 0
!
      IF (TYPE .EQ. 'gl') THEN
        DO i=1,inum
          irel      = iorg(index) - 1 + i
          weight(i) = wgl(irel)
          abscis(i) = xgl(irel)
        END DO
        go to 400
      END IF
!
      IF (TYPE .EQ. 'gh') THEN
        DO i=1,inum
          irel      = iorg(index) - 1 + i
          weight(i) = wgh(irel)
          abscis(i) = xgh(irel)
        END DO
        go to 400
      END IF
!
      ier = 1
  400 RETURN
!
      END SUBROUTINE d01bbf

      FUNCTION dfhx (beta)
!---------------------------------------------------------------------
! -- approximation to the function defined by Janev & Presnyakov,
! -- j. phys. b 13, 4233 (1980).
!---------------------------------------------------------------------


      IMPLICIT  NONE

      REAL(DP) dd,beta,beta1,a,dfhx
      INTEGER ia

      DIMENSION dd(38)
      DATA      dd/                                                 &
        .10450_DP, .121_DP, .138_DP, .157_DP, .175_DP, .200_DP,     &
        .229_DP, .260_DP, .300_DP, .339_DP,                         &
        .367_DP, .389_DP, .402_DP, .410_DP, .398_DP, .376_DP,       &
        .346_DP, .317_DP, .285_DP, .255_DP,                         &
        .227_DP, .205_DP, .185_DP, .168_DP, .152_DP, .138_DP,       &
        .124_DP, .110_DP, .099_DP, .089_DP,                         &
        .079_DP, .070_DP, .062_DP, .054_DP, .047_DP, .041_DP,       &
        .035_DP, .02898_DP/

      beta1 = 1.0_DP / beta
      IF (beta1 .LE. 0.2_DP  )  go to 110
      IF (beta1 .GE. 1.0e3_DP)  go to 120
      a    = 10.0_DP * LOG10 (beta1) + 8.0
      ia   = MIN0 (37, INT (a))
      dfhx = dd(ia)+(a-FLOAT (ia))*(dd(ia+1)-dd(ia))
      RETURN

  110 dfhx = 0.5 * beta * (1.0 - 1.0/(8.0 * beta * SQRT (beta)))    &
                        * EXP (-SQRT (2.0 * beta))
      RETURN

  120 dfhx = 4.0 * beta * LOG (1.4 * beta1)
      RETURN

      END FUNCTION dfhx


      SUBROUTINE decomp (ndim, n, a, cond, ipvt, work)
!--------------------------------------------------------------------------------
!
!
!     decomposes a real matrix by gaussian elimination
!     and estimates the condition of the matrix.
!
!     use solveq to compute solutions to linear systems.
!
!     input..
!
!        ndim = declared row dimension of the array containing  a.
!
!        n = order of the matrix.
!
!        a = matrix to be triangularized.
!
!     output..
!
!        a  contains an upper triangular matrix  u  and a permuted
!          version of a lower triangular matrix  i-l  so that
!          (permutation matrix)*a = l*u .
!
!        cond = an estimate of the condition of  a .
!           for the linear system  a*x = b, changes in  a  and  b
!           may cause changes  cond  times as large in  x .
!           if  cond+1.0 .eq. cond , a is singular to working
!           precision.  cond  is set to  1.0e+32  if exact
!           singularity is detected.
!
!        ipvt = the pivot vector.
!           ipvt(k) = the index of the k-th pivot row
!           ipvt(n) = (-1)**(number of interchanges)
!
!     work space..  the vector  work  must be declared and included
!                   in the call.  its input contents are ignored.
!                   its output contents are usually unimportant.
!
!     the determinant of a can be obtained on output by
!        det(a) = ipvt(n) * a(1,1) * a(2,2) * ... * a(n,n).
!--------------------------------------------------------------------------------



      IMPLICIT NONE

      INTEGER(I4B)  ipvt(n), ndim, n
      INTEGER(I4B)  nm1, i, j, k, kp1, kb, km1, m
      REAL(DP)      a(ndim,n), cond, work(n)
      REAL(DP)      ek, t, anorm, ynorm, znorm
!
      ipvt(n) = 1
      IF (n .EQ. 1)  go to 80
      nm1 = n - 1
!
!     compute 1-norm of a
!
      anorm = 0.0
      DO j=1,n
        t = 0.0
        DO i=1,n
          t = t + ABS (a(i,j))
        END DO
        IF (t .GT. anorm)  anorm = t
      END DO
!
!     gaussian elimination with partial pivoting
!
      DO 35 k=1,nm1
         kp1 = k+1
!
!        find pivot
!
         m = k
         DO i=kp1,n
           IF (ABS (a(i,k)) .GT. ABS (a(m,k))) m = i
         END DO
         ipvt(k) = m
         IF (m .NE. k) ipvt(n) = -ipvt(n)
         t = a(m,k)
         a(m,k) = a(k,k)
         a(k,k) = t
!
!        skip step if pivot is zero
!
         IF (t .EQ. 0.0)  go to 35
!
!        compute multipliers
!
         DO i=kp1,n
           a(i,k) = -a(i,k)/t
         END DO
!
!        interchange and eliminate by columns
!
         DO 30 j=kp1,n
           t = a(m,j)
           a(m,j) = a(k,j)
           a(k,j) = t
           IF (t .EQ. 0.0)  go to 30
           DO i=kp1,n
             a(i,j) = a(i,j) + a(i,k)*t
           END DO
   30    CONTINUE
   35 CONTINUE
!
!     cond = (1-norm of a)*(an estimate of 1-norm of a-inverse)
!     estimate obtained by one step of inverse iteration for the
!     small singular vector.  this involves solving two systems
!     of equations, (a-transpose)*y = e  and  a*z = y  where  e
!     is a vector of +1 or -1 chosen to cause growth in y.
!     estimate = (1-norm of z)/(1-norm of y)
!
!     solve (a-transpose)*y = e
!
      DO k=1,n
        t = 0.0
        IF (k .EQ. 1)  go to 45
        km1 = k-1
        DO i=1,km1
          t = t + a(i,k)*work(i)
        END DO
   45   ek = 1.0
        IF (     t .LT. 0.0)  ek = -1.0
        IF (a(k,k) .EQ. 0.0)  go to 90
        work(k) = -(ek + t)/a(k,k)
      END DO
!
      DO kb=1,nm1
        k   = n - kb
        t   = 0.0
        kp1 = k+1
        DO i=kp1,n
          t = t + a(i,k) * work(k)
        END DO
        work(k) = t
        m       = ipvt(k)
        IF (m .NE. k) THEN
          t       = work(m)
          work(m) = work(k)
          work(k) = t
        END IF
      END DO
!
      ynorm = 0.0
!
      DO i=1,n
         ynorm = ynorm + ABS (work(i))
      END DO
!
!     solve a*z = y
!
      CALL solveq(ndim, n, a, work, ipvt)
!
      znorm = 0.0
!
      DO i=1,n
        znorm = znorm + ABS (work(i))
      END DO
!
!     estimate condition
!
      cond = anorm*znorm/ynorm
      IF (cond .LT. 1.0)  cond = 1.0
      RETURN
!
!     1-by-1
!
   80 cond = 1.0
      IF (a(1,1) .NE. 0.0)  RETURN
!
!     exact singularity
!
   90 cond = 1.0e+32
      RETURN
!
      END SUBROUTINE decomp


     FUNCTION fsgxne (vb, te, zne)
! ----------------------------------------------------------------------
! this subprogram evaluates local inverse mean free path for electron impact
! ionization from fitted results of freeman and jones, clm-r137,culham (1974).
!     input:
!             vb  - speed of impinging neutral (cm/sec)
!             te  - target electron temperature (ev)
!             zne - target electron density (cm**-3)
! ----------------------------------------------------------------------
       IMPLICIT NONE

       REAL(DP) vb,te,zne,alogt,expo,fsgxne

      REAL(DP)   cfione(7)
      DATA       cfione /-3.173850e+01_DP,  1.143818e+01_DP, -3.833998_DP,   &
                        7.046692e-01_DP, -7.431486e-02_DP,  4.153749e-03_DP, &
                       -9.486967e-05_DP/

      alogt = zeroc
      IF (te .GT. 1.0_DP    )  alogt = LOG (te)
      IF (te .GT. 1.0e+05_DP)  alogt = 11.51_DP
      expo = (((((cfione(7) *alogt + cfione(6))*alogt + cfione(5))*alogt     &
                + cfione(4))*alogt + cfione(3))*alogt + cfione(2))*alogt     &
                + cfione(1)
      fsgxne = EXP (expo) * zne / vb
      RETURN

      END  FUNCTION fsgxne 

      FUNCTION fsgxni (atw, eova, zni, zzi)
! ----------------------------------------------------------------------
!  this subprogram calculates inverse mean free path due to proton and
!     impurity impact ionization.  proton impact ionization cross
!     sections from fitted results of freeman and jones,clm-r137, culham
!     (1974).  impurity impact ionization cross sections from r.e. olson
!     et al., Phys. Rev. Lett. 41, 163 (1978).

!     input:
!             atw  - atomi! weight of target ion
!             eova - relative energy of impinging neutral (ev/amu)
!             zni  - density of target ion (cm**-3)
!             zzi  - average charge state of target ion
! ----------------------------------------------------------------------

       IMPLICIT NONE

       REAL(DP) fsgxni,atw,eova,zni,zzi,aloge,expo,sigi,ekev
       REAL(DP) cfionp(7)
       DATA      cfionp /-4.203309e+01_DP,  3.557321_DP    , -1.045134_DP,          &
                         3.139238e-01_DP, -7.454475e-02_DP,  8.459113e-03_DP,       &
                        -3.495444e-04_DP/

      IF (atw .LE. 3.01_DP) THEN
        aloge = LOG10 (eova)
        aloge = aloge * 2.302585093_DP - 6.907755279_DP
        IF (aloge .LE. -2.30258_DP) THEN
          sigi = 0.0
        ELSE
          expo = (((((cfionp(7) *aloge + cfionp(6))*aloge        &
                   + cfionp(5))*aloge + cfionp(4))*aloge         &
                   + cfionp(3))*aloge + cfionp(2))*aloge         &
                   + cfionp(1)
          sigi = EXP (expo)
        END IF
        fsgxni = sigi*zni
      ELSE
        ekev   = 1.0e-3_DP*eova
        fsgxni = 1.0e-17_DP*zni*46.0_DP*zzi*(32.0_DP*zzi/ekev)*   &
                    (1.0_DP - EXP (-ekev/(32.0_DP*zzi)))
      END IF
      RETURN

      END FUNCTION fsgxni


      FUNCTION fsgxncx (atw, e, zni)
! ----------------------------------------------------------------------
! this subprogram calculates inverse mean free path due to charge exchange
! from the fitted results of freeman and jones, clm-r137, culham (1974).
!     input:
!             atw - atomi! weight of target ion
!             e   - relative energy of impinging neutral (ev/amu)
!             zni - density of target ion (cm**-3)
! ----------------------------------------------------------------------
      IMPLICIT NONE

      REAL(DP) atw,e,zni,aloge,sigcx,fsgxncx

      IF (atw .GT. 3.01_DP) THEN
        sigcx = zeroc
      ELSE
        aloge = LOG10 (e)
        sigcx = 0.6937e-14_DP * (1.0_DP - 0.155_DP*aloge)**2 /     &
                             (1.0_DP + 0.1112e-14_DP*e**3.3)
      END IF
      fsgxncx = sigcx*zni
      RETURN

      END  FUNCTION fsgxncx


      SUBROUTINE adassgxn (ebkev, ebfac, ibion, mb, nebin, vbeam, &
                           debin,sgxnmi, atw_beam)
!
!
! ----------------------------------------------------------------------
!
! ADASSGXN calculates the effective neutral beam stopping cross sections
! using the JET Atomi! Data and Analysis Structure (ADAS). The cross
! sections are returned in array sgxn for each beam, beam energy
! component, and FREYA flux zone.
!
! The plasma ions are assumed to be fully stripped. Only ion species
! H, He, B, Be, C, O, N, and Ne are available from ADAS. Any other ion
! will terminate the code. This can be improved later.
!
! The routine was made to be as compatible as possible with the existing
! code. A call is made to the original cross section package in order to
! determine relative deposition fractions only. The HEXNB package is
! totally avoided.
!
! Reference:   Finkenthal, Ph.D. Thesis, UC-Berkeley, 1994
!
! Created  :   06-jul-1995    by  D. Finkenthal
! ----------------------------------------------------------------------
!
!     input:
!
!          ebkev(mb)      - full energy of mb-th neutral beamline
!                           (keV)
!          ebfa!          - factor defining upper bound on energy bin
!                           range,  .gt. 1.0
!          ibion          - index of beam ion species
!                           (e.g. atwb = atw(ibion))
!                           if ibion = -1 beam is dt mixture
!          atw_beam         (use atw_beam for mass in this case)
!          mb             - number of beamlines modeled
!          mfm1           - number of flux zones minus one
!          nebin          - number of energy bins (rotation case only)
!          vbeam(ie,mb)   - speed of ie-th energy group of the mb-th
!                           beamline (cm/sec)
!          zne(mfm1)      - local electron density (cm**-3)
!          zni(mfm1,nion) - local density of nion-th ion species
!                           (cm**-3)
!          zte(mfm1)      - local electron temperature (KeV)
!          zti0(mfm1)     - local ion temperature (KeV)
!          zzi(mfm1,nion) - local average charge state of nion-th ion
!                           species
!     input from common /io/:
!          ncrt
!     input from common /ions/:
!          namei, atw
!     input from common /numbrs/:
!          nprim,nimp,nion
!
!     output:
!          sgxn(i,j,k,l)
!               i - mode index
!                   = 1, fraction of reactions producing electrons;
!                   = 2, fraction of reactions producing species 1 ion;
!                   = 3, fraction of reactions producing species 2 ion;
!                   = 4, total inverse mean free path;
!               j - FREYA zone index
!               k - beam index
!                   k = 3*(ib-1) + ie, where ib is the beam index and ie
!                                      is the beam energy group index
!               l - index for relative energy.  bins are equispaced in
!                   delta(energy) between 0 and max(ebkev(ib))*ebfac,
!                   with bin width given by
!                             delta(energy) = ebmax*ebfac/nebin.
!                   ebfac .ge. 1.0 and nebin are user supplied namelist
!                   variables.
!          sgxnmi(ie,mb)
!                 - minimum inverse mean free path for ie-th energy
!                   group of mb-th beamline.  calculated for stationary
!                   target case only.  otherwise, calculated for each
!                   neutral trajectory in subroutine INJECT.
!          debin
!                 - width of energy bins (keV/amu)
!                   Used only for rotatation case.
! ----------------------------------------------------------------------
!


       USE nf_param,                                  ONLY : kion,kbe,ksge,kb,kcmp1

       USE  zonal_data,                              ONLY : mfm1,zne,zni,zte,zzi,zti0 => zti

      IMPLICIT NONE
!
!     argument list:
      REAL(DP) ebkev(kb),vbeam(ke,kb)

      REAL(DP) sgxnmi(ke,kb)
 
      REAL(DP) ebfac,debin,atw_beam
      INTEGER(I4B) ibion,mb,nebin
!
!     local storage:
      REAL(DP)  sgxeff(3), atwb,tiev, &
                ecol,ebmax,ebin,sgxnscale
      !cnz(kz,20), zeffx(kz)

      REAL(DP),ALLOCATABLE,DIMENSION(:,:)      :: cnz
      REAL(DP),ALLOCATABLE,DIMENSION(:)        :: zeffx

      INTEGER(I4B) init,izbm,i,j,k,ierr,ie,ib,ind,jreff
!      SAVE init, izatom, atwb, izbm

      DATA init /0/
!
      IF(.NOT. ALLOCATED(cnz))ALLOCATE(cnz(mfm1,20))
      cnz(:,:) = 0.0_DP
      IF(.NOT. ALLOCATED(zeffx))ALLOCATE(zeffx(mfm1))
      zeffx(:) = 0.0_DP
      IF (init .EQ. 0) THEN
!
! ----------------------------------------------------------------------
! determine atomic  number of primary ions
! ----------------------------------------------------------------------
!

      DO i=1,nprim
        izatom(i) = 1
        IF (namep(i) .EQ. 'he') &
        izatom(i) = 2
      END DO
!
! ----------------------------------------------------------------------
! determine atomic number of impurity ions
! ----------------------------------------------------------------------
!
 
      IF (nimp .EQ. 0)  go to 3430
      DO i=1,nimp
        k = nprim + i
        izatom(k) = 0
        IF (namei(i) .EQ. 'he')  izatom(k) =  2
        IF (namei(i) .EQ. 'c' )  izatom(k) =  6
        IF (namei(i) .EQ. 'o' )  izatom(k) =  8
        IF (namei(i) .EQ. 'si')  izatom(k) = 14
        IF (namei(i) .EQ. 'ar')  izatom(k) = 18
        IF (namei(i) .EQ. 'cr')  izatom(k) = 24
        IF (namei(i) .EQ. 'fe')  izatom(k) = 26
        IF (namei(i) .EQ. 'ni')  izatom(k) = 28
        IF (namei(i) .EQ. 'kr')  izatom(k) = 36
        IF (namei(i) .EQ. 'mo')  izatom(k) = 42
        IF (namei(i) .EQ. 'w' )  izatom(k) = 74
      END DO
!
! Get beam species
!
 3430 IF (ibion .GT. 0) THEN
        atwb = atw(ibion)
        izbm = izatom(ibion)
      ELSE
        atwb = atw_beam
        izbm = 1
      END IF



!
! ----------------------------------------------------------------------
! Check to make sure that the impurities requested are compatable with
! ADAS (i.e., H, He, B, Be, C, O, N, or Ne). Terminate with error if an
! invalid impurity is listed.
! ----------------------------------------------------------------------
!

      DO k=1,nion
        IF (izatom(k) .GT. 10) THEN
            IF(myid == master)WRITE (ncrt, 1010)
            lerrno = 202 + iomaxerr
            CALL terminate(lerrno,nlog)
        END IF
 1010 FORMAT &
         (' *** Execution Terminated:'                                 / &
          '     The ADAS database only contains atomic cross-sections' / &
          '     for fully-stripped H, He, B, Be, C, O, N, and Ne ions.'/ &
          '     You must restrict your choice of impurity species to'  / &
          '     these ions.')
      END DO
!
      init = 1
      END IF    ! init
!
! ----------------------------------------------------------------------
! Get original cross sections-Used to determine relative deposition
! fractions only. Total cross section is determined using ADAS.
! ----------------------------------------------------------------------
!

      CALL nbsgold (ebkev, ebfac, ibion, mb, nebin, vbeam, &
                    debin, sgxnmi, atw_beam)


!
! ----------------------------------------------------------------------
! Set up the cnz arrays (concentrations of plasma ions) and Zeffx
! zni and zne are the (FREYA zone) densities of electron and ions:
! ----------------------------------------------------------------------
!
      DO i=1, mfm1
        DO k = 1, nion
          cnz(i,k) = 0.0
          IF (izatom(k) .EQ.  1)  cnz(i, 1) = zni(i,k)/zne(i)
          IF (izatom(k) .EQ.  2)  cnz(i, 2) = zni(i,k)/zne(i)
          IF (izatom(k) .EQ.  4)  cnz(i, 4) = zni(i,k)/zne(i)
          IF (izatom(k) .EQ.  5)  cnz(i, 5) = zni(i,k)/zne(i)
          IF (izatom(k) .EQ.  6)  cnz(i, 6) = zni(i,k)/zne(i)
          IF (izatom(k) .EQ.  7)  cnz(i, 7) = zni(i,k)/zne(i)
          IF (izatom(k) .EQ.  8)  cnz(i, 8) = zni(i,k)/zne(i)
          IF (izatom(k) .EQ. 10)  cnz(i,10) = zni(i,k)/zne(i)
        END DO

!
! The Zeff(i=1,...mfm1) array has been stored in zzi(i,nion+1)
        zeffx(i) = zzi(i,nion+1)
!
      END DO


!
! ----------------------------------------------------------------------
!     stationary plasma case
! ----------------------------------------------------------------------
!
      IF (nebin .EQ. 0) THEN
        ierr = 0
        DO    i=1,mfm1
          DO ib=1,mb
            tiev = 1.0e3_DP * zti0(i)
            ecol = 1.0e3_DP * ebkev(ib) / (atwb)

            CALL adasqh6 (tiev,ecol,izbm,0,zne(i),zeffx(i), &
                          cnz(i,2),cnz(i,4),cnz(i,5),cnz(i,6), &
                          cnz(i,7),cnz(i,8),cnz(i,10), &
                          sgxeff,ierr)
            IF (ierr .EQ. 1) THEN
               IF(myid == master)WRITE (ncrt, FMT='("subroutine ADASSGXN: problem #1")')
               lerrno = 222 + iomaxerr
               CALL terminate(lerrno,nlog)
            END IF
            DO ie = 1,ke
              ind  = ke*(ib-1) + ie
              sgxn(4,i,ind,1) = zne(i)*sgxeff(ie)
              sgxnmi(ie,ib)   = MAX (sgxnmi(ie,ib),sgxn(4,i,ind,1))
            END DO
          END DO
        END DO
!
      ELSE
!
! ----------------------------------------------------------------------
!     rotating plasma case
! ----------------------------------------------------------------------
!

        ebmax = zeroc
        DO ib=1,mb
          ebmax = MAX (ebmax,ebkev(ib))
        END DO
        ebmax = ebmax/atwb
        debin = ebmax*ebfac/FLOAT (nebin)
!
        DO i=1,mfm1
          tiev  = 1.0e3*zti0(i)
          jreff = 0
          DO j=1,nebin
            ebin = FLOAT (j) * debin
            ecol =     1.0e3 *  ebin
!
!  Only call ADAS if beam energy (ecol) is greater than 5 keV/amu.
!  Skip over the energy bins that are less than 5 keV/amu. These
!  low energy bins will be calculated using the old cross sections
!  scaled to the adas data at some reference energy.
!
            IF (ecol .GE. 5000.0) THEN
 
              CALL adasqh6 (tiev,ecol,izbm,1,zne(i),zeffx(i), &
                            cnz(i,2),cnz(i,4),cnz(i,5),cnz(i,6), &
                            cnz(i,7),cnz(i,8),cnz(i,10), &
                            sgxeff,ierr)
              IF (ierr .EQ. 1) THEN
                 IF(myid == master)WRITE (ncrt, FMT='("subroutine ADASSGXN: problem #1")')
                 lerrno = 223 + iomaxerr
                 CALL terminate(lerrno,nlog)
             END IF
!
! Calculate a scaling factor from the first available (ge. 5.0)
!   ADAS bin energy This will be used to scale the old cross
!   sections to fit in with ADAS for the low energy bins.
!
              IF (jreff .EQ. 0) THEN
                jreff     = j
                sgxnscale = zne(i)*sgxeff(1)/sgxn(4,i,1,j)
              END IF
              DO k=1,ke*mb
                sgxn(4,i,k,j) = zne(i)*sgxeff(1)
              END DO
            END IF
          END DO

!
! Now scale the old cross sections to match the first available ADAS
!   datapoint for the energy bins below 5.0 keV/amu
!
          IF (jreff .GT. 1) THEN
            DO  j=1,jreff-1
             DO k=1,ke*mb
               sgxn(4,i,k,j) = sgxn(4,i,k,j)*sgxnscale
             END DO
            END DO
          END IF
        END DO
      END IF
      RETURN
!

      END SUBROUTINE adassgxn



      SUBROUTINE adasqh6 (ti, ecol, bmz, iflag, ne, zeff, conche, &
                          concbe, concb, concc, concn, conco, concne, &
                          qrat, ierr)

!
!  Calculate and return the rate coefficient for beam stopping
!  This routine is a modified version of L.D. Horton's (JET)
!  QHIOCH6 routine used to access and evaluate the ADAS beam
!  stopping cross sections from ion-specific effective data files.
!
!   changes from QHIOCH6:
!      1) The need for the NAG library has been eliminated using
!         the spline routines SPLINE and SEVAL.
!      2) Only the full-energy component can be calculated by
!         setting iflag=1. This allows lower energies to be
!         read from the ion-specific files.
!
!                                  D.F. Finkenthal 6-JUL-95
!
!   changes from QHIOCH5:
!      1) modified to force rereading of input files when beam
!         species has changed from the last call (using ipass)
!
!                                  L.D. Horton    2 June 92
!
!   changes from QHIOCH4:
!      1) modified to accept beam species as input and to then
!         calculate the beam stopping for either hydrogen or
!         helium beams
!
!                                  L.D. Horton   16 August 91
!
!   changes from QHIOCH3:
!      1) incorporated Fritsch's new calculation for excitation by
!         Be4+ and He2+
!      2) improved calculation of Maxwellian averages in bundled-n
!         code
!      3) proper inclusion of beam energy dependence
!          - stopping rate is now read from a matrix on beam energy
!            and target density; only the temperature dependence is
!            done separately
!      4) included boron, nitrogen and neon as input concentrations
!          - the code skips all zero concentrations
!
!                                  L.D. Horton   31 July 91
!
!   also:  back to Wilhelm's 3 energies-at-a-time so that all spline
! fits can be done at once. Only if the beam energy has changed by
! more than 1% will the density splines be redone.  Since only one
! energy is used per bank, this means that the spline will be done
! only twice per call of ATTS4.
!
!                                  L.D. Horton   14 August 91
!
! ti        : REAL   : ion temperature in eV
! ecol      : REAL   : collision (=beam) energy in eV/amu
! bmz       : INTEGER: beam nuclear charge
! iflag     : INTEGER: flag for full energy calculation only
! ne        : REAL   : electron density in m**-3
! conche    : REAL   : relative Helium concentration   ( 0 < CHE < 1)
! concbe    : REAL   : relative Beryllium concentration
! concb     : REAL   : relative Boron concentration
! concc     : REAL   : relative Carbon concentration
! concn     : REAL   : relative Nitrogen concentration
! conco     : REAL   : relative Oxygen concentration
! concne    : REAL   : relative Neon concentration
! qrat(3)   : REAL   : requested cross section m**2 for full, half,
!                      and third energy beam components
!
!




      IMPLICIT NONE

! argument list:
  REAL(DP) ti,ecol,zeff, ne, qrat(3),zeff8
  REAL(DP)     conche, concbe, concb, concc, concn, conco, concne

  INTEGER(I2B) get_next_io_unit,nunadas
! ipass : file read switch.  Reread files if beam species has changed
!
      INTEGER    ipass
      CHARACTER  dsn(8,2)*35
      DATA       ipass/0/
      DATA       dsn  / '/data/adas/h_h1.dat'  , &
                        '/data/adas/h_he2.dat' , &
                        '/data/adas/h_be4.dat' , &
                        '/data/adas/h_b5.dat'  , &
                        '/data/adas/h_c6.dat'  , &
                        '/data/adas/h_n7.dat'  , &
                        '/data/adas/h_o8.dat'  , &
                        '/data/adas/h_ne10.dat', &
                        '/data/adas/he_h1.dat' , &
                        '/data/adas/he_he2.dat', &
                        '/data/adas/he_be4.dat', &
                        '/data/adas/he_b5.dat' , &
                        '/data/adas/he_c6.dat' , &
                        '/data/adas/he_n7.dat' , &
                        '/data/adas/he_o8.dat' , &
                        '/data/adas/he_ne10.dat'/
!
! Physics Constants
!

      REAL(DP) ,PARAMETER :: ev = 1.6022e-12
!      PARAMETER (amu = 1.6605e-24)
!      PARAMETER (eV  = 1.6022e-12)
      REAL(DP)      amu

!
! Local variables
!

      INTEGER    bmz, iflag
      INTEGER,SAVE :: nebeam
 


      INTEGER(I4B)    nsp, maxe, maxn, maxt, ierr,nchars
      PARAMETER (nsp = 8)              ! 8 different ion species
      PARAMETER (maxe = 15, maxn = 10, maxt = 10)
      INTEGER(I4B)    isp, z(nsp), neb(nsp), ndens(nsp), ntemp(nsp), i, j, k
      REAL(DP)     eb(maxe,nsp), dens(maxn,nsp), temp(maxt,nsp)
      REAL(DP)     tref(nsp), ebref(nsp), denref(nsp), svref(nsp)
      REAL(DP)     sven(maxe,maxn,nsp), svt(maxt,nsp)
      CHARACTER  line*80
      DATA       z /1, 2, 4, 5, 6, 7, 8, 10/

      REAL(DP)     ti8, ne8, qrat8(3)
      REAL(DP),    SAVE  :: ecol8
      REAL(DP)     conc(nsp)
!
      INTEGER(I4B)     ifail ! apparently this was intended as an error flag on return
                         ! from pline_coeff. But such a return value does not exist
                         ! in spline_coeff. Hence it is simply commented out here HSJ

      REAL(DP)     be(maxe,maxn,nsp),ce(maxe,maxn,nsp),de(maxe,maxn,nsp)
      REAL(DP)     bt(maxt,nsp),ct(maxt,nsp),dt(maxt,nsp)
      REAL(DP)     bn(maxn,nsp,3),cn(maxn,nsp,3),dn(maxn,nsp,3)
      REAL(DP)     svintn(maxn,nsp,3),svintt,svtot(nsp),svtcor(nsp)
!
      REAL(DP) zeffm1,vbeam
!
      amu =  AMU_Value *1000._DP            ! grams

      ierr    = 0
      ti8     = ti
      zeff8   = zeff
      ne8     = ne * 1.0e-13
!***  conc(1) = conch
      conc(2) = conche
      conc(3) = concbe
      conc(4) = concb
      conc(5) = concc
      conc(6) = concn
      conc(7) = conco
      conc(8) = concne
!
! open and read input file only once
!
      IF (ipass .NE. bmz) THEN
        IF(myid == master)THEN
           WRITE (ncrt, 1100)
           WRITE (nlog, 1100)
        ENDIF
!
        ecol8 = ecol
!

        DO isp=1,nsp   ! loop over species
          nchars  =LEN_TRIM(ADJUSTL(adas_xsct_path))
          PRINT *,'file =',TRIM(ADJUSTL(adas_xsct_path))//dsn(isp,bmz)
          nunadas = get_next_io_unit ()
          OPEN (unit = nunadas, &
                file = adas_xsct_path(1:nchars)//dsn(isp,bmz), &
                                               status = 'OLD')
          READ (nunadas,  1000) z(isp),svref(isp)
          READ (nunadas, '(a)') line
          READ (nunadas,  1001) neb(isp),ndens(isp),tref(isp)
          READ (nunadas, '(a)') line
          READ (nunadas,  1002) (eb(j,isp),j=1,neb(isp))
          READ (nunadas,  1002) (dens(k,isp),k=1,ndens(isp))
          READ (nunadas, '(a)') line

          DO k=1,ndens(isp)
            READ (nunadas, 1002) (sven(j,k,isp),j=1,neb(isp))
          END DO

          READ  (nunadas, '(a)') line
          READ  (nunadas,  1003) ntemp(isp),ebref(isp),denref(isp)
          READ  (nunadas, '(a)') line
          READ  (nunadas,  1002) (temp(j,isp),j=1,ntemp(isp))
          READ  (nunadas, '(a)') line
          READ  (nunadas,  1002) (svt(j,isp),j=1,ntemp(isp))
!          CALL giveupus(nunadas)
          CLOSE (nunadas)
!
          DO k = 1, ndens(isp)
            dens(k,isp) = dens(k,isp) * 1.0e-13
          END DO
!
! spline the data on energy and temperature
!
          DO k=1,ndens(isp)
            ifail = 0
            CALL spline_coef (neb(isp),eb(1,isp),sven(1,k,isp), &
                         be(1,k,isp),ce(1,k,isp),de(1,k,isp))
!            IF (ifail .NE. 0) THEN
!              WRITE (6, '(a)')  ' spline error #1 in ADASQH6'
!              ierr = 1
!              RETURN
!            END IF
          END DO
          ifail = 0
          CALL spline_coef (ntemp(isp),temp(1,isp),svt(1,isp), &
                       bt(1,isp),ct(1,isp),dt(1,isp))
!          IF (ifail .NE. 0) THEN
!            WRITE (6, '(a)')  ' spline error #2 in ADASQH6'
!            ierr = 1
!            RETURN
!          END IF
!
! Determine if only the full energy component is requested
! Default is that all three beam energy components are to be determined
! Only the full energy component is interesting for He beams
!

          nebeam = 3
          IF (iflag .GT. 0)  nebeam = 1
          IF (  bmz .EQ. 2)  nebeam = 1
!
! spline on density for each requested energy component (1 or all 3)
!
          DO   i=1,nebeam
            DO k=1,ndens(isp)
              svintn(k,isp,i) = seval(neb(isp),ecol8/i,eb(1,isp), &
                   sven(1,k,isp),be(1,k,isp),ce(1,k,isp),de(1,k,isp))
!
            END DO
!
            ifail = 0
            CALL spline_coef(ndens(isp),dens(1,isp),svintn(1,isp,i), &
                        bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
!
!            IF (ifail .NE. 0) THEN
!              WRITE (6, '(a)')  ' spline error #3 in ADASQH6'
!              ierr = 1
!              RETURN
!            END IF
          END DO
        END DO
!
 1000   FORMAT (i5, 8x, d9.3)
 1001   FORMAT (2i5, 7x, d9.3)
 1002   FORMAT (8(1x, d9.3) / 8(1x, d9.3))
 1003   FORMAT (i5, 7x, d9.3, 7x, d9.3)
 1100   FORMAT (/ &
       ' Using ADAS, the effective beam stopping cross sections:' / &
       ' ***** Opening and Reading the Atomic Data Tables *****'  /)
!
        ipass = bmz
      ELSE
!
!  redo density splines only if beam energy has changed by more than 1%
!
        IF (ABS (ecol8-ecol) / ecol8 .GT. 0.01) THEN
          DO isp=1,nsp
            ecol8 = ecol
            DO i=1,nebeam
              DO k=1,ndens(isp)
                svintn(k,isp,i) = seval(neb(isp),ecol8/i,eb(1,isp), &
                      sven(1,k,isp),be(1,k,isp),ce(1,k,isp),de(1,k,isp))
              END DO
              ifail = 0
              CALL spline_coef (ndens(isp),dens(1,isp),svintn(1,isp,i), &
                           bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
!              IF (ifail .NE. 0) THEN
!                WRITE (6, '(a)')  ' spline error #4 in ADASQH6'
!                ierr = 1
!                RETURN
!              END IF
            END DO
          END DO
        END IF
      END IF
!
! calculate correction to requested temperature
!
      DO isp=1,nsp
        IF (ti8 .LE. temp(1,isp)) THEN
          svtcor(isp) = svt(1,isp)/svref(isp)
        ELSE
          svintt      = seval(ntemp(isp),ti8,temp(1,isp), &
                        svt(1,isp),bt(1,isp),ct(1,isp),dt(1,isp))
          svtcor(isp) = svintt/svref(isp)
        END IF
      END DO
!
! scale the input concentrations to match the required zeff
!
           zeffm1 = 0.0
           DO isp=2,nsp
             zeffm1 = zeffm1 + z(isp)*(z(isp)-1)*conc(isp)
           END DO
           IF (zeffm1 .GT. 1.0e-5) THEN
             DO isp=2,nsp
               conc(isp) = (zeff8-1.0) / zeffm1 * conc(isp)
             END DO
           END IF
           conc(1) = 1.0
           DO isp=2,nsp
             conc(1) = conc(1) - z(isp)*conc(isp)
           END DO
!
! evaluate at three energy components
!
      DO i=1,nebeam
!
! interpolate to requested density
!
        DO isp=1,nsp
          IF (ne8 .LE. dens(1,isp)) THEN
            svtot(isp) = seval(ndens(isp),dens(1,isp),dens(1,isp), &
                  svintn(1,isp,i),bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
          ELSE
            svtot(isp) = seval(ndens(isp),ne8,dens(1,isp), &
                  svintn(1,isp,i),bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
          END IF
        END DO

!
!  construct the total stopping cross section
!
        qrat8(i) = 0.0
        DO isp=1,nsp
          qrat8(i) = qrat8(i) + svtot(isp)*svtcor(isp)*z(isp)*conc(isp)
        END DO
!
!  divide by the beam speed to get a cross section and convert to m**2
!
        vbeam   = SQRT (2.0 * ecol8 / i * ev / amu)
        qrat(i) = qrat8(i) / vbeam
      END DO
      RETURN
!
      END SUBROUTINE  adasqh6


      FUNCTION polyf (c, n, x)

!------------------------------------------------------------------
!  -- evaluate the polynomial c(1)+c(2)*x+...+c(n)*x**(n-1).
!------------------------------------------------------------------
      USE cpub_dat

      IMPLICIT NONE

      INTEGER(I4B) j,n

      REAL(DP) c(n),polyf,x,polysv
      REAL(SP) cpuaa,cpubb

      CALL SECOND (cpuaa)
      polyf = c(n)
      DO j=1,n-1
        polyf  = polyf*x+c(n-j)
        polysv = polyf
      END DO
      CALL SECOND (cpubb)
      cpuplf = cpuplf + cpubb - cpuaa
      RETURN

      END  FUNCTION polyf

      FUNCTION poly3f (a, ndim1, ndim2, ndim3,        &
                            x1, n1, x2, n2, x3, n3)
!--------------------------------------------------------------
! --- Modified by J. Mandrekas, for the SuperCode
! --- New version with new fits, 03/27/92
!-------------------------------------------------------------
      IMPLICIT NONE

      INTEGER(i4B) ndim1, ndim2, ndim3, n1,n2, n3, i, i1, i2, i3


      REAL(DP)     a(ndim1,ndim2,ndim3),poly3f

      REAL(DP)     x1, x2, x3, y1(4), y2(4), y3(4)



      y1(1) = 1.0
      DO i = 2, n1
        y1(i) = y1(i-1)*x1
      END DO

      y2(1) = 1.0
      DO i = 2, n2
        y2(i) = y2(i-1)*x2
      END DO

      y3(1) = 1.0
      DO i = 2, n3
        y3(i) = y3(i-1)*x3
      END DO

      poly3f = 0.0
      DO i1 = 1, n1
        DO i2 = 1, n2
          DO i3 = 1, n3
            poly3f = poly3f + a(i1,i2,i3)*y1(i1)*y2(i2)*y3(i3)
          END DO
        END DO
      END DO
      RETURN

      END  FUNCTION poly3f


      FUNCTION sdaccf (ni, nj, z, er)
!---------------------------------------------------------
! -- dacc cross section.
!     Janev & Presnyakov, j. phys. b 13, 4233 (1980).
!---------------------------------------------------------

      USE hexnb_data,                      ONLY : ms,mc,     & ! param
                                                  f,ar         ! b3
      IMPLICIT  NONE
      REAL(DP) omega,z,er,alam,beta,sdaccf
      INTEGER(I4B) ni,nj

      omega  = 0.5_DP * (1.0_DP/ ((FLOAT (ni))**2) - 1.0_DP / ((FLOAT (nj))**2))
      alam   = SQRT (f(ni,nj)/(2.0_DP*omega))
      beta   = z*alam*omega/(er/24982.0_DP)
      sdaccf = (1.7595e-16_DP)*(z*alam/omega)*dfhx(beta)
      RETURN

      END FUNCTION sdaccf 


      FUNCTION scxzf (i, z, er)
!----------------------------------------------------------------------
! -- cross section for electron loss via collision WITH impurity.
!----------------------------------------------------------------------

      USE hexnb_data,                           ONLY : kdene,kdeni,kdenz,  &
                                                       ksvi,ksvz,ksve,krad,&
                                                       ngh,ngl               ! b1
      IMPLICIT  NONE

      REAL(DP) z,er,ansq,xtilde,scxzf
      INTEGER(I4b) i,ni,id

      ni = nfhx(i)
      id = 2
      IF (ni .EQ. 1)  id = 1

      ansq   = (FLOAT (ni))**2
      xtilde = ansq * er / (31250.0_DP * z)
      scxzf  = zeroc
      IF (xtilde .LT. 1.0_DP)  scxzf = 1.465e-11_DP * xtilde * (ansq*z*z/er)
      RETURN

      END FUNCTION scxzf


       FUNCTION scxif (i, er)
!------------------------------------------------------------------------
! --  estimated cross sections for charge exchange
!------------------------------------------------------------------------

      USE hexnb_data,                              ONLY : kdene,kdeni,kdenz,   &
                                                          ksvi,ksvz,ksve,krad, &
                                                          ngh,ngl                 ! b1
      IMPLICIT NONE 

      REAL(DP) er,aer,siif1,scxif,xtilde,ansq

      INTEGER(I4B) i


 
      IF (i .GE. 2)  go to 200

! --- 1s

      aer = LOG (er)

!     charge exchange
!       below 1.0e6 eV, use Riviere's fit
!       above 1.0e6 eV, use power-law extrapolation.

      IF (er .GT. 1.0e6_DP)  go to 105
      siif1 = (0.6937e-14_DP) * (1.0_DP - 0.06732_DP   *aer)**2 &
                          / (1.0_DP + 0.1112e-14_DP*er  **3.3)
      go to 110
  105 siif1 = (4.8363e-22_DP)/((1.0e-6_DP)*er)**3.3
  110 scxif = siif1
      RETURN

  200 CONTINUE
      ansq   = (nfhx(i))**2
      xtilde = ansq*er / 9260.0_DP
      scxif  = zeroc
      IF (xtilde .LT. 1.0_DP)  scxif = 1.465e-11_DP * xtilde * (ansq/er)
      RETURN

      END FUNCTION scxif



       FUNCTION seef (i, j, er)
!------------------------------------------------------------------------
! --- cross sections for excitation from i to j due to electron impact
!------------------------------------------------------------------------

         USE hexnb_data,                            ONLY : ms,mc,        & ! param
                                                           kdene,kdeni,  &
                                                           kdenz,ksvi,   &
                                                           ksvz,ksve,    &
                                                           krad,ngh,ngl, &  ! b1
                                                           nouthx,istart,& 
                                                           ihxbug,       &  ! b2
                                                           en,dg,ae,be,  &
                                                           de1,de2,ge1,  &
                                                           ge2              ! b4
         

      IMPLICIT NONE
      
       REAL(DP) ryd,er,aer,emin12,emax12,emin13,emax13,seef,zero
       REAL(DP)ae13(9),ae12(8)
       INTEGER(I4B)i,j,ni,nj,id

! --- cross sections for excitation from i to j due to electron impact

      DATA ryd /13.6_DP/

! --- data for 1s to 2s (e impact):

      DATA emin12/10.2_DP/, emax12/1140.5_DP/
    
      DATA ae12/                                                                  &
                -1.5138513657e-17_DP, 1.8231019028e-16_DP,-1.9274654811e-16_DP,   &
                8.9291530932e-17_DP,-2.2108041618e-17_DP, 3.0448025264e-18_DP,    &
                -2.2039298214e-19_DP, 6.5491238735e-21_DP/

! --- data for 1s to 2p (e impact)

      DATA emin13/10.2_DP/, emax13/747.4_DP/
    
      DATA ae13/                                                               &
               -2.1078372920e-14_DP, 3.7548968459e-14_DP,-2.8788463637e-14_DP, &
               1.2394689683e-14_DP,-3.2696005067e-15_DP, 5.4068250737e-16_DP,  &
              -5.4759059881e-17_DP, 3.1084887079e-18_DP,-7.5815695055e-20_DP/

      zero = zeroc
      IF (i .GE. j)  go to 300
      seef = 0.0
      aer  = LOG (er)
      IF ((i .GE. 3) .OR. (j .GE. 4))  go to 250
      go to (210, 240), i

  210 go to (300, 220, 230), j


! --- 1s to 2s:

  220 IF (er .LE. emin12)  RETURN
      IF (er .LE. 11.1  )  go to 221
      IF (er .GE. emax12)  go to 222
      seef = polyf (ae12, 8 ,aer)
      go to 223

! --- linear interpolation for region just above threshold:

  221 seef = ((er - 10.2) / 0.9) * 1.6297e-17
      go to 223

! --- asymptotic:

  222 seef = (5.312e-16) / er
  223 RETURN

! --- 1s to 2p:

  230 IF (er .LE. emin13)  RETURN
      IF (er .GE. emax13)  go to 232
      seef = polyf (ae13, 9, aer)
      go to 233

! --- asymptotic:

  232 seef = (2.65571e-15/er)*(aer-2.120681)
  233 RETURN

! --- 2s to 2p:

  240 seef = (8.617e-15/er)*LOG (1.14e4*er)
      RETURN

  250 ni = nfhx(i)
      nj = j - 1

! --- from Vriens & Smeets (eq.(14)):

      seef = (1.75947e-16*ryd)*(ae(ni,nj)*LOG (0.5 * er/ryd+de1(ni,nj)) &
          + be(ni,nj))/(er+ge1(ni,nj))
      seef = MAX (seef, zero)
      id   = 2
      IF (ni .EQ. 1)  id = 1
      RETURN

  300 ihxbug = 6
      IF (nouthx .GT. 0 .AND. myid == master) THEN
        WRITE (nouthx, 3939)  i, j
 3939   FORMAT (' ERROR in SEEF: i4 and j are in error i4 = ', i3 / &
                '                                       j = ', i3)
      END IF
      RETURN

      END FUNCTION seef



      FUNCTION sezf (i, j, z, er)
!-----------------------------------------------------------------------
! -- cross sections for excitation from i to j due to ion impact
!-----------------------------------------------------------------------

      USE hexnb_data,                       ONLY : kdene,kdeni,kdenz, &
                                                   ksvi,ksvz,ksve,    &
                                                   krad,ngh,ngl,      &  ! b1
                                                   nouthx,istart,     &
                                                   ihxbug                ! b2
      IMPLICIT  NONE



      REAL(DP)  min12,emax12,ai12(11),emin12,emin13,emax13,        &
                ai13(12),zero,z,er,sezf,sezfsave,zz1,aer
      INTEGER(I4B) i,j,ni,id,nj
! --- data for 1s to 2s (p impact):

      DATA      emin12/1.0e3_DP/, emax12/1.0e6_DP/
      DATA      ai12/                                                    &
       -3.3812171941e+05_DP, 3.5976726159e+05_DP,-1.7037301759e+05_DP,   &
        4.7282644608e+04_DP,-8.5167211591e+03_DP, 1.0405267504e+03_DP,   &
       -8.7343818297e+01_DP, 4.9754307346e+00_DP,-1.8412279745e-01_DP,   &
        3.9984716057e-03_DP,-3.8707719188e-05_DP/

! --- data for 1s to 2p (p impact):

      DATA      emin13/1.0e3_DP/, emax13/1.0e6_DP/
      DATA      ai13/                                                    &
       -1.3490069287e+06_DP, 1.4573274888e+06_DP,-7.0815407013e+05_DP,   &
        2.0432442357e+05_DP,-3.8900004212e+04_DP, 5.1319758650e+03_DP,   &
       -4.7884757131e+02_DP, 3.1608287540e+01_DP,-1.4469255104e+00_DP,   &
        4.3759824250e-02_DP,-7.8716881911e-04_DP, 6.3824434435e-06_DP/

      zero = zeroc
      IF (i .GE. j) THEN
        ihxbug = 4
        IF (nouthx .GT. izero .AND. myid == master) THEN
          WRITE  (nouthx, 100)  i, j
  100     FORMAT (' ERROR in SEZF: i4 = ', i3, '    j = ', i3)
        END IF
        sezf = zeroc
        RETURN
      END IF

      IF (z .GT. 1.01_DP)  go to 300
      aer = LOG (er)
      IF ((i .GE. 3) .OR. (j .GE. 4))  go to 250
      go to (210, 240), i

  210 go to (400, 220, 230), j

! --- 1s to 2s:

  220 IF (er .GE. emin12)  go to 221
      sezf = 2.8875e-18*(er/emin12)**0.7093
      go to 223

  221 IF (er .GE. emax12)  go to 222
      sezf = EXP (polyf(ai12,11,aer))
      go to 223

  222 sezf     = 1.9564e-12 / er
  223 sezfsave = sezf
      RETURN

! --- 1s to 2p:

  230 IF (er .GE. emin13)  go to 231
      sezf = 1.4053e-17*(er/emin13)**0.95695
      go to 233

  231 IF (er .GE. emax13)  go to 232
      sezf = EXP (polyf(ai13,12,aer))
      go to 233

  232 sezf     = (6.7085e-13/er) * LOG (5701.79*er)
  233 sezfsave = sezf
      RETURN

! --- 2s to 2p:

  240 sezf     = (1.584e-10/er) * LOG (0.62*er)
      sezf     = MAX (sezf, zero)
      sezfsave = sezf
      RETURN

  250 ni       = nfhx(i)
      id       = 2
      IF (ni .EQ. 1)  id = 1
      nj       = j-1
      zz1      = 1.0
      sezf     = sdaccf(ni,nj,zz1,er)
      sezfsave = sezf
      RETURN

! --- impurity scattering:

  300 IF (i .EQ. 1 .AND. j .EQ. 2)  go to 312
      IF (i .EQ. 1 .AND. j .EQ. 3)  go to 313
      sezf     = 0.0
      sezfsave = sezf
      IF ((i .EQ. 2) .AND. (j .EQ. 3))  RETURN
      ni       = nfhx(i)
      nj       = j - 1
      id       = 2
      IF (ni .EQ. 1)  id = 1
      sezf     = sdaccf(ni,nj,z,er)
      sezfsave = sezf
      RETURN

  312 sezf     = 0.25 * sdaccf(1,2,z,er)
      sezfsave = sezf
      RETURN

  313 sezf     = 0.75 * sdaccf(1,2,z,er)
      sezfsave = sezf
      RETURN

  400 ihxbug = 5
      IF (nouthx .GT. 0 .AND. myid == master) THEN
        WRITE (nouthx, '(/ a /)')  ' ERROR in SEZF: state "j" is 0'
      END IF
      sezf = zeroc
      RETURN

      END FUNCTION sezf



      FUNCTION sizf (i, z, er)
!-----------------------------------------------------------------------
! -- cross section for electron loss via collision with impurity.
!-----------------------------------------------------------------------

      USE hexnb_data,                        ONLY :  kdene,kdeni,kdenz,  &
                                                     ksvi,ksvz,ksve,krad,&
                                                     ngh,ngl               ! b1
      IMPLICIT NONE

      INTEGER(I4B)  ni,i,id,ansq,sizf
      REAL(DP) z,er
!
      ni = nfhx(i)
      id = 2
      IF (ni .EQ. 1)  id = 1
!
      ansq = (FLOAT (ni))**2
      sizf = (1.465e-11_DP)*(ansq*z*z/er)*(1.0_DP-EXP (-ansq*er/(31250.0_DP*z)))
      RETURN
!
      END FUNCTION sizf


      FUNCTION sief (i, er)
!----------------------------------------------------------------------
! -- cross sections for ionization due to e impact
!----------------------------------------------------------------------

      USE hexnb_data,                           ONLY : kdene,kdeni, &
                                                       kdenz,ksvi,  &
                                                       ksvz,ksve,   &
                                                       krad,ngh,    &
                                                       ngl            ! b1
      IMPLICIT NONE

      REAL(DP)     ryd,sief,amin1,amax1,emin2,emax2,          &
                   sbea,ansq,ery,aer,emin1,emax1,er
      REAL(DP)     ae1(7),ae2(9)

      INTEGER(I4B) i
      DATA         ryd/13.6_DP/


! --- data for e-impact ionization of 1s

      DATA      emin1/13.723_DP/, emax1/620.0_DP/

      DATA      ae1/                                                          &
               2.3036652148e-15_DP,-3.4298197580e-15_DP, 1.9465132340e-15_DP, &
              -5.4457519508e-16_DP, 8.0929651995e-17_DP,-6.1527279210e-18_DP, &
               1.8889193736e-19_DP/

! --- data for e-impact ionization of 2s

      DATA emin2/3.4_DP/, emax2/386.0_DP/

      DATA      ae2/                                                           &
                8.4162370468e-15_DP,-2.3545910382e-14_DP, 2.4728342689e-14_DP, &
               -1.2718639429e-14_DP, 3.6382722533e-15_DP,-6.0044378263e-16_DP, &
                5.5114285405e-17_DP,-2.4253587346e-18_DP, 3.0573876445e-20_DP/
!
! --- bea cross section (ansq = n**2, ery = e/ryd):

      sbea(ansq,ery) =  (3.519e-16_DP)*(5.0_DP * ansq/3.0_DP -1.0_DP/ery-2.0_DP/ &
                        (3.0_DP*ansq*ery*ery)) /(ery+3.25_DP/ansq)

      sief = zeroc
      IF (i .GE. 3)  go to 130
      aer  = LOG (er)
      go to (110, 120), i

! --- 1s

  110 IF (er .LE. emin1)  RETURN
      IF (er .GE. emax1)  go to 112
      sief = polyf(ae1,7,aer)
      go to 113
  112 sief = (1.3563e-15_DP/er)*(aer+1.823647_DP)
  113 RETURN

! --- 2s

  120 IF (er .LE. emin2)  RETURN
      IF (er .GE. emax2)  go to 122
      sief = polyf(ae2,9,aer)
      go to 123
  122 sief = (8.195137e-15_DP/er) * (aer-0.9445204_DP)
  123 RETURN

! --- 2p and higher

  130 ansq = (FLOAT (i-1))**2
      IF (er .LE. ryd/ansq)  RETURN
      sief = sbea(ansq,er/ryd)
      RETURN

      END FUNCTION sief


      SUBROUTINE sigfit (ebeam, dene, te, denz, sig)
! ----------------------------------------------------------------------
!.....author:
!       Charles D. Boley
!       Lawrence Livermore National Laboratory
!       L-574
!       Livermore, CA 94550
!       (415) 423-7365
!
!     version of 3/11/92
!
!.....This subroutine computes the neutral beam stopping cross section,
!     for a hydrogenic beam, as a function of beam energy, electron
!     density, electron temperature, and impurity density.  Four
!     impurities (He, C, O, Fe) are allowed.
!
!     The code employs a fit based on results of the author's beam
!     penetration code.  The latter code contains the detailed atomic
!     cross sections recommended at the IAEA Specialists' Meeting on
!     Required Atomic Data Base for Neutral Beam Penetration in Large
!     Tokamaks (Vienna, 4/89).  The fit applies for beam energies from
!     10 keV/amu to 10 MeV/amu and for all densities and electron
!     temperatures of interest.
!
!     The fit is independent of the ion temperature and hence does not
!     distinguish among hydrogen ion isotopes in the plasma.  (The actual
!     cross sections were evaluated at Ti = 15 keV.)  The ion temperature
!     has little effect provided that it is small compared to the beam
!     energy per amu.  (This condition normally is a triviality, except
!     perhaps for the 1/3 energy component of a deuterium beam.)  The
!     following table gives an idea of the variation with ion temperature
!     at low beam energy, for injection into a pure H plasma with
!     dene=1.e14 and Te=10 keV (energies in keV; cross sections in this
!     table in units of 10**-16 cm**2):
!
!     ebeam  sig(Ti=1)  sig(Ti=15)  sig(Ti=30)  sig(fit)
!       10     11.43       9.96        8.91      11.18
!       30      6.06       5.12        4.71       5.35
!       50      3.57       3.56        3.39       3.75
!      100      2.09       2.11        2.09       2.28
!
!     The fit was evaluated with B(perpendicular component) = 5 T.
!     Variations of the magnetic field have an insignificant effect on
!     the stopping cross section.
!
!     Accuracy of fit: rms error about 2.5%
!                      max error about 12%
!
!     Input parameters --
!       ebeam: beam energy per amu (keV/amu) -- 10. .le. ebeam .le. 1.e4
!       dene:  electron density (cm**-3) -- 1.e12 .le. dene .le. 1.e15
!       te:    electron temperature (keV) -- te .ge. 1.
!       denz:  impurity densities (cm**-3) -- array of dimension 4
!         denz(1): He
!         denz(2): C
!         denz(3): O
!         denz(4): Fe
!         Note: denz(i)=0. is permissible
!     Output parameter --
!       sig:   beam stopping cross section (cm**2)
!
!     Modified by John Mandrekas (GIT) for the SUPERCODE
!     New version with fits valid from lower beam energies (10 keV)
!     03/27/92, John Mandrekas, GIT
!------------------------------------------------------------------------------------
      USE hexnb_data,                            ONLY : mzz, mt1, mt2, mt3, & ! param
                                                        A1,A2,nth1,nth2,    &
                                                        nth3,ntz1,ntz2,     &
                                                        ntz3
      IMPLICIT NONE
!

      INTEGER(I4B) i, iinit
!
!      INTEGER(I4B)       nth1, nth2, nth3, ntz1(mzz), ntz2(mzz), ntz3(mzz)
!      REAL(DP)           A1(mt1,mt2,mt3), A2(mt1,mt2,mt3,mzz)
!      COMMON /cfit/ A1, A2, nth1, nth2, nth3, ntz1, ntz2, ntz3
!
      REAL(DP) s1, ebeam, te, sig, sum, dene, &
             denz(mzz), s2(mzz), az(mzz), xx(3)
      DATA az/2.,6.,8.,26./
      DATA iinit/0/



!
      IF (iinit.EQ.0) THEN
        iinit = 1
        CALL initfit
      END IF
!
      xx(1) = LOG (ebeam)
      xx(2) = LOG (dene/1.e13)
      xx(3) = LOG (te)
      s1 = poly3f(A1,mt1,mt2,mt3,xx(1),nth1,xx(2),nth2,xx(3),nth3)
      DO i = 1, 4
         IF (denz(i).GT.0.) THEN
            s2(i)=poly3f(A2(1,1,1,i),mt1,mt2,mt3,xx(1),ntz1(i), &
                  xx(2),ntz2(i),xx(3),ntz3(i))
         ELSE
           s2(i)=0.
         END IF
      END DO
!
      sum=0.0
      DO i = 1, 4
         sum = sum + (denz(i)*az(i)*(az(i)-1.)/dene)*s2(i)
      END DO
!
      sig = (1.e-16/ebeam)*EXP(s1)*(1.+sum)
      RETURN
!
      END SUBROUTINE sigfit 

       FUNCTION siif (i, er)
!---------------------------------------------------------------------------
! --  cross sections for charge exchange and ion-impact ionization
!---------------------------------------------------------------------------

      USE hexnb_data,                                 ONLY : kdene,kdeni, & 
                                                             kdenz,ksvi,  &
                                                             ksvz,ksve,   &
                                                             krad,ngh,    &
                                                             ngl              ! b1
      IMPLICIT NONE
!

!
      REAL(DP)       aer,er,siif,ansq,siif1,siif2

      INTEGER(I4B) i 
!
! --- data for ion-impact ionization of h(1s)
!
      REAL(DP)  ai1(6)
      DATA      ai1/ -4.410635e+02_DP,  1.846170e+02_DP, -3.429509e+01_DP, &
                      3.217263e+00_DP, -1.512004e-01_DP,  2.825854e-03_DP/
!
! --- data for ion-impact ionization of h(1s) (Freeman & Jones)
!
      IF (i .GE. 2)  go to 200
!
! --- 1s:
!
      aer = LOG (er)
!
!     charge exchange:
!       below 1.0e6 eV, use Riviere's fit
!       above 1.0e6 eV, use power-law extrapolation
!
      IF (er .GT. 1.0e6_DP)  go to 105
      siif1 = (0.6937e-14_DP) * (1.0_DP - 0.06732_DP    * aer)**2 &
                           / (1.0_DP + 0.1112e-14_DP * er**3.3)
      go to 110
  105 siif1 = (4.8363e-22_DP) / ((1.0e-6_DP)*er)**3.3
!
!     ion-impact ionization.
!     at low energies or high energies, use Riviere's fits
!     at intermediate energies, use fit to rkj's curve.
!
!     Freeman and Jones option
!
  110 IF (er .LT.    807.4_DP)  go to 113
      IF (er .GT. 154400.0_DP)  go to 114
      siif2 = EXP (polyf(ai1,6,aer))
      go to 120
!
  113 siif2 = EXP (-80.206_DP +8.156_DP *aer-0.3784_DP *aer*aer)
      go to 120
  114 siif2 = (1.56346e-12_DP/er)*(aer-1.792160_DP)
!
  120 siif  = siif1 + siif2
      RETURN
!
  200 ansq  = (nfhx(i))**2
      siif  = (1.465e-11_DP)*(ansq/er)*(1.0_DP-EXP (-ansq*er/9260.0_DP))
      RETURN
!
      END FUNCTION siif

      SUBROUTINE solveq (ndim, n, a, b, ipvt)
!-----------------------------------------------------------------------
! -- used in conjunctin with decomp
!
!   solution of linear system, a*x = b .
!   do not use if decomp has detected singularity.
!
!   input..
!
!     ndim = declared row dimension of array containing a .
!
!     n = order of matrix.
!
!     a = triangularized matrix obtained from decomp .
!
!     b = right hand side vector.
!
!     ipvt = pivot vector obtained from decomp .
!
!   output..
!
!     b = solution vector, x .
!
!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(I4B)  ndim, n, ipvt(n), kb, km1, nm1, kp1, i, k, m
      REAL(DP)   a(ndim,n), b(n), t
!
!     forward elimination
!
      IF (n .EQ. 1)  go to 50
      nm1 = n-1
      DO 20 k=1, nm1
         kp1 = k+1
         m = ipvt(k)
         t = b(m)
         b(m) = b(k)
         b(k) = t
         DO 10 i=kp1, n
             b(i) = b(i) + a(i,k)*t
   10    CONTINUE
   20 CONTINUE
!
!     back substitution
!
      DO 40 kb=1,nm1
         km1 = n-kb
         k = km1+1
         b(k) = b(k)/a(k,k)
         t = -b(k)
         DO 30 i=1, km1
             b(i) = b(i) + a(i,k)*t
   30    CONTINUE
   40 CONTINUE
   50 b(1) = b(1)/a(1,1)
      RETURN
!
      END SUBROUTINE solveq




     SUBROUTINE wrap_xboley (atw,ebkev,ibion,mb,nebin,ebfac, &
                             debin,nion,sgxnmi,atw_beam)
!
! ----------------------------------------------------------------------
!
!  OUTPUT
!  debin      energy bin width if nebin .ne. 0
!
! introduce a wrapper that will allow inclusion of onetwo include
! files without affecting xboley subroutine
! ----------------------------------------------------- HSJ ------------
!
       USE nf_param,                                 ONLY : kj,ke,kb,kcmp1,maxp,nap,kion
       USE  zonal_data,                              ONLY : mfm1,zne,zni,zte,zzi

      IMPLICIT NONE

!     argument list:
      REAL(DP) ebfac,debin,atw_beam
      REAL(DP) ebkev(*),sgxnmi(ke,*),atw(*)
      !zne(*),zni(kz,*),zte(*),zzi(kz,*),
      INTEGER(I4B) ibion,mb,nebin,nion

!     local storage:
!      REAL(DP) sgxnd(kj,ke,kb) ! this is wrong in onetwo 
      REAL(DP) sgxnd(mfm1,ke,kb) 
      ! dnz(kz,4)

      ! dnz(kz,4)
      REAL(DP),ALLOCATABLE,DIMENSION(:,:)        :: dnz
      REAL(DP) ebin(1),ebmax                  ! EBIN(1) FOR XBOLEY ROUTINE
      INTEGER(I4B)  ie,ib,i,j,kbb,mb_1,je,jemax,kz

      kz = mfm1
      IF(.NOT. ALLOCATED(dnz))ALLOCATE(dnz(mfm1,4))
      dnz(:,:) = 0.0_DP

!
      IF (nebin .EQ. 0) THEN   ! no plasma rotation case
          CALL xboley (atw, zzi, ebkev, ibion, mb, mfm1, zne, zte, zni, &
                       dnz, nion, sgxnd, sgxnmi, atw_beam, kz, nw, nh, &
                       maxp, nap, kion, ke, kb)
          DO ie=1,ke
             DO ib=1,mb
                DO j=1,mfm1
                   DO i=1,kcmp1
                      IF (i .EQ. 4) THEN
                          sgxn(i,j,3*(ib-1)+ie,1)=sgxnd(j,ie,ib)
                      ELSE
                          sgxn(i,j,3*(ib-1)+ie,1)=0.0_DP
                      END IF
                    END DO
                END DO
             END DO
           END DO
       ELSE
!
!        include toroidal rotation; use existing method of energy bins
!
!            get the energy bin structure
!
             IF (nebin .NE. 0) THEN
                ebmax = 0.0
                DO ib=1,mb ! use largest beam energy to define the bins
                  ebmax = MAX (ebmax, ebkev(ib))
                END DO
                ebmax = ebmax / atw_beam
                debin = ebmax * ebfac / FLOAT (nebin)
                jemax=nebin
              END IF
!
! have to break out individual beams because of energy bin structure
!
             DO ib=1,mb
                IF (nebin .EQ. 0) THEN ! use to eliminate nebin=0 branch
                  debin = ebkev(ib)
                  jemax = 1
                END IF
                DO je=1,jemax
                   ebin(1) = FLOAT (je) * debin * atw_beam
                   mb_1 = 1    ! force xboley to work on 1 beam at atime
                   kbb  = 1

                   CALL xboley (atw, zzi, ebin, ibion, mb_1, mfm1, zne, &
                                zte, zni, dnz, nion, sgxnd, sgxnmi, &
                                atw_beam, kz, nw, nh, maxp, nap, kion, &
                                ke, kbb)
                   DO ie=1,ke
                     DO j=1,mfm1
                       DO i=1,kcmp1
                         IF (i .EQ. 4) THEN
                           sgxn(i,j,3*(ib-1)+ie,je)=sgxnd(j,ie,ib)
                         ELSE
                           sgxn(i,j,3*(ib-1)+ie,1)=0.0_DP
                         END IF
                       END DO
                     END DO
                   END DO
                END DO
             END DO
      END IF

      RETURN

!
      END       SUBROUTINE wrap_xboley

      SUBROUTINE initfit
! -------------------------------------------------------------------------------
! --- New version, with different fits valid from lower beam energies
! --- Received from C. Boley, LLNL
! --- Modified for use by the SuperCode, John Mandrekas, 03/27/92
!-------------------------------------------------------------------------------

      USE hexnb_data,                            ONLY : mzz, mt1, mt2, mt3, & ! param
                                                        A1,A2,nth1,nth2,    &
                                                        nth3,ntz1,ntz2,     &
                                                        ntz3
      IMPLICIT NONE
!
!      INTEGER mz, mt1, mt2, mt3
!      PARAMETER(mz=4, mt1=4, mt2=3, mt3=2)
      INTEGER i
!
!      INTEGER       nth1, nth2, nth3, ntz1(mz), ntz2(mz), ntz3(mz)
!      REAL*8        A1(mt1,mt2,mt3), A2(mt1,mt2,mt3,mz)
!      COMMON /cfit/ A1, A2, nth1, nth2, nth3, ntz1, ntz2, ntz3
!
!..   pure plasma
!
      nth1 = 3
      nth2 = 3
      nth3 = 2
      A1(1,1,1)= 3.95e+00
      A1(1,1,2)= 1.60e-02
      A1(1,2,1)=-3.84e-02
      A1(1,2,2)=-5.98e-03
      A1(1,3,1)=-3.10e-03
      A1(1,3,2)=-1.09e-03
      A1(2,1,1)= 3.67e-01
      A1(2,1,2)=-2.15e-02
      A1(2,2,1)= 3.07e-02
      A1(2,2,2)= 1.78e-03
      A1(2,3,1)= 3.16e-03
      A1(2,3,2)= 3.47e-04
      A1(3,1,1)=-9.95e-03
      A1(3,1,2)= 6.19e-04
      A1(3,2,1)=-2.36e-03
      A1(3,2,2)=-1.67e-04
      A1(3,3,1)=-1.31e-04
      A1(3,3,2)=-2.28e-05
!
      DO i = 1, 4
         ntz1(i) = 4
         ntz2(i) = 2
         ntz3(i) = 2
      END DO
!
!..   He
!
      i=1
      A2(1,1,1,i)=-1.76e+00
      A2(1,1,2,i)=-2.90e-01
      A2(1,2,1,i)= 5.43e-02
      A2(1,2,2,i)= 1.04e-02
      A2(2,1,1,i)= 7.49e-01
      A2(2,1,2,i)= 1.57e-01
      A2(2,2,1,i)=-3.74e-02
      A2(2,2,2,i)=-6.70e-03
      A2(3,1,1,i)=-7.28e-02
      A2(3,1,2,i)=-2.45e-02
      A2(3,2,1,i)= 6.11e-03
      A2(3,2,2,i)= 1.14e-03
      A2(4,1,1,i)= 2.05e-03
      A2(4,1,2,i)= 1.34e-03
      A2(4,2,1,i)=-2.82e-04
      A2(4,2,2,i)=-5.42e-05
!
!..   C
!
      i=2
      A2(1,1,1,i)=-1.89e-01
      A2(1,1,2,i)=-3.22e-02
      A2(1,2,1,i)= 5.43e-02
      A2(1,2,2,i)= 3.83e-03
      A2(2,1,1,i)=-1.41e-03
      A2(2,1,2,i)= 8.98e-03
      A2(2,2,1,i)=-3.34e-02
      A2(2,2,2,i)=-1.97e-03
      A2(3,1,1,i)= 3.10e-02
      A2(3,1,2,i)= 7.43e-04
      A2(3,2,1,i)= 5.08e-03
      A2(3,2,2,i)= 1.96e-04
      A2(4,1,1,i)=-2.54e-03
      A2(4,1,2,i)=-4.21e-05
      A2(4,2,1,i)=-2.20e-04
      A2(4,2,2,i)= 8.32e-07
!
!..   O
!
      i=3
      A2(1,1,1,i)=-1.07e-01
      A2(1,1,2,i)=-3.36e-02
      A2(1,2,1,i)= 4.90e-02
      A2(1,2,2,i)= 3.77e-03
      A2(2,1,1,i)=-3.72e-02
      A2(2,1,2,i)= 1.11e-02
      A2(2,2,1,i)=-2.89e-02
      A2(2,2,2,i)=-1.84e-03
      A2(3,1,1,i)= 3.34e-02
      A2(3,1,2,i)= 7.20e-06
      A2(3,2,1,i)= 4.14e-03
      A2(3,2,2,i)= 1.64e-04
      A2(4,1,1,i)=-2.50e-03
      A2(4,1,2,i)= 1.08e-05
      A2(4,2,1,i)=-1.65e-04
      A2(4,2,2,i)= 2.45e-06
!
!..   Fe
!
      i=4
      A2(1,1,1,i)= 5.46e-02
      A2(1,1,2,i)=-3.89e-02
      A2(1,2,1,i)= 1.71e-02
      A2(1,2,2,i)=-7.35e-04
      A2(2,1,1,i)=-8.97e-02
      A2(2,1,2,i)= 2.36e-02
      A2(2,2,1,i)=-7.46e-03
      A2(2,2,2,i)= 9.94e-04
      A2(3,1,1,i)= 2.96e-02
      A2(3,1,2,i)=-4.71e-03
      A2(3,2,1,i)= 2.52e-04
      A2(3,2,2,i)=-3.05e-04
      A2(4,1,1,i)=-1.75e-03
      A2(4,1,2,i)= 3.63e-04
      A2(4,2,1,i)= 4.10e-05
      A2(4,2,2,i)= 2.37e-05
!
      RETURN
!
      END SUBROUTINE initfit


      SUBROUTINE hexnb (istarx, iexcix, ilorenx, mstatx, ncontx, &
                        er0x, tex, tix, numix, amix, denix, &
                        numzx, izx, amzx, izstrx, denzx, bperpx, &
                        ldene, ldeni, ldenz, lsvi, lsvz, lsve, lrad, &
                        lngh, lngl, louthx, lcor, &
                        nsigmav, lambda, hexfrac, ihxerr)
!-------------------------------------------------------------------------------
!
!     version 21
!
! --- author:
!       c. d. boley
!       pppl 1984
! --- modified by: (vax / modular)
!       r. m. weiland
!       pppl july 1985
! --- ref.:
!       c. d. boley, r. k. janev, and d. e. post, Phys. Rev. Letters
!       52, 534 (1984).
! --- calling sequence:
!
!     call hexnb (istart, iexcit, ilorent, mstate, ncont,
!    .            er0, tex, tix, numi, ami, deni,
!    .            numz, iz, amz, izstrp, denz, bperp,
!    .            kdene, kdeni, kdenz, ksvi, ksvz, ksve, krad,
!    .            ngh, ngl, nouthx, ncorin, nsigmav, lambda)
!
! --- input parameters:
!       istart = 1 to initialize rad rates & fine pts       0 otherwise
!       iexcit if =0 then ignore contribution of excited states
!              if =1 then include contributionof excited states.
!       ilorent =0 or 1       whether or not to calculate the maximum
!                        principal quantum number ns of the populated
!                        states as given by the lorentz ionization
!                        limit (1  = => yes)
!       mstate< parameter ms
!               for ilorent = 0, use ns = mstate+1
!       note: the operationalsignificance of ns in the code is to
!             set an upper limit such that any excitations to levels n
!             higher than ns are counted as "ionizations".
!       ncont< parameter mc
!               upper bound to number of continuum states
!       er0    beam energy per amu (ev)
!       te     electron temperature (ev)
!       ti     ion temperature (ev)
!       numi   number of hydrogenic ion species
!       ami    masses of the hydrogenic ions (amu)
!       deni   densities of the hydrogenic ions (cm**-3)
!       numz   number of impurity species
!              if numz gt 0, the coronal data file 'coronb' is read.
!       iz     atomic numbers of the impurities
!       izstrp =0 for coronal equilibrium values of <z> and <zsq>
!              =1 for fully stripped impurities
!       amz    atomic mass number of the impurities [iff izstrp#0]
!       denz   densities of the impurities (cm**-3)
!       bperp  magnetic field perpendicular
!              to beamline (tesla (for ilorent = 1))
!              kdene  =0: discard electron reactions
!                     =1: include electron reactions (default)
!              kdeni  =0: discard ion (hydrogen) reactions
!                     =1: include ion reactions (default)
!              kdenz  =0: discard impurity reactions
!                     =1: include impurity reactions (default)
!              ksvi   =0: simple multiplication instead of full sigma-v
!                         integral for ions
!                     =1: full sigma-v integral for ions
!              ksvz   =0: simple multiplication instead of full sigma-v
!                         integral for impurities
!                     =1: full sigma-v integral for impurities
!              ksve   =0: sigma-v integrals for electron reactions
!                         involving 1s, 2s, and 2p       rates for other
!                         electron reactions from Vriens & Smeets.
!                     =1: full sigma-v integrals for electron reactions
!              krad   =0: discard hydrogen line radiation
!                     =1: include hydrogen line radiation .
!              ngh        order of gauss-hermite integrations for
!                         ion and impurity reactions (default 24)
!              ngl        order of gauss-laguerre integrations for
!                         electron reactions
!                         permissible values: 10,16,20,24
!       nouthx unit number for messages(to unit nouthx if > 0)
!              (also nouthx>0 gives detailed diagnostic output)
!       ncorin unit number for reading coronb ce data file
!
! --- output:
!       nsigmav n<sigma-v> [sec**-1]: all processes included
!       lambda  mean free path (cm)
!       hexfrac fraction of 3rd excited state
!       ihxerr: error return
!         0     no error
!         1     radiation file coronb not found
!         2     D01BBF error in setting up ngh integration arrays
!         3     D01BBF error in setting up ngl integration arrays
!         4     SEZF   error: states "i" and    "j" are the same!
!         5     SEZF   error: state  "j" is 0
!         6     SEEF   error: states "i" and/or "j" are in error
!         7     CIEF   error
!         8     CEEF   error
!         9     EIGRF  error: eigenvalue determination is incorrect
!         10    input parameter error
! --- note: this routine uses IMSL libraries
!
! --- foreign files required (SOME INFORMATION BELOW IS ARCHAIC):
!
!       hx2:[wieland.hex]coronb.:  this file contains all the coronal
!                                  equilibrium data for the various impurities
!       dsk0:[imsl]imslibs/lib -- the IMSL replacements for the NAG codes.
!       nag:    f02aff et. al.
!       imsl:   eigrf  et. al.
!       hx2:[wieland.lib]wielib  for diagnostic matrix solvers
!                                DECOMP and SOLVEQ for ax = b
!       subroutine SECOND -- a CPU timing routine of your choice
!
! --- recommended namelist values:
!
!       kdene   = 1
!       kdeni   = 1
!       kdenz   = 1
!       ksvi    = 0
!       ksvz    = 0
!       ksve    = 0
!       krad    = 1
!       ngh     = 10
!       ngl     = 10
!       iexcit  = 1
!       ilorent = 0
!       mstate  = 4
!       ncont   = 30
!
!       these result in a cpu time of approx. 7 sec for a 25 point
!       radial profile of n*sigma(rho) with reasonably good accuracy.
!--------------------------------------------------------------------------------------------

      USE cpub_dat

      USE hexnb_data,                                   ONLY :   er0, v0, te, ti, ami,    &
                                                                 deni, amz,denz, zcor,    &
                                                                 zsqcor, dene, xfrac,ns,  &
                                                                 nc, numi, numz, iz,      &
                                                                 izstrp,kdene,kdeni,kdenz,&
                                                                 ksvi,ksvz,ksve,krad,ngh, &
                                                                 ngl,nouthx,istart,ihxbug,&  
                                                                 iexcit,ilorent,mstate,   &
                                                                 ncont,mi,mz,ms,mc
      IMPLICIT NONE



!

!     argument list:
      REAL(DP)     amix(mi),denix(mi),denzx(mz), amzx(mz) 

      REAL(DP)     er0x,tex,tix,bperpx, nsigmav, lambda,hexfrac

      INTEGER(I4B) istarx,iexcix,mstax,ncontx,numix,numzx,              &
                   izstrx(mz),ldene,ldeni,ldenz,lsvi,lsvz,lsve,lrad,    &
                   lngh,lngl,louthx,lcor,ihxerr,izx(mz),kz,ilorenx,     &
                   mstatx

!      local storage:
       REAL(DP) bperp,xeig

       INTEGER(I4B) i,ncor,ki,iwatch
!

!


      istart = istarx
!
! --- move vbls in hexnb commons
!
      iexcit  = iexcix
      ilorent = ilorenx
      mstate  = mstatx
      ncont   = ncontx
      er0     = er0x
      te      = tex
      ti      = tix
      numi    = numix
      IF (numi .GT. mi)  go to 30
      DO i=1,mi
        ami (i) = amix(i)
        deni(i) = denix(i)
      END DO
      IF (numz .GT. mz)  go to 30
      numz = numzx
      DO i=1,mz
        iz(i)     = izx(i)
        izstrp(i) = izstrx(i)
        amz(i)    = amzx(i)
        denz(i)   = denzx(i)
      END DO
      bperp   = bperpx
      kdene   = ldene
      kdeni   = ldeni
      kdenz   = ldenz
      ksvi    = lsvi
      ksvz    = lsvz
      ksve    = lsve
      krad    = lrad
      ngh     = lngh
      ngl     = lngl
      nouthx  = louthx
      ncor    = lcor
      ihxbug  = izero
      nsigmav = zeroc
      lambda  = zeroc
      IF (istart .EQ. izero)  go to 5
!
      IF (mstate+1 .GT. ms .OR. ncont+1 .GT. mc)  go to 30
!***  ncor = 31
      CALL hradin (ncor, numz, iz, izstrp, amz)
      IF (ihxbug .GT. 0)  go to 20
      CALL hxinit
!
    5 v0 = 1.3841e6 * SQRT (er0)
      CALL hxradi (te, zcor, zsqcor, numz, iz, izstrp, iwatch)
      dene = 0.0
!
      DO ki=1,numi
        dene = dene + deni(ki)
      END DO
!
      DO kz=1,numz
        dene = dene + zcor(kz) * denz(kz)
      END DO
!
! --- determine calculational mode: include excitations or not?
!
      go to (21, 22),  iexcit + 1
!
!     no excitations included
!
   21 nc = 1
!***  if (nouthx .gt. 0)  write (nouthx, 1000)
!1000 format (/ ' no excitations')
      go to 23
!
!     include excitations
!
   22 nc = ncont + 1
!***  if (nouthx .gt. 0)  write (nouthx,1004) ncont
!1004 format (/ ' excitations, with continuum at n = ',i3)
!
   23 CALL lorent (v0, bperp, nc, ns)
!***  if (nouthx .gt. 0)  write (nouthx, 1005) ns
!1005 format (/ ' max princ qn= ',i3)
!
      CALL hxsvi
      IF (ihxbug .GT. 0)  go to 20
!
      CALL hxsve
      IF (ihxbug .GT. 0)  go to 20
!
      CALL matri
!
      CALL eigen (xeig)
      IF (ihxbug .GT. 0)  go to 20
!
      lambda  = xeig
      nsigmav = v0 / xeig
      hexfrac = xfrac
!
   20 ihxerr  = ihxbug
      istart  = 0
      RETURN
!
   30 ihxbug = 10
      IF (nouthx .GT. 0 .AND. myid == master)  WRITE (nouthx, 40) numi,numz,mstate,ncont
   40 FORMAT (' **** inconsistency between input vbls and upper', &
              ' limits as defined by parameters:' / &
                3x, 'numi= '  , i4                / &
                3x, 'numz= '  , i4                / &
                3x, 'mstate= ', i4                / &
                3x, 'ncont= ' , i4)
      go to 20
!
      END      SUBROUTINE hexnb


      SUBROUTINE xboley (atw, zzi, ebkev, ibion, mb, mfm1, zne, zte, &
                         zni, dnz, nion, sgxn, sgxnmi, atw_beam, &
                         kz, ki, kj, maxp, nap, kion, ke, kb)
!------------------------------------------------------------------------------------
!
!   This subroutine calculates the relevant cross section quantities
!   -sgxn- and -sgxnmi-, needed by the NFREYA code, using fits that are
!   based on the recent work by Janev, Boley and Post. See routine
!   sigfit for more details.
!   created  on 01-05-90 by J. Mandrekas
!   modified on 01-08-91 by J. Mandrekas to treat He3 as He4 for the
!   ARIES-III calculations
!   modified on 02-20-92 by J. Mandrekas for SUN/486/UNICOS platforms
!-------------------------------------------------------------------------------------
 
      IMPLICIT NONE
!
!     argument list:
      INTEGER(I4B) kion,kz,kb,ke,kj,maxp,nap,mfm1,mb,nion,ki,ibion

      REAL(DP)  atw(kion), zzi(kz,kion), ebkev(kb), zne(kz), zte(kz), &
                zni(kz,kion), sgxn(kz,ke,kb), sgxnmi(ke,kb), denz(4), &
                dnz(kz,4),atw_beam


!     local storage:
      INTEGER(I4B) i,k,jhe3,jhe4,jcarb,joxy,jfe,dhe3,dhe4,iz,ib,ie,j

      REAL(DP) sgxnm,sig,ebeam,dene,te


      DATA      jhe3 /0/, jhe4 /0/, jcarb /0/, joxy /0/, jfe /0/, &
                dhe3 /0.0/, dhe4 /0.0/
!
!     zero-out the impurity densities
!
      DO   i=1,mfm1
        DO k=1,4
          dnz(i,k) = 0.0_DP
        END DO
      END DO
!
! --- recognize impurity species (including He)
!
      DO i = 1, nion
        IF (atw(i) .EQ.  3.0 .AND. zzi(1,i) .EQ.  2.0)  jhe3  = i
        IF (atw(i) .EQ.  4.0 .AND. zzi(1,i) .EQ.  2.0)  jhe4  = i
        IF (atw(i) .EQ. 12.0 .AND. zzi(1,i) .EQ.  6.0)  jcarb = i
        IF (atw(i) .EQ. 16.0 .AND. zzi(1,i) .EQ.  8.0)  joxy  = i
        IF (atw(i) .EQ. 56.0 .AND. zzi(1,i) .EQ. 26.0)  jfe   = i
      END DO
      DO iz = 1, mfm1
        DO k = 1, nion                  ! this loop is bogus HSJ
          IF (jhe3  .NE. 0)  dhe3      = zni(iz,jhe3)
          IF (jhe4  .NE. 0)  dhe4      = zni(iz,jhe4)
                             dnz(iz,1) = dhe3 + dhe4
          IF (jcarb .NE. 0)  dnz(iz,2) = zni(iz,jcarb)
          IF (joxy  .NE. 0)  dnz(iz,3) = zni(iz,joxy)
          IF (jfe   .NE. 0)  dnz(iz,4) = zni(iz,jfe)
        END DO
      END DO
!
      DO iz = 1, mfm1
        denz(1) = dnz(iz,1)
        denz(2) = dnz(iz,2)
        denz(3) = dnz(iz,3)
        denz(4) = dnz(iz,4)
        te   = zte(iz)
        dene = zne(iz)
        DO   ib=1,mb
          DO ie=1,3
            ebeam = ebkev(ib) / (ie*atw_beam)
            CALL sigfit (ebeam, dene, te, denz, sig)
            sgxn(iz,ie,ib) = dene * sig
          END DO
        END DO
      END DO

      DO ib = 1, mb
        IF (ib .LE. 1 .OR. ebkev(ib) .NE. ebkev(1)) THEN
          DO j = 1, 3
            sgxnm = sgxn(1,j,ib)
            DO i = 2, mfm1
              IF (sgxnm .LT. sgxn(i,j,ib))  sgxnm = sgxn(i,j,ib)
            END DO
            sgxnmi(j,ib) = 1.0 / sgxnm
          END DO
        ELSE
          DO j = 1, 3
            sgxnmi(j,ib) = sgxnmi(j,1)
          END DO
        END IF
      END DO
      RETURN
!
      END       SUBROUTINE xboley

    SUBROUTINE distribute_xsct
!----------------------------------------------------------------------
! -- broadcast xsct data to other cpus
!----------------------------------------------------------------------
      

    USE  neutral_beams,        ONLY : de_tk


    IMPLICIT NONE
    INTEGER nsend


#ifdef USEMPI
           CALL MPI_BCAST(de_tk,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
 
           nsend = ke*kb

           CALL MPI_BCAST(hxfrac,nsend,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
 
           CALL MPI_BCAST(sgxnmi,nsend,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)

           nsend = kcmp1*mfm1*kbe*ksge

           CALL MPI_BCAST(sgxn,nsend,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpiierr)
          
#endif    
      RETURN

      END SUBROUTINE distribute_xsct



   END MODULE xsct
