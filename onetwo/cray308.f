      subroutine out
c
      USE param
      USE fusion
      USE io 
      USE ions
      USE neut
      USE nub
      USE nub2
      USE solcon
      USE soln
      USE mhdpar 
      USE extra
      USE rf
      USE yoka
      USE numbrs
      USE mesh
      USE sourc                          ! totrf
      USE machin
      USE tfact
      USE geom
      USE flags
      USE tordlrot
      USE constnts
      USE tcoef
      USE bd_condtn
      USE mixcom
      USE rhog
      USE flx
      USE flxav
      USE neo2d
      USE gpsi
      USE zeffcom
      USE pelcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray308.f,v 1.71 2014/05/08 17:48:12 stjohn Exp $"/
c
c --- this subroutine prints out results at selected times
c

      include 'co2.i'

      include 'imsl.i'


      include 'sxrcom.i'

c
      dimension  fluxc(kk), ena(kion), enx(kj), w3neo(kj), xkangsav(kj)
      character  blank*1, jp*2, aions*6, aneuts*6
c
c add bootstrap current totals
c
      dimension  tjbni(kion), tjbti(kion)
c
c set ihead = 1 so heading will be printed the next time that info is called
c
*%%%% write (6, *) '%%%% entering subroutine OUT'
      ihead = 1
c
c ----------------------------------------------------------------------
c                        quantities at time point n-1/2
c ----------------------------------------------------------------------
c
c calculate time and time point number
c

 

      t     = n
      timet = time
      if (n .eq. 0)  go to 110
      t     =     t - (1.0-theta)
      timet = timet - (1.0-theta)*dtt
  110 if (itimav .eq. 1)      t = n
      if (itimav .eq. 1)  timet = time
c
c  calculate fluxes at center
c
*%%%% write (6, *) '%%%% in subroutine OUT, calling RHOMSH'
      call rhomsh (timet)
      do 120 k=1,nk-1 -iangrot
  120 fluxc(k) = 0.0
      if (iangrot .eq. 1)  fluxc(nk) = 0.0

c
c print out fluxes
c
      nwrt = nion + 2
      call header (nout, timet, t)
      if (nion .eq. 1)  write (nout, 1071) (i,i=1,nion),
     .                 (itranflag(i),i = 1,nwrt),itranflag(nwrt+2)
      if (nion .eq. 2)  write (nout, 1072) (i,i=1,nion),
     .                 (itranflag(i),i = 1,nwrt),itranflag(nwrt+2)
      if (nion .eq. 3)  write (nout, 1073) (i,i=1,nion),
     .                 (itranflag(i),i = 1,nwrt),itranflag(nwrt+2)
      if (nion .eq. 4)  write (nout, 1074) (i,i=1,nion),
     .                 (itranflag(i),i = 1,nwrt),itranflag(nwrt+2)
      if (nion .eq. 5)  write (nout, 1075) (i,i=1,nion),
     .                 (itranflag(i),i = 1,nwrt),itranflag(nwrt+2)
      j  = 1
      jp = '  '
      write (nout, 1080) j, jp,  r(j), (fluxc(i), i=1,nwrt),
     .                                  fluxc(nwrt+2)
      jp = '.5'
c
      do 150 j=1,nj-1
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj-1)  go to 150
        write (nout, 1080)
     .    j, jp, ra(j), (flux(i,j), i = 1,nwrt), flux(nwrt+2,j)
  150 continue
c
      j  = nj
      jp = '  '
      write (nout, 1080) j, jp, r (j),(fluxb(i), i=1,nwrt),fluxb(nwrt+2)
c
c print out components of fluxes
c
      if (jflux .eq. 0)  go to 210
c      write (6, *) '%%%% in subroutine OUT, calling PFLUX'
      call pflux (timet, t)
c
c print out transport coefficients
c
  210 if (jcoef .eq. 0)  go to 220
c       write (6, *) '%%%% in subroutine OUT, calling PCOEF'
      call pcoef (timet, t)



c
c print out components of particle and energy sources
c
  220 if (jsourc .eq. 0)  go to 230
c
c --- save xkangrot (it is converted from mesh center to mesh in psourc
c --- (xkeneo is also converted from mesh center to mesh in psourc)
c
      if (itran(nk) .ne. 0)  call copya (xkangrot, xkangsav, nj)
c      write (6, *) '%%%% in subroutine OUT, calling PSOURC'

      call psourc (timet, t)

c      write (6, *) '%%%% in subroutine OUT, returning from PSOURC'
c
c --- restore xkangrot
c
      if (itran(nk) .ne. 0)  call copya (xkangsav, xkangrot, nj)
c
c print out Fred Marcus output 
c
  230 continue
      if (ifred .eq. 0)  go to 240
*%%%% write (6, *) '%%%% in subroutine OUT, calling FRED'
      call fred

c
c --- print out angular momentum information
c

  240 if (iangrot .ne. 0 )  
     .                 call printang (nout, timet, t, jprt)

c
c print table of torque densities
c
*%%%% write (6, *) '%%%% in subroutine OUT, calling PTORQUE'
      call ptorque (time, t)

c
c print out detailed fusion information
c
      if (jtfus .eq. 0 .or. ifus .eq. 0)  go to 400
*%%%% write (6, *) '%%%% in subroutine OUT, calling PTHERF'
      call ptherf (timet, t)

c
c print particle and energy confinement time information
c
  400 if (enalp(1) .gt. 0.0)  call fustable
c
      nunit = nout
      iunit = nunit
      call header (iunit, timet, t)
      write (iunit, 1155) entot, dentot, stot, taup,
     .                    eetot, deetot, qetot, tauee,
     .                    etot, detot, qtot, taue, volume, entaue,
     .                    sfarea, ec, dec, qcen, tauec, angmtot,
     .                    dangmtot, storquet, tauangt
      write (ncrt , 1155) entot, dentot, stot, taup,
     .                    eetot, deetot, qetot, tauee,
     .                    etot, detot, qtot, taue, volume, entaue,
     .                    sfarea, ec, dec, qcen, tauec, angmtot,
     .                    dangmtot, storquet, tauangt
c
c print out neutral data
c


      if (nneu .eq. 0)  go to 440
      call header (nout, timet, t)
      if (  nneu .eq. 1  )  write (nout, 1090)  namep(in), in
      if (  nneu .eq. 2  )  write (nout, 1090) (namep(i ), i, i=1,2)
      if (widths .ne. 0.0)  write (nout, 1091)  widths
      write (nout, 1092)
c
      do 420 j=1,njs
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj .and. j .ne. njs)  go to 420
        write (nout,1094) j,r(j),tn(j,1),ennw(j,1),ennv(j,1),volsn(j,1),
     .    vneut(j,1), tn(j,2),ennw(j,2),ennv(j,2),volsn(j,2),vneut(j,2)
 420  continue
c
      if (idiagn .eq. 0)  go to 440
      if ( nengn .eq. 0)  go to 440
      m1 = 1
      if (nneu .eq. 1)  m1 = in
      imslmd = 'nuspec'
      call nuspec (enn(1,m1),enn(1,2),rtandn,rhdn,rmajor,
     .             englstn,nengn,f1w,f2w,fvn,
     .             spflux(1,m1),spflux(1,2))
      write (nout, 1096)  rtandn,rhdn,f1w,f2w,fvn
      do i=1,nengn
        write (nout, 1098)  englstn(i), spflux(i,1), spflux(i,2)
      end do
c
c ----------------------------------------------------------------------
c                        quantities at time point n
c ----------------------------------------------------------------------
c
c calculate integration factors and get new rho mesh
c
  440 t = n
c
      call rhomsh(time)
      do j=1,nj
        fact(j) = 4.0 * pi**2 * rmajor * hcap(j)
      end do
      cconst = 1.0
      call trap1 (r, fact, cconst, nj, volume)
      do j=1,nj
        fact(j) = fact(j) / volume
      end do
c
c print out densities for electrons, ions, and neutrals
c
      nunit = nout
      call header (nunit, time, t)
      write (nunit, 1015)  (nameu(i), i=1,nion)
      write (nunit, 1017)
      jlines = 0
c
      do 510 j=1,nj
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 510
        jlines = jlines+1
        write (nunit, 1020) j, r(j), zeff(j), ene(j), (en(j,i),i=1,nion)
  510 continue
c
      do j=1,nj
        enx(j) = ene(j) * zeff(j)
      end do
c
      call trapv (r,enx,fact,nj,zeffa)
      zeffa = zeffa/eneav
      do 525 i=1,nion
  525 call trapv (r,en(1,i),fact,nj,ena(i))
      write (nunit, 1022)  zeffa, eneav, (ena(i),i=1,nion)
      if (codeid .ne. 'onedee')  go to 528
      jlines = jlines + 2
      write (nunit, 1023)  enebar
  528 write (nunit, 1024)  entotn, enitn, snaddt
c
      if (jlines .gt. 15)  call header (nunit, time, t)
      write (nunit, 1016)  namep(1), namep(2)
c
      do 531 j=1,nj
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 531
        write (nunit, 1021)  j,r(j),enbeam(j),enalp(j),enn(j,1),enn(j,2)
  531 continue
c
      call trapv (r, enbeam  , fact, nj, enba )
      call trapv (r, enalp   , fact, nj, enaa )
      call trapv (r, enn(1,1), fact, nj, enn1a)
      call trapv (r, enn(1,2), fact, nj, enn2a)
      write (nunit, 1025)  enba, enaa, enn1a, enn2a
c
c print out density and charge of impurity ions
c
      if (nimp .eq. 0)  go to 600
c
      do i=1,nimp
        k = nprim + i
        call header (nout, time, t)
        write (nout, 1030)  namei(i)
        do 540 j=1,nj
          j1prt  = ((j-1)/jprt)*jprt
          if (j1prt .ne. j-1 .and. j .ne. nj)  go to 540
          rotzsq = SQRT (zsq(j,k))
          write (nout,1040) j, r(j), en(j,k), z(j,k), rotzsq
  540   continue
      end do
c
c print out temperatures, magnetic field, current density, and electric field
c
  600 nunit = nout
      call header (nunit, time, t)
      write (nunit, 1045)
      write (ncrt , 1045)
c
      do 610 j=1,nj
        j1prt = (j-1) / jprt * jprt
****    bp    = 0.0
****    if (r(j) .ne. 0.0)  bp = rbp(j)/(r(j)*fcap(j)*gcap(j)*hcap(j))
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 610
        if (wmix .eq. 0.0)  write (nunit, 1060)  j, r(j), te(j), ti(j),
     .    tn(j,1), tn(j,2), bpol(j), curden(j), etor(j), q(j)
        if (wmix .ne. 0.0)  write (nunit, 1060)  j, r(j), te(j), ti(j),
     .    tn(j,1), tn(j,2), bpol(j), curden(j), etor(j), q(j), psis(j)
        if (   j .eq. 1  )  write (ncrt , 1060)  j, r(j), te(j), ti(j),
     .    tn(j,1), tn(j,2), bpol(j), curden(j), etor(j), q(j), psis(j)
  610 continue
c
c ----------------------------------------------------------------------
c calculate and print out volume-average temperatures
c ----------------------------------------------------------------------
c
      call trapv (r, te, fact, nj, teav)
****  teav = teav / volume
      call trapv (r, ti, fact, nj, tiav)
****  tiav = tiav / volume
      write (nunit, 1062)  teav, tiav
      write (ncrt , 1062)  teav, tiav
      write (nunit, 1064)
c
c ----------------------------------------------------------------------
c print out current density and current drive profiles
c ----------------------------------------------------------------------
c
      icprt = ibeam + ibcur + iboot
      do 630 k=1,krf
         if(irfcur(k) .ne. 0.0)icprt=icprt+1
  630 icprt = icprt+irf(k)
      if (icprt .eq. 0)  go to 650
      write (nunit, 1053)
c
      do 640 j=1,nj
        j1prt = (j-1) / jprt * jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 640
        write (nunit, 1056)  j, r(j), curden(j), curohm(j), curboot(j),
     .              curbi(j), curbe(j), curbet(j), currf(j), curdri(j),
     .              curpar_soln(j)
c
  640 continue
c
      call trapv (r, curden , hcap, nj, totcurden)
      totcurden  = 2.0 * pi * totcurden
      call trapv (r, curpar_soln , bsqncap, nj, totcurpar_dmc)
      totcurpar_dmc  = 2.0 * pi * totcurpar_dmc
      call trapv (r, cur_tor_ps_soln , hcap, nj, totcurpar_ps_dmc)
      totcurpar_ps_dmc = 2.0 * pi * totcurpar_ps_dmc
      totcurpar_dmc  =       totcurpar_dmc  + totcurpar_ps_dmc
c
      if (curtype .eq. 0) then
c
c     NOTE THAT THESE INTEGRALS,WITH DA=2.*PI*HCAP*r ,ARE CORRECT ONLY
c     for current densitites of the form ,< JX *R0/R>,a fact we are
c     neglecting here and elsewhere   HSJ
c
           call trapv (r, curohm , hcap, nj, totohm)
           call trapv (r, curboot, hcap, nj, totboot)
           call trapv (r, curbi  , hcap, nj, totbi)
           call trapv (r, curbe  , hcap, nj, totbe)

           call trapv (r, currf  , hcap, nj, totrf)
           call trapv (r, curdri , hcap, nj, totdri)
           call trapv (r, curbet , hcap, nj, totbet)
c
           totohm  = 2.0 * pi * totohm
           totboot = 2.0 * pi * totboot
           totbet  = 2.0 * pi * totbet
           totbi   = 2.0 * pi * totbi
           totbe   = 2.0 * pi * totbe
           totrf   = 2.0 * pi * totrf
           totdri  = 2.0 * pi * totdri
           totbeam = totdri - totrf

      else ! diamagnetic correction, added 8/25/98 HSJ
c
c          additional pressure gradient driven term added to bootstrap
c
 
           call trapv (r, curohm , bsqncap, nj, totohm_dmc)
           call trapv (r, curboot, bsqncap, nj, totboot_dmc)
           call trapv (r, curbi  , bsqncap, nj, totbi_dmc)
           call trapv (r, curbe  , bsqncap, nj, totbe_dmc)
           call trapv (r, currf  , bsqncap, nj, totrf_dmc)
           call trapv (r, curdri , bsqncap, nj, totdri_dmc)
           call trapv (r, curbet , bsqncap, nj, totbet_dmc)
           totohm  = 2.0 * pi * totohm_dmc
           totboot = 2.0 * pi * totboot_dmc + totcurpar_ps_dmc
           totbet  = 2.0 * pi * totbet_dmc
           totbi   = 2.0 * pi * totbi_dmc
           totbe   = 2.0 * pi * totbe_dmc
           totrf   = 2.0 * pi * totrf_dmc
           totdri  = 2.0 * pi * totdri_dmc
           totbeam = totdri - totrf
c      print *,'tcurden,totohm,totboot,totcurpar_ps_dmc,totbeam,totrf =',
c     .        totcurden,totohm,totboot,totcurpar_ps_dmc,totbeam,totrf
      end if
c
      write (nunit, 1054)  totcur(1), totohm, totboot, totbi, totbe,
     .                     totbet, totrf, totdri
      do l=1,krf
      if (rf_ext_curtot(l) .ne. 0.0)
     .  write (nunit, 1058)  extcurrf_id(l), rf_ext_curtot(l)
 1058  format (' currf contains an externally input contribution' /
     .         ' identified in  inone as ', a /
     .         ' with total current of ', 1pe14.4,' amps')
      end do

c
c      the profile to be printed here is extcurrf_curr, together with
c      other RF current sources.
c
c ----------------------------------------------------------------------
c     Write out current profiles to the "nclout" file
c --------------------------------------------------------------- daniel
c
      inclout = 0
      if (inclout .eq. 1) then
        call nclout (nunncl, time, t, k, r, nj, jprt, jhirsh, curden,
     .               curohm, curboot, curbi, curbe, curbet, currf,
     .               curdri, totcur(1), totohm, totboot, totbi, totbe,
     .               totbet, totrf, totdri)
      end if
c
c ----------------------------------------------------------------------
c     print out several average quantities
c ----------------------------------------------------------------------
c
  650 write (nunit, 1051) totcur(1), voltag, betae, betae0,
     .                    betai, betai0, betab, betab0, betaa, betaa0,
     .                    beta, beta0, betap, xlica
c
c     print out breakdown of bootstrap current from Hirshman
c
      if (jhirsh .eq. 0)  go to 3100
c
c --- convert trapped/passing fraction from half mesh to full mesh
c     for printing
c
      call mescon (ftfc, dr, nj)
      ftfc(1) = 0.0
      if (icprt .eq. 0)  go to 3100
      nunit   = nout
      call header (nunit, time, t)
      if (jhirsh .gt. 30)  write (nunit, 3018) jhirsh
 3018 format (' WARNING: jhirsh =', I5                         /
     .        ' This is a single ion species fluid model.'     /
     .        ' The following breakdown into various gradient' /
     .        ' terms is therefore artificial.')
c
      if (jhirsh .eq. 2 .or. jhirsh .eq. 22 .or. jhirsh .eq. 89)
     .write  (nunit, 3019)
 3019 format (' The fast ion density gradient term is an ad hoc'   /
     .        ' approximation which assumes that fast ions can be' /
     .        ' treated as a thermal species.')
      write (nunit, 3015)  (nameu(i), i=1,nion)
      write (nunit, 3017)
      jlines = 0
c
      do 3510 j=1,nj
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 3510
        jlines = jlines + 1
        write  (nunit, 3020)  j,r(j),ftfc(j),xjbne(j),xjbnf(j),
     .                       (xjbni(j,i),i = 1,nion)
 3020   format (i5, f10.2, 1p10e12.3)
 3510 continue
c
      call trapv (r, xjbne, hcap, nj, tjbne)
      call trapv (r, xjbnf, hcap, nj, tjbnf)
      call trapv (r, xjbte, hcap, nj, tjbte)
      tjbne = 2.0 * pi * tjbne
      tjbnf = 2.0 * pi * tjbnf
      tjbte = 2.0 * pi * tjbte
      do i=1,nion
        call trapv (r,xjbni(1,i),hcap,nj,tjbni(i))
        call trapv (r,xjbti(1,i),hcap,nj,tjbti(i))
        tjbni(i) = 2.0 * pi * tjbni(i)
        tjbti(i) = 2.0 * pi * tjbti(i)
      end do
      write (nunit, 3054)  tjbne,tjbnf,(tjbni(i),i=1,nion)
      write (nunit, 3016) (nameu(i),i=1,nion)
      write (nunit, 3017)
      jlines = 0
      do 3511 j=1,nj
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 3511
        jlines = jlines+1
        write  (nunit, 3021) j,r(j),ftfc(j),xjbte(j),
     .                      (xjbti(j,i),i = 1,nion)
 3021   format (i5, f10.2, 1p9e12.3)
 3511 continue
      write (nunit, 3055)  tjbte, (tjbti(i), i=1,nion)
c
c ----------------------------------------------------------------------
c print out number of neutrons produced by D-D fusion, if deuterium
c is present as a primary ion.  also print out neutron rates due to
c knock-on and beam deuterium-plasma processes.
c ----------------------------------------------------------------------
c
 3100 if (iddfus .eq. 0)  go to 670
      call header (nout, time, t)
      write (nout, 1055)
      ksymp1 = nbeams*3 + 1
      do 660 j=1,nj
        j1prt = ((j-1)/jprt) * jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 660
        if (iddfus .eq. 1)  enddf =      en(j,id )
        if (iddfus .eq. 2)  enddf = fd * en(j,idt)
        if (iddfus .eq. 3)  enddf =      en(j,id )
        if (iddfus .eq. 4)  enddf = 0.0 ! tritium case - no printout yet
        write (nout, 1056) j, r(j), enddf, ti(j), ddfusn(j), ddknck(j),
     .                     ddbeam(j), beam_thermalddn(j,ksymp1),
     .                     beam_thermalddp(j,ksymp1)
  660 continue
      write (nout, 1057) ddnthm,ddknct,ddbmt,ddntot,beam_thermal_ddntot,
     .                   beam_thermal_ddptot
c
c ----------------------------------------------------------------------
c print out co2 and zeff information
c ----------------------------------------------------------------------
c
  670 if (jco2 .eq. 0 .and. jzeff .eq. 0)  go to 710
      if (iddfus .eq. 0)  call header (nout, time, t)
      if (  jco2 .eq. 0)  go to 690
      write (nout, 1500)
      do 680 i=1,nco2
  680 write (nout, 1510)  rtco2(i), patco2(i), denco2(i)
  690 if (jzeff .eq. 0)  go to 710
      write (nout,1520) zeffwl,zeffdl,zefrad
      do 700 i=1,8
        i2 = i + 8
        write (nout, 1530)  i, rtzeff(i ), phzeff(i),
     .                     i2, rtzeff(i2), phzeff(i2)
  700 continue
c
c ----------------------------------------------------------------------
c print summary profiles to qikone file
c ----------------------------------------------------------------------
c
  710 call header (nqik, time, t)
      call header (nout, time, t)
      write (nqik, 1400)
      write (nout, 1400)
c
      aions  = ' ions '
      aneuts = ' neuts'
      write  (nqik, 1410) (nameu(i), aions , i=1,nion),
     .                    (namep(i), aneuts, i=1,2   )
      write  (nout, 1410) (nameu(i), aions , i=1,nion),
     .                    (namep(i), aneuts, i=1,2   )
 1410 format (/ 4x,'j',4x,'r',3x,'r/a',6x,'zeff',7x,'ene',8x,'enbeam',
     .          6x, 'enalp', 2x, 6(4x, a2, a6))
c
      blank = ' '
      write  (nqik, 1420)  (blank, i=1,nion+5)
      write  (nout, 1420)  (blank, i=1,nion+5)
 1420 format (7x, '(cm)', 16x, 10(2x, a1, '(1/cm**3)'))
      do 780 j=1,nj
        j1prt = (j-1) / jprt * jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 780
        write (nqik, 1430) j, r(j), roa(j), zeff(j), ene(j),
     .                     enbeam(j), enalp(j),
     .                    (en(j,i), i=1,nion), enn(j,1), enn(j,2)
        write (nout, 1430) j, r(j), roa(j), zeff(j), ene(j),
     .                     enbeam(j), enalp(j),
     .                    (en(j,i), i=1,nion), enn(j,1), enn(j,2)
  780 continue
      write (nqik, 1440)
      write (nout, 1440)
      do 790 j=1,nj
        j1prt = (j-1) / jprt * jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 790
        wbw       = wbeam(j)*1.60217733e-16
        wbjcm3(j) = wbw
        waw       = walp(j)*1.60217733e-16
        wew       = 1.5*ene(j)*te(j)*1.60217733e-16
        wejcm3(j) = wew
        ensum     = 0.0
        do 795 i=1,nion
  795   ensum     = ensum + en(j,i)
        wiw       = 1.5*ensum*ti(j)*1.60217733e-16
        wijcm3(j) = wiw
        write (nqik, 1450)  j,r(j),roa(j),te(j),ti(j),wew,wiw,wbw,waw
        write (nout, 1450)  j,r(j),roa(j),te(j),ti(j),wew,wiw,wbw,waw
  790 continue
      write (nqik, 1460)
****  write (nout, 1460)
      do 800 j=1,nj
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 800
        et    = etor(j) + esaw(j)
        psirx =   1.0e3 * psir(j)
        write (nqik, 1470) j, r(j), roa(j), q(j), etor(j), et,
     .                     curden(j), curohm(j), curboot(j), curdri(j),
     .                     bpol(j), psirx
****    write (nout, 1470) j, r(j), roa(j), q(j), etor(j), et,
**** .                     curden(j), curohm(j), curboot(j), curdri(j),
**** .                     bpol(j), psirx
  800 continue
      write (nqik, 1471)  totcur(1), totohm, totboot, totdri
****  write (nout, 1471)  totcur(1), totohm, totboot, totdri
c
c ----------------------------------------------------------------------
c print transport summary page
c ----------------------------------------------------------------------
c
****  if (ilastp .eq. 0 .or. iyoka .eq. 0)  go to 810
      if (ilastp .eq. 0)  go to 810    ! print, independent of iyoka HSJ
      call psumry (nout)
      call psumry (nqik)
c
c ----------------------------------------------------------------------
c print tweak information
c ----------------------------------------------------------------------
c
  810 if (ttweak .ne. 0.0 .or.   w2mix .ne. 0.0 .or. w3mix .ne. 0.0
     .                    .or. (ilastp .ne. 0  .and. iyoka .ne. 0) )
     .                    call header (nout, time, t)
      if (ttweak .ne. 0.0 .or.   w2mix .ne. 0.0 .or. w3mix .ne. 0.0
     .                    .or. (ilastp .ne. 0  .and. iyoka .ne. 0) )
     .                    call header (nqik, time, t)
      do 830 nunit=nout,nqik
      if (ttweak .eq. 0.0)  go to 820
      w33 = wneo(3,3)
      if (w3typt .ne. 0.0)  w33 = w3typt
      write (nunit, 2050)  ttweak, w33, zeffc(1), zeffb(1)
      write (nunit, 2052)  fusnin, ddntot, ticin, ti(1), voltin, voltag,
     .                     qcin, q(1)
c
  820 if (w2mix .eq. 0.0 .and. w3mix .eq. 0.0)  go to 830
      write (nunit, 2060) timmix, rsx, rmixx, epste, epsti
      write (nunit, 2062) dtemix, dtecal, s3mix, s3cal, s71mix, s71cal,
     .                    s18mix, s18cal,
     .                    dtimix, dtical, fusmix, fuscal, trmix, trcal
      write (nunit, 2070)
      do 825 j=1,nj
        j1prt = (j-1) / jprt * jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 825
        write (nunit, 2072) j, r(j), roa(j), tem(j), tep(j),
     .                      tim(j), tip(j), ddfusm(j), ddfusp(j)
  825 continue
      write (nunit, 2075)  sxrm(3,1),sxrp(3,1),sxrm(2,3),sxrp(2,3),
     .                     sxrm(2,2),sxrp(2,2),ddntm,ddntp
  830 continue
c
c ----------------------------------------------------------------------
c print out sxr information
c ----------------------------------------------------------------------
c
      if (jsxr .ne. 0) then
        call header (nqik, time, t)
        write (nqik, 1210)
        k1 = 1
        k2 = 1
        if (jsxr .eq. 2) k2 = 2
        write (nqik, 1220) (namar(k),namar(k),k=k1,k2)
        if (jsxr .eq. 1)  write (nqik, 1230)
        if (jsxr .eq. 2)  write (nqik, 1235)
        do 840 i=1,ndiode(k1)
  840   write (nqik, 1240) (idiode(i,k),roamin(i,k),sxr(i,k),sxcal(i,k),
     .                      k = k1,k2)
        k1 = k2+1
        k2 = k1+1
        write (nqik, 1220) (namar(k),namar(k),k=k1,k2)
        write (nqik, 1235)
        do 850 i=1,ndiode(k1)
  850   write (nqik, 1240) (idiode(i,k),roamin(i,k),sxr(i,k),sxcal(i,k),
     .                      k = k1,k2)
      end if
c
****  call torayprt (iotoray, nj, r, ene, en, kj, te, ti,
**** .               bpol, curden, q, time, t, npsi, nion)
c
c ----------------------------------------------------------------------
c Write transport bank data to yokfil in namelist format (OBSOLETE)
c ----------------------------------------------------------------------
c


      if (ilastp .eq. 0 .or. iyoka .eq. 0)  return
c
      we     =  eetot  * 1.0e-3
      wi     = (etot-eetot-eatot-ebtot) * 1.0e-3
      wfi    =  ebtot  * 1.0e-3
      wtot   =  etot   * 1.0e-3
      tau    =  taue   * 1.0e3
      taupn  =  taupin * 1.0e3
      volumm =  volume * 1.0e-6
      circcm =  circum
      axjka  =  totcur(1) * 1.0e-3
      btorkg =  btor   * 1.0e-3
      pbtor  =  ptor   * 1.0e-6
      ptotmw =  qtot   * 1.0e-6
      pradmw =  pradt  * 1.0e-6
      betatp =  beta   * 100.0
      betaep =  betae  * 100.0
      betaip =  betai  * 100.0
      betabp =  betab  * 100.0
      taueth = (we+wi) * tau / wtot
c
****  write (nyok, 1268)
****  write (nyok, 1270) ishot,itime,rminor,kappa,volumm,
**** .                   xmagn1,ymagn1,rgeom,circcm,btorkg,axjka,qstar,
**** .                   voltag,voltoh,ddntot,enebar,zfrac
****  write (nyok, 1272) enepk,tepk,tipk,curpk,
**** .                   pbtor,sthru,pfil,pradmw,
**** .                   poh,pbe,pbi,prf,ptotmw
****  write (nyok, 1274) betatp,betaep,betaip,betabp,betap,
**** .                   tau,taueth,taupn,wtot,we,wi,wfi
c
      do 910 j=1,nj
      enx  (j) = 1.0e3 * psir(j)
      w3neo(j) = 0.0
  910 if (xkineo(j) .ne. 0.0)
     .w3neo(j) = wneo(3,3) * xkiinv(j) / xkineo(j)
      j11 = nj / 10
c
****  write (nyok, 1280) 'R =',(r(j),j=1,nj,j11)
****  write (nyok, 1281) 'HCAP =',(hcap(j),j=1,nj,j11)
****  write (nyok, 1282) 'PSIR =',(enx(j),j=1,nj,j11)
****  write (nyok, 1282) 'ENE =',(ene(j),j=1,nj,j11)
****  do 920 i=1,nion
**920 write (nyok, 1285)  i,(en(j,i),j=1,nj,j11)
****  write (nyok, 1282) 'ENB =',(enbeam(j),j=1,nj,j11)
****  write (nyok, 1281) 'TE =',(te(j),j=1,nj,j11)
****  write (nyok, 1281) 'TI =',(ti(j),j=1,nj,j11)
****  write (nyok, 1282) 'CURDEN =',(curden(j),j=1,nj,j11)
****  write (nyok, 1281) 'Q =',(q(j),j=1,nj,j11)
****  write (nyok, 1281) 'ZEFF =',(zeff(j),j=1,nj,j11)
****  write (nyok, 1282) 'XKE =',(xkeinv(j),j=1,nj,j11)
****  write (nyok, 1282) 'XKI =',(xkiinv(j),j=1,nj,j11)
****  write (nyok, 1282) 'W3NEO =',(w3neo(j),j=1,nj,j11)
****  write (nyok, 1282) 'WEJCM3 =',(wejcm3(j),j=1,nj,j11)
****  write (nyok, 1282) 'WIJCM3 =',(wijcm3(j),j=1,nj,j11)
****  write (nyok, 1282) 'WBJCM3 =',(wbjcm3(j),j=1,nj,j11)
****  write (nyok, 1282) 'SLENE =',(slene(j),j=1,nj,j11)
****  write (nyok, 1282) 'SLTE =',(slte(j),j=1,nj,j11)
****  write (nyok, 1282) 'SLTI =',(slti(j),j=1,nj,j11)
****  write (nyok, 1282) 'SLPRES =',(slpres(j),j=1,nj,j11)
****  write (nyok, 1282) 'SHEARP =',(shearp(j),j=1,nj,j11)
****  write (nyok, '('' $'')
c

      return
c
 1015 format (5x,'densities for electrons and ions' //
     .        4x,'j',7x,'r',7x,'zeff',7x,'ene  ',5(5x,a2,' ions'))
 1017 format (10x,'(cm)',11x,2(3x,'(1/cm**3)'))
 1016 format (5x,'densities: fast ions, neutrals'
     .        4x,'(neutrals at time n-1/2)' //
     .        4x,'j',7x,'r',7x,'enbeam',7x,'enalp',5x,a2,' neuts',
     .        4x,a2,' neuts' /
     .       10x,'(cm)  ',4(3x,'(1/cm**3)'))
 1020 format (1x,i4,f10.2,f9.3,1x,1p9e12.3)
 1021 format (1x,i4,f10.2,1x,1p9e12.3)
 1022 format (/8x,'average',f9.3,1x,1p9e12.3)
 1025 format (/8x,'average',1x,1p9e12.3)
 1023 format (/3x,'line average',10x,1pe12.3)
 1024 format (/ '  entotn =',1pe11.3, '  enitn  =',1pe11.3,
     .          '  snaddt  =',1p2e13.3 ///)
 1030 format (5x,'density and charge for impurity element  ',a2 //
     .        4x,'j',7x,'r',8x,'en',10x,' z ',4x,'sqrt(zsq)'     /
     .       10x,'(cm)',4x,'(1/cm**3)')
 1040 format (1x,i4,f10.2,1pe12.3,0p2f10.2)
 1045 format (10x,'temperatures, magnetic field, current density, ',
     .       'electric field, safety factor, and helical flux'  //
     .       4x,'j',7x,'r',9x,'te',8x,'ti',7x,'tn1',7x,'tn2',
     .       6x,'bp',8x,'curden',6x,'etor',8x,'q',8x,'psis'     /
     .       10x,'(cm)',1x,4(5x,'(keV)'),5x,'(G)',5x,
     .       '(A/cm**2)',4x,'(v/cm)',15x,'(G-cm2)')
 1051 format (/// '  totcur =', 1pe11.3, ' a'                   /
     .            '  voltag =', 1pe11.3, ' v   (at n-1/2)'      //
     .            '  betae  =', 1pe11.3, 8x, 'betae0 =',1pe11.3 /
     .            '  betai  =', 1pe11.3, 8x, 'betai0 =',1pe11.3 /
     .            '  betab  =', 1pe11.3, 8x, 'betab0 =',1pe11.3 /
     .            '  betaa  =', 1pe11.3, 8x, 'betaa0 =',1pe11.3 /
     .            '  beta   =', 1pe11.3, 8x, 'beta0  =',1pe11.3 //
     .            '  betap  =', 1pe11.3                         //
     .            '  li     =', 1pe11.3)
 1053 format (//5x,'current densities (A/cm**2) and total currents (a)'/
     .          5x,'(curbi, curbe,curbet, currf and curohm at time '   /
     .          5x,'point n-1/2,curden at time point n)'              //
     .          4x,'j',7x,'r',6x,'curden',6x,'curohm',5x,'curboot',7x,
     .         'curbi',7x,'curbe',7x,'curbet',7x,'currf',5x,
     .         'curdrive',5x,'curpar' /
     .          10x,'(cm)')
 1054 format (/ 10x, 'total', 1p8e12.3)
 1055 format (  10x,'Neutrons produced by D-D fusion',
     .           5x,'[= 1/2 # of D-D fusions]'         //
     .           4x,'j',7x,'r',8x,'enD',9x,'ti',5x,'D-D neutrons',
     .            '  Knock-ons',3x,'D-beam neutrons'   /
     .          10x,'(cm)',4x,'(1/cm**3)',5x,'(keV)',4x,'(1/cm**3-s)',
     .           1x,'(1/cm**3-s)',1x,'(1/cm**3-s)'     /)
 1056 format (1x,i4,f10.2,1p9e12.3)
 1057 format (/6x,'total # of DD thermal neutrons:',1pe12.3,' /s',
     .       /6x,'total # of knock-on   neutrons:',1pe12.3,' /s',
     .       /6x,'total # of beam-deut  neutrons:',1pe12.3,' /s',
     .       /6x,'total neutron production  rate:',1pe12.3,' /s',
     .       /6x,'total # of beam-th_deut neutrons new :',1pe12.3,' /s',
     .      /6x, 'total # of beam-th_deut protons  new :',1pe12.3,' /s')
 1060 format (1x,i4,f10.2,4f10.4,1p6e11.3)
 1062 format (/1x,'volume-average',2f10.4)
 1064 format (/ '  note:  tn and etor are evaluated at time point ',
     .            'n-1/2.')
*1065 format (1h1,'time =',3pf8.2,' ms,   time point =',i4     //)
*1066 format (1h1,'time =',3pf8.2,' ms,   time point =',0pf6.1 //)
 1071 format (37x,'fluxes'//4x,'j',7x,'r',1x,  (7x,'ion #',i1),7x,
     .       'elec.e.',6x,'ion e.'9x,'angmtm'/10x,'(cm)',2x,  (2x,
     .       '(1/cm**2-s)'),2x,2('(keV/cm**2-s)'),1x,'(gm/sec**2)' /
     .         8x,4(10x,a))
 1072 format (37x,'fluxes'//4x,'j',7x,'r',1x, 2(7x,'ion #',i1),7x,
     .       'elec.e.',6x,'ion e.'9x,'angmtm'/10x,'(cm)',2x, 2(2x,
     .       '(1/cm**2-s)'),2x,2('(keV/cm**2-s)'),1x,'(gm/sec**2)' /
     .       8x,5(10x,a))
 1073 format (37x,'fluxes'//4x,'j',7x,'r',1x, 3(7x,'ion #',i1),7x,
     .       'elec.e.',6x,'ion e.'9x,'angmtm'/10x,'(cm)',2x, 3(2x,
     .       '(1/cm**2-s)'),2x,2('(keV/cm**2-s)'),1x,'(gm/sec**2)' /
     .       8x,6(10x,a))
 1074 format (37x,'fluxes'//4x,'j',7x,'r',1x, 4(7x,'ion #',i1),7x,
     .       'elec.e.',6x,'ion e.'9x,'angmtm'/10x,'(cm)',2x, 4(2x,
     .       '(1/cm**2-s)'),2x,2('(keV/cm**2-s)'),1x,'(gm/sec**2)' /
     .       8x,7(10x,a))
 1075 format (37x,'fluxes'//4x,'j',7x,'r',1x, 5(7x,'ion #',i1),7x,
     .       'elec.e.',6x,'ion e.'9x,'angmtm'/10x,'(cm)',2x, 5(2x,
     .       '(1/cm**2-s)'),2x,2('(keV/cm**2-s)'),1x,'(gm/sec**2)' /
     .       8x,8(10x,a))
 1080 format (i4, a2, f9.2, 1x, 1p8e13.3)
 1090 format (10x,'Neutral profiles for species  ',a2,' = ',i2,
     .        5x,a2,' = ',i2)
 1091 format (10x,'Scrape-off layer width =',f8.2,' cm')
 1092 format (/4x,'j',7x,'r',8x,'tn1',6x,'ennw1',7x,'ennv1',
     .  7x,'volsn1',7x,'vz',8x,'tn2',4x,'ennw2',7x,'ennv2',7x,'volsn2',
     .  6x,'vz' /
     .  10x,'(cm)',6x,'(keV)',2(3x,'(1/cm**3)'), '  (1/cm**3-s)',
     .  4x,'cm/sec',4x,'(keV)',2(3x,'(1/cm**3)'), '  (1/cm**3-s)',
     .   '  cm/sec')
 1094 format (i5, f10.2, f10.4, 1p4e12.3, 0pf10.4, 1p4e11.3)
 1096 format (// 10x, 'Specific flux of neutrals ',
     .'(1/cm**2-s-keV-sterad)'/10x,'rtandn = ',f8.2,' cm,  rhdn=',
     . f8.2,' cm      (f1w,f2w,fvn = ',3f6.3,')' //
     . 11x, 'energy', 7x, 'spflux1', 8x, 'spflux2' / 12x, '(keV)')
 1098 format (7x, f10.4, 3x, 1pe12.3, 3x, 1pe12.3)
 1155 format (
     . //5x,'entot  = ',1pe11.3 / 5x,'dentot = ',1pe11.3,' 1/s',
     .  /5x,'stot   = ',1pe11.3,' 1/s',13x,'taup   = ',1pe11.3,' s',
     . //5x,'eetot  = ',1pe11.3,' j'/5x,'deetot = ',1pe11.3,' w'
     .  /5x,'qetot  = ',1pe11.3,' w',
     .  15x,'tauee  = ',1pe11.3,' s'
     . //5x,'etot   = ',1pe11.3,' j' /5x,'detot  = ',1pe11.3,' w',
     .  /5x,'qtot   = ',1pe11.3,' w',15x,'taue   = ',1pe11.3,' s',
     . //5x,'volume = ',1pe11.3,' cm**3',11x,'entaue = ',1pe11.3,
     .          ' s/cm**3'//5x,'surface area (cm**2) = ',1pe11.3,
     . //5x,'ec     = ',1pe11.3,' j/cm**3'
     .  /5x,'dec    = ',1pe11.3,' W/cm**3',
     .  /5x,'qcen   = ',1pe11.3,
     . ' W/cm**3',9x,'tauec  = ',1pe11.3,' s' //
     .  5x,'angmtot (kg*m**2/sec) = ',1pe11.3 /
     .  5x,'dangmtot/dt (nt-m) = ',1pe11.3 /
     .  5x,'storquet (nt-m) = ',1pe11.3 /
     .  5x,'tauang (sec) = ',1pe11.3)
 1210 format (5x, 'sxr signals (arbitrary units)')
 1220 format (// 2(18x,a8,4x,a8))
 1230 format ('  diode  roamin',5x,'sxr',8x,'sxcal' /)
 1235 format (2( '  diode  roamin', 5x, 'sxr', 8x, 'sxcal', 2x))
 1240 format (2(i6, f8.3, 1p2e12.3))
*1268 format ('1YOKA DATA')
*1270 format (' $TRBANK'                             /
**** .        ' ISHOT = ',i8,',  ITIME=',i5, ','       /
**** .        ' RMINOR = ',f6.1,',  ELONG=',f7.2, ','  /
**** .        ' VOLUME = ',f7.2, ','                   /
**** .        ' RAX = ',f7.1,',  ZAX=',f7.1, ','       /
**** .        ' RGEOM = ',f7.1, ','                    /
**** .        ' CIRCUM = ',f7.1, ','                   /
**** .        ' BTOR = ',f7.1,',  TOTCUR=',f7.0, ','   /
**** .        ' QSTAR = ',f7.2, ','                    /
**** .        ' VOLTAG = ',f7.2,',  VOLTOH=',f7.2, ',' /
**** .        ' FUSRAT = ',1pe11.3, ','                /
**** .        ' ENEBAR = ',1pe11.3,',  ZFRAC=',0pf7.2, ',')
*1272 format (' ENEPK = ',f7.2,',  TEPK=',f7.2,',  TIPK=',f7.2,
**** .      ',  CURPK = ',f7.2, ','                    /
**** .        ' PBTOR = ',f7.2,',  STHRU=',f7.2,
**** .       ',  PBCX = ',f7.2,',  PRAD=',f7.2, ','    /
**** .         ' POHM = ',f7.2,',  PBE=',f7.2,',  PBI=',f7.2,
**** .        ',  PRF = ',f7.2,',  PTOT=',f7.2, ',')
*1274 format (' BETAT = ',f7.2, ','                                  /
**** .       ' BETAE = ',f7.2,',  BETAI=',f7.2,',  BETAB=',f7.2, ',' /
**** .       ' BETAP = ',f7.2, ','                                   /
**** .       ' TAUE = ',f7.1,',  TAUETH='f7.1,',  TAUP=',f7.1, ','   /
**** .       ' WTOT = ',f7.1, ','                                    /
**** .       ' WE = ',f7.1,',  WI=',f7.1,',  WB=',f7.1, ',')
*1280 format (' ',a8,1x,6(f6.2, ',')           /
**** .       (10x,6(f6.2, ',')))
*1281 format (' ',a8,1x,6(f6.3, ',')           /
**** .       (10x,6(f6.3, ',')))
*1282 format (' ',a8,1x,6(1pe10.3, ',')        /
**** .       (10x,6(1pe10.3, ',')))
*1285 format (' EN',i1,' =',4x,6(1pe10.3, ',') /
**** .       (10x,6(1pe10.3, ',')))
 1400 format (20x,'densities, temperatures, currents, and related',
     .       ' quantities')
 1430 format (1x,i4,f6.1,f5.1,f10.3,1x,1p9e12.3)
 1440 format (/4x,'j',4x,'r',3x,'r/a',6x,'te',7x,'ti',8x,'we',10x,'wi',
     .       8x,'wbeam',7x,'walp' /
     .       7x,'(cm)',5x,2(4x,'(keV)'),4(3x,'(J/cm**3)'))
 1450 format (1x,i4,f6.1,f5.1,2f9.3,1p4e12.3)
 1460 format (/4x,'j',4x,'r',3x,'r/a',6x,'q',6x,'etor',
     .       7x,'et+esaw',5x,'curden',6x,'curohm',5x,'curboot',
     .       5x,'curdrive',7x,'bpol',8x,'psir' /
     .       7x,'(cm)',21x,2(6x,'(V/cm)'),2x,3(3x,'(A/cm**2)'),
     .       6x,'(G)',6x,'(G-cm**2)')
 1470 format (1x,i4,f6.1,f5.1,1f9.3,1x,1p8e12.3)
 1471 format (/39x,'total currents:    ',1p4e12.3)
 1500 format (//10x,'Line-averaged electron densities measured by',
     .       ' CO2 interferometer array'
     .       // '  rtan(cm)', '  path(cm)', '  density(cm-3)')
 1510 format (1x,f8.1,f10.1,1pe12.3)
 1520 format (// 10x,
     .       'Visible continuum emission as measured by Z-effective',
     .       ' profile diagnostic'
     .       / 12x, 'zeffwl = ',f7.1,' A,  zeffdl=',f5.1,' A,  zefrad=',
     .       1pe12.3,' cm'
     .       // '  chn.',4x,'rtan(cm)',6x,'photons/s',13x,
     .       'chn.',4x,'rtan(cm)',5x,'photons/s')
 1530 format (2x,i3,4x,f8.1,4x,1pe12.3,13x,i3,4x,0pf8.1,4x,1pe12.3)
 2050 format (// '  ttweak =',1pe11.3,' s'/ '  w33    =',1pe11.3,
     .       8x,'zeffc  =',1pe11.3,8x,'zeffb  =',1pe11.3)
 2052 format (/ '  fusnin =',1pe11.3,' s-1',4x,'fusn   =',1pe11.3,' s-1'
     .        / '  ticin  =',1pe11.3,' keV',4x,'ti(1)  =',1pe11.3,' keV'
     .        / '  voltin =',1pe11.3,' v',6x,'voltag =',1pe11.3,' v'
     .        / '  qcin   =',1pe11.3,8x,'q(1)   =',1pe11.3)
 2060 format (// '  timmix =',1pe11.3,' s',6x,'rs/a   =',1pe11.3,
     .        8x,'rmix/a =',1pe11.3
     .        /'  epste  =',1pe11.3,' keV',4x,'epsti  =',1pe11.3,' keV')
 2062 format (/ '  dtemix =',1pe11.3,' keV',4x,'dtecal =',1pe11.3,' keV'
     .        / '  s3mix  =',1pe11.3,8x,'s3cal  =',1pe11.3
     .        / '  s71mix =',1pe11.3,8x,'s71cal =',1pe11.3
     .        / '  s18mix =',1pe11.3,8x,'s18cal =',1pe11.3
     .        / '  dtimix =',1pe11.3,' keV',4x,'dtical =',1pe11.3,' keV'
     .        / '  fusmix =',1pe11.3,8x,'fuscal =',1pe11.3
     .        / '  trmix  =',1pe11.3,8x,'trcal  =',1pe11.3)
 2070 format (/ 16x, 'temperatures before and after mixing' //
     .           4x,'j',4x,'r',3x,'r/a',5x,'tem',6x,'tep',
     .           6x,'tim',6x,'tip',
     .           5x,'ddfusm',6x,'ddfusp' /
     .           7x,'(cm)',5x,4(4x,'(keV)'),2(' (1/s-cm**3)'))
 2072 format (i5, f6.1, f5.1, 4f9.3, 1p2e12.3)
 2075 format (/ '  sxr(3,1)   =', 1p2e12.3 /
     .          '  sxr(2,3)   =', 1p2e12.3 /
     .          '  sxr(2,2)   =', 1p2e12.3 /
     .          '  fusn       =', 1p2e12.3)
 3015 format (// 5x,'bootstrap current from density gradient' //
     .           4x,'j',7x,'r',8x,'ft/fc',8x,'elec ',4x,'fast ions',
     .           4x,a2,' ions',4(5x,a2,' ions'))
 3016 format (// 5x,'bootstrap current from temperature gradient' //
     .           4x,'j',7x,'r',8x,'ft/fc',8x,'elec ',5(5x,a2,' ions'))
 3017 format (10x,'(cm)',13x,2(3x,'(a/cm**2)'))
 3054 format (/ 22x, 'total', 1p9e12.3)
 3055 format (/ 22x, 'total', 1p8e12.3)
c
      end

      subroutine pcoef (timet, t)
c
      USE param
      USE io 
      USE ions
      USE numbrs
      USE mesh
      USE sourc    ! eta is defined here
      USE tfact
      USE geom
      USE paleocl,only : include_paleo,paleocl_tr,paleocl_output
      USE tordlrot
      USE tcoef
      USE soln,      ONLY : ene
      implicit  integer (i-n), real*8 (a-h, o-z)
c

      include 'rebut.i'

c
      character*2 jp
c
c print out diffusion coefficients
c
      call header (nout, timet, t)
      write (nout, 1082)
      jp = '.5'
      do 360 j=1,nj-1
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj-1)  go to 360
      write (nout, 1080) j,jp,ra(j),dneo(j),dtyp(j),
     .                   dsaw(j),disl(j),drebut(j)*wrebut,d(1,1,j)
  360 continue
c
c print out electron thermal conductivities
c
      call header (nout, timet, t)
      write (nout, 1084)
      do 370 j=1,nj-1
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj-1)  go to 370
c
c --- xkerebut includes neoclassical so we don't want to count
c --- neoclassical twice in sum
c
      dumy = xkerebut(j) * wrebut
      sum  = xkeneo(j) + xketyp(j) + xkesaw(j) + xkeisl(j) + xkebal(j)
     .     + dumy
      sum  = sum - wrebut * xkeneo(j)
      write (nout, 1080) j,jp,ra(j),xkeneo(j),xketyp(j),xkesaw(j),
     .                   xkeisl(j),xkebal(j),dumy,sum
  370 continue

c
c --- print paleoclassical results
c 
      if(include_paleo .gt. 0)then
          call header (nout, timet, t)
          call paleocl_output(nout,jprt)
          write(nout,1105)
          do  j=1,nj-1
               !enea =0.5*(ene(j)+ene(j-1))
               enea =0.5*(ene(j+1)+ene(j)) ! corrected 6/7/2012 HSJ
               j1prt = ((j-1)/jprt)*jprt
               if (j1prt .ne. j-1 .and. j .ne. nj-1) go to 371
                  write (nout, 1104) j,jp,ra(j),xkeneo(j),
     .                       chie_paleo(j)*enea   ! ,paleo_flux(j)
 371              continue
          enddo
      endif


c
c print out ion thermal conductivities
c
      call header (nout, timet, t)
      write (nout, 1086)
c
      do 380 j=1,nj-1
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj-1)  go to 380
      dumy = xkirebut(j)*wrebut
      sum  = xkineo(j) + xkityp(j) + xkisaw(j) + xkiisl(j) + xkibal(j)
     .     + dumy
      sum  = sum - wrebut*xkineo(j)           ! xkirebut includes xkineo
      write (nout, 1080) j,jp,ra(j),xkineo(j),xkityp(j),xkisaw(j),
     .                   xkiisl(j),xkibal(j),dumy,sum
  380 continue
c
c print out trapped particle fraction, resistivity, resistive skin time
c and collision frequencies
c     tres = 4 * pi * L**2/eta*c**2
c
      call header (nout,timet,t)
      write (nout, 1090)
      write (nout, 1100) (nameu(k),k=1,nion)
      write (nout, 1102)
      do 400 j=1,nj-1
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj-1)  go to 400
      tres  = 0.0
      if (eta(j) .ne. 0.0)  tres = 1.396e-20*r(nj)**2/eta(j)
      write (nout, 1080) j,jp,ra(j),ftrap(j),eta(j)*8.98755e11,tres,
     .                   xnuse(j),(xnus(k,j), k = 1,nion)

  400 continue
c
c     print out effect of trapped electrons on the ohmic power, ftr
c     since eta is in CGS (sec**-1) and curohm is in amps/cm**2,
c     the  following multiplicative constant, 9.0e7,
c     give pohm1 in watts/cm**3
c
      do 500 j=1,nj
      pohm1(j) = eta(j)*curohm(j)**2
  500 pohm2(j) = (1.0 - ftrap(j))*pohm1(j)
      call trapv (r,pohm1,hcap,nj,pohm1t)
      call trapv (r,pohm2,hcap,nj,pohm2t)
      ftr = pohm1t / pohm2t
      write  (nout, 1081)  ftr
 1081 format (/// ' multiplicative effect of trapped electrons',
     .            ' on ohmic power =' , 1pe10.2)
c
 1080 format (i4,a2,f9.2,1x,1p8e13.3,1pe12.3)
 1082 format (30x,'ion diffusion coefficients (cm**2/s)' //
     .       4x,'j',7x,'r',7x,'dneo(1,1)',
     .       8x,'dtyp',9x,'dsaw',9x,'disl',9x,'drebut',9x,'sum'
     .       / 10x,'(cm)')
 1084 format (30x,'electron thermal conductivities (1/cm-s)' //
     .       4x,'j',7x,'r',9x,'xkeneo',7x,'xketyp',
     .       7x,'xkesaw',7x,'xkeisl',7x,'xkebal',8x,'xkerebut',8x,'sum'
     .       / 10x,'(cm)')
 1086 format (30x,'ion thermal conductivities (1/cm-s)' //
     .       4x,'j',7x,'r',9x,'xkineo',7x,'xkityp',7x,'xkisaw',
     .       7x,'xkiisl',7x,'xkibal',7x,'xkirebut',7x,'sum' /10x,'(cm)')
 1090 format (20x,'trapped electron fraction, resistivity, ',
     .       'resistive skin time and collision frequencies'//)
 1100 format (4x,'j',8x,'r',8x,'ftrap',9x,'eta',9x,'tres',
     .        9x,'xnuse',7x,5(a2,'-xnus',6x))
c 1102 format (11x,'(cm)',18x,'(ohm-cm)',10x,'(s)',10x,'(s)')
 1102 format (11x,'(cm)',18x,'(ohm-cm)',10x,'(s)',10x,' ') ! xnuse is unitless HSJ 10/1/10
 1104 format (i4,a2,f9.4,1x,3(1pe12.4))
 1105 format (30x,'electron paleoclassical conductivities (1/cm-s)' //
     .       4x,'j',7x,'r',9x,'xkeneo',7x,'xke_paleo')      ! ,7x ,
                   ! .       ' paleo flux  #/(cm2 sec)')
      return
c
      end

      subroutine pflux (timet, t)
c
      USE param
      USE io 
      USE numbrs
      USE mesh
      USE geom
      USE tordlrot
      USE tcoef
      USE flx
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c this subroutine outputs flux components
c
c      include 'param.i'
c      include 'geom.i'
c      include 'io.i'
c      include 'numbrs.i'
c      include 'flx.i'
c      include 'mesh.i'
c      include 'tcoef.i'
c      include 'tordlrot.i'
c
      character*2 jp
      dimension   fluxp(kk), etorp(kk)
c
c write out energy fluxes (keV/cm**2-s)
c

      call visc_flux ! HSJ 5/5/14 load qomegapi and flxangce
                     ! for initial time.
      call header (nout,timet,t)
      write (nout,1010)
      write (nout,1016)
      write (nout,1015)
      write (nout,1020)
      jp = '.5'
      do 110 j=1,nj-1
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj-1)  go to 110
      conde = 0.0
      condi = 0.0
      do 105 k=1,nk
      conde = conde - d(nion+1,k,j)*dudr(k,j)
  105 condi = condi - d(nion+2,k,j)*dudr(k,j)
      conve = flux(nion+1,j) - conde
      convi = flux(nion+2,j) - condi-iangrot*angrcple*(omegapi(j)
     .      + flxangce(j))
      write (nout, 1070) j,jp,ra(j),conde,conve,flux(nion+1,j),
     .                   condi,convi,omegapi(j),flxangce(j),
     .                                          flux(nion+2,j)

 
  110 continue
c
c  write out header information
c
      do 150 i=1,nion+3
      call header (nout,timet,t)
      if (i .gt. nion)  go to 112
      write (nout, 1040) i
      call flxtitle(nout,i,nk)
      go to 118
 112  if (i .ne. nion+1)  go to 114
      write (nout, 1080)
      call flxtitle(nout,i,nk)
      go to 118
 114  if (i .ne. nion+2)  go to 116
      write (nout, 1100)
      call flxtitle(nout,i,nk)
      go to 118
 116  write (nout, 1110)
      call flxtitle(nout,i,nk)
c
c  write out flux components
c
  118 do 150 j=1,nj-1
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj-1)  go to 150
      sum   = 0.0
c
      do 130 k=1,nk
        if (i .eq. nion+3)  go to 120
        fluxp(k) = -d(i,k,j)*dudr(k,j)
        sum      = sum + fluxp(k)
        go to 130
  120   etorp(k) = 1.0e-8 * ra(j) * d(i,k,j) * dudr(k,j)
     .                    * 0.5 * (hcap(j)+hcap(j+1))
        sum = sum + etorp(k)
  130 continue
c
      if (i .eq. nion+3)  go to 140
      write (nout,1070) j,jp,ra(j), (fluxp(k), k=1,nk),sum
      go to 150
  140 write (nout,1070) j,jp,ra(j), (etorp(k), k=1,nk),sum
  150 continue
c
*1000 format (1h1, 'time = ', 3pf8.2,' ms,   time point =',0pf6.1//)
 1010 format (30x, 'energy fluxes (keV/cm**2-s)'//)
 1015 format (28x, 'electron energy', 27x, 'ion energy')
 1016 format (' NOTE: in this table conde and condi are calculated'  /
     .        ' from models (i.e., neoc, rlw, etc.) conve and convi' /
     .        ' are defined as energy flux - cond . The energy'      /
     .        ' flux in turn is given in the previous table '        /
     .        ' and was determined either from models (simulation)'  /
     .        ' or from div flux = source (analysis)')
 1020 format ('  j',9x,'r',8x,'cond',9x,'conv',9x,'sum',10x,'cond',9x,
     .        'conv',9x,'omegapi',6x,'cvctvrot',6x,'sum' / 10x,'(cm)')
 1040 format (20x, 'components of particle flux ',
     .             '(1/cm**2-s) for ion#',i3//)
 1070 format (i4,a2,f9.2,1p9e13.3)
 1080 format (20x,'components of electron energy flux (keV/(cm**2*s)',
     .           ' due to conduction'//)
 1100 format (20x,'components of ion energy flux (keV/(cm**2*s)',
     .           ' due to conduction'//)
 1110 format (20x,'terms in generalized ohm''s law (V/cm)'//)
      return
c
      end







      subroutine psourc (timet, t)
c

c
c --- this subroutine prints out details of the particle and energy
c --- and angular momentum sources
c
c --- Revisions: 08/30/94 G.M. Staebler - added ANGROT and SPERP to
c                                         output plotfile trpltfil
c
c ----------------------------------------------------------------------
c
      USE param
      USE fusion
      USE ions
      USE io 
      USE glf23
      USE neut
      USE nub
      USE nub2
      USE solcon
      USE soln
      USE mhdpar
      USE rf
      USE extra
      USE yoka
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE tfact
      USE geom
      USE flags
      USE tordlrot
      USE constnts
      USE tcoef
      USE bd_condtn
      USE mixcom
      USE rhog
      USE ifs
      USE staebler
      USE flx
      USE weiland
      USE pelcom
      USE P_Nfreya_12_interface,              ONLY : use_P_Nfreya
c      USE flx,                                ONLY : anal_eng_flux_e,
c     .                                               anal_eng_flux_i

      USE flx ! ifort doesnt like above one 

      implicit  integer (i-n), real*8 (a-h, o-z)




      include 'exptlprof.i'

      include 'rebut.i'
      include 'storage.i'
      include 'shay.i'


      complex*16   zet, zz, zp
      dimension psum(kj), pdum(kj), qdum(kj)
      dimension dudt(kj), dflux(kj),dpidt_tot(kj)
      dimension ss(kj)

      dimension srcint(kj),swork(kj)
      dimension energe(kj), energi(kj)

      dimension pdpedt(kj), psawe (kj), pbeame_rot(kj),
     .          pdelte(kj), pexch (kj), pohmic(kj), pione(kj),
     .          prad  (kj), pbeame(kj), prfe  (kj), pfuse(kj),
     .          ptfuse(kj), pbfuse(kj), pe2d  (kj), pmag (kj)

      dimension pdpidt(kj), psawi (kj), pbeami_rot(kj),
     .          pdelti(kj), pneut (kj), qneut (kj),
     .          pioni (kj), pcx   (kj), pbeami(kj), prfi  (kj),
     .          pfusi (kj), pi2d  (kj), ptfusi(kj), pbfusi(kj)
      dimension tauetr(kj),tauitr(kj),tauir(kj),taupe(kj),tauelc(kj)
      dimension pbeamf(ke,kb),pbeams(ke,kb)
      dimension fpe(ke,kb),fpi(ke,kb),fpcx(ke,kb)
      dimension pomegapi(kj),pomegale(kj),amassden(kj),pangce(kj)
      equivalence (pomegapi(1), xdum(1))
      equivalence (pomegale(1), ydum(1))
c
      real*8     realjp
      character  jp*2, moderun*10
c
      data       five_halfs_te,five_halfs_ti /2.5,2.5/

c
c  calculate normalization factor for averages
c
*%%%% write (6, *) '%%%% entering subroutine PSOURC'
c
      zero              = 0.0
      one_third         = 1.0 / 3.0
      ten_to_minus_10th = 1.0e-10
      ten_to_minus_20th = 1.0e-20
      if(no_te_convection .eq. -1)five_halfs_te =1.5
      if(no_ti_convection .eq. -1)five_halfs_ti =1.5
      if(no_te_convection .eq. 1)five_halfs_te = 0.0
      if(no_ti_convection .eq. 1)five_halfs_ti = 0.0
c
      do 10 j=1,nj
   10   fact(j) = 4.0 * pi**2 * rmajor*hcap(j)    ! d(volume)/drho
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling TRAP'
      cconst = 1.0
      call trap1 (r, fact, cconst, nj, volume)     ! get total volume
      do 20 j=1,nj
   20 fact(j) = fact(j) / volume                  ! fractional volume
c
c ----------------------------------------------------------------------
c calculate corrected fluxes at time point theta.
c calculate energy fluxes (keV/cm**2-s)
c ----------------------------------------------------------------------
c

      onemt = 1.0 - theta
      do 130 j=1,nj-1
      do 120 k=1,nk                         ! nk = nprim+nimp+3+iangrot
c
c    for impurities itran(k) is set to zero in file "cray101.f"
c    itran(k) = 0 for analysis mode dependent variables
c
      du        = theta*(u(k,j+1)-u(k,j))
     .          + onemt*(usave(k,j+1)-usave(k,j))
      dudr(k,j) = du/dr(j)              ! at time point theta
 120  continue                          ! get all gradients first HSJ
      do k=1,nk                         ! now apply these new gradientsa
         ! flux(k,j) = 0.0              ! HSJ corrected 8/9/10
         if (itran(k) .eq.1) then       ! see below for analysis mode
           flux(k,j) = 0.0              ! HSJ corected 8/9/10
           do i=1,nk
             flux(k,j) = flux(k,j)
     .        - d(k,i,j)*dudr(i,j)      ! flux is on half grid
           end do
           if (k .le. nprim) then       ! HSJ 2/14/96
               ena=0.5*(theta*(u(k,j+1)+u(k,j))
     .                   +onemt*(usave(k,j+1)+usave(k,j)))
               flux(k,j)=flux(k,j) +ena*vpinch(j)
           end if
         end if
      end do
  130 continue

c
c --- nion = nprim+nimp,nion+1=>te,nion+2=>ti,nion+3=>rbp,nion+4=>rotation
c --- if itran(nion+1),itran(nion+2) =0 then conde and/or condi
c --- are not used so the following is o.k. ... HSJ
c --- conde,condi,conve,convi on half grid
c
 
      do j=1,nj-1
        k = nion + 1                  ! for te
        conde(j) = flux(k,j)          ! used only if itran(nion+1) = 1
        temm = theta*u(k,j)  +onemt*usave(k,j  ) ! te at r(j  ), theta
        tepp = theta*u(k,j+1)+onemt*usave(k,j+1) ! te at r(j+1), theta
        conve(j)  = five_halfs_te * (MAX (fluxe(j), zero) * temm +
     .                     MIN (fluxe(j), zero) * tepp)
        if (no_te_convection .eq. 1) conve(j) = 0.0
c        flux(k,j) =   conde(j) + conve(j)      ! disable 4/7/14 HSJ
        IF(itran(k) == 0)THEN
            conde(j) = flux(k,j) - conve(j)      ! enable 4/7/14 HSJ
        ELSE
            flux(k,j) =   conde(j) + conve(j)  
        ENDIF
c
        k = nion + 2            ! for ti
        condi(j)  = flux(k,j)   ! used only if itran(nion+2) = 1
        timm      = theta*u(k,j  ) + onemt*usave(k,j  )
        tipp      = theta*u(k,j+1) + onemt*usave(k,j+1)
        convi(j)  = five_halfs_ti * (MAX (fluxi(j), zero)*timm +
     .                     MIN (fluxi(j), zero)*tipp)
        if (no_ti_convection .eq. 1)  convi(j) = 0.0
c         flux(k,j) = condi(j) + convi(j)      ! disable 4/7/14 HSJ
        IF(itran(k) == 0)THEN
            condi(j) = flux(k,j)-convi(j) -omegapi(j) - flxangce(j)  
        ELSE
            flux(k,j) = condi(j) + convi(j) 
        ENDIF

c
        if (iangrot .eq. 1) then
c
c ---     as above, if itran(nk) = 0 then only the convective part
c ---     (fluxangc, calculated here) will be used later in this routine
c
          fluxangv(j) = flux(nk,j)
          angrota     = (theta*(u(nk,j)+u(nk,j+1))+onemt*
     .                  (usave(nk,j)+usave(nk,j+1)))*0.5
          flxmass     = 0.0
          do 145 k=1,nprim
            if (k .gt. nprim)  go to 145
            flxmass = flxmass+atw(k)*flux(k,j)
  145     continue
          fluxangc(j) = 0.5 * (r2capi(j)
     .                       + r2capi(j+1)) * angrota * flxmass * xmassp
          flux(nk,j)  = fluxangc(j)+fluxangv(j)
        else
          fluxangv(j) = 0.0
          fluxangc(j) = 0.0
          xkangrot(j) = 0.0
          flxangce(j) = 0.0
          omegapi (j) = 0.0
          vischeat(j) = 0.0
          qangc   (j) = 0.0
        end if
      end do


c
c extrapolate for end points
c
      do k=1,nk
*%%%%   write (6, *) '%%%% in subroutine PSOURC, calling EXTRAP'
        call extrap (ra(nj-2), ra(nj-1), r(nj),
     .               flux(k,nj-2), flux(k,nj-1), fluxb(k))
      end do
c
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling EXTRAP (5 times)'
      call extrap (ra(nj-2),ra(nj-1),r(nj),conde(nj-2),
     .             conde(nj-1),condeb)
      call extrap (ra(nj-2),ra(nj-1),r(nj),conve(nj-2),
     .             conve(nj-1),conveb)
      call extrap (ra(nj-2),ra(nj-1),r(nj),condi(nj-2),
     .             condi(nj-1),condib)
      call extrap (ra(nj-2),ra(nj-1),r(nj),convi(nj-2),
     .             convi(nj-1),convib)
      call extrap (ra(nj-2),ra(nj-1),r(nj),fluxangc(nj-2),
     .             fluxangc(nj-1),fluxangc(nj))
c
c --- if te, ti are in simulation mode then calculate the divergence of
c --- conduction term, qconde, qcondi
c
      if (itran(nion+1) .eq. 1)     ! te
     .  call divflx (conde,condeb,hcap,r,ra,drr,1,nj,1,qconde)
      if (itran(nion+2) .eq. 1)     ! ti
     .  call divflx (condi,condib,hcap,r,ra,drr,1,nj,1,qcondi)


c
c --- calculate divergence of convection terms for te,ti,and rotation,
c --- for analysis or simulation mode
c
      call divflx (conve,conveb,hcap,r,ra,drr,1,nj,1,qconve)
      call divflx (convi,convib,hcap,r,ra,drr,1,nj,1,qconvi)
      call divflx (fluxangc,fluxangc(nj),hcap,r,ra,drr,1,nj,1,qangc)
c
c ----------------------------------------------------------------------
c print out components of particle sources and
c particle balance table for each primary ion species
c ----------------------------------------------------------------------
c
      do 390 k=1,nprim

      call header (nout, timet, t)
      i = 0
      if (k .le. 2)  i = ineut(k)
      write (nout, 1110) k, namep(k), itran(k),i
c
c     sbfus is beam-thermal dt fusion rate (other,non d-t, reactions
c     are neglected). we can have d_f(t_th,n)he and t_f(d_th,n)he  HSJ

      do 320 j=1,njs
c      sfust = 0.0
      sfust = -stfus(j)              !HSJ 08/17/04
c
      if (ibion .gt. 0) then        !ibion = -1 if beam is dt mixture
                                    !and ifus < 0.
                                    !otherwise ibion is such that 
                                    !fast ion species is same as 
                                    !thermal species namep(ibion)
c
c         beam, if present, must be dt mixture,(ifus >0):
c
c          if (k .eq. idt)  sfust = sfust - 2.0*(stfus(j) + sbfus(j))
c         sfust is set to -stfus above so only subtract once:  HSJ 08/17/04
          if (k .eq. idt)  sfust = sfust - stfus(j) -2.0* sbfus(j)
c
c         beam, if present, could be d or t . if beam is t then there
c         is thermal depletion of t. if beam is d then there is
c         thermal depletion of d:
c
c          if (k .eq. id .and. ibion .eq. it)
c     .                 sfust = sfust -     (stfus(j) + sbfus(j))
          if (k .eq. id .and. ibion .eq. it)
     .                 sfust = sfust -  sbfus(j) 
c          if (k .eq. it .and. ibion .eq. id)
c     .                  sfust = sfust -     (stfus(j) + sbfus(j))
          if (k .eq. it .and. ibion .eq. id)
     .                  sfust = sfust - sbfus(j)
c
      else   ! beam is (single fluid) dt mixture,thermal species
c              is d and t (two separate fluids)
c              first get t_f(d_th,n)he4:
c
          jk=3*nbeams+1
c          if (k .eq. id)  sfust = sfust - (stfus(j) +
c     .                            beam_thermaldth_tf(j,jk))
         if (k .eq. id)  sfust = sfust - beam_thermaldth_tf(j,jk)
c
c         now get d_f(t_th,n)he4
c
c          if (k .eq. it)  sfust = sfust - (stfus(j) +
c     .                               beam_thermaltth_df(j,jk))
          if (k .eq. it)  sfust = sfust - beam_thermaltth_df(j,jk)
      end if
c
c     sfus is total dt fusion rate,gives source of he
c     from beam-thermal and thermal-thermal dt reactions:      HSJ
c
c      if (k .eq. ihe)  sfust =              sfust    + sfus (j)
      if (k .eq. ihe)  sfust =  sfus (j)
      ss(j)    = sfust
      sother   = s2d(j,k)
      if (itimav .eq. 1)  sother = sother + ssaw(j,k)
      swork(j) = sother
      j1prt    = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj .and. j .ne. njs)  go to 320
      sbm    = 0.0
      if (k .eq. ibion) sbm = sbeam(j)
      if (ibion .lt. 0) then ! beam slowing down source is dt fluid HSJ
         if (k .eq. id)  sbm = sbeam(j)*fdbeam
         if (k .eq. it)  sbm = sbeam(j)*(1.-fdbeam)
      end if
      sion1  = 0.0
      srcom1 = 0.0
      scx1   = 0.0
      sbcx1  = 0.0
      if (k .gt. 2)  go to 310
      sion1  = sion(j,k)
      srcom1 = srecom(j,k)
      scx1   = scx(j,k)
      sbcx1  = sbcx(j,k)
  310 sum    = s(k,j)
      if (itimav .eq. 1)  sum = sum + ssaw(j,k)
      write  (nout, 1155)  j, r(j), roa(j), sion1, srcom1, scx1, sbm,
     .                     sbcx1, sfust, sother, sum
 1155 format (1x, i4, f6.1, f5.1, 1p9e12.3)
  320 continue



c
      if (k .gt. 2)  go to 330
c      write (6, *) '%%%% in subroutine PSOURC, calling TRAPV (4 times)'
      call trapv (r,sion(1,k),fact,nj,sion1)
      call trapv (r,srecom(1,k),fact,nj,srcom1)
      call trapv (r,scx(1,k),fact,nj,scx1)
      call trapv (r,sbcx(1,k),fact,nj,sbcx1)
c      write (6, *) '%%%% in subroutine PSOURC, done TRAPV (4 times)'
  330 sbm = 0.0
      if (k .eq. ibion)
     .call trapv (r,sbeam,fact,nj,sbm)
      if (ibion .lt. 0 .and. k .eq. id) then
                 call trapv (r,sbeam,fact,nj,sbm)
                 sbm=fdbeam*sbm
      else if (ibion .lt. 0 .and. k .eq. it) then
                 call trapv (r,sbeam,fact,nj,sbm)
                 sbm=(1.-fdbeam)*sbm
      end if
      call trapv (r,swork,fact,nj,sother)
      call trapv (r,ss,fact,nj,sfust)
      do 340 j=1,nj
  340 swork(j) = s(k,j)
      call trapv (r,swork,fact,nj,sum)
      write (nout, 1152) sion1,srcom1,scx1,sbm,sbcx1,sfust,sother,sum
c
c now compute particle balance table for ion k.
c also compute the particle confinement time.
c
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling DIVFLX'
      call divflx (flux, fluxb, hcap, r, ra, drr, k, nj, kk, dflux)
c
      do j=1,nj
        dudt(j) = 0.0
        ss(j) = s(k,j)
        if (   dtt .ne. 0.0 .and. itran(k) .ne. 0) then
                dudt(j) = (u(k,j)-usave(k,j))/dtt
        else if(itran(k) .ne. 0)then
                dudt(j)=ss(j)-dflux(j)                  ! HSJ 1/13/96
        end if
        if (itimav .eq. 1  )  dudt(j) = (uav(k,j)-uav0(k,j))*dtsumi
        if (itimav .eq. 1  )  ss  (j) =  ss(j) + ssaw(j,k)
        swork(j) = ss(j) - dudt(j)           ! div flux
        if (itran(k) .eq. 0)
     .    dflux(j) = swork(j)         ! analysis mode get dflux this way
        if (k .eq. in)  dinv(j) = swork(j)   ! HSJ 2/9/96
      end do
c
      dflux(nj) = ss(nj) - dudt(nj)
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling TRAP3'
c     srcint(j)=volume integral of(div flux)= #/sec ions leaving region j
      call trap3 (r, swork, nj, hcap, volfac, srcint)
      sum = 0.0
      call zeroa (swork, nj)
c
      do j=2,nj
        sum = sum + trapf(j,r,en(1,k),hcap,volfac)
        if (srcint(j) .ne. 0.0)  swork(j) = sum / srcint(j)
      end do
c
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling EXTRAP'
      call extrap (r(2),r(3),r(1),swork(2),swork(3),swork(1))
      write (nout, 1090)
c
      do 380 j=1,nj
        dudtsv(k,j) = dudt(j) ! save for output in subroutine ITER_DBASE
        j1prt       = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 380
        err         = ss(j) - dudt(j) - dflux(j)
        write (nout, 1155) j, r(j), roa(j), dudt(j), dflux(j),
     .                     ss(j), err, swork(j)
  380 continue
c
      call trapv (r,dudt,fact,nj,dua)
      call trapv (r,dflux,fact,nj,dfa)
      call trapv (r,ss,fact,nj,ssa)
      err=ssa - dua - dfa
      write (nout, 1152)  dua, dfa, ssa, err




  390 continue   ! end loop over primary ions






c
c --- print a balance table for toroidal momentum
c --- note: if itran(nk) = 1 we use the flux calculated at the beginning
c --- of this routine. if itran(nk) = 0 we recalculate the flux from its
c --- divergence below.
c
      iangne0 : if (iangrot .ne. 0) then 
c
c --- calculate the divergence of the angular momentum flux
c ---   = (1.0 / (hcap*r)) * (d/dr)(hcap*r*flux)
c
        if (itran(nk) .eq. 1)
     .  call divflx(flux,fluxb,hcap,r,ra,drr,nk,nj,kk,dflux)
c
c --- get the total time rate of change of angular momentum
c --- terms from d(omega)/dt and d(density)/dt
c
        do 3000 j=1,nj
          amassden(j) = 0.0
          dudt(j)     = 0.0
            dudtpp    = 0.0
            do 3100 k=1,nion
              dens1   = u(k,j)
              dens2   = usave(k,j)
              if (itimav .eq. 1)  dens1 = uav(k,j)
              if (itimav .eq. 1)  dens2 = uav0(k,j)
              amassden(j) = amassden(j)+atw(k)*dens1
              if (k .gt. nprim)  go to 3100
              dudtp = 0.0
              if (   dtt .ne. 0.0)  dudtp = (dens1-dens2)/dtt
              if (itimav .eq. 1  )  dudtp = (dens1-dens2)*dtsumi
              dudtp  = dudtp*atw(k)
              dudtpp = dudtpp+dudtp
 3100       continue
            amassden(j) = amassden(j)*xmassp
            aomega1     = u(nk,j)
            aomega2     = usave(nk,j)
            if (itimav .eq. 1)  aomega1 = uav (nk,j)
            if (itimav .eq. 1)  aomega1 = uav0(nk,j)
            dudtpp = dudtpp*aomega1*r2capi(j)*xmassp
            if (   dtt .ne. 0.0)  dudtp = (aomega1-aomega2)/dtt
            if (itimav .eq. 1  )  dudtp = (aomega1-aomega2)*dtsumi
            dudt(j)     = dudtpp+amassden(j)*r2capi(j)*dudtp
            avgangmt(j) = amassden(j)*r2capi(j)*aomega1
            ss(j)       = s(nk,j)   ! total torque including beam
                                    ! uses beam torque only in sub printang
        if (itran(nk) .eq. 1)  go to 3000
        dflux(j) = ss(j)-dudt(j)
        qangv(j) = dflux(j)-qangc(j)
 3000   swork(j) = ss(j)-dudt(j)
c
        dflux(nj) = swork(nj)
        call trap3 (r, swork, nj, hcap, volfac, srcint)
        sum       = 0.0
        call zeroa (swork,nj)
c
        do j=2,nj
          sum = sum + trapf (j, r, avgangmt, hcap, volfac)
          if (srcint(j) .ne. 0.0)  swork(j) = sum / srcint(j)
        end do
c
        call extrap (r(2),r(3),r(1),swork(2),swork(3),swork(1))
        write (nout, 1091)
c
        do 3300 j=1,nj
          j1prt = ((j-1)/jprt)*jprt
          if (j1prt .ne. j-1 .and. j .ne. nj)  go to 3300
          err   = ss(j) - dudt(j) - dflux(j)
          write (nout, 1155) j,r(j),roa(j),dudt(j),dflux(j),
     .                       ss(j),err,swork(j)
 3300   continue
c
        call trapv (r,dudt,fact,nj,dua)
        call trapv (r,dflux,fact,nj,dfa)
        call trapv (r,ss,fact,nj,ssa)
        err = ssa - dua - dfa
        write (nout, 1152)  dua, dfa, ssa, err
c
c --- if itran(nk) = 0 we need to recalculate the viscous and total flux
c
      if (itran(nk) .eq. 0) then
*%%%%   write (6, *) '%%%% in subroutine PSOURC, calling FLXCAL'
        call flxcal (qangv,drr,hcap,nj,r,ra,fluxangv)
        do j=1,nj-1
          flux(nk,j) = fluxangv(j)+fluxangc(j)
          call extrap (ra(nj-2), ra(nj-1), r(nj),
     .                 flux(nk,nj-2), flux(nk,nj-1), fluxb(nk))
          denomtr    = dudr(nk,j)*rmajor**2*amassden(j)
          if (denomtr .ne. 0.0)  xkangrot(j) = -fluxangv(j)/denomtr
        end do
      end if


 

c
c --- at this point we have the total, viscous, and convective angular momentum
c --- flux for both settings of itran(nk). Now get subsidiary
c --- quantities and make corrections to the ion energy flux. (these
c --- corrections will be valid only if itran(nion+2) = 1, since flux(nion+2,j)
c --- is recalculated in the ion energy table output below if
c --- itran(nion+2) = 0)
c --- quantities related to ion energy flux
c

          if (angrcple .ne. 0.0)call visc_flux

      end if  iangne0 




c
c ----------------------------------------------------------------------
c print out components of electron particle sources
c and electron particle balance table
c ----------------------------------------------------------------------
c
c      write (6, *) '%%%% in subroutine PSOURC, calling HEADER #2'
      call header (nout, timet, t)
*%%%% write (6, *) '%%%% in subroutine PSOURC, done  HEADER #2'
      write (nout, 1115)
c
      do 430 j=1,njs
        sion1  = 0.0
        srcom1 = 0.0
        simp   = 0.0
        s2ds   = 0.0
        sbm    = sbion(j)
        do 410 k=1,nprim
        if (k .le. 2) sion1  =  sion1 + z(j,k)*sion  (j,k)
        if (k .le. 2) srcom1 = srcom1 + z(j,k)*srecom(j,k)
  410   s2ds = s2ds + z(j,k)*s2d(j,k)
        if (nimp .eq. 0)  go to 421
c
        do 420 i=1,nimp
        k    = nprim + i
        simp = simp + dzdtim(j,i)*en(j,k)
  420   s2ds = s2ds + z(j,k)*s2d(j,k)
c
  421   swork(j) = simp
        ss(j)    = s2ds
        ssawt    = 0.0
        if (itimav .eq. 1)  ssawt = ssawe(j)
        sum      = sione(j) + ssawt
        j1prt    = ((j-1)/jprt) * jprt
        if (j1prt .ne. j-1 .and. j .ne. nj .and. j .ne. njs)  go to 430
        write (nout, 1155) j,r(j),roa(j),sion1,srcom1,
     .                     simp,sbm,s2ds,ssawt,sum
  430 continue
c
      call trapv (r,sion(1,1),fact,nj,sion1)
      call trapv (r,srecom(1,1),fact,nj,srcom1)
      sion2  = 0.0
      srcom2 = 0.0
      if (nprim .ge. 2)  call trapv (r,sion(1,2),fact,nj,sion2)
      if (nprim .ge. 2)  call trapv (r,srecom(1,2),fact,nj,srcom2)
      sion1  = z(1,1) * sion1  + z(1,2) * sion2
      srcom1 = z(1,1) * srcom1 + z(1,2) * srcom2
      call trapv (r,swork,fact,nj,simp)
      call trapv (r,sbion,fact,nj,sbm)
      call trapv (r,ss,fact,nj,s2ds)
      if (itimav .eq. 1)  call trapv (r,ssawe,fact,nj,ssawt)
      call trapv (r,sione,fact,nj,sum)
      sum = sum + ssawt
      write (nout, 1152) sion1,srcom1,simp,sbm,s2ds,ssawt,sum
c
c now compute electron particle balance and confinement time
c
      call extrap (ra(nj-2),ra(nj-1),r(nj),
     .             fluxe(nj-2),fluxe(nj-1),fluxeb)
      call divflx (fluxe,fluxeb,hcap,r,ra,drr,1,nj,1,dflux)
      do 460 j=1,nj
      dudt(j) = 0.0
      ss(j)   = sione(j)
      if (   dtt .ne. 0.0)  dudt(j) = (ene(j)-enesav(j))/dtt
      if (itimav .eq. 1  )  dudt(j) = (eneav1(j)-eneav0(j))*dtsumi
      if (itimav .eq. 1  )  ss  (j) = ss(j) + ssawe(j)
  460 swork(j)  = ss(j) - dudt(j)
      dflux(nj) = ss(nj) - dudt(nj)
      call trap3(r,swork,nj,hcap,volfac,srcint)
      sum = 0.0
      call zeroa (taupe,nj)
      do 470 j=2,nj
      sum = sum + trapf(j,r,ene,hcap,volfac)
  470 if (srcint(j) .ne. 0.0) taupe(j) = sum/srcint(j)
      call extrap (r(2),r(3),r(1),taupe(2),taupe(3),taupe(1))
      write (nout,1090)
c
      do 480 j=1,nj
        j1prt = ((j-1)/jprt) * jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 480
        err   = ss(j) - dudt(j) - dflux(j)
        write (nout, 1155)  j,r(j),roa(j),dudt(j),dflux(j),
     .                      ss(j),err,taupe(j)
 480  continue
c




      call trapv (r,dudt,fact,nj,dua)
      call trapv (r,dflux,fact,nj,dfa)
      call trapv (r,ss,fact,nj,ssa)
      err = ssa - dua - dfa
      write (nout, 1152)  dua, dfa, ssa, err
c
c ----------------------------------------------------------------------
c print out components of electron energy source
c ----------------------------------------------------------------------
c
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling HEADER #3'
      call header (nout, timet, t)
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling HEADER #4'
      call header (nqik, timet, t)
      k       =  nion + 1                    ! k points to te
      moderun = 'simulation'
      if (itran(k) .eq. 0)  moderun = 'analysis'
      write (nout, 1120) moderun
      write (nout, 1121)
      write (nqik, 1120) moderun
      write (nqik, 1121)
      if (itran(k) .eq. 1)
     .  call divflx(flux, fluxb, hcap, r, ra, drr, k, nj, kk, dflux)
      do 520 j=1,nj
      dpedt(j) = 0.0
      ss(j)    = s(k,j)
****  if (n .eq. 0)  go to 510
      denfdt   = 0.0
      if (dtt .ne. 0.0)
     .  denfdt = ((enbeam(j)-enbs(j))+2.0*(enalp(j)-enasav(j)))/dtt
      if (itimav .eq. 1)
     .  denfdt = ((enbav1(j)-enbav0(j))+2.0*(enaav(j)-enaav0(j)))*dtsumi
      ss(j)    = ss(j) + 1.5*te(j)*denfdt
      if (itimav .eq. 1 .and. w2mix .gt. 0.0) ss(j) =
     .                                        ss(j) + qmag(j) + qsawe(j)
      if (itimav .eq. 1 .and. w2mix .gt. 0.0) qohm(j) =
     .                                        qohm(j) + qmag(j)
      if (dtt .ne. 0.0)
     .  dpedt(j) = 1.5*(ene(j)*u(k,j)-enesav(j)*usave(k,j))/dtt
      if (itimav .eq. 1)
     .  dpedt(j) = 1.5*(eneav1(j)*uav(k,j)-eneav0(j)*uav0(k,j))*dtsumi
      if (itran(k) .eq. 1)  go to 510
      qconde(j) = -qdelt(j)-qexch(j)+qohm(j)-qione(j)-qrad(j)+qbeame(j)
     .  + qrfe(j)+qfuse(j)+qe2d(j)
     .  -omegale(j) -dpedt(j) -qconve(j)  
c     Note qconde is equivalent to 
c         s(nion+1,j)  -dpedt(j) -qconve(j) 
c               here s includes  -qdelt(j)
c      s(nion+1,j) is equal to ssdum:

   

      if (itimav .eq. 1 .and. w2mix .gt. 0.0)  qconde(j) =
     .                                         qconde(j) + qsawe(j)
      if (                    w2mix .lt. 0.0)  qconde(j) =
     .                                         qconde(j) + qsawe(j)
      !dflux(j)  = ss(j) - dpedt(j)
      ! dflux(j) = ss(j) 
      !dflux(j)  = qconde(j) + dpedt(j) + qconve(j)
      !dflux(j)   = qconde(j) + dpedt(j)
      !dflux(j)  = ss(j) - dpedt(j)
      dflux(j)   = ss(j) -  qdelt(J)
  510 if (j .ne. nj)  go to 515
      !dflux(j)  = ss(j) - dpedt(j) 
      !dflux(j)  = ss(j) 
      !dflux(j)  = qconde(j) + dpedt(j) + qconve(j) 
      !dflux(j)   = qconde(j) + dpedt(j)
      !dflux(j)  = ss(j) - dpedt(j)
      dflux(j)   = ss(j) -  qdelt(j)
  515 j1prt     = ((j-1)/jprt) * jprt
      if (j1prt .ne. j-1 .and. j .ne. nj)  go to 520
      qd        = -qdelt(j)
      qg        = -qexch(j)
      qie       = -qione(j)
      qr        = -qrad(j)
      dpedtc(j) = dpedt(j)        ! copy to be saved in sourc.i
      dpedtj    = dpedt (j) * joupkev
      qcondj    = qconde(j) * joupkev
      qconvj    = qconve(j) * joupkev
      qd        = qd * joupkev
      qg        = qg * joupkev
      qohmj     = qohm(j) * joupkev
      qie       = qie * joupkev
      qr        = qr * joupkev
      qomegaj   = -omegale(j)*joupkev
      write (nout, 1155) j,r(j),roa(j),dpedtj,qcondj,qconvj,qd,
     .                   qg,qohmj,qie,qr,qomegaj
      write (nqik, 1155) j,r(j),roa(j),dpedtj,qcondj,qconvj,qd,
     .                   qg,qohmj,qie,qr,qomegaj

  520 continue

c

       CALL flxcal (dflux,drr,hcap,nj,r,ra,anal_eng_flux_e)

      write (nout, 1125)
      write (nout, 1121)
      write (nqik, 1125)
      write (nqik, 1121)
      do 540 j=1,nj
        err    = ss(j) - dpedt(j) - dflux(j)
        j1prt  = ((j-1)/jprt) * jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 540
        qbeamj = qbeame(j) * joupkev
        qrfej  = qrfe(j) * joupkev
        qe2dj  = qe2d(j) * joupkev
        qtfusj = qtfuse(j) * joupkev
        qbfusj = qbfuse(j) * joupkev
        qmagj  = qmag(j) * joupkev
        qsawej = qsawe(j) * joupkev
        if (itimav .eq. 0 .and. w2mix .gt. 0.0)  qmagj  = 0.0
        if (itimav .eq. 0 .and. w2mix .gt. 0.0)  qsawej = 0.0
        err    = err * joupkev
        write (nout, 1155)  j,r(j),roa(j),qbeamj,qrfej,
     .                      qe2dj,qtfusj,qbfusj,qsawej,err,qmagj
        write (nqik, 1155)  j,r(j),roa(j),qbeamj,qrfej,
     .                      qe2dj,qtfusj,qbfusj,qsawej,err,qmagj
  540 continue
c
      do l=1,krf
        if (extqerf(l) .ne. 0.0) then
          write  (nout, 1156)  extqerf_id(l), rf_ext_qetot(l)
          write  (nqik, 1156)  extqerf_id(l), rf_ext_qetot(l)
 1156     format (' qrfe contains the external electron',
     .            ' heat source ', a /
     .            ' with total input heating of ',1pe12.4,' watts')
        end if
      end do
c
c     print out components of integrated electron power sources
c
      const = 4.0*pisq*rmajor*joupkev
      call trap3(r,dpedt,nj,hcap,const,pdpedt)
      call trap3(r,qconde,nj,hcap,const,pconde)
      call trap3(r,qconve,nj,hcap,const,pconve)
      call trap3(r,qdelt,nj,hcap,-const,pdelte)
      call trap3(r,qexch,nj,hcap,-const,pexch)
      call trap3(r,qohm,nj,hcap,const,pohmic)
      call trap3(r,qione,nj,hcap,-const,pione)
      call trap3(r,qrad,nj,hcap,-const,prad)
      call trap3(r,qbeame,nj,hcap,const,pbeame)
      call trap3(r,qbeame_rot,nj,hcap,const,pbeame_rot)
      call trap3(r,qrfe,nj,hcap,const,prfe)
      call trap3(r,qe2d,nj,hcap,const,pe2d)
      call trap3(r,qfuse,nj,hcap,const,pfuse)
      call trap3(r,qtfuse,nj,hcap,const,ptfuse)
      call trap3(r,qbfuse,nj,hcap,const,pbfuse)
      call trap3(r,qmag,nj,hcap,const,pmag)
      call trap3(r,qsawe,nj,hcap,const,psawe)
      call trap3(r,ss,nj,hcap,const,psum)
      call trap3(r,omegale,nj,hcap,const,pomegale)
      write (nout, 1122)
      write (nout, 1121)
      write (nqik, 1122)
      write (nqik, 1121)
      do 542 j=2,nj
      j1prt = ((j-1)/jprt) * jprt
      if (j1prt .ne. j-1 .and. j .ne. nj)  go to 542
      write (nout,1155) j, r(j), roa(j), pdpedt(j), pconde(j), pconve(j)
     .                 ,pdelte(j), pexch(j),pohmic(j), pione(j), prad(j)
     .                 ,-pomegale(j)
      write (nqik,1155) j, r(j), roa(j), pdpedt(j), pconde(j), pconve(j)
     .                 ,pdelte(j), pexch(j),pohmic(j), pione(j), prad(j)
     .                 ,-pomegale(j)
  542 continue
      pradt = ABS (prad(nj))
c
      write (nout, 1123)
      write (nout, 1121)
      write (nqik, 1123)
      write (nqik, 1121)
c


 

      do 544 j=2,nj
        err = psum(j) - pdpedt(j) - pconde(j) - pconve(j)
        j1prt  = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 544
        psawej = psawe(j)
        pmagj  = pmag(j)
        if (itimav .eq. 0 .and. w2mix .gt. 0.0)  psawej = 0.0
        if (itimav .eq. 0 .and. w2mix .gt. 0.0)  pmagj  = 0.0
        write (nout, 1155) j, r(j), roa(j), pbeame(j), prfe(j), pe2d(j),
     .                     ptfuse(j), pbfuse(j), psawej, err, pmagj
        write (nqik, 1155) j, r(j), roa(j), pbeame(j), prfe(j), pe2d(j),
     .                     ptfuse(j), pbfuse(j), psawej, err, pmagj
  544 continue
c
      do j=2,nj
        j1prt = ((j-1)/jprt)*jprt
        if (.not. (j1prt .ne. j-1 .and. j .ne. nj)) then
          write (nout, 1155) j, r(j), roa(j), pbeame(j), pbeame_rot(j)
        end if
      end do
c
c ----------------------------------------------------------------------
c print out components of ion energy source
c ----------------------------------------------------------------------
c
      k = nion + 2    ! k points to ti
      moderun = 'simulation'
      if (itran(k) .eq. 0)  moderun = 'analysis'
      call header (nout, timet, t)
      call header (nqik, timet, t)
      write (nout, 1140) moderun
      write (nout, 1121)
      write (nqik, 1140) moderun
      write (nqik, 1121)
      if (itran(k) .eq. 1) then

        call divflx (flux, fluxb, hcap, r, ra, drr, k, nj, kk, dflux)
      end if
c

      do 620 j=1,nj
c

      dpidt_tot(j) = 0.0

      ss(j) = s(k,j)

      if (itimav .eq. 1 .and. w3mix .gt. 0.0)  ss(j) = ss(j)+qsawi(j)

      sum = 0.0
c

      do 605 i=1,nion
      dnt = 0.0
      if (dtt .ne. 0.0)
     .  dnt = (u(i,j)*u(k,j)-usave(i,j)*usave(k,j))/dtt
      if (itimav .eq. 1)
     .  dnt = (uav(i,j)*uav(k,j)-uav0(i,j)*uav0(k,j))*dtsumi
  605 sum   = sum + dnt

c
      dpidt_tot(j) = 1.5 * sum
c
c --- add terms due to toroidal rotation
c
      if ((iangrot .ne. 0) .and. (angrcple .ne. 0.0)) then



        dtavg = 0.0
        if (   dtt .ne. 0.0)  dtavg = 1.0 / dtt
        if (itimav .eq. 1  )  dtavg = dtsumi
        sum = 0.0

        do 606 i=1,nion
          avg1 = u(i,j)
          avg0 = usave(i,j)
          if (itimav .eq. 1)  avg1 = uav(i,j)
          if (itimav .eq. 1)  avg0 = uav0(i,j)
          if (i .gt. nprim)  go to 606
          dndt = (avg1-avg0)*dtavg*atw(i)
          sum  = sum + dndt
  606   continue

        avg1 = u(nk,j)
        avg0 = usave(nk,j)
        if (itimav .eq. 1)  avg1 = uav (nk,j)
        if (itimav .eq. 1)  avg0 = uav0(nk,j)
        dnt         = (avg1-avg0)*dtavg
        sum         = sum * 0.5 * r2capi(j)*avg1**2*xmassp
        angterm     = amassden(j)*xmassp*r2capi(j)*avg1*dnt
        wdnidt(j)   = sum*angrcple*kevperg
        aniwdwdt(j) = angterm*angrcple*kevperg
        dpidt_tot(j)    = dpidt_tot(j) + wdnidt(j) + aniwdwdt(j)


c
c --- 6.2415064e8 converts from erg/cm**3*sec to keV/(cm**3*sec)
c
      else
        wdnidt  (j) = 0.0
        aniwdwdt(j) = 0.0
      end if
c
      if (itran(k) .eq. 1)  go to 610
      qcondi(j) = qdelt(j)+qexch(j)+qioni(j)-qcx(j)+qbeami(j)
     .          + qrfi(j)+qfusi(j)+qi2d(j)-dpidt_tot(j)-qconvi(j)
     .          + (omegale(j)-qomegapi(j)-qangce(j)+sprcxe(j)
     .          + spreimpe(j)+sprcxree(j))*iangrot
c    NOTE:   qcondi is equivalent to 
c                = s(nion+2,j) -dpidt_tot(j) -qconvi(j)
c                             -qomegapi(j)-qangce(j)    



 


      if (itimav .eq. 1 .and. w3mix .gt. 0.0)  qcondi(j) =
     .                                         qcondi(j) + qsawi(j)
      if (                    w3mix .lt. 0.0)  qcondi(j) =
     .                                         qcondi(j) + qsawi(j)

      !dflux(j) = ss(j) - dpidt_tot(j)
      !dflux(j) = ss(j)
      !dflux(j) = qcondi(j) + dpidt_tot(j) + qconvi(j) ! 88888999
      !dflux(j) = qcondi(j) + dpidt_tot(j) ! 88888999
      !dflux(j) = ss(j) - dpidt_tot(j)
      !dflux(j) = ss(j) - dpidt_tot(j) + qdelt(j)
      dflux(j)  = ss(j) + qdelt(j)
  610 if (j .ne. nj)  go to 615
      !dflux(j)  = ss(j) - dpidt_tot(j)
      !dflux(j)  = ss(j)
      !dflux(j)  = qcondi(j) + dpidt_tot(j) + qconvi(j) ! 88888999
      !dflux(j)  = qcondi(j)+ dpidt_tot(j) ! 88888999
      !dflux(j)  = ss(j) - dpidt_tot(j)
       dflux(j)  = ss(j) + qdelt(j)
  615 qneut(j) = qioni(j) - qcx(j)
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj)  go to 620
      qcv       = qconvi(j) 
      qc        = -qcx(j)
      qpi       = qomegapi(j) * joupkev
      qomegaj   = omegale(j)  * joupkev
      dpidtc(j) = dpidt_tot(j)             ! copy to be saved in sourc.i
      dpidtj    = dpidt_tot(j)    * joupkev
      qcondj    = qcondi(j)   * joupkev
      qcv       = qcv         * joupkev
      qdeltj    = qdelt(j)    * joupkev
      qexchj    = qexch(j)    * joupkev
      qionij    = qioni(j)    * joupkev
      qc        = qc          * joupkev
      write (nout, 1155) j,r(j),roa(j),dpidtj,qcondj,qcv,qdeltj,
     .                   qexchj,qionij,qc ,qpi,qomegaj
      write (nqik, 1155) j,r(j),roa(j),dpidtj,qcondj,qcv,qdeltj,
     .                   qexchj,qionij,qc,qpi,qomegaj
  620 continue
      


       do j=1,nj
c          dflux(j)  = ss(j) + qdelt(j)
           dflux(j)  = ss(j) -dpidt_tot(j)
       ENDDO
       CALL flxcal (dflux,drr,hcap,nj,r,ra,anal_eng_flux_i)
 
c

      write (nout, 1135)
      write (nout, 1121)
      write (nqik, 1135)
      write (nqik, 1121)
c

      do 640 j=1,nj
      qdum(j) = 0.0
      do 630 jb=1,nbeams
      do 630 ic=1,3
  630 qdum(j) = qdum(j) + qb(j,ic,jb)
      qdum(j) = 0.62415064e16*qdum(j) - qbeame(j) - qbeami(j)
      err=ss(j) - dpidt_tot(j) - dflux(j)
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj)  go to 640
      qbeamj = qbeami(j) * joupkev
      qrfij  = qrfi(j) * joupkev
      qi2dj  = qi2d(j) * joupkev
      qtfusj = qtfusi(j) * joupkev
      qbfusj = qbfusi(j) * joupkev
      qsawij = qsawi(j) * joupkev
      if (itimav .eq. 0 .and. w3mix .gt. 0.0)  qsawij = 0.0
      err   =     err * joupkev
      qdumj = qdum(j) * joupkev
      write (nout, 1155)  j,r(j),roa(j),qbeamj,qrfij,qi2dj,
     .                    qtfusj,qbfusj,qsawij,err,qdumj
      write (nqik, 1155)  j,r(j),roa(j),qbeamj,qrfij,qi2dj,
     .                    qtfusj,qbfusj,qsawij,err,qdumj
  640 continue
c
      do l=1,krf
        if (extqirf(l) .ne. 0.0) then
          write  (nout, 1157)  extqirf_id(l), rf_ext_qitot(l)
          write  (nqik, 1157)  extqirf_id(l), rf_ext_qitot(l)
 1157     format (' qrfi contains the external ion heat source ', a /
     .            ' with total input heating of ', 1pe12.4, ' watts')
        end if
      end do


 


c
c print out components of integrated ion power sources

c
      call trap3 (r,dpidt_tot,nj,hcap,const,pdpidt)
      call trap3 (r,qcondi,nj,hcap,const,pcondi)
      call trap3 (r,qioni,nj,hcap,const,pioni)
      call trap3 (r,qcx,nj,hcap,-const,pcx)
      call trap3 (r,qbeami,nj,hcap,const,pbeami)
      call trap3 (r,qbeami_rot,nj,hcap,const,pbeami_rot)
      call trap3 (r,qrfi,nj,hcap,const,prfi)
      call trap3 (r,qfusi,nj,hcap,const,pfusi)
      call trap3 (r,qi2d,nj,hcap,const,pi2d)
      call trap3 (r,qtfusi,nj,hcap,const,ptfusi)
      call trap3 (r,qbfusi,nj,hcap,const,pbfusi)
      call trap3 (r,qsawi,nj,hcap,const,psawi)
      call trap3 (r,ss,nj,hcap,const,psum)
      call trap3 (r,qdum,nj,hcap,const,pdum)
      call trap3 (r,qomegapi,nj,hcap,const,pomegapi)
      call trap3 (r,qangce,nj,hcap,const,pangce)
c
      write (nout, 1144)
      write (nout, 1121)
      write (nqik, 1144)
      write (nqik, 1121)
      sum6      = 0.0
      do 660 j=2,nj
      sum6      = sum6 + trapf (j,r,qconvi,hcap,const)
      pconvi(j) = sum6
      pdelti(j) = -pdelte(j)
      pexchi    = -pexch(j)
      pneut(j)  = pioni(j) + pcx(j)
      j1prt     = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj)  go to 660
      write (nout,1155) j, r(j), roa(j), pdpidt(j), pcondi(j), pconvi(j)
     .      ,pdelti(j), pexchi, pioni(j), pcx(j),pomegapi(j),pomegale(j)
      write (nqik,1155) j, r(j), roa(j), pdpidt(j), pcondi(j), pconvi(j)
     .      ,pdelti(j), pexchi, pioni(j), pcx(j),pomegapi(j),pomegale(j)
  660 continue
c
      write (nout,1145)
      write (nout,1121)
      write (nqik,1145)
      write (nqik,1121)
      do 670 j=2,nj
      err=psum(j)-pdpidt(j)-pcondi(j)-pconvi(j)-pomegapi(j)-pangce(j)
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj)  go to 670
      psawij = psawi(j)
      if (itimav .eq. 0 .and. w3mix .gt. 0.0)  psawij = 0.0
      write (nout, 1155) j, r(j), roa(j), pbeami(j), prfi(j), pi2d(j),
     .                   ptfusi(j), pbfusi(j), psawij, err, pdum(j)
      write (nqik, 1155) j, r(j), roa(j), pbeami(j), prfi(j), pi2d(j),
     .                   ptfusi(j), pbfusi(j), psawij, err, pdum(j)
  670 continue
c


      do j=2,nj
      j1prt = ((j-1)/jprt)*jprt
      if (.not. (j1prt .ne. j-1 .and. j .ne. nj)) then
        write (nout, 1155) j, r(j), roa(j), pbeami(j),
     .                     pbeami_rot(j),pdum(j)
        end if
      end do
c
c     save quantities for yoka data bank
c
      poh =           pohmic(nj)  * 1.0e-6
      pbe =           pbeame(nj)  * 1.0e-6
      pbi =           pbeami(nj)  * 1.0e-6
      prf = (prfe(nj) + prfi(nj)) * 1.0e-6 

c
c ----------------------------------------------------------------------
c --- print out energy terms related to angular rotation
c ----------------------------------------------------------------------
c
      if ((iangrot .ne. 0)  .and. (angrcple .ne. 0.0)) then
         do 718 j=1,nj-1
  718    omegdgam(j) = qomegapi(j) * joupkev + vischeat(j)
*%%%%    write (6, *) '%%%% in subroutine PSOURC, calling EXTRAP #1'
         call extrap (ra(nj-2),ra(nj-1),r(nj),vischeat(nj-2),
     .               vischeat(nj-1),vischeat(nj))
*%%%%    write (6, *) '%%%% in subroutine PSOURC, calling EXTRAP #2'
         call extrap (ra(nj-2),ra(nj-1),r(nj),omegdgam(nj-2),
     .               omegdgam(nj-1),omegdgam(nj))
*%%%%    write (6, *) '%%%% in subroutine PSOURC, calling HEADER #7'
         call header (nout, timet, t)
         write  (nout, 1190)
 1190    format (1x, 33(1h-), 'ION ENERGY SOURCES DUE TO ',
     .                        'ANGULAR ROTATION, (W/CM**3)', 33(1H-))
         write (nout,1191)
 1191    format (4x,'j',4x,'r',3x,'r/a',4x,'wdnidt',5x,'niwdwdt',
     .    5x,'omegdgam',4x,'vischeat',5x,'qangce',
     .    4x,'omegale',6x,'th cx  ',4x,'rec +fcx',
     .    4x,'e-impact')
        write (nout, 1121)
        do 721 j=1,nj
          j1prt = ((j-1)/jprt)*jprt
          if (j1prt .ne. j-1 .and. j .ne. nj)  go to 721
          wdnidtj  = wdnidt(j)*joupkev
          andwdtj  = aniwdwdt(j)*joupkev
          qangcj   = qangce(j)*joupkev
          sprcxej  = sprcxe(j)*joupkev
          sprcxrej = sprcxree(j)*joupkev
          spreimpj = spreimpe(j)*joupkev
          qomegaj  = omegale(j)*joupkev
          write (nout, 1155)  j,r(j),roa(j),wdnidtj,andwdtj,omegdgam(j),
     .              vischeat(j),qangcj,qomegaj,sprcxej,sprcxrej,spreimpj
  721   continue
c
c --- calculate integrated sources (watts) and vol avg. (w/cm**3)
c
         const2 = 4.0 * pisq*rmajor
         const3 = const2*joupkev
         call trap3(r,wdnidt,nj,hcap,const3,pwdnidt)
         call trap3(r,aniwdwdt,nj,hcap,const3,pniwdwdt)
         call trap3(r,qangce,nj,hcap,const3,pqangce)
         call trap3(r,sprcxe,nj,hcap,const3,psprcxe)
         call trap3(r,sprcxree,nj,hcap,const3,psprcxee)
         call trap3(r,spreimpe,nj,hcap,const3,pspreimp)
         call trap3(r,vischeat,nj,hcap,const2,pvscheat)
         call trap3(r,omegdgam,nj,hcap,const2,pomdgam)
          avwdnidt = pwdnidt(nj)/volume
          avwdwdt  = pniwdwdt(nj)/volume
          avmdgam  = pomdgam(nj)/volume
          avscheat = pvscheat(nj)/volume
          avangce  = pqangce(nj)/volume
          avomegle = pomegale(nj)/volume
          avprcxe  = psprcxe(nj)/volume
          avprcxee = psprcxee(nj)/volume
          avpreimp = pspreimp(nj)/volume
          write (nout, 1197) avwdnidt,avwdwdt,avmdgam,
     .      avscheat,avangce,avomegle,avprcxe,avprcxee,avpreimp
 1197   format ('  VOLUME AVERAGE',1p9e12.3)



 
c
c --- print out integrated sources
c
*%%%%   write (6, *) '%%%% in subroutine PSOURC, calling HEADER #8'
        call header (nout, timet, t)
        write  (nout, 1194)
 1194   format (1x,33(1h-),'INTEGRATED ION ENERGY SOURCES DUE TO',
     .                    ' ANGULAR ROTATION, (W)',33(1H-))
        write (nout,1196)
 1196    format (4x,'j',4x,'r',3x,'r/a',4x,'pwdnidt',4x,'pniwdwdt',
     .    4x,'pomegdgam', '  pvischeat',5x,'pqangce',
     .    3x,'pomegale',5x,'p-th cx ',3x,'p-rec+fcx',
     .    3x,'p-eimpact')
        write (nout,1121)
        do 722 j=1,nj
          j1prt = ((j-1)/jprt)*jprt
          if (j1prt .ne. j-1 .and. j .ne. nj)  go to 722
          write (nout, 1155)  j,r(j),roa(j),pwdnidt(j),pniwdwdt(j),
     .    pomdgam(j),pvscheat(j),pqangce(j),pomegale(j),psprcxe(j),
     .    psprcxee(j),pspreimp(j)
  722   continue
      else
        do 723 j=1,nj
          pwdnidt (j) = 0.0
          pniwdwdt(j) = 0.0
          pqangce (j) = 0.0
          psprcxe (j) = 0.0
          psprcxee(j) = 0.0
          pspreimp(j) = 0.0
          pvscheat(j) = 0.0
          pomdgam (j) = 0.0
  723     omegdgam(j) = 0.0
      end if
c
c ----------------------------------------------------------------------
c print out charge balance, i.e., Faraday's law
c ----------------------------------------------------------------------
c
      k       = nion+3
      moderun = 'simulation'
      if (itran(k) .eq. 0)  moderun = 'analysis'
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling HEADER #9'
      call header (nout, timet, t)
      write (nout, 1143)  twkfar, moderun
      call divflx(flux, fluxb, hcap, r, ra, drr, k, nj, kk, dflux)
      do 720 j=1,nj
      dudt(j)  = 0.0
      dflux(j) = r(j)*dflux(j)
      ss(j)    = s(nk-iangrot,j)*r(j)
      if (j .eq. 1)  go to 715
****  if (n .eq. 0)  go to 710
      if (dtt .ne. 0.0)
     .  dudt(j) = (u(k,j)-usave(k,j))/
     .            (r(j)*fcap(j)*gcap(j)*hcap(j)**2*dtt*twkfar)
      if (itimav   .eq. 1)  dudt(j) = (uav(k,j)-uav0(k,j))*dtsumi /
     .                          (r(j)*fcap(j)*gcap(j)*hcap(j)**2*twkfar)
      if (itran(k) .eq. 1)  go to 710
      dflux(j) = ss(j) - dudt(j)
  710 if (j .ne. nj)  go to 715
      dflux(j) = ss(j) - dudt(j)
  715 j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj)  go to 720
      err=ss(j) - dudt(j) - dflux(j)
      write (nout, 1155)  j,r(j),roa(j),dudt(j),dflux(j),ss(j),err
  720 continue
c
      call trapv (r,dudt,fact,nj,dua)
      call trapv (r,dflux,fact,nj,dfa)
      call trapv (r,ss,fact,nj,ssa)
      err=ssa - dua - dfa
      write (nout, 1152)  dua, dfa, ssa, err
c
c ----------------------------------------------------------------------
c compute various sums of neutral beam injection sources
c (pfil and sthru are for Yoka data bank)
c ----------------------------------------------------------------------
c
      pbap   = 0.0
      pbsap  = 0.0
      pwall  = 0.0
      pborb  = 0.0
      pbplaf = 0.0
      pbplas = 0.0
      pbel   = 0.0
      pbion  = 0.0
      pfil   = 0.0
      sthru  = 0.0


      GO TO 850 


      ! obsolete code :
      beam_model:    IF(use_P_Nfreya)THEN
        !   CALL P_Nfreya_output(timet, t)
      !IF (ibeam .lt. 2)  go to 850
      ELSEIF (ibeam .GE. 2) THEN beam_model
      do 810 jb=1,nbeams
      do 810 ic=1,3
      pbap   = pbap  + pbeam(ic,jb)
      pbsap  = pbsap + fap(ic,jb)*pbeam(ic,jb)
      pwall  = pwall + fwall(ic,jb)*pbeam(ic,jb)
      pborb  = pborb + forb(ic,jb)*pbeam(ic,jb)
      pbeamf(ic,jb) = pbeam(ic,jb)*
     .                (1.0-fap(ic,jb)-fwall(ic,jb)-forb(ic,jb))
      pbplaf = pbplaf + pbeamf(ic,jb)
      fpe(ic,jb)  = 0.0
      fpi(ic,jb)  = 0.0
      fpcx(ic,jb) = 0.0
      do 805 j=1,nj
      psum(j) = qbbe(j,ic,jb)
  805 pdum(j) = qbbi(j,ic,jb)
      call trapv (r,psum,hcap,nj,fpe(ic,jb))
      call trapv (r,pdum,hcap,nj,fpi(ic,jb))
      call trapv (r,qb(1,ic,jb),hcap,nj,pbeams(ic,jb))
      pbeams(ic,jb) = volfac*pbeams(ic,jb)
      if (pbeams(ic,jb) .le. 0.0)  go to 810
      fpe(ic,jb)  = volfac*fpe(ic,jb)/pbeams(ic,jb)
      fpi(ic,jb)  = volfac*fpi(ic,jb)/pbeams(ic,jb)
      fpcx(ic,jb) = 1.0-fpe(ic,jb)-fpi(ic,jb)
      pbplas      = pbplas + pbeams(ic,jb)
      pbel        = pbel   + fpe(ic,jb)*pbeams(ic,jb)
      pbion       = pbion  + fpi(ic,jb)*pbeams(ic,jb)
      pfil        = pfil   + fpcx(ic,jb)*pbeams(ic,jb)
  810 continue
      fsap   = 0.0
      fw     = 0.0
      florb  = 0.0
      if (pbap .eq. 0.0)  go to 812
      fsap   = pbsap/pbap
      fw     = pwall/pbap
      florb  = pborb/pbap
  812 fpbe   = 0.0
      fpbi   = 0.0
      fpbcx  = 0.0
      if (pbplas .eq. 0.0)  go to 814
      fpbe   = pbel/pbplas
      fpbi   = pbion/pbplas
      fpbcx  = pfil/pbplas
  814 ptor   = pbap - pbsap
      if (ptor .gt. 0.0)  sthru = (pwall/ptor) * 100.0
      pfil   = pfil * 1.0e-6
c
c print out neutral beam injection sources
c
      do 830 jb=1,nbeams
      do 830 ic=1,3
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling HEADER #10'
      call header (nout, timet, t)
      write (nout, 1146)  jb, ic
      do 820 j=1,nj
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj)  go to 820
      write (nout,1153) j,r(j),hibr(j,ic,jb),hdep(j,ic,jb),
     .   zeta(j,ic,jb), qbsav(j,ic,jb), qb(j,ic,jb),
     .   fbe(j,ic,jb), fbi(j,ic,jb), taupb(j,ic,jb), taueb(j,ic,jb)
  820 continue
      write (nout,8110)         ebeam(ic,jb)
      write (nout,8120)         bion(ic,jb)
      write (nout,8130)         bneut(ic,jb)
      write (nout,8140) pbap,   pbeam(ic,jb)
      write (nout,8150) fsap,   fap(ic,jb)
      write (nout,8160) fw,     fwall(ic,jb)
      write (nout,8161) florb,  forb(ic,jb)
      write (nout,8162) pbplaf, pbeamf(ic,jb)
      write (nout,8163) pbplas, pbeams(ic,jb)
      write (nout,8164) fpbe,   fpe(ic,jb)
      write (nout,8166) fpbi,   fpi(ic,jb)
      write (nout,8168) fpbcx,  fpcx(ic,jb)
      write (nout,8170)
      fxorb = forb(ic,jb)/(1.0-fap(ic,jb)-fwall(ic,jb))
      write (nout,8180) fb11(ic,jb), wb11(ic,jb),
     .    fb10(ic,jb), wb10(ic,jb), fb01(ic,jb), wb01(ic,jb),
     .    fb00(ic,jb), wb00(ic,jb), fxorb, fber(ic,jb)
  830 continue
c
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling HEADER #11'
      call header (nqik, timet, t)
      write (nqik,8200)
      write (nqik,8110)        (( ebeam(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8120)        ((  bion(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8130)        (( bneut(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8140) pbap,  (( pbeam(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8150) fsap,  ((   fap(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8160) fw,    (( fwall(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8161) florb, ((  forb(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8162) pbplaf,((pbeamf(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8163) pbplas,((pbeams(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8164) fpbe,  ((   fpe(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8166) fpbi,  ((   fpi(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8168) fpbcx, ((  fpcx(ic,jb),ic=1,3),jb=1,nbeams)
      write (nqik,8170)
      do 840 jb=1,nbeams
      write (nqik,8210) jb,(ebeam(ic,jb),ic=1,3)
      do 840 j=1,nj
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj)  go to 840
      write (nqik,8220) j,r(j),roa(j),
     .         (qb(j,ic,jb),fbe(j,ic,jb),fbi(j,ic,jb),ic = 1,3)
  840 continue
c  850 continue
      ENDIF beam_model




 850   CONTINUE
    
       CALL Nfreya_output(timet,t)


c
c ----------------------------------------------------------------------
c calculate particle diffusion coefficient, electron thermal
c conductivity, and ion thermal conductivity if appropriate
c variable in =1 if neutral corresponds to primary ion species 1
c              2                                               2
c variable in is from neut.i INCLUDE file
c note that dinv was calculated above for the case itran(in) = 0
c (i.e., analysis mode)
c ----------------------------------------------------------------------
c
      k = in
      if (itran(k) .eq. 0)  go to 864
      do 862 j=1,nj-1
  862 dinv(j) = d(k,k,j)
     .+d(k,nion+1,j)*(te(j)-te(j+1))/(en(j,k)-en(j+1,k))
      go to 868
  864 do 865 j=1,nj
  865 pdum(j) = dinv(j)
      call zeroa (dinv,nj)
      call flxcal (pdum,drr,hcap,nj,r,ra,swork)
      do 866 j=1,nj-1
      grad = (theta*(u(k,j+1)-u(k,j))+onemt*(usave(k,j+1)-usave(k,j)))
     .       /dr(j)
      if (grad .eq. 0.0)  go to 866
      ena = 0.5 * theta*(u(in,j)+u(in,j+1))
     .    + 0.5 * onemt*(usave(in,j)+usave(in,j+1))
      dinv(j) = (-swork(j)+ena*vpinch(j))/grad
  866 continue
  868 continue
c
      k = nion + 1    ! te
      call zeroa (chietrinv,nj)
      if (itran(k) .eq. 0)  go to 874
      do 872 j=1,nj-1
          detemp = te(j+1) -te(j)  !HSJ 3/26/03
          if(detemp .ne. 0.0)then
            xkeinv(j) = d(k,k,j)
     .      +d(k,in,j)*(en(j,in)-en(j+1,in))/detemp
          else
             xkeinv(j) =0.0
          endif
          grad = (te(j+1)-te(j))/dr(j)
          gconde(j) = -d(k,k,j)*grad*joupkev
          gconve(j) =conve(j)*joupkev
 872      gconvde(j) = gconde(j) + gconve(j)
      go to 878
  874 call zeroa (xkeinv,nj)
      call zeroa (gconde,nj-1)
      call flxcal (qconde,drr,hcap,nj,r,ra,gconde)
      call zeroa (gconve,nj-1)
      call flxcal (qconve,drr,hcap,nj,r,ra,gconve)
      do 876 j=1,nj-1
        grad = (theta*(u(k,j+1)-u(k,j)) + onemt*(usave(k,j+1)
     .         -usave(k,j)))/dr(j) ! space- and time-centered derivative
        gconde(j) =  gconde(j) * joupkev    ! conductive heat flux
        gconve(j) =  gconve(j) * joupkev    !  in watts/cm**2
        gconvde(j) = gconde(j) + gconve(j)
        if (grad .eq. 0.0)  go to 876
        ene_theta=ene(j)
        if (itimav .eq. 0)
     .             ene_theta = theta*ene(j)+onemt*enesav(j)
        if (eneav .ne. 0.0)
     .             chietrinv(j) = -gconvde(j)/(grad*eneav*joupkev)
        xkeinv(j) = -gconde(j) / (grad*joupkev)
  876 continue
  878 continue
c
      k = nion + 2    ! ti
      call zeroa (chiitrinv,nj)  ! transport chi,(sum of cond+conv)
      if (itran(k) .eq. 0)  go to 884
      do 882 j=1,nj-1
          ditemp = ti(j+1) -ti(j)  !HSJ 3/26/03
          if(ditemp .ne. 0.0)then
               xkiinv(j) = d(k,k,j)
     .          +d(k,in,j)*(en(j,in)-en(j+1,in))/ditemp
          else
               xkiinv(j) =0.0
          endif
          grad = (ti(j+1)-ti(j))/dr(j)
          gcondi(j) = -d(k,k,j)*grad*joupkev
          gconvi(j) =convi(j)*joupkev
 882      gconvdi(j) = gcondi(j) + gconvi(j)

      go to 888
c     determine the power balance chi-ion if Ti is in analysis mode:
  884 call zeroa (xkiinv, nj)      ! see notes, vol 2, pg 32 for details
      call zeroa (gcondi,nj-1)
      call flxcal (qcondi,drr,hcap,nj,r,ra,gcondi)
      call zeroa (gconvi,nj-1)
      call flxcal (qconvi,drr,hcap,nj,r,ra,gconvi)
      do 886 j=1,nj-1
        grad = (theta*(u(k,j+1)-u(k,j)) +
     .                  onemt*(usave(k,j+1)-usave(k,j)))/dr(j)
        gcondi(j) =  gcondi(j) * joupkev    ! conductive heat flux
c                                           in watts/cm**2
        gconvi(j) =  gconvi(j) * joupkev
        gconvdi(j)=gconvi(j)+gcondi(j)
        if (grad .eq. 0.0)  go to 886
           enps = 0.0
           do jprim=1,nprim
             if (itimav .eq. 0) then
               enps = enps + theta* u(jprim,j)
     .              + onemt*usave(jprim,j)
             else
               enps = enps + u(jprim,j)
             end if
           end do
           xkiinv(j) = -gcondi(j) / (grad*joupkev)
           chiitrinv(j) = -gconvdi(j)/ (joupkev*grad*enps)
  886 continue
  888 continue
c
c --- print electron and ion (diagonal component) conductive heat flux
c
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling HEADER #12'
      call header (nout, timet, t)
      write  (nout, 1179)  itranflag(nion+1), itranflag(nion+2)
 1179 format (20x, 'electron and ion DIAGONAL conductive ',
     .             'heat flux w/cm**2' //
     .              4x,'j',9x,'r',7x,'conde', 8x,'condi' /
     .             22x,a,10x,a)
          jp = '.5'
          do j=1,nj-1
            j1prt = ((j-1)/jprt)*jprt
            if (j1prt .ne. j-1 .and. j .ne. nj-1)  go to 1177
            write  (nout, 1178)  j, jp, ra(j), gconde(j), gcondi(j)
 1178       format (1x, i4, a2, f9.2, 2(1pe12.3))
 1177       continue
          end do
c
c print out diffusivity, thermal conductivities, etc.
c on half grid (in rho) at central time
c
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling HEADER #13'
      call header (nout, timet, t)
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling HEADER #14'
      call header (nqik, timet, t)
      write (nout,1160) in,namep(1),itranflag(in),itranflag(nion+1),
     .            itranflag(nion+2),itranflag(nion+1),itranflag(nion+2)
      write (nqik,1160) in,namep(1),itranflag(in),itranflag(nion+1),
     .            itranflag(nion+2),itranflag(nion+1),itranflag(nion+2)
      jp = '.5'
      do 751 j=1,nj-1
      eeneav(j) = 0.5 * theta*(ene(j)+ene(j+1))
     .          + 0.5 * onemt*(enesav(j)+enesav(j+1))
      xndinv(j) = dinv(j)*eeneav(j)
      chieinv(j) =  xkeinv(j)/eeneav(j)
      gr2a = (grho2_mesh(j+1)+grho2_mesh(j))*0.5
      gr1a = (grho1_mesh(j+1)+grho1_mesh(j))*0.5
      if(gr2a .eq. 0.0)call stop('gr2a problem',0)
      chieinv_transp(j) = chieinv(j)*gr1a/gr2a
c     Transp uses the chienv of Onetwo times <(grad rho)**2>/<grad rho>
      enpav(j) = 0.0
      do 752 jprim=1,nprim
      enpav(j) = enpav(j) + 0.5 * (en(j,jprim)+en(j+1,jprim))
  752 continue
      chiinv(j) = 0.0
      if (enpav(j) .ne. 0.0)  chiinv(j) = xkiinv(j)/enpav(j)
      chiiinv_transp(j) = chiinv(j)*gr1a/gr2a
c     Transp uses the chiinv of Onetwo times <(grad rho)**2>/<grad rho>
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj-1)  go to 751
      write  (nout, 1200)  j, jp, ra(j), dinv(j), chieinv(j), chiinv(j),
     .                     xkeinv(j),xkiinv(j),xnuse(j),xnus(1,j),
     .                     ftrap(j),eta(j)*8.98755179e11
      write  (nqik, 1200)  j, jp, ra(j), dinv(j), chieinv(j), chiinv(j),
     .                     xkeinv(j),xkiinv(j),xnuse(j),xnus(1,j),
     .                     ftrap(j),eta(j)*8.98755179e11
 1200 format (1x, i4, a2, f9.2, 9(1pe12.3))
  751 continue
c
c ----------------------------------------------------------------------
c --- calculate single fluid chi,xeffsf,on the half grid
c ----------------------------------------------------------------------
c
      do j=1,nj-1
          enea = eeneav(j)
          ensum = enpav(j)
          teprime = (theta*(u(nion+1,j+1)-u(nion+1,j)) +
     .        onemt*(usave(nion+1,j+1)-usave(nion+1,j)))/dr(j)
          tiprime = (theta*(u(nion+2,j+1)-u(nion+2,j)) +
     .        onemt*(usave(nion+2,j+1)-usave(nion+2,j)))/dr(j)
          denom =enea*teprime+ensum*tiprime
          if (denom .ne. 0.0) then
          xeffsf(j) = (chieinv(j)*enea*teprime+chiinv(j)*ensum*tiprime)
     .              / denom
          else
             xeffsf(j) = 0.0
          end if
      end do
c
c  convert xndinv, xkeinv, xkiinv, and xkineo from mesh centers
c  to mesh points for output to plot file
c
      call mescon (xndinv  , dr, nj)
      call mescon (xkeinv  , dr, nj)
      call mescon (xkiinv  , dr, nj)
      call mescon (xkineo  , dr, nj)
      call mescon (chiineo , dr, nj)
      call mescon (xkeneo  , dr, nj)
      call mescon (xkangrot, dr, nj)
c
c  calculate and print diagnostic transport coefficients
c
      write (nout, 1201)
      write (nqik, 1201)
c
      do 215 j=1,nj-1
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj-1)  go to 215
        zet    = CMPLX (one_third, zetaim(j))
        call zfunc(zet,zz,zp)
        xkecha = xkepp(j)*
     .           ABS ((1.0 + zet*zz) / (1.0 + (0.0,1.0)*zetaim(j)*zz))
        xkirat = 0.0
        if (xkineo(j) .ne. 0.0)  xkirat = wneo(3,3)*xkiinv(j)/xkineo(j)
        write (nout, 1200)  j, jp, ra(j),xkecar(j),xkeohk(j),
     .                      xkepp(j),xkedom(j),xkirat,
     .                      shearp(j),slene(j),slte(j),chiwneo(j)
        write (nqik, 1200)  j, jp, ra(j),xkecar(j),xkeohk(j),
     .                      xkepp(j),xkedom(j),xkirat,
     .                      shearp(j),slene(j),slte(j),chiwneo(j)
  215 continue
c
c     convert scale lengths from mesh centers to mesh points
c     for output to the yokfil file (which is OBSOLETE)
c
      call mescon (slene , dr, nj)
      call mescon (slte  , dr, nj)
      call mescon (slti  , dr, nj)
      call mescon (slpres, dr, nj)
      call mescon (shearp, dr, nj)
c
c ----------------------------------------------------------------------
c --- more diagnostic coefficients
c --- central time, on the full grid (in rho)
c ----------------------------------------------------------------------
c
      do 3005 j=1,nprim
        call mescon (etaim   (1,j), dr, nj)
 3005   call mescon (xkwmatdm(1,j), dr, nj)
      if (nprim .ge. 2)
     .call mescon (xkwmatdm(1,nprim+1),dr,nj)
      call mescon (xkematdm,dr,nj)
      call mescon (vionzgrd,dr,nj)
      call mescon (chiwneo,dr,nj)
      call mescon (chiwmatm,dr,nj)
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling HEADER #15'
      call header (nout, timet, t)
      write (nout, 3010)
 3010 format (2x,40(1h-),'ELECTROSTATIC DRIFT WAVE RELATED QUANTITIES',
     .       40(1H-))
      do 3020 j=1,nj
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 3020
        eeneav(j) = theta*ene(j)+onemt*enesav(j)
        chieinvd = xkeinv(j)/eeneav(j)
        enpav(j) = 0.0
        do k=1,nprim
          enpav(j) = enpav(j)+en(j,k)
        end do
        chiinvd = xkiinv(j)/enpav(j)
        if ((nprim .eq. 1) .and. (j .eq. 1))
     .  write  (nout, 3015)
 3015   format (4x,'j',4x,'r',3x,'r/a',6x,'chie',3x,
     .         'chiematdm',6x,'chiwneo', '  chiimatdm',3x,'chi-torrot',
     .          3x,'vionzgrd',5x,'eta-ion' /
     .          7x,'(cm)',5x,5(1x,'(cm**2/s)',2x),3x,'1/sec')
        if ((nprim .ge. 2) .and. (j .eq. 1))
     .  write  (nout, 3016)  namep(1), namep(2)
 3016   format (4x,'j',4x,'r',3x,'r/a',1x,'chiematdm',6x,
     .         'chiwneo', '  chiimatdm',3x,'chiimatdm',2x,
     .         'chiimatdmtot',1x,'chi-torrot',3x,'vionzgrd',
     .           5x,'eta-ion',5x,'chiwneo+md' /
     .          42x,a7,5x,a7 /
     .           7x,'(cm)',5x,6(1x,'(cm**2/s)',2x),3x,'1/sec',
     .          15x,'(cm**2/s)')
      if (nprim .eq. 1)
     .  write (nout, 1155)  j, r(j), roa(j), chieinvd, xkematdm(j),
     .  chiinvd, xkwmatdm(j,1), xkangrot(j), vionzgrd(j), etaim(j,1)
      if (nprim .ge. 2)
     .  write (nout, 1155)  j, r(j), roa(j), xkematdm(j), chiwneo(j),
     .  xkwmatdm(j,1), xkwmatdm(j,2), xkwmatdm(j,nprim+1), xkangrot(j),
     .  vionzgrd(j), etaim(j,1), chiwmatm(j)
 3020 continue
c
      ifsflag = INT (fs(1))     ! flow shear turbulence suppression flag
c
c ----------------------------------------------------------------------
c print out Rebut-Lallia and single fluid models
c we convert everything to the full grid first, to be consistent
c with plot output (to unit ntrplt; see below)
c ----------------------------------------------------------------------
c
      call mescon (crit_grad,dr,nj)
      call mescon (grad_te,dr,nj)
      grad_te(1) = 0.0                              ! boundary condition
      call mescon (xchierl,dr,nj)
      call mescon (xchiirl,dr,nj)
      call mescon (xeffrl,dr,nj)
      call mescon (xeffsf,dr,nj)
      call print_rebut(r,nj,jprt,nout,t,ifsflag)
c
c ----------------------------------------------------------------------
c print out Hsieh and single-fluid models
c we convert everything to the full grid first, to be consistent
c with plot output (to unit ntrplt, see below)
c ----------------------------------------------------------------------
c
      if (wshay .ne. 0.0) then
c
c --- convert to full grid
c
        call mescon (schie,dr,nj)
        call mescon (schii,dr,nj)
        call mescon (scheff,dr,nj)
        call mescon (xkeshay,dr,nj)
        call mescon (xkishay,dr,nj)
      end if
c
      if (wshay .ne. 0.0 .and. .not. scsmult) then
c
c --- standard print routine for Hsieh model,used when constant
c --- multipliers are known a priori (i.e., from values given in inone)
c
      call print_shay (nout,jprt,r,timeshay,t,xkeshay,xkishay,
     .                 xkeinv,xkiinv,xkeneo,xkineo,
     .                 psir,rmajorvec,grad_te,nj,te,smult,skimult,
     .                 snexp,sbpexp,srexp,sbigrexp,stexp,sdtdrexp,
     .                 sdenscale,irinshay,iroutshay,ifsflag)
      else if (wshay .ne. 0.0 .and. scsmult) then
c
c --- special least squares calculation and print routine used when
c --- the constant multipliers in Hsieh model are to be found
c
            call shay_multiplier (srincev,sroutcev,nout,r,timeshay,
     .                 t,xkeshay,xkishay,xkeinv,xkiinv,xkeneo,xkineo,
     .                 psir,rmajorvec,xdum,ydum,zdum,sdum,vdum,
     .                 nj,smult,skimult,smultstder,skimultstder,te,
     .                 slim95,skilim95,irinshaycev,iroutshaycev,skimass,
     .                 snexp,sbpexp,srexp,sbigrexp,stexp,sdtdrexp,
     .                 sdenscale,irinshay,iroutshay,smulta,skimulta,
     .                 maxtimes,ktimes,ifsflag)
      end if
c
c --- Staebler-Hinton shear flow turbulence suppression factors
c
      if (ifsflag .gt. 0)  then
        call print_staebler (nout, t)
c
c       convert sperp to full grid for plotting
c
        call mescon (sperp, dr, nj)
      end if
c
c --- phenomenological shear flow turbulence suppression factors
c
      if (irgc .ne. 0)  call print_rgcmodel (nout, t)
c
c --- Weiland model output  HSJ 2/22/96
c
      if (include_weiland .eq. 1) then
        call mescon (xchie_weiland,dr,nj) ! convert half to full grid
        call mescon (xchii_weiland,dr,nj)
        call mescon (xeffwl,dr,nj)
        call mescon (xke_weiland,dr,nj)
        call mescon (xki_weiland,dr,nj)
        call mescon (d_weiland,dr,nj)
        call print_weiland(r,nj,jprt,nout,t,ifsflag)
      end if
c
c --- IFS model output  HSJ 2/11/97
c
      if (include_ifs .eq. 1) then
c
c       convert the ip_chi2 results to  the full grid for
c       output and plotting purposes . Extrapolation done by
c       mescon on the end points may cause trouble:
c
        call mescon (chi_i_ifs, dr, nj)
        chi_i_ifs(1 ) = MAX (zero, chi_i_ifs(1 ))
        chi_i_ifs(nj) = MAX (zero, chi_i_ifs(nj))
        call mescon (chi_e_ifs, dr, nj)
        chi_e_ifs(1 ) = MAX (zero, chi_e_ifs(1 ))
        chi_e_ifs(nj) = MAX (zero, chi_e_ifs(nj))
        call mescon (RLTcrit_ifs , dr, nj)
        call mescon (RLTcritz_ifs, dr, nj)
        call print_ifs (r, nj, jprt, nout, t, ifsflag)
      end if
c
      call mescon (xchietot,dr,nj)       ! convert half to full grid
      call mescon (xchiitot,dr,nj)
      call mescon (xketot,dr,nj)
      call mescon (xkitot,dr,nj)
      call mescon (xdchitot,dr,nj)
      call mescon (xkangtot,dr,nj)      
      call mescon (xkangrot,dr,nj)         !new line  HSJ 10/16/01
c
c ----------------------------------------------------------------------
c calculate and print out confinement times
c ----------------------------------------------------------------------
c
      const1 = 1.5 * const
      do j=1,nj
        energe(j) = theta*ene(j)*u(nion+1,j)
     .            + onemt*enesav(j)*usave(nion+1,j)
        energi(j) = 0.0
        do k=1,nion
          energi(j) = energi(j) + theta*u(k,j)*u(nion+2,j)
     .              + onemt*usave(k,j)*usave(nion+2,j)
        end do
      end do
c
c ----------------------------------------------------------------------
c --- get the local angular momentum conf. time.
c --- (theta weighting of torque terms is neglected)
c ----------------------------------------------------------------------
c
      if (iangrot .ne. 0) then
        do 795 j=1,nj
          avg0 = 0.0
          avg1 = 0.0
          do 796 k=1,nion
            avg0 = avg0+usave(k,j)*usave(nk,j)*atw(k)
            avg1 = avg1+u(k,j)*u(nk,j)*atw(k)
  796     continue
          avg0 = avg0*r2capi0(j)*xmassp
          avg1 = avg1*r2capi(j)*xmassp
          avgangmt(j) = onemt*avg0+theta*avg1
          avgangdt(j) = 0.0
          if (dtt .ne. 0.0)  avgangdt(j) = (avg1-avg0)/dtt
          xxxx = storqueb(j) - avgangdt(j)
          tauang(j) = 0.0
          if (xxxx .ne. 0.0)  tauang(j) = avgangmt(j)/xxxx
  795   continue
      end if
      sume      = 0.0
      sumi      = 0.0
      suma      = 0.0
      sumb      = 0.0
      ataue(1)  = tauec
      tauer(1)  = taueec
      tauelc(1) = tauec
c
      do j=2,nj
c
c  local confinement time tauelc
c
        wnumer = 1.5*(energe(j)+energi(j))+walp(j)+wbeam(j)
        qdenom = qohm(j)+qbeame(j)+qbeami(j)+qfuse(j)+qfusi(j)+qrfe(j)+
     .           qrfi(j)-dpedt(j) - dpidt_tot(j)
        tauelc(j) = wnumer / MAX (qdenom, ten_to_minus_20th)
        sume      = sume + trapf(j,r,energe,hcap,const1)
        sumi      = sumi + trapf(j,r,energi,hcap,const1)
        suma      = suma + trapf(j,r,walp,hcap,const)
        sumb      = sumb + trapf(j,r,wbeam,hcap,const)
        ataue(j)  = (sume + sumi + suma + sumb)
     .           / (pohmic(j)+pbeame(j)+pbeami(j)+pfuse(j)+pfusi(j)
     .            + prfe(j)+prfi(j)-pdpedt(j)-pdpidt(j))
        tauetr(j) = sume / MAX (pconde(j)+pconve(j), ten_to_minus_10th)
        tauitr(j) = sumi / MAX (pcondi(j)+pconvi(j), ten_to_minus_10th)
        denom     = MAX (pohmic(j)+pbeame(j)+pfuse(j)+prfe(j),
     .                   ten_to_minus_10th)
        tauer(j)  = sume / denom
        denom     = MAX (pdelti(j)+pbeami(j)+pfusi(j)+prfi(j),
     .                   ten_to_minus_10th)
        tauir(j)  = sumi / denom
      end do
c
      call extrap (r(2),r(3),r(1),tauetr(2),tauetr(3),tauetr(1))
      call extrap (r(2),r(3),r(1),tauitr(2),tauitr(3),tauitr(1))
      call extrap (r(2),r(3),r(1),tauir(2),tauir(3),tauir(1))
      write (nout,1170)
      write (nqik,1170)
      do 792 j=1, nj
      j1prt = ((j-1)/jprt)*jprt
      if (j1prt .ne. j-1 .and. j .ne. nj)  go to 792
      tres  = 0.0
      if (eta(j) .ne. 0.0)  tres = 1.396e-20*r(nj)**2/eta(j)
      write (nout, 1155) j,r(j),roa(j),ataue(j),tauetr(j),tauitr(j),
     .                   tauer(j),tauir(j),taupe(j),tauelc(j),tres,
     .                   tauang(j)
      write (nqik, 1155) j,r(j),roa(j),ataue(j),tauetr(j),tauitr(j),
     .                   tauer(j),tauir(j),taupe(j),tauelc(j),tres,
     .                   tauang(j)
  792 continue
c
c print out debug variable list calculated in diffus
c
*%%%% write (6, *) '%%%% in subroutine PSOURC, calling HEADER #16'
      call header (nout, timet, t)
      write (nout, 1161)
      jp = '.5'
c
      do 793 j=1,nj-1
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 793
        realjp = 0.5
        write (nout, 1155)  j, realjp, ra(j), (ydebug(j,i), i=1,8)
  793 continue
c
      write (nout, 1162)
c
      do 794 j=1,nj-1
        j1prt = ((j-1)/jprt)*jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 794
        realjp = 0.5
        write (nout, 1155)  j, realjp, ra(j), (ydebug(j,i), i=1,8)
  794 continue
c
c ----------------------------------------------------------------------
c output for plotting (since psourc is called by out which in turn
c is called whenever printed output is to be generated,
c ie at the prtlst times, this has the effect of giving us plots
c at the prtlst times as well,the pltlst of times for
c plotting has no effect here) HSJ
c ----------------------------------------------------------------------
c

 

      call exptlprofiles (timet)       ! get exptl profiles if available

      write (ntrplt, '(a)')  ' **** continue ****'
    
      iflag = 4
      write (ntrplt, '(i10)') iflag

      write (ntrplt, 9100) timet

c 9100 format (1p6e12.6)
 9100 format(6(1pe14.6))         !needs to match other formats for ntrplt
 9101 format (6(2x,i5))
      do 910 j=1,nj
      pdum(j) = 0.0
  910 if (wneo(3,3) .ne. 0.0)  pdum(j) = xkineo(j)/wneo(3,3)

      write (ntrplt, 9100) (r(i), i=1, nj)
      write (ntrplt, 9100) (pdum(i), i=1, nj)
      if (wneo(3,3) .ne. 0.0) then
            const = 1.0 / wneo(3,3)
            call copya(chiineo,pdum,nj)
            call multpl1(pdum,nj,const)
      end if

      write (ntrplt, 9100) (pdum(i)  , i=1, nj)  ! chiineo w/o wneo(3,3)
      write (ntrplt, 9100) (xkiinv(i), i=1, nj)
      write (ntrplt, 9100) (chiinv(i), i=1, nj)
      write (ntrplt, 9100) (chiiinv_transp(i), i=1, nj)
      write (ntrplt, 9100) (chieinv_transp(i), i=1, nj)
      write (ntrplt, 9100) (xkeinv(i), i=1, nj)
      write (ntrplt, 9100) (chieinv(i), i=1, nj)
      write (ntrplt, 9100) (xchierl(i), i=1, nj)
      write (ntrplt, 9100) (xchiirl(i), i=1, nj)
      write (ntrplt, 9100) (xeffrl(i), i=1, nj)
      write (ntrplt, 9100) (xeffsf(i), i=1, nj)
      write (ntrplt, 9100) (crit_grad(i), i=1, nj)
cHSJ      write (ntrplt, 9100) (grad_te(i), i=1, nj)
      write (ntrplt, 9100) (rmajorvec(i), i=1, nj)
      write (ntrplt, 9100) (hcap(i), i=1, nj)
      write (ntrplt, 9100) (gcap(i), i=1, nj)
      write (ntrplt, 9100) (fcap(i), i=1, nj)
      write (ntrplt, 9100) (rmajor,  i=1, nj) ! treat as constant vector
c     write (ntrplt, 9100) (r(i),    i=1, nj)
      write (ntrplt, 9100) (chietrinv(i), i=1,nj)
      write (ntrplt, 9100) (chiitrinv(i), i=1,nj)
      write (ntrplt, 9100) (gconde(i), i=1,nj)
      write (ntrplt, 9100) (gconve(i), i=1,nj)
      write (ntrplt, 9100) (gconvde(i), i=1,nj)
      write (ntrplt, 9100) (gconvi(i), i=1,nj)
      write (ntrplt, 9100) (gcondi(i), i=1,nj)
      write (ntrplt, 9100) (gconvdi(i), i=1,nj)
      write (ntrplt, 9100) (eta(j) * 8.98755179e11, j=1,nj) ! ohm-cm
c
c --- toroidal rotation frequency
c
      if (iangrot .eq. 0)  go to 920
      write (ntrplt, 9100)  (angrot(i), i=1,nj)
c
c --- perpendicular velocity shear from Staebler-Hinton model
c
 920  if (ifsflag .eq. 0)  go to 930
      write (ntrplt, 9100)  (sperp(i),i=1,nj-1)
c
c --- Hsieh model data (wshay and scsmult were output in "cray101.f")
c
 930  if (wshay .gt. 0.0) then
        ktimes = MAX (ktimes, 1)
        write (ntrplt, 9100) smulta(ktimes),   smultstder(ktimes),
     .                       skimulta(ktimes), skimultstder(ktimes)
        write (ntrplt, 9100) slim95(ktimes),skilim95(ktimes),
     .                       snexp,sbpexp,srexp
        write (ntrplt, 9100) sbigrexp,stexp,sdtdrexp
        write (ntrplt, 9101) irinshay,iroutshay,irinshaycev,
     .                       iroutshaycev
        write (ntrplt, 9100) (xkeshay(i), i=1,nj)
        write (ntrplt, 9100) (xkishay(i), i=1,nj)
        write (ntrplt, 9100) (  schie(i), i=1,nj)
        write (ntrplt, 9100) (  schii(i), i=1,nj)
      end if

c
c --- Weiland model data (wweiland and include_weiland were output
c                         in file "cray101.f")
c
      if (include_weiland .gt. 0.0) then
        write (ntrplt, 9100) (xke_weiland(i), i=1,nj)
        write (ntrplt, 9100) (xki_weiland(i), i=1,nj)
        write (ntrplt, 9100) (xchie_weiland(i), i=1,nj)
        write (ntrplt, 9100) (xchii_weiland(i), i=1,nj)
        write (ntrplt, 9100) (d_weiland(i), i=1,nj)
        write (ntrplt, 9100) (xkangwl(i), i=1,nj)
        write (ntrplt, 9100) (xeffwl(i), i=1,nj)
      end if
c
c --- IFS model data (dorl_kotch and include_ifs were output
c                         in file "cray101.f")
c

      if (include_ifs .gt. 0.0) then
        write (ntrplt, 9100) (xke_ifs(i), i=1,nj)
        write (ntrplt, 9100) (xki_ifs(i), i=1,nj)
        write (ntrplt, 9100) (chi_e_ifs(i), i=1,nj)
        write (ntrplt, 9100) (chi_i_ifs(i), i=1,nj)
        write (ntrplt, 9100) (d_ifs(i), i=1,nj)
        write (ntrplt, 9100) (xkang_ifs(i), i=1,nj)
        rhod_ifs(1) = rhod_max_ifs    ! use rhod(1) to save this value
        write (ntrplt, 9100) (rhod_ifs(i), i=1,nj)
        rhod_ifs(1) = 0.0             ! restore the original value
        write (ntrplt, 9100) (rlt_ifs(i), i=1,nj)
        write (ntrplt, 9100) (rln_ifs(i), i=1,nj)
        write (ntrplt, 9100) (rlne_ifs(i), i=1,nj)
        write (ntrplt, 9100) (shat_ifs(i), i=1,nj)
        write (ntrplt, 9100) (RLTcrit_ifs(i), i=1,nj)
        write (ntrplt, 9100) (RLTcritz_ifs(i), i=1,nj)
        write (ntrplt, 9100) (g_ifs(i), i=1,nj)
      end if



      write (ntrplt, 9100) (xndinv(i), i=1,nj)
      write (ntrplt, 9100) (qohm  (i), i=1,nj)
      write (ntrplt, 9100) (qdelt (i), i=1,nj)
      write (ntrplt, 9100) (qrad  (i), i=1,nj)
      write (ntrplt, 9100) (qconde(i), i=1,nj)
      write (ntrplt, 9100) (qconve(i), i=1,nj)
      call copya (qsawe,swork,nj)
      if (itimav .eq. 0 .and. w2mix .gt. 0.0)  call zeroa (swork,nj)
      write (ntrplt, 9100) (swork (i), i=1,nj)
      write (ntrplt, 9100) (dpedt (i), i=1,nj)
      write (ntrplt, 9100) (qbeame(i), i=1,nj)
      write (ntrplt, 9100) (qrfe  (i), i=1,nj)
      write (ntrplt, 9100) (qfuse (i), i=1,nj)
      write (ntrplt, 9100) (qneut (i), i=1,nj)
      write (ntrplt, 9100) (qcondi(i), i=1,nj)
      write (ntrplt, 9100) (qconvi(i), i=1,nj)
c
      call copya (qsawi, swork, nj)
      if (itimav .eq. 0 .and. w3mix .gt. 0.0)  call zeroa (swork, nj)
      write (ntrplt, 9100) (swork (i), i=1,nj)
      write (ntrplt, 9100) (dpidt_tot (i), i=1,nj)
      write (ntrplt, 9100) (qbeami(i), i=1,nj)
      write (ntrplt, 9100) (qrfi  (i), i=1,nj)
      write (ntrplt, 9100) (qfusi (i), i=1,nj)
      write (ntrplt, 9100) (pohmic(i), i=1,nj)
      write (ntrplt, 9100) (pdelte(i), i=1,nj)
      write (ntrplt, 9100) (prad  (i), i=1,nj)
      write (ntrplt, 9100) (pconde(i), i=1,nj)
      write (ntrplt, 9100) (pconve(i), i=1,nj)

c
      call copya (psawe, swork, nj)
      if (itimav .eq. 0 .and. w2mix .gt. 0.0)  call zeroa (swork, nj)
      write (ntrplt, 9100) (swork (i), i=1,nj)
      write (ntrplt, 9100) (pdpedt(i), i=1,nj)
      write (ntrplt, 9100) (pbeame(i), i=1,nj)
      write (ntrplt, 9100) (prfe  (i), i=1,nj)
      write (ntrplt, 9100) (pfuse (i), i=1,nj)
      write (ntrplt, 9100) (pdelti(i), i=1,nj)
      write (ntrplt, 9100) (pneut (i), i=1,nj)
      write (ntrplt, 9100) (pcondi(i), i=1,nj)
      write (ntrplt, 9100) (pconvi(i), i=1,nj)
      call copya (psawi,swork,nj)
      if (itimav .eq. 0 .and. w3mix .gt. 0.0)  call zeroa (swork,nj)
      write (ntrplt, 9100) (swork  (i), i=1,nj)
      write (ntrplt, 9100) (pdpidt (i), i=1,nj)
      write (ntrplt, 9100) (pbeami (i), i=1,nj)
      write (ntrplt, 9100) (prfi   (i), i=1,nj)
      write (ntrplt, 9100) (pfusi  (i), i=1,nj)
      write (ntrplt, 9100) (ataue  (i), i=1,nj)
      write (ntrplt, 9100) (tauetr (i), i=1,nj)
      write (ntrplt, 9100) (tauitr (i), i=1,nj)
      write (ntrplt, 9100) (tauelc (i), i=1,nj)
      write (ntrplt, 9100) (tauer  (i), i=1,nj)
      write (ntrplt, 9100) (tauir  (i), i=1,nj)
      write (ntrplt, 9100) (taupe  (i), i=1,nj)
      write (ntrplt, 9100) (te     (i), i=1,nj)
      write (ntrplt, 9100) (ti     (i), i=1,nj)
      write (ntrplt, 9100) (exptlne(i), i=1,nj)
      write (ntrplt, 9100) (exptlte(i), i=1,nj)
      write (ntrplt, 9100) (exptlti(i), i=1,nj)

c
c --- toroidal rotation data
c

      if (iangrot .eq. 1)
     .write (ntrplt, 9100)  (exptlangrot(i), i=1,nj)
      write (ntrplt, 9100)  (exptlcur   (i), i=1,nj)
c
      m1 = 1
      if (nneu .eq. 1)  m1 = in
c
      write (ntrplt, 9100) (tn(i,m1), i=1,nj)
      write (ntrplt, 9100) (  ene(i), i=1,nj)
      write (ntrplt, 9100) ((en(i,j), j=1,nprim)     , i=1,nj)
      write (ntrplt, 9100) ((en(i,j), j=nprim+1,nion), i=1,nj)
      write (ntrplt, 9100) (    enn(i,m1), i=1,nj)
      write (ntrplt, 9100) ( spflux(k,m1), k=1,nengn)
      write (ntrplt, 9100) (englstn(k)   , k=1,nengn)
c
 1090 format (/// 10x,'particle balance and confinement time' //
     .        4x,'j',4x,'r',3x,'r/a',6x,'den/dt',5x,'div.flux',6x,
     .        'source',7x,'rhs-lhs',7x,'taup' /
     .        7x,'(cm)',8x,4('(1/cm**3-s)',2x),4x,'(s)')
 1091 format (///10x,'momentum balance and confinement time' //
     .       4x,'j',4x,'r',3x,'r/a',6x,'dang/dt',5x,'div.flux',6x,
     .       'source',7x,'rhs-lhs',7x,'tauang' /
     .       7x,'(cm)',8x,4('(g/cm*s**2)',2x),4x,'(s)')
 1110 format (20x,'particle sources (1/cm**3-s) for species #',i2,
     .       ', name: ',a2,4x,'itenp =',i2,4x,'ineut =', i2 /
     .       20x,'Note: sother includes  1-1/2d and',
     .       ' sawtooth sources' //
     .       4x,'j',4x,'r',3x,'r/a',6x,'sion',8x,'srecom',
     .       5x,'scx',9x,'sbeam',8x,'sbcx',8x,'sfusion',5x,'sother',
     .       6x,'sum'/7x,'(cm)')
 1115 format (20x,'particle sources (1/cm**3-s) for electrons' //
     .       4x,'j',4x,'r',3x,'r/a',4x,'sion-neut',5x,'srecom',
     .       6x,'sion-imp',7x,'sbeam',9x,'s2d',
     .       9x,'ssawe',9x,'sum'/7x,'(cm)')
 1120 format (20x,'electron energy sources (W/cm**3)','  mode =',a //
     . 4x,'j',4x,'r',3x,'r/a',3x,'1.5*dpe/dt', '  qconde',6x,'qconve',
     . 5x,' qdelt',8x,'qexch',6x,'qohm',7x,' qione',
     .       6x,' qrad',9x,'omegale')
 1121 format (7x,'(cm)')
 1122 format (/20x,'integrated electron energy sources (W)' //
     .        4x,'j',4x,'r',3x,'r/a',1x,'int 1.5*dpe/dt',3x,'pconde',
     .        3x,'pconve',5x,' pdelt',7x,'pexch',7x,'pohm',7x,
     .       ' pione',6x,' prad',8x,'pomegale')
 1123 format (/ 4x,'j',4x,'r',3x,'r/a',6x,'pbeame',7x,'prfe',9x,'pe2d',
     .          6x,'ptfuse',6x,'pbfuse',7x,'psawe',5x,'rhs-lhs',
     .          5x,'pmag')
 1125 format (/4x,'j',4x,'r',3x,'r/a',6x,'qbeame',7x,'qrfe',9x,'qe2d',
     .         6x,'qtfuse',6x,'qbfuse',7x,'qsawe',5x,'rhs-lhs',
     .         6x,'qmag')
 1135 format (/4x,'j',4x,'r',3x,'r/a',3x,'qbeami',6x,'qrfi',
     .         8x,'qi2d',7x,'qtfusi',6x,'qbfusi',7x,'qsawi',
     .         6x,'rhs-lhs',6x,'qbcx')
 1140 format (20x,'ion energy sources (W/cm**3)','  mode =',a //
     .   4x,'j',4x,'r',3x,'r/a',3x,'1.5*dpi/dt', '  qcondi',6x,'qconvi',
     .   6x,'qdelt',8x,'qexch',6x,'qioni',8x,' qcx',6x,'qomegapi',
     .   6x,'omegale')
 1143 format (8x,'charge balance, i.e., Faraday''s law (G/cm-s)',
     .      5x,'twkfar = ',f7.1,'  mode =',a //
     .      4x,'j',4x,'r',3x,'r/a',6x,'dbp/dt',4x,'r*div.flux',
     .      5x,'r*source',8x,'rhs-lhs'/7x,'(cm)')
 1144 format (/20x,'integrated ion energy sources (W)' //
     .    4x,'j',4x,'r',3x,'r/a',1x,'int 1.5*dpi/dt',3x,'pcondi',
     .    3x,'pconvi',6x,'pdelt',8x,'pexch',6x,'pioni',
     .    6x,' pcx',9x,'pomegapi',5x,'pomegale')
 1145 format (/4x,'j',4x,'r',3x,'r/a',5x,'pbeami',8x,'prfi',
     .       9x,'pi2d',6x,'ptfusi',6x,'pbfusi',7x,'psawi',
     .       5x,'rhs-lhs',7x,'pbcx')
 1146 format (20x,'neutral beam injection sources, beam number ',i2,
     .       /30x,'component',i6,
     .       //15x,2( '  normalized'),4x,'average',
     .       5x,'fast ion',6x,'delayed',8x,'energy fraction',
     .       5x,'p. slowing',3x,'e. slowing'
     .       /4x,'j',7x,'r',1x,2(5x,'hot ion'), '  pitch angle',
     .       2(3x,'e. source*'),8x,'deposited in',
     .       4x,2(4x,'down time')
     .       /10x,'(cm)',3x,'birth rate',3x,'dep. rate',5x,'cosine',1x,
     .       2(4x,'(W/cm**3)'),4x,'electrons',6x,'ions',
     .       2(8x,'(s)',2x))
 1152 format (/7x,'average',1x,1p9e12.3)
 1153 format (1x,i4,f9.2,3f12.4,1x,1p2e13.3,0p2f12.4,1p2e13.3)
*1154 format (// 20x,'components of RF heating power (W)' // 19x,
**** .      'powabs',5x,'powcol',9x,' powe',8x,'powd' / 14x,1p9e13.3)
 1160 format (20x,'selected transport coefficients'
     .        //4x,'j',9x,'r',7x,'d(',i1,')',
     .        8x,'chie',8x,'chii',9x,'xke',9x,'xki',
     .        8x,'xnuse',6x,a2,'-xnus',6x,'ftrap',5x,'eta' /
     .        12x, '(cm)', 3(3x,'(cm**2/s)'),2(4x,'(1/cm-s)'),40x,
     .        '(ohm-cm)' /
     .        20x,a,10x,a,6x,a,6x,a,6x,a)
 1161 format (4x,'j',9x,'r',7x,'ydebug(1-8)')
 1162 format (4x,'j',9x,'r',7x,'ydebug(8-16)')
 1170 format (/20x,'confinement times (s)' //
     .       4x,'j',4x,'r',3x,'r/a',6x,'taue',7x,'tauetr',5x,'tauitr',
     .       6x,'tauer',7x,'tauir',7x,'taupe',7x,'tauelc',6x,'tres',
     .       6x,'tangmtmloc' /
     .       7x,'(cm)')
 1201 format (/20x,'diagnostic transport coefficients' //
     .          4x,'j',9x,'r',7x,'xkecar',7x,'xkeohk',
     .          7x,'xkepp ',7x,'xkedom',
     .     7x,' xki/ ', '  shearp',7x,'slene',8x,'slte',7x,'chiwneo' /
     .    12x,'(cm)',5x,4('(1/cm-s)',5x),' xkineo',15x,'(cm)',9x,'(cm)',
     .          5x,'cm**2/sec')

c the following  format block also appers in sub P_Nfreya_output
 8110 format (// ' particle energy (keV)', 22x, 'total',6(4x,f8.3))
 8120 format (' ion beam intensity (part./s)',24x,1p6e12.3)
 8130 format (' neutral beam intensity (part./s) to ap.',13x,1p6e12.3)
 8140 format (' neutral beam power (W) to aperture',6x,1p7e12.3)
 8150 format (' fraction stopped by aperture',9x,7(3x,f9.4))
 8160 format (' fraction incident on wall (shinethrough)',
     .          f9.4,6(3x,f9.4))
 8161 format (' fraction lost on orbits',14x,7(3x,f9.4))
 8162 format (' neutral beam power (W) in plasma ',7x,1p7e12.3)
 8163 format (' slowed  beam power (W) in plasma*',7x,1p7e12.3)
 8164 format (' fraction deposited in electrons',6x,7(3x,f9.4))
 8166 format (' fraction deposited in ions',11x,7(3x,f9.4))
 8168 format (' fraction lost by fast ion charge ex.',1x,
     .       7(3x,f9.4))
 8170 format (/' * excludes aperture, shinethrough, and orbit',
     .       ' losses, but includes loss due to fast ion charge ex.')
 8180 format (/ ' fraction and width of various type orbits' //
     .      15x, 'type',16x,'fraction',7x,'width'/51x,'(cm)' /
     .       3x, 'passing and axis-circling',3x,f12.4,f12.2  /
     .       3x, 'passing and not circling',4x,f12.4,f12.2   /
     .       3x, 'trapped and axis-circling',3x,f12.4,f12.2  /
     .       3x, 'trapped and not circling',4x,f12.4,f12.2   /
     .       3x, 'lost on orbit',15x,f12.4                   /
     .       3x, 'error detected',14x,f12.4)
 8200 format (   20x, 'neutral beam injection sources')
 8210 format (// 20x, 'fast ion power source and fraction deposited',
     .               ' in electrons and ions    (beam',i3,')' //
     .                 3x,'energy (keV):',12x,3(f8.3,24x) /
     .                 4x,'j',4x,'r',3x,'r/a',5x,
     .                 3(1x,'W/cm**3',5x,'elec',4x,'ions',7x))
 8220 format (1x,i4,f6.1,f5.1,3(2x,1pe12.3,0p2f8.3,2x))
*%%%%      write (6, *) '%%%% exiting  subroutine PSOURC'


      return
c
      end


