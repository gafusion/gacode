      subroutine iter_dbase129
c   NOTE: this is version 1.29 of cray204.f it was added as
c         a separate iterdb routine to statisfy old users of the
c         routine. It is invoked by using the switch iterdb = 1
c         The routine in cray204.f is now invoked by using iterdb = 2.
c
c
c
c ----------------------------------------------------------------------
c --- collects and writes out (or reads in) data for ITER database
c ----------------------------------------------------------------------
c
c --- INPUT
c
c INCLUDE file constnts.i
c     psimks           conversion factor convers psi from kg/cm**2
c                      to volt*second/radian
c
c INCLUDE file iterdb.i:
c     irwflag          read/write flag, included to make the reading
c                      of the file created by this subroutine
c                      consistent with the writing of the file.
c                      irwflag = 0  ==>  WRITE the data
c                      irwflag = 1  ==>  READ  the data
c                      if this subroutine is used in read mode note that
c                      the relevant arrays that are read in are defined
c                      in INCLUDE files as described below.
c                      NOTE IF IRWFLAG=1, THE I/O UNIT NUMBER ASSIGNED
c                      TO THE FILE FOR READING PURPOSES, NITERDB, IS CURRENTLY
c                      HARD CODED AT 11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     iterdb           if set to 1 then this subroutine will execute.
c                      Otherwise this subroutine returns immediately without
c                      performing any calculations. (by doing it this way
c                      the necessary INCLUDE files can all be brought in
c                      at this subroutine, rather than at the calling level)
c
c     iterdbfilename   name of file which will be created or appended to
c     iterdsc          switch used to output description of parameters
c     niterdb          Fortran I/O unit number for output OR INPUT - SEE ABOVE
c
c INCLUDE file contour.i:
c     nplasbdry
c     rplasbdry(j)
c     zplasbdry(j)     j=1,2..nplasbdry describes the plasma boundary
c
c INCLUDE file extra.i:
c     q(j)             j=1,2,...nj, the safety factor on rho grid
c     betap
c     beta
c
c INCLUDE file etc.i:
c     tocur            total plasma current, amps
c
c INCLUDE file io.i
c     ncrt
c     nout             unit numbers for diagnostic output
c INCLUDE file psig.i:
c     rho(j)           j=1,2,..npsi the rho values corresponding
c                      to the psi grid in meters
c     qpsi(j)          the corresponding q values
c     psivolp(j)       volume, meter**3, enclosed by flux surface
c     grho1(j)         < ABS (grad rho)>
c     grho2(j)         <(grad rho)**2>
c     rmajavnpsi(j)    j=1,2,..npsi average major radius
c                      at elevation of magnetic axis
c     rminavnpsi(j)    same for minor radius
c     psival(j)        j=1,2..npsi
c
c     triangnpsi(j)    triangularity
c     pindentnpsi(j)   indentation
c     sfareanpsi(j)    surface area of flux surfaces
c     cxareanpsi(j)    cross-sectional area
c
c     btor             mag field at R0, in vacuum, tesla
c     elong(j)         elongation
c
c INCLUDE file geom.i:
c     fcap(j)
c     gcap(j)
c     hcap(j)          the geometric factors
c     ali              plasma inductance
c
c INCLUDE file ions.i:
c     namep(i)         i=1,..nprim
c     namei(i)         i=1,..nimp
c     namen(i)         i=1,..nneu
c     zeff(j)          j=1,2...nj
c     nameb
c
c INCLUDE file nub.i:
c     ibion            index (into namep) of beam species
c     sbion
c
c INCLUDE file nub2.i:
c     enbeam(j)        j=1,2..nj  fast ion density due to beam
c
c INCLUDE file neut.i:
c     enn(j,k)         j=1,2..nj(njs),k=1,2 neutral density
c     ennw             neutral density due to wall source
c     ennv                                    volume
c     volsn            volume source of neutrals
c
c INCLUDE file machin.i:
c     rmajor
c     kappa
c     btor
c
c INCLUDE file mesh.i
c     r(j)             j=1,2...nj,rho grid
c
c INCLUDE file solcon.i
c     time            current time, sec
c
c INCLUDE file numbrs.i
c     nj              transport grid size
c     nprim
c     nneu
c     nimp
c     nion
c
c INCLUDE file rhog.i
c     psir(j)         j=1,2...nj psi on r (i.e., rho) grid
c
c INCLUDE file soln.i
c     te(j)           electron and ion temperatures,keV
c     ti(j)
c     ene(j)          electron density
c     en(j,k)         j=1,2..nj,k=1..nion,where nion=nprim+nimp is the
c                     number of primary plus impurity ion species
c     curden(j)       current  density
c     rbp(j)          rho*bp field
c
c INCLUDE file sourc.i:
c     curboot(j)
c     curohm(j)
c     currf(j)
c     curbe(i)
c     curbi(j)
c     qrad(j)         radiated power
c     totohm
c     totbeam
c     totboot
c     totrf           total currents, amps
c     dpedtc(j), dpidtc(j)
c     qconde,qcondi
c     qconve,qconvi
c     qdelt
c     qexch
c     qione,qioni
c     qbeame,qbeami
c     qrfe, qrfi
c     qe2d, qi2d
c     qtfuse,qtfusi
c     qmag
c     qsawe,qsawi
c     qbfusi
c     sbeam
c
c INCLUDE file tcoef.i:
c     cheinv(j)
c     chiinv(j)        electron and ion thermal cond.
c     xkineo(j)        ion neoclassical thermal conductivity
c
c INCLUDE file tordlrot.i:
c     angrot(j)       angular rotation speed, radian/second
c     storqueb(j)     beam torque density nt-m/m**3 !!error the value is in dyne cm /cm**3 HSJ 8/3/01  !! corrected 03/16/04 HSJ
c                   
c
c INCLUDE file  yoka.i
c     ishot           shot number
c
c --- OUTPUT
c
c     ALL OUTPUT IS TO FILE POINTED TO BY ITERDBFILENAME
c     THIS ROUTINE DOES NOT (AND SHOULD NOT IN THE FUTURE)
c     PERFORM ANY CALCULATIONS WHICH DEFINE NEW PARAMETERS
c     BECAUSE THIS IS ONLY AN OUTPUT ROUTINE!
c
c ----------------------------------------------------- 11/22/94 --- HSJ
c


c
      USE param
      USE ions
      USE neut
      USE solcon
      USE soln
      USE io
      USE contour
      USE mhdpar
      USE nub  
      USE nub2
      USE extra
      USE yoka
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE geom
      USE tordlrot
      USE constnts
      USE tcoef
      
      USE psig
      USE rhog
      USE flxav
      USE iterdbmd
      USE etc
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: iterdb_v1.29.f,v 1.5 2011/12/09 20:07:57 stjohn Exp $"/

      include 'storage.i'

c
      dimension  aspline(kpsi), bspline(kpsi), cspline(kpsi),
     .           cs2spline(kpsi,3), dspline(kpsi), espline(kpsi),
     .           fspline(kpsi), bpar(4), work(kj)
      logical    monotonic
      character  starflag*2, headerline*132
c
      integer    LENGTH
c
c      equivalence (aspline(1)    , xdum(1))
c      equivalence (bspline(1)    , xdum(kpsi+1))
c      equivalence (cspline(1)    , wdum(1))
c      equivalence (dspline(1)    , vdum(1))
c      equivalence (espline(1)    , ydum(1))
c      equivalence (fspline(1)    , ydum(kpsi+1))
c      equivalence (work(1)       , zdum(1))
c      equivalence (cs2spline(1,1), sdum(1))
c      equivalence (bpar(1)       , rdum(1))
c
      data icall/0/, convert/1.60217733e-10/ ! keV/sec/cm*3 > watts/m**3
      save icall   , convert
c
      if (iterdb .ne. 1)  return
      iterdsc = 1                   ! always write descriptor
c
  7   format (         a   )
  8   format (5(2x,    a  ))        ! common character write/read format
  9   format (5(2x,   i6  ))        ! common integer   write/read format
 10   format (5(2x,1pe14.4))        ! common floating  write/read format
c
c OPENing method depends on whether read or write, and whether or not first time
c
      call getioun(niterdb,niterdb) ! need a free ioc every call

      if (irwflag .eq. 0) then      ! WRITE to file, might be existing
        if (icall .eq. 0) then      !   first time, so delete, then open
          call DESTROY (iterdbfilename)
          open (unit = niterdb, file = iterdbfilename, status = 'NEW')
        else                        !   open to append to existing file
          call OPEN_APPEND (niterdb, iterdbfilename)
        end if
      else                          ! READ from existing file
        open (unit = niterdb, file = iterdbfilename, status = 'OLD',
     .        err = 2)
        go to 3
    2   write  (ncrt, 4)  TRIM(iterdbfilename)
        write  (nout, 4)  TRIM(iterdbfilename)
    4   format (                                                     /
     .      ' ERROR: subroutine ITER_DBASE has encountered an error' /
     .      '        the ITER database file "', a, '" cannot be opened')
        call giveupus(niterdb)
        call STOP ('subroutine ITER_DBASE: cannot open database', 171)
    3   continue
      end if
c
c write a header line each time routine is called:
c   header lines are identified by     **                  in first two columns
c   other comment lines (if any) have  *b  (where b=blank) in first two columns
c



      if (irwflag .eq. 0) then
          write  (niterdb, 1)  ishot, time
    1     format ('** shot number = ', i6, '  time = ', 1pe14.6)

c
c --- check that psir is monotonic; it may not be in certain cases where
c --- the current profile was evolved
c
        call check_monotonic (psir, nj, monotonic, 1)
        if (.not. monotonic) then    ! psir is not monotonic
          write  (niterdb, 5) time
    5     format (' the psi grid, calculated from the poloidal' /
     .            ' B field evolution, is not monotonic'        /
     .            ' therefore the data at this time (', 1pe12.6,
     .            ') was not calculated')
          go to 2000
        end if
      else
        read (niterdb, 7) headerline
      end if
c
c --- some scalar quantities:
c --- integer and character parameters
c
      if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 250)
 250      format ('*  ishot : shot number')
          write  (niterdb, 9) ishot
          if (iterdsc .ne. 0)
     .    write  (niterdb, 300)
 300      format ('*  nj : the size of the vectors printed',
     .                   ' in this file')
          write  (niterdb, 9) nj
          if (iterdsc .ne. 0)
     .    write  (niterdb, 310)
 310      format ('*  nion : the number of ion species')
          write  (niterdb, 9) nion
          if (iterdsc .ne. 0)
     .    write  (niterdb, 320)
 320      format ('*  nprim : the number of primary ion species')
          write  (niterdb, 9) nprim
          if (iterdsc .ne. 0)
     .    write  (niterdb, 330)
 330      format ('*  nimp : the number of impurity ion species')
          write  (niterdb, 9) nimp
          if (iterdsc .ne. 0)
     .    write  (niterdb, 340)
 340      format ('*  nneu : the number of neutral ion species')
          write  (niterdb, 9) nneu
          if (iterdsc .ne. 0)
     .    write  (niterdb, 350)
 350      format ('*  ibion : index of beam species,-1 means',
     .            ' beam is dt mixture')
          write  (niterdb, 9) ibion
          if (iterdsc .ne. 0)
     .    write  (niterdb, 360)
 360      format ('*  namep : name(s) of primary ion species')
          write  (niterdb, 8) (namep(i),i=1,nprim)
          if (iterdsc .ne. 0)
     .    write  (niterdb, 370)
 370      format ('*  namei : name(s) of impurity ion species')
          write  (niterdb, 8) (namei(i),i=1,nimp)
          if (iterdsc .ne. 0)
     .    write  (niterdb, 380)
 380      format ('*  namen : name(s) of neutral ion species')
          write  (niterdb, 8) (namen(i), i=1,nneu)
       else
          read   (niterdb, 7) starflag
          read   (niterdb, 9) ishot
          read   (niterdb, 7) starflag
          read   (niterdb, 9) nj
          read   (niterdb, 7) starflag
          read   (niterdb, 9) nion
          read   (niterdb, 7) starflag
          read   (niterdb, 9) nprim
          read   (niterdb, 7) starflag
          read   (niterdb, 9) nimp
          read   (niterdb, 7) starflag
          read   (niterdb, 9) nneu
          read   (niterdb, 7) starflag
          read   (niterdb, 9) ibion
          read   (niterdb, 7) starflag
          read   (niterdb, 8) (namep(i),i=1,nprim)
          read   (niterdb, 7) starflag
          read   (niterdb, 8) (namei(i),i=1,nimp)
          read   (niterdb, 7) starflag
          read   (niterdb, 8) (namen(i),i=1,nneu)
       end if
c
c --- real (i.e., floating point) parameters
c
      if (irwflag .eq. 0) then
c
      rgeom    = rmajavnpsi(1)
      rmag     = rmajavnpsi(npsi)
      deltao   = triangnpsi(1)
      areao    = cxareanpsi(1)
      pindento = pindentnpsi(1)
      volo     = psivolp(1)
      te0      = te(1)
      ti0      = ti(1)
c
c      if (irwflag .eq. 0) then
         if (iterdsc .ne. 0)
     .   write  (niterdb, 490)
 490     format ('*  time : time at which data is printed')
         write  (niterdb, 10) time
         if (iterdsc .ne. 0)
     .   write  (niterdb, 500)
 500     format ('*  Rgeom : major radius of geometric',
     .             ' center at elevation of magnetic axis, meters')
         write  (niterdb, 10) rgeom*0.01
         if (iterdsc .ne. 0)
     .   write  (niterdb, 501)
 501     format ('*  Rmag : major radius of magnetic axis, meters')
         write  (niterdb, 10) rmag*0.01
         if (iterdsc .ne. 0)
     .   write  (niterdb, 502)
 502     format ('*  R0 : major radius of vacuum btor ref',
     .             ' location, meters')
         write  (niterdb, 10) rmajor*0.01
         if (iterdsc .ne. 0)
     .   write  (niterdb, 503)
 503     format ('*  kappa : plasma elongation')
         write  (niterdb, 10) kappa
         if (iterdsc .ne. 0)
     .   write  (niterdb, 504)
 504     format ('*  delta : plasma triangularity')
         write  (niterdb, 10) deltao
         if (iterdsc .ne. 0)
     .   write  (niterdb, 505)
 505     format ('*  pindent : plasma indentation')
         write  (niterdb, 10) pindento
         if (iterdsc .ne. 0)
     .   write  (niterdb, 506)
 506     format ('*  volo : plasma volume, meters**3')
         write  (niterdb, 10) volo*1.0e-6
         if (iterdsc .ne. 0)
     .   write  (niterdb, 507)
 507     format ('*  cxareao : plasma cross-sectional area, meters**2')
         write  (niterdb, 10) areao*1.0e-4
         if (iterdsc .ne. 0)
     .   write  (niterdb, 508)
 508     format ('*  Btor : vacuum toroidal field at rmajor, tesla')
         write  (niterdb, 10) btor*1.0e-4
         if (iterdsc .ne. 0)
     .   write  (niterdb, 509)
 509     format ('*  total, ohmic, bootstrap, beam and RF',
     .             ' currents, amps')
         write  (niterdb, 10) tocur, totohm, totboot, totbeam, totrf
         if (iterdsc .ne. 0)
     .   write  (niterdb, 510)
 510     format ('*  betap : poloidal beta')
         write  (niterdb, 10) betap
         if (iterdsc .ne. 0)
     .   write  (niterdb, 511)
 511     format ('*  beta : toroidal beta')
         write  (niterdb, 10) beta
         if (iterdsc .ne. 0)
     .   write  (niterdb, 512)
 512     format ('*  ali : plasma inductance')
         write  (niterdb, 10) ali
         if (iterdsc .ne. 0)
     .   write  (niterdb, 513)
 513     format ('*  te0 : central electron temperature')
         write  (niterdb, 10) te0
         if (iterdsc .ne. 0)
     .   write  (niterdb, 514)
 514     format ('*  ti0 : central ion temperature')
         write  (niterdb, 10) ti0
      else
         read   (niterdb,  7) starflag
         read   (niterdb, 10) time
         read   (niterdb,  7) starflag
         read   (niterdb, 10) rgeom
         read   (niterdb,  7) starflag
         read   (niterdb, 10) rmag
         read   (niterdb,  7) starflag
         read   (niterdb, 10) rmajor
         read   (niterdb,  7) starflag
         read   (niterdb, 10) kappa
         read   (niterdb,  7) starflag
         read   (niterdb, 10) deltao
         read   (niterdb,  7) starflag
         read   (niterdb, 10) pindento
         read   (niterdb,  7) starflag
         read   (niterdb, 10) volo
         read   (niterdb,  7) starflag
         read   (niterdb, 10) areao
         read   (niterdb,  7) starflag
         read   (niterdb, 10) btor
         read   (niterdb,  7) starflag
         read   (niterdb, 10) tocur, totohm, totboot, totbeam, totrf
         read   (niterdb,  7) starflag
         read   (niterdb, 10) betap
         read   (niterdb,  7) starflag
         read   (niterdb, 10) beta
         read   (niterdb,  7) starflag
         read   (niterdb, 10) ali
         read   (niterdb,  7) starflag
         read   (niterdb, 10) te0
         read   (niterdb,  7) starflag
         read   (niterdb, 10) ti0
      end if
c
c --- psir grid
c
      if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 1010)
 1010     format ('*  psi on rho grid, volt*second/radian')
          do j=1,nj
            work(j) = psimks*psir(j)
          end do
          write  (niterdb, 10) (work(j), j=1,nj)
      else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (psir(j), j=1,nj)
      end if
c
c --- rho grid
c
      if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 1020)
 1020     format ('*  rho grid, meters')
          do j=1,nj
            work(j) = r(j) * 1.0e-02
          end do
          write (niterdb, 10) (work(j), j=1,nj)
      else
          read  (niterdb,  7) starflag
          read  (niterdb, 10) (r(j)   , j=1,nj)
      end if
c
c --- fcap
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 1051)
 1051     format ('*  fcap, (i.e., f(psilim)/f(psi))')
          write  (niterdb, 10) (fcap(j), j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (fcap(j), j=1,nj)
       end if
c
c --- gcap
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 1052)
 1052     format ('*  gcap, (i.e., <(grad rho)**2*(R0/R)**2>)')
          write  (niterdb, 10) (gcap(j), j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (gcap(j), j=1,nj)
       end if
c
c --- hcap
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 1053)
 1053     format ('*  hcap, (i.e., (dvolume/drho)/(4*pi*pi*R0*rho))')
          write  (niterdb, 10) (hcap(j), j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (hcap(j), j=1,nj)
       end if
c
c --- te (in keV)
c
      if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 1030)
 1030     format ('*  electron temperature, keV')
          write  (niterdb, 10) (te(j), j=1,nj)
      else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (te(j), j=1,nj)
      end if
c
c --- ti (in keV)
c
      if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 1040)
 1040     format ('*  ion temperatue, keV')
          write  (niterdb, 10) (ti(j), j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (ti(j), j=1,nj)
       end if
c
c --- safety factor
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 1050)
 1050     format ('*  q (i.e., safety factor) profile')
          write  (niterdb, 10) (ABS (q(j)), j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (     q(j) , j=1,nj)
       end if
c
c --- electron density
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3054)
 3054     format ('*  electron density, #/meter**3')
          write  (niterdb, 10) (ene(j)*1.0e+6, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (ene(j)       , j=1,nj)
       end if
c
c --- primary, impurity densities
c
       jp = 0
       ji = 0
       do jj=1,nion
           if (jj .le. nprim)  jp = jp + 1
           if (jj .gt. nprim)  ji = ji + 1
           if (irwflag .eq. 0) then
               if (iterdsc .ne. 0 .and. jj .le. nprim)
     .         write  (niterdb, 3055) namep(jp)
 3055          format ('*  primary ion density,',
     .                 ' #/meter**3, species: ', a)
               if (iterdsc .ne. 0 .and. jj .gt. nprim)
     .         write  (niterdb, 3056) namei(ji)
 3056          format ('*  impurity ion density,',
     .                 ' #/meter**3, species: ', a)
               write  (niterdb, 10) (en(j,jj)*1.0e+6, j=1,nj)
           else
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (en(j,jj)       , j=1,nj)
           end if
        end do
c
c --- primary ion species sources
c
       do jj=1,nprim
           if (irwflag .eq. 0) then
               if (iterdsc .ne. 0)
     .         write  (niterdb, 4155) namep(jj)
 4155          format ('*  sion : source due to ionization,',
     .                 ' #/(meter**3*second), species: ', a)
               write  (niterdb, 10) (sion(j,jj)*1.0e+6, j=1,nj)
c
               if (iterdsc .ne. 0)
     .         write (niterdb, 4156) namep(jj)
 4156          format ('*  srecom : source due to recombination,',
     .                 ' #/(meter**3*second), species: ', a)
               write  (niterdb, 10) (srecom(j,jj)*1.0e+6, j=1,nj)
c
               if (iterdsc .ne. 0)
     .         write  (niterdb, 4157) namep(jj)
 4157          format ('*  scx : source due to cx thermal neut.,',
     .                 ' #/(meter**3*second), species: ', a)
               write  (niterdb, 10) (scx(j,jj)*1.0e+6, j=1,nj)
c
               if (iterdsc .ne. 0)
     .         write  (niterdb, 4158) namep(jj)
 4158          format ('*  sbcx : sink due to cx with beam neut.,',
     .                 ' #/(meter**3*second), species: ', a)
               write  (niterdb, 10) (sbcx(j,jj)*1.0e+6, j=1,nj)
c
               if (iterdsc .ne. 0)
     .         write  (niterdb, 4159) namep(jj)
 4159          format ('*  s : total source rate,',
     .                 ' #/(meter**3*second), species: ', a)
               write  (niterdb, 10) (s(jj,j)*1.0e+6, j=1,nj)
c
               if (iterdsc .ne. 0)
     .         write  (niterdb, 4160) namep(jj)
 4160          format ('*  dudt : s dot,',
     .                 ' #/(meter**3*second), species: ', a)
               write  (niterdb, 10) (dudtsv(jj,j)*1.0e+6, j=1,nj)
c
           else
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (sion(j,jj)  , j=1,nj)
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (srecom(j,jj), j=1,nj)
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (scx(j,jj)   , j=1,nj)
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (sbcx(j,jj)  , j=1,nj)
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (s(jj,j)     , j=1,nj)
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (dudtsv(jj,j), j=1,nj)
           end if
        end do
c
c --- fast ion density
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 4000) nameb
 4000     format ('*  fast ion density, #/meter**3, species: ', a)
          write  (niterdb, 10) (enbeam(j)*1.0e+6, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (enbeam(j)       , j=1,nj)
       end if
c
c --- neutral densities
c
       do jn=1,nneu
           if (irwflag .eq. 0) then
               write  (niterdb, 3057) namen(jn)
 3057          format ('*  neutral density, #/meter**3, species: ', a)
               write  (niterdb, 10) (enn(j,jn)*1.0e+6, j=1,nj)
           else
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (enn(j,jn)       , j=1,nj)
           end if
        end do
c
c --- neutral densities
c
       do jn=1,nneu
           if (irwflag .eq. 0) then
               write  (niterdb, 3200) namen(jn)
 3200          format ('*  neutral density from wall source,',
     .                 ' #/meter**3, species: ', a)
               write  (niterdb, 10) (ennw(j,jn)*1.0e+6, j=1,nj)
           else
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (ennw(j,jn)       , j=1,nj)
           end if
        end do
c
c --- neutral densities
c
       do jn=1,nneu
           if (irwflag .eq. 0) then
               write  (niterdb, 3210) namen(jn)
 3210          format ('*  neutral density from volume source,',
     .                 ' #/meter**3, species: ', a)
               write  (niterdb, 10) (ennv(j,jn)*1.0e+6, j=1,nj)
           else
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (ennv(j,jn)       , j=1,nj)
           end if
        end do
c
c --- neutral source
c
       do jn=1,nneu
           if (irwflag .eq. 0) then
               write  (niterdb, 3420) namen(jn)
 3420          format ('*  volume source of neutrals,',
     .                 ' #/(meter**3*second), species: ', a)
               write  (niterdb, 10) (volsn(j,jn)*1.0e+6, j=1,nj)
           else
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (volsn(j,jn)       , j=1,nj)
           end if
        end do
c
c --- electron source due to beams
c
           if (irwflag .eq. 0) then
               write  (niterdb, 3230)
 3230          format ('*  sbion : beam electron source,',
     .                 ' #/(meter**3*second)')
               write  (niterdb, 10) (sbion(j)*1.0e+6, j=1,nj)
           else
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (sbion(j)       , j=1,nj)
           end if
c
c --- thermal ion source due to beams
c
           if (irwflag .eq. 0) then
               write  (niterdb, 3086)
 3086          format ('*  sbion : beam thermal ion source,',
     .                 ' #/(meter**3*second)')
               write  (niterdb, 10) (sbeam(j)*1.0e+6, j=1,nj)
           else
               read   (niterdb,  7) starflag
               read   (niterdb, 10) (sbeam(j)       , j=1,nj)
           end if
c
c --- current density
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3058)
 3058     format ('*  total current density, amps/meter**2')
          write  (niterdb, 10) (curden(j)*1.0e4, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (curden(j)      , j=1,nj)
       end if
c
c --- ohmic current density
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3059)
 3059   format ('*  ohmic current density, amps/meter**2')
        write  (niterdb, 10) (curohm(j)*1.0e4, j=1,nj)
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (curohm(j)      , j=1,nj)
      end if
c
c --- bootstrap current density
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3060)
 3060   format ('*  bootstrap current density, amps/meter**2')
        write  (niterdb, 10) (curboot(j)*1.0e4, j=1,nj)
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (curboot(j)      , j=1,nj)
      end if
c
c --- beam current density
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3061)
 3061   format ('*  beam-driven current density, amps/meter**2')
        write  (niterdb, 10) (curdbeam(j)*1.0e4, j=1,nj)
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (curdbeam(j)      , j=1,nj)
      end if
c
c --- RF current density
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3070)
 3070   format ('*  RF current density, amps/meter**2')
        write  (niterdb, 10) (currf(j)*1.0e4, j=1,nj)
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (currf(j)      , j=1,nj)
      end if
c
c --- rho*bp0*fcap*gcap*hcap
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3080)
 3080   format ('*  rho*bp0*fcap*gcap*hcap, tesla*meters')
        write  (niterdb, 10) (rbp(j)*1.0e-6, j=1,nj)
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (rbp(j)       , j=1,nj)
      end if
c
c --- zeff profile
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3110)
 3110   format ('*  zeff profile')
        write  (niterdb, 10) (zeff(j), j=1,nj)
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (zeff(j), j=1,nj)
      end if
c
c --- angular rotation speed profile
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3120)
 3120   format ('*  angular rotation speed profile, rad/sec')
        write  (niterdb, 10) (angrot(j), j=1,nj)
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (angrot(j), j=1,nj)
      end if
c
c --- thermal diff. profiles, electron and ion
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3130)
 3130   format ('*  electron thermal diffusivity, meters**2/sec',
     .            ' on half grid')
        write  (niterdb, 10) (chieinv(j)*1.0e-4, j=1,nj)
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (chieinv(j)       , j=1,nj)
      end if
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3140)
 3140   format ('*  ion thermal diffusivity, meters**2/second',
     .            ' on half grid')
        write  (niterdb, 10) (chiinv(j)*1.0e-4, j=1,nj)
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (chiinv(j)       , j=1,nj)
      end if
c
c --- ion neoclassical thermal conductivity
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0) write (niterdb, 3145)
 3145     format ('*  ion neoclassical thermal conductivity,',
     .            ' 1/(meter*second), on half grid')
          write  (niterdb, 10) (xkineo(j) * 100.0, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (xkineo(j)      , j=1,nj)
       end if
c
c --- d(electron energy)/dt profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3150)
 3150     format ('*  wdot, electrons, watts/meter**3')
          write  (niterdb, 10) (dpedtc(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (dpedtc(j)        , j=1,nj)
       end if
c
c --- d(ion energy)/dt profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3160)
 3160     format ('*  wdot, ions, watts/meter**3')
          write  (niterdb, 10) (dpidtc(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (dpidtc(j)        , j=1,nj)
       end if
c
c --- electron conduction profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3170)
 3170     format ('*  electron conduction, watts/meter**3')
          write  (niterdb, 10) (qconde(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qconde(j)        , j=1,nj)
       end if
c
c --- ion conduction profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3180)
 3180     format ('*  ion conduction, watts/meter**3')
          write  (niterdb, 10) (qcondi(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qcondi(j)        , j=1,nj)
       end if
c
c --- electron convection profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3190)
 3190     format ('*  electron convection, watts/meter**3')
          write  (niterdb, 10) (qconve(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qconve(j)        , j=1,nj)
       end if
c
c --- ion convection profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3410)
 3410     format ('*  ion convection, watts/meter**3')
          write  (niterdb, 10) (qconvi(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qconvi(j)        , j=1,nj)
       end if
c
c --- beam electron profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3215)
 3215     format ('*  power to elec. from beam, watts/meter**3')
          write  (niterdb, 10) (qbeame(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qbeame(j)        , j=1,nj)
       end if
c
c --- electron ion equilibration profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3220)
 3220     format ('*  qdelt : electron-ion equilibration,',
     .            ' watts/meter**3')
          write  (niterdb, 10) (-qdelt(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) ( qdelt(j)        , j=1,nj)
       end if
c
c --- beam ion profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3450)
 3450     format ('*  power to ions from beam, watts/meter**3')
          write  (niterdb, 10) (qbeami(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qbeami(j)        , j=1,nj)
       end if
c
c --- RF electron heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3240)
 3240     format ('*  qrfe, RF electron heating, watts/meter**3')
          write  (niterdb, 10) (qrfe(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qrfe(j)        , j=1,nj)
       end if
c
c --- RF ion heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3250)
 3250     format ('*  qrfi, -RF ion heating, watts/meter**3')
          write  (niterdb, 10) (qrfi(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qrfi(j)        , j=1,nj)
       end if
c
c --- qione heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3260)
 3260     format ('*  qione, electron power density due to',
     .              ' recombination and impact ionization,',
     .              ' watts/meter**3')
          write  (niterdb, 10) (-qione(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) ( qione(j)        , j=1,nj)
       end if
c
c --- qioni heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3270)
 3270     format ('*  qioni, ion power density due to',
     .              ' recombination and impact ionization,',
     .              ' watts/meter**3')
          write  (niterdb, 10) (qioni(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qioni(j)        , j=1,nj)
       end if
c
c --- qxc, ion heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3275)
 3275     format ('*  qcx, ion power density due to',
     .              ' neutral-ion charge exchange, watts/meter**3')
          write  (niterdb, 10) (-qcx(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) ( qcx(j)        , j=1,nj)
       end if
c
c --- 2d electron heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3280)
 3280     format ('*  2d electron heating, watts/meter**3')
          write  (niterdb, 10) (qe2d(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qe2d(j)        , j=1,nj)
       end if
c
c --- 2d ion heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3290)
 3290     format ('*  2d ion heating, watts/meter**3')
          write  (niterdb, 10) (qi2d(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qi2d(j)        , j=1,nj)
       end if
c
c --- fusion electron heating  profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3300)
 3300     format ('*  fusion electron heating, watts/meter**3')
          write  (niterdb, 10) (qfuse(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qfuse(j)        , j=1,nj)
       end if
c
c --- fusion ion heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3310)
 3310     format ('*  fusion ion heating, watts/meter**3')
          write  (niterdb, 10) (qfusi(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qfusi(j)        , j=1,nj)
       end if
c
c --- beam fusion electron heating profile
c --- (fraction of beam fusion energy deposited on thermal electron distribution
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3312)
 3312     format ('*  beam fusion electron heating, watts/meter**3')
          write  (niterdb, 10) (qbfuse(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qbfuse(j)        , j=1,nj)
       end if
c
c --- beam fusion ion heating profile
c --- (fraction of beam fusion energy deposited on thermal ion distribution)
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3315)
 3315     format ('*  beam fusion ion heating, watts/meter**3')
          write  (niterdb, 10) (qbfusi(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qbfusi(j)        , j=1,nj)
       end if
c
c --- mag electron heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3320)
 3320     format ('*  qmag electron heating, watts/meter**3')
          write  (niterdb, 10) (qmag(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qmag(j)        , j=1,nj)
       end if
c
c --- sawtooth electron heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3330)
 3330     format ('*  sawtooth electron heating, watts/meter**3')
          write  (niterdb, 10) (qsawe(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qsawe(j)        , j=1,nj)
       end if
c
c --- sawtooth ion heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3340)
 3340     format ('*  sawtooth ion  heating, watts/meter**3')
          write  (niterdb, 10) (qsawi(j)*convert,j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qsawi(j)        , j=1,nj)
       end if
c
c --- radiated power density
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3345)
 3345     format ('*  radiated power density, watts/meter**3')
          write  (niterdb, 10) (-ABS (qrad(j))*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (      qrad(j)         , j=1,nj)
       end if
c
c --- qohm,ohmic heating profile
c
       if (irwflag .eq. 0) then
          if (iterdsc .ne. 0)
     .    write  (niterdb, 3350)
 3350     format ('*  (electron) ohmic power density, watts/meter**3')
          write  (niterdb, 10) (qohm(j)*convert, j=1,nj)
       else
          read   (niterdb,  7) starflag
          read   (niterdb, 10) (qohm(j)        , j=1,nj)
       end if
c
c --- convert vectors defined on the npsi (i.e., MHD) grid to corresponding
c --- quantities on the rho (i.e., transport) grid
c
      tension = 0.0               ! don't use tension option of tspline
      tmax    = 0.0               ! max allowed tension
      bpar(1) = 0.0               ! set boundary conditions on spline
      bpar(2) = 0.0
      bpar(3) = 0.0
      bpar(4) = 0.0
c
c --- take care of a roundoff problem
c --- changed to relative measure HSJ 8/25/98
c
      if  (psir(nj) .ne. 0.0) then
        If (ABS((psival(1)-psir(nj))/(psir(nj))) .lt. 1.0e-10)
     .                                     psival(1) = psir(nj)
      else
        If (ABS(psival(1)-psir(nj)) .lt. 1.0e-10)
     .                                     psival(1) = psir(nj)
      end if
c
c --- average major radius
c
       if (irwflag .eq. 0) then
              call tspline (psival,rmajavnpsi,npsi,bpar,cs2spline,kpsi,
     .                      ier,tension,aspline,bspline,cspline,dspline,
     .                      espline,fspline,tmax,psir,work,nj)
             if (iterdsc .ne. 0)
     .       write  (niterdb, 1100)
 1100        format ('*  average major radius of each flux surface,',
     .               ' meters, evaluated at elevation of magnetic axis')
             do j=1,nj
               work(j) = 1.0e-2*work(j)     ! convert to meters
             end do
             write  (niterdb, 10) (work(j)      , j=1,nj)
        else
             read   (niterdb,  7) starflag
             read   (niterdb, 10) (rmajavnpsi(j), j=1,nj)
        end if
c
c --- average minor radius
c
        if (irwflag .eq. 0) then
           call tspline (psival,rminavnpsi,npsi,bpar,cs2spline,kpsi,
     .                   ier,tension,aspline,bspline,cspline,dspline,
     .                   espline,fspline,tmax,psir,work,nj)
           if (iterdsc .ne. 0)
     .     write  (niterdb, 1110)
 1110      format ('*  average minor radius of each flux surface, '
     .                'meters, evaluated at elevation of magnetic axis')
           do j=1,nj
             work(j) = 1.0e-2*work(j)       ! convert to meters
           end do
           write (niterdb, 10) (work(j)      , j=1,nj)
         else
           read  (niterdb,  7) starflag
           read  (niterdb, 10) (rminavnpsi(j), j=1,nj)
         end if
c
c --- volume
c
         if (irwflag .eq. 0) then
             call tspline (psival,psivolp,npsi,bpar,cs2spline,kpsi,
     .                     ier,tension,aspline,bspline,cspline,dspline,
     .                     espline,fspline,tmax,psir,work,nj)
             if (iterdsc .ne. 0)
     .       write  (niterdb, 1120)
 1120        format ('*  volume of each flux surface, meters**3')
             do j=1,nj
                work(j) = work(j) * 1.0e-6         ! convert to meter**3
             end do
             write (niterdb, 10) (work(j)   , j=1,nj)
        else
             read  (niterdb,  7) starflag
             read  (niterdb, 10) (psivolp(j), j=1,nj)
        end if
c
c ---  elongation
c
        if (irwflag .eq. 0) then
             call tspline (psival,elongx,npsi,bpar,cs2spline,kpsi,
     .                     ier,tension,aspline,bspline,cspline,dspline,
     .                     espline,fspline,tmax,psir,work,nj)
             if (iterdsc .ne. 0)
     .       write  (niterdb, 1130)
 1130        format ('*  elongation of each flux surface')
             write  (niterdb, 10) (work(j)  , j=1,nj)
          else
             read   (niterdb,  7) starflag
             read   (niterdb, 10) (elongx(j), j=1,nj)
          end if
c
c --- triangularity
c
          if (irwflag .eq. 0) then
            call tspline (psival,triangnpsi,npsi,bpar,cs2spline,kpsi,
     .                    ier,tension,aspline,bspline,cspline,dspline,
     .                    espline,fspline,tmax,psir,work,nj)
            if (iterdsc .ne. 0)
     .      write  (niterdb, 1140)
 1140       format ('*  triangularity of each flux surface')
            write  (niterdb, 10) (work(j)      , j=1,nj)
          else
            read   (niterdb,  7) starflag
            read   (niterdb, 10) (triangnpsi(j), j=1,nj)
          end if
c
c --- indentation
c
          if (irwflag .eq. 0) then
            call tspline (psival,pindentnpsi,npsi,bpar,cs2spline,kpsi,
     .                    ier,tension,aspline,bspline,cspline,dspline,
     .                    espline,fspline,tmax,psir,work,nj)
            if (iterdsc .ne. 0)
     .      write  (niterdb, 1150)
 1150       format ('*  indentation of each flux surface')
            write  (niterdb, 10) (work(j)       , j=1,nj)
          else
            read   (niterdb,  7) starflag
            read   (niterdb, 10) (pindentnpsi(j), j=1,nj)
          end if
c
c --- surface area
c
          if (irwflag .eq. 0) then
c
****        call tspline (psival,sfareanpsi,npsi,bpar,cs2spline,kpsi,
****  .                   ier,tension,aspline,bspline,cspline,dspline,
****  .                   espline,fspline,tmax,psir,work,nj)
c
c           change definition of surface area calc
c           from using sfareanpsi to using
c
c               4 * pi * pi * R0 * h * r * ABS (grad rho)
c
c           to be consistent with other laboratories:
c
            call tspline (psival,grho1npsi,npsi,bpar,cs2spline,kpsi,
     .                    ier,tension,aspline,bspline,cspline,dspline,
     .                    espline,fspline,tmax,psir,work,nj)
            if (iterdsc .ne. 0)
     .      write  (niterdb, 1160)
 1160       format ('*  surface area of each flux surface, meter**2' /
     .              '*  this is 4*pi*pi*R0*hcap*rho*<ABS(grad rho)>;',
     .                ' see file "iter1.ps" for explanation')
            do j=1,nj
              work(j) = work(j) * 1.0e-4         ! convert to meter**2
            end do
****        write (niterdb, 10) (work(j)      , j=1,nj)
            write (niterdb, 10) (39.478 * rmajor * hcap(j) * r(j)
     .                    * ABS (work(j))     , j=1,nj)
          else
            read  (niterdb,  7)  starflag
            read  (niterdb,  7)  starflag
            read  (niterdb, 10) (sfareanpsi(j), j=1,nj)
          end if
c
c --- cross-sectional area
c
           if (irwflag .eq. 0) then
             call tspline (psival,cxareanpsi,npsi,bpar,cs2spline,kpsi,
     .                     ier,tension,aspline,bspline,cspline,dspline,
     .                     espline,fspline,tmax,psir,work,nj)
             if (iterdsc .ne. 0)
     .       write  (niterdb, 1165)
 1165        format ('*  cross-sectional area of each flux'
     .                 ' surface, meters**2')
             do j=1,nj
               work(j) = 1.0e-4 * work(j)
             end do
             write (niterdb, 10) (      work(j), j=1,nj)
           else
             read  (niterdb,  7) starflag
             read  (niterdb, 10) (cxareanpsi(j), j=1,nj)
           end if
c
c --- flux surface average grad rho
c
            if (irwflag .eq. 0) then
              call tspline (psival,grho1npsi,npsi,bpar,cs2spline,kpsi,
     .                      ier,tension,aspline,bspline,cspline,dspline,
     .                      espline,fspline,tmax,psir,work,nj)
              if (iterdsc .ne. 0)
     .        write  (niterdb, 1170)
 1170         format ('*  flux surface average absolute grad rho')
              write  (niterdb, 10) (work(j)     , j=1,nj)
            else
              read   (niterdb,  7) starflag
              read   (niterdb, 10) (grho1npsi(j), j=1,nj)
            end if
c
c --- flux surface average (grad rho)**2
c
            if (irwflag .eq. 0) then
                 call tspline (psival,grho2npsi,npsi,bpar,cs2spline,
     .                         kpsi,ier,tension,aspline,bspline,cspline,
     .                         dspline,espline,fspline,tmax,psir,
     .                         work,nj)
                 if (iterdsc .ne. 0)
     .           write  (niterdb, 1180)
 1180            format ('*  flux surface average (grad rho)**2')
                 write  (niterdb, 10) (work(j)     , j=1,nj)
            else
                 read   (niterdb,  7) starflag
                 read   (niterdb, 10) (grho2npsi(j), j=1,nj)
            end if
c
c --- plasma boundary
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 1200)
 1200   format ('*  nplasbdry : number of points on plasma boundary')
        write  (niterdb, 9) nplasbdry
        if (iterdsc .ne. 0)
     .  write  (niterdb, 1210)
 1210   format ('*  r points for plasma boundary, meters')
        write  (niterdb, 10) (rplasbdry(j), j=1,nplasbdry)
        if (iterdsc .ne. 0)
     .  write  (niterdb, 1215)
 1215   format ('*  z points for plasma boundary, meters')
        write  (niterdb, 10) (zplasbdry(j), j=1,nplasbdry)
      else
        read   (niterdb,  7)  starflag
        read   (niterdb,  9)  nplasbdry
        read   (niterdb,  7)  starflag
        read   (niterdb, 10) (rplasbdry(j), j=1,nplasbdry)
        read   (niterdb,  7)  starflag
        read   (niterdb, 10) (zplasbdry(j), j=1,nplasbdry)
      end if

c
c --- nbeam torque density - new input added here so that backward 
c     compatibility is maintained  
c
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3125)
 3125   format ('* beam   torque density, nt-m/m**3')
        write  (niterdb, 10) (0.1*storqueb(j), j=1,nj) !0.1 multiplier added 03/16/04 HSJ
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (storqueb(j), j=1,nj)
      end if

c  12/8/09 add pressb and press HSJ:
      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3225)
 3225   format ('* beam  pressure, nt/m**2')
        write  (niterdb, 10) (0.1*pressb(j), j=1,nj) !dyne/cm^2 o nt/m^2
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (pressb(j), j=1,nj)
        pressb(:) = pressb*10.
      end if

      if (irwflag .eq. 0) then
        if (iterdsc .ne. 0)
     .  write  (niterdb, 3325)
 3325   format ('* total  pressure, nt/m**2')
        write  (niterdb, 10) (0.1*press(j), j=1,nj) !dyne/cm^2 o nt/m^2
      else
        read   (niterdb,  7) starflag
        read   (niterdb, 10) (press(j), j=1,nj)
        press(:) = press*10.
      end if


c
 2000 call giveupus(niterdb)
      close (unit = niterdb, status = 'KEEP')
c
c     subsequent calls will use existing file since ICALL is SAVEd
c
      icall = icall + 1
      return
c
      end
