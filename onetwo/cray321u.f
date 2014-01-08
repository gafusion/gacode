
      subroutine adasqh6 (ti, ecol, bmz, iflag, ne, zeff, conche,
     .                    concbe, concb, concc, concn, conco, concne,
     .                    qrat, ierr)

c
c  Calculate and return the rate coefficient for beam stopping
c  This routine is a modified version of L.D. Horton's (JET)
c  QHIOCH6 routine used to access and evaluate the ADAS beam
c  stopping cross sections from ion-specific effective data files.
c
c   changes from QHIOCH6:
c      1) The need for the NAG library has been eliminated using
c         the spline routines SPLINE and SEVAL.
c      2) Only the full-energy component can be calculated by
c         setting iflag=1. This allows lower energies to be
c         read from the ion-specific files.
c
c                                  D.F. Finkenthal 6-JUL-95
c
c   changes from QHIOCH5:
c      1) modified to force rereading of input files when beam
c         species has changed from the last call (using ipass)
c
c                                  L.D. Horton    2 June 92
c
c   changes from QHIOCH4:
c      1) modified to accept beam species as input and to then
c         calculate the beam stopping for either hydrogen or
c         helium beams
c
c                                  L.D. Horton   16 August 91
c
c   changes from QHIOCH3:
c      1) incorporated Fritsch's new calculation for excitation by
c         Be4+ and He2+
c      2) improved calculation of Maxwellian averages in bundled-n
c         code
c      3) proper inclusion of beam energy dependence
c          - stopping rate is now read from a matrix on beam energy
c            and target density; only the temperature dependence is
c            done separately
c      4) included boron, nitrogen and neon as input concentrations
c          - the code skips all zero concentrations
c
c                                  L.D. Horton   31 July 91
c
c   also:  back to Wilhelm's 3 energies-at-a-time so that all spline
c fits can be done at once. Only if the beam energy has changed by
c more than 1% will the density splines be redone.  Since only one
c energy is used per bank, this means that the spline will be done
c only twice per call of ATTS4.
c
c                                  L.D. Horton   14 August 91
c
c ti        : REAL   : ion temperature in eV
c ecol      : REAL   : collision (=beam) energy in eV/amu
c bmz       : INTEGER: beam nuclear charge
c iflag     : INTEGER: flag for full energy calculation only
c ne        : REAL   : electron density in m**-3
c conche    : REAL   : relative Helium concentration   ( 0 < CHE < 1)
c concbe    : REAL   : relative Beryllium concentration
c concb     : REAL   : relative Boron concentration
c concc     : REAL   : relative Carbon concentration
c concn     : REAL   : relative Nitrogen concentration
c conco     : REAL   : relative Oxygen concentration
c concne    : REAL   : relative Neon concentration
c qrat(3)   : REAL   : requested cross section m**2 for full, half,
c                      and third energy beam components
c
c

      USE param
      USE io 
      USE ext_prog_info, only : nchars_12,onetwo_xsct
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'io.i'

c
c ipass : file read switch.  Reread files if beam species has changed
c
      integer    ipass
      character  dsn(8,2)*35
      data       ipass/0/
      data       dsn  / '/data/adas/h_h1.dat'  ,
     .                  '/data/adas/h_he2.dat' ,
     .                  '/data/adas/h_be4.dat' ,
     .                  '/data/adas/h_b5.dat'  ,
     .                  '/data/adas/h_c6.dat'  ,
     .                  '/data/adas/h_n7.dat'  ,
     .                  '/data/adas/h_o8.dat'  ,
     .                  '/data/adas/h_ne10.dat',
     .                  '/data/adas/he_h1.dat' ,
     .                  '/data/adas/he_he2.dat',
     .                  '/data/adas/he_be4.dat',
     .                  '/data/adas/he_b5.dat' ,
     .                  '/data/adas/he_c6.dat' ,
     .                  '/data/adas/he_n7.dat' ,
     .                  '/data/adas/he_o8.dat' ,
     .                  '/data/adas/he_ne10.dat'/
c
c Physics Constants
c
      real*8     amu, eV
      parameter (amu = 1.6605e-24)
      parameter (eV  = 1.6022e-12)
c
c Local variables
c
c      integer    bmz, iflag, nebeam
      integer    bmz, iflag
      integer,save :: nebeam
      real*8     ti, ecol, ne, qrat(3)
      real*8     conche, concbe, concb, concc, concn, conco, concne
****  real*8     conch
      integer    nsp, maxe, maxn, maxt, ierr
      parameter (nsp = 8)              ! 8 different ion species
      parameter (maxe = 15, maxn = 10, maxt = 10)
      integer    isp, z(nsp), neb(nsp), ndens(nsp), ntemp(nsp), i, j, k
      real*8     eb(maxe,nsp), dens(maxn,nsp), temp(maxt,nsp)
      real*8     tref(nsp), ebref(nsp), denref(nsp), svref(nsp)
      real*8     sven(maxe,maxn,nsp), svt(maxt,nsp), seval
      character  line*80
      data       z /1, 2, 4, 5, 6, 7, 8, 10/
c
c      real*8     ti8, ecol8, ne8, qrat8(3)
      real*8     ti8, ne8, qrat8(3)
      real*8,    save  :: ecol8
      real*8     conc(nsp)
c
      integer    ifail
      real*8     be(maxe,maxn,nsp),ce(maxe,maxn,nsp),de(maxe,maxn,nsp)
      real*8     bt(maxt,nsp),ct(maxt,nsp),dt(maxt,nsp)
      real*8     bn(maxn,nsp,3),cn(maxn,nsp,3),dn(maxn,nsp,3)
      real*8     svintn(maxn,nsp,3),svintt,svtot(nsp),svtcor(nsp)
c
      real*8 zeffm1
      real*8 vbeam
c
      ierr    = 0
      ti8     = ti
      zeff8   = zeff
      ne8     = ne * 1.0e-13
****  conc(1) = conch
      conc(2) = conche
      conc(3) = concbe
      conc(4) = concb
      conc(5) = concc
      conc(6) = concn
      conc(7) = conco
      conc(8) = concne
c
c open and read input file only once
c
      if (ipass .ne. bmz) then
        write (ncrt, 1100)
        write (nout, 1100)
c
        ecol8 = ecol
c
        do isp=1,nsp
          call getioun(nunadas,nunadas)
          print *,'file =',onetwo_xsct(1:nchars_12)//dsn(isp,bmz) !jmp.den
          open (unit = nunadas,
     .          file = onetwo_xsct(1:nchars_12)//dsn(isp,bmz), 
     .                                         status = 'OLD')
          read (nunadas,  1000) z(isp),svref(isp)
          read (nunadas, '(a)') line
          read (nunadas,  1001) neb(isp),ndens(isp),tref(isp)
          read (nunadas, '(a)') line
          read (nunadas,  1002) (eb(j,isp),j=1,neb(isp))
          read (nunadas,  1002) (dens(k,isp),k=1,ndens(isp))
          read (nunadas, '(a)') line
          do k=1,ndens(isp)
            read (nunadas, 1002) (sven(j,k,isp),j=1,neb(isp))
          end do
          read  (nunadas, '(a)') line
          read  (nunadas,  1003) ntemp(isp),ebref(isp),denref(isp)
          read  (nunadas, '(a)') line
          read  (nunadas,  1002) (temp(j,isp),j=1,ntemp(isp))
          read  (nunadas, '(a)') line
          read  (nunadas,  1002) (svt(j,isp),j=1,ntemp(isp))
          call giveupus(nunadas)
          close (nunadas)
c
          do k = 1, ndens(isp)
            dens(k,isp) = dens(k,isp) * 1.0e-13
          end do
c
c spline the data on energy and temperature
c
          do k=1,ndens(isp)
            ifail = 0
            call spline_12 (neb(isp),eb(1,isp),sven(1,k,isp),
     .                   be(1,k,isp),ce(1,k,isp),de(1,k,isp))
            if (ifail .ne. 0) then
              write (6, '(a)')  ' spline error #1 in ADASQH6'
              ierr = 1
              return
            end if
          end do
          ifail = 0
          call spline_12 (ntemp(isp),temp(1,isp),svt(1,isp),
     .                 bt(1,isp),ct(1,isp),dt(1,isp))
          if (ifail .ne. 0) then
            write (6, '(a)')  ' spline error #2 in ADASQH6'
            ierr = 1
            return
          end if
c
c Determine if only the full energy component is requested
c Default is that all three beam energy components are to be determined
c Only the full energy component is interesting for He beams
c
          nebeam = 3
          if (iflag .gt. 0)  nebeam = 1
          if (  bmz .eq. 2)  nebeam = 1
c
c spline on density for each requested energy component (1 or all 3)
c
          do   i=1,nebeam
            do k=1,ndens(isp)
              svintn(k,isp,i) = seval(neb(isp),ecol8/i,eb(1,isp),
     .             sven(1,k,isp),be(1,k,isp),ce(1,k,isp),de(1,k,isp))
c
            end do
c
            ifail = 0
            call spline_12(ndens(isp),dens(1,isp),svintn(1,isp,i),
     .                  bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
c
            if (ifail .ne. 0) then
              write (6, '(a)')  ' spline error #3 in ADASQH6'
              ierr = 1
              return
            end if
          end do
        end do
c
 1000   format (i5, 8x, d9.3)
 1001   format (2i5, 7x, d9.3)
 1002   format (8(1x, d9.3) / 8(1x, d9.3))
 1003   format (i5, 7x, d9.3, 7x, d9.3)
 1100   format (/
     . ' Using ADAS, the effective beam stopping cross sections:' /
     . ' ***** Opening and Reading the Atomic Data Tables *****'  /)
c
        ipass = bmz
      else
c
c  redo density splines only if beam energy has changed by more than 1%
c
        if (ABS (ecol8-ecol) / ecol8 .gt. 0.01) then
          do isp=1,nsp
            ecol8 = ecol
            do i=1,nebeam
              do k=1,ndens(isp)
                svintn(k,isp,i) = seval(neb(isp),ecol8/i,eb(1,isp),
     .                sven(1,k,isp),be(1,k,isp),ce(1,k,isp),de(1,k,isp))
              end do
              ifail = 0
              call spline_12 (ndens(isp),dens(1,isp),svintn(1,isp,i),
     .                     bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
              if (ifail .ne. 0) then
                write (6, '(a)')  ' spline error #4 in ADASQH6'
                ierr = 1
                return
              end if
            end do
          end do
        end if
      end if
c
c calculate correction to requested temperature
c
      do isp=1,nsp
        if (ti8 .le. temp(1,isp)) then
          svtcor(isp) = svt(1,isp)/svref(isp)
        else
          svintt      = seval(ntemp(isp),ti8,temp(1,isp),
     .                  svt(1,isp),bt(1,isp),ct(1,isp),dt(1,isp))
          svtcor(isp) = svintt/svref(isp)
        end if
      end do
c
c scale the input concentrations to match the required zeff
c
           zeffm1 = 0.0
           do isp=2,nsp
             zeffm1 = zeffm1 + z(isp)*(z(isp)-1)*conc(isp)
           end do
           if (zeffm1 .gt. 1.0e-5) then
             do isp=2,nsp
               conc(isp) = (zeff8-1.0) / zeffm1 * conc(isp)
             end do
           end if
           conc(1) = 1.0
           do isp=2,nsp
             conc(1) = conc(1) - z(isp)*conc(isp)
           end do
c
c evaluate at three energy components
c
      do i=1,nebeam
c
c interpolate to requested density
c
        do isp=1,nsp
          if (ne8 .le. dens(1,isp)) then
            svtot(isp) = seval(ndens(isp),dens(1,isp),dens(1,isp),
     .            svintn(1,isp,i),bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
          else
            svtot(isp) = seval(ndens(isp),ne8,dens(1,isp),
     .            svintn(1,isp,i),bn(1,isp,i),cn(1,isp,i),dn(1,isp,i))
          end if
        end do
c
c  construct the total stopping cross section
c
        qrat8(i) = 0.0
        do isp=1,nsp
          qrat8(i) = qrat8(i) + svtot(isp)*svtcor(isp)*z(isp)*conc(isp)
        end do
c
c  divide by the beam speed to get a cross section and convert to m**2
c
        vbeam   = SQRT (2.0 * ecol8 / i * ev / amu)
        qrat(i) = qrat8(i) / vbeam
      end do
      return
c
      end
      subroutine adassgxn (ebkev, ebfac, ibion, mb, mfm1, nebin, vbeam,
     .                     zne, zni, zte, zti0, zzi, debin, sgxn,
     .                     sgxnmi, atw_beam)
c

c
c ----------------------------------------------------------------------
c
c ADASSGXN calculates the effective neutral beam stopping cross sections
c using the JET Atomic Data and Analysis Structure (ADAS). The cross
c sections are returned in array sgxn for each beam, beam energy
c component, and FREYA flux zone.
c
c The plasma ions are assumed to be fully stripped. Only ion species
c H, He, B, Be, C, O, N, and Ne are available from ADAS. Any other ion
c will terminate the code. This can be improved later.
c
c The routine was made to be as compatible as possible with the existing
c code. A call is made to the original cross section package in order to
c determine relative deposition fractions only. The HEXNB package is
c totally avoided.
c
c Reference:   Finkenthal, Ph.D. Thesis, UC-Berkeley, 1994
c
c Created  :   06-jul-1995    by  D. Finkenthal
c ----------------------------------------------------------------------
c
c     input:
c
c          ebkev(mb)      - full energy of mb-th neutral beamline
c                           (keV)
c          ebfac          - factor defining upper bound on energy bin
c                           range,  .gt. 1.0
c          ibion          - index of beam ion species
c                           (e.g. atwb = atw(ibion))
c                           if ibion = -1 beam is dt mixture
c          atw_beam         (use atw_beam for mass in this case)
c          mb             - number of beamlines modeled
c          mfm1           - number of flux zones minus one
c          nebin          - number of energy bins (rotation case only)
c          vbeam(ie,mb)   - speed of ie-th energy group of the mb-th
c                           beamline (cm/sec)
c          zne(mfm1)      - local electron density (cm**-3)
c          zni(mfm1,nion) - local density of nion-th ion species
c                           (cm**-3)
c          zte(mfm1)      - local electron temperature (KeV)
c          zti0(mfm1)     - local ion temperature (KeV)
c          zzi(mfm1,nion) - local average charge state of nion-th ion
c                           species
c     input from common /io/:
c          ncrt,nout
c     input from common /ions/:
c          namei, atw
c     input from common /numbrs/:
c          nprim,nimp,nion
c
c     output:
c          sgxn(i,j,k,l)
c               i - mode index
c                   = 1, fraction of reactions producing electrons;
c                   = 2, fraction of reactions producing species 1 ion;
c                   = 3, fraction of reactions producing species 2 ion;
c                   = 4, total inverse mean free path;
c               j - FREYA zone index
c               k - beam index
c                   k = 3*(ib-1) + ie, where ib is the beam index and ie
c                                      is the beam energy group index
c               l - index for relative energy.  bins are equispaced in
c                   delta(energy) between 0 and max(ebkev(ib))*ebfac,
c                   with bin width given by
c                             delta(energy) = ebmax*ebfac/nebin.
c                   ebfac .ge. 1.0 and nebin are user supplied namelist
c                   variables.
c          sgxnmi(ie,mb)
c                 - minimum inverse mean free path for ie-th energy
c                   group of mb-th beamline.  calculated for stationary
c                   target case only.  otherwise, calculated for each
c                   neutral trajectory in subroutine INJECT.
c          debin
c                 - width of energy bins (keV/amu)
c                   Used only for rotatation case.
c ----------------------------------------------------------------------
c

      USE param
      USE ions
      USE io 
      USE numbrs
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character rcs_id*63
      save      rcs_id
      data      rcs_id /
     ."$Id: cray321u.f,v 1.6 2013/12/12 18:37:01 stjohn Exp $"/
c      include 'param.i'
c      include 'io.i'
c      include 'numbrs.i'
c      include 'ions.i'
c
      dimension ebkev(kb),vbeam(ke,kb),zne(kz),zni(kz,kion)
      dimension zte(kz),zti0(kz),zzi(kz,kion)
      dimension sgxn(kcmp1,kz,kbe,ksge),sgxnmi(ke,kb)
c
      real*8 sgxeff(3), cnz(kz,20), zeffx(kz), izatom(kion)
      save init, izatom, atwb, izbm
      data init /0/
c
      if (init .eq. 0) then
c
c ----------------------------------------------------------------------
c determine atomic number of primary ions
c ----------------------------------------------------------------------
c
      do i=1,nprim
        izatom(i) = 1
        if (namep(i) .eq. 'he')
     .  izatom(i) = 2
      end do
c
c ----------------------------------------------------------------------
c determine atomic number of impurity ions
c ----------------------------------------------------------------------
c
      if (nimp .eq. 0)  go to 3430
      do i=1,nimp
        k = nprim + i
        izatom(k) = 0
        if (namei(i) .eq. 'he')  izatom(k) =  2
        if (namei(i) .eq. 'c' )  izatom(k) =  6
        if (namei(i) .eq. 'o' )  izatom(k) =  8
        if (namei(i) .eq. 'si')  izatom(k) = 14
        if (namei(i) .eq. 'ar')  izatom(k) = 18
        if (namei(i) .eq. 'cr')  izatom(k) = 24
        if (namei(i) .eq. 'fe')  izatom(k) = 26
        if (namei(i) .eq. 'ni')  izatom(k) = 28
        if (namei(i) .eq. 'kr')  izatom(k) = 36
        if (namei(i) .eq. 'mo')  izatom(k) = 42
        if (namei(i) .eq. 'w' )  izatom(k) = 74
      end do
c
c Get beam species
c
 3430 if (ibion .gt. 0) then
        atwb = atw(ibion)
        izbm = izatom(ibion)
      else
        atwb = atw_beam
        izbm = 1
      end if
c
c ----------------------------------------------------------------------
c Check to make sure that the impurities requested are compatable with
c ADAS (i.e., H, He, B, Be, C, O, N, or Ne). Terminate with error if an
c invalid impurity is listed.
c ----------------------------------------------------------------------
c
      do k=1,nion
        if (izatom(k) .gt. 10) then
          write (ncrt, 1010)
          write (nout, 1010)
          call STOP ('subroutine ADASSGXN: unallowed impurity', 181)
        end if
 1010 format
     .   (' *** Execution Terminated:'                                 /
     .    '     The ADAS database only contains atomic cross-sections' /
     .    '     for fully-stripped H, He, B, Be, C, O, N, and Ne ions.'/
     .    '     You must restrict your choice of impurity species to'  /
     .    '     these ions.')
      end do
c
      init = 1
      end if    ! init
c
c ----------------------------------------------------------------------
c Get original cross sections-Used to determine relative deposition
c fractions only. Total cross section is determined using ADAS.
c ----------------------------------------------------------------------
c
      call nbsgold (ebkev, ebfac, ibion, mb, mfm1, nebin, vbeam,
     .              zne, zni, zte, zzi, debin, sgxn, sgxnmi,
     .                                               atw_beam)
c
c ----------------------------------------------------------------------
c Set up the cnz arrays (concentrations of plasma ions) and Zeffx
c zni and zne are the (FREYA zone) densities of electron and ions:
c ----------------------------------------------------------------------
c
      do i=1, mfm1
        do k = 1, nion
          cnz(i,k) = 0.0
          if (izatom(k) .eq.  1)  cnz(i, 1) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  2)  cnz(i, 2) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  4)  cnz(i, 4) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  5)  cnz(i, 5) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  6)  cnz(i, 6) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  7)  cnz(i, 7) = zni(i,k)/zne(i)
          if (izatom(k) .eq.  8)  cnz(i, 8) = zni(i,k)/zne(i)
          if (izatom(k) .eq. 10)  cnz(i,10) = zni(i,k)/zne(i)
        end do
c
c The Zeff(i=1,...mfm1) array has been stored in zzi(i,nion+1)
        zeffx(i) = zzi(i,nion+1)
c
      end do
c
c ----------------------------------------------------------------------
c     stationary plasma case
c ----------------------------------------------------------------------
c
      if (nebin .eq. 0) then
        ierr = 0
        do    i=1,mfm1
          do ib=1,mb
            tiev = 1.0e3 * zti0(i)
            ecol = 1.0e3 * ebkev(ib) / (atwb)
            call adasqh6 (tiev,ecol,izbm,0,zne(i),zeffx(i),
     .                    cnz(i,2),cnz(i,4),cnz(i,5),cnz(i,6),
     .                    cnz(i,7),cnz(i,8),cnz(i,10),
     .                    sgxeff,ierr)
            if (ierr .eq. 1) then
              call STOP ('subroutine ADASSGXN: problem #1', 183)
            end if
            do ie = 1,ke
              ind  = ke*(ib-1) + ie
              sgxn(4,i,ind,1) = zne(i)*sgxeff(ie)
              sgxnmi(ie,ib)   = MAX (sgxnmi(ie,ib),sgxn(4,i,ind,1))
            end do
          end do
        end do
c
      else
c
c ----------------------------------------------------------------------
c     rotating plasma case
c ----------------------------------------------------------------------
c
        ebmax = 0.0
        do ib=1,mb
          ebmax = MAX (ebmax,ebkev(ib))
        end do
        ebmax = ebmax/atwb
        debin = ebmax*ebfac/FLOAT (nebin)
c
        do i=1,mfm1
          tiev  = 1.0e3*zti0(i)
          jreff = 0
          do j=1,nebin
            ebin = FLOAT (j) * debin
            ecol =     1.0e3 *  ebin
c
c  Only call ADAS if beam energy (ecol) is greater than 5 keV/amu.
c  Skip over the energy bins that are less than 5 keV/amu. These
c  low energy bins will be calculated using the old cross sections
c  scaled to the adas data at some reference energy.
c
            if (ecol .ge. 5000.0) then
              call adasqh6 (tiev,ecol,izbm,1,zne(i),zeffx(i),
     .                      cnz(i,2),cnz(i,4),cnz(i,5),cnz(i,6),
     .                      cnz(i,7),cnz(i,8),cnz(i,10),
     .                      sgxeff,ierr)
              if (ierr .eq. 1) then
                call STOP ('subroutine ADASSGXN: problem #2', 184)
              end if
c
c Calculate a scaling factor from the first available (ge. 5.0)
c   ADAS bin energy This will be used to scale the old cross
c   sections to fit in with ADAS for the low energy bins.
c
              if (jreff .eq. 0) then
                jreff     = j
                sgxnscale = zne(i)*sgxeff(1)/sgxn(4,i,1,j)
              end if
              do k=1,ke*mb
                sgxn(4,i,k,j) = zne(i)*sgxeff(1)
              end do
            end if
          end do
c
c Now scale the old cross sections to match the first available ADAS
c   datapoint for the energy bins below 5.0 keV/amu
c
          if (jreff .gt. 1) then
            do  j=1,jreff-1
             do k=1,ke*mb
               sgxn(4,i,k,j) = sgxn(4,i,k,j)*sgxnscale
             end do
            end do
          end if
        end do
      end if
      return
c
      end



      subroutine beam_init(task,taskl,beam_restart_file_length,
     .                     beamon,beamoff,btime,nbeams,nsourc,
     .                     time0,timmax,bctime,beam_end,beam_cycles)
c--------------------------------------------------------HSJ---------------
c    INPUT
c    task        string set to 'ignore file' or 'use_file'
c                'use_file' is not implemented yet however
c    beam_restart_file_length
c    beamon(i)  the time at which source 1 of beam i first turns on
c    beamoff(i) the time interval during  which beam i is on
c    btime(i)   the time interval during which beam i is off
c    nbeams
c    nsourc
c    time0 
c    timmax 
c    bctime()
c    beam_end(i)  the time beyond which  the beam is off (and never turns on
c                 again)
c    beam_cycles(i) alternative to beam_end ,beam_end takes precedence
c----------------------------------------------------------------------------
      USE param
      USE io 
      USE nub4
      implicit none
      character *(*) task
      integer taskl, l, ie, i, j, k,nj,kbe1
      integer beam_restart_file_length,nbeams,nsourc
      real *8 timmax_beam,dtp,dt,
     .        beamon(*),beamoff(*),bctime(*),beam_end(*),
     .        timenew,time0,dttt,btime(*),taus,timmax,
     .        beam_cycles(*)

c      include 'param.i'    !kbs,kj,
c      include 'nub4.i'
c      include 'io.i'
      include 'beam_plot_dat.i'

       if(task(1:taskl) == 'use file')then
          !read in everything
          read(nb_strt)timmax_beam, time0_beam, dtp, dt
          read(nb_strt)beam_thermal_cutoff
          read(nb_strt)method,beam_pulse_control
          read(nb_strt)n_pulse,nj,kj1,kbs1,kbe1,kb1,kt1
          read(nb_strt)nbeams,nsourc,nbe,nbeam_pts
          read(nb_strt)taus
          read(nb_strt)pbeamOn
          read(nb_strt)pbeamOff
          read(nb_strt)pssv
          read(nb_strt)pssvoff
          read(nb_strt)enbeam_part
          read(nb_strt)enbeam_part_nl
          read(nb_strt)wenbeam_part
          read(nb_strt)Qfi_part
          read(nb_strt)Rfi_part
          read(nb_strt)Rfe_part
          read(nb_strt)enbeam_tot  !zeroed below,not needed?
          read(nb_strt)enbeam_tot_nl
          read(nb_strt)Qfi_tot  !zeroed below,not needed?
          read(nb_strt)Rfi_tot  !zeroed below,not needed?
          read(nb_strt)Rfe_tot  !zeroed below,not needed?
          read(nb_strt)beam_intensity
          read(nb_strt)tau0_vlj
          read(nb_strt)tau0_vljoff
          read(nb_strt)vcrit
          read(nb_strt)vbeam
          read(nb_strt)vthi
          read(nb_strt)waveform(1),timplot(1)
          nbeam_pts=1
C         ERROR CHECKING THAT COULD NOT BE DONE UNTIL FILE IS READ IN:
c         the ending time value in file beam_restart_file must be the
c         same as the starting time value to be used in this run:
          if( ABS(timmax_beam - time0) .gt. 1.d-5)then
                     print *,' error, time0 in inone and time0 in'
                     print *, beam_restart_file(1:
     .                                   beam_restart_file_length)
                     print *,' must match '
                     print *,'     when beam_init_restart_file = 0'
                     print *,'     and beam_restart_file name is ',
     .                              ' specified'
                     print *,' to get self consistent results you must'
                     print *,' either change the time value in the'
                     print *,' namelist or regenerate the '
                     print *,' beam_restart_file '
                     CALL EXIT(1)
          endif

       else  
           do l=1,kb
             do ie=1,ke
                do i = 1,kbs
                   do j=1,kj
                      enbeam_part(j,i,ie,l) = 0.0d0
c                     enbeam_part(j,i,ie,l) = 0.0d0
                      wenbeam_part(j,i,ie,l) = 0.0d0
                      Qfi_part(j,i,ie,l) = 0.0d0
                      Rfi_part(j,i,ie,l) = 0.0d0
                      Rfe_part(j,i,ie,l) = 0.0d0
                      do k=1,kt
                         tau0_vlj(k,j,i,ie,l)=-1.0
                         tau0_vljoff(k,j,i,ie,l)=-1.0
                         pssv(k,j,i,ie,l)= .false.
                         pssvoff(k,j,i,ie,l)= .false.
                      enddo
                   enddo
                enddo
             enddo
          enddo

c      set beam_end from beam_cycles if beam_end was not specified
c      in namelist inone:
       do j= 1,nbeams
          if(beam_end(j) .le. 0.0 )then
             if (beam_cycles(j) .gt. 0.0)then
                beam_end(j) = beam_cycles(j)*(btime(j)+beamoff(j)) +
     .                                                 beamon(j)
             endif
          endif
       enddo
          tauslow = .150             !rough estimate of slowing down time
          call pulse_start_time (nbeams,beamon,time0,timmax,btime,
     .                     beamoff,tauslow,nsourc,source2_phase,
     .                     kt,kj,kbs,kb,bstime,beam_end,
     .                     pbeamOn,pbeamOff,bctime,n_pulse)

       endif

       if(task /= 'use file')then
           timenew=time0
       else             ! create a small first step(so plots will start right)
           dttt=0.0000001
           timenew = time0 + dttt
       endif

       return
       end



      subroutine beam_int(tau0_vu,tau0_vl,vb_l,vc_l,taus_l,enn_l,
     .                    emzrat_l,enf_int,enf_int_nl,
     .                    wenf_int,Qfi_int,Rfe_int,Rfi_int,
     .                                               no_fast_losses)
c-----------------------------------------------------------------------
c     integrate the beam distribution function to get the fast ion
c     stored energy density and related quantities.
c--------------------------------------------------------------HSJ------
      implicit none
      real *8 taus_l,vsqfunc,vfsq,
     .        pcx, cxinta,enf_int,wenf_int,tau0_vu,tau0_vl,
     .        vb_l,vc_l,enn_l,wenf,enf_int_nl,Qfi_int,vf,vf3,
     .        b,vc_l3,vb_l3,Rfe_int,Rfi_int,emzrat_l,enrgy_factor
      integer j
      include 'gauss_info.i'
      logical no_fast_losses

      common /cx_bblock/ vc_l3,enrgy_factor ! use vc_l3 here

c     generate Guass-legendre weights which will be used to
c     evaluate the integrals:
      call gauleg (tau0_vl,tau0_vu,tau0_x, tau0_w, ngauss)



         vb_l3 = vb_l*vb_l*vb_l
         wenf_int=0.0d0
         Qfi_int=0.0d0
         enf_int = 0.0d0
         Rfe_int = 0.0d0
         Rfi_int = 0.0d0
         enf_int_nl = (tau0_vu-tau0_vl)
         do j=1,ngauss                                  !evaluate the integrals
            wenf =  vsqfunc(taus_l,tau0_x(j))
            vfsq = vc_l*vc_l*wenf
            vf  = SQRT(vfsq) 
            vf3 = vf * vfsq
            b = (vb_l3 + vc_l3)/(vf3 + vc_l3)
            b = b**(emzrat_l/3.d0)
            if(no_fast_losses)then                      !charge exchange,fusion,etc
               wenf_int = wenf_int + wenf*tau0_w(j)
               Qfi_int  = Qfi_int  + tau0_w(j)/vf 
               Rfe_int  = Rfe_int  + tau0_w(j)*b*vf
               Rfi_int  = Rfi_int  + tau0_w(j)*b/vfsq
            else
               pcx = cxinta(vf,vb_l,taus_l,enn_l)            !charge exchange part
               enf_int   = enf_int  + pcx*tau0_w(j)
               wenf_int  = wenf_int + wenf*pcx*tau0_w(j)
               Qfi_int   = Qfi_int  + pcx*tau0_w(j)/vf
               Rfe_int   = Rfe_int  + pcx*b*vf*tau0_w(j)
               Rfi_int   = Rfi_int  + pcx*b*tau0_w(j)/vfsq
            endif
            
         enddo
      return
      end



        subroutine beam_plot(timenew,s0,bon,boff,ie)
        implicit none
        include "beam_plot_dat.i"
        integer ie
        real *8
     .         timenew,s0,bon,boff,tmeps,tpeps,eps2

c        print *,"nbeam_pts =",nbeam_pts

        eps = 1.e-12  
        eps2 =0.5*eps
        tmeps = timenew - eps
        tpeps = timenew + eps
        if(nbeam_pts .gt. 1 .and. nbeam_pts  .lt. nplt_max-1) then
           if(bon .gt. 0.0d0)then
                if(ABS(timplot(nbeam_pts)-timenew) .gt. eps2)then
                    !here a new pulse is just turning on
                    !before the turnon we have:
                    nbeam_pts = nbeam_pts + 1
                    timplot(nbeam_pts)= tmeps
                    waveform(nbeam_pts) = waveform(nbeam_pts -1)
                    !and after  the turnon we have:
                    nbeam_pts = nbeam_pts + 1
                    timplot(nbeam_pts)= timenew
                    waveform(nbeam_pts) = waveform(nbeam_pts-1) + s0
                else  ! dont change the time,just add up the pulses
                    !here different pulses have the same switching times
                    waveform(nbeam_pts) = waveform(nbeam_pts) + s0
                endif
           else if(boff .gt. 0.0d0)then
                if(ABS(timplot(nbeam_pts)-timenew) .gt. eps2)then
                    !here a new pulse is just turning off
                    !before the turnon we have:
                    nbeam_pts = nbeam_pts + 1
                    timplot(nbeam_pts)  = timenew 
                    waveform(nbeam_pts)=waveform(nbeam_pts -1)
                    !and after  the turnoff  we have:
                    nbeam_pts = nbeam_pts + 1
                    timplot(nbeam_pts)  = tpeps
                    waveform(nbeam_pts) = waveform(nbeam_pts-1) - s0
                else  ! dont change the time,just subtract out the pulses
                    !here different pulses have the same switching times
                    waveform(nbeam_pts) = waveform(nbeam_pts) - s0
                endif
           else
             call STOP ('subroutine BEAM_PLOT: unallowed bon,boff', 0)
           endif
        else                                ! initialization assumes beams 
           nbeam_pts = 1                    ! are initially off
           timplot(nbeam_pts) = tmeps
           waveform(nbeam_pts) = 0.0
           nbeam_pts = nbeam_pts + 1
           timplot(nbeam_pts)= tpeps
           waveform(nbeam_pts) = waveform(nbeam_pts-1) + s0
        endif
        if(s0 .lt. 0.0)then
           print *,"s0 < 0.0 not valid"
           call STOP ('subroutine BEAM_PLOT: unallowed no S0', 0)
        endif

        return
        end


      subroutine  beam_prop(atwb,eb,bpow,vbeam,bion,bneut,bntot)
c----------------------------------------------------------------------
      USE param, only : ke,kb
      USE nub, only: atw_beam,nbeams,ebkev,neg_ion_source,fbcur,bcur
      implicit none
C INPUT
c  from nub: atw_beam,nbeams,ebkev,neg_ion_source,fbcur,bcur


c OUTPUT through argument list, see description on freya:
c  atwb
c  bntot
c  eb
c  bion
c  bneut
c  bpow
c  vbeam


c LOCAL 
c ebx
c beff
c ebev
c  --------------------------------------------------------------------
      real *8 atwb,bntot,eb(ke,kb),bpow(ke,kb),ebx,beff,
     .        ebev,vbeam(ke,kb),bion(ke,kb),bneut(ke,kb)
      integer ib,ie
      atwb  = atw_beam
      bntot = 0.0
      do 20 ib=1,nbeams
      do 20 ie=1,3
      eb(ie,ib) = ebkev(ib)/ie
      ebx       = eb(ie,ib)/atwb
      if (neg_ion_source(ib) .gt. 0) then
        beff = 0.98           ! arbitrary efficiency HSJ
      else
        call logint (ebx, beff)
      end if
      ebev = 1.0e3*ebx
      vbeam(ie,ib) = 1.384e6 * SQRT (ebev)
      bion(ie,ib)  = 0.625e19*fbcur(ie,ib)*bcur(ib)
      bneut(ie,ib) = ie*beff*bion(ie,ib)
      bpow(ie,ib)  = eb(ie,ib)*bneut(ie,ib)/0.625e16
   20 bntot        = bntot + bneut(ie,ib)
  
      return
      end





      subroutine beam_time_dependance(taus,dt,vcrit,ebeam,vbeam,vthi,
     .                     time0,time,timmax,pbeamOn,pbeamOff,
     .                     enbeam_part,enbeam_tot,enbeam_part_nl,
     .                     enbeam_tot_nl,wenbeam_part,wenbeam_tot,
     .                     Qfi_part,Qfi_tot,Rfi_part,Rfe_part,
     .                     Rfi_tot,Rfe_tot,
     .                     kj,kbs,kbe,kb,beam_time_init,
     .                     kt,beam_thermal_cutoff,
     .                     tau0_vlj,tau0_vljoff,nbeams,nsourc,
     .                     neg_ion_source,beam_thermal_speed,
     .                     source_strength,nj,pssv,
     .                     pssvoff,method,therm_frac,
     .                     no_fast_losses,enntot,mass_beam,emzrat,Pf0,
     .                     nf_tot_source,pwf_tot_source)

c-----------------------------------------------------------------------------------------------
c     subroutine stages the beam integrals to be performed for the case where the beams
c     are represented by arbitrary square waves. The inherent approximation here is that
c     taus remains constant ( energy and particle confinement times are varying only slowly
c     during the fast ion slowing down time).    


c     This subroutine uses some F90 array initialization
c --- INPUT by Argument list
c     taus(j)      fast ion Spitzer slowing down time at grid point j
c     dt           current time step in seconds
c     vcrit(j,l)   critcal velocity, beam l, rho grid point j
c     beam_thermal_cutoff      
c     vthi         thermal ion speed
c     nbeams           number of beamlines
c     nsourc           number of sources per beamline
c     pssv,pssvoff logical arrays, must be initialized to false before the first time
c                  that this routine is called. Should not be changed by the user
c                  after the initialization.
c     enntot(j)       j=1...nj neutral density 1/cm**3
c     mass_beam(k) k=1,kb in grams
c     emzrat       passed through to sub eval_beam_integrals


c     pbeamOn(kt,kbs,kb) contains the power level switching times in secs, for source kbs,
c                    beamline kb, time index kt
c                    power level switching includes intitial turnon, 
c                    final turnoff, sources that are  on before the beginning
c                    of the simulation time (ie time0), and sources that 
c                    switch from one power level to another at some time instant.
c                    This subroutine  assumes that if the source switches from 
c                    power level s1 to level s2 at time ts  then the entry ts 
c                    appears in pbeamOn to bring in source level s2 and the entry ts 
c                    also appears in pbeamOff to turn off source level s1.
c     pbeamOff(kt,kbs,kb) same for pulse off times
c
c     source_strength(kt,kj,kbs,kbe,kb)
c     neg_ion_source(kb)   =0  for standard 3 energy source
c                          =1  for full energy neg. ion source




c---OUTPUT
c
c   Qfi_part     energy deposition to ions,kev/(cm3*sec)
c   Qfi_tot      Qfi_part summed over all beams,sources,energies
c   Rfi_tot
c   Rfe_tot
c   nf_tot_source(J), J=1,2..NJ 
C                 the particle source rate,#/(cm**3sec), at this time summed over beams
c                 energies and pulses.
c   pwf_tot_source(J), J=1,2..NJ 
c                 the fast ion power,kev/cm**3/sec, at this time, summed over beams
c                 energies and pulses.
c ----------------------------------------------------------------------------HSJ-2/8/02--------
      implicit none
      integer i,j,nj,ie,k,l,kj,kbs,kbe,kb,kt,nbeams,nsourc,
     .        nbtmax,nbe, method,
     .        neg_ion_source(kb),beam_time_init,beam_thermal_cutoff
      logical*1 done,pssv(kt,kj,kbs,kbe,kb),
     .        pssvoff(kt,kj,kbs,kbe,kb),no_fast_losses
      real *8 
     . taus(nj),dt,timestep,vcrit(kj,kb),vbeam(kbe,kb),time0,timmax,
     . enbeam_part(kj,kbs,kbe,kb),vcutoff, tau0_vlj(kt,kj,kbs,kbe,kb),
     . tau0_vljoff(kt,kj,kbs,kbe,kb),timenew,timetol,enbeam_tot(kj),
     . tau0_vb,tau0_vcutoff,source_strength(kt,kj,kbs,kbe,kb),
     . vthi(kj),pbeamOn(kt,kbs,kb),pbeamOff(kt,kbs,kb),
     . time,timeon,timeoff,tau0_vu,tau0_vl,contrib,contrib_en,
     . ftau0,vu,vl,dtslowb,dtnew,bon,boff,dtslowbon,dtslowboff,
     . enf_contrib,wenf_contrib,wenbeam_tot(kj),vb3,vc3,vmin3,vth3,
     . wenbeam_part(kj,kbs,kbe,kb),contrib_wen,ebeam(kbe,kb),
     . vthcut3,therm_frac,one_third,vmax3,previous_en,
     . previous_wen,enntot(kj),mass_beam(kb),enf_contrib_nl,
     . previous_en_nl,enbeam_tot_nl(kj),enbeam_part_nl(kj,kbs,kbe,kb),
     . contrib_en_nl,previous_Qfi,contrib_Qfi,Qfi_part(kj,kbs,kbe,kb),
     . Qfi_tot(kj),Qfi_contrib,emzrat(kj),Rfi_contrib,Rfe_contrib,Pf0,
     . previous_Rfi,previous_Rfe,contrib_Rfi,contrib_Rfe,
     . Rfi_part(kj,kbs,kbe,kb),Rfe_part(kj,kbs,kbe,kb),Rfi_tot(kj),
     . Rfe_tot(kj),beam_thermal_speed,eval_vthcut3,tau0_vcut,
     . tau0_cut_avg,nf_tot_source(kj),pwf_tot_source(kj)


      beam_time_init = 1
      timenew=time+dt
      one_third=1.d0/3.d0
      timetol =5.e-9           !make consistent ( slightly larger than tol in chekdt)
      nbtmax = kt
      bon =-1.
      boff = -1.
      nf_tot_source =0.0d0   !array 
      pwf_tot_source =0.0d0  !array

      PRINT *,' '
      print *,'***********starting beam_time_dependance****************'
      print *,'time,timmax,timenew,dt',time,timmax,timenew,dt



      do l=1,nbeams                  !loop over  all beamlines
         if( beam_thermal_cutoff .eq. -1)
     .     tau0_vcut = tau0_cut_avg(taus,vcrit,kj,kb,nj,l)
         do i=1,nsourc               !loop over nsourc  sources per beamline
            if( timmax .lt. pbeamOn(1,i,l) )go to 100 ! this source is not used in current run
            nbe =3 
            if(neg_ion_source(l) .eq. 1)nbe =1
            do ie = 1,nbe         !loop over energies for each beam line
               vb3 = vbeam(ie,l)**3

               do j =1,nj        !loop over transport grid
                 vc3 = vcrit(j,l)**3
                 vth3 = vthi(j)**3
                 vthcut3 = eval_vthcut3(j)
                 if( beam_thermal_cutoff .eq. -1) then
                      tau0_vcutoff = tau0_vcut
                 else
                      tau0_vcutoff = one_third*taus(j)*
     .                                         LOG((vthcut3+vc3)/vc3)
                 endif
                 !ion with initial speed vbeam will slow to vthcut speed in time tau0_vb;
                 tau0_vb = one_third*taus(j)*LOG((vb3+vc3)/vc3)
                 do k=1,nbtmax    ! loop over input  time values in pbeamOn(k,i,l)
                    contrib_en=0.0d0
                    contrib_en_nl =0.0d0
                    contrib_wen=0.0d0
                    contrib_Qfi = 0.0d0
                    contrib_Rfi = 0.0d0
                    contrib_Rfe = 0.0d0
                    previous_en_nl = enbeam_part_nl(j,i,ie,l )
                    previous_en    = enbeam_part(j,i,ie,l )
                    previous_wen   = wenbeam_part(j,i,ie,l )
                    previous_Qfi   = Qfi_part(j,i,ie,l )
                    previous_Rfi   = Rfi_part(j,i,ie,l )
                    previous_Rfe   = Rfe_part(j,i,ie,l )
                    if(pbeamOn(k,i,l) - timenew  .ge. timetol)go to 200

                    if(.not. pssvoff(k,j,i,ie,l) ) then  ! if pulse k is contributing
                     if(ABS(timenew - pbeamOn(k,i,l)) .lt. timetol)then
                       tau0_vlj(k,j,i,ie,l) = tau0_vb  !lower limit of integral at this time
                       tau0_vu = tau0_vb               !equals upper limit
                       if(j .lt.2) then
                        print *, 'beam just ON  SECTION,l,i,ie =',l,i,ie
                        write(*,'("source_strength,k,j,i,ie,l ",
     .                     1pe14.6,2x,5(2x,i4))')
     .                   source_strength(k,j,i,ie,l),k,j,i,ie,l
                        print *,'tau0_vu',tau0_vu
                    print*, ' j up to 1,timenew,k,i,j,ie,l =',timenew,k,
     .                         i,j,ie,l
                       endif
                       bon = 1.0
                       boff = 0.0
                       call beam_plot(timenew,
     .                             Source_strength(k,j,i,ie,l),
     .                               bon,boff,ie)
                       nf_tot_source(j) = nf_tot_source(j) + 
     .                                  Source_strength(k,j,i,ie,l)
                       pwf_tot_source(j)  = pwf_tot_source(j) + 
     .                         Source_strength(k,j,i,ie,l)*ebeam(ie,l) !kev/(cm**3sec)
                 else  if ( timenew .gt. pbeamOn(k,i,l) .and.
     .                  timenew .le. pbeamOff(k,i,l) ) then !use .ge. instead of .gt. to set initial limits
                       done=.false.
                       !get oontribution from this time step due to source that was
                       !turned on before the current time:
                       timeon = timenew-pbeamOn(k,i,l)
                       timestep = MIN(dt,timeon)
                       if( .not. pssv(k,j,i,ie,l))then
                            if(tau0_vlj(k,j,i,ie,l) .ge. 0.0)then
                                 !set upper limit to previous lower limit:
                                 tau0_vu = tau0_vlj(k,j,i,ie,l)
                            endif
                            !set lower limit of integration
                            vmin3 = (vb3+vc3)*exp(-3.d0*timeon/
     .                                           taus(j))-vc3
                            vmin3=MAX(vmin3,vthcut3)
                            if(ABS(vmin3 - vthcut3 ) .gt. 1.e-12)then
                                 tau0_vl= one_third*taus(j)*
     .                                           LOG((vmin3+vc3)/vc3)
                            else
                                 tau0_vl = tau0_vcutoff
                            endif
                    if(j .eq. 1 .and. k .eq. 1)
     .              print *,'taus(1),LOG(vmin3+vc3)/vc3)',taus(1),
     .                 LOG((vmin3+vc3)/vc3)
                      tau0_vlj(k,j,i,ie,l)=tau0_vl
                    if(tau0_vu .gt. tau0_vl)then
                            call eval_beam_integrals(tau0_vu,
     .                            tau0_vl,Source_strength(k,j,i,ie,l),
     .                            no_fast_losses,taus(j),vcrit(j,l),
     .                            vbeam(ie,l),ebeam(ie,l),enntot(j),
     .                            mass_beam(l),emzrat(j),Pf0,
     .                            enf_contrib,enf_contrib_nl,
     .                            wenf_contrib,Qfi_contrib,
     .                            Rfi_contrib,Rfe_contrib)
                    if(j*k*i*l .eq. 1 )then
       write(*,'("tau0_vu,tau0_vl,enf_contrib_nl,tnew,ie = ",
     .                                       4(1pe14.6,2x),i5)')
     .              tau0_vu,tau0_vl,enf_contrib_nl,timenew,ie
       write(*,'("tau0_vu,tau0_vl,vmin3,vc3,ie = ",
     .                                       4(1pe14.6,2x),i5)')
     .              tau0_vu,tau0_vl,vmin3,vc3,ie
                    endif
                    nf_tot_source(j) = nf_tot_source(j) + 
     .                               Source_strength(k,j,i,ie,l)
                       pwf_tot_source(j)  = pwf_tot_source(j) + 
     .                         Source_strength(k,j,i,ie,l)*ebeam(ie,l)      !kev/(cm**3sec)
                  else
                     print *,'tau0_vu .lt. tau0_vl'
                     enf_contrib =0.0
                      enf_contrib_nl =0.0
                      wenf_contrib =0.0
                      Qfi_contrib =0.0
                      Rfe_contrib =0.0
                      Rfi_contrib =0.0


                  endif
                            contrib_en = enf_contrib
                            contrib_en_nl = enf_contrib_nl
                            contrib_wen = wenf_contrib
                            contrib_Qfi = Qfi_contrib
                            contrib_Rfe = Rfe_contrib
                            contrib_Rfi = Rfi_contrib
                        endif                              !end .not. pssv branch
                        if(ABS(timenew -pbeamOff(k,i,l)) .lt. timetol)
     .                          pssv(k,j,i,ie,l)=.true.
                        !add contribution from this time step to total from previous time step:
                        enbeam_part(j,i,ie,l)  = enbeam_part(j,i,ie,l) 
     .                                              + contrib_en
                        enbeam_part_nl(j,i,ie,l) = 
     .                         enbeam_part_nl(j,i,ie,l) + contrib_en_nl
                        wenbeam_part(j,i,ie,l)  = wenbeam_part(j,i,ie,l) 
     .                                        + contrib_wen
                        Qfi_part(j,i,ie,l) = Qfi_part(j,i,ie,l)  
     .                                        + contrib_Qfi
                        Rfi_part(j,i,ie,l) = Rfi_part(j,i,ie,l)  
     .                                        + contrib_Rfi
                        Rfe_part(j,i,ie,l) = Rfe_part(j,i,ie,l)  
     .                                        + contrib_Rfe
                        if(j .lt. 2)then
                        print *,'timenew =',timenew
                        print *,'tau0_vlj(k,j,i,ie,l) =',
     .                      tau0_vlj(k,j,i,ie,l)
                         print *,'enbeam_part(j,i,ie,l)= ',
     .                        enbeam_part(j,i,ie,l),j,i,ie,l
                         print *,'contrib_en = ',contrib_en
                         print *,'pbeamOn  k,j,i,ie,l =',k,j,i,ie,l
                         print *,'contrib_en_nl=',contrib_en_nl
                         print *,'enf_contrib_nl=',enf_contrib_nl
                         print *,'previous_en_nl=',previous_en_nl
                        print *,'pssv(k,j,i,ie,l)  ',pssv(k,j,i,ie,l)
                        if(pssv(k,j,i,ie,l))then
                           print *,'contrib_en_nl=',contrib_en_nl
                           print*,'enbeam_part_nl(j,i,ie,l )=',
     .                                  enbeam_part_nl(j,i,ie,l )
                        endif
                        print *,' '
                        endif
                     endif      !end  time .ge. pbeamOn branch



c
c------------------------beamOff section starts here----------------------------------
c




c                   if( ABS(timenew - pbeamOff(k,i,l)) .lt. timetol )then !beam is turned off at this time  HSJ 2/22/02 replaced with following
                   if( ABS(time - pbeamOff(k,i,l)) .lt. timetol )then !beam is turned off at this time
c                      previous_en    = enbeam_part(j,i,ie,l )
c                      previous_en_nl = enbeam_part_nl(j,i,ie,l )
c                      previous_wen   = wenbeam_part(j,i,ie,l )
cc                     previous_Qfi   = Qfi_part(j,i,ie,l )
c                      previous_Rfi   = Rfi_part(j,i,ie,l )
c                      previous_Rfe   = Rfe_part(j,i,ie,l )
                       tau0_vljoff(k,j,i,ie,l)=tau0_vlj(k,j,i,ie,l)
                       tau0_vlj(k,j,i,ie,l) = tau0_vb ! both lower and upper limits require monitoring
                       !if beam was on long enough to be fully developed then
                       !tau0_vljoff = tau0_vcutoff otherwise tau0_vljoff > tau0_vcutoff
                       !since beam is just turned off no contribution is calculated:
                       contrib_en=0.0d0
                       contrib_en_nl =0.0d0
                       contrib_wen=0.0d0
                       contrib_Qfi = 0.0d0
                       contrib_Rfi = 0.0d0
                       contrib_Rfe = 0.0d0
                       if(j .lt. 2)then
                       print *,
     .                 'beam just OFF section enbeam_part(j,i,ie,l)= ',
     .                                enbeam_part(j,i,ie,l)
                       print*, 'timenew,k,i=',timenew,k,i
                        print *,'tau0_vljoff(k,j,i,ie,l) =',
     .                      tau0_vljoff(k,j,i,ie,l)
                       endif
                       bon = 0.0
                       boff = 1.0
                       call beam_plot(time,
     .                              Source_strength(k,j,i,ie,l),
     .                               bon,boff,ie)
                  else if(timenew .gt. pbeamOff(k,i,l))then
                       !get oontribution from this time step due to source that was
                       !turned off before the current time:
                       timeoff = timenew-pbeamOff(k,i,l)
                       timeon  = timenew-pbeamOn(k,i,l)
                       timestep = MIN(dt,timeoff)
                       contrib_en = 0.0d0
                       contrib_wen = 0.0d0
                       contrib_en_nl = 0.0d0
                       contrib_Qfi =0.0d0
                       contrib_Rfi = 0.0d0
                       contrib_Rfe = 0.0d0
                       !set lower limit of integration
                       vmin3 = (vb3+vc3)*exp(-3.d0*timeon/taus(j))-vc3
                       vmin3=MAX(vmin3,vthcut3)
                       if(ABS(vmin3 - vthcut3 ) .gt. 1.e-12)then
                            tau0_vl= one_third*taus(j)*
     .                                           LOG((vmin3+vc3)/vc3)
                       else
                            tau0_vl = tau0_vcutoff
                       endif
                       tau0_vljoff(k,j,i,ie,l) = tau0_vl
                       !decrease upper limit of integration due to dt time elapse
                       vmax3 = (vb3+vc3)*exp(-3.d0*timeoff/taus(j))-vc3
                       vmax3=MAX(vmax3,vthcut3)
                       if(ABS(vmax3 - vthcut3 ) .gt. 1.e-12)then
                            tau0_vu= one_third*taus(j)*
     .                                           LOG((vmax3+vc3)/vc3)
                       else
                            tau0_vu = tau0_vcutoff
                       endif
                       tau0_vlj(k,j,i,ie,l) = tau0_vu
                       !calculate the amount that was lost during this time step
                       if(tau0_vu .gt. tau0_vl)then
                            call eval_beam_integrals(tau0_vu,
     .                             tau0_vl,Source_strength(k,j,i,ie,l),
     .                             no_fast_losses,taus(j),vcrit(j,l),
     .                             vbeam(ie,l),ebeam(ie,l),enntot(j),
     .                             mass_beam(l),emzrat(j),Pf0,
     .                             enf_contrib,enf_contrib_nl,
     .                             wenf_contrib,Qfi_contrib,
     .                             Rfi_contrib,Rfe_contrib)
                                   contrib_en =  enf_contrib
                                   contrib_en_nl = enf_contrib_nl
                                   contrib_wen = wenf_contrib
                                   contrib_Qfi = Qfi_contrib
                                   contrib_Rfe =Rfe_contrib
                                   contrib_Rfi =Rfi_contrib
                            if(j .lt. 2)then
                                   print *,'pbeamOff timenew,k,i =',
     .                                                   timenew,k,i
                                   print *,'Otau0_vu, tau0_vl =',
     .                                    tau0_vu,tau0_vl
                            endif
                        else 
                            pssvoff(k,j,i,ie,l) = .true. !pulse k can no longer contributute
                            contrib_en = 0.0d0
                            contrib_en_nl = 0.0d0
                            contrib_wen = 0.0d0
                            contrib_Qfi = 0.0d0
                            contrib_Rfi = 0.0d0
                            contrib_Rfe = 0.0d0
                            if(j .lt. 2)then
        write(6,'("k,i,ie,timenew,dt,pssvoff",3(2x,i3),
     .                                   2(2x,1pe14.6),2x,l1)')
     .                k,i,ie, timenew,dt,pssvoff(k,j,i,ie,l)
                            print *,'O0tau0_vu, tau0_vl =',
     .                                    tau0_vu,tau0_vl
                            endif
                         endif


                        !subtract  contribution from this time step from total 
                        !from previous time step(note that contrib is negative 
                        !here so we are subtracting)

                        enbeam_part(j,i,ie,l ) = enbeam_part(j,i,ie,l )
     .                                        +contrib_en -previous_en
                        enbeam_part_nl(j,i,ie,l ) =
     .                        enbeam_part_nl(j,i,ie,l )+contrib_en_nl
     .                                                -previous_en_nl
                        if(j .lt. 2)then
                   print *,'enbeam_part_nl =',enbeam_part_nl(j,i,ie,l )
                   print *,'previous_en_nl=', previous_en_nl
            print *,'new enbeam_part_nl =',enbeam_part_nl(j,i,ie,l )
                      endif
                        previous_en = contrib_en
                        previous_en_nl = contrib_en_nl
                        wenbeam_part(j,i,ie,l) = wenbeam_part(j,i,ie,l) 
     .                         +contrib_wen -previous_wen
                        previous_wen = contrib_wen
                        Qfi_part(j,i,ie,l) = Qfi_part(j,i,ie,l)
     .                         +contrib_Qfi - previous_Qfi
                        Rfi_part(j,i,ie,l) = Rfi_part(j,i,ie,l)
     .                         +contrib_Rfi - previous_Rfi
                        Rfe_part(j,i,ie,l) = Rfe_part(j,i,ie,l)
     .                         +contrib_Rfe - previous_Rfe
                        previous_Qfi   = contrib_Qfi
                        previous_Rfi   = contrib_Rfi
                        previous_Rfe   = contrib_Rfe
                        if(j .lt. 2)then
                        print *,' '
                        print *,'timenew =',timenew
                        print *,'tau0_vlj(k,j,i,ie,l) =',
     .                      tau0_vlj(k,j,i,ie,l)
                        print *,'tau0_vljoff(k,j,i,ie,l) =',
     .                      tau0_vljoff(k,j,i,ie,l)
                         print *,'enbeam_part(jo,i,ie,l)= ',
     .                           enbeam_part(j,i,ie,l),j,i,ie,l
                         print *,'contrib_en = ',contrib_en
                         print *,'pbeamOff section,k,i =',k,i
                         print *,'contrib_en_nl=',contrib_en_nl
                         print *,'enf_contrib_nl=',enf_contrib_nl
                         print *,'previous_en_nl=',previous_en_nl
                        print *,'pssv(k,j,i,ie,l)  ',pssv(k,j,i,ie,l)
                        if(pssv(k,j,i,ie,l))then
                           print *,'contrib_en_nl=',contrib_en_nl
                           print*,'enbeam_part_nl(j,i,ie,l )=',
     .                                  enbeam_part_nl(j,i,ie,l )
                        endif
                        endif
                     endif      !end  time .ge. beamOff branch
                  if(done) go to 200 !skip the beamOff test since beam must be on first

                  endif   !pssvoff(k,j,i,ie,l) = false branch


c--------------------------beam off section ends here--------------------------------


 
 300             continue
c               print *,'bottom of k loop, k,j,i,ie,l = ',k,j,i,ie,l
               enddo              !end  k loop over times
 200           enbeam_tot(j) = enbeam_tot(j)+enbeam_part(j,i,ie,l)
               enbeam_tot_nl(j) = enbeam_tot_nl(j)+
     .                                     enbeam_part_nl(j,i,ie,l)
               wenbeam_tot(j) = wenbeam_tot(j)+wenbeam_part(j,i,ie,l)
               Qfi_tot(j)     = Qfi_tot(j) + Qfi_part(j,i,ie,l)
               Rfi_tot(j)     = Rfi_tot(j) + Rfi_part(j,i,ie,l)
               Rfe_tot(j)     = Rfe_tot(j) + Rfe_part(j,i,ie,l)
c               print *,'end loop on spatial index j'
               if(j .lt. 2)then
               print *,'enbeam_part(j,i,ie,l) =',
     .                 enbeam_part(j,i,ie,l),j,i,ie,l
               print *,'j=1 only,with i,ie,l =',i,ie,l
               endif
              enddo               !end  j loop over transport grid
              print *,'end loop on beam energies'
            enddo                 !end  ie loop on beam energies
 100          print *,'end loop on sources per beam line'
            enddo                  !end  i loop on sources per beamline
           print *,'end loop on beam lines'
       enddo                      !end  l loop on beamlines

         print *,'on return ', 
     .     'enbeam_part_nl(1,1,1,1),enbeam_part_nl(1,1,1,2)',
     .     ' enbeam_tot_nl(1)= ',
     .     enbeam_part_nl(1,1,1,1),enbeam_part_nl(1,1,1,2),
     .     enbeam_tot_nl(1)
         print *,'***********ending beam_time_dependance***************'
         print *, ' '
         return 
         end

      real*8 function bkef (vpar, tpar, emzpar)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c This routine evaluates the function Ke, which is related to the
c transfer of parallel momentum from fast ions slowing down on electrons.
c The function is defined by Callen et al., IAEA Tokyo, Vol. I, 645 (1974).
c ----------------------------------------------------------------------
c
      external       bkefun
      common /gecom/ vcvo, tstcx, emzrat
c
      vcvo   = vpar
      tstcx  = tpar
      emzrat = emzpar
      a      = 0.0
      b      = 1.0
      ep     = 0.1
      m4     = 3
      bkef   = asimp (a, b, ep, m, m4, bkefun)
      return
c
      end

      real*8 function bkefun (y)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c evaluate argument of integral defining Ke, the fraction of initial
c fast ion parallel momentum collisionally transferred to thermal electrons.
c see comments under function gefun HSJ
c ----------------------------------------------------------------------
c
      common /gecom/ vcvo, tstcx, emzrat
c
      bkefun = 0.0
      if (y .gt. 0.0) then
        v3       = vcvo**3
        arg      = (1.0 + v3)/(y**3+v3)
        alogarg  = LOG (arg)
        pcxlog   = -tstcx*alogarg*0.33333333333334
****    pcx      = (arg)**(-tstcx/3.0)
****    b        = (y**3*arg)**(emzrat/3.0)
****    bkefun   = y**3*pcx*b/(y**3+v3)
        alogy3v3 = LOG (y**3+v3)
        alog3y   = 3.0 * LOG (y)
        blog     = (alog3y+alogarg)*emzrat*0.33333333333334
        bkeflog  = alog3y+pcxlog+blog-alogy3v3
        if (bkeflog .lt. -30.0) then
          bkefun = 0.0
        else
          bkefun = EXP (bkeflog)
        end if
      end if
      return
c
      end

      subroutine bproc (hicm, kb, ke, kj, mb, nj, sb, sbcx, sbion)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --------------------------------------------------------------------------
c subroutine calculates ion and neutral sources due to neutral beam injection
c ---------------------------------------------------------------------------
c
c     inputs:
c          hicm(j,ie,ib,icm) - hot ion creation mode array.  contains
c                              fraction of hot ion birth rate creating
c                              electrons, species 1 and species 2
c                              neutrals.
c          kcm               - maximum number of hot ion creation modes
c          ke                - maximum number of beam energy groups
c          kj                - maximum number of zones
c          mb                - number of beams
c          nj                - number of radial points
c          sb(j,ie,ib)       - source rate of hot ions (#/cm**3-sec)
c
c     outputs:
c          sbcx(j,i)         - hot ion source rate at r(j) due to charge
c                              exchange with species i (#/cm**3-sec)
c          sbion(j)          - hot ion source rate at r(j) producing
c                              electrons (#/cm**3-sec)
c ----------------------------------------------------------------------
c
      dimension hicm(kj,ke,kb,*),sb(kj,ke,*)
      dimension sbcx(kj,2),sbion(*)
c
      do j=1,nj
        sbcx (j,1) = 0.0
        sbcx (j,2) = 0.0
        sbion(j  ) = 0.0
      end do
c
      do    ib=1,mb
        do  ie=1,ke
          do j=1,nj  ! NO CHANGE HERE DUE TO 'DT' BEAMS AND
c                      'D' AND 'T' THERMALS HSJ
            sbcx(j,1) = sbcx(j,1) + hicm(j,ie,ib,2)*sb(j,ie,ib)
            sbcx(j,2) = sbcx(j,2) + hicm(j,ie,ib,3)*sb(j,ie,ib)
            sbion(j)  = sbion(j)
     .               + (1.0-hicm(j,ie,ib,2)-hicm(j,ie,ib,3))*sb(j,ie,ib)
          end do
        end do
      end do
      return
c
      end

      real*8 function ceef (i, j, te)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- e-impact excitation rate coefficient from Vriens & Smeets
c     used only for excitations to final state with n > 2 (i.e. j > 3).
c
      parameter  (ms = 21, mc = 35)
      common /b1/ kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      common /b2/ nouthx, istart, ihxbug
      common /b4/ en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc),
     .            de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
c
      data        ryd/13.6/
c
      ceef1 (ni, nj, te) = 1.6e-7 * SQRT (te)
     . * EXP (-(ryd / te) * (1.0 / FLOAT (ni)**2 - 1.0 / FLOAT (nj)**2))
     . * (ae(ni,nj) * LOG (0.3 * te / ryd + de2(ni,nj)) + be(ni,nj))
     . / (te + ge2(ni,nj) * LOG (1.0 + (FLOAT (ni)**3) * te / ryd))
c
      if (i .ge. j .or. j .le. 3) then
        ihxbug = 8
        if (nouthx .gt. 0) then
          write (nouthx, 10)  i, j
   10     format (/ ' ERROR in CEEF: i = ', i3, '    j = ', i3 /)
        end if
        ceef = 0.0
      else
        ni   = nfhx(i)
        nj   = j - 1
        ceef = ceef1 (ni, nj, te)
        id   = 2
        if (ni .eq. 1)  id = 1
      end if
c
      return
c
      end

      real*8 function cief (i, te)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- e-impact ionization rate coefficient from Vriens & Smeets,
c     Phys. Rev. a 22, 940 (1980).
c     used only for 2p (i = 3) or for n > 2 (i=n+1).
c
      parameter  (ms = 21, mc = 35)
      common /b1/ kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      common /b2/ nouthx,istart,ihxbug
      common /b4/ en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc),
     .            de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
c
      ent(i) = en(i) / te
c
      if (i .le. 2) then
        ihxbug = 7
        if (nouthx .gt. 0) then
          write  (nouthx, 10) i
   10     format (/ ' ERROR in CIEF: i = ', i3 /)
        end if
        cief = 0.0
      else
        cief = 9.56e-6 * EXP (-ent(i))
     .                 / (te * SQRT (te) * (ent(i)**2.33
     .                 + 4.38 * ent(i)**1.72 + 1.32 * ent(i)))
      end if
c
      return
c
      end
c




      real*8 function cxevala (vf)
c
      implicit  none
      real *8 
     .     vf,vfsq,vfcb,vcrit3,erel,enrgy_factor,sigcx,cxrat,vc_l3
      common /cx_bblock/ vc_l3,enrgy_factor
c
c ----------------------------------------------------------------------
c --- evaluates the function which is integrated by cxint
c     cxeval = (v**2/(v**3+vc**3))(sigcx*v)
c --- where sigcx*v is the charge exchange rate in cm**3/sec
c --- sigcx*v is evaluated by function cxr.
c ----------------------------------------------------------------------
c
c
      vfsq =vf*vf
      vfcb = vfsq*vf                              !( v in cm/sec)
      erel   = enrgy_factor*vfsq                  ! in keV (assumes vrel = vbeam)
      sigcx  = cxrat(erel)                         ! so debug can look at it
      cxevala = (vfsq/(vfcb+vc_l3))*sigcx
      return
c
      end







      real*8 function cxinta (vlow_lim,vup_lim,taus_l,enn_l)
c
      implicit  none
      real *8 vlow_lim,vup_lim,answ,taus_l,enn_l
c
c ----------------------------------------------------------------------
c --- function calculates the exponential integrating factor due to
c --- loss of fast ions by charge exchange.
c ----------------------------------------------------------------------
c
c
      external cxevala
c
      call qgaus2a(cxevala,vlow_lim,vup_lim,answ)
      answ  = answ * enn_l                ! enn_l is neutral density, #/cm**3
      answ  = taus_l * answ
      cxinta = DEXP (-answ)
      return
c
      end

      real*8 function cxrat (x)
c
      implicit  none
      real *8 x,dum,tc
c
c ----------------------------------------------------------------------
c this function calculates the charge exchange rate for hydrogen atoms
c interacting with protons in units of cm**3/s.
c x is in units of keV for 1.0e-3 .le. x .le. 100, the
c the formula is taken from the paper by r.l. freeman and e.m. jones
c clm-r 137 culham laboratory 1974.
c for x .lt. 1.0e-3, a rate coefficient derived from an analytic average
c over the approximate cross section  sigma = 0.6937e-14*(1.0-0.155*LOG10
c (e/1ev))**2 is used.  this cross section is an approximation to that
c given by riviere, nuclear fusion 11,363(1971).
c ----------------------------------------------------------------------
c
      if (x .lt.   1.0e-3)  go to 10
      if (x .gt. 100     )  go to 20
      tc  = DLOG (x)+6.9077553D0
      dum = 0.2530205D-4-tc*0.8230751D-6
      tc  = -0.1841757D2+tc*(0.528295D0-tc*
     .      (0.2200477D0-tc*(0.9750192D-1-tc*
     .      (0.1749183D-1-tc*(0.4954298D-3+tc*(0.2174910D-3-tc*dum))))))
      tc  = DEXP (tc)
      cxrat = tc
      return
c
   10 tc  = 0.50654D0 - 6.7316D-2 * DLOG  (x)
      cxrat =           3.3340D-7 * DSQRT (x) * (tc*tc + 7.454D-3)
      return
c
   20 cxrat = 0.D0
      return
c
      end



      subroutine d01bbf (type, inum, weight, abscis, ier)
c
      USE d01bff_dat

      implicit  integer (i-n), real*8 (a-h, o-z)
c
      character*2 type
      dimension   weight(inum), abscis(inum)
      common /b2/ nouthx, istart, ihxbug
      dimension   num (4)
      data        num  /10, 16, 20, 24/
      dimension   iorg(4)
      data        iorg / 1, 11, 27, 47/
c
c      common / cgh / wgh(70), xgh(70)
c

c

c
c      common / cgl / wgl(70), xgl(70)
c

c


      do i=1,4
        index = i
        if (inum .eq. num(i))  go to 60
      end do
      ier = 1
      return
c
   60 ier = 0
c
      if (type .eq. 'gl') then
        do i=1,inum
          irel      = iorg(index) - 1 + i
          weight(i) = wgl(irel)
          abscis(i) = xgl(irel)
        end do
        go to 400
      end if
c
      if (type .eq. 'gh') then
        do i=1,inum
          irel      = iorg(index) - 1 + i
          weight(i) = wgh(irel)
          abscis(i) = xgh(irel)
        end do
        go to 400
      end if
c
      ier = 1
  400 return
c
      end

      subroutine decomp (ndim, n, a, cond, ipvt, work)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     decomposes a real matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     use solveq to compute solutions to linear systems.
c
c     input..
c
c        ndim = declared row dimension of the array containing  a.
c
c        n = order of the matrix.
c
c        a = matrix to be triangularized.
c
c     output..
c
c        a  contains an upper triangular matrix  u  and a permuted
c          version of a lower triangular matrix  i-l  so that
c          (permutation matrix)*a = l*u .
c
c        cond = an estimate of the condition of  a .
c           for the linear system  a*x = b, changes in  a  and  b
c           may cause changes  cond  times as large in  x .
c           if  cond+1.0 .eq. cond , a is singular to working
c           precision.  cond  is set to  1.0e+32  if exact
c           singularity is detected.
c
c        ipvt = the pivot vector.
c           ipvt(k) = the index of the k-th pivot row
c           ipvt(n) = (-1)**(number of interchanges)
c
c     work space..  the vector  work  must be declared and included
c                   in the call.  its input contents are ignored.
c                   its output contents are usually unimportant.
c
c     the determinant of a can be obtained on output by
c        det(a) = ipvt(n) * a(1,1) * a(2,2) * ... * a(n,n).
c
      integer  ipvt(n), ndim, n
      integer  nm1, i, j, k, kp1, kb, km1, m
      real*8   a(ndim,n), cond, work(n)
      real*8   ek, t, anorm, ynorm, znorm
c
      ipvt(n) = 1
      if (n .eq. 1)  go to 80
      nm1 = n - 1
c
c     compute 1-norm of a
c
      anorm = 0.0
      do j=1,n
        t = 0.0
        do i=1,n
          t = t + ABS (a(i,j))
        end do
        if (t .gt. anorm)  anorm = t
      end do
c
c     gaussian elimination with partial pivoting
c
      do 35 k=1,nm1
         kp1 = k+1
c
c        find pivot
c
         m = k
         do i=kp1,n
           if (ABS (a(i,k)) .gt. ABS (a(m,k))) m = i
         end do
         ipvt(k) = m
         if (m .ne. k) ipvt(n) = -ipvt(n)
         t = a(m,k)
         a(m,k) = a(k,k)
         a(k,k) = t
c
c        skip step if pivot is zero
c
         if (t .eq. 0.0)  go to 35
c
c        compute multipliers
c
         do i=kp1,n
           a(i,k) = -a(i,k)/t
         end do
c
c        interchange and eliminate by columns
c
         do 30 j=kp1,n
           t = a(m,j)
           a(m,j) = a(k,j)
           a(k,j) = t
           if (t .eq. 0.0)  go to 30
           do i=kp1,n
             a(i,j) = a(i,j) + a(i,k)*t
           end do
   30    continue
   35 continue
c
c     cond = (1-norm of a)*(an estimate of 1-norm of a-inverse)
c     estimate obtained by one step of inverse iteration for the
c     small singular vector.  this involves solving two systems
c     of equations, (a-transpose)*y = e  and  a*z = y  where  e
c     is a vector of +1 or -1 chosen to cause growth in y.
c     estimate = (1-norm of z)/(1-norm of y)
c
c     solve (a-transpose)*y = e
c
      do k=1,n
        t = 0.0
        if (k .eq. 1)  go to 45
        km1 = k-1
        do i=1,km1
          t = t + a(i,k)*work(i)
        end do
   45   ek = 1.0
        if (     t .lt. 0.0)  ek = -1.0
        if (a(k,k) .eq. 0.0)  go to 90
        work(k) = -(ek + t)/a(k,k)
      end do
c
      do kb=1,nm1
        k   = n - kb
        t   = 0.0
        kp1 = k+1
        do i=kp1,n
          t = t + a(i,k) * work(k)
        end do
        work(k) = t
        m       = ipvt(k)
        if (m .ne. k) then
          t       = work(m)
          work(m) = work(k)
          work(k) = t
        end if
      end do
c
      ynorm = 0.0
c
      do i=1,n
         ynorm = ynorm + ABS (work(i))
      end do
c
c     solve a*z = y
c
      call solveq(ndim, n, a, work, ipvt)
c
      znorm = 0.0
c
      do i=1,n
        znorm = znorm + ABS (work(i))
      end do
c
c     estimate condition
c
      cond = anorm*znorm/ynorm
      if (cond .lt. 1.0)  cond = 1.0
      return
c
c     1-by-1
c
   80 cond = 1.0
      if (a(1,1) .ne. 0.0)  return
c
c     exact singularity
c
   90 cond = 1.0e+32
      return
c
      end

      real*8 function dfhx (beta)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- approximation to the function defined by Janev & Presnyakov,
c     j. phys. b 13, 4233 (1980).
c
      dimension dd(38)
      data      dd/
     . .10450, .121, .138, .157, .175, .200, .229, .260, .300, .339,
     .   .367, .389, .402, .410, .398, .376, .346, .317, .285, .255,
     .   .227, .205, .185, .168, .152, .138, .124, .110, .099, .089,
     .   .079, .070, .062, .054, .047, .041, .035, .02898/
c
      beta1 = 1.0 / beta
      if (beta1 .le. 0.2  )  go to 110
      if (beta1 .ge. 1.0e3)  go to 120
      a    = 10.0 * LOG10 (beta1) + 8.0
      ia   = MIN0 (37, INT (a))
      dfhx = dd(ia)+(a-FLOAT (ia))*(dd(ia+1)-dd(ia))
      return
c
  110 dfhx = 0.5 * beta * (1.0 - 1.0/(8.0 * beta * SQRT (beta)))
     .                  * EXP (-SQRT (2.0 * beta))
      return
c
  120 dfhx = 4.0 * beta * LOG (1.4 * beta1)
      return
c
      end

      subroutine eigen (xeig)
      USE cpub_dat
      USE replace_imsl,                    ONLY : my_eigrf
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      parameter    (ms = 21, mc = 35, mz = 1, mi = 2)
c
c      common /cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8,
c     .              cpuii, cpuiz, cpuplf
      common /b0  / er0, v0, te, ti, ami(mi), deni(mi), amz(mz),
     .              denz(mz), zcor(mz), zsqcor(mz), dene,
     .              ns, nc, numi, numz, iz(mz), izstrp(mz)
      common /b2  / nouthx,istart,ihxbug
      common /b7  / q(ms+1,ms+1)
      common /b9  / eigvr(ms+1,ms+1),eigvl(ms+1),cj(ms+1)
      common /b10 / xfrac
c
      dimension qq(ms+1,ms+1), work2(3*(ms+1))
      dimension w(2*(ms+1)),z(2*(ms+1),ms+1)
      dimension cjmat(ms+1,ms+1), workdc(ms+1), ipvt(ms+1)
      real*8    infac(ms+1,ms+1), lexp(ms+1), in(ms+1)
c
      data      xtest /1.0/
c
c     on return from eigrf, w contains the complex eigenvalues
c     as w = [wr(1),wi(1), ... ,wr(ns),wi(ns)]
c     on return from eigrf, z contains the complex eigenvectors
c     as z(col1) =  [zr(1),zi(1), ... ,zr(ns),zi(ns)], etc.,
c     each column representing an eigenvector
c
      call SECOND (cpua)
c
      do 10 i=1,ns
      do 10 j=1,ns
        qq(i,j) = q(i,j)
   10 continue
      if (nouthx .gt. 0) then
        write  (nouthx, 11)
   11   format ( 'q _ matrix')
        do i=1,ns
          write  (nouthx, 12)  (q(i,j), j=1,ns)
   12     format (1x, 1p7e11.3)
        end do
      end if
c
****  ifail1 = 0
****  call f02aff (qq,ms+1,ns,eigr,eigi,intger,ifail1)
****  if (ifail1 .eq. 0)  go to 15
c
      ijob = 1

      call my_eigrf (qq, ns, ms+1, ijob, w, z, ms+1, work2, ifail1) ! 888889999
      if (ifail1 .eq. 0)  go to 15
      ihxbug = 9 
      if (nouthx .gt. 0) then
        write  (nouthx, 3939) ifail1
 3939   format (' ERROR in EIGRF: ifail1 = ', i4)
      end if
      return
c
   15 do j=1,ns
        do i=1,ns
          eigvr(i,j) = z(2*i-1,j)
        end do
        eigvl(j) = w(j*2-1)
      end do
c
      eigmin = 1.0e30
      do i=1,2*ns,2
        eigmin = MIN (eigmin, ABS (w(i)))
      end do
c
****  do 20 i=1,ns
****    eigmin = MIN (eigmin, ABS (eigr(i)))
***20 continue
c
      xeig = 0.0
      if (eigmin .ne. 0.0)  xeig = v0 / eigmin
c
c     determine fraction of 3rd excited state (approximate)
c
      sum = 0
      do i=1,ns
        sum = sum + eigvr(i,1)
      end do
c
      xfrac = eigvr(3,1) / sum
      if (nouthx .eq. 0)  return
c
      write  (nouthx, 1000)  xeig
 1000 format (' length (eigenvalue):  ', 1pe11.4)
      write  (nouthx, 1010)  (i, eigvl(i), i=1,ns)
 1010 format ('    i= ', i2, '         eigenvalue= ', 1pe11.4)
      do j=1,ns
        write  (nouthx, 1020)  j, (eigvr(i,j), i=1,ns)
 1020   format (' eigenvector[', i2, ']= ' / (10x, 1p5e11.4))
      end do
c
c     apply initial condition constraint
c
      do 1200 i=1,ns
      do 1200 j=1,ns
        cjmat(i,j) = eigvr(i,j)
 1200 continue
c
      cj(1) = 1.0
      do i=2,ns
        cj(i) = 0.0
      end do
c
      call decomp (ms+1, ns, cjmat, cond, ipvt, workdc)
      call solveq (ms+1, ns, cjmat,   cj, ipvt)
c
      write  (nouthx, 1320)  (cj(i), i=1,ns)
 1320 format (' bc constants: ' / (10x,1p5e11.4))
c
      do j=1,ns
        do i=1,ns
          lexp(i)    = eigvl(i) / v0
          infac(j,i) = cj(i) * eigvr(j,i)
        end do
        write  (nouthx, 1400)  j, (infac(j,i), lexp(i), i=1, ns)
 1400   format ( ' *** i(',i2,') =    ' /
     .         (12x, 1pe15.7, ' * EXP (', 1pe15.7, ' * x)') )
      end do
c
      do j=1,ns
        in(j) = 0.0
        do i=1,ns
          in(j) = in(j) + infac(j,i) * EXP (lexp(i)*xtest)
        end do
      end do
c
      write  (nouthx, 1450)
 1450 format (/// ' qnj * ij / v0 ...........' /)
c
      do j=1,ns
        sum = 0.0
        do i=1,ns
          sum = sum + q(j,i) * in(i)
        end do
        write  (nouthx, 1500)  j, sum/v0
 1500   format (15x,i3,5x,1pe15.7)
      end do
c
      call SECOND (cpub)
      cpu7 = cpu7 + cpub - cpua
c
      return
c
      end

      real*8 function encapf (vcvo, tstcx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  This routine evaluates the function N, which is related to the rate
c     at which fast ions slow down on electrons and thermal ions.  This
c     function is defined by Callen et al., IAEA Tokyo, Vol. I, 645 (1974).
c     It is assumed that taus/taucx is independent of the ion speed.
c     This allows calculation of Pcx,the probability against charge
c     exchange, as
c           Pcx = [(v0**3+vc**3)/(v**3+vc**3)]**(-taus/(3.0*taucx))
c     Subsequently the integral defining N can be evaluated with the
c     result that
c         N = (taucx/taus)*(1.0-EXP (-tauf/taus))
c
c --- input
c  tstcx            ratio of taus/taucx
c                   for fast alpha particle this ratio is taken as 0.0
c
c  vcvo             SQRT (Ec/E0)
c
c ------------------------------------------------------------------ HSJ
c
      v3   = vcvo**3
      tfts = LOG (1.0 + 1.0/v3) / 3.0
      if (tstcx .le. 0.01) then
        encapf = tfts
      else
        encapf = (1.0 - EXP (-tfts*tstcx)) / tstcx
      end if
      return
c
      end


      subroutine eval_beam_integrals(tau0_vu,tau0_vl,Sbeam0,
     .                         no_fast_losses,taus_l,vc_l,vb_l,
     .                         eb_l,enn_l,mass_beam_l,emzrat_l,
     .                         Pf0,enf_contrib,
     .                         enf_contrib_nl,wenf_contrib,
     .                         Qfi_contrib,Rfi_contrib,Rfe_contrib)
c-------------------------------------------------------------------HSJ-2/23/00
c      This subroutine calculates quantities determined by
c      quadrature from the analytic fast ion distribution function.
c      1) the fast ion density:                           enf_contrib
c      2) fast ion energy:                                wenf_contrib
c      2) depostion of energy to electron:                
c      3)                        ions
c      4) depostion of parallel momentum 
c      5) torque density for toroidal rotation 
c      4)

c---INPUT
c     enn_l         local neutral density 1/cm**3
c     mass_beam_l   in grams
c     eb_l          beam energy,kev
c     emzrat_l      mi<Z>/(mf[Z])
c     Sbeam0        Beam intensity, #/cm**3/sec
c---OUTPUT
c
c    enf_contrib        fast ion density 1/cm**3
c    Qfi_contrib        kev/(cm**3 sec) energy deposited on ions
c    Qfe                kev/(cm**3 sec) energy deposited on electrons
c                       Qfe is determined by  wenf. Hence it is not 
c                       necessary to calculate it in this subroutine
c    wenf_contrib       kev/cm**3 stored fast ion energy density
c                       toroidal momentum deposited on ions
c    Rfi_contrib        momentum rate density (g cm/sec)/(cm**3 sec)     ions
c    Rfe_contrib                                                         electrons
c
c--------------------------------------------------------------------------------

      implicit none
      real *8
     .     enf_contrib,tau0_vu,tau0_vl, Sbeam0,wenf_contrib,
     .     gaussint,taus_l,vc_l,vb_l,enf_int,wenf_int,v2,enn_l,
     .     vc_l3,mass_beam_l,enrgy_factor,enf_int_nl,Qfi_contrib,
     .     eb_l,enf_contrib_nl,Qfi_int,Rfe_int,Rfi_int,Rfe_contrib,
     .     Rfi_contrib,emzrat_l,Pf0
      logical no_fast_losses
      common /cx_bblock/ vc_l3,enrgy_factor



      vc_l3 = vc_l*vc_l*vc_l
      v2=(vc_l/vb_l)**2
      enrgy_factor = 3.1207D8 * mass_beam_l  ! converts 0.5*mass_beam*vb**2 to kev

      call beam_int(tau0_vu,tau0_vl, vb_l,vc_l,taus_l,enn_l,emzrat_l,
     .             enf_int,enf_int_nl,wenf_int,Qfi_int,Rfe_int,Rfi_int,
     .                                               no_fast_losses)


      enf_contrib    = Sbeam0*enf_int                             !1/(cm**3) with losses
      enf_contrib_nl = Sbeam0*enf_int_nl                          !1/(cm**3) no losses
      wenf_contrib   = Sbeam0*eb_l*v2*wenf_int                    !kev/cm**3
      Qfi_contrib    = Sbeam0*eb_l*v2*vc_l*2.d0*Qfi_int/taus_l    !kev/(cm**3 sec)
      Rfe_contrib    = Sbeam0*Pf0*Rfe_int/(taus_l*vb_l)           !g/(cm**2 sec**2)
      Rfi_contrib    = Sbeam0*Pf0*Rfi_int*vc_l3*
     .                           (1.d0+emzrat_l)/(taus_l*vb_l)    !g/(cm**2 sec**2)

      return
      end


      real *8 function eval_vthcut3(j)
      USE param
      ! create a single routine where vthcut3 is evaluated:

      USE nub4
         implicit none
         real *8 vthcut3
         integer*4 j
c         include 'param.i'
c         include 'nub4.i'
      
         if(beam_thermal_cutoff == -1)then
            vthcut3 = beam_thermal_speed**3
         else if( beam_thermal_cutoff == 0)then
            vthcut3 = (therm_frac*vthi(j))**3
         else if(beam_thermal_cutoff == 1)then
            vthcut3 = vthi(j)**3
         endif         
         eval_vthcut3 = vthcut3
         
      return
      end




      subroutine freya (mb, atw, codeid, elong, mi, mj, nion, psi, r,
     .                  rin, rmax, z, zax, zmin, zmax, bpow, eb, iexcit,
     .                  enbeams,time,time0)
c
c
c ----------------------------------------------------------------------
c
c  this subroutine calculates the particle and energy sources
c     due to neutral beam injection.
c
c  the input quantities are:
c
c  mb              Number of neutral beam injectors 
c  anglev(ib)      Vertical angle (degrees) between optical axis
c                    and horizontal plane; a positive value indicates
c                    particles move upward
c  angleh(ib)      Horizontal angle (degrees) between optical axis and
c                    vertical plane passing through pivot point and
c                    toroidal axis; a zero value denotes perpendicular
c                    injection, while a positive value indicates par-
c                    ticles move in the co-current direction
c  nsourc          Number of sources per beamline.
c                    If 1, source is centered on beamline axis.
c                    If nsourc = 2, distinguish between the beamline
c                     axis and the source centerline (optical axis).
c                     The two sources are assumed to be mirror images
c                     through the beamline axis.
c                    In either case, the exit grid plane is perpendicula
c                    to the beamline axis, and contains the source
c                    exit grid center(s).
c                    If nsourc = 2, the alignment of the sources w.r.t.
c                    the beamline axis is specified through bhofset,
c                    bvofset, and bleni (described further below).
c  bvofset(ib)     Vertical offset from beamline axis to center
c                    of each source (cm; used only for nsourc = 2)
c  bhofset(ib)     Horizontal offset from beamline axis to center
c                    of each source (cm; used only for nsourc = 2)
c  bleni(ib)       Length along source centerline (source optical axis)
c                    source to intersection point with the beamline axis
c  sfrac1(ib)      Fraction of source current per beamline coming
c                    from upper source (used only for nsourc = 2)
c  bcur(ib)        Total current (a) in ion beam (used only if bptor
c                    is zero)
c  bptor(ib)       Total power (w) through aperture into torus; when
c                    nonzero, bptor takes precedence over bcur
c  bshape(ib)      Beam shape
c                    'circ': circular
c                    'rect': rectangular
c  bheigh(ib)      Height of source (cm)
c  bwidth(ib)      Width of source (cm); diameter for
c                     circular source.
c  bhfoc(ib)       Horizontal focal length of source (cm)
c  bvfoc(ib)       Vertical focal length of source (cm)
c  bhdiv(ib)       Horizontal divergence of source (degrees)
c  bvdiv(ib)       Vertical divergence of source (degrees)
c  b1ins(i)        ratio of R*B to Rax*Bax along horizontal
c                    chord inside and through the magnetic axis vs.
c                    a uniform mesh (rinsid(i)) in major radius;
c                    b1ins is needed only if iborb .gt. 0
c  b1ots(i)        ratio of R*B to Rax*Bax along horizontal
c                    chord outside and through the magnetic axis vs.
c                    a uniform mesh (rotsid(i)) in major radius;
c                    b1ots is needed only if iborb = 1
c  b2ins(i)        ratio of B to Btor along horizontal
c                    chord inside and through the magnetic axis vs.
c                    a uniform mesh (rinsid(i)) in major radius;
c                    b2ins is needed only if iborb = 1
c  b2ots(i)        ratio of B to Btor along horizontal
c                    chord outside and through the magnetic axis vs.
c                    a uniform mesh (rotsid(i)) in major radius;
c                    b2ots is needed only if iborb = 1
c  ebkev(ib)       Maximum particle energy in source (keV)
c  enbeams         beam density,passed to prep_mcgo
c  fbcur(ie,ib)    Fraction of current at energy ebkeV/ie
c  iborb           Flag for modeling orbit effects on beam-
c                    generated fast ions
c                    3, use mcgo for fast ion slowing down,initial
c                       conditions for the fast ions are generated
c                       using freya.
c                    2, model prompt loss with STAMBAUGH model
c                    1, model orbit effects
c                    0, do not model orbit effects
c  npart           Number of particles followed into plasma
c                    (suggest 10000)
c  npskip          Ratio of number of particles followed into plasma
c                    to number of source particles (suggest 1)
c  naptr           Total number of apertures encountered by a particle
c                    as is moves from the source into the plasma chamber
c                    Maximum is specified by parameter nap ( = 4).
c                    First set of apertures encountered by the particle
c                    are assumed centered on the source axis, and subseq
c                    apertures are centered on the beamline axis; the
c                    distinction is made through ashape.
c  ashape(iap,ib)  Aperture shape.
c                   Prefix 's-' indicates source axis centered.
c                   Prefix 'b-' indicates beamline axis centered.
c                     's-circ'          'b-circ'
c                     's-rect'          'b-rect'
c                     's-vert'          'b-vert'
c                     's-horiz'         'b-horiz'
c                                       'b-d3d'
c                   (circ = circular aperture, rect=rectagular,
c                    vert = limits vertical height of source particles,
c                    horiz = limits horizontal height of source particles,
c                    d3d= special DIII-D polygonal aperture)
c                    'circ': circular
c                    'rect': rectangular
c  aheigh(iap,ib)  Height of aperture (cm)
c  awidth(iap,ib)  Width  of aperture (cm); diameter for circular aperture
c  alen(iap,ib)    Length from source to aperture for 's-type' aperatur
c                    and from exit grid plane along beamline axis for
c                    'b-type' apertures.
c  blenp(ib)       Distance along beamline axis from source exit
c                    plane to the fiducial "pivot" point.
c  rpivot(ib)      Radial position of pivot (cm)
c  zpivot(ib)      Axial position of pivot (cm)
c     atw(k)          atomic mass of ion species k
c     codeid          flux surface geometry
c                       "onedee"     : elliptical
c                       anything else: nonelliptical
c     elong           elongation (height/width) of elliptical
c                       cross-section plasma
c     ibion           ion index of beam species
c     mf              number of flux zones
c     mi              number of radial mesh points for nonelliptical
c                       plasma
c     mj              number of axial mesh points for nonelliptical
c                       plasma
c     nion            number of ion species
c     pinsid(i)       poloidal magnetic flux (G-cm2) along horizontal
c                       chord inside and through the magnetic axis vs.
c                       a uniform mesh (rinsid(i)) in major radius;
c                       pinsid is needed only if iborb .gt. 0
c     potsid(i)       poloidal magnetic flux (G-cm2) along horizontal
c                       chord outside and through the magnetic axis vs.
c                       a uniform mesh (rotsid(i)) in major radius;
c                       potsid(1) and potsid(mf) are needed if
c                       codeid .ne. 'onedee'; the entire potsid array is
c                       needed if iborb = 1
c     psi(i,j)        poloidal magnetic flux (G-cm2) at mesh point i,j
c                       (needed only if codid .ne. 'onedee')
c     csgn            sign of the scalar product of plasma current and
c                       toroidal unit vector phi, j dot phi, for the
c                       right-handed triad (R,phi,Z) where Z points
c                       upward in the lab (DIII-D) frame.
c     psivol(i)       volume (cm3) of flux zone i; depending upon
c                       whether codeid = 'onedee' or not, psivol is
c                       chosen such that either r or SQRT (psi-psiax)
c                       varies a constant amount from one flux surface
c                       to the next
c     rinsid(i)       major radius (cm) of uniform mesh along horizontal
c                       chord inside and through the magnetic axis;
c                       rinsid is needed only if iborb = 1
c     rotsid(i)       major radius (cm) of uniform mesh along horizontal
c                       chord outside and through the magnetic axis;
c                       rotsid(1) and rotsid(mf) are needed if
c                       codeid = 'onedee'; the entire rotsid array is
c                       needed if iborb = 1
c     r(i)            radial mesh point i for nonelliptical plasma (cm)
c     rin             major radius of inside of vacuum vessel (cm)
c                     for 2D equilibria this number must coincide with
c                     rmhdgrid(1).
c     rmax            maximum radial position of plasma (cm)
c     rpivot(ib)      radial position of pivot point (cm)
c     sfrac1(ib)      fraction of source current per beamline coming
c                       from upper source (used only for nsourc = 2)
c     sofset(ib)      vertical offset from optical axis to center
c                       of each source (cm; used only for nsourc = 2)
c     z(j)            axial mesh point j for nonelliptical plasma (cm)
c     zax             axial position of magnetic axis (cm) for
c                       elliptical plasma
c     zmin            minimum axial position of plasma (cm)
c     zmax            maximum axial position of plasma (cm)
c     zpivot(ib)      axial position of pivot point (cm)
c     zne(i)          electron density in zone i (cm-3)
c     zni(i,k)        density of ion species k in zone i (cm-3)
c     zte(i)          electron temperature in zone i (keV)
c     zzi(i,k)        charge number of ion species k in zone i
c     iexcit          hexnb switch
c                     0, (default) do not account for atomic excitation
c                        in neutral stopping calculation
c                    >0, use hexnb to account for atomic excitation
c     inubpat        switch controlling 2-d beam deposition calculation
c                    0, (default) bypass 2-d deposition
c                    1, determine beam deposition on (r,z) grid, defined
c                       by npat, and output neutral beam deposition
c                       information to file 'beamdep'
c     npat           dimensions of (r,z) grid for neutral beam
c                    deposition calculation.  used if inubpat = 1.0
c                    npat(1) = number of new 'r' elements
c                    npat(2) = number of new 'z' elements
c                    default:  npat(1) = mi, npat(2) = mj
c
c  the output quantities are:
c
c     bion(ie,ib)     intensity of ion beam (particles/s)
c     bneut(ie,ib)    intensity of neutral beam (particles/s)
c     bpow(ie,ib)     beam power to aperture (w)
c     eb(ie,ib)       particle energy (keV)
c     fap(ie,ib)      fraction of beam stopped by aperture
c     fwall(ie,ib)    fraction of beam incident on wall (shinethrough)
c     forb(ie,ib)     fraction of beam lost on orbits
c     ftrapfi(i,ie,ib)fraction of trapped fast ions in each zone
c     fb11(ie,ib)     fraction of ions passing and axis-encircling
c     fb10(ie,ib)     fraction of ions passing and not encircling
c     fb01(ie,ib)     fraction of ions trapped and axis-encircling
c     fb00(ie,ib)     fraction of ions trapped and not encircling
c     fber(ie,ib)     fraction of ions trapped for which error was detected
c     hibrz(i,ie,ib)  normalized hot ion birth rate
c     hdepz(i,ie,ib)  normalized hot ion deposition rate
c     hicmz(i,ie,ib,imd)
c                     hot ion creation mode (i.e. ionization or cx)
c     npts            number of birth points to be plotted
c     xpts(ii)        x coordinate of birth point
c     ypts(ii)        y coordinate of birth point
c     zpts(ii)        z coordinate of birth point
c     vx(ii)          x component of birth velocity
c     vy(ii)          y component of birth velocity
c     vz(ii)          z component of birth velocity
c     wb11(ie,ib)     orbit width of fb11 ions (cm)
c     wb10(ie,ib)     orbit width of fb10 ions (cm)
c     wb01(ie,ib)     orbit width of fb01 ions (cm)
c     wb00(ie,ib)     orbit width of fb00 ions (cm)
c     zetaz(i,ie,ib)  average pitch angle cosine of deposited hot ions
c     angmpz(i,ie,ib) average toroidal angular momentum of a single ion
c                     born in zone i.  no orbit smearing is done.  for a
c                     pencil beam, angmpz will be constant across all
c                     zones.  the orbit smearing of the toroidal
c                     momentum is accomplished (in postnub) by
c                     multiplying angmpf by the smeared deposition rate
c                     of fast ions.
c olossc(i,ie,ib) the cold return current flowing through each   ! rs
c                 flux zone ro neutralize the orbit loss current ! rs
c                 (A/cm**3)                                      ! rs
c                 assume routine olossa (activated by iborb = 2) ! rs
c                 has supplied tmin(nchi,izone) and chi(nchi)    ! rs
c                 tmin(ev) is the energy of the loss boundary    ! rs
c                 chi(degrees) is array of pitch angles (180-90) ! rs
c                 nchi is the dimension of chi                   ! rs
c                 izone is the label of the flux zone            ! rs
c                 izone = 1 is the first zone at magnetic axis   ! rs
c                 izone = mfm1(mf-1) is zone inside the boundary ! rs
c ----------------------------------------------------------------------
c
      USE param
      USE io 
      USE mhdpar
      USE nub   ! hicmz,etc
      USE nub2
      USE verbose
      USE nub4
      USE constnts
      USE mcgo
      implicit  integer (i-n), real*8 (a-h, o-z)
D     include 'mpif.h'

c
      character*8 codeid
      character *256 spawn_command
      parameter  (nxx2 = nw*2, nyx2 = nh*2)
      dimension   atw(kion),psi(nw,nh),r(nw),z(nh)
      dimension   bpow(ke,kb),eb(ke,kb)
      dimension   cangv(kb),cangh(kb),sangv(kb),sangh(kb),thetp(kb),
     .            thetpp(kb),costp(kb),sintp(kb),costpp(kb),sintpp(kb),
     .            iatype(nap,kb)
c     .            vbeam(ke,kb)    vbeam is now part of nub4.i HSJ
      dimension   sgxn(kcmp1,kz,kbe,ksge),sgxnloc(kbe),sgxnmi(ke,kb),
     .            hxfrac(ke,kb),enbeams(*)
      dimension  mlost1(ke,kb),mlost2(ke,kb)
      dimension   fin(kf),fot(kf),zin(kf),zot(kf)
      dimension   wt(kz),zetorb(kz),nmbrz(kz)
      dimension   idebug(5)
      dimension   rzpat(nxx2,nyx2,ke,kb)

D      real *8 ,DIMENSION(:,:),ALLOCATABLE::   t_fap,  t_fwall
D      real *8 ,DIMENSION(:,:),ALLOCATABLE::   t_forb, t_b11
D      real *8 ,DIMENSION(:,:),ALLOCATABLE::   t_fb10, t_fb01
D      real *8 ,DIMENSION(:,:),ALLOCATABLE::   t_fb00, t_fber
D      real *8 ,DIMENSION(:,:),ALLOCATABLE::   t_wb11, t_wb10
D      real *8 ,DIMENSION(:,:),ALLOCATABLE::   t_wb01, t_wb00
D      real *8 ,DIMENSION(:,:),ALLOCATABLE::   t_fb11
D      real *8 ,DIMENSION(:,:,:),ALLOCATABLE:: t_ftrapfi, t_hibrz
D      real *8 ,DIMENSION(:,:,:),ALLOCATABLE:: t_hdepz, t_angmpz
D      real *8 ,DIMENSION(:,:,:),ALLOCATABLE:: t_zetaz, t_olossc
D      real *8 ,DIMENSION(:,:,:,:),ALLOCATABLE:: t_hicmz
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_nmbrz
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_xpts, t_ypts, t_zpts
D      real *8 ,DIMENSION(:),ALLOCATABLE:: t_vx, t_vy, t_vz
      data pio180 /0.017453293/



c
      external LENGTH
c
c ----------------------------------------------------------------------
c general FREYA initialization
c ----------------------------------------------------------------------
c
       numprocs =1
       myid = 0 
       master = 0

D      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr) !get processor id
D      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr) !get num processors
D      print *,"freya, myid = ", myid
D      print *,"freya, numprocs = ",numprocs
D      if(numprocs .gt. 1)npskip =1 

D      allocate(t_hibrz(kz,ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_hibrz",myid,istat) 
D      t_hibrz(:,:,:) =0.0 

D      allocate(t_hdepz(kz,ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_hdepz",myid,istat)
D      t_hdepz(:,:,:) =0.0 

D      allocate(t_angmpz(kz,ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_angmpz",myid,istat)
D      t_angmpz(:,:,:) = 0.0 

D      allocate(t_ftrapfi(kz,ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_ftrapfi",myid,istat)
D      t_ftrapfi(:,:,:) = 0.0 

D      allocate(t_zetaz(kz,ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_zetaz",myid,istat)
D      t_zetaz(:,:,:) = 0.0 

D      allocate(t_olossc(kz,ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_olossc",myid,istat)
D      t_olossc(:,:,:) = 0.0 

D      allocate(t_hicmz(kz,ke,kb,kcm),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_hicmzz",myid,istat)
D      t_hicmz(:,:,:,:) = 0.0 

D      allocate(t_fber(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_fber",myid,istat)
D      t_fber(:,:) = 0.0
 
D      allocate(t_nmbrz(kz),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_nmbrz",myid,istat)
D      t_nmbrz(:) = 0.0 

D      allocate(t_fb11(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_fb11",myid,istat)
D      t_fb11(:,:) = 0.0
 
D      allocate(t_wb11(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_wb11",myid,istat)
D      t_wb11(:,:) = 0.0
 
D      allocate(t_fb10(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_fb10",myid,istat)
D      t_fb10(:,:) = 0.0
 
D      allocate(t_wb10(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_wb10",myid,istat)
D      t_wb10(:,:) = 0.0

D      allocate(t_fb01(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_fb01",myid,istat)
D      t_fb01(:,:) = 0.0
 
D      allocate(t_wb01(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_wb01",myid,istat)
D      t_wb01(:,:) = 0.0
 
D      allocate(t_fb00(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_fb00",myid,istat)
D      t_fb00(:,:) = 0.0
 
D      allocate(t_wb00(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_wb00",myid,istat)
D      t_wb00(:,:) = 0.0

D      allocate(t_fap(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_fap",myid,istat)
D      t_fap(:,:) = 0.0

D      allocate(t_fwall(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_fwall",myid,istat)
D      t_fwall(:,:) = 0.0

D      allocate(t_forb(ke,kb),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_forb",myid,istat)
D      t_forb(:,:) = 0.0

D      allocate(t_xpts(maxp),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_xpts",myid,istat)
D      t_xpts(:) = 0.0

D      allocate(t_ypts(maxp),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_ypts",myid,istat)
D      t_ypts(:) = 0.0

D      allocate(t_zpts(maxp),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_zpts",myid,istat)
D      t_zpts(:) = 0.0

D      allocate(t_vx(maxp),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_vx",myid,istat)
D      t_vx(:) = 0.0

D      allocate(t_vy(maxp),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_vy",myid,istat)
D      t_vy(:) = 0.0

D      allocate(t_vz(maxp),STAT = istat)
D      if(istat .ne. 0)
D    .          call allocate_error("t_vz",myid,istat)
D      t_vz(:) = 0.0


      zero          = 0.0
      one           = 1.0
      one_millionth = 1.0e-6
      maxpp = maxp                                    !for MPI use, each process supplies maxpp
      if(numprocs .gt. 1) maxpp = maxpp/numprocs      !elements for subsequent plotting
c
c  set up some data for optional (r,z) deposition calculation
c
      if (freyavb .gt. 0 .and. myid .eq. 0)  
     .       write (ncrt, '('' starting FREYA,time = '',1pe18.8)')time
      if (inubpat .eq. 1) ! not veryfied
     .  call setrz (npat,r(1),r(mi),z(1),z(mj),drpat,dzpat,nrpat,nzpat)
c
c  zero orbit output parameters
c  if output is desired, these need to be set accordingly
c
      norb = 0
****  norb = nout
      losseval = 1  ! if iborb=2 do setup calcs in GET_ORBLOSS   ! rsHSJ
      do 10 i=1,5
   10 idebug(i) = 0.0

      npart_full  =0.0 !array
      npart_half  =0.0 !array
      npart_third =0.0 !array

      
c
c     initialize flux surface quantities
c
      mfm1   = mf-1
      rmajor = rotsid(1)
      drin   = (rinsid(1)-rinsid(mf))/mfm1
      drot   = (rotsid(mf)-rotsid(1))/mfm1
      drutp  = SQRT (potsid(mf)-potsid(1))/mfm1
      drin   =  MAX (drin , one_millionth)
      drot   =  MAX (drot , one_millionth)
      drutp  =  MAX (drutp, one_millionth)
      elong  =  MAX (elong, one_millionth)
      drini  = 1.0 / drin
      droti  = 1.0 / drot
      drutpi = 1.0 / drutp
      elongi = 1.0 / elong
c
      if (codeid .ne. 'onedee') then
        mim1  = mi - 1
        mjm1  = mj - 1
        dr    = (r(mi)-r(1))/mim1
        dz    = (z(mj)-z(1))/mjm1
        dri   = 1.0 / dr
        dzi   = 1.0 / dz
      end if
c
c     calculate sines and cosines of various angles
c
      do ib=1,mb
        cangv(ib) = COS (anglev(ib)*pio180)
        cangh(ib) = COS (angleh(ib)*pio180)
        sangv(ib) = SIN (anglev(ib)*pio180)
        sangh(ib) = SIN (angleh(ib)*pio180)
      end do
c
c  calculate beam power; account for neutralizer efficiency
c
      call beam_prop(atwb,eb,bpow,vbeam,bion,bneut,bntot)

c
c     calculate macroscopic cross sections (modified for MPI - HSJ)
c
      call nbsgxn (ebkev,fe_tk,ibion,mb,mfm1,ne_tk,vbeam,
     .             zne,zni,zte,zzi,de_tk,hxfrac,sgxn,sgxnmi,atw_beam)

c
c     calculate total plasma volume
c
      volume = 0.0
      do 30 i=1,mfm1
   30 volume = volume + psivol(i)
c
c ----------------------------------------------------------------------
c     begin loop over beams
c ----------------------------------------------------------------------
c
      iskip = 1 + (npart-1)/maxp
      ic    = 0
      npts  = 0



         mlost1(:,:) =0
         mlost2(:,:) = 0 
      do 200 ib=1,mb

         if(time_dep_beam .eq.1) then
c           note that we assume that all sources of beam ib are on
c           or off together. (ie source2_phase .ne. 0.0  is not 
c           implemented at this time)
c           The logic in sub source dictated that Freya is called.
c           We now have to decide what beam(s) are actually on.
c           Skip those beams that are currently off:
            print *,'beamon , source2_p,ib',
     .                     beamOn(ib),source2_phase(ib),ib
            if(time .lt. beamon(ib)+source2_phase(ib))   go to 200
c            if(time .gt. beam_end(ib)+source2_phase(ib)) go to 200
c              skip beamlines/sources  that aren't currently on:
c               do nsrc = 1, nsourc
c                   insrc = nsrc
                   insrc = 1              ! source phase not implemented
                   do npls =1,kt
                      npulse(ib) = npls
                      if( time .ge. pbeamOn(npls,insrc,ib) .and. 
     .                  time .le. pbeamOff(npls,insrc,ib)) go to 202
                   enddo
c               enddo
               npulse(ib) = 0
               print *,'skipping beam #',ib,' at time',time
               go to 200  !source is not currently on,skip this beam

           endif



 202        if (freyavb .gt. 0 .and. myid .eq. master )
     .    write (ncrt, '('' master  starting injection for beam #'',
     .    i2, '' time  = '',1pe14.6)') ib,time
D         if (freyavb .gt. 0 .and. myid .ne. master )
D    .    write (ncrt, '('' slave starting injection for beam #'',
D    .    i2, '' time  = '',1pe14.6)') ib,time
c
        if (iborb .eq. 3 .and. myid .eq. master ) then ! MCGO may be run..
c                              ..open the file for saving deposition
          if (read_mcgo_file(ib) .eq. 0) then
            call DESTROY (mcgo_input_file2(ib))
            call getioun(nmcgo2,nmcgo2)
            open (unit = nmcgo2, file = mcgo_input_file2(ib),
     .          status = 'NEW' , form = 'UNFORMATTED')
            write(nmcgo2) ebkev(ib)
          else
            go to 200 ! iborb = 3 and read_mcgo_file(ib) = 1 ===>
c                       MCGO results to be read in
          end if
        end if
c
c       Determine aperture designators
c
        do i=1,nap

c            removed the following 4/1/2011 HSJ
c          if (nashape(i,ib) .eq. 's-circ' ) iatype(i,ib) = 1
c          if (nashape(i,ib) .eq. 's-rect' ) iatype(i,ib) = 2
c          if (nashape(i,ib) .eq. 's-vert' ) iatype(i,ib) = 3
c          if (nashape(i,ib) .eq. 's-horiz') iatype(i,ib) = 4
c          if (nashape(i,ib) .eq. 'b-circ' ) iatype(i,ib) = 5
c          if (nashape(i,ib) .eq. 'b-rect' ) iatype(i,ib) = 6
c          if (nashape(i,ib) .eq. 'b-vert' ) iatype(i,ib) = 7
c          if (nashape(i,ib) .eq. 'b-horiz') iatype(i,ib) = 8
c          if (nashape(i,ib) .eq. 'b-d3d'  ) iatype(i,ib) = 9

c replaced above by Bob Harvey's corrections 4/1/2011 HSJ:
c There is also a correction in sub rotate
         IF (nsourc.eq.1) then
             if(nashape(i,ib).eq.'s-circ')  iatype(i,ib)=5
             if(nashape(i,ib).eq.'s-rect')  iatype(i,ib)=6
             if(nashape(i,ib).eq.'s-vert')  iatype(i,ib)=7
             if(nashape(i,ib).eq.'s-horiz') iatype(i,ib)=8
          ELSE IF (nsourc.gt.1) then
             if(nashape(i,ib).eq.'s-circ')  iatype(i,ib)=1
             if(nashape(i,ib).eq.'s-rect')  iatype(i,ib)=2
             if(nashape(i,ib).eq.'s-vert')  iatype(i,ib)=3
             if(nashape(i,ib).eq.'s-horiz') iatype(i,ib)=4
          ENDIF
          IF (nashape(i,ib) .EQ. 'b-circ' ) iatype(i,ib) = 5
          IF (nashape(i,ib) .EQ. 'b-rect' ) iatype(i,ib) = 6
          IF (nashape(i,ib) .EQ. 'b-vert' ) iatype(i,ib) = 7
          IF (nashape(i,ib) .EQ. 'b-horiz') iatype(i,ib) = 8
          IF (nashape(i,ib) .EQ. 'b-d3d'  ) iatype(i,ib) = 9
        end do
c




c
c       Some angles for subroutine ROTATE
c       NOTE that indentation is messed up from here to 200 CONTINUE
c
      thetp(ib)  =
     .  ATAN2 (bvofset(ib), SQRT (bleni(ib)**2-bvofset(ib)**2))
      costp(ib)  = COS (thetp(ib))
      sintp(ib)  = SIN (thetp(ib))
      thetpp(ib) =
     .  ATAN2 (bhofset(ib), SQRT (bleni(ib)**2-bvofset(ib)**2))
      costpp(ib) = COS (thetpp(ib))
      sintpp(ib) = SIN (thetpp(ib))






c
c ----------------------------------------------------------------------
c begin loop over beam energy components
c ----------------------------------------------------------------------



      do 201 ie=1,3
D      if (freyavb .gt. 0 .and. myid .eq. master)
D    .  write (ncrt, '('' master starting energy component'', i2)') ie
D      if (freyavb .gt. 0 .and. myid .ne. master)
D    .  write (ncrt, '('' slave starting energy component'', i2)') ie
      fap(ie,ib)   = 0.0
      fwall(ie,ib) = 0.0
      forb(ie,ib)  = 0.0
      fb11(ie,ib)  = 0.0
      fb10(ie,ib)  = 0.0
      fb01(ie,ib)  = 0.0
      fb00(ie,ib)  = 0.0
      fber(ie,ib)  = 0.0
      wb11(ie,ib)  = 0.0
      wb10(ie,ib)  = 0.0
      wb01(ie,ib)  = 0.0
      wb00(ie,ib)  = 0.0

      do 110 i=1,mfm1
         ftrapfi(i,ie,ib) = 0.0
         hibrz(i,ie,ib)   = 0.0
         hdepz(i,ie,ib)   = 0.0
         angmpz(i,ie,ib)  = 0.0
         zetaz(i,ie,ib)   = 0.0
         hicmz(i,ie,ib,1) = 0.0
         hicmz(i,ie,ib,2) = 0.0
         hicmz(i,ie,ib,3) = 0.0
         olossc(i,ie,ib)  = 0.0                                        ! rs
  110 continue

      do 310 i=1,mfm1
  310 nmbrz(i) = 0
c
c --- nmbrz(i) counts ions born in zone i
c
c ----------------------------------------------------------------------
c begin loop over particles
c ----------------------------------------------------------------------
c
      npar = (bneut(ie,ib)/bntot)*npart      !#neutrals to launch for this
****  if (npar .eq. 0)  go to 200            !energy component
      if (npar .eq. 0)  go to 201
      nparx  = 0
      newpar = 0


      
      do 180 ipar = myid+1,npar,numprocs       !particle loop starts here

c
      if (MOD (ipar-1,npskip) .eq. 0)  newpar = 1
      if (newpar .eq. 0)  go to 120
c
c  generate neutral particle at beam source
c
      call sorspt(nbshape,bheigh,bwidth,bhfoc,bvfoc,bhdiv,bvdiv,ib,ie,
     .            isourc,ke,nsourc,sfrac1,vbeam,x0,y0,z0,vx0,vy0,vz0)
c
c  transform coordinates and advance particle to pivot point
c

      call rotate(naptr,iatype,aheigh,awidth,alen,bhofset,bvofset,cangv,
     .            cangh,ib,isourc,costp,sintp,costpp,sintpp,blenp,
     .            nsourc,sangv,sangh,rpivot,zpivot,mlost,x0,y0,z0,vx0,
     .            vy0,vz0,mlost1(ie,ib),mlost2(ie,ib))
c

c  skip injection if particle is lost at aperture
c
  120 if (mlost .ne. 0)  go to 160

c
c  inject particle into plasma, i.e., follow particle from pivot
c     point into or through the plasma
c

      call inject(atwb, codeid, de_tk,drutpi,droti,dri,ds_tk,dzi,
     .            elongi,ib,ie,kbe,ke,kz,nw,mfm1,mim1,mjm1,ne_tk,newpar,
     .            nout,potsid(1),psi,r,rmajor,rin,rmax,sgxn,sgxnloc,
     .            sgxnmi,x0,y0,z0,vx0,vy0,vz0,vbeam,z,zangrot,zax,zmin,
     .            zmax,izone,pzone,rzone,rpos,xpos,ypos,zpos,myid,
     .            tenter,smax,texit)



c
c  skip birth data if:  particle missed plasma
c
      if (izone .ge. mf)  go to 170
c
c  accumulate hot ion birth rate and creation mode
c
 
      hibrz(izone,ie,ib)   = hibrz(izone,ie,ib) + 1.0
      hicmz(izone,ie,ib,1) = hicmz(izone,ie,ib,1) + sgxnloc(1)
      hicmz(izone,ie,ib,2) = hicmz(izone,ie,ib,2) + sgxnloc(2)
      hicmz(izone,ie,ib,3) = hicmz(izone,ie,ib,3) + sgxnloc(3)
D     t_hibrz(izone,ie,ib) = hibrz(izone,ie,ib)
D     t_hicmz(izone,ie,ib,1) = hicmz(izone,ie,ib,1)
D     t_hicmz(izone,ie,ib,2) = hicmz(izone,ie,ib,2)
D     t_hicmz(izone,ie,ib,3) = hicmz(izone,ie,ib,3)
c
c  calculate and accumulate pitch angle at birth point
c  calculate angular momentum deposited by each monte carlo ion
c
c  NOTE AS OF 10/18/94 THE LINES MARKED WITH ! rs WERE ADDED. SOME OF
c  THESE CHANGES WILL AFFECT THE NEUTRAL BEAM RESULTS (PRESUMABLY ONLY
c  SLIGHTLY) EVEN WHEN THE STAMBAUGH ORBIT LOSS MODEL IS NOT USED. ... HSJ
c    should replace zetai this with  full formula HSJ
c
****  zetai  = csgn*(xpos*vy0-ypos*vx0)/(rpos*vbeam(ie,ib))   ! ORIGINAL
      vplane = SQRT (vx0**2+vy0**2)                           ! rs
      zetai  = csgn*(xpos*vy0-ypos*vx0)/(rpos*vplane)         ! rs
      zetai  = MIN (zetai,  one)
      zetai  = MAX (zetai, -one)
****  vtroid = csgn*zetai*vbeam(ie,ib)                        ! ORIGINAL
      vtroid = csgn*zetai*vplane                              ! rs
      vrad   = -SQRT (1.0 - zetai**2) * vplane                ! rs
      angmtm = rpos*vtroid*atwb*1.673e-24
      angmpz(izone,ie,ib) = angmpz(izone,ie,ib) + angmtm
****  if (iborb .eq. 0)                                       ! ORIGINAL
      if (iborb .ne. 1)                                       ! rs
     .  zetaz(izone,ie,ib) =  zetaz(izone,ie,ib) + zetai
D      t_zetaz(izone,ie,ib) =  zetaz(izone,ie,ib)
c
c  save occasional birth point for subsequent plotting
c
      ic = ic + 1
      if (MOD (ic-1, iskip) .eq. 0 ) then
        npts       = npts + 1
        npts=min(npts,maxpp)            !dont let the arrays overflow (numprocs =1)
        xpts(npts) = xpos               ! or save maxpp pts from each process (numprocs > 1)
        ypts(npts) = ypos
        zpts(npts) = zpos
        vx(npts)   = vx0
        vy(npts)   = vy0
        vz(npts)   = vz0
        pitch_a(npts) = zetai
D       if(numprocs .gt. 1)then
D          t_xpts(myid*maxpp+npts) = xpts(npts)
D          t_ypts(myid*maxpp+npts) = ypts(npts)
D          t_zpts(myid*maxpp+npts) = zpts(npts)
D          t_vx(myid*maxpp+npts) = vx(npts)
D          t_vy(myid*maxpp+npts) = vy(npts)
D          t_vz(myid*maxpp+npts) = vz(npts)
D       endif
      end if
c
c    save results for mcgo:
c
      if (iborb .eq. 3 .and. myid .eq. master) then
         if (read_mcgo_file(ib) .eq. 0) then
             if (ie .eq. 1)  npart_full(ib)=npart_full(ib)+1
             if (ie .eq. 2)  npart_half(ib)=npart_half(ib)+1
             if (ie .eq. 3)  npart_third(ib)=npart_third(ib)+1
             ppos = ATAN2 (ypos,xpos) ! MCGO converts to degrees
             write(nmcgo2) rpos
             write(nmcgo2) zpos 
             write(nmcgo2) ppos
             write(nmcgo2) zetai
             write(nmcgo2) pzone
             write(nmcgo2) ie
         end if
      end if
c
c  calculate (r,z) grid location
c
      if (inubpat .gt. 0 .and. myid .eq. master) then
        i = (rpos-r(1))/drpat + 1.0
        j = (zpos-z(1))/dzpat + 1.0
        rzpat(i,j,ie,ib) = rzpat(i,j,ie,ib) + 1.0
      end if
c
c ----------------------------------------------------------------------
c
c  calculate the beam orbit loss using OLOSSA subroutine         ! rs
c
      if (iborb .ne. 2)  go to 150                               ! rs
c
c  ignore ions born inside mag. axis, assume they are confined   ! rs
c
      if (rpos .le. rotsid(1))  go to 140                        ! rsHSJ
      call get_orbloss (losseval, izone, rpos, zpos, vtroid,
     .                  vrad, vz0, vbeam(ie,ib), phi, elossb)
      if (phi .gt. 180 .or. phi .lt. 90)  go to 140              ! rs
c
c  if the ion energy is greater than elossb, it is lost          ! rs
c
      if (eb(ie,ib) .lt. elossb)  go to 140                      ! rs
c
c  take away the angular momentum and zetai contributions        ! rs
c
      angmpz(izone,ie,ib) = angmpz(izone,ie,ib)-angmtm           ! rs
D     t_angmpz(izone,ie,ib) = angmpz(izone,ie,ib)
      zetaz (izone,ie,ib) = zetaz (izone,ie,ib)-zetai            ! rs
D     t_zetaz (izone,ie,ib) = zetaz (izone,ie,ib)
c
c  accumulate the radial current that flows in to neutralize the ! rs
c  electron left behind by the lost ion                          ! rs
c  accumulate this current as a number of particles passing      ! rs
c  through each zone                                             ! rs
c
      do  jl=izone,mfm1                                          ! rs
         olossc(jl,ie,ib) = olossc(jl,ie,ib) + 1.0               ! rs
D         t_olossc(jl,ie,ib) = olossc(jl,ie,ib)
      enddo
      go to 175                                                  ! rs
c
c  here we have a confined ion                                   ! rs
c
  140 hdepz(izone,ie,ib) = hdepz(izone,ie,ib) + 1.0              ! rs
D     t_hdepz(izone,ie,ib) = hdepz(izone,ie,ib)
      go to 180                                                  ! rs
c
c ----------------------------------------------------------------------
c
c  calculate orbit widths and orbit loss
c
****  if (iborb .ne. 0)
  150 if (iborb .eq. 1) then                                     ! rs
        nparx  = nparx + 1
        call orbit_12(atwb,b1ins,b1ots,b2ins,b2ots, codeid,ic,idebug,
     .                iskip,
     .             izone,mf,norb,pinsid,potsid,pzone,rinsid,rotsid,
     .             rzone,rpos,vbeam(ie,ib),zetai,zpos,fin,fot,zin,zot,
     .             ipass,iaxis,ier,izp,wid,wt,zetorb)
        if (ier .ne. 0) then
          fber(ie,ib) = fber(ie,ib) + 1.0
D         t_fber(ie,ib) = fber(ie,ib)
        else if (izp .gt. mfm1) then
          go to 175
        else
          nmbrz(izp) = nmbrz(izp) + 1
D         t_nmbrz(izp) = nmbrz(izp)
          if      (ipass .eq. 1 .and. iaxis .eq. 1) then
            fb11(ie,ib) = fb11(ie,ib) + 1.0
            wb11(ie,ib) = wb11(ie,ib) + wid
D           t_fb11(ie,ib) = fb11(ie,ib)
D           t_wb11(ie,ib) = wb11(ie,ib)
          else if (ipass .eq. 1 .and. iaxis .eq. 0) then
            fb10(ie,ib) = fb10(ie,ib) + 1.0
            wb10(ie,ib) = wb10(ie,ib) + wid
D           t_fb10(ie,ib) = fb10(ie,ib)
D           t_wb10(ie,ib) = wb10(ie,ib)
          else if (ipass .eq. 0 .and. iaxis .eq. 1) then
            fb01(ie,ib) = fb01(ie,ib) + 1.0
            wb01(ie,ib) = wb01(ie,ib) + wid
            ftrapfi(izp,ie,ib) = ftrapfi(izp,ie,ib) + 1.0
D           t_fb01(ie,ib) = fb01(ie,ib)
D           t_wb01(ie,ib) = wb01(ie,ib)
D           t_ftrapfi(izp,ie,ib) = ftrapfi(izp,ie,ib)
          else if (ipass .eq. 0 .and. iaxis .eq. 0) then
            fb00(ie,ib) = fb00(ie,ib) + 1.0
            wb00(ie,ib) = wb00(ie,ib) + wid
            ftrapfi(izp,ie,ib) = ftrapfi(izp,ie,ib)+1.0
D           t_fb00(ie,ib) = fb00(ie,ib)
D           t_wb00(ie,ib) = wb00(ie,ib)
D           t_ftrapfi(izp,ie,ib) = ftrapfi(izp,ie,ib)
          end if
        end if
c
c  accumulate hot ion deposition rate and average pitch angle cosine
c
        do i=1,mfm1
          hdepz(i,ie,ib) = hdepz(i,ie,ib) + wt(i)
D         t_hdepz(i,ie,ib) = hdepz(i,ie,ib)
c         angmpz(i,ie,ib) = angmpz(i,ie,ib)+wt(i)*angmtm
          zetaz(i,ie,ib) = zetaz(i,ie,ib) + wt(i)*zetorb(i)
D         t_zetaz(i,ie,ib) = zetaz(i,ie,ib)
        end do
      end if                !iborb=1
      go to 180
c
c  accumulate particles that are not deposited in plasma
c
  160 fap(ie,ib)   = fap(ie,ib) + 1.0
D     t_fap(ie,ib)   = fap(ie,ib)
      go to 180
  170 fwall(ie,ib) = fwall(ie,ib) + 1.0
D     t_fwall(ie,ib) = fwall(ie,ib)
      go to 180
  175 forb(ie,ib)  = forb(ie,ib) + 1.0
D     t_forb(ie,ib)  = forb(ie,ib)
c
c ----------------------------------------------------------------------
c end loop over particles
c ----------------------------------------------------------------------
c
  180 newpar = 0

D     call MPI_Barrier( MPI_COMM_WORLD,ierr) !wait until all processors are done

D     if(numprocs .gt. 1)then
C     hibrz(:,ie,ib)=0.0   !explicit zeronot necessary
D     call MPI_Allreduce(t_hibrz(1,ie,ib),hibrz(1,ie,ib),kz,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)



C     angmpz(:,ie,ib)=0.0   !explicit zeronot necessary
D     call MPI_Allreduce(t_angmpz(1,ie,ib),angmpz(1,ie,ib),kz,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

C      hicmz(:,ie,ib,1) = 0.0   !explicit zeronot necessary
D     call MPI_Allreduce(t_hicmz(1,ie,ib,1),hicmz(1,ie,ib,1),kz,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

C      hicmz(:,ie,ib,2) = 0.0   !explicit zeronot necessary
D     call MPI_Allreduce(t_hicmz(1,ie,ib,2),hicmz(1,ie,ib,2),kz,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

C      hicmz(:,ie,ib,3) = 0.0   !explicit zeronot necessary
D     call MPI_Allreduce(t_hicmz(1,ie,ib,3),hicmz(1,ie,ib,3),kz,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

C     hdepz(:,ie,ib)=0.0   !explicit zeronot necessary
D     call MPI_Allreduce(t_hdepz(1,ie,ib),hdepz(1,ie,ib),kz,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)


D     call MPI_Allreduce(t_nmbrz,nmbrz,kz,MPI_Integer,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_ftrapfi(1,ie,ib),ftrapfi(1,ie,ib),kz,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_olossc(1,ie,ib),olossc(1,ie,ib),kz,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_fap(ie,ib),fap(ie,ib),1,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_fwall(ie,ib),fwall(ie,ib),1,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_forb(ie,ib),forb(ie,ib),1,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_fb10(ie,ib),fb10(ie,ib),1,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)


D     call MPI_Allreduce(t_fb11(ie,ib),fb11(ie,ib),1,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_fb00(ie,ib),fb00(ie,ib),1,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_fb01(ie,ib),fb01(ie,ib),1,
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D      endif   !numprocs  > 1 

c
c  normalize average pitch angle cosine.  normalize momentum and birth
c     mode in each shell to a single particle.
c
      do 182 i=1,mfm1
          xnorm = hibrz(i,ie,ib)
          if (xnorm .ne. 0.0) then
              angmpz(i,ie,ib)  = angmpz(i,ie,ib)/xnorm
              hicmz(i,ie,ib,1) = hicmz(i,ie,ib,1)/xnorm
              hicmz(i,ie,ib,2) = hicmz(i,ie,ib,2)/xnorm
              hicmz(i,ie,ib,3) = hicmz(i,ie,ib,3)/xnorm
          end if
          if (iborb .ne. 0)  xnorm = hdepz(i,ie,ib)
          if (xnorm .ne. 0.0) xnorm = 1.0/xnorm
****      angmpz(i,ie,ib) = xnorm*angmpz(i,ie,ib)
  182 zetaz(i,ie,ib) = xnorm*zetaz(i,ie,ib)
c
c  get the fraction of trapped ions in each zone.  if orbit effects are
c     not turned on or if itrapfi = 0 then this effect is not included.
c
      do i=1,mfm1
        nmbrz(i) = MAX0 (nmbrz(i),1) 
        ftrapfi(i,ie,ib) = ftrapfi(i,ie,ib)/nmbrz(i)
      end do
c
c  normalize loss fractions and hot ion birth rate
c
      fap(ie,ib)   = fap(ie,ib)/npar
      fwall(ie,ib) = fwall(ie,ib)/npar
      forb(ie,ib)  = forb(ie,ib)/npar
      xloss1       = fap(ie,ib) + fwall(ie,ib)
      xloss2       = fap(ie,ib) + fwall(ie,ib) + forb(ie,ib)
****  if (xloss1 .ge. 1.0)  go to 200
      if (xloss1 .ge. 1.0)  go to 201
c
      do i=1,mfm1
        hibrz(i,ie,ib) =
     .  hibrz(i,ie,ib)*volume/((1.0-xloss1)*npar*psivol(i))
        if      ( iborb .eq. 0  ) then
          hdepz(i,ie,ib) = hibrz(i,ie,ib)
        else if (xloss2 .lt. 1.0) then
          hdepz(i,ie,ib) =
     .    hdepz(i,ie,ib)*volume/((1.0-xloss2)*npar*psivol(i))
        end if
      if (iborb .eq. 2)  olossc(i,ie,ib) = olossc(i,ie,ib) / npar   ! rs
      end do
c
c  normalize orbit widths and fractions
c
      if (iborb  .eq. 1) then
        if (fb11(ie,ib) .ne. 0.0)  wb11(ie,ib) = wb11(ie,ib)/fb11(ie,ib)
        if (fb10(ie,ib) .ne. 0.0)  wb10(ie,ib) = wb10(ie,ib)/fb10(ie,ib)
        if (fb01(ie,ib) .ne. 0.0)  wb01(ie,ib) = wb01(ie,ib)/fb01(ie,ib)
        if (fb00(ie,ib) .ne. 0.0)  wb00(ie,ib) = wb00(ie,ib)/fb00(ie,ib)
        fb11(ie,ib) = fb11(ie,ib)/nparx
        fb10(ie,ib) = fb10(ie,ib)/nparx
        fb01(ie,ib) = fb01(ie,ib)/nparx
        fb00(ie,ib) = fb00(ie,ib)/nparx
        fber(ie,ib) = fber(ie,ib)/nparx
      end if
c
      if (iborb .eq. 3 .and. read_mcgo_file(ib) .eq. 0 )
     .  write(nmcgo2)(hibrz(i,ie,ib),i=1,mfm1)

  201 continue                ! end loop over energy components


      if (iborb .eq. 3 .and. read_mcgo_file(ib) .eq. 0 ) then
        call giveupus(nmcgo2)
        close (unit = nmcgo2) ! close the fast ion birth file
      end if
  200 continue                ! end loop over beams
c      print *,'mlost1,mlost2 =',mlost1(1:3,1:2),mlost2(1:3,1:2)
 

c
c  renormalize currents and powers to bptor
c
      do ib=1,mb
        if (bptor(ib) .gt. 0.0) then
          bptorx = 0.0
          do 210 ie=1,3
  210     bptorx = bptorx + (1.0-fap(ie,ib))*bpow(ie,ib)
          if (bptorx .gt. 0.0) then
            xnorm = bptor(ib)/bptorx
            bcur(ib) = xnorm*bcur(ib)
            do 220 ie=1,3                                        ! rs
              bion (ie,ib) = xnorm*bion(ie,ib)
              bneut(ie,ib) = xnorm*bneut(ie,ib)
              bpow (ie,ib) = xnorm*bpow(ie,ib)
              if (iborb .ne. 2)  go to 220                       ! rs
              do 219 j=1,mfm1                                    ! rs
**219         olossc(j,ie,ib) = olossc(j,ie,ib)*xnorm*bneut(ie,ib) ! rs
**** .                        *1.6e-19/psivol(j)                 ! rs
  219         olossc(j,ie,ib) = olossc(j,ie,ib)*bneut(ie,ib)     ! rsHSJ
     .                        *1.6e-19/psivol(j)                 ! rs
c             return olossc in amperes per cm**3                 ! rs
  220         continue                                           ! rs
          end if
        end if
      end do
c
      if (iborb .eq. 3) then
c
         do ib=1,mb
c
c          prepare an input file for MCGO for each beam in ONETWO:
c          (the values of bnuet, etc. from above are required in
c          prep_mcgo.)
c
           if (read_mcgo_file(ib) .eq. 0)
     .         call prep_mcgo (ib, bpow, ke, enbeams)
           if (spawn_mcgo .gt. 0) then ! spawn_mcgo is set in sub INIT
                 spawn_command = mcgo_path(1:LENGTH(mcgo_path      )) //
     .       '/mcgo < ' 
     .      // mcgo_input_file(ib)(1:LENGTH(mcgo_input_file))
             if (freyavb .gt. 0)
     .           write(ncrt,'(''spawning mcgo:'',/,
     .           2x,a,2x,
     .           /,''  with beam number '',i3,//)')
     .           spawn_command(1:LENGTH(spawn_command)),ib
c
c           call mcgo with the two files,
c           mcgo_input_12(ib) and mcgo_input_file2(ib), for each beam:
c
            write (ncrt, '(a)') 'mcgo_path = ' // mcgo_path
            if (ISHELL (spawn_command(1:LENGTH(spawn_command))) .ne. 0)
     .      call STOP ('subroutine FREYA: failure of spawned MCGO', 270)
c
           end if
c
c          read the file produced by MCGO unless
c          we are just generating MCGO input files:
c          the files to be read were created by the above spawn or they
c          existed a priori and were specified in inone:
c
           if (spawn_mcgo .ne. 0 .or. read_mcgo_file(ib) .ne. 0)
     .         call read_mcgo_files (ib, eb, ke, bpow)
         end do
      end if
c
c     calculate neutral beam density on (r,z) grid
c
      if (inubpat .eq. 1 .and. codeid .ne. 'onedee' .and. iborb .ne. 3
     .     .and. myid .eq. master )
     .  call nbdep2d (psi, mi, mj, r, z, potsid, mf, rzpat, nrpat,
     .                nzpat, ke, mb, sgxn, vbeam, hxfrac)


D     call MPI_Allreduce(t_xpts,xpts,maxp,                     !(inplace reduction not until MPI 2.0)
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)


D     call MPI_Allreduce(t_ypts,ypts,maxp,                     !(inplace reduction not until MPI 2.0)
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)
     
D     call MPI_Allreduce(t_zpts,zpts,maxp,                     !(inplace reduction not until MPI 2.0)
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_vx,vx,maxp,                         !(inplace reduction not until MPI 2.0)
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_vy,vy,maxp,                         !(inplace reduction not until MPI 2.0)
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

D     call MPI_Allreduce(t_vz,vz,maxp,                         !(inplace reduction not until MPI 2.0)
D    &                       MPI_Double_Precision,
D    &                       MPI_SUM,MPI_COMM_WORLD,ierr)

c
c      get rid of arrays used for MPI:
D      deallocate(t_hibrz,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_hibrz",myid,istat) 


D      deallocate(t_hdepz,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_hdepz",myid,istat)


D      deallocate(t_angmpz,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_angmpz",myid,istat)


D      deallocate(t_ftrapfi,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_ftrapfi",myid,istat)


D      deallocate(t_zetaz,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_zetaz",myid,istat)


D      deallocate(t_olossc,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_olossc",myid,istat)


D      deallocate(t_hicmz,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_hicmzz",myid,istat)

D      deallocate(t_fber,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_fber",myid,istat)

 
D      deallocate(t_nmbrz,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_nmbrz",myid,istat)


D      deallocate(t_fb11,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_fb11",myid,istat)

 
D      deallocate(t_wb11,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_wb11",myid,istat)

 
D      deallocate(t_fb10,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_fb10",myid,istat)

 
D      deallocate(t_wb10,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_wb10",myid,istat)

D      deallocate(t_fb01,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_fb01",myid,istat)

 
D      deallocate(t_wb01,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_wb01",myid,istat)

 
D      deallocate(t_fb00,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_fb00",myid,istat)

 
D      deallocate(t_wb00,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_wb00",myid,istat)


D      deallocate(t_fap,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_fap",myid,istat)

D      deallocate(t_fwall,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_fwall",myid,istat)

D      deallocate(t_forb,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_forb",myid,istat)


D      deallocate(t_xpts,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_xpts",myid,istat)


D      deallocate(t_ypts,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_ypts",myid,istat)

D      deallocate(t_zpts,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_zpts",myid,istat)

D      deallocate(t_vx,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_vx",myid,istat)


D      deallocate(t_vy,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_vy",myid,istat)

D      deallocate(t_vz,STAT = istat)
D      if(istat .ne. 0)
D    .          call deallocate_error("t_vz",myid,istat)


      if (freyavb .gt. 0 )  write (ncrt, '('' done with FREYA'')')


      return

      end

      real*8 function fsgxncx (atw, e, zni)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subprogram calculates inverse mean free path due to charge exchange
c from the fitted results of freeman and jones, clm-r137, culham (1974).
c
c     input:
c             atw - atomic weight of target ion
c             e   - relative energy of impinging neutral (ev/amu)
c             zni - density of target ion (cm**-3)
c ----------------------------------------------------------------------
c
      if (atw .gt. 3.01) then
        sigcx = 0.0
      else
        aloge = LOG10 (e)
        sigcx = 0.6937e-14 * (1.0 - 0.155*aloge)**2 /
     .                       (1.0 + 0.1112e-14*e**3.3)
      end if
      fsgxncx = sigcx*zni
      return
c
      end

      real*8 function fsgxne (vb, te, zne)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subprogram evaluates local inverse mean free path for electron impact
c ionization from fitted results of freeman and jones, clm-r137,culham (1974).
c
c     input:
c             vb  - speed of impinging neutral (cm/sec)
c             te  - target electron temperature (ev)
c             zne - target electron density (cm**-3)
c ----------------------------------------------------------------------
c
      dimension cfione(7)
      data      cfione /-3.173850e+01,  1.143818e+01, -3.833998    ,
     .                   7.046692e-01, -7.431486e-02,  4.153749e-03,
     .                  -9.486967e-05/
c
      alogt = 0.0
      if (te .gt. 1.0    )  alogt = LOG (te)
      if (te .gt. 1.0e+05)  alogt = 11.51
      expo = (((((cfione(7) *alogt + cfione(6))*alogt + cfione(5))*alogt
     .          + cfione(4))*alogt + cfione(3))*alogt + cfione(2))*alogt
     .          + cfione(1)
      fsgxne = EXP (expo) * zne / vb
      return
c
      end

      real*8 function fsgxni (atw, eova, zni, zzi)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  this subprogram calculates inverse mean free path due to proton and
c     impurity impact ionization.  proton impact ionization cross
c     sections from fitted results of freeman and jones,clm-r137, culham
c     (1974).  impurity impact ionization cross sections from r.e. olson
c     et al., Phys. Rev. Lett. 41, 163 (1978).
c
c     input:
c             atw  - atomic weight of target ion
c             eova - relative energy of impinging neutral (ev/amu)
c             zni  - density of target ion (cm**-3)
c             zzi  - average charge state of target ion
c ----------------------------------------------------------------------
c
      dimension cfionp(7)
      data      cfionp /-4.203309e+01,  3.557321    , -1.045134,
     .                   3.139238e-01, -7.454475e-02,  8.459113e-03,
     .                  -3.495444e-04/
c
      if (atw .le. 3.01) then
        aloge = LOG10 (eova)
        aloge = aloge * 2.302585093 - 6.907755279
        if (aloge .le. -2.30258) then
          sigi = 0.0
        else
          expo = (((((cfionp(7) *aloge + cfionp(6))*aloge
     .              + cfionp(5))*aloge + cfionp(4))*aloge
     .              + cfionp(3))*aloge + cfionp(2))*aloge
     .              + cfionp(1)
          sigi = EXP (expo)
        end if
        fsgxni = sigi*zni
      else
        ekev   = 1.0e-3*eova
        fsgxni = 1.0e-17*zni*46.0*zzi*(32.0*zzi/ekev)*
     .              (1.0 - EXP (-ekev/(32.0*zzi)))
      end if
      return
c
      end



      real *8 function ftau0(a,ts,vc)
      implicit none
      real *8 a,vc,ts
        ftau0 = 0.3333333*ts*
     .                 LOG((a**3+vc**3)/vc**3)
        return
        end

      real*8 function gef (vpar, tpar)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  This routine evaluates the function Ge, which is related to the
c     transfer of energy from fast ions slowing down on electrons.  The
c     function is defined by Callen et al., IAEA Tokyo, Vol. I, 645
c     (1974).
c ----------------------------------------------------------------------
c
      external       gefun
      common /gecom/ vcvo, tstcx, emzrat
      data           root13 /0.577350/, fact4 /0.302300/
c
****  root13 = SQRT (1.0/3.0)
****  fact4  = root13 * ATAN (root13)
c
      if (tpar .le. 0.01) then
        rooty = 1.0 / vpar
        y     = rooty**2
        term2 = LOG ((1.0-rooty+y)/(1.0 + rooty)**2)/(6.0 * y)
        term3 = root13* ATAN (root13*(2.0*rooty-1.0))/y
        term4 = fact4/y
        gef   = 2.0 * (0.5-term2-term3-term4)
      else
        vcvo  = vpar
        tstcx = tpar
        a     = 0.0
        b     = 1.0
        ep    = 0.1
        m4    = 3
        gef   = asimp(a,b,ep,m,m4,gefun)
      end if
      return
c
      end

      real*8 function gefun (y)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  gefun is the argument of the Ge integral given by Callen.  the
c     substitution y = v/v0 has been made and the ratio of slowing down
c     time to charge exchange time, tstcx, is assumed to be independent
c     of speed.  the integral determing probability against charge
c     exchange can then be done analytically and leads to the expression
c     for pcx given here.  gefun is integrated using subroutine ASIMP.
c     for small densities (such as occur in large h mode gradients) tstcx
c     can become quite large (>1000) which causes numerical problems.
c     hence this routine was modified 4/26/89 by HSJ
c ----------------------------------------------------------------------
c
      common /gecom/ vcvo, tstcx, emzrat
      data           alog2 /0.693147181/
c
      if (y .le. 0.0) then
        gefun = 0.0
      else
        v3       = vcvo**3
        arg      = (1.0 + v3)/(y**3+v3)
        pcxlog   = -tstcx*LOG (arg)*0.33333333333334
        alogy3v3 = LOG (y**3+v3)
        gefunlog = alog2 + 4.0 * LOG (y) + pcxlog - alogy3v3
        if (gefunlog .lt. -30.0) then
          gefun = 0.0
        else
          gefun = EXP (gefunlog)
        end if
      end if
****  pcx = (arg)**(-tstcx/3.0)
****  gefun = 2.0 * y**4*pcx/(y**3+v3)
      return
c
      end

      subroutine getsgxn (e, iz, ns, ktk, ib, ie, sgxn, nbins,
     .                    debin, kbe, kz, nout, sgxnloc)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subroutine locates n*sigma in the lookup table, sgxn
c ----------------------------------------------------------------------
c
c --- input through argument list:
c          debin          - energy bin width     (keV/amu)
c          e(ns)          - energy to be located (keV/amu)
c          ib             - beam line index
c          ie             - beam energy group index
c          iz(ns)         - flux zone index
c          nbins          - total number of energy bins
c          ns             - total number of points to be evaluated
c          sgxn(i,iz,j,k) - array of neutral stopping data
c                           i  - data type index
c                                =1, fraction of interactions producing
c                                    electrons;
c                                =2, fraction of interactions producing
c                                    neutral species 1;
c                                =3, fraction of interactions producing
c                                    neutral species 2;
c                                =4, inverse mean free path (cm**-1)
c                           iz - flux zone index
c                           j  - beam/energy index, j = 3*(ib-1)+ie
c                           k  - energy bin index
c
c --- output through argument list:
c          sgxnloc(i)     - local neutral stopping data (cm**-1)
c                           i - data type index (see above)
c                           note:  iftns>1, only the inverse mean free
c                           path is evaluated (i = 4) and ns data points
c                           are passed.  this facilitates faster
c                           execution when evaluating the max imfp along
c                           the collisionless neutral trajectory (used
c                           in rotating discharges only).  otherwise
c                           more detailed information is passed in the
c                           first four elements only (ns = 1).
c ----------------------------------------------------------------------
c
c   Created:    Glenn Sager?               Date:  ??-???-????
c
c   changes from original version:
c      1) The code will now terminate if the maximum rotational energy
c         bin is exceeded. Detailed Error message provide to both the
c         OUTONE file (nout) and the standard output device (ncrt = 6).
c      2) Fixed some math in array index calculations.
c      3) Today is the 50th Aniversary of Hiroshima. NEVER AGAIN!
c
c                                  D.F. Finkenthal 6-AUG-95
c
c ----------------------------------------------------------------------
c
      dimension  e(*), iz(*), sgxn(4,kz,kbe,*), sgxnloc(*), imaxa(200)
      integer,save :: imax 
c
      ncrt = 6
c
      ind = 3 * (ib - 1) + ie
      if (ns .gt. 1) then
        do i=1,ns
          ibin       = e(i) / debin + 1.0
          imaxa(i)   = ibin
          ibin       = MIN0 (ibin, nbins)
          sgxnloc(i) = sgxn (4, iz(i), ind, ibin)
        end do
        imax = maxaf (imaxa, 1, ns)
      else       
        ibin       = e(1) / debin + 1.0
        imax       = MAX0 (ibin, imax )
        ibin       = MIN0 (ibin, nbins)
        sgxnloc(1) = sgxn (1, iz(1), ind, ibin)
        sgxnloc(2) = sgxn (2, iz(1), ind, ibin)
        sgxnloc(3) = sgxn (3, iz(1), ind, ibin)
        sgxnloc(4) = sgxn (4, iz(1), ind, ibin)
 
      end if

c
      if (imax .gt. nbins) then
        emax = amaxaf (e, 1, ns)
        write (nout, 110)  nbins * debin, emax
        write (ncrt, 110)  nbins * debin, emax
        call STOP ('subroutine GETSGXN: max bin energy exceeded', 51)
      end if
c
 110  format (/
     . ' ERROR in GETSGXN: maximum rotational energy bin exceeded'     /
     . ' The highest calculated energy bin is now ', f10.2, ' keV/amu' /
     . ' but GETSGXN tried to look up a value of  ', f10.2, ' keV/amu' /
     . ' The FE_TK input parameter must be increased to overcome',
     . ' this problem!'                                                /
     . ' ONETWO will be terminated to avoid further problems.')
c
      return
c
      end

      subroutine hexnb (istarx, iexcix, ilorenx, mstatx, ncontx,
     .                  er0x, tex, tix, numix, amix, denix,
     .                  numzx, izx, amzx, izstrx, denzx, bperpx,
     .                  ldene, ldeni, ldenz, lsvi, lsvz, lsve, lrad,
     .                  lngh, lngl, louthx, lcor,
     .                  nsigmav, lambda, hexfrac, ihxerr)

      USE cpub_dat
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     version 21
c
c --- author:
c       c. d. boley
c       pppl 1984
c --- modified by: (vax / modular)
c       r. m. weiland
c       pppl july 1985
c --- ref.:
c       c. d. boley, r. k. janev, and d. e. post, Phys. Rev. Letters
c       52, 534 (1984).
c --- calling sequence:
c
c     call hexnb (istart, iexcit, ilorent, mstate, ncont,
c    .            er0, tex, tix, numi, ami, deni,
c    .            numz, iz, amz, izstrp, denz, bperp,
c    .            kdene, kdeni, kdenz, ksvi, ksvz, ksve, krad,
c    .            ngh, ngl, nouthx, ncorin, nsigmav, lambda)
c
c --- input parameters:
c       istart = 1 to initialize rad rates & fine pts       0 otherwise
c       iexcit if =0 then ignore contribution of excited states
c              if =1 then include contributionof excited states.
c       ilorent =0 or 1       whether or not to calculate the maximum
c                        principal quantum number ns of the populated
c                        states as given by the lorentz ionization
c                        limit (1  = => yes)
c       mstate< parameter ms
c               for ilorent = 0, use ns = mstate+1
c       note: the operational significance of ns in the code is to
c             set an upper limit such that any excitations to levels n
c             higher than ns are counted as "ionizations".
c       ncont< parameter mc
c               upper bound to number of continuum states
c       er0    beam energy per amu (ev)
c       te     electron temperature (ev)
c       ti     ion temperature (ev)
c       numi   number of hydrogenic ion species
c       ami    masses of the hydrogenic ions (amu)
c       deni   densities of the hydrogenic ions (cm**-3)
c       numz   number of impurity species
c              if numz gt 0, the coronal data file 'coronb' is read.
c       iz     atomic numbers of the impurities
c       izstrp =0 for coronal equilibrium values of <z> and <zsq>
c              =1 for fully stripped impurities
c       amz    atomic mass number of the impurities [iff izstrp#0]
c       denz   densities of the impurities (cm**-3)
c       bperp  magnetic field perpendicular
c              to beamline (tesla (for ilorent = 1))
c              kdene  =0: discard electron reactions
c                     =1: include electron reactions (default)
c              kdeni  =0: discard ion (hydrogen) reactions
c                     =1: include ion reactions (default)
c              kdenz  =0: discard impurity reactions
c                     =1: include impurity reactions (default)
c              ksvi   =0: simple multiplication instead of full sigma-v
c                         integral for ions
c                     =1: full sigma-v integral for ions
c              ksvz   =0: simple multiplication instead of full sigma-v
c                         integral for impurities
c                     =1: full sigma-v integral for impurities
c              ksve   =0: sigma-v integrals for electron reactions
c                         involving 1s, 2s, and 2p       rates for other
c                         electron reactions from Vriens & Smeets.
c                     =1: full sigma-v integrals for electron reactions
c              krad   =0: discard hydrogen line radiation
c                     =1: include hydrogen line radiation .
c              ngh        order of gauss-hermite integrations for
c                         ion and impurity reactions (default 24)
c              ngl        order of gauss-laguerre integrations for
c                         electron reactions
c                         permissible values: 10,16,20,24
c       nouthx unit number for messages(to unit nouthx if > 0)
c              (also nouthx>0 gives detailed diagnostic output)
c       ncorin unit number for reading coronb ce data file
c
c --- output:
c       nsigmav n<sigma-v> [sec**-1]: all processes included
c       lambda  mean free path (cm)
c       hexfrac fraction of 3rd excited state
c       ihxerr: error return
c         0     no error
c         1     radiation file coronb not found
c         2     D01BBF error in setting up ngh integration arrays
c         3     D01BBF error in setting up ngl integration arrays
c         4     SEZF   error: states "i" and    "j" are the same!
c         5     SEZF   error: state  "j" is 0
c         6     SEEF   error: states "i" and/or "j" are in error
c         7     CIEF   error
c         8     CEEF   error
c         9     EIGRF  error: eigenvalue determination is incorrect
c         10    input parameter error
c --- note: this routine uses IMSL libraries
c
c --- foreign files required (SOME INFORMATION BELOW IS ARCHAIC):
c
c       hx2:[wieland.hex]coronb.:  this file contains all the coronal
c                                  equilibrium data for the various impurities
c       dsk0:[imsl]imslibs/lib -- the IMSL replacements for the NAG codes.
c       nag:    f02aff et. al.
c       imsl:   eigrf  et. al.
c       hx2:[wieland.lib]wielib  for diagnostic matrix solvers
c                                DECOMP and SOLVEQ for ax = b
c       subroutine SECOND -- a CPU timing routine of your choice
c
c --- recommended namelist values:
c
c       kdene   = 1
c       kdeni   = 1
c       kdenz   = 1
c       ksvi    = 0
c       ksvz    = 0
c       ksve    = 0
c       krad    = 1
c       ngh     = 10
c       ngl     = 10
c       iexcit  = 1
c       ilorent = 0
c       mstate  = 4
c       ncont   = 30
c
c       these result in a cpu time of approx. 7 sec for a 25 point
c       radial profile of n*sigma(rho) with reasonably good accuracy.
c
      parameter (ms = 21, mc = 35)
      parameter (mz =  1, mi =  2)
c
      dimension     amix(mi),denix(mi),izx(mz),denzx(mz),izstrx(mz),
     .              amzx(mz)
c
c      common /cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8,
c     .              cpuii, cpuiz, cpuplf
c
      common /b0  / er0, v0, te, ti, ami(mi), deni(mi), amz(mz),
     .              denz(mz), zcor(mz), zsqcor(mz), dene,
     .              ns, nc, numi, numz, iz(mz), izstrp(mz)
      common /b1  / kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      common /b2  / nouthx,istart,ihxbug
      common /b8  / iexcit,ilorent,mstate,ncont
      common /b10 / xfrac
c
      real*8        nsigmav, lambda
c

c
      istart = istarx
c
c --- move vbls in hexnb commons
c
      iexcit  = iexcix
      ilorent = ilorenx
      mstate  = mstatx
      ncont   = ncontx
      er0     = er0x
      te      = tex
      ti      = tix
      numi    = numix
      if (numi .gt. mi)  go to 30
      do i=1,mi
        ami (i) = amix(i)
        deni(i) = denix(i)
      end do
      if (numz .gt. mz)  go to 30
      numz = numzx
      do i=1,mz
        iz(i)     = izx(i)
        izstrp(i) = izstrx(i)
        amz(i)    = amzx(i)
        denz(i)   = denzx(i)
      end do
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
      ihxbug  = 0
      nsigmav = 0.0
      lambda  = 0.0
      if (istart .eq. 0)  go to 5
c
      if (mstate+1 .gt. ms .or. ncont+1 .gt. mc)  go to 30
****  ncor = 31
      call hradin (ncor, numz, iz, izstrp, amz)
      if (ihxbug .gt. 0)  go to 20
      call hxinit
c
    5 v0 = 1.3841e6 * SQRT (er0)
      call hxradi (te, zcor, zsqcor, numz, iz, izstrp, iwatch)
      dene = 0.0
c
      do ki=1,numi
        dene = dene + deni(ki)
      end do
c
      do kz=1,numz
        dene = dene + zcor(kz) * denz(kz)
      end do
c
c --- determine calculational mode: include excitations or not?
c
      go to (21, 22),  iexcit + 1
c
c     no excitations included
c
   21 nc = 1
****  if (nouthx .gt. 0)  write (nouthx, 1000)
*1000 format (/ ' no excitations')
      go to 23
c
c     include excitations
c
   22 nc = ncont + 1
****  if (nouthx .gt. 0)  write (nouthx,1004) ncont
*1004 format (/ ' excitations, with continuum at n = ',i3)
c
   23 call lorent (v0, bperp, nc, ns)
****  if (nouthx .gt. 0)  write (nouthx, 1005) ns
*1005 format (/ ' max princ qn= ',i3)
c
      call hxsvi
      if (ihxbug .gt. 0)  go to 20
c
      call hxsve
      if (ihxbug .gt. 0)  go to 20
c
      call matri
c
      call eigen (xeig)
      if (ihxbug .gt. 0)  go to 20
c
      lambda  = xeig
      nsigmav = v0 / xeig
      hexfrac = xfrac
c
   20 ihxerr  = ihxbug
      istart  = 0
      return
c
   30 ihxbug = 10
      if (nouthx .gt. 0)  write (nouthx, 40) numi,numz,mstate,ncont
   40 format (' **** inconsistency between input vbls and upper',
     .        ' limits as defined by parameters:' /
     .          3x, 'numi= '  , i4                /
     .          3x, 'numz= '  , i4                /
     .          3x, 'mstate= ', i4                /
     .          3x, 'ncont= ' , i4)
      go to 20
c
      end

      subroutine hradin (ncor, numz, nz, izstrp, amz)
c
      USE ext_prog_info, only : nchars_12,onetwo_xsct
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     if numz ne 0, read the post radiation tables (file 'coronb').
c
c     arrangement of arad(j,k,i) for given i:
c
c    1    2    3    4    5    6    7    8    9    10   11   12   13   14
c 1  t1   t2   t3   t4   0.   t1   t2   t3   t4   0.   t1   t2   t3   t4
c 2  t2   t3   t4   t5   0.   t2   t3   t4   t5   0.   t2   t3   t4   t5
c 3 a(0) a(0) a(0) a(0)  0.  b(0) b(0) b(0) b(0)  0.  c(0) c(0) c(0) c(0
c 4  .    .    .    .    .    .    .    .    .    .    .    .    .    .
c 5  .    .    .    .    .    .    .    .    .    .    .    .    .    .
c 6  .    .    .    .    .    .    .    .    .    .    .    .    .    .
c 7  .    .    .    .    .    .    .    .    .    .    .    .    .    .
c 8 a(5) a(6) a(5) a(5)  0.  b(5) b(5) b(5) b(5)  0.  c(5) c(5) c(5) c(5
c ----------------------------------------------------------------------
c
      parameter      (mz = 1)
      common /b2/     nouthx,istart,ihxbug
      common /locrad/ arad(8,15,mz)
      dimension       amz(*), nz(*), izstrp(*), dum(8)
      character*2     spec, cdum, cmod
      character*4     loop, loop2
c

c
      if (numz .eq. 0)  return
      izsum = 0
      do i=1,numz
        izsum = izsum + izstrp(i)
      end do
      if (izsum .eq. numz)  return
      if ( numz .gt. mz  )  go to 900
c
      do imp=1,numz
        do k=1,15
          do j=1,8
            arad(j,k,imp) = 0.0
          end do
        end do
      end do
c
      call getioun(ncor,ncor)
      open (unit = ncor, status = 'OLD', err = 901, iostat = iostat,
     .      file = onetwo_xsct(1:nchars_12) // '/data/coronb')
c
      do imp=1,numz
c
c       search for desired element in data table
c
        loop = 'cont'
        do while (loop .eq. 'cont')
c
c         at first record of data set
c
          read (ncor, 1010, err=902, end=902)  spec, cmod, iz, ia
          if (iz .eq. nz(imp)) then
c
c           found correct data set
c
            loop = 'exit'
c
          else
c
c           search for first record of next data set
c
            loop2 = 'cont'
            do while (loop2 .eq. 'cont')
              read (ncor, 1020, err=902, end=902) (dum(j), j=1,8)
              read (ncor, 1015, err=902, end=902)  cdum, cmod
              if (cdum .ne. spec)  loop2 = 'exit'
            end do
            backspace (unit = ncor)
c
          end if
        end do
c
c       found desired element
c
        amz(imp) = ia
c
c       read data for this species.
c       three normal exit conditions:
c         (1)  read last species in file (read to end of file);
c         (2)  number of vectors in data set exceeds array dimension;
c         (3)  first record of next species encountered;
c
        k    = 1
        loop = 'cont'
        do while (loop .eq. 'cont')
          read (ncor, 1020, err=903, end=903) (arad(j,k,imp), j=1,8)
          read (ncor, 1015, err=903, end=120)  cdum, cmod
          k = k + 1
          if (cmod .eq. 'zb')  k = MAX0 (k,  6)
          if (cmod .eq. 'zs')  k = MAX0 (k, 11)
          if (k .gt. 15 .or. cdum .ne. spec)  loop = 'exit'
        end do
c
  120   rewind (unit = ncor)
      end do
      call giveupus(ncor)
      close (unit = ncor)
      return
c
c ----------------------------------------------------------------------
c error exits
c ----------------------------------------------------------------------
c
c     number of species exceeds array dimension
c
  900 if (nouthx .gt. 0)  write (nouthx, 9000) numz
      go to 999
c
c     file does not exist
c
  901 if (nouthx .gt. 0)  write (nouthx, 9010)  iostat
 9010 format (' ERROR opening file coronb, iostatus = ', i8)
      go to 999
c
c     could not find desired element
c
  902 if (nouthx .gt. 0)  write (nouthx, 9020)  nz(imp), ncor
 9020 format (' element #', i10, ' not found on unit ', i5)
      go to 999
c
c     error reading species data set
c
  903 if (nouthx .gt. 0)  write (nouthx, 9030)  nz(imp), ncor
c
  999 ihxbug = 1
      call giveupus(ncor)
      close (unit = ncor)
c
c ----------------------------------------------------------------------
c format statements
c ----------------------------------------------------------------------
c
 1010 format (a2, 10x, a2, 7x, 2i4)
 1015 format (a2, 10x, a2)
 1020 format (2e10.3, 3e15.6 / 3e15.6)
 9000 format (' ERROR in subroutine HRADIN' /
     .     7x, 'number of species = ', i8, 'exceeds array dimensions')
 9030 format (' ERROR reading element number', i10, ' on unit ', i5)
      return
c
      end

      real*8 function hxgf (n, y)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      g0(r) =   0.9935+0.2328/r-0.1296/(r*r)
      g1(r) = -(0.6282-0.5598/r+0.5299/(r*r))/r
      g2(r) =  (0.3887-1.1810/r+1.4700/(r*r))/(r*r)
c
      rn = FLOAT (n)
      if (n .eq. 1)  go to 11
      if (n .eq. 2)  go to 12
      hxgf = g0(rn)+g1(rn)/y+g2(rn)/(y*y)
      return
c
   11 hxgf = 1.1330-0.4059/y+.07014/(y*y)
      return
c
   12 hxgf = 1.0785-0.2319/y+.02947/(y*y)
      return
c
      end

      subroutine hxinit
c
      USE cpub_dat
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      parameter    (ms = 21, mc = 35)
c      common /cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8,
c     .              cpuii, cpuiz, cpuplf
      common /b1  / kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      common /b3  / f(ms,mc),ar(ms+1,ms+1)
      common /b4  / en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc),
     .              de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
      common /b8  / iexcit,ilorent,mstate,ncont
      dimension     be1(mc)
      data          ryd/13.6/
c
c --- total radiation rate from level n2 to level n1:
c
      rrate(n2,n1) = (8.0323e9) *
     .             (((1.0/FLOAT (n1))**2-(1.0/FLOAT (n2))**2) *
     .                   (FLOAT (n1)/FLOAT (n2)))**2*f(n1,n2)
c
      call SECOND (cpua)
c
c --- tabulate oscillator strengths
c
      call hxosc
c
      dg(1) = 1.0
      dg(2) = 1.0
      dg(3) = 3.0
      en(1) = ryd
      en(2) = ryd / 4.0
      en(3) = ryd / 4.0
c
      do i=4,mc+1
        dg(i) = (FLOAT (i-1))**2
        en(i) =  ryd / dg(i)
      end do
c
      do n1=1,mc
        an1     = FLOAT (n1)
        be1(n1) = 1.4 * LOG (an1) / an1
     .          - 0.70 / an1 - 0.51 / (an1*an1) + 1.16 / (an1**3)
     .          - 0.55 / (an1**4)
      end do
c
      do 9 n1=1,mstate
        an1  = FLOAT (n1)
        en1  = ryd / (an1 * an1)
        do 9 n2=n1+1,mc
          an2        = FLOAT (n2)
          en12       = en1 - ryd / (an2*an2)
          ae(n1,n2)  = 2.0 * ryd * f(n1,n2) / en12
          be(n1,n2)  = 4.0 * ryd * ryd * (1.0 / (en12*en12)
     .               + 4.0 * en1 / (3.0 * en12**3)
     .               + be1(n1) * en1 * en1 / (en12**4)) / (an2**3)
          de1(n1,n2) = EXP (-be(n1,n2) / ae(n1,n2)) - 0.4 * en12 / ryd
          s          = an2 - an1
          de2(n1,n2) = EXP (-be(n1,n2) / ae(n1,n2))
     .               + 0.06 * s * s / (an1 * an1 * an2)
          ge1(n1,n2) = ryd * (8.0 + 23.0 * (s/an1)**2)
     .               / (8.0 + 1.1 * an2 * s + 0.8 / (s * s)
     .               + 0.4 * (s - 1.0) * SQRT (an2 * an2 * an2 / s))
          ge2(n1,n2) = ryd * (3.0 + 11.0 * (s / an1)**2)
     .               / (6.0 + 1.6 * an2 * s + 0.3 / (s * s)
     .               + 0.8 * (s - 0.6) * SQRT (an2 * an2 * an2 / s))
    9 continue
c
c --- radiation rates
c
      do   i=1,mstate+1
        do j=1,mstate+1
          ar(i,j) = 0.0
        end do
      end do
c
      if (krad .eq. 0)  go to 1000
c
      do 210 j=4,mstate+1
      nj = j-1
      ar(j,1) = rrate(nj,1)
      do 210 i=4,j-1
      ni = i-1
      ar(j,i) = rrate(nj,ni)
  210 continue
c
c --- 2s to 1s:
c
      ar(2,1) = 0.0
c
c --- 2p to 1s:
c
      ar(3,1) = (4.0/3.0)*rrate(2,1)
c
c --- 2p to 2s:
c
      ar(3,2) = 0.0
c
      do 220 j=4,mstate+1
      nj   = j - 1
      ajsq = (FLOAT (nj))**2
c
c --- total radiation to 2s+2p
c
      tot = rrate(nj,2)
c
c --- fraction to 2s:
c
      frac2s  = 12.0*(ajsq-1.0)*(ajsq-4.0)
      frac2s  = frac2s/(frac2s+ajsq*(ajsq-4.0)+32.0*ajsq*(ajsq-1.0))
c
      ar(j,2) = frac2s*tot
      ar(j,3) = (1.0-frac2s)*tot
  220 continue
c
 1000 call SECOND (cpub)
      cpu1 = cpub - cpua
c
      return
c
      end

      subroutine hxosc
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- calculate oscillator strengths
c
      parameter (ms = 21, mc = 35)
      common /b3/f(ms,mc),ar(ms+1,ms+1)
      common /b8/iexcit,ilorent,mstate,ncont
c
      pi    = ACOS (-1.0)
      const = 32.0 / (SQRT (27.0) * pi)
c
      do 10 i=1,mstate
      ai = FLOAT (i)
      do 10 j=i+1,mc
      aj = FLOAT (j)
      y  = 1.0 - (ai/aj)**2
      f(i,j) = (const*ai/((aj*y)**3))*hxgf(i,y)
   10 continue
      return
c
      end

      subroutine hxradi (teev, z, zsq, numimp, iz, izstrp, iwatch)
c
      USE cpub_dat
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- evaluate coronal z, zsq, and radiation.
c
      parameter (mz = 1, mi = 2)
c      common /cpub  / cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8,
c     .                cpuii, cpuiz, cpuplf
      common /locrad/ arad(8,15,mz)
      dimension rad(mz),z(mz),zsq(mz),b(3),iz(mz),izstrp(mz)
      data erad /1.0/
c
      call SECOND (cpua)
      iwatch = 0
      ihigh  = 0
      tekev  = teev*1.0e-3
      do 60 imp=1,numimp
      if (izstrp(imp) .eq. 1)  go to 51
c
c --- find temperature region. test for less than 5 intervals.
c
      do 10 j=1,5
        if (arad(2,j,imp) .eq. 0.0)  go to 10
        jj = j
        if (tekev .ge. arad(1,j,imp)  .and.
     .      tekev .lt. arad(2,j,imp))  go to 30
   10 continue
c
c --- temperature out of range
c
      iwatch = 1
      if (tekev .ge. arad(1,1,imp))  go to 20
c
c --- temperature too low. return zero power
c
****  rad(imp) = 0.0
****  z  (imp) = 0.0
****  zsq(imp) = 0.0
c
c --- interpolation for low te (cdb):
c
      t1 = LOG10 (arad(1,1,imp))
      do 15 j=1,3
        jp   = 5*j - 4
        b(j) = arad(3,jp,imp) + t1*(arad(4,jp,imp) + t1*(arad(5,jp,imp)
     .       + t1*(arad(6,jp,imp) + t1*(arad(7,jp,imp)
     .       + t1*arad(8,jp,imp)))))
   15 continue
      factor   = (tekev/arad(1,1,imp))**erad
      rad(imp) = factor*(10.0**b(1))
      z(imp)   = factor*b(2)
      zsq(imp) = factor*b(3)
      go to 50
c
c --- temperature too high. compute value for tmax, then extrapolate.
c
   20 tl    = LOG10 (arad(2,jj,imp))
      ihigh = 1
   30 if (ihigh .eq. 0)
     .tl    = LOG10 (tekev)
c
c --- polynomial fits  bb = sum(a(ik)*log(te)**k) etc.
c
      bb = 0.0
      cc = 0.0
      dd = 0.0
c
      do kk=1,6
        k  = 7-kk
        bb = bb*tl + arad((k+2),jj,imp)
        cc = cc*tl + arad((k+2),(jj+5),imp)
        dd = dd*tl + arad((k+2),(jj+10),imp)
      end do
c
      rad(imp) = 0.0
      z(imp)   = cc
      zsq(imp) = dd
      if (bb .lt. -38.0)  go to 50
      rad(imp) = 10.0**bb
c
c --- for te > tmax extrapolate as for pure Bremsstrahlung
c
      if (ihigh .eq. 1)  rad(imp) =
     .                   rad(imp) * SQRT (tekev / arad(2,4,imp))
      ihigh = 0
   50 continue
      go to 60
c
c --- fully stripped option [ izstrp(imp) = 1 ]
c
   51 z(imp)   = iz(imp)
      zsq(imp) = iz(imp)**2
   60 continue
c
      call SECOND (cpub)
      cpu2 = cpu2 + cpub - cpua
c
      return
c
      end

      subroutine hxsve

      USE cpub_dat
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- calculate the rate coefficients for collisions with electrons.
c
      parameter   (ms = 21, mc = 35, mz = 1, mi = 2)
c      common /cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8,
c     .             cpuii, cpuiz, cpuplf
      common /b0  / er0, v0, te, ti, ami(mi), deni(mi), amz(mz),
     .              denz(mz), zcor(mz), zsqcor(mz), dene,
     .              ns, nc, numi, numz, iz(mz), izstrp(mz)
      common /b1  / kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      common /b2  / nouthx,istart,ihxbug
      common /b4  / en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc),
     .              de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
      common /b6  / cii(ms+1,mi),cei(ms+1,ms+2,mi),ciz(ms+1,mz),
     .              cez(ms+1,ms+2,mz),cie(ms+1),cee(ms+1,ms+2),
     .              ccxi(ms+1,mi),ccxz(ms+1,mz)
      dimension et(mc+1)
      dimension xgl(64),wgl(64)
      character*2 typent
****  external d01bax
c
      sinh1(arg) = SINH (arg)/arg
c
      call SECOND (cpua)
c
c --- initialize:
c
      if (istart .eq. 0)  go to 5
      pi   = ACOS (-1.0)
      sqpi = SQRT ( pi )
c
c     e-impact thresholds (temp.)
c
      et(1) = (13.723/13.6)*en(1)
      do i=2,mc+1
        et(i) = en(i)
      end do
c
      ifail  =  0
      typent = 'gl'
      call d01bbf (typent, ngl, wgl, xgl, ifail)
****  call d01bbf (d01bax, 0.0, 1.0, 0, ngl, wgl, xgl, ifail)
      if ( ifail .eq. 0) go to 5
      if (nouthx .gt. 0) then
        write  (nouthx, 3939) ifail
 3939   format (' ERROR in HXSVE calling D01BBF: ifail= ', i4)
      end if
      ihxbug = 3
      return
c
    5 nsp1 = ns+1
c
      do 205 i=1,ns
      cie(i) = 0.0
      do 205 j=i+1,nsp1
      cee(i,j) = 0.0
  205 continue
c
      ve = 5.931e7 * SQRT (te)
      if (kdene .eq. 0)  go to 1000
      go to (220, 210)  ksve + 1
c
c --- Gauss-Laguerre integrations:
c
  210 facte = (2.0/sqpi)*ve * EXP (-v0*v0/(ve*ve))
      do 211 kk=1,ngl
      do 212 i=1,ns
      cie(i) = cie(i)
     .  +wgl(kk)*(xgl(kk)+et(i)/te)
     .  *sinh1(2.0*v0 * SQRT (xgl(kk)+et(i)/te)/ve)
     .  *sief(i,et(i)+te*xgl(kk))
     .  *facte * EXP (-et(i)/te)
      do 213 j=i+1,ns
      cee(i,j) = cee(i,j)
     .         + wgl(kk)*(xgl(kk)+(et(i)-et(j))/te)
     .          *sinh1(2.0*v0 * SQRT (xgl(kk)+(et(i)-et(j))/te)/ve)
     .          *seef(i,j,et(i)-et(j)+te*xgl(kk))
     .          *facte * EXP (-(et(i)-et(j))/te)
  213 continue
      do 214 j=nsp1,nc
      cee(i,nsp1) = cee(i,nsp1)
     .  +wgl(kk)*(xgl(kk)+(et(i)-et(j))/te)
     .  *sinh1(2.0*v0 * SQRT (xgl(kk)+(et(i)-et(j))/te)/ve)
     .  *seef(i,j,et(i)-et(j)+te*xgl(kk))
     .  *facte * EXP (-(et(i)-et(j))/te)
  214 continue
  212 continue
  211 continue
      go to 1000
c
c --- Gauss-Laguerre averages for ionization of 1s, ionization of 2s,
c     and excitations among 1s, 2s, and 2p:
c
  220 facte = (2.0/sqpi)*ve * EXP (-v0*v0/(ve*ve))
      do 221 kk=1,ngl
      do 222 i=1,MIN0 (2,ns)
      cie(i) = cie(i)
     .  +wgl(kk)*(xgl(kk)+et(i)/te)
     .  *sinh1(2.0*v0 * SQRT (xgl(kk)+et(i)/te)/ve)
     .  *sief(i,et(i)+te*xgl(kk))
     .  *facte * EXP (-et(i)/te)
      do 222 j=i+1,MIN0 (3,ns)
      cee(i,j) = cee(i,j)
     .  +wgl(kk)*(xgl(kk)+(et(i)-et(j))/te)
     .  *sinh1(2.0*v0 * SQRT (xgl(kk)+(et(i)-et(j))/te)/ve)
     .  *seef(i,j,et(i)-et(j)+te*xgl(kk))
     .  *facte * EXP (-(et(i)-et(j))/te)
  222 continue
  221 continue
c
c --- rates from Vriens & Smeets for other reactions:
c     ionizations--
c
      do 223 i=3,ns
      cie(i) = cief(i,te)
  223 continue
c
c     excitations--
c
      do 224 i=1,ns
        do 225 j=i+1,ns
          if (j .le. 3)  go to 225
          cee(i,j) = ceef(i,j,te)
  225   continue
        do 226 j=nsp1,nc
          cee(i,nsp1) = cee(i,nsp1)+ceef(i,j,te)
  226   continue
  224 continue
c
 1000 call SECOND (cpub)
      cpu5 = cpu5 + cpub - cpua
c
      return
c
      end

      subroutine hxsvi
c
      USE cpub_dat
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- calculate the rate coefficients for collisions with ions and impurities
c
      parameter    (ms = 21, mc = 35, mz = 1, mi = 2)
c      common /cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8,
c     .              cpuii, cpuiz, cpuplf
      common /b0  / er0, v0, te, ti, ami(mi), deni(mi), amz(mz),
     .              denz(mz), zcor(mz), zsqcor(mz), dene,
     .              ns, nc, numi, numz, iz(mz), izstrp(mz)
      common /b1  / kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      common /b2  / nouthx,istart,ihxbug
      common /b6  / cii(ms+1,mi),cei(ms+1,ms+2,mi),ciz(ms+1,mz),
     .              cez(ms+1,ms+2,mz),cie(ms+1),cee(ms+1,ms+2),
     .              ccxi(ms+1,mi),ccxz(ms+1,mz)
      dimension     vi(mi),vz(mz)
      dimension     xgh(64),wgh(64),w1(32)
      character*2   typent
****  external      d01baw
      data          ister0 /1/, zz1 /1.0/
c
c --- initialize:
c
      call SECOND (cpua)
c
      if (istart .eq. 0)  go to 5
      pi   = ACOS (-1.0) 
      sqpi = SQRT (pi)
c
      ngh1 = (ngh+1)/2
      do i=1,ngh1
        w1(i) = 2.0
      end do
      if (2*ngh1-ngh .eq. 1)  w1(ngh1) = 1.0
      ifail  = 0
      typent = 'gh'
      call d01bbf (typent, ngh, wgh, xgh, ifail)
****  call d01bbf (d01baw, 0.0, 1.0, 0, ngh, wgh, xgh, ifail)
      if (ifail  .eq. 0) go to 5
      if (nouthx .gt. 0) then
        write (nouthx, 3939) ifail
 3939   format (' ERROR in HXSVI calling D01BBF: ifail = ', i4)
      end if
      ihxbug = 2
      return
c
    5 continue
      nsp1 = ns+1
      if (er0 .ne. er0old) ister0 = 1
c
c --- ion reactions-----------------------------------------------
c
      call SECOND (cpua1)
c
c --- first zero out arrays:
c
      if ( (ister0 .eq. 1 .or. istart .eq. 1) .or. ksvi .eq. 1) then
        do 10 ki=1,numi
        do 10 i=1,ns
        cii(i,ki) = 0.0
        ccxi(i,ki) = 0.0
        do 10 j=i+1,nsp1
        cei(i,j,ki) = 0.0
   10   continue
      end if
c
      do ki=1,numi
        vi(ki) = 1.3841e6 * SQRT (ti/ami(ki))
      end do
c
      if (kdeni .eq. 0)  go to 100
      if ( ksvi .eq. 0)  go to 40
c
c --- gauss-hermite integrations:
c
      do 20 kk=1,ngh1
      do 20 ki=1,numi
      xp = (xgh(kk)+v0/vi(ki))**2
      xm = (xgh(kk)-v0/vi(ki))**2
      s = 1.0
      if (xgh(kk)-v0/vi(ki) .lt. 0.0)  s = -1.0
      ei = ti/ami(ki)
      do 20 i=1,ns
      cii(i,ki) = cii(i,ki)+wgh(kk)*w1(kk)
     .  *(xp*siif(i,ei*xp)-s*xm*siif(i,ei*xm))
      ccxi(i,ki) = ccxi(i,ki)+wgh(kk)*w1(kk)
     .  *(xp*scxif(i,ei*xp)-s*xm*scxif(i,ei*xm))
      do j=i+1,ns
        cei(i,j,ki) = cei(i,j,ki)+wgh(kk)*w1(kk)
     .              *(xp*sezf(i,j,zz1,ei*xp)-s*xm*sezf(i,j,zz1,ei*xm))
      end do
      do j=nsp1,nc
        cei(i,nsp1,ki) = cei(i,nsp1,ki)+wgh(kk)*w1(kk)
     .                *(xp*sezf(i,j,zz1,ei*xp)-s*xm*sezf(i,j,zz1,ei*xm))
      end do
   20 continue
c
      do 30 ki=1,numi
      facti = vi(ki)*vi(ki)/(2.0*sqpi*v0)
      do 30 i=1,ns
      cii(i,ki) = facti*cii(i,ki)
      ccxi(i,ki) = facti*ccxi(i,ki)
      do 30 j=i+1,nsp1
      cei(i,j,ki) = facti*cei(i,j,ki)
   30 continue
      go to 100
c
c --- simple multiplication instead of maxwellian averaging:
c ---      cii and cei are independent of plasma parameters in this
c ---      approximation       so as long as er0 doesn't change, use
c ---      the old values.
c ---      ister0 = 1 ==> restart everything
c ---             = 0 ==> cruise in no-update mode
c
c --- redo everything if er0 changes
c
   40 if ((ister0 .eq. 1 .or. istart .eq. 1) .or. ksvi .eq. 1) then
      do 41 ki=1,numi
      do 41 i=1,ns
      cii (i,ki) =  siif(i,er0)*v0
      ccxi(i,ki) = scxif(i,er0)*v0
      do 42 j=i+1,ns
      cei(i,j,ki) = sezf(i,j,zz1,er0)*v0
   42 continue
      do j=nsp1,nc
        cei(i,nsp1,ki) = cei(i,nsp1,ki)+sezf(i,j,zz1,er0)*v0
      end do
   41 continue
      end if
  100 continue
      call SECOND (cpua2)
      cpuii = cpuii + cpua2 - cpua1
c
c --- impurity reactions------------------------------------------
c
c --- zero out arrays:
c
      do 110 kz=1,numz
      vz(kz) = 1.3841e6 * SQRT (ti/amz(kz))
      do 110 i=1,ns
      ciz(i,kz) = 0.0
      ccxz(i,kz) = 0.0
      do 110 j=i+1,nsp1
      cez(i,j,kz) = 0.0
  110 continue
c
      if (kdenz .eq. 0)  go to 1000
      if (ksvz  .eq. 0)  go to  140
c
c --- gauss-hermite integrations:
c
      do 120 kk=1,ngh1
      do 120 kz=1,numz
      xp = (xgh(kk)+v0/vz(kz))**2
      xm = (xgh(kk)-v0/vz(kz))**2
      s = 1.0
      if (xgh(kk)-v0/vz(kz) .lt. 0.0)  s = -1.0
      ez = ti/amz(kz)
      zz = zcor(kz)
      do 120 i=1,ns
      ciz(i,kz) = ciz(i,kz)+wgh(kk)*w1(kk)
     .  *(xp*sizf(i,zz,ez*xp)-s*xm*sizf(i,zz,ez*xm))
      ccxz(i,kz) = ccxz(i,kz)+wgh(kk)*w1(kk)
     .  *(xp*scxzf(i,zz,ez*xp)-s*xm*scxzf(i,zz,ez*xm))
      do j=i+1,ns
        cez(i,j,kz) = cez(i,j,kz)+wgh(kk)*w1(kk)
     .              *(xp*sezf(i,j,zz,ez*xp)-s*xm*sezf(i,j,zz,ez*xm))
      end do
      do j=nsp1,nc
        cez(i,nsp1,kz) = cez(i,nsp1,kz)+wgh(kk)*w1(kk)
     .                 *(xp*sezf(i,j,zz,ez*xp)-s*xm*sezf(i,j,zz,ez*xm))
      end do
  120 continue
c
      do kz=1,numz
        factz = vz(kz)**2/(2.0*sqpi*v0)
        do i=1,ns
          ciz(i,kz) = factz*ciz(i,kz)
          ccxz(i,kz) = factz*ccxz(i,kz)
          do j=i+1,nsp1
            cez(i,j,kz) = factz*cez(i,j,kz)
          end do
        end do
      end do
c
      go to 1000
c
c --- simple multiplications:
c
  140 do kz=1,numz
        do i=1,ns
          ciz (i,kz) = sizf(i,zcor(kz),er0)*v0
          ccxz(i,kz) = scxzf(i,zcor(kz),er0)*v0
          do j=i+1,ns
            cez(i,j,kz) = sezf(i,j,zcor(kz),er0)*v0
          end do
          do j=nsp1,nc
            cez(i,nsp1,kz) = cez(i,nsp1,kz)+sezf(i,j,zcor(kz),er0)*v0
          end do
        end do
      end do
c
 1000 call SECOND (cpub)
      cpu4   = cpu4  + cpub - cpua
      cpuiz  = cpuiz + cpub - cpua2
      ister0 = 0
      er0old = er0
      return
c
      end

      integer function imax (f, mfm1)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c This function computes the index of the maximum element of the vector
c   f excluding the first and last elements.
c ----------------------------------------------------------------------
c
      dimension f(*)
c
      fmax = f(2)
      imax = 2
      do i=3,mfm1
        if (f(i) .ge. fmax) then
          fmax = f(i)
          imax = i
        end if
      end do
      return
c
      end

      subroutine initfit
c
c ----------------------------------------------------------------------
c --- New version, with different fits valid from lower beam energies
c --- Received from C. Boley, LLNL
c --- Modified for use by the SuperCode, John Mandrekas, 03/27/92
c
      implicit none
c
      integer mz, mt1, mt2, mt3
      parameter(mz=4, mt1=4, mt2=3, mt3=2)
      integer i
c
      integer       nth1, nth2, nth3, ntz1(mz), ntz2(mz), ntz3(mz)
      real*8        A1(mt1,mt2,mt3), A2(mt1,mt2,mt3,mz)
      common /cfit/ A1, A2, nth1, nth2, nth3, ntz1, ntz2, ntz3
c
c..   pure plasma
c
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
c
      do i = 1, 4
         ntz1(i) = 4
         ntz2(i) = 2
         ntz3(i) = 2
      end do
c
c..   He
c
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
c
c..   C
c
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
c
c..   O
c
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
c
c..   Fe
c
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
c
      return
c
      end

      subroutine inject (atw, codeid, debin,drutpi,droti,dri,ds1,dzi,
     .                   elongi,ib,ie,kbe,ke,kz,ki,mfm1,mim1,mjm1,nebin,
     .                   newpar,nout,psiax,psi,r,rmajor,rin,rmax,sgxn,
     .                   sgxnloc,sgxnmi,x0,y0,z0,vx0,vy0,vz0,vbeam,z,
     .                   zangrot,zax,zmin,zmax,izone,pzone,rzone,rpos,
     .                   xpos,ypos,zpos,myid,tenter,smax,texit)
c
      USE constnts,only : kevperg,xmassp
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  this subroutine follows the particle from the pivot point into,
c     through, or around the plasma.  note:  bilinear interpolation is
c     used to obtain psi as a function of track length along a
c     collisionless neutral trajectory (rotating discharges only).
c     bicubic spline interpolation was tested, but found to provide no
c     appreciable increase in accuracy.
c ----------------------------------------------------------------------
c
      external  RANDOM12                         ! random number generator
      dimension cvec(200)
      data      cvec
     .     /  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0,
     .       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
     .       21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0,
     .       31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 40.0,
     .       41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0,
     .       51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0, 60.0,
     .       61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0,
     .       71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0,
     .       81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0, 90.0,
     .       91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0,100.0,
     .      101.0,102.0,103.0,104.0,105.0,106.0,107.0,108.0,109.0,110.0,
     .      111.0,112.0,113.0,114.0,115.0,116.0,117.0,118.0,119.0,120.0,
     .      121.0,122.0,123.0,124.0,125.0,126.0,127.0,128.0,129.0,130.0,
     .      131.0,132.0,133.0,134.0,135.0,136.0,137.0,138.0,139.0,140.0,
     .      141.0,142.0,143.0,144.0,145.0,146.0,147.0,148.0,149.0,150.0,
     .      151.0,152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0,160.0,
     .      161.0,162.0,163.0,164.0,165.0,166.0,167.0,168.0,169.0,170.0,
     .      171.0,172.0,173.0,174.0,175.0,176.0,177.0,178.0,179.0,180.0,
     .      181.0,182.0,183.0,184.0,185.0,186.0,187.0,188.0,189.0,190.0,
     .      191.0,192.0,193.0,194.0,195.0,196.0,197.0,198.0,199.0,200.0/
c
      parameter (ktk = 100)
      character  codeid*8
      dimension  psi(ki,*), sgxn(4,kz,6,*), sgxnloc(*), sgxnmi(ke,*),
     .           r(*), vbeam(ke,*), z(*), zangrot(*), e1(ktk), i1(ktk),
     .           sgxntab(ktk)
      data       seed0    /0.0/
c


C     the following assumes that data from previous
c     particle is saved. this ok only if passed through
c     argument list. HSJ

      if (newpar .eq. 0)  go to 100
c
c calculate times for particle to enter and exit toroidal box surrounding plasma
c
      call timtor (rin,rmax,x0,y0,z0,vx0,vy0,vz0,zmin,zmax,tenter,texit)
      if (tenter .le. -1.0e10)  go to 140
c
c advance particle to edge of box
c
      x0 = x0 + vx0*tenter
      y0 = y0 + vy0*tenter
      z0 = z0 + vz0*tenter
c

      if (nebin .ne. 0) then
c
c ----------------------------------------------------------------------
c follow collisionless neutral trajectory to obtain minimum mean free path.
c this is required to account for toroidal rotation.
c
c ALSO:
c    Only include those energy groups that fall in a valid flux zone,
c    (izone1 .le. mfm1) for the search in subroutine GETSGXN. The array
c    of energy values e1 is now ne1 elements long rather than n1.
c                                            Daniel Finkenthal    9-6-95
c ----------------------------------------------------------------------
c
      smax = vbeam(ie,ib) * (texit-tenter)
      n1   = 2.0 + smax/ds1
      if (n1 .gt. ktk) then
        n1  = ktk
        ds1 = smax / FLOAT (n1-1)
        write (nout, 200) ds1
      end if
      dt1 = smax / (FLOAT (n1-1) * vbeam(ie,ib))
      ne1 = 0
      if (codeid .eq. 'onedee') then
        do i=1,n1
          delt   = (cvec(i) - 1.0) * dt1
          x1     = x0 + delt*vx0
          y1     = y0 + delt*vy0
          z1     = z0 + delt*vz0
          r1     = SQRT (x1**2 + y1**2)
          ir     = (r1-r(1))*dri + 1.0
          iz     = (z1-z(1))*dzi + 1.0
          ir     = MIN0 (ir,mim1)
          iz     = MIN0 (iz,mjm1)
          if (ir .le. 0 .or. iz .le. 0)go to 10  !HSJ 1/28/2000
          p1     = SQRT ((r1-rmajor)**2+(elongi*(z1-zax))**2)
          izone1 = p1*droti + 1.0
          if (izone1 .le. mfm1) then
            ne1     = ne1 + 1
            i1(ne1) = izone1
            usq     = ( x1**2 + y1**2 ) * zangrot(i1(ne1))**2
            vdotu   = (x1*vy0 - y1*vx0) * zangrot(i1(ne1))
            vrel1   = vbeam(ie,ib)**2 + usq - 2.0*vdotu
            e1(ne1) = 0.5 * xmassp * vrel1 * kevperg
            e1(ne1) = ABS (e1(ne1))
          end if
        end do
      else
        do i=1,n1
          delt  = (cvec(i) - 1.0) * dt1
          x1    = x0 + delt*vx0
          y1    = y0 + delt*vy0
          z1    = z0 + delt*vz0
          r1    = SQRT (x1**2+y1**2)
          ir    = (r1-r(1))*dri + 1.0
          iz    = (z1-z(1))*dzi + 1.0
          ir    = MIN0 (ir,mim1)
          iz    = MIN0 (iz,mjm1)
          if (ir .le. 0 .or. iz .le. 0)go to 10  !HSJ 1/28/2000
          area1 = (r1-r(ir))*(z1-z(iz))
          area2 = (r(ir+1)-r1)*(z1-z(iz))
          area3 = (r(ir+1)-r1)*(z(iz+1)-z1)
          area4 = (r1-r(ir))*(z(iz+1)-z1)
          p1    = (area3*psi(ir,iz)   + area4*psi(ir+1,iz)
     .          +  area2*psi(ir,iz+1) + area1*psi(ir+1,iz+1))*dri*dzi
          p1     =  MAX (p1,psiax)
          izone1 = SQRT (p1-psiax)*drutpi + 1.0
          if (izone1 .le. mfm1) then
            ne1     = ne1 + 1
            i1(ne1) = izone1
            usq     = (x1**2+y1**2)*zangrot(i1(ne1))**2
            vdotu   = (x1*vy0-y1*vx0)*zangrot(i1(ne1))
            vrel1   = vbeam(ie,ib)**2 + usq - 2.0*vdotu
            e1(ne1) = 0.5 * xmassp * vrel1 * kevperg
            e1(ne1) = ABS (e1(ne1))
          end if
        end do
      end if
 10   call getsgxn (e1, i1, ne1, ktk, ib, ie, sgxn, nebin, debin, kbe,
     .              kz, nout, sgxntab)
      xnorm         = amaxaf (sgxntab, 1, ne1)
      sgxnmi(ie,ib) = 1.0 / xnorm
      end if
c
c ----------------------------------------------------------------------
c inject neutral into plasma
c ----------------------------------------------------------------------
c
c  set coordinates and time for entering box
c
  100 xpos      = x0
      ypos      = y0
      zpos      = z0
      tt        = tenter
      izone     = mfm1 + 1     ! initially neutral is outside the plasma
      smin_step = 0.1                           ! 0.1 cm min step or
      smin_step = MIN (smin_step, smax/1000.0)  ! make scale-independent
      smin_time = smin_step / SQRT (vx0**2 + vy0**2 + vz0**2)
c
c  follow particle into plasma
c
* 110 dfac  = -LOG (RANF   (     ))

  110 dfac  = -LOG (RANDOM12 (seed0))
c
c     if neutral is not yet in the plasma (izone ge mf) then take steps
c     of minimum size 1 mm until we enter the plasma. we could find the
c     exact plasma boundary but that would be overkill. with 1mm step
c     size we certainly are within any physics scales we could resolve
c     near the plasma edge. this avoids wasting a lot of steps outside
c     the plasma and will be a significant savings if the cross sections
c     are large. dfac is 1/exponentially distributed, so for large steps
c     don't modify tstep. --------------------------- 27 Oct 95 ---- HSJ
c
      if (izone .le. mfm1) then                     !  inside the plasma
        tstep =      dfac * sgxnmi(ie,ib) / vbeam(ie,ib)
      else
        tstep = MAX (dfac * sgxnmi(ie,ib) / vbeam(ie,ib),
     .               smin_time)                     ! outside the plasma
      end if
c
      tt    = tt  + tstep
      if (tt .ge. texit)  go to 140
      xpos = xpos + vx0*tstep
      ypos = ypos + vy0*tstep
      zpos = zpos + vz0*tstep
      rpos = SQRT (xpos**2 + ypos**2)
c
c  determine zone in which particle collides for 'onedee' geometry
c
      if (codeid .eq. 'onedee') then
        rzone2 = (rpos-rmajor)**2 + (elongi*(zpos-zax))**2
        rzone  = SQRT (rzone2)
        izone  = rzone*droti + 1.0
      else
c
c  determine zone in which particle collides for general geometry;
c     use bilinear interpolation away from magnetic axis,
c     and biquadratic interpolation near the axis.
c
        i     = (rpos-r(1))*dri+1.0
        j     = (zpos-z(1))*dzi+1.0
        i     = MIN0 (i,mim1)
        j     = MIN0 (j,mjm1)
        psix  = MIN  (psi(i,j),psi(i+1,j),psi(i,j+1),psi(i+1,j+1))
        ptest = (psix-psiax)*(drutpi/mfm1)**2
        if (ptest .ge. 0.02) then
          area1 = (rpos-r(i))*(zpos-z(j))
          area2 = (r(i+1)-rpos)*(zpos-z(j))
          area3 = (r(i+1)-rpos)*(z(j+1)-zpos)
          area4 = (rpos-r(i))*(z(j+1)-zpos)
          pzone = (area3*psi(i,j) + area4*psi(i+1,j)
     .          + area1*psi(i+1,j+1) + area2*psi(i,j+1))*dri*dzi
        else
          call pfit(psi(i-1,j-1),r(i-1),z(j-1),rpos,zpos,ki,pzone,dum,
     .              dum)
        end if
        pzone =  MAX (pzone,psiax)
        izone = SQRT (pzone-psiax)*drutpi + 1.0
      end if
c
c     if particle has psuedo-collision, continue following particle.
c     if particle has real collision, return.
c
      if (izone .gt. mfm1)  go to 110 ! the particle is inside the box..
c                                     ..but still outside the plasma
      if (nebin .ne. 0   ) then
        usq   = (rpos*zangrot(izone))**2
        vdotu = (xpos*vy0-ypos*vx0)*zangrot(izone)
        vrel2 = vbeam(ie,ib)**2 + usq - 2.0*vdotu
        eova  = 0.5 * xmassp*vrel2*kevperg
        eova  = ABS (eova)
        call getsgxn (eova, izone, 1, ktk, ib, ie, sgxn, nebin,
     .                debin, kbe, kz, nout, sgxnloc)
      else
        index      = ke*(ib-1) + ie
        sgxnloc(1) = sgxn(1,izone,index,1)
        sgxnloc(2) = sgxn(2,izone,index,1)
        sgxnloc(3) = sgxn(3,izone,index,1)
        sgxnloc(4) = sgxn(4,izone,index,1)
      end if
*     if (RANF   (     ) .gt. sgxnloc(4) * sgxnmi(ie,ib))  go to 110
      if (RANDOM12 (seed0) .gt. sgxnloc(4) * sgxnmi(ie,ib))  go to 110
      return
c
c  set flag to indicate that particle hit wall
c
  140 izone = mfm1 + 1
  200 format (' WARNING from subroutine INJECT:'                  /
     .        '         maximum number of grid elements exceeded' /
     .        '         increasing ds1 to ', e10.3, ' cm')
      return
c
      end

      subroutine logint (x, y)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c interpolates y(x) quadratically and logarithmically
c ----------------------------------------------------------------------
c
      dimension  xdat(15), ydat(15)
c
      data xdat/4.0,6.0,8.0,10.0,20.0,30.0,40.0,60.0,80.0,100.0,
     .          200.0,300.0,400.0,600.0,800.0/
      data ydat/8.95e-01,8.75e-01,8.70e-01,8.65e-01,8.20e-01,
     .          7.25e-01,6.25e-01,4.40e-01,2.90e-01,1.90e-01,2.40e-02,
     .          5.25e-03,1.20e-03,1.60e-04,5.40e-05/
c
      mdat  = 15
      mdatm = mdat - 1
c
      do i0=2,mdatm
        if (xdat(i0) .ge. x)  go to 11
      end do
c
   11 if (i0 .gt. mdatm) i0 = mdatm
      im = i0-1
      ip = i0+1
      ylogp = LOG (ydat(ip))
      ylog0 = LOG (ydat(i0))
      ylogm = LOG (ydat(im))
      dm = x-xdat(im)
      d0 = x-xdat(i0)
      dp = x-xdat(ip)
      d0m = xdat(i0)-xdat(im)
      dp0 = xdat(ip)-xdat(i0)
      dpm = xdat(ip)-xdat(im)
      dm0 = -d0m
      d0p = -dp0
      dmp = -dpm
      facm = d0*dp/(dm0*dmp)
      fac0 = dm*dp/(d0m*d0p)
      facp = dm*d0/(dpm*dp0)
      ylog = facm*ylogm+fac0*ylog0+facp*ylogp
      y = EXP (ylog)
c
      return
c
      end

      subroutine lorent (v0, bperp, nc, ns)
c
      USE cpub_dat
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      parameter    (ms = 21)
c      common /cpub/ cpu1, cpu2, cpu3, cpu4, cpu5, cpu6, cpu7, cpu8,
c     .              cpuii, cpuiz, cpuplf
      common /b5  / al(ms+1)
      common /b8  / iexcit, ilorent, mstate, ncont
      data          w0 /4.1341e16/, almin /1.0e-10/, almax /1.0e15/,
     .              expl /4.0/
      dimension     sl(30)
c
      data sl/
     .   1.00000000e+00, 2.50000000e-01, 2.77777778e-02, 1.73611111e-03,
     .   6.94444444e-05, 1.92901235e-06, 3.93675989e-08, 6.15118733e-10,
     .   7.59405843e-12, 7.59405843e-14, 6.27608135e-16, 4.35838982e-18,
     .   2.57892889e-20, 1.31578005e-22, 5.84791131e-25, 2.28434036e-27,
     .   7.90429189e-30, 2.43959626e-32, 6.75788439e-35, 1.68947110e-37,
     .   10*0.0/
c
      fl (n, an3, el) =
     .   ((4.0 / (an3*el))**(2*n-1)) * EXP (-2.0 / (3.0*an3*el)) / an3
c
c     input:
c       v0
c       bperp
c       nc
c     output:
c       ns
c       al(i),i = 1,ns
c
      call SECOND (cpua)
c
      one   = 1.0
      eight = 8.0
c
      if (nc .ge. 2)  go to 1
      ns    = 1
      al(1) = 0.0
      go to 30
c
    1 if (ilorent .eq. 0) then
        ns = mstate + 1
        do il=1,ns
          al(il) = 0.0
        end do
        go to 30
      end if
c
      if (bperp .eq. 0.0)  go to 3
c
c --- ep = ABS (v x b)/(electric field at first Bohr radius)
c
      ep    = v0*bperp / 5.1417e13
      ncrit = 0.5 / (ep**(1.0/expl))
      if (ncrit .le. 1)  ns = 1
      if (ncrit .ge. 2)  ns = ncrit + 1
      if (ns .le. (mstate+1))  go to 4
c
    3 ns = mstate + 1
c
    4 do 10 i=1,ns
        if (i .ge. 4)  go to 13
        go to (11, 12, 12), i
   11   al(1) = fl (1, one  , ep)
        go to 10
   12   al(i) = fl (2, eight, ep) * 0.5
        go to 10
   13   ni = i-1
        ani3  = (FLOAT (ni))**3
        al(i) = fl (ni, ani3, ep) * sl(ni)
   10 continue
c
c --- normalize rates; zero out those which are too small
c
      do i=1,ns
        al(i) = al(i)*w0
        if (al(i) .le. almin)
     .  al(i) = 0.0
      end do
c
c --- if some rates are too high, truncate system
c
      do i=1,ns
        i1 = i
        if (al(i) .gt. almax)  go to 21
      end do
      go to 30
c
   21 ns = i1 - 1
c
   30 call SECOND (cpub)
      cpu3 = cpu3 + cpub - cpua
      return
c
      end

      subroutine matri
c
      USE cpub_dat
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- form the matrix of rate coefficients which occurs
c     in the system of rate equations.
c     evaluate deexcitation via detailed balance.
c --- vectorized
c
      parameter    (ms = 21, mc = 35, mz = 1, mi = 2)
c      common /cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8,
c     .              cpuii, cpuiz, cpuplf
      common /b0  / er0, v0, te, ti, ami(mi), deni(mi), amz(mz),
     .              denz(mz), zcor(mz), zsqcor(mz), dene,
     .              ns, nc, numi, numz, iz(mz), izstrp(mz)
      common /b3  / f(ms,mc),ar(ms+1,ms+1)
      common /b4  / en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc),
     .              de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
      common /b5  / al(ms+1)
      common /b6  / cii(ms+1,mi),cei(ms+1,ms+2,mi),ciz(ms+1,mz),
     .              cez(ms+1,ms+2,mz),cie(ms+1),cee(ms+1,ms+2),
     .              ccxi(ms+1,mi),ccxz(ms+1,mz)
      common /b7  / q(ms+1,ms+1)
      dimension     dge(ms+1),dgi(ms+1)
c
cc*** dump c arrays
c      write (60,3310) cii
c3310  format (/ ' cii= ' / (1x,1p5e11.3))
c      write (60,3315) ccxii
c3315  format (/ ' ccxii= ' / (1x,1p5e11.3))
c      write (60,3320) cei
c3320  format (/ ' cei= ' / (1x,1p5e11.3))
c      write (60,3330) ciz
c3330  format (/ ' ciz= ' / (1x,1p5e11.3))
c      write (60,3335) ccxz
c3335  format (/ ' ccxz= ' / (1x,1p5e11.3))
c      write (60,3340) cez
c3340  format (/ ' cez= ' / (1x,1p5e11.3))
c      write (60,3350) cie
c3350  format (/ ' cie= ' / (1x,1p5e11.3))
c      write (60,3360) cee
c3360  format (/ ' cee= ' / (1x,1p5e11.3))
c
      call SECOND (cpua)
c
      nsp1 = ns+1
      do 10 i=1,ns
      do 10 j=1,ns
        q(i,j) = 0.0
   10 continue
c
c --- detailed-balance factors (note that en's are positive):
c
      do i=1,ns
        dge(i) = dg(i) * EXP (en(i)/te)
        dgi(i) = dg(i) * EXP (en(i)/ti)
      end do
c
c --- rates due to collisions with electrons
c     radiation rates       lorentz rates:
c
      do 20 j=1,ns
      q(j,j) = -cie(j)*dene -al(j)
      do 21 jp=j+1,nsp1
      q(j,j) = q(j,j) -cee(j,jp)*dene
   21 continue
      do 22 jp=1,j-1
      q(j,j) = q(j,j)
     .  -(dge(jp)/dge(j))*cee(jp,j)*dene
     .  -ar(j,jp)
   22 continue
      do 23 i=j+1,ns
      q(i,j) = q(i,j)
     .  +cee(j,i)*dene
   23 continue
      do 24 i=1,j-1
      q(i,j) = q(i,j)
     .  +(dge(i)/dge(j))*cee(i,j)*dene
     .  +ar(j,i)
   24 continue
   20 continue
c
c --- add rates due to collisions with ions:
c
      do 30 ki=1,numi
      do 30 j=1,ns
      q(j,j) = q(j,j)
     .  -cii(j,ki)*deni(ki)
      do 31 jp=j+1,nsp1
      q(j,j) = q(j,j) - cei(j,jp,ki)*deni(ki)
   31 continue
      do 32 jp=1,j-1
      q(j,j) = q(j,j) - (dgi(jp)/dgi(j))*cei(jp,j,ki)*deni(ki)
   32 continue
      do 33 i=j+1,ns
      q(i,j) = q(i,j) + cei(j,i,ki)*deni(ki)
   33 continue
      do 34 i=1,j-1
      q(i,j) = q(i,j) + (dgi(i)/dgi(j))*cei(i,j,ki)*deni(ki)
   34 continue
   30 continue
c
c --- add rates due to collisions with impurities:
c
      do 40 kz=1,numz
      do 40 j=1,ns
      q(j,j) = q(j,j)
     .  -ciz(j,kz)*denz(kz)
      do 41 jp=j+1,nsp1
      q(j,j) = q(j,j)
     .  -cez(j,jp,kz)*denz(kz)
   41 continue
      do 42 jp=1,j-1
      q(j,j) = q(j,j)
     .  -(dgi(jp)/dgi(j))*cez(jp,j,kz)*denz(kz)
   42 continue
      do 43 i=j+1,ns
      q(i,j) = q(i,j)
     .  +cez(j,i,kz)*denz(kz)
   43 continue
      do 44 i=1,j-1
      q(i,j) = q(i,j)
     .  +(dgi(i)/dgi(j))*cez(i,j,kz)*denz(kz)
   44 continue
   40 continue
c
      call SECOND (cpub)
      cpu6 = cpu6 + cpub - cpua
      return
c
      end

      subroutine nbdep2d (psi, mi, mj, r, z, potsid, mf, rzpat, nrpat,
     .                    nzpat, me, mb, sgxn, vbeam, hxfrac)
c
      USE param
      USE mhdpar
      USE io 
      USE replace_imsl, ONLY : my_ibcieu
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c calculates neutral beam deposition on (r,z) grid.  grid size
c determined by nrpat and nzpat, the number of equally spaced
c elements in the r and z axes.  outputs the 3rd excited state
c component, rzhex, to file 'beamdep' for postprocessing.
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'io.i'
c
      parameter (kzm1 = kz-1, nxny = nw*nh, nxx2 = nw*2, nyx2 = nh*2,
     .           kwork = 3*(nh-1)+nh)
      dimension  hxfrac(ke,kb),psi(nw,nh),r(nw),z(nh),potsid(kf),
     .           rzpat(nxx2,nyx2,ke,kb),sgxn(kcmp1,kz,kbe,*),
     .           vbeam(ke,kb)
      dimension  bpar1(4),bpar2(4),capsig(kz),cspln1(kzm1,3),
     .           inc(nyx2),isupp(4,ke,kb),pc(kz),
     .           psikp(nxny),psipat(nxx2,nyx2),rr(nxx2),
     .           rznub(nxx2,nyx2,ke,kb),rzsig(nxny),splnwk(kwork),
     .           zz(nyx2)
c
      logical    exist
c
c     set (r,z) grid modification option
c
      igrid = 0
      if (nrpat .ne. mi .or. nzpat .ne. mj)  igrid = 1 ! per Jeff Moller
c
c     zero out arrays
c
      do i=1,4
        bpar1(i) = 0.0
        bpar2(i) = 0.0
      end do
c
c     set up integer array to allow vectorization further on
c
      ndum = MAX (nrpat, nzpat)
      if (ndum .gt. nyx2)
     .  call STOP ('subroutine NBDEP2D: Jeff Moller bailout', 68)
c
      do i=1,ndum
        inc(i) = i-1
      end do
c
c     get psi on zone centers
c
      mfm1 = mf-1
      do i=1,mfm1
        pc(i) = 0.5 * (potsid(i) + potsid(i+1))
      end do
c
c     if user defined new (r,z) grid, interpolate psi(i,j) onto it
c
      if (igrid .eq. 1) then
        dr      = (r(mi)-r(1)) / (nrpat-1)
        rr(1)   = r(1)
        do i=2,nrpat
          rr(i) = rr(1)+inc(i)*dr
        end do
        dz      = (z(mj)-z(1))/(nzpat-1)
        zz(1)   = z(1)
        do i=2,nzpat
          zz(i) = zz(1)+inc(i)*dz
        end do
        call my_ibcieu (psi, nw, r, mi, z, mj, rr, nrpat, zz, nzpat,
     .               psipat, nxx2, splnwk, ier)
      end if
c
c     begin loop over beam and beam energy
c
      do   ib=1,mb
        do ie=1,me
c
c         get an estimate of the support of rzpat
c
          call support (rzpat, nrpat, nzpat, ie, ib, isupp, ifail)
c
c         get spline fits to macroscopic neutral beam attenuation cross
c         sections as a function of psi
c
          ind = ke*(ib-1) + ie
          do i=1,mfm1
            capsig(i) = sgxn(4,i,ind,1)
          end do
c
          call icsicu1 (pc, capsig, mfm1, bpar1, cspln1, kzm1, ier)
c
c         loop over (r,z) points
c
          do   i=isupp(1,ie,ib),isupp(2,ie,ib)
            do j=isupp(3,ie,ib),isupp(4,ie,ib)
              if (rzpat(i,j,ie,ib) .ne. 0) then
                n = n+1
                if (igrid .eq. 0) then
                  psikp(n) = psi   (i,j)
                else
                  psikp(n) = psipat(i,j)
                end if
              end if
            end do
          end do
c
          call icsevu1 (pc, capsig, mfm1, cspln1, kzm1,
     .                  psikp, rzsig, n,ier)
          n = 0
          do   i=isupp(1,ie,ib),isupp(2,ie,ib)
            do j=isupp(3,ie,ib),isupp(4,ie,ib)
              if (rzpat(i,j,ie,ib) .ne. 0.0) then
                n = n + 1
                rznub(i,j,ie,ib) = rzpat(i,j,ie,ib)
     .                           / (vbeam(ie,ib) * rzsig(n))
              end if
            end do
          end do
        end do
      end do
c
c     output subset of rzhex to file 'beamdep'
c
      nout = nbdep
      inquire (unit = nout, exist = exist)
      if (exist) then
        call DESTROY ('beamdep')
      else
        call getioun(nbdep,nbdep)
        nout = nbdep
      end if
      open  (unit = nout, file = 'beamdep', status = 'UNKNOWN')
      write (nout, 1000)  nrpat, nzpat, me, mb
      if (igrid .eq. 0) then
        write (nout, 1010)   r(1),  r(mi)   ,  z(1),  z(mj)
      else
        write (nout, 1010)  rr(1), rr(nrpat), zz(1), zz(nzpat)
      end if
c
      do   ib=1,mb
        do ie=1,me
          write (nout, 1000) (isupp(i,ie,ib),i=1,4)
          write (nout, 1010) hxfrac(ie,ib)
          write (nout, 1010)
     .    ((rznub(i,j,ie,ib), i=isupp(1,ie,ib),isupp(2,ie,ib)),
     .                        j=isupp(3,ie,ib),isupp(4,ie,ib))
        end do
      end do
c
 1000 format (4(5x,i4))
 1010 format (7(2x,e16.7))
      call giveupus(nout)
      close (unit = nout)
      return
c
      end

      subroutine nbsgold (ebkev,ebfac,ibion,mb,mfm1,nebin,vbeam,zne,zni,
     .                    zte,zzi,debin,sgxn,sgxnmi,atw_beam)
c
c
c ----------------------------------------------------------------------
c
c     This subroutine calculates the neutral beam attenuation array sgxn
c     using the old cross sections based on Freeman & Jones (1972).
c     No excited state effects are included. These cross sections have
c     been found to be outdated and are not to used except for purposes
c     of comparison to the newer ADAS data.
c
c     Created:     12-JUL-1995     Daniel Finkenthal
c
c     input:
c
c          ebkev(mb)      - maximum energy of mb-th neutral beamline
c                           (keV)
c          ebfac          - factor defining upper bound on energy bin
c                           range,  > 1.0
c          ibion          - index of beam ion species
c                           (e.g., atwb = atw(ibion))
c          mb             - number of beamlines modeled
c          mfm1           - number of flux zones minus one
c          nebin          - number of energy bins (rotation case only)
c          vbeam(ie,mb)   - speed of ie-th energy group of the mb-th
c                           beamline (cm/sec)
c          zne(mfm1)      - local electron density (cm**-3)
c          zni(mfm1,nion) - local density of nion-th ion species
c                           (cm**-3)
c          zzi(mfm1,nion) - local average charge state of nion-th ion
c                           species
c     input from common /io/:
c          ncrt,nout,nouthx,ncorin
c     input from common /ions/:
c          namei, atw
c     input from common /numbrs/:
c          nprim,nimp,nion
c
c     output:
c          sgxn(i,j,k,l)
c               i - mode index
c                   = 1, fraction of reactions producing electrons;
c                   = 2, fraction of reactions producing species 1 ion;
c                   = 3, fraction of reactions producing species 2 ion;
c                   = 4, total inverse mean free path;
c               j - FREYA zone index
c               k - beam index
c                   k = 3*(ib-1) + ie, where ib is the beam index and ie
c                                      is the beam energy group index
c               l - index for relative energy.  bins are equispaced in
c                   delta(energy) between 0 and max(ebkev(ib))*ebfac,
c                   with bin width given by
c                             delta(energy) = ebmax*ebfac/nebin.
c                   ebfac .ge. 1.0 and nebin are user supplied namelist
c                   variables.
c          sgxnmi(ie,mb)
c                 - minimum inverse mean free path for ie-th energy
c                   group of mb-th beamline.  calculated for stationary
c                   target case only.  otherwise, calculated for each
c                   neutral trajectory in subroutine INJECT.
c          debin
c                 - width of energy bins (keV/amu)
c ----------------------------------------------------------------------
c
      USE param
      USE ions
      USE io
      USE numbrs 
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'io.i'
c      include 'numbrs.i'
c      include 'ions.i'
c
      dimension ebkev(*),vbeam(ke,*),zne(*),zni(kz,*),zte(*),zzi(kz,*)
      dimension sgxn(kcmp1,kz,kbe,*),sgxnmi(ke,*)
      dimension sgxne(kbe)
c
c ----------------------------------------------------------------------
c stationary plasma case
c ----------------------------------------------------------------------
c
      if (nebin .eq. 0) then
        do 220 ib=1,mb
        do 220 ie=1,ke
        sgxnmi(ie,ib) = 0.0
        ind  = ke*(ib-1) + ie
        eova = 1.0e3*ebkev(ib)/(FLOAT (ie)*atw_beam)
        vbin = vbeam(ie,ib)
        do 220 i=1,mfm1
        teev = 1.0e3 * zte(i)
        sgxne(    ind  ) = fsgxne(vbin,teev,zne(i))
        sgxn (1,i,ind,1) = 0.0
        do k=1,nion
          sgxn(1,i,ind,1) = sgxn(1,i,ind,1)
     .                    + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
          if ((k .le. nprim) .and. (k .le. 2.0))
     .      sgxn(k+1,i,ind,1) = fsgxncx(atw(k),eova,zni(i,k))
        end do
        sgxn(1,i,ind,1) = sgxn(1,i,ind,1) + sgxne(ind)
        sgxn(4,i,ind,1) = sgxn(1,i,ind,1) + sgxn(2,i,ind,1)
     .                  + sgxn(3,i,ind,1)
        sgxn(1,i,ind,1) = sgxn(1,i,ind,1)/sgxn(4,i,ind,1)
        sgxn(2,i,ind,1) = sgxn(2,i,ind,1)/sgxn(4,i,ind,1)
        sgxn(3,i,ind,1) = sgxn(3,i,ind,1)/sgxn(4,i,ind,1)
        sgxnmi(ie,ib) = MAX (sgxnmi(ie,ib),sgxn(4,i,ind,1))
  220   continue
        do   ib=1,mb
          do ie=1,ke
            sgxnmi(ie,ib) = 1.0 / sgxnmi(ie,ib)
          end do
        end do
      else
c
c ----------------------------------------------------------------------
c rotating plasma case
c ----------------------------------------------------------------------
c
        ebmax = 0.0
        do ib=1,mb
           ebmax = MAX (ebmax,ebkev(ib))
        end do
        ebmax = ebmax/atw_beam
        debin = ebmax*ebfac/FLOAT (nebin)
c
        do 340 i=1,mfm1
        teev = 1.0e3*zte(i)
c       electron impact ionization independent of rotation speed
        do 320 ib=1,mb
        do 320 ie=1,ke
        ind = ke*(ib-1) + ie
        sgxne(ind) = fsgxne(vbeam(ie,ib),teev,zne(i))
  320   continue
        do 340 j=1,nebin
        ebin = FLOAT (j)*debin
        eova = 1.0e3*ebin
        sgxn(1,i,1,j) = 0.0
        do 330 k=1,nion
        sgxn(1,i,1,j) = sgxn(1,i,1,j)
     .                + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
        if ((k .le. nprim) .and. (k .le. 2))
     .    sgxn(k+1,i,1,j) = fsgxncx(atw(k),eova,zni(i,k))
  330   continue

        do 340 k=ind,1,-1       ! (ind=ke*(ib-1) +ie,beam energy index)
          sgxn(1,i,k,j) = sgxn(1,i,1,j) + sgxne(k)
          sgxn(4,i,k,j) = sgxn(1,i,k,j) + sgxn(2,i,1,j) + sgxn(3,i,1,j)
          sgxn(1,i,k,j) = sgxn(1,i,k,j)/sgxn(4,i,k,j)
          sgxn(2,i,k,j) = sgxn(2,i,1,j)/sgxn(4,i,k,j)
          sgxn(3,i,k,j) = sgxn(3,i,1,j)/sgxn(4,i,k,j)
  340   continue
      end if

      return
c
      end

      subroutine nbsgxn (ebkev, ebfac, ibion, mb, mfm1, nebin, vbeam,
     .                   zne, zni, zte, zzi, debin, hxfrac, sgxn,
     .                   sgxnmi,atw_beam)
c

c
c ----------------------------------------------------------------------
c
c  this subroutine calculates the neutral beam attenuation array, sgxn.
c
c     input:
c
c          ebkev(mb)      - maximum energy of mb-th neutral beamline
c                           (keV)
c          ebfac          - factor defining upper bound on energy bin
c                           range,  > 1.0
c          atw_beam         mass no. of beam
c
c          mb             - number of beamlines modeled
c          mfm1           - number of flux zones minus one
c          nebin          - number of energy bins (rotation case only)
c          vbeam(ie,mb)   - speed of ie-th energy group of the mb-th
c                           beamline (cm/sec)
c          zne(mfm1)      - local electron density (cm**-3)
c          zni(mfm1,nion) - local density of nion-th ion species
c                           (cm**-3)
c          zzi(mfm1,nion) - local average charge state of nion-th ion
c                           species
c     input from common /io/:
c          ncrt,nout,nouthx,ncorin
c     input from common /ions/:
c          namei, atw
c     input from common /nub3/:
c          iexcit,ilorent,mstate,ncont,kdene,kdeni,kdenz,ksvi,ksvz,ksve,
c          krad,ngh,ngl,znipm,atwpm,iz,zti,izstrp
c     input from common /numbrs/:
c          nprim,nimp,nion
c
c     output:
c          sgxn(i,j,k,l)
c               i - mode index
c                   = 1, fraction of reactions producing electrons;
c                   = 2, fraction of reactions producing species 1 ion;
c                   = 3, fraction of reactions producing species 2 ion;
c                   = 4, total inverse mean free path;
c               j - FREYA zone index
c               k - beam index
c                   k = 3*(ib-1) + ie, where ib is the beam index and ie
c                                      is the beam energy group index
c               l - index for relative energy.  bins are equispaced in
c                   delta(energy) between 0 and max(ebkev(ib))*ebfac,
c                   with bin width given by
c                             delta(energy) = ebmax*ebfac/nebin.
c                   ebfac .ge. 1.0 and nebin are user supplied namelist
c                   variables.
c          sgxnmi(ie,mb)
c                 - minimum inverse mean free path for ie-th energy
c                   group of mb-th beamline.  calculated for stationary
c                   target case only.  otherwise, calculated for each
c                   neutral trajectory in subroutine INJECT.
c          hxfrac(ie,mb)
c                 - fractional density of n = 3 excited state neutral.
c                   calculated in stationary target case only.
c          debin
c                 - width of energy bins (keV/amu)
c ----------------------------------------------------------------------
c
      USE param
      USE ions
      USE io 
      USE numbrs
      USE nub3
      implicit  integer (i-n), real*8 (a-h, o-z)
D     include 'mpif.h'
c      include 'param.i'
c      include 'nub3.i'
c      include 'io.i'
c      include 'numbrs.i'
c      include 'ions.i'
c
      parameter (ksrc = 3)
      dimension  ebkev(*),vbeam(ke,*),zne(*),zni(kz,*),zte(*),zzi(kz,*)
      dimension  hxfrac(ke,*),sgxn(kcmp1,kz,kbe,*),sgxnmi(ke,*)
      dimension  sgxne(kbe)
      myid =0
      master = myid
D     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr) !get processor id
D     call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr) !get num processors 
c
c ----------------------------------------------------------------------
c Go right to ADAS routines and return if  iexcit = 5
c ----------------------------------------------------------------------
c
       
      if (iexcit .eq. 5) then
        if(myid .eq. master )
     .   call adassgxn (ebkev,ebfac,ibion,mb,mfm1,nebin,vbeam,
     .                 zne,zni,zte,zti,zzi,debin,sgxn,sgxnmi,
     .                                              atw_beam)

D        if(numprocs .gt. 1)then
D             ksg = kcmp1*kz*kbe*ksge
D             call MPI_BCAST(sgxn,ksg,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D             call MPI_BCAST(sgxnmi,kbe,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D             call MPI_BCAST(hxfrac,kbe,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D             call MPI_BCAST(debin,1,MPI_Double_Precision,
D    .                      master,MPI_COMM_WORLD,ierr)
D        endif

        return
      end if
c
c ---------------------------------------------------- HSJ-2/5/98 ------
c --- Boley's parameterization including MSI effects:
c ----------------------------------------------------------------------
c
      if (iexcit .eq. 6) then
        call wrap_xboley (atw,zzi,ebkev,ibion,mb,mfm1,nebin,ebfac,
     .                     debin,zne,zte,zni,nion,sgxn,sgxnmi,atw_beam)
        return
      end if
c
c ----------------------------------------------------------------------
c HEXNB initialization
c ----------------------------------------------------------------------
c
      if (iexcit .ne. 0) then
        if (nouthx .gt. 0) then
          call getioun(nouthx,nouthx)
          open (unit = nouthx, file = 'hexnbout', status = 'UNKNOWN')
        end if
        istart = 1
c
c       separate primary and impurity ions
c
        do j=1,nion
          if (j .le. nprim) then
            atwpm(j) = atw(j)
          else
            k        = j - nprim
            atwim(k) = atw(j)
          end if
        end do
c
c       get atomic number of impurities
c
        do i=1,nimp
          if (namei(i) .eq. 'he')  iz(i) =  2
          if (namei(i) .eq. 'c ')  iz(i) =  6
          if (namei(i) .eq. 'o ')  iz(i) =  8
          if (namei(i) .eq. 'si')  iz(i) = 14
          if (namei(i) .eq. 'ar')  iz(i) = 18
          if (namei(i) .eq. 'cr')  iz(i) = 24
          if (namei(i) .eq. 'fe')  iz(i) = 26
          if (namei(i) .eq. 'ni')  iz(i) = 28
          if (namei(i) .eq. 'kr')  iz(i) = 36
          if (namei(i) .eq. 'mo')  iz(i) = 42
          if (namei(i) .eq. 'w ')  iz(i) = 74
        end do
c
c       no Lorentz ionization limit
c
        bperp = 0.0
      end if
c
c ----------------------------------------------------------------------
c stationary plasma case
c ----------------------------------------------------------------------
c
      if (nebin .eq. 0) then
        do   ib=1,mb
          do ie=1,ke
            sgxnmi(ie,ib) = 0.0
            ind  = ke*(ib-1) + ie
            eova = 1.0e3*ebkev(ib)/(FLOAT (ie)*atw_beam)
            vbin = vbeam(ie,ib)
            do i=1,mfm1
              teev = 1.0e3*zte(i)
              if (iexcit .ne. 0) then
                tiev = 1.0e3*zti(i)
                do j=1,nion
                  if (j .le. nprim) then
                    znipm(j) = zni(i,j)
                  else
                    k        = j - nprim
                    zniim(k) = zni(i,j)
                  end if
                end do
              end if
              sgxne(ind) = fsgxne(vbin,teev,zne(i))
              sgxn(1,i,ind,1) = 0.0
              do k=1,nion
                sgxn(1,i,ind,1) = sgxn(1,i,ind,1)
     .                          + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
                if ((k .le. nprim) .and. (k .le. 2.0))
     .            sgxn(k+1,i,ind,1) = fsgxncx(atw(k),eova,zni(i,k))
              end do
              sgxn(1,i,ind,1) = sgxn(1,i,ind,1) + sgxne(ind)
              sgxn(4,i,ind,1) = sgxn(1,i,ind,1) + sgxn(2,i,ind,1)
     .                        + sgxn(3,i,ind,1)
              sgxn(1,i,ind,1) = sgxn(1,i,ind,1)/sgxn(4,i,ind,1)
              sgxn(2,i,ind,1) = sgxn(2,i,ind,1)/sgxn(4,i,ind,1)
              sgxn(3,i,ind,1) = sgxn(3,i,ind,1)/sgxn(4,i,ind,1)
              if (iexcit .ne. 0) then
                call hexnb (istart, 1, ilorent, mstate, ncont, eova,
     .                      teev, tiev, nprim, atwpm, znipm, nimp, iz,
     .                      atwim, izstrp, zniim, bperp, kdene, kdeni,
     .                      kdenz, ksvi, ksvz, ksve, krad, ngh, ngl,
     .                      nouthx, ncorin, rerate, rmfp, hexfrac,
     .                      ihxerr)
                if (ihxerr .ne. 0) write (ncrt, 1000) ihxerr
                if (ihxerr .eq. 1) then
                  write (ncrt, 1010)
                  write (nout, 1010)
                  call STOP ('subroutine NBSGXN: problem #1', 45)
                end if
                sgxn(4,i,ind,1) = 1.0/rmfp
                if (i .eq. 1) hxfrac(ie,ib) = hexfrac
              end if
              sgxnmi(ie,ib) = MAX (sgxnmi(ie,ib),sgxn(4,i,ind,1))
            end do
          end do
        end do
        do   ib=1,mb
          do ie=1,ke
            sgxnmi(ie,ib) = 1.0 / sgxnmi(ie,ib)
          end do
        end do
      else
c
c ----------------------------------------------------------------------
c rotating plasma case
c ----------------------------------------------------------------------
c
        ebmax = 0.0
        do ib=1,mb
          ebmax = MAX (ebmax,ebkev(ib))
        end do
        ebmax = ebmax/atw_beam
        debin = ebmax*ebfac/FLOAT (nebin)
c
        do 340 i=1,mfm1
        teev = 1.0e3 * zte(i)
        if (iexcit .ne. 0) then
          tiev = 1.0e3 * zti(i)
          do j=1,nion
            if (j .le. nprim) then
              znipm(j      ) = zni(i,j)
            else
              zniim(j-nprim) = zni(i,j)
            end if
          end do
        end if
c
c       electron impact ionization independent of rotation speed
c
        do   ib=1,mb
          do ie=1,ke
            ind        = ke*(ib-1) + ie
            sgxne(ind) = fsgxne(vbeam(ie,ib),teev,zne(i))
          end do
        end do
c
        do 340 j=1,nebin
        ebin = FLOAT (j) * debin
        eova =     1.0e3 * ebin
        sgxn(1,i,1,j) = 0.0
        do k=1,nion
          sgxn(1,i,1,j) = sgxn(1,i,1,j)
     .                  + fsgxni(atw(k),eova,zni(i,k),zzi(i,k))
          if ((k .le. nprim) .and. (k .le. 2.0))
     .      sgxn(k+1,i,1,j) = fsgxncx(atw(k),eova,zni(i,k))
        end do
        if (iexcit .ne. 0) then
          call hexnb (istart, 1, ilorent, mstate, ncont, eova,
     .                teev, tiev, nprim, atwpm, znipm, nimp, iz,
     .                atwim, izstrp, zniim, bperp, kdene, kdeni,
     .                kdenz, ksvi, ksvz, ksve, krad, ngh, ngl, nouthx,
     .                ncorin, rerate, rmfp, hexfrac, ihxerr)
          if (ihxerr .ne. 0) write (ncrt, 1000) ihxerr
          if (ihxerr .eq. 1) then
            write (ncrt, 1010)
            write (nout, 1010)
            call STOP ('subroutine NBSGXN: problem #2', 46)
          end if
        end if
        do 340 k=ind,1,-1
          sgxn(1,i,k,j) = sgxn(1,i,1,j) + sgxne(k)
          sgxn(4,i,k,j) = sgxn(1,i,k,j) + sgxn(2,i,1,j) + sgxn(3,i,1,j)
          sgxn(1,i,k,j) = sgxn(1,i,k,j)/sgxn(4,i,k,j)
          sgxn(2,i,k,j) = sgxn(2,i,1,j)/sgxn(4,i,k,j)
          sgxn(3,i,k,j) = sgxn(3,i,1,j)/sgxn(4,i,k,j)
          if (iexcit .ne. 0)  sgxn(4,i,k,j) = 1.0/rmfp
  340   continue
      end if
c
 1000 format (' subroutine NBSGXN reports a HEXNB return code of ', i5)
 1010 format (' ERROR: execution terminated - file "coronb" not found')
c
      return
c
      end

      subroutine newgrid (xold, yold, nold, xnew, ynew, nnew)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c convert to new grid
c ----------------------------------------------------------------------
c
      dimension  xold(*), yold(*), xnew(*), ynew(*)
      data       tol/1.0e-20/
c
      iold = 1
      inew = 1
      sex  = (yold(nold)-yold(nold-1))/(xold(nold)-xold(nold-1))
c
c ----------------------------------------------------------------------
c interpolate to new grid that is compute ynew(xnew(inew))
c ----------------------------------------------------------------------
c
 2000 if (inew .gt. nnew)  return
 2100 if (iold .gt. nold)  go to 2200
      if (xnew(inew) .gt. xold(1))  go to 2105
      ynew(inew) = yold(1)
      inew       = inew+1
      go to 2000
c
 2105 del  = reldif(xold(iold),xnew(inew))
      if (del .gt. tol)  go to 2110
      ynew(inew) = yold(iold)
      inew = inew+1
      go to 2000
 2110 if (xold(iold) .gt. xnew(inew))  go to 2120
      iold = iold+1
      go to 2100
 2120 ipm1 = iold-1
      s    = (yold(iold)-yold(ipm1))/(xold(iold)-xold(ipm1))
      ynew(inew) = yold(ipm1)+(xnew(inew)-xold(ipm1))*s
      inew = inew+1
      go to 2000
 2200 ynew(inew) = yold(nold)+sex*(xnew(inew)-xold(nold))
      if (ynew(inew) .lt.  0.0)
     .ynew(inew) = 0.0
      inew       = inew + 1
      go to 2000
c
      end

      integer function nfhx (i)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      nfhx = i - 1
      if (i .le. 2)
     .nfhx = i
      return
c
      end

      subroutine nubplt
c
      USE param
      USE io 
      USE neut
      USE nub
      USE nub2
      USE solcon
      USE soln
      USE numbrs
      USE mesh
      USE transp,only : use_nubeam
      USE machin
      USE geom
      USE tordlrot
      USE oloss_dat
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     record debug information
c


      include 'storage.i'

c
      integer mcgo_run
      dimension    xp(kstore),yp(kstore),work(kstore)
      equivalence (xdum(1),xp(1)), (ydum(1),yp(1)), (zdum(1),work(1))
c
c ----------------------------------------------------------------------
c if 'onedee', unit nscr is not opened and closed in this routine
c ----------------------------------------------------------------------
c


C     not yet set up for nubeam
      if(use_nubeam)return

      if (codeid .eq. 'onedee')  go to 1600
      if (beamon(1) .gt. timmax)  return
      mcgo_run=0
      if (iborb .eq. 3)  mcgo_run = 1
 1600 write (nbplt, '(a)') '****continue****'
      write (nbplt,  8005)  time,ene(1),te(1),rmax,rmin,npts,mcgo_run
      write (nbplt,  8010) (xpts(i),   i=1,npts)
      write (nbplt,  8010) (ypts(i),   i=1,npts)
      write (nbplt,  8010) (zpts(i),   i=1,npts)
      write (nbplt,  8010) (rpts(i),   i=1,npts)
!     new output , onetwo v3.94
      write (nbplt,  8010) (vx(i),     i=1,npts)
      write (nbplt,  8010) (vy(i),     i=1,npts)
      write (nbplt,  8010) (vz(i),     i=1,npts)
      write (nbplt,  8010) (pitch_a(i),i=1,npts)

      if (iborb .gt. 1) then
        write (nbplt, 8010)  (storque  (i), i=1,nj)
        write (nbplt, 8010)  (storqueb (i), i=1,nj)
        write (nbplt, 8010)  (spbolt   (i), i=1,nj)
        write (nbplt, 8010)  (pitchang (i), i=1,npitch)
        write (nbplt, 8010)  (pitchrmaj(i), i=1,mf)
        do i=1,mf
          write (nbplt, 8010)  (pitchlos(j,i), j=1,npitch)
        end do
      end if
c
      do jb=1,nbeams
        write (nbplt, 8010) (bneut(i,jb),i=1,3)
        write (nbplt, 8010) (ebeam(i,jb),i=1,3)
        write (nbplt, 8010) (pbeam(i,jb),i=1,3)
        write (nbplt, 8010) (fap(i,jb),i=1,3)
        write (nbplt, 8010) (fwall(i,jb),i=1,3)
        write (nbplt, 8010) (forb(i,jb),i=1,3)
        do ie=1,3
          write (nbplt, 8010) (hibr  (j,ie,jb), j=1,nj)
          write (nbplt, 8010) (hdep  (j,ie,jb), j=1,nj)
          write (nbplt, 8010) (angmpf(j,ie,jb), j=1,nj)
        end do
      end do
c
      write (nbplt, 8010)  (r(j), j=1,nj)
      if (codeid .eq. 'onedee')  return
      call getioun(nscr,nscr)
      open (unit = nscr, file = 'scratch', status = 'OLD')
      do i=1,mf
        read  (nscr , 8100) mp
        write (nbplt, 8100) mp
        if (mp .eq. 0)  go to 2030
        read  (nscr , 8120)  (xp(j), yp(j), j=1,mp)
        write (nbplt, 8120)  (xp(j), yp(j), j=1,mp)
      end do
c
 2030 call giveupus(nscr)
      close (unit = nscr)
      return
c
 8100 format (i6)
 8120 format (6e12.5)
 8005 format (5e12.3,2i10)
 8010 format (6e12.3)
c
      end

      subroutine orbit_12 (atw,b1ins,b1ots,b2ins,b2ots, codeid,
     .                  ic,idebug,iskip,izi,mf,norb,
     .                  pinsid,potsid,pzone,rinsid,rotsid,rzone,ri,v,
     .                  zetai,zi,finsid,fotsid,zinsid,zotsid,
     .                  ipass,iaxis,ier,izp,wid,wt,zeta)
c
      USE constnts, only : xmassp
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  This subroutine calculates the following orbit parameters for a
c     fast ion:
c       ipass: 1, particle is passing
c              0, particle is trapped
c       iaxis: 1, particle circles the axis
c              0, particle does not circle the axis
c       ier  : 1, routine failed to obtain orbit
c              0, routine obtained reasonable orbit
c       izp  : outermost zone of orbit
c       wid  : width of orbit
c       wt   : weight giving approximate fraction of time spent
c              by fast ion in each zone
c       zeta : approximate pitch angle cosine in each zone along
c              orbit
c
      character  codeid*8, flagm*1, flagp*1
      dimension  b1ins(*), b1ots(*), b2ins(*), b2ots(*)
      dimension  idebug(*)
      dimension  pinsid(*), potsid(*), rinsid(*), rotsid(*)
      dimension  finsid(*), fotsid(*), zinsid(*), zotsid(*)
      dimension  wt(*), zeta(*)
      dimension  rin(3), rot(3)
      data eoverc / 1.600e-20 /
c
c     initialize parameters
c
      zero   = 0.0
      mfm1   = mf-1
      drini  = mfm1/(rinsid(1)-rinsid(mf))
      droti  = mfm1/(rotsid(mf)-rotsid(1))
      drutpi = mfm1 / SQRT (potsid(mf)-potsid(1))
      flagm  = ' '
      flagp  = ' '
      ier    = 0
      izm    = izi
      izp    = izi
      do 10 i=1,3
      rin(i) = 0.0
   10 rot(i) = 0.0
      rm     = 0.0
      rp     = 0.0
      zetam  = 0.0
      zetap  = 0.0
      wid    = 0.0
      do 20 i=1,mfm1
      wt(i)  = 0.0
   20 zeta(i) = 0.0
      wt(izi) = 1.0
      zeta(izi) = zetai
c
c calculate ion mass times c/e; initialize rmajor and psiax
c
      xmass  = atw*xmassp/eoverc
      rmajor = rotsid(1)
      psiax  = potsid(1)
c
c  calculate poloidal magnetic flux at the fast ion birth point
c
      if (codeid .eq. 'onedee') then
        psii = yinter(droti,mfm1,0.0,rzone,potsid)
      else
        psii = pzone
      end if
c
c  approximate ratios of magnetic fields at the fast ion birth point
c
      if (ri .ge. rmajor) then
        b1i = yinter(droti,mfm1,rmajor,ri,b1ots)
        b2i = yinter(droti,mfm1,rmajor,ri,b2ots)
      else
        b1i = yinter(drini,mfm1,rmajor,ri,b1ins)
        b2i = yinter(drini,mfm1,rmajor,ri,b2ins)
      end if
c
c  calculate constants for orbit
c
      psiref = xmass*v*ri
      pang   = psiref*zetai/b2i - psii
      c1     = (1.0-zetai**2)*ri/b1i
      c2     = ri/psiref
      c3     = pang
c
c  determine whether fast ion is passing or trapped
c
      if (ABS (zetai) .eq. 1.0) then
        ipass = 1
      else if (zetai .eq. 0.0) then
        ipass = 0
      else
        rtip = (1.0-zetai**2)*ri/b1i
        if (rtip .ge. rmajor) then
          b1tip  = yinter(droti,mfm1,rmajor,rtip,b1ots)
          rtip   = rtip*b1tip
          psitip = yinter(droti,mfm1,rmajor,rtip,potsid)
        else
          b1tip  = yinter(drini,mfm1,rmajor,rtip,b1ins)
          rtip   = rtip*b1tip
          psitip = yinter(drini,mfm1,rmajor,rtip,pinsid)
        end if
        pangt = -psitip
        if (pang .gt. pangt) then
          ipass = 1
        else
          ipass = 0
        end if
      end if
c
c  determine whether ion circles the magnetic axis or not
c
      iaxis  = 1
      zetax2 = 1.0 - (1.0-zetai**2)*ri/(rmajor*b1i)
      if (zetax2 .le. 0.0)  go to 80
      zetax  = SQRT (zetax2)
      if (ipass .eq. 0 .or. zetai .lt. 0.0) zetax = -zetax
      pangx = xmass*v*zetax*rmajor - psiax
      if (ipass .eq. 1 .and. pang .lt. pangx)  go to 90
      if (ipass .eq. 0 .and. pang .gt. pangx)  go to 90
   80 iaxis = 0
c
c  find roots on inside and outside of axis
c
   90 do i=1,mf
        finsid(i) = rinsid(i)**2 - c1*rinsid(i)*b1ins(i)
     .              - (c2*b2ins(i)*(c3+pinsid(i)))**2
        fotsid(i) = rotsid(i)**2 - c1*rotsid(i)*b1ots(i)
     .              - (c2*b2ots(i)*(c3+potsid(i)))**2
      end do
c
      nin = 0
      not = 0
c
      do i=1,mfm1
        if (finsid(i) .eq. 0.0 .or. finsid(i)*finsid(i+1) .lt. 0.0) then
          nin = nin+1
          if(nin .le. 3)
     .     rin(nin) = (rinsid(i+1)*finsid(i)-rinsid(i)*finsid(i+1))
     .              /(finsid(i)-finsid(i+1))
        end if
        if (fotsid(i) .eq. 0.0 .or. fotsid(i)*fotsid(i+1) .lt. 0.0) then
          not = not+1
          if(not .le. 3)
     .     rot(not) = (rotsid(i+1)*fotsid(i)-rotsid(i)*fotsid(i+1))
     .              /(fotsid(i)-fotsid(i+1))
        end if
      end do
c
c  match roots to orbit extrema
c
      if (nin .gt. 3 .or. not .gt. 3) then
        flagp = '*'
        go to 340
      end if
      if (ipass .eq. 0 .and. iaxis .eq. 0) then
        if (not .lt. 2)  go to 210
        rp = rot(2)
        rm = rot(1)
      else if (ipass .eq. 1 .and. iaxis .eq. 1) then
        if (zetai .ge. 0.0) then
          if (nin .eq. 1) then
            if (not .eq. 0)  go to 210
            rp = rot(1)
            rm = rin(1)
          else if (nin .eq. 2) then
            if (not .le. 1)  go to 210
            rp = rot(2)
            rm = rin(2)
            im = (rmajor-rm)*drini + 2.0
            if (finsid(im-2) .lt. 0.0) then
              rm = rfine(b1ins,b2ins,c1,c2,c3,finsid,im,1,mfm1,pinsid,
     .                   rinsid,rmajor)
            end if
          else if (nin .eq. 0) then
            if (not .le. 1)  go to 210
            rp = rot(2)
            im = imax(finsid,mfm1)
            if (finsid(im+1) .gt. finsid(im-1)) im = im+1
            rm = rfine(b1ins,b2ins,c1,c2,c3,finsid,im,1,mfm1,pinsid,
     .                 rinsid,rmajor)
            if (rm .eq. 0.0) then
              flagm = '*'
              go to 340
            end if
          else
            if (not .eq. 0)  go to 210
            rp = rot(1)
            rm = rin(3)
          end if
        else
          if (nin .eq. 0)  go to 210
          rp = rin(1)
          rm = rot(1)
          ip = (rmajor-rp)*drini + 1.0
          if (finsid(ip+2) .lt. 0.0) then
            rp = rfine(b1ins,b2ins,c1,c2,c3,finsid,ip,-1,mfm1,pinsid,
     .                 rinsid,rmajor)
          end if
        end if
      else if (ipass .eq. 1 .and. iaxis .eq. 0) then
        if (zetai .ge. 0.0) then
          if (not .lt. 2)  go to 210
          rp = rot(2)
          rm = rot(1)
        else
          if (nin .lt. 2)  go to 210
          rp = rin(2)
          rm = rin(1)
        end if
      else if (ipass .eq. 0 .and. iaxis .eq. 1) then
        if (not .eq. 0)  go to 210
        rp = rot(1)
        rm = rin(1)
      end if
      go to 220
  210 izp = mf
      go to 350
c
c  calculate poloidal magnetic flux and index of outermost zone of orbit
c
  220 if (rp .ge. rmajor) then
        psip = yinter(droti,mfm1,rmajor,rp,potsid)
        b2p = yinter(droti,mfm1,rmajor,rp,b2ots)
      else
        psip = yinter(drini,mfm1,rmajor,rp,pinsid)
        b2p = yinter(drini,mfm1,rmajor,rp,b2ins)
      end if
      if (codeid .eq. 'onedee') then
        xzp = ABS (rp-rmajor)*droti + 1.0
      else
        xzp = SQRT (psip-psiax)*drutpi + 1.0
      end if
      izp = xzp
      izp = MIN0 (izp,mf)
      zetap = (pang+psip)*b2p/(xmass*v*rp)
      isp = 1
      if (zetai .lt. 0. .and. ipass .eq. 1) isp = -1
      if (isp*zetap .lt. -0.01 .or. izp .lt. izi-1) then
        flagp = '*'
        go to 340
      end if
c
c  determine whether ion is confined or lost; return if it is lost
c
      if (izp .gt. mfm1)  go to 350
c        if(rm .eq. 0)then
c           print *,ipass,iaxis
c           print *,'rin =',rin
c           print *,'rot =',rot
c           print *,'zetai,nin =',zetai,nin
c        endif
      if(rm .le. 0.00) go to 350  !ctr injection lost orbit ?? 
c
c  calculate poloidal magnetic flux and index of innermost zone of orbit
c
      if (rm .lt. rmajor) then
        psim = yinter(drini,mfm1,rmajor,rm,pinsid)
        b2m = yinter(drini,mfm1,rmajor,rm,b2ins)
      else
        psim = yinter(droti,mfm1,rmajor,rm,potsid)
        b2m = yinter(droti,mfm1,rmajor,rm,b2ots)
      end if
      if (codeid .eq. 'onedee') then
        xzm = ABS (rmajor-rm)*droti + 1.0
      else
        xzm = SQRT (psim-psiax)*drutpi + 1.0
      end if
      izm = xzm
      zetam = (pang+psim)*b2m/(xmass*v*rm)
      ism = -1
      if (zetai .ge. 0. .and. ipass .eq. 1) ism = 1
      if (ism*zetam .lt. -0.01 .or. izm .gt. izi+1) then
        flagm = '*'
        go to 340
      end if
c
c  calculate orbit width
c
      wid = ABS (rp-rmajor) - ABS (rm-rmajor)
      wid = MAX (wid, zero)
c
c  calculate weights by zone; assume ion spends equal time in each
c     zone it fully traverses and proportionately less time in the
c     innermost and outermost zones of the orbit
c
      if (izp .eq. izm)  go to 350
      wtm = 1.0 - (xzm - izm)
      wtp = xzp - izp
      wtx = 1.0/(wtm+wtp+izp-izm-1.0)
      wt(izm) = wtm*wtx
      wt(izp) = wtp*wtx
      if (izp .ne. izm+1) then
      do 310 i=izm+1,izp-1
  310 wt(i) = wtx
      end if
c
c  calculate pitch angle cosines by zone; assume cosine varies
c     linearly with zone index
c
      zdif      = zetap - zetam
      zeta(izm) = zetam + 0.5 * wt(izm)*zdif
      zeta(izp) = zetap - 0.5 * wt(izp)*zdif
c
      if (izp .ne. izm+1) then
        do i=izm+1,izp-1
          zeta(i) = zeta(i-1) + wt(i)*zdif
        end do
      end if
c
      go to 350
c
c  set error flag if error was detected
c
  340 izp = izi
      ier = 1
c
c  print out selected results
c
  350 if (norb .eq. 0)  return
      if (  ic .eq. 1)  write (norb, 1010)
      if (MOD (ic-1,iskip) .ne. 0)  return
      write (norb, 1020)  ic, ipass, iaxis, ism, isp, nin, not, ri, zi,
     .                    zetai, rm, flagm, zetam, rp, flagp, zetap,
     .                    wid, izi, izm, izp
      if (ic .ne. idebug(1) .and. ic .ne. idebug(2) .and.
     .    ic .ne. idebug(3) .and. ic .ne. idebug(4) .and.
     .                            ic .ne. idebug(5))  return
      do i=1,mf
        zinsid(i) = (pang+pinsid(i))*b2ins(i)/(xmass*v*rinsid(i))
        zotsid(i) = (pang+potsid(i))*b2ots(i)/(xmass*v*rotsid(i))
      end do
      do ii=mf,11,-10
        write (norb,2010) (i, i=ii,ii-9,-1)
        write (norb,2020) (rinsid(i), i=ii,ii-9,-1)
        write (norb,2030) (zinsid(i), i=ii,ii-9,-1)
        write (norb,2040) (finsid(i), i=ii,ii-9,-1)
      end do
      do ii=1,mfm1,10
        write (norb,2010) (i, i=ii,ii+9)
        write (norb,2020) (rotsid(i), i=ii,ii+9)
        write (norb,2030) (zotsid(i), i=ii,ii+9)
        write (norb,2040) (fotsid(i), i=ii,ii+9)
      end do
      return
c
 1010 format (/' orbit parameters' //
     .  4x,'ic  ip  ix ism isp nin not',7x,'ri',7x,'zi',4x,'zetai',
     .  7x,'rm',4x,'zetam',7x,'rp',4x,'zetap',5x,'wid','  izi izm izp')
 1020 format (i6,6i4,2f9.2,f9.4,2(f9.2,a1,f8.4),f8.2,1x,3i4)
 2010 format (/ ' i   :',10i12)
 2020 format (  ' r   :',10f12.2)
 2030 format (  ' zeta:',10f12.4)
 2040 format (  ' f   :',1p10e12.4)
c
      end

      subroutine pfit (p, x, y, xv, yv, nx, f, dfdx, dfdy)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension p(nx,*),x(*),y(*)
      dimension cx(4),cy(4),cxp(4),cyp(4), c2xp(4),c2yp(4)
c
c     bi-quadratic interpolator
c
c     input
c     1.  p     - effectivly a 4x4 matrix of function values
c                 but for generality the first dimension of
c                 p is nx.
c     2.  x     - associates a grid to the first dimension of p
c     3.  y     - associates a grid to the second dimension of p
c     4.  xv    - location at which bi-quadratic function is evaluated
c     5.  yv    - location at which bi-quadratic function is evaluated
c     6.  nx    - first dimension of p in calling program
c
c     output
c     1.  f     - interpolated value of p at (xv,yv)
c     2.  dfdx  - interpolated value of x-partial of p at (xv,yv)
c     3.  dfdy  - interpolated value of y-partial of p at (xv,yv)
c
c ----------------------------------------------------------------------
c
      a1    = (x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))
      a2    = (x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))
      a3    = (x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))
      a4    = (x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3))
c
      cx(1) = (xv-x(2))*(xv-x(3))*(xv-x(4))/a1
      cx(2) = (xv-x(1))*(xv-x(3))*(xv-x(4))/a2
      cx(3) = (xv-x(1))*(xv-x(2))*(xv-x(4))/a3
      cx(4) = (xv-x(1))*(xv-x(2))*(xv-x(3))/a4
c
      cxp(1) = ((xv-x(3))*(xv-x(4))
     .   +      (xv-x(2))*(xv-x(3))
     .   +      (xv-x(2))*(xv-x(4)))/a1
      cxp(2) = ((xv-x(3))*(xv-x(4))
     .   +      (xv-x(1))*(xv-x(3))
     .   +      (xv-x(1))*(xv-x(4)))/a2
      cxp(3) = ((xv-x(1))*(xv-x(2))
     .   +      (xv-x(1))*(xv-x(4))
     .   +      (xv-x(2))*(xv-x(4)))/a3
      cxp(4) = ((xv-x(1))*(xv-x(2))
     .   +      (xv-x(1))*(xv-x(3))
     .   +      (xv-x(2))*(xv-x(3)))/a4
c
      c2xp(1) = 2.0 * (3.0*xv-x(2)-x(3)-x(4))/a1
      c2xp(2) = 2.0 * (3.0*xv-x(1)-x(3)-x(4))/a2
      c2xp(3) = 2.0 * (3.0*xv-x(1)-x(2)-x(4))/a3
      c2xp(4) = 2.0 * (3.0*xv-x(1)-x(2)-x(3))/a4
c
      b1      = (y(1)-y(2))*(y(1)-y(3))*(y(1)-y(4))
      b2      = (y(2)-y(1))*(y(2)-y(3))*(y(2)-y(4))
      b3      = (y(3)-y(1))*(y(3)-y(2))*(y(3)-y(4))
      b4      = (y(4)-y(1))*(y(4)-y(2))*(y(4)-y(3))
c
      cy(1)   = (yv-y(2))*(yv-y(3))*(yv-y(4))/b1
      cy(2)   = (yv-y(1))*(yv-y(3))*(yv-y(4))/b2
      cy(3)   = (yv-y(1))*(yv-y(2))*(yv-y(4))/b3
      cy(4)   = (yv-y(1))*(yv-y(2))*(yv-y(3))/b4
c
      cyp(1) = ((yv-y(3))*(yv-y(4))
     .   +      (yv-y(2))*(yv-y(3))
     .   +      (yv-y(2))*(yv-y(4)))/b1
      cyp(2) = ((yv-y(3))*(yv-y(4))
     .   +      (yv-y(1))*(yv-y(3))
     .   +      (yv-y(1))*(yv-y(4)))/b2
      cyp(3) = ((yv-y(1))*(yv-y(2))
     .   +      (yv-y(1))*(yv-y(4))
     .   +      (yv-y(2))*(yv-y(4)))/b3
      cyp(4) = ((yv-y(1))*(yv-y(2))
     .   +      (yv-y(1))*(yv-y(3))
     .   +      (yv-y(2))*(yv-y(3)))/b4
      c2yp(1) = 2.0 * (3.0*yv-y(2)-y(3)-y(4))/b1
      c2yp(2) = 2.0 * (3.0*yv-y(1)-y(3)-y(4))/b2
      c2yp(3) = 2.0 * (3.0*yv-y(1)-y(2)-y(4))/b3
      c2yp(4) = 2.0 * (3.0*yv-y(1)-y(2)-y(3))/b4
c
      f      = 0.0
      dfdx   = 0.0
      dfdy   = 0.0
      d2fdx2 = 0.0
      d2fdy2 = 0.0
c
      do 10 i=1,4
      do 10 j=1,4
        f      = f + cx(i)*cy(j)*p(i,j)
        dfdx   = dfdx+cxp(i)*cy(j)*p(i,j)
        dfdy   = dfdy+cx(i)*cyp(j)*p(i,j)
        d2fdx2 = d2fdx2+c2xp(i)*cy(j)*p(i,j)
        d2fdy2 = d2fdy2+c2yp(j)*cx(i)*p(i,j)
   10 continue
      if (nx .gt. 0)  return
      dfdx = dfdx/d2fdx2
      dfdy = dfdy/d2fdy2
      return
c
      end

      subroutine polfit (ideg, npnts, x, y, coeff, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension x(101),y(101),coeff(7),sums(203), a(49),b(7),wkarea(7)
c
      ier    = 0
      idegp1 = ideg+1
      if ( ideg .ge. npnts+1)  ier = 1
      if (npnts .gt. 101    )  ier = 1
      if (  ier .eq. 1      )  return
c
      nsums = 2*ideg
      do i=0,nsums
        s = 0.0
        do j=1,npnts
          s = s+x(j)**i
        end do
        sums(i+1) = s
      end do
      k = 0
      do i=1,idegp1
        jstart = i
        jend   = i+idegp1-1
        do j=jstart,jend
          k      = k+1
          a(k)   = sums(j)
        end do
      end do
      do i=0,ideg
        s = 0.0
        do j=1,npnts
          s = s+y(j)*(x(j)**i)
        end do
        b(i+1) = s
      end do
c
c ----------------------------------------------------------------------
c solve set of simultaneous equations
c ----------------------------------------------------------------------
c
      call leqt1fl (a, 1, idegp1, idegp1, b, 0, wkarea, ier)
      do i=1,idegp1
        coeff(i) = b(i)
      end do
      return
c
      end

      real*8 function polval (xval, coeff, ideg)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension coeff(*)
c
      yval = coeff(ideg+1)
      do i=ideg+1,2,-1
        yval = yval*xval+coeff(i-1)
      end do
      polval = yval
      return
c
      end

      real*8 function poly3f (a, ndim1, ndim2, ndim3,
     .                        x1, n1, x2, n2, x3, n3)
c
c --- Modified by J. Mandrekas, for the SuperCode
c --- New version with new fits, 03/27/92
c
      implicit none
c
      integer ndim1, ndim2, ndim3, n1,n2, n3, i, i1, i2, i3
      real*8  a(ndim1,ndim2,ndim3)
      real*8  x1, x2, x3, y1(4), y2(4), y3(4)
c
      y1(1) = 1.0
      do i = 2, n1
        y1(i) = y1(i-1)*x1
      end do
c
      y2(1) = 1.0
      do i = 2, n2
        y2(i) = y2(i-1)*x2
      end do
c
      y3(1) = 1.0
      do i = 2, n3
        y3(i) = y3(i-1)*x3
      end do
c
      poly3f = 0.0
      do i1 = 1, n1
        do i2 = 1, n2
          do i3 = 1, n3
            poly3f = poly3f + a(i1,i2,i3)*y1(i1)*y2(i2)*y3(i3)
          end do
        end do
      end do
      return
c
      end

      real*8 function polyf (c, n, x)
c
      USE cpub_dat
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- evaluate the polynomial c(1)+c(2)*x+...+c(n)*x**(n-1).
c
****  real*8 c
c
c      common /cpub/ cpu1,cpu2,cpu3,cpu4,cpu5,cpu6,cpu7,cpu8,
c     .             cpuii, cpuiz, cpuplf
      dimension c(n)
c
      call SECOND (cpuaa)
      polyf = c(n)
      do j=1,n-1
        polyf  = polyf*x+c(n-j)
        polysv = polyf
      end do
      call SECOND (cpubb)
      cpuplf = cpuplf + cpubb - cpuaa
      return
c
      end

      subroutine postnub
c
      USE param
      USE solcon
      USE soln
      USE io 
      USE ions
      USE neut
      USE nub
      USE nub2
      USE mhdpar
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE tordlrot
      USE constnts
      USE rhog
      USE mcgo
      USE flxav
      USE gpsi
      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'storage.i'

c
c ----------------------------------------------------------------------
c this subroutine does the smoothing for the FREYA output
c ----------------------------------------------------------------------
c
      dimension  gfact(kf)
      dimension  ideg(4), nuse(4), nupdate(4)
      dimension  spbrt(kj),drho(kj)
      equivalence (ydum(1),drho(1))                              ! rsHSJ
      equivalence (xdum(1),spbrt(1))
c
c
c --- zero the FREYA calculated total beam torque:
c
      do j=1,nj
        spbrt (j) = 0.0
        spbolt(j) = 0.0 ! accumulates torque from orbit loss     ! rs
      end do
c
c ----------------------------------------------------------------------
c obtain r values (freyr) corresponding to psi values used in FREYA
c ----------------------------------------------------------------------
c
      call newgrid(psir,r,nj,psif,freyr,mf)
      call newgrid(r,hcap,nj,freyr,gfact,mf)
c
c ----------------------------------------------------------------------
c convert freyr to mid point values
c ----------------------------------------------------------------------
c
      do i=1,mfm1
        freyr(i) = 0.5 * (freyr(i) + freyr(i+1))
      end do
c
c ----------------------------------------------------------------------
c transform angmpz to transport coordinate system
c ----------------------------------------------------------------------
c
      do jbeam=1,nbeams
        do ie=1,ke
          call multpl1 (angmpz(1,ie,jbeam),mfm1,csgn)
        end do
      end do
c
c ----------------------------------------------------------------------
c smooth hibr and hdep
c note that if mcgo is used (iborb = 3) then hdep is not defined
c (only hibr is set in this case)
c ----------------------------------------------------------------------
c --- smoothing can lead to non physical profiles! as of 9/22/87 the
c --- smoothing can be turned on or off by input switch hdepsmth.
c --- see more hdepsmoothing further below
      if (hdepsmth .gt. 0.0)  go to 2111
      nsmooth = 2
      do i=1,4
        nuse(i)    = 10
        nupdate(i) = 10
        ideg(i)    = 2
      end do
c
      do 2110 jbeam=1,nbeams
      do 2110    ie=1,3
      do i=1,nsmooth
        call smooth (freyr,hibrz(1,ie,jbeam),mfm1,ideg(i),nuse(i),
     .               nupdate(i))
        call smooth (freyr,hdepz(1,ie,jbeam),mfm1,ideg(i),nuse(i),
     .               nupdate(i))
      end do
 2110 continue
c
c ----------------------------------------------------------------------
c convert hibr, hdep, zeta, ftrapfi and angmp from FREYA grid to transport grid
c ----------------------------------------------------------------------
c
c compute drho increments needed by orbit loss torque calc       ! rs
c                                                                ! rs
 2111 do 2130 j=1,nj            ! rs
      if (j .ne. 1)  go to 2125                                  ! rs
      drho(j) = r(2) / 2.0                                       ! rsHSJ
      go to 2130                                                 ! rs
 2125 if (j .ne. nj)  go to 2127                                 ! rs
      drho(j) = (r(nj)-r(nj-1)) / 2.0                            ! rs
      go to 2130                                                 ! rs
 2127 drho(j) = (r(j+1)-r(j-1)) / 2.0                            ! rs
 2130 continue                                                   ! rs
c                                                                ! rs
      do 2140 jbeam=1,nbeams
         if (iborb .eq. 3)
     .      call newgrid (freyr,beam_net_cur_mcgo(1,jbeam),mfm1,r,
     .              beam_net_cur_12_mcgo(1,jbeam),nj)
      do 2140 ie=1,3
      if(neg_ion_source(jbeam) .eq. 1 .and. ie .gt. 1)go to 2140
      call newgrid (freyr,hibrz(1,ie,jbeam),mfm1,r,
     .              hibr(1,ie,jbeam),nj) ! convert from zone to rho grid
c
      call newgrid (freyr,hdepz(1,ie,jbeam),mfm1,r,
     .              hdep(1,ie,jbeam),nj)
c
      if (mcgo_info_available .gt. 0) then ! get MCGO on transport grid
           call newgrid (freyr,enere_mcgo(1,ie,jbeam),mfm1,r,
     .              enere_12_mcgo(1,ie,jbeam),nj)
           call newgrid (freyr,eneri_mcgo(1,ie,jbeam),mfm1,r,
     .              eneri_12_mcgo(1,ie,jbeam),nj)
           call newgrid (freyr,density_mcgo(1,ie,jbeam),mfm1,r,
     .              density_12_mcgo(1,ie,jbeam),nj)
           call newgrid (freyr,prespar_mcgo(1,ie,jbeam),mfm1,r,
     .              prespar_12_mcgo(1,ie,jbeam),nj)
           call newgrid (freyr,presprp_mcgo(1,ie,jbeam),mfm1,r,
     .              presprp_12_mcgo(1,ie,jbeam),nj)
           call newgrid (freyr,dtneut_mcgo(1,ie,jbeam),mfm1,r,
     .              dtneut_12_mcgo(1,ie,jbeam),nj)
           call newgrid (freyr,ddneut_mcgo(1,ie,jbeam),mfm1,r,
     .              ddneut_12_mcgo(1,ie,jbeam),nj)
           call newgrid (freyr,curdens_mcgo(1,ie,jbeam),mfm1,r,
     .              curdens_12_mcgo(1,ie,jbeam),nj)
           call newgrid (freyr,fbth_mcgo(1,ie,jbeam),mfm1,r,
     .              fbth_12_mcgo(1,ie,jbeam),nj)
           call newgrid (freyr,fene_mcgo(1,ie,jbeam),mfm1,r,
     .              fbe_12_mcgo(1,ie,jbeam),nj)
           call newgrid (freyr,feni_mcgo(1,ie,jbeam),mfm1,r,
     .              fbi_12_mcgo(1,ie,jbeam),nj)
      end if
c
c --- do some smoothing HSJ
c
      if (hdepsmth .gt. 10 .and. hdepsmth .lt. 999. ) then
            npass = hdepsmth - 10
            nhdep = 6  !note that this introduces a grid dependence
            do ij=1,npass
                 hibr(1,ie,jbeam)=( 2.*hibr(1,ie,jbeam)+
     .                               hibr(2,ie,jbeam)+
     .                               hibr(3,ie,jbeam)+
     .                               hibr(4,ie,jbeam)+
     .                               hibr(5,ie,jbeam))/nhdep
                 hibr(2,ie,jbeam)=( hibr(1,ie,jbeam)+
     .                               2.*hibr(2,ie,jbeam)+
     .                               hibr(3,ie,jbeam)+
     .                               hibr(4,ie,jbeam)+
     .                               hibr(5,ie,jbeam))/nhdep
                 hdep(1,ie,jbeam)=( 2.*hdep(1,ie,jbeam)+
     .                               hdep(2,ie,jbeam)+
     .                               hdep(3,ie,jbeam)+
     .                               hdep(4,ie,jbeam)+
     .                               hdep(5,ie,jbeam))/nhdep
                 hdep(2,ie,jbeam)=( hdep(1,ie,jbeam)+
     .                               2.*hdep(2,ie,jbeam)+
     .                               hdep(3,ie,jbeam)+
     .                               hdep(4,ie,jbeam)+
     .                               hdep(5,ie,jbeam))/nhdep
               do j=3,nj-2
                  hibr(j,ie,jbeam)=(hibr(j-2,ie,jbeam)+
     .                               hibr(j-1,ie,jbeam)+
     .                               2.*hibr(j  ,ie,jbeam)+
     .                               hibr(j+1,ie,jbeam)+
     .                               hibr(j+2,ie,jbeam))/nhdep
                  hdep(j,ie,jbeam)=(hdep(j-2,ie,jbeam)+
     .                               hdep(j-1,ie,jbeam)+
     .                               2.*hdep(j  ,ie,jbeam)+
     .                               hdep(j+1,ie,jbeam)+
     .                               hdep(j+2,ie,jbeam))/nhdep
               end do
                 hibr(nj-1,ie,jbeam)= (hibr(nj-4,ie,jbeam)+
     .                               hibr(nj-3,ie,jbeam)+
     .                               hibr(nj-2,ie,jbeam)+
     .                               2.*hibr(nj-1,ie,jbeam)+
     .                               hibr(nj,ie,jbeam))/nhdep
                 hibr(nj,ie,jbeam)=(hibr(nj-4,ie,jbeam)+
     .                               hibr(nj-3,ie,jbeam)+
     .                               hibr(nj-2,ie,jbeam)+
     .                               hibr(nj-1,ie,jbeam)+
     .                               2.*hibr(nj,ie,jbeam))/nhdep
                 hdep(nj-1,ie,jbeam)=(hdep(nj-4,ie,jbeam)+
     .                               hdep(nj-3,ie,jbeam)+
     .                               hdep(nj-2,ie,jbeam)+
     .                               2.*hdep(nj-1,ie,jbeam)+
     .                               hdep(nj,ie,jbeam))/nhdep
                 hdep(nj,ie,jbeam)=(hdep(nj-4,ie,jbeam)+
     .                               hdep(nj-3,ie,jbeam)+
     .                               hdep(nj-2,ie,jbeam)+
     .                               hdep(nj-1,ie,jbeam)+
     .                               2.*hdep(nj,ie,jbeam))/nhdep
            end do
c
c renormalize the smoothed profiles to give the plasma volume when integrated
c
            vhdep =0.0
            vhibr =0.0
            call trapv (r, hdep(1,ie,jbeam), hcap, nj, vhdep)
            call trapv (r, hibr(1,ie,jbeam), hcap, nj, vhibr)
            vhdep=vhdep*volfac
            vhibr=vhibr*volfac
c it is possible for some beam/energy components to be zero
c because freya assignes zero particles to such a channel.
c (this could happend for example of the beam powers were greatly different)
c In this case renormalization should be skipped because the hdep,hibr are 
c identically zero:
         if(vhdep*vhibr .ne. 0.0)then
            if(iborb .ne. 3) then   !hdep is zero in this case, so vhdep =0
               do j=1,nj
                 hdep(j,ie,jbeam)=hdep(j,ie,jbeam)*(volume/vhdep)
                 hibr(j,ie,jbeam)=hibr(j,ie,jbeam)*(volume/vhibr)
               end do
            endif
         endif
      else if(hdepsmth .gt. 999)then
            npass = hdepsmth - 999 
            nhdep = 6  !note that this introduces a grid dependence
            do ij=1,npass
                 njsm = nj/2
                 if(ij .gt. npass-20)njsm = njsm+1
                 njsm = min(njsm,nj)
                 hibr(1,ie,jbeam)=( 2.*hibr(1,ie,jbeam)+
     .                               hibr(2,ie,jbeam)+
     .                               hibr(3,ie,jbeam)+
     .                               hibr(4,ie,jbeam)+
     .                               hibr(5,ie,jbeam))/nhdep
                 hibr(2,ie,jbeam)=( hibr(1,ie,jbeam)+
     .                               2.*hibr(2,ie,jbeam)+
     .                               hibr(3,ie,jbeam)+
     .                               hibr(4,ie,jbeam)+
     .                               hibr(5,ie,jbeam))/nhdep
                 hdep(1,ie,jbeam)=( 2.*hdep(1,ie,jbeam)+
     .                               hdep(2,ie,jbeam)+
     .                               hdep(3,ie,jbeam)+
     .                               hdep(4,ie,jbeam)+
     .                               hdep(5,ie,jbeam))/nhdep
                 hdep(2,ie,jbeam)=( hdep(1,ie,jbeam)+
     .                               2.*hdep(2,ie,jbeam)+
     .                               hdep(3,ie,jbeam)+
     .                               hdep(4,ie,jbeam)+
     .                               hdep(5,ie,jbeam))/nhdep
               do j=3,njsm
                  hibr(j,ie,jbeam)=(hibr(j-2,ie,jbeam)+
     .                               hibr(j-1,ie,jbeam)+
     .                               2.*hibr(j  ,ie,jbeam)+
     .                               hibr(j+1,ie,jbeam)+
     .                               hibr(j+2,ie,jbeam))/nhdep
                  hdep(j,ie,jbeam)=(hdep(j-2,ie,jbeam)+
     .                               hdep(j-1,ie,jbeam)+
     .                               2.*hdep(j  ,ie,jbeam)+
     .                               hdep(j+1,ie,jbeam)+
     .                               hdep(j+2,ie,jbeam))/nhdep
               end do

            end do
c
c renormalize the smoothed profiles to give the plasma volume when integrated
c
            call trapv (r, hdep(1,ie,jbeam), hcap, nj, vhdep)
            call trapv (r, hibr(1,ie,jbeam), hcap, nj, vhibr)
            vhdep=vhdep*volfac
            vhibr=vhibr*volfac
            do j=1,nj
               if(iborb .ne. 3) !hdep is zero in this case, so vhdep =0
     ,           hdep(j,ie,jbeam)=hdep(j,ie,jbeam)*(volume/vhdep)
                 hibr(j,ie,jbeam)=hibr(j,ie,jbeam)*(volume/vhibr)
            end do
      end if
c
c --- end smoothing
c
      call newgrid(freyr,zetaz(1,ie,jbeam),mfm1,r,
     .             zeta(1,ie,jbeam),nj)
      call newgrid(freyr,hicmz(1,ie,jbeam,1),mfm1,r,
     .             hicm(1,ie,jbeam,1),nj)
      call newgrid(freyr,hicmz(1,ie,jbeam,2),mfm1,r,
     .             hicm(1,ie,jbeam,2),nj)
      call newgrid(freyr,hicmz(1,ie,jbeam,3),mfm1,r,
     .             hicm(1,ie,jbeam,3),nj)
c
c --- negative zeta is possible. Undo the zeros set in newgrid if necessary:
c
      if (zetaz(mfm1,ie,jbeam) .lt. 0.0) then
        do j=nj,1,-1
          jj = j
          if (r(jj) .lt. freyr(mfm1))  go to 510
        end do
        write  (ncrt, 600)
        write  (nout, 600)
        write  (nout, 600)
  600   format (' subroutine POSTNUB detected error in zetaz intrp.' /)
        call STOP ('subroutine POSTNUB: unspecified problem', 47)
  510   jj = jj + 1
        slope = (zetaz(mfm1,ie,jbeam)-zetaz(mfm1-1,ie,jbeam))
     .        / (freyr(mfm1)-freyr(mfm1-1))
        do j=jj,nj
          drrj = r(jj)-freyr(mfm1)
          zeta(j,ie,jbeam) = zetaz(mfm1,ie,jbeam)+drrj*slope
        end do
      end if
c
      call newgrid(freyr,ftrapfi(1,ie,jbeam),mfm1,r,
     .             ftrapfit(1,ie,jbeam),nj)
      call newgrid(freyr,angmpz(1,ie,jbeam),mfm1,r,
     .             angmpf(1,ie,jbeam),nj)
      call newgrid(freyr,olossc(1,ie,jbeam),mfm1,r,              ! rs
     .             oloss(1,ie,jbeam),nj)                         ! rs
 2140 continue
c
c ----------------------------------------------------------------------
c renormalize hot ion birth rate and deposition rate
c ----------------------------------------------------------------------
c
      constv = 4.0 * pi**2 * rmajor
      do 2260 jb=1,nbeams
      do 2260 ie=1,3
      call trapv (r,hibr(1,ie,jb),hcap,nj,xbir)
      call trapv (r,hdep(1,ie,jb),hcap,nj,xdep)
      if (xbir .ne. 0.0) xbir = volume/(xbir*constv)
      if (xdep .ne. 0.0) xdep = volume/(xdep*constv)
      do 2250 j=1,nj
      hibr(j,ie,jb) = hibr(j,ie,jb)*xbir
      hdep(j,ie,jb) = hdep(j,ie,jb)*xdep
 2250 continue
 2260 continue
c
c ----------------------------------------------------------------------
c calculate particle, energy, and parallel momentum sources in plasma
c    due to neutral beam(s)
c      sb  is in particles/cm3-s
c      qb  is in w/cm3
c      spb is in g/cm2-s2
c      spbr (toroidal angular momentum source) is in g/cm-sec**2
c      spbolr (toroidal torque density from orbit loss current, gm/cm-sec   ! rs
c      pb0 (average initial parallel momentum per ion) is in gm/cm-sec
c ----------------------------------------------------------------------
c

      do 2280 jb=1,nbeams
      do 2280 ie=1,3
      xloss = fap(ie,jb) + fwall(ie,jb) + forb(ie,jb)
      vbeam = 1.384e6 * SQRT (1.0e3*ebeam(ie,jb)/atw_beam)
      do 2270 j=1,nj
      qb(j,ie,jb)   = (1.0-xloss)*pbeam(ie,jb)*hdep(j,ie,jb)/volume
      sb(j,ie,jb)   = 0.625e16*qb(j,ie,jb)/ebeam(ie,jb)
c
c --- sbpure(j,ie,jb) is set to the beam depostion sb here. sbpure is
c --- saved so that the modification due to time and charge exchange
c --- effects which are folded into sb in slow2 are not introduced.
c
      sbpure(j,ie,jb) = sb(j,ie,jb)                              ! HSJ
      spb(j,ie,jb)  = atw_beam*xmassp*vbeam*zeta(j,ie,jb)*sb(j,ie,jb)
****  spbr(j,ie,jb) = spb(j,ie,jb)*rcap(j)
      spbr(j,ie,jb) = angmpf(j,ie,jb)*sb(j,ie,jb)
      spbrt(j)      = spbrt(j) + spbr(j,ie,jb)
      bpflav = 0.0                                               ! rsHSJ
****  if (j .ne. 1)                                              ! rsHSJ
**** .  bpflav = rbp(j)/(fcap(j)*hcap(j)*gcap(j)*r(j)) / 1000.0  ! rsHSJ
      bpflav = bprmaj(j)                           ! bprmaj is in kgauss
      spbolr(j,ie,jb) = -oloss(j,ie,jb)*bpflav*rcap(j)*drho(j)   ! rsHSJ
      spbolr(j,ie,jb) = spbolr(j,ie,jb)*100.0                    ! rsHSJ
c
c --- dv = (4*pi**2) * R0 * hcap * r * dr = constv * hcap * r * dr ! HSJ
c --- assume oloss(amps),bpflav is flux averaged bp in (kgauss)  ! rs
c --- rcap(j) is (time-evolved) flux-averaged major radius (cm), ! rsHSJ
c --- drho has units  (cm),rbp is in gauss-cm                    ! rsHSJ
c --- spbolr (gm-cm/sec**2)                                      ! rs
c --- accumulate total orbit loss current torque                 ! rs
c       for routine source                                       ! rs
c
      spbolt(j)    = spbolt(j)+spbolr(j,ie,jb)                   ! rs
      pb0(j,ie,jb) = angmpf(j,ie,jb)*rcap(j)/r2capi(j)
 2270 continue
 2280 continue
c
c --- calculate the beam and orbit loss current torque in nt-m:  ! rs
c --- this is the prompt beam torque as determined directly from FREYA
c
      call trapv (r,spbrt,hcap,nj,stbeamt)
      stbeamt = stbeamt*1.0e-07*constv
      call trapv (r,spbolt,hcap,nj,stoloss)                      ! rs
      stoloss = stoloss*1.0E-07                                  ! rs
c
c ----------------------------------------------------------------------
c convert flux values to kgauss-cm**2
c ----------------------------------------------------------------------
c
      if (codeid .eq. 'onedee')  go to 2320
      do 2310 i=1,nw
      do 2310 j=1,nh
 2310   p(i,j) = 1.0e-3*p(i,j)
 2320 continue
      return
c
      end

      subroutine prenub (enbeams, mhdmethd)
c ----------------------------------------------------------------------
c
c this subroutine calculates certain flux surface information needed by FREYA
c if orbit effects are to be considered, then calculate the
c magnetic flux (pinsid and potsid) vs. major radius (rinsid and rotsid)
c through the magnetic axis.  if orbit effects are not being
c modeled, then simply store the flux and radius values at the
c magnetic axis and limiter in the appropriate arrays.  all flux
c and radius values are in CGS units here.
c
c --- input:
c
c      mhdmethd             character variable if mhdmethd='tdem'
c                           then we need to trace the plasma surface
c                           otherwise this switch is not used
c
c      INCLUDE file param:
c          kj
c          kf,kz
c
c      INCLUDE file bicube (used in 1-1/2-d cases only):
c          cspln(n2cspln,nw,nh2)
c          wnoperm(nwork)
c          pds(6)
c
c      INCLUDE file constnts:
c          pi
c          u0                          4 * pi * e-07
c
c      INCLUDE file flxav:
c          npsi
c          xmagn1
c          ymagn1        1-1/2d case only
c
c       INCLUDE file mhdpar:
c          nw
c          nh
c          kpsi
c
c       INCLUDE file mhdbcdtn:
c          mhdmultp
c
c      INCLUDE file mesh:
c          r(nj)                     used in 1d case only
c
c      INCLUDE file geom:
c          codeid
c          gcap(j)           used only if codeid = onedee
c
c      INCLUDE file mhdgrid (used in 1-1/2d cases only):
c          rmhdgrid
c          zmhdgrid
c
c      INCLUDE file io
c          ncrt
c          nout
c          nscr
c
c      INCLUDE file ions:
c          z(kj,kion)
c
c      INCLUDE file machin(1-d case only):
c          rmajor
c          rminor
c          btor
c          kappa
c
c      INCLUDE file numbrs:
c          nion
c          nj
c
c      INCLUDE file nub
c          mf,mfm1      #flux surfaces and #zones (mfm1 = mf-1) for FREYA grid
c                            (maximum possible mf is given by parameter kf)
c          ne_tk
c
c      INCLUDE file nub2:
c          iborb
c
c      INCLUDE file psig:
c          fpsi(npsi)         f(psi) defined over psival grid (1-1/2d only)
c          psival(npsi)
c          psivolp(npsi)       volume of flux zones on psival grid
c
c      INCLUDE file rhog:
c          psir(nj)                input in 1-1/2d case
c
c      INCLUDE file small:
c          p(nw,nh)         kgauss-cm**2 on input,changed on output (1-1/2d)
c
c      INCLUDE file soln:
c          en(j,i)           ion density,ion i,grid pt. j
c          ene(j)            electron density
c          rbp(nj)           used in 1-d case only
c          te(nj)
c          ti(nj)
c
c      INCLUDE file storage:
c          xdum(mf)
c          ydum(mf)
c          zdum(mf)
c          wdum(mf)         temporary storage vectors
c
c      INCLUDE file tordlrot:
c          angrot               toroidal rotation (on nj rho mesh)
c
c      INCLUDE file limiter:
c          xlimiter(1..nlimiter+2)
c          ylimiter(1...nlimiter+2)
c
c --- output---------------------------------------------------------------
c      INCLUDE file nub:
c          psif(mf)
c          psivol(mfm1)     volume of flux zones
c          zangrot(mfm1)        angrot on mfm1 grid
c          zne(mfm1)
c          zni(mfm1,i)
c          zte(mfm1)
c          zzi(mfm1,i)
c
c      INCLUDE file nub2:
c          b1ins(mf)
c          b2ins(mf)
c          b1ots(mf)
c          b2ots(mf)
c          pinsid(mf)
c          potsid(mf)
c          rinsid(mf)
c          rotsid(mf)
c
c      INCLUDE file nub3:
c          zti(mfm1)
c
c      INCLUDE file rhog:
c          psir(nj)              output for 1-d case
c
c      INCLUDE file small:
c          p(nw,nh)   for 1-1/2d case p is changed to gauss-cm**2 on output
c                       p(nw,nh) is changed back to kgauss-cm**2 in postnub
c
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE fusion
      USE ions
      USE io 
      USE solcon
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE nub
      USE nub2
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE tordlrot
      USE constnts
      USE nub3
      USE soln2d
      USE psig
      USE rhog
      USE bicube
      USE flxav
      USE neo2dp
      USE gpsi
      USE mhdbcdtn
      USE replace_imsl,                    ONLY : my_ibcccu,my_dbcevl1


      implicit  integer (i-n), real*8 (a-h, o-z)

      include 'storage.i'

c
      dimension     work(kf),zenbeam(kf),enbeams(*),zenbeamold(kf)
      dimension     powork(kj), rowork(kj) ,rxbp(kj)
      character*(*) mhdmethd
      data          icall_prenub /0/
c
c --- equivalences for 1-d run:
c
      equivalence (rxbp(1),xdum(1)),(rowork(1),ydum(1))
      equivalence (work(1),zdum(1)),(powork(1),wdum(1))
      equivalence (zenbeam(1),vdum(1))
c
      icall_prenub = icall_prenub + 1


      if (codeid .eq. 'onedee') then
c
c ----------------------------------------------------------------------
c calculate psif, an array of flux zone boundaries used in FREYA.
c for a 1-D run the zones are of equal width in minor radius.
c ----------------------------------------------------------------------
c

          do 10 j=1,nj
   10     psir(j) = r(j)
          dpsi = r(nj)/mfm1
          do 20 i=1,mf
   20     psif(i) = (i-1)*dpsi
c
c ----------------------------------------------------------------------
c calculate volumes of flux zones for elliptical cross sections
c ----------------------------------------------------------------------
c
          fact = 2.0 * pi**2 * rmajor*kappa*(rminor/mfm1)**2
          do 30 i=1,mfm1
   30     psivol(i) = fact*(i**2-(i-1)**2)
c
c ----------------------------------------------------------------------
c do calculation for elliptical cross sections
c ----------------------------------------------------------------------
c
          rotsid(1) = rmajor
          rotsid(mf) = rmajor + rminor
          if (iborb .ne. 0) then
              rkappa = SQRT (kappa)
              rowork(1) = 0.0
              rxbp(1) = 0.0
              work(1) = 1.0
              do 40 j=2,nj
              rowork(j) = r(j)/rkappa
              rxbp(j) = rmajor*rbp(j)/(r(j)*gcap(j))
   40         work(j) = SQRT (1.0 + (rkappa*rxbp(j)/(rmajor*btor))**2)
              call trap2(r,rxbp,powork,nj)
              drot = rminor/mfm1
              do 50 i=1,mfm1
   50         rotsid(i) = (i-1)*drot
              rotsid(mf) = rminor
              call intrp (0,1,rowork,powork,nj,rotsid,potsid,
     .                   mf)
              call intrp (0,1,rowork,work,nj,rotsid,b2ots,mf)
              do 60 i=1,mf
              rinsid(i) = rmajor - rotsid(i)
              rotsid(i) = rmajor + rotsid(i)
              pinsid(i) = potsid(i)
              b2ins(i) = b2ots(i)
              b1ins(i) = b2ots(i)
   60         b1ots(i) = b2ots(i)
          end if
      else        ! calculations for 1-1/2 d cases
c
c ----------------------------------------------------------------------
c calculate psif, an array of flux zone boundaries used in FREYA.
c     for a 1-1/2-D run the zones are of equal width in the square
c     root of the flux.  psif is in units of cm for 1-D runs and
c     kgauss-cm**2 for 1-1/2-D runs.
c     NOTE: square root makes flux zones near axis too small to get
c           accurate contours. square of flux makes zones near axis
c           too large. Therefore the best choice still is linear in psi.
c           the dropoff in the density and increase in volume near the
c           plasma boundary approximately lead to constant deposition.
c     NOTE: the square root spacinf is also used in FREYA so if it
c           is changed here it must be changed in FREYA also!
c     ALSO: The same is true for subroutine INJECT.
c ----------------------------------------------------------------------
c
          drutp = SQRT (psir(nj)-psir(1))/mfm1
          do 100 i=1,mf
  100     psif(i) = psir(1) + ((i-1)*drutp)**2
c
****      dpsi = (psir(nj)-psir(1))/mfm1
****      do 100 i=1,mf
**100     psif(i) = psir(1)+(i-1)*dpsi
c
c ----------------------------------------------------------------------
c  calculate volumes of flux zones for nonelliptical cross sections:
c           flux is in units of kgauss-cm**2; dimensions are in cm;
c ----------------------------------------------------------------------
c
          call revers (psif, work, mf) 
          call intrp (1, 1, psival, psivolp, npsi, psif, zdum, mf)
          do 140 i=1,mfm1
  140     psivol(i) = zdum(mf-i)-zdum(mf-i+1)
          call revers (psif, work, mf) 
c
c ----------------------------------------------------------------------
c convert flux values to gauss-cm**2
c ----------------------------------------------------------------------
c
          do 160 i=1,nw
          do 160 j=1,nh
  160     p(i,j) = 1.0e3*p  (i,j)
          psiax  = 1.0e3*psif( 1)
          psilim = 1.0e3*psif(mf)
c
c ----------------------------------------------------------------------
c get bicubic representation of psi ( = p):
c (reverse the sign of psi so that CNTOUR will work)
c ----------------------------------------------------------------------
c
          cconst = -1.0
          call multpl1 (p,nwh,cconst)
          call my_ibcccu  (p,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                     wnoperm,ier)
          call multpl1 (p,nwh,cconst)
c
c ----------------------------------------------------------------------
c generate the psi contours corresponding to the FREYA grid. these are
c copied into the neutral beam plotting file in subroutine NUBPLT.
c ----------------------------------------------------------------------
c
          call getioun(nscr,nscr)
          open (unit = nscr, file = 'scratch', status = 'UNKNOWN')
          bperr     = 0.05
          iconvg    = 0
          arcl      = 2.0    ! 2 cm arclength increment
          taxis     = 5.0
          tlim      = 30.0
          a         = (tlim-taxis)/(1000.0*(psif(1)-psif(mf)))
          bincp     = taxis
          delta_psi = 1000.0 * (-psif(mf-1)+psif(mf))
          do j=1,mfm1
             i = mf-j+1
             psi_psif = psif(i)
             if (ifixshap .eq. 1) then
c
c               the actual contour to be traced has the original (at
c               time t=time0) values of psi associated with the (R,Z)
c               grid IF the equilibrium is not evolved. Under this
c               condition the values of psi in psif will not map
c               back to the (R,Z) contours which are fixed in time.
c               Here we assume that the poloidal flux fraction inside each
c               (R,Z) contour remains invariant as a fucntion of time.
c               get normalized poloidal flux fraction at current time,psifn:
c
                psifn = (psi_psif -psif(1))/(psif(mf)-psif(1))
c
c               get psi_psif, the psi value that represents the same
c               normalized poloidal flux fraction as psifn but at the time the
c               eqdsk was read. Then get the (R,Z) contour points that
c               correspondt to psi_psif on the eqdsk  and ASSUME that
c               this contour corresponds to
c               the current value of psi, psif(i)
c
                psi_psif = psifn*(psir_at_t0(nj)-psir_at_t0(1)) +
     .                                           psir_at_t0(1)
             end if
c
             ptrace = -psi_psif * 1000.0
c
c --- prevent searching for plasma boundary by moving psilim in slightly:
c
c            if (j .eq. 1)  ptrace = -1005.0*psif(i)
             if (j .eq. 1)  ptrace = ptrace +0.01*delta_psi
c
             iauto  = 0
             iconvg = 0
             dang   = a*(ptrace+1000.0*psif(1))+bincp
             drx    = 0.0
             dry    = 0.0
             if (j .eq. 1 .and. mhdmethd .ne. 'tdem')  go to 501
             if (j .eq. 1 .and. mhdmethd .eq. 'tdem') then
               iauto  = 1
               iconvg = 0
             end if

             call cntour (xmagn1,ymagn1,ptrace,rcmin,rcmax,zcmin,
     .                    zcmax,zrcmin,zrcmax,rzcmin,rzcmax,dang,arcl,
     .                    bperr,drx,dry,100.0*xlimiter(nlimiter+1),
     .                    100.0*xlimiter(nlimiter+2),
     .                    100.0*ylimiter(nlimiter+1),
     .                    100.0*ylimiter(nlimiter+2),
     .                    iauto,iautoc,wdum,zdum,ndum,rmhdgrid,nw,
     .                    zmhdgrid,nh,cspln,n2cspln,nh2,iounit,kstore,
     .                    ierr,xdum,iconvg,delta_psi)

            if (ierr .ne. 0)
     .      call STOP ('subroutine PRENUB: returned CNTOUR error', 263)
            go to 502
  501       continue

              call fixedcntour (rplasbdry,zplasbdry,nplasbdry,
     .                        wdum,zdum,ndum,
     .                        rcmin,rcmax,zcmin,zcmax,
     .                        rzcmin,rzcmax,zrcmin,zrcmax,
     .                        rmhdgrid,zmhdgrid,nw,nh,
     .                        xdum,cspln,n2cspln,nh2,pds)
         call limiter_check(rcmin/100.,rcmax/100.,zcmin/100.,zcmax/100.,
     .                    xlimiter,ylimiter,nlimiter)
  502        write (nscr, '( i6   )')  ndum
             write (nscr, '(6e12.5)') (wdum(i),zdum(i),i=1,ndum)
c             get the flux zone volumes,udum((1) = total plasma volume,

             call volcalc (wdum,zdum,ndum,xmagn1,ymagn1,udum(j),area)
          end do
          udum(mf)   = 0.0     ! udum(mf) = volume at magnetic axis =0.0:
          volume     = 0.0
          do j=1,mfm1
            psivol(mfm1-j+1) = udum(j)-udum(j+1)
            volume = volume+psivol(mfm1-j+1) ! stored in machin.i
          end do
          ndum = 0
          write (nscr, '(i6)') ndum
          call giveupus(nscr)
          close (unit = nscr)
c
c ----------------------------------------------------------------------
c do calculation for nonelliptical cross sections
c get psi and grad psi for these grids. IBCCCU sets up bicubic
c spline array cspln, for psi values given in array p:
c ----------------------------------------------------------------------
c
          potsid(1)  = psiax
          potsid(mf) = psilim
          if (iborb .ne. 0) then
              call my_ibcccu (p,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                     wnoperm,ier)
c
c --- first get inside and outside major radius at elevation of magnetic axis
c
              isigncur = -1
              if (psiax .gt. psilim)  isigncur = 1
              npts = 3
              call getrmaj (ymagn1, npts, ierr, xdum(1), xdum(npts+1),
     .                      isigncur, psilim)
              rinsid(mf) = xdum(1)
              rotsid(mf) = xdum(npts)
c
c --- next get inboard and outboard (uniform) spatial grids
c --- at magnetic axis elevation
c
              drot = (rotsid(mf)-xmagn1)/mfm1
              drin = (xmagn1-rinsid(mf))/mfm1
              do 200 i=1,mfm1
                rotsid(i) = xmagn1+(i-1)*drot
  200         rinsid(i) = xmagn1+(1-i)*drin
              do 210 i=1,mf
                call my_dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                       rotsid(i),ymagn1,pds,ier,2)
                potsid(i) = pds(1)     ! psi     at (rotsid(i),ymagn1)
                b2ots(i) = pds(2)      ! dpsi/dr at (rotsid(i),ymagn1)
                call my_dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                       rinsid(i),ymagn1,pds,ier,2)
                pinsid(i) = pds(1)
  210         b2ins(i) = pds(2)
              b2ins(1) = 0.0
              b2ots(1) = 0.0
c
c --- given f(psi) ( = fpsi) on the psival grid,get f(psi) on the rotsid
c --- ( = fpsio) and rinsid (=fpsii) grids. The MHD grids are numbered
c --- from the plasma edge inward so we have to temporarily reverse them:
c
              do 410 i=1,npsi
  410         psival(i) = 1.0e3*psival(i)
              call revers(psival,work,npsi)
              call revers(fpsi,work,npsi)
              call intrp (1,1,psival,fpsi,npsi,potsid,fpsio,mf)
              call intrp (1,1,psival,fpsi,npsi,pinsid,fpsii,mf)
              do 420 i=1,npsi
  420         psival(i) = 1.0e-3*psival(i)
              call revers(psival,work,npsi)
              call revers(fpsi,work,npsi)
c
c --- next define the magnetic field ratios required in orbit:
c
              do 440 i=1,mf
              b1ins(i) = SQRT (fpsii(i)**2+b2ins(i)**2)
              b1ots(i) = SQRT (fpsio(i)**2+b2ots(i)**2)
              b2ins(i) = b1ins(i) / ABS (fpsii(i))
  440         b2ots(i) = b1ots(i) / ABS (fpsio(i))
              do 450 i=mf,1,-1
              b1ins(i) = b1ins(i)/b1ins(1)
  450         b1ots(i) = b1ots(i)/b1ots(1)
          end if
      end if
c
c ----------------------------------------------------------------------
c convert densities, temperatures, and angular rotation speed from
c   transport mesh points to FREYA zones for 1-d and 1-1/2d cases:
c ----------------------------------------------------------------------
c
      call tozone(psir,ene,nj,psif,zne,mfm1,nout,ncrt)
      call tozone(psir,te,nj,psif,zte,mfm1,nout,ncrt)
      call tozone(psir,ti,nj,psif,zti,mfm1,nout,ncrt)
      if (ne_tk .ne. 0)
     .call tozone(psir,angrot,nj,psif,zangrot,mfm1,nout,ncrt)
      do 600 i=1,nion
      call tozone(psir,en(1,i),nj,psif,zni(1,i),mfm1,nout,ncrt)
  600 call tozone(psir,z(1,i),nj,psif,zzi(1,i),mfm1,nout,ncrt)
c
c ----------------------------------------------------------------------
c Assume that the stored fast ion density presents a target with
c similar cross sections as the corresponding thermal density
c ----------------------------------------------------------------------
c
      if (fast_ion_target .gt. 0) then
         call tozone(psir,enbeams,nj,psif,zenbeam,mfm1,nout,ncrt)
c
c        add the beam density to the thermal density
c        that will be used by FREYA
c
         brelax = 0.5
         if (icall_prenub .gt. 1) then
         do i=1,mfm1
           zenbeam(i)    = brelax*zenbeam(i)+(1.-brelax)*zenbeamold(i)
           zenbeamold(i) = zenbeam(i)
         end do
         else
           zenbeamold(i) = 0.0
         end if
         if (ibion .gt. 0) then
           call adda  (zenbeam, zni(1,ibion), mfm1)
         else
           call addac (zenbeam, zni(1,id), mfm1,     fdbeam)
           call addac (zenbeam, zni(1,it), mfm1, 1.0-fdbeam)
         end if
      end if
c
c Put the Zeff array in an extra row of the zzi flux-zone array (dff)
c
      call tozone (psir,zeff,nj,psif,zzi(1,nion+1),mfm1,nout,ncrt)
c
c ----------------------------------------------------------------------
c transfrom to FREYA (lab) coordinate system.
c psi is assumed to correspond to plasma current in the positive
c (FREYA) toroidal system, so carry the sign in var:csgn.
c ----------------------------------------------------------------------
c
      csgn = mhdmultp
      if (ne_tk .ne. 0)  call multpl1 (zangrot,mfm1,csgn)

      return
c
      end

      subroutine prenub_pre (enbeams)
c
c ----------------------------------------------------------------------
c     use this subroutine to get mhdmethd into PRENUB without
c     including mhdcom in source.
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE mhdpar
      USE mhdcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c

c
      call prenub (enbeams, mhdmethd)

      return
c
      end





      subroutine pulse_start_time(nbeams,beamon,time0,timmax,
     .           btime,beamoff,tauslow,nsourc,source2_phase,
     .           kt,kj,kbs,kb,bstime,beam_end,
     .           pbeamOn,pbeamOff,bctime,n_pulse)
c-----------------------------------------------------------------------HSJ----
c  subroutine returns the starting time,bstime,required for multiple
c  beam times and the source switching times pbeamOn,pbeamOff.
c --- INPUT ---
c   time0          the start time of the analysis
c   nbeams         the number of active beams
c   beamon(i)      i=1,..nbeams  the start times of each beam
c   beam_end(i)                  end time
c   tauslow        the longest slowing down time taken over beams and regions
c                  the effect of tauslow is as follows. 
c                   a) if beamon(ib) > time0-tauslow then beamon(ib) is not
c                       changed from its input value.
c                   b) if beamon(ib) < time0-tauslow then beamon(ib)
c                      is reset to the largest value consistent with
c                      being one whole slowing down time below time0.
c                      For pulsed beams this is done by assuming that
c                      on input beamon(ib) gives the correct phase of the
c                      pulse relative to time0. Note that this means that
c                      you cannot simply set beamon(ib) =-99 UNLESS the
c                      beam has a pulse length longer than timmax-time0.
c                   
c   btime(i)    i=1,..neams the pulse on time
c   beamoff(i)  i=1,..nbeams (sec) the pulse off time
c   nsourc        number of sources per beamline (typically 2)
c   source2_phase(i) i=1,..neams is the phase in seconds of the turnon
c                  time of source 2 relative to source one for beam i.
c
c
c --- OUTPUT ---
c   bstime         the time at which the beam calculations are to start
c
c   beamon(ib)     ib=1,2,..nbeams the inital beam on time for
c                  source #1 of each beam. (Other sources of each
c                  beam have an inital turnon time of beamon(ib)+source2_phase(ib))
c                  This routine may 
c                  change beamon only if beamon < time0 on input.
c                 ( Even then beamon will be changed only if beamon .lt. time0-tuaslow)
c   pbeamOn(kt,kbs,kb)  the beam on/off pulse times in sec. These times
c   pbeamOff(kt,kbs,kb) are indexd by pulse number (kt),source number (kbs)
c                       and beam number (kb). pbeamOn(1,1,ib) is identical
c                       to beamon(ib)
c   n_pulse             the maximum number of pulses to be followed
c   timmax
c   time0           these times will be reset if the beam was on
c                   before time0
c----------------------------------------------------------------------HSJ-3/25/00--
      implicit none
      integer l,nbeams,nsourc,i,k,kt,kj,kbs,kb,kmax,n_pulse
      real *8 bstime,time0,beamon(kb),timestart,timmax,
     .        tauslow,tbon,tboff,beamoff(kb),btime(kb),
     .        phase,source2_phase(kb),first_pulse_off,
     .        pbeamOn(kt,kbs,kb), pbeamOff(kt,kbs,kb),timmax_save,
     .        bctime(*),beam_end(*)
      include "beam_plot_dat.i"
      bstime = timmax + 100.d0   !arbitrary initilization
      timestart = time0 -tauslow ! start at least  1 slowing down time back from time0
      do i=1,nsourc
            do l=1,nbeams
               tbon = beamon(l)      !beamon is time beam is first turned on
                                     !(NOT  length of time pulse is on)
               if(i == 2 )tbon = beamon(l) + source2_phase(l)
              
               if(tbon .lt. timestart)then
                   !here we may want to change beamon(l). 
                   !the start time, tbon, is back far enough
                   !to be asymtotic but the pulse duration may
                   !modify the situation:
                   first_pulse_off = tbon+btime(l)
                   if(first_pulse_off .gt. time0)then
                      !we are more than a slowing down time back from
                      !time0 and the pulse is not turned off before
                      !time0 so we can reset the start time to a slowing
                      !down time before time0:
                      tbon = timestart
                      if(i == 2)tbon = beamon(l) + source2_phase(l)
                   else
                      !here at least one full pulse exists
                      !before time0. We want to eliminate as many full
                      !cycles as possible while still maintaining a start
                      !time at least a slowing down time back of time0:
                      do while(tbon < timestart)
                        tboff = tbon + btime(l)
                        tbon  = tboff+beamoff(l)
                      enddo                     
                      tbon  = tbon - beamoff(l) - btime(l)
                   endif
                else if(tbon .lt. time0)then 
                      !here the beam is turned on before time0 but
                      !is not on long enough to be asymtotic(according
                      !to our estimate of tauslow).
                      !no change to tbon is necessary.
                else
                      !tbon .gt. time0,no change necessary
                endif

                bstime = MIN(bstime,tbon)         ! get the furthest back beam line
c                print *,'bstime,i= ',bstime,i
                if(i .eq. 1)beamon(l) = tbon
c                print *,'beamon(l),l = ',beamon(l),l
             enddo                                ! end loop over beams,l=1,...nb
      enddo                                       ! end i =1,..nsourc

c      if(bstime .gt. time0)bstime = time0    !why this ?? 


c     set up pbeamOn and pbeamOff
      pbeamOn =  timmax+100.d0      ! prevents unnecessary looping in beam_time
      pbeamOff = timmax+100.d0
      timmax_save = timmax

      if(bstime .lt. time0)then     !we are going to generate the startup file
           timmax = time0
           time0  = bstime  
           bctime(1) = time0        !use intial/bc at original time0
      endif





     
      do i=1,nsourc
         do l=1,nbeams
               k=0
               tbon = beamon(l)
               if(i == 2 )tbon = beamon(l) + source2_phase(l)
               do while (tbon <=  timmax_save)
                  tboff = tbon + btime(l)
                  if(tbon  >=  time0 .and. k < kt )then
                     k=k+1
                       if( tbon .le. beam_end(l))
     .                          pbeamOn(k,i,l)=tbon 
                       if( tboff .le. beam_end(l))     
     .                          pbeamOff(k,i,l)=tboff
                  elseif (k .eq. kt) then
                      call STOP('subroutine PULSE_Start_TIME:
     .  recompile with larger kt error #?', 0)
                  endif
                  tbon  = tboff+beamoff(l)
               enddo
         enddo
      enddo
      n_pulse = k   ! set to kt is ok, freya checks this
      do l=1,nbeams
         beamon(l) = pbeamon(1,1,l)
         print *,'beamon(l),l = ', beamon(l),l
         do i=1,nsourc
            do k=1,n_pulse
               print *,k,i,l
               print *,pbeamOn(k,i,l),pbeamOff(k,i,l)
            enddo
         enddo
      enddo


          s0_start = 0.0d0
         

c     following  used globally (beam_plot_dat.i)
      eps = ABS (.00001 * bstime )
      nbeam_pts=1
      timplot(nbeam_pts) =  bstime-eps
      waveform(nbeam_pts) = s0_start
      return
      end


      subroutine qgaus2a (func, a, b, ss)
c
      implicit  none
      integer j,nquad
      parameter (nquad = 5)
      real *8 x(nquad),w(nquad),xm,xr,ss,dx,a,b,func
c
      external    func
      data x
     .     /.1488743389,.4333953941,.6794095682,.8650633666,.9739065285/
      data w
     .     /.2955242247,.2692667193,.2190863625,.1494513491,.0666713443/
c
      xm = 0.5D0*(b+a)
      xr = 0.5D0*(b-a)
      ss = 0
      do j=1,nquad
        dx = xr*x(j)
        ss = ss+w(j)*(func(xm+dx)+func(xm-dx))
      end do
      ss = xr * ss
      return
c
      end

      real*8 function ranorm (fseed)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     generate one normal (0,1) pseudo-random number using a method
c     based on central limit theorem and power residue method.
c     output becomes truly normal as k (below) goes to infinity.
c     general formula for the deviate is:
c
c         y  =  (sum of x(i),i=1,k) -k/2.0) / SQRT (k/12.0)
c
c     where x(i) are are uniformly distributed on 0,1.
c     method is borrowed from IBM; they use k = 12
c
      external    RANDOM12                       ! random number generator
      parameter  (k = 12, rtkd12 = 1.0)
      data       seed0 /0.0/
c
****  rtkd12 = SQRT (k / 12.0)
      a      = 0.0
c
c --- note:  this loop is vector hazard
c
      do i=1,k
****    y = RANF (fseed)
*       y = RANF (     )
        y = RANDOM12 (seed0)
        a = a + y
      end do
      ranorm = (a - 0.5 * k) / rtkd12
      return
c
      end

      subroutine read_mcgo_files (ib, eb, ireb, bpow)
c
c ----------------------------------------------------------------------
c     This subroutine reads the MCGO-generated output
c --------------------------------------------------------9/9/98--HSJ---
c
      USE param
      USE io 
      USE nub
      USE nub2
      USE machin
      USE mcgo
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'machin.i'
c      include 'mcgo.i'
c      include 'nub.i'
c      include 'nub2.i'
c      include 'io.i'
c
      integer   onetwo_beam_mcgo
      dimension eb(ireb,*),bpow(ireb,*)
      real*8    zero
      real*8    kevjou
      data      kevjou /1.6022e-16 /
c
      zero = 0.0
c
****  do ib=1,nbeams
c
         if (read_mcgo_file(ib) .eq. 1) then
c
c          here the file name can be anything set by the user in inone:
           mcgo_output = mcgo_output_file(ib)
c
         else
c
c          input from a just spawned mcgo is expected to be in these files.
c          The user has not been given the opportunity to change the names.
c          the following must match mcgo file names. Links that automate
c          this are not in place at present!!! see mcgo include file
c          onetwo.i and subroutine mcgo::prep_onetwo
c
           if (ib .eq. 1)  mcgo_output = 'mcgo_12_output' // '_beam1'
           if (ib .eq. 2)  mcgo_output = 'mcgo_12_output' // '_beam2'
         end if
c
          call getioun(nmcgo3,nmcgo3)
          open (unit = nmcgo3, file = mcgo_output, status = 'OLD',
     .          form = 'UNFORMATTED', err = 100)
          read(nmcgo3)mfm1_mcgo,kep1_mcgo,onetwo_beam_mcgo
          read(nmcgo3)(rho_mcgo(j,ib),j=1,mf)
          read(nmcgo3)(psiz_mcgo(j,ib),j=1,mf)
          read(nmcgo3)(beam_net_cur_mcgo(j,ib),
     .                                   j=1,mf) ! mf contains total
          mf_mcgo=mfm1_mcgo+1
c
c         NOTE THAT INDEX MF CONTAINS TOTALS NORMALLY:
c
          read(nmcgo3)((prespar_mcgo(j,i,ib),j=1,mf_mcgo),
     .                                       i=1,ke) ! joules/cm3
          read(nmcgo3)((presprp_mcgo(j,i,ib),j=1,mf_mcgo),i=1,ke)
          read(nmcgo3)((density_mcgo(j,i,ib),j=1,mf_mcgo),
     .                                                    i=1,kep1_mcgo)
          read(nmcgo3)((curdens_mcgo(j,i,ib),j=1,mf_mcgo),
     .                                                   i=1,kep1_mcgo)
          read(nmcgo3)((dtneut_mcgo(j,i,ib),j=1,mf_mcgo),i=1,kep1_mcgo)
          read(nmcgo3)((ddneut_mcgo(j,i,ib),j=1,mf_mcgo),i=1,kep1_mcgo)
          read(nmcgo3)((enere_mcgo(j,i,ib),j=1,mf_mcgo),
     .                                     i=1,kep1_mcgo) ! watts/cm3
          read(nmcgo3)((eneri_mcgo(j,i,ib),j=1,mf_mcgo),i=1,kep1_mcgo)
          read(nmcgo3)((eng_cx_mcgo(j,i,ib),j=1,mf),i=1,kep1_mcgo)
          read(nmcgo3)((eng_forb_mcgo(j,i,ib),j=1,mf),i=1,kep1_mcgo)
          read(nmcgo3)((hibrz(j,ie,ib),j=1,mf),ie=1,ke)
          read(nmcgo3)((fbth_mcgo(j,ie,ib),j=1,mf),ie=1,ke)
c
c         reading bpow also loads pbeam since they refer to the same
c         storage locations:
c
          read(nmcgo3)(bpow(ie,ib),ie=1,ke)
          read(nmcgo3)(fap(ie,ib),ie=1,ke)
          read(nmcgo3)(fwall(ie,ib),ie=1,ke)
          read(nmcgo3)(bneut(ie,ib),ie=1,ke)
          read(nmcgo3)(forbit_mcgo(ie,ib),ie=1,ke)
          read(nmcgo3)(fpcx_mcgo(ie,ib),ie=1,ke)
          read(nmcgo3)((fene_mcgo(j,ie,ib),j=1,mf_mcgo),ie=1,ke)
          read(nmcgo3)((feni_mcgo(j,ie,ib),j=1,mf_mcgo),ie=1,ke)
          read(nmcgo3)(eb(ie,ib),ie=1,ke)
c
c         read in some initial points used in mcgo for plotting in nubplt:
c
          read(nmcgo3)(nparts_12(ie),ie=1,ke)
          k_12=0
          do j=1,ke
            k_12 = k_12+nparts_12(j)
          end do
          npts=k_12    ! in nub2.i
          read(nmcgo3)(rpts(j),j=1,k_12)
          read(nmcgo3)(zpts(j),j=1,k_12)
c
c         read phinit (which is in degrees) into ypts:
c
          read (nmcgo3) (ypts(j), j=1,k_12)
          do j=1,k_12
             angle_mcgo=ypts(j)*1.745329e-02
             xpts(j)=rpts(j)*cos(angle_mcgo) !cm
             ypts(j)=rpts(j)*sin(angle_mcgo)
          end do
c
      call giveupus(nmcgo3)
      close (unit = nmcgo3)
c
c     load some onetwo variables:
c
        call copya(psiz_mcgo(1,ib),psif,mf)
            do i=1,ke
               forb(i,ib)=forbit_mcgo(i,ib) !orbit losses assumed prompt
c
c              hdepz is normalized to the power that would go into
c              the thermal ion and electron distributions if no
c              power were lost due to charge exchange.
c              qb is derived from hdep (which comes from the zonal
c              hdepz) . The integral of qb,see variable pbeams,
c              is the fast ion power  in the plasma and it must
c              equal the power to the aperature,pbeam (or bpow),
c              times (1-fap-fwall -forbit) ,with a
c              correction due to beam slowing down effects if the
c              beam is not in steady state(as determined by the slow1-slow2
c              set of subroutines). If Mcgo is not run then the ions
c              are assumed to slow down according to the depostion profile
c              so the thermal source ends up being proportional to the
c              fast ion depostion shape . With Mcgo the thermal source will
c              have a different shape due to crossing of flux surfaces
c              during slowing down. Hence it is not possible to
c              fit the Mcgo results cleanly into the Freya framework.
c              The charge exchange  part,given by fpcx_mcgo,
c              is not used directly
c              in onetwo due to the Callen formulation. But the
c              the effect is implicit in enere_mcgo and eneri_mcgo.
c              Here we define and hdepz from the Mcgo calcs,which
c              is equivalent to the hdepz normally determined in Freya.
c              But the thermal deposition is based on enere_mcgo and
c              eneri_mcgo and eng_cx_mcgo,not the qb derived from hdepz!!
c
               hcon = bneut(i,ib)*(1.-fap(i,ib)-fwall(i,ib)-forb(i,ib))*
     .                eb(i,ib)*kevjou/volume
               if (hcon .ne. 0) then
                   do j=1,mf
c                    give hdep the thermal energy deposition shape:
c
****                 hdepz(j,i,ib)=(enere_mcgo(j,i,ib)+
**** .                 eneri_mcgo(j,i,ib)+eng_cx_mcgo(j,i,ib))/hcon
c
c                    or give hdep the fast ion denergy density shape:
c
                      hdepz(j,i,ib)=1.5*(0.666*presprp_mcgo(j,i,ib)+
     .                         0.333*prespar_mcgo(j,i,ib))/hcon
                   end do
               else
                  call multpl1 (hdepz(1,i,ib), mf, zero)
               end if
            end do
****  end do
c
      mcgo_info_available = 1
      return
c
  100 call STOP ('subroutine READ_MCGO_FILE: cannot read file', 272)
      end

      real*8 function reldif (x1, x2)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      del    = ABS (x1-x2)
      xabs   = ABS (x1   )
      if (xabs .gt. 1.0e-20)  del = del / xabs
      reldif = del
      return
c
      end

      subroutine revers (x ,work, n)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension  work(*), x(*)
c
      do i=1,n
        work(i) = x(i)
      end do
      k = 0
      do i=n,1,-1
        k    = k + 1
        x(k) = work(i)
      end do
      return
c
      end

      real*8 function rfine (b1ins, b2ins, c1, c2, c3, finsid, i, ileft,
     .                       mfm1, pinsid, rinsid, rmajor)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  This function uses a fine mesh to calculate:
c   a. the leftmost root of finsid within the interval from rinsid(i)
c      to rinsid(i-1) when ileft = 1 or
c   b. the rightmost root of finsid within the interval from rinsid(i+1)
c      to rinsid(i) when ileft = -1.
c
      dimension  b1ins(*), b2ins(*), finsid(*), pinsid(*), rinsid(*)
c
      r1 = rinsid(i)
      f1 = finsid(i)
      drin = rinsid(i-1)-rinsid(i)
      drini = 1.0/drin
      do 10 j=1,10
      r2 = r1 + ileft*0.1*drin
      psi2 = yinter(drini,mfm1,rmajor,r2,pinsid)
      b12 = yinter(drini,mfm1,rmajor,r2,b1ins)
      b22 = yinter(drini,mfm1,rmajor,r2,b2ins)
      f2 = r2**2 - c1*r2*b12 - (c2*b22*(c3+psi2))**2
      if (f2 .gt. 0.0)  go to 20
      r1 = r2
      f1 = f2
   10 continue
      rfine = 0.0
      return
   20 rfine = (r1*f2-r2*f1)/(f2-f1)
      return
c
      end

      subroutine rotate (naptr,iatype,aheigh,awidth,alen,bhofset,
     .                   bvofset,cangv,cangh,ib,isourc,costp,sintp,
     .                   costpp,sintpp,blenp,nsourc,sangv,sangh,
     .                   rpivot,zpivot,mlost,x0,y0,z0,vx0,vy0,vz0,
     .                   mlost1,mlost2)
c
      USE param,        ONLY : nap
      implicit  integer (i-n), real*8 (a-h, o-z)
      real,save :: iap
c
c  this subroutine advances a particle from source to the pivot point,
c  and transforms coordinates.
c
c     translation(s):  Particle to each aperture aligned along the souce
c                      axis.  Return with mlost = 1 if particle hits aperture.
c     rotations:       about y-axis so z-axis is aligned vertically;
c                      about new z-axis so x-axis parallel to beamline axis.
c     translation:     coordinates to intersection of beamline axis with
c                      source exit plane.
c     translation:     particles to each of beamline axis centered apertures,
c                      checking for mlost = 1 condition.
c     translation:     particle and then coordinate axis to pivot point.
c     rotation:        x-axis through pivot point and torus center.
c     translation:     origin to torus center.
c
c      parameter (nap = 10) 1%^&#@ 
      dimension iatype(nap,*), aheigh(nap,*), awidth(nap,*),
     .          alen(nap,*), bhofset(*), bvofset(*), blenp(*), cangv(*),
     .          cangh(*), sangv(*), sangh(*), rpivot(*), zpivot(*),
     .          costp(*), sintp(*), costpp(*),  sintpp(*)
c
c      if (nsourc .eq. 1)  go to 19 replaced with follwoing 4/1/2011 HSJ

c
cBH110314:  IF (nsourc .EQ. 1)  go to 19
cBH110314:  Define iap for nsourc=1 case
      IF(nsourc.eq.1)  THEN
         iap = 0
         go to 19
      ENDIF
c
c     Move particle to each of source-axis centered apertures,
c     and test for particle passage.
c
      alen1 = 0.0
      alen2 = 0.0
      do 10 i=1,naptr
c
c     iatype .le. 4 are source-centered apertures
c

      if (iatype(i,ib) .gt. 4)  go to 11
      alen1 = alen(i,ib)-alen2
      alen2 = alen(i,ib)
      tpvt = -alen1/vx0
      x0 = x0+vx0*tpvt
      y0 = y0+vy0*tpvt
      z0 = z0+vz0*tpvt
c
      mlost = 0
      go to (12,13,14,15),  iatype(i,ib)
   12 if ((y0**2+z0**2) .le. (0.5 * awidth(i,ib))**2)  go to 10
      go to 16
   13 if (ABS (y0) .le. 0.5 * awidth(i,ib) .and.
     .    ABS (z0) .le. 0.5 * aheigh(i,ib))  go to 10
      go to 16
   14 if (ABS (z0) .le. 0.5 * aheigh(i,ib))  go to 10
      go to 16
   15 if (ABS (y0) .le. 0.5 * awidth(i,ib))  go to 10
   16 mlost = 1
      mlost1 = mlost1 +1
      return
   10 iap = i
   11 continue


c
c     Source center is taken to lie at vertical distance
c     bvofset and horiz. distance bhofset from beam line axis.
c     Source centerline intesects beam line axis at distance
c     bleni from source.  Coords. are changed from source
c     coordinates to those with x-axis along beamline axis
c     directed in +R direction, z in the vertical direction, and
c     origin in source exit plane.
c
c     Rotate about y-axis of source coordinate system so x-axis
c     is vertical and in exit plane defined by beam axis
c     and source centers.
      temp = x0*costp(ib)-isourc*z0*sintp(ib)
      z0 = isourc*x0*sintp(ib)+z0*costp(ib)
      x0 = temp
      temp = vx0*costp(ib)-isourc*vz0*sintp(ib)
      vz0 = isourc*vx0*sintp(ib)+vz0*costp(ib)
      vx0 = temp
c
c     Rotate about z-axis of source system so y-axis lies
c     in source exit plane.
      temp = x0*costpp(ib)-isourc*y0*sintpp(ib)
      y0 = isourc*x0*sintpp(ib)+y0*costpp(ib)
      x0 = temp
      temp = vx0*costpp(ib)-isourc*vy0*sintpp(ib)
      vy0 = isourc*vx0*sintpp(ib)+vy0*costpp(ib)
      vx0 = temp
c
c     Translate coordinate axes to beamline axis.
      x0 = x0
      y0 = y0+isourc*bhofset(ib)
      z0 = z0+isourc*bvofset(ib)
c
c     Translate particle to beamline centered apertures and set
c     mlost = 1 if particle lost.
c
   19 continue
      if (iap .eq. naptr)  go to 29
      alen1 = 0.0
      alen2 = ABS (x0)
      do 20 i=iap+1,naptr
      alen1 = alen(i,ib)-alen2
      if (alen1 .lt. 0.0)
     .  call STOP ('subroutine ROTATE: unspecified problem', 48)
      alen2 = alen(i,ib)
      tpvt = -alen1/vx0
      x0 = x0+vx0*tpvt
      y0 = y0+vy0*tpvt
      z0 = z0+vz0*tpvt
      mlost = 0
      go to (22,23,24,25,26),  iatype(i,ib)-4
   22 if ((y0**2+z0**2) .le. (0.5 * awidth(i,ib))**2)  go to 20
      go to 27
   23 if (ABS (y0) .le. 0.5 * awidth(i,ib) .and.
     .    ABS (z0) .le. 0.5 * aheigh(i,ib))  go to 20
      go to 27
   24 if (ABS (z0) .le. 0.5 * aheigh(i,ib))  go to 20
      go to 27
   25 if (ABS (y0) .le. 0.5 * awidth(i,ib))  go to 20
      go to 27
c     DIII-D  special case:
   26 if (ABS (z0) .gt. 24.75)  go to 27
      if (ABS (y0) .le. 13.45)  go to 20
      if (ABS (y0) .gt. 21.7)  go to 27
      if (ABS (z0) .gt. 24.75*((ABS (y0)-13.45)/(-8.25)+1))  go to 27
      go to 20
   27 mlost = 1
      mlost2 = mlost2+1 
      return
   20 continue
   29 continue
c
c  translate particle and then coordinate axes to pivot point
c
      tpvt = -(blenp(ib)+x0)/vx0
      x0 = 0.0
      y0 = y0+vy0*tpvt
      z0 = z0+vz0*tpvt
c
      if (sangv(ib) .eq. 0.0)  go to 30
      zcos = cangv(ib)
      zsin = sangv(ib)
      temp = x0*zcos+z0*zsin
      z0 = -x0*zsin+z0*zcos
      x0 = temp
      temp = vx0*zcos+vz0*zsin
      vz0 = -vx0*zsin+vz0*zcos
      vx0 = temp
   30 continue
c
      if (sangh(ib) .eq. 0.0)  go to 40
      zcos = cangh(ib)
      zsin = sangh(ib)
      temp = x0*zcos+y0*zsin
      y0 = -x0*zsin+y0*zcos
      x0 = temp
      temp = vx0*zcos+vy0*zsin
      vy0 = -vx0*zsin+vy0*zcos
      vx0 = temp
   40 continue
c
      x0 = x0 + rpivot(ib)
      z0 = z0 + zpivot(ib)
      return
c
      end

      real*8 function scxif (i, er)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- estimated cross sections for charge exchange
c
      common /b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
c
      if (i .ge. 2)  go to 200
c
c --- 1s
c
      aer = LOG (er)
c
c     charge exchange
c       below 1.0e6 eV, use Riviere's fit
c       above 1.0e6 eV, use power-law extrapolation.
c
      if (er .gt. 1.0e6)  go to 105
      siif1 = (0.6937e-14) * (1.0 - 0.06732   *aer)**2
     .                     / (1.0 + 0.1112e-14*er  **3.3)
      go to 110
  105 siif1 = (4.8363e-22)/((1.0e-6)*er)**3.3
  110 scxif = siif1
      return
c
  200 continue
      ansq   = (nfhx(i))**2
      xtilde = ansq*er / 9260.0
      scxif  = 0.0
      if (xtilde .lt. 1.0)  scxif = 1.465e-11 * xtilde * (ansq/er)
      return
c
      end

      real*8 function scxzf (i, z, er)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- cross section for electron loss via collision with impurity.
c
      common /b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
c
      ni = nfhx(i)
      id = 2
      if (ni .eq. 1)  id = 1
c
      ansq   = (FLOAT (ni))**2
      xtilde = ansq * er / (31250.0 * z)
      scxzf  = 0.0
      if (xtilde .lt. 1.0)  scxzf = 1.465e-11 * xtilde * (ansq*z*z/er)
      return
c
      end

      real*8 function sdaccf (ni, nj, z, er)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- dacc cross section.
c     Janev & Presnyakov, j. phys. b 13, 4233 (1980).
c
      parameter (ms = 21, mc = 35)
      common /b3/ f(ms,mc), ar(ms+1,ms+1)
c
      omega  = 0.5 * (1.0 / ((FLOAT (ni))**2) - 1.0 / ((FLOAT (nj))**2))
      alam   = SQRT (f(ni,nj)/(2.0*omega))
      beta   = z*alam*omega/(er/24982.0)
      sdaccf = (1.7595e-16)*(z*alam/omega)*dfhx(beta)
      return
c
      end

      real*8 function seef (i, j, er)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- cross sections for excitation from i to j due to electron impact
c
      parameter (ms = 21, mc = 35)
      common /b1/kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      common /b2/nouthx,istart,ihxbug
      common /b4/en(mc+1),dg(mc+1),ae(ms,mc),be(ms,mc),
     .          de1(ms,mc),de2(ms,mc),ge1(ms,mc),ge2(ms,mc)
      data ryd /13.6/
c
      real*8 polyf
c
c --- data for 1s to 2s (e impact):
c
      data emin12/10.2/, emax12/1140.5/
      dimension ae12(8)
      data ae12/
     . -1.5138513657e-17, 1.8231019028e-16,-1.9274654811e-16,
     .  8.9291530932e-17,-2.2108041618e-17, 3.0448025264e-18,
     . -2.2039298214e-19, 6.5491238735e-21/
c
c --- data for 1s to 2p (e impact)
c
      data emin13/10.2/, emax13/747.4/
      dimension ae13(9)
      data ae13/
     . -2.1078372920e-14, 3.7548968459e-14,-2.8788463637e-14,
     .  1.2394689683e-14,-3.2696005067e-15, 5.4068250737e-16,
     . -5.4759059881e-17, 3.1084887079e-18,-7.5815695055e-20/
c
      zero = 0.0
      if (i .ge. j)  go to 300
      seef = 0.0
      aer  = LOG (er)
      if ((i .ge. 3) .or. (j .ge. 4))  go to 250
      go to (210, 240), i
c
  210 go to (300, 220, 230), j
c
c --- 1s to 2s:
c
  220 if (er .le. emin12)  return
      if (er .le. 11.1  )  go to 221
      if (er .ge. emax12)  go to 222
      seef = polyf (ae12, 8 ,aer)
      go to 223
c
c --- linear interpolation for region just above threshold:
c
  221 seef = ((er - 10.2) / 0.9) * 1.6297e-17
      go to 223
c
c --- asymptotic:
c
  222 seef = (5.312e-16) / er
  223 return
c
c --- 1s to 2p:
c
  230 if (er .le. emin13)  return
      if (er .ge. emax13)  go to 232
      seef = polyf (ae13, 9, aer)
      go to 233
c
c --- asymptotic:
c
  232 seef = (2.65571e-15/er)*(aer-2.120681)
  233 return
c
c --- 2s to 2p:
c
  240 seef = (8.617e-15/er)*LOG (1.14e4*er)
      return
c
  250 ni = nfhx(i)
      nj = j - 1
c
c --- from Vriens & Smeets (eq.(14)):
c
      seef = (1.75947e-16*ryd)*(ae(ni,nj)*LOG (0.5 * er/ryd+de1(ni,nj))
     .     + be(ni,nj))/(er+ge1(ni,nj))
      seef = MAX (seef, zero)
      id   = 2
      if (ni .eq. 1)  id = 1
      return
c
  300 ihxbug = 6
      if (nouthx .gt. 0) then
        write (nouthx, 3939)  i, j
 3939   format (' ERROR in SEEF: i4 and j are in error i4 = ', i3 /
     .          '                                       j = ', i3)
      end if
      return
c
      end

      subroutine setrz (ndim, rmin, rmax, zmin, zmax, dr, dz, nr, nz)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this subroutine sets up (r,z) grid parameters for FREYA when the
c optional two-dimensional deposition calculation is used.
c default parameters correspond to original eqdsk grid.
c ----------------------------------------------------------------------
c
      dimension ndim(2)
c
      nr =  ndim(1)
      nz =  ndim(2)
      dr = (rmax-rmin)/(nr-1)
      dz = (zmax-zmin)/(nz-1)
      return
c
      end

      real*8 function sezf (i, j, z, er)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- cross sections for excitation from i to j due to ion impact
c
      common /b1/ kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      common /b2/ nouthx, istart, ihxbug
c
      real*8 polyf
c
c --- data for 1s to 2s (p impact):
c
      data      emin12/1.0e3/, emax12/1.0e6/
      dimension ai12(11)
      data      ai12/
     .  -3.3812171941e+05, 3.5976726159e+05,-1.7037301759e+05,
     .   4.7282644608e+04,-8.5167211591e+03, 1.0405267504e+03,
     .  -8.7343818297e+01, 4.9754307346e+00,-1.8412279745e-01,
     .   3.9984716057e-03,-3.8707719188e-05/
c
c --- data for 1s to 2p (p impact):
c
      data      emin13/1.0e3/, emax13/1.0e6/
      dimension ai13(12)
      data      ai13/
     .  -1.3490069287e+06, 1.4573274888e+06,-7.0815407013e+05,
     .   2.0432442357e+05,-3.8900004212e+04, 5.1319758650e+03,
     .  -4.7884757131e+02, 3.1608287540e+01,-1.4469255104e+00,
     .   4.3759824250e-02,-7.8716881911e-04, 6.3824434435e-06/
c
      zero = 0.0
      if (i .ge. j) then
        ihxbug = 4
        if (nouthx .gt. 0) then
          write  (nouthx, 100)  i, j
  100     format (' ERROR in SEZF: i4 = ', i3, '    j = ', i3)
        end if
        sezf = 0.0
        return
      end if
c
      if (z .gt. 1.01)  go to 300
      aer = LOG (er)
      if ((i .ge. 3) .or. (j .ge. 4))  go to 250
      go to (210, 240), i
c
  210 go to (400, 220, 230), j
c
c --- 1s to 2s:
c
  220 if (er .ge. emin12)  go to 221
      sezf = 2.8875e-18*(er/emin12)**0.7093
      go to 223
c
  221 if (er .ge. emax12)  go to 222
      sezf = EXP (polyf(ai12,11,aer))
      go to 223
c
  222 sezf     = 1.9564e-12 / er
  223 sezfsave = sezf
      return
c
c --- 1s to 2p:
c
  230 if (er .ge. emin13)  go to 231
      sezf = 1.4053e-17*(er/emin13)**0.95695
      go to 233
c
  231 if (er .ge. emax13)  go to 232
      sezf = EXP (polyf(ai13,12,aer))
      go to 233
c
  232 sezf     = (6.7085e-13/er) * LOG (5701.79*er)
  233 sezfsave = sezf
      return
c
c --- 2s to 2p:
c
  240 sezf     = (1.584e-10/er) * LOG (0.62*er)
      sezf     = MAX (sezf, zero)
      sezfsave = sezf
      return
c
  250 ni       = nfhx(i)
      id       = 2
      if (ni .eq. 1)  id = 1
      nj       = j-1
      zz1      = 1.0
      sezf     = sdaccf(ni,nj,zz1,er)
      sezfsave = sezf
      return
c
c --- impurity scattering:
c
  300 if (i .eq. 1 .and. j .eq. 2)  go to 312
      if (i .eq. 1 .and. j .eq. 3)  go to 313
      sezf     = 0.0
      sezfsave = sezf
      if ((i .eq. 2) .and. (j .eq. 3))  return
      ni       = nfhx(i)
      nj       = j - 1
      id       = 2
      if (ni .eq. 1)  id = 1
      sezf     = sdaccf(ni,nj,z,er)
      sezfsave = sezf
      return
c
  312 sezf     = 0.25 * sdaccf(1,2,z,er)
      sezfsave = sezf
      return
c
  313 sezf     = 0.75 * sdaccf(1,2,z,er)
      sezfsave = sezf
      return
c
  400 ihxbug = 5
      if (nouthx .gt. 0) then
        write (nouthx, '(/ a /)')  ' ERROR in SEZF: state "j" is 0'
      end if
      sezf = 0.0
      return
c
      end

      real*8 function sief (i, er)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- cross sections for ionization due to e impact
c
      common /b1/ kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
      data        ryd/13.6/
      real*8      polyf
c
c --- data for e-impact ionization of 1s
c
      data      emin1/13.723/, emax1/620.0/
      dimension ae1(7)
      data      ae1/
     .          2.3036652148e-15,-3.4298197580e-15, 1.9465132340e-15,
     .         -5.4457519508e-16, 8.0929651995e-17,-6.1527279210e-18,
     .          1.8889193736e-19/
c
c --- data for e-impact ionization of 2s
c
      data emin2/3.4/, emax2/386.0/
      dimension ae2(9)
      data      ae2/
     .          8.4162370468e-15,-2.3545910382e-14, 2.4728342689e-14,
     .         -1.2718639429e-14, 3.6382722533e-15,-6.0044378263e-16,
     .          5.5114285405e-17,-2.4253587346e-18, 3.0573876445e-20/
c
c --- bea cross section (ansq = n**2, ery = e/ryd):
c
      sbea(ansq,ery) =
     .  (3.519e-16)*(5.0 * ansq/3.0-1.0/ery-2.0/(3.0*ansq*ery*ery))
     .  /(ery+3.25/ansq)
c
      sief = 0.0
      if (i .ge. 3)  go to 130
      aer  = LOG (er)
      go to (110, 120), i
c
c --- 1s
c
  110 if (er .le. emin1)  return
      if (er .ge. emax1)  go to 112
      sief = polyf(ae1,7,aer)
      go to 113
  112 sief = (1.3563e-15/er)*(aer+1.823647)
  113 return
c
c --- 2s
c
  120 if (er .le. emin2)  return
      if (er .ge. emax2)  go to 122
      sief = polyf(ae2,9,aer)
      go to 123
  122 sief = (8.195137e-15/er) * (aer-0.9445204)
  123 return
c
c --- 2p and higher
c
  130 ansq = (FLOAT (i-1))**2
      if (er .le. ryd/ansq)  return
      sief = sbea(ansq,er/ryd)
      return
c
      end

      subroutine sigfit (ebeam, dene, te, denz, sig)
c
c ----------------------------------------------------------------------
c.....author:
c       Charles D. Boley
c       Lawrence Livermore National Laboratory
c       L-574
c       Livermore, CA 94550
c       (415) 423-7365
c
c     version of 3/11/92
c
c.....This subroutine computes the neutral beam stopping cross section,
c     for a hydrogenic beam, as a function of beam energy, electron
c     density, electron temperature, and impurity density.  Four
c     impurities (He, C, O, Fe) are allowed.
c
c     The code employs a fit based on results of the author's beam
c     penetration code.  The latter code contains the detailed atomic
c     cross sections recommended at the IAEA Specialists' Meeting on
c     Required Atomic Data Base for Neutral Beam Penetration in Large
c     Tokamaks (Vienna, 4/89).  The fit applies for beam energies from
c     10 keV/amu to 10 MeV/amu and for all densities and electron
c     temperatures of interest.
c
c     The fit is independent of the ion temperature and hence does not
c     distinguish among hydrogen ion isotopes in the plasma.  (The actual
c     cross sections were evaluated at Ti = 15 keV.)  The ion temperature
c     has little effect provided that it is small compared to the beam
c     energy per amu.  (This condition normally is a triviality, except
c     perhaps for the 1/3 energy component of a deuterium beam.)  The
c     following table gives an idea of the variation with ion temperature
c     at low beam energy, for injection into a pure H plasma with
c     dene=1.e14 and Te=10 keV (energies in keV; cross sections in this
c     table in units of 10**-16 cm**2):
c
c     ebeam  sig(Ti=1)  sig(Ti=15)  sig(Ti=30)  sig(fit)
c       10     11.43       9.96        8.91      11.18
c       30      6.06       5.12        4.71       5.35
c       50      3.57       3.56        3.39       3.75
c      100      2.09       2.11        2.09       2.28
c
c     The fit was evaluated with B(perpendicular component) = 5 T.
c     Variations of the magnetic field have an insignificant effect on
c     the stopping cross section.
c
c     Accuracy of fit: rms error about 2.5%
c                      max error about 12%
c
c     Input parameters --
c       ebeam: beam energy per amu (keV/amu) -- 10. .le. ebeam .le. 1.e4
c       dene:  electron density (cm**-3) -- 1.e12 .le. dene .le. 1.e15
c       te:    electron temperature (keV) -- te .ge. 1.
c       denz:  impurity densities (cm**-3) -- array of dimension 4
c         denz(1): He
c         denz(2): C
c         denz(3): O
c         denz(4): Fe
c         Note: denz(i)=0. is permissible
c     Output parameter --
c       sig:   beam stopping cross section (cm**2)
c
c     Modified by John Mandrekas (GIT) for the SUPERCODE
c     New version with fits valid from lower beam energies (10 keV)
c     03/27/92, John Mandrekas, GIT
c
      implicit none
c
      integer mz, mt1, mt2, mt3
      parameter(mz=4, mt1=4, mt2=3, mt3=2)
      integer i, iinit
c
      integer       nth1, nth2, nth3, ntz1(mz), ntz2(mz), ntz3(mz)
      real*8        A1(mt1,mt2,mt3), A2(mt1,mt2,mt3,mz)
      common /cfit/ A1, A2, nth1, nth2, nth3, ntz1, ntz2, ntz3
c
      real*8 s1, ebeam, te, sig, sum, dene, poly3f,
     .       denz(mz), s2(mz), az(mz), xx(3)
      data az/2.,6.,8.,26./
      data iinit/0/
c
      if (iinit.eq.0) then
        iinit = 1
        call initfit
      end if
c
      xx(1) = LOG (ebeam)
      xx(2) = LOG (dene/1.e13)
      xx(3) = LOG (te)
      s1 = poly3f(A1,mt1,mt2,mt3,xx(1),nth1,xx(2),nth2,xx(3),nth3)
      do i = 1, 4
         if (denz(i).gt.0.) then
            s2(i)=poly3f(A2(1,1,1,i),mt1,mt2,mt3,xx(1),ntz1(i),
     .            xx(2),ntz2(i),xx(3),ntz3(i))
         else
           s2(i)=0.
         end if
      end do
c
      sum=0.0
      do i = 1, 4
         sum = sum + (denz(i)*az(i)*(az(i)-1.)/dene)*s2(i)
      end do
c
      sig = (1.e-16/ebeam)*exp(s1)*(1.+sum)
      return
c
      end

      real*8 function siif (i, er)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- cross sections for charge exchange and ion-impact ionization
c
      real*8      polyf
      common /b1/ kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
c
c --- data for ion-impact ionization of h(1s)
c
      dimension ai1(6)
      data      ai1/ -4.410635e+02,  1.846170e+02, -3.429509e+01,
     .                3.217263e+00, -1.512004e-01,  2.825854e-03/
c
c --- data for ion-impact ionization of h(1s) (Freeman & Jones)
c
      if (i .ge. 2)  go to 200
c
c --- 1s:
c
      aer = LOG (er)
c
c     charge exchange:
c       below 1.0e6 eV, use Riviere's fit
c       above 1.0e6 eV, use power-law extrapolation
c
      if (er .gt. 1.0e6)  go to 105
      siif1 = (0.6937e-14) * (1.0 - 0.06732    * aer)**2
     .                     / (1.0 + 0.1112e-14 * er**3.3)
      go to 110
  105 siif1 = (4.8363e-22) / ((1.0e-6)*er)**3.3
c
c     ion-impact ionization.
c     at low energies or high energies, use Riviere's fits
c     at intermediate energies, use fit to rkj's curve.
c
c     Freeman and Jones option
c
  110 if (er .lt.    807.4)  go to 113
      if (er .gt. 154400.0)  go to 114
      siif2 = EXP (polyf(ai1,6,aer))
      go to 120
c
  113 siif2 = EXP (-80.206+8.156*aer-0.3784*aer*aer)
      go to 120
  114 siif2 = (1.56346e-12/er)*(aer-1.792160)
c
  120 siif  = siif1 + siif2
      return
c
  200 ansq  = (nfhx(i))**2
      siif  = (1.465e-11)*(ansq/er)*(1.0-EXP (-ansq*er/9260.0))
      return
c
      end

      real*8 function sizf (i, z, er)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- cross section for electron loss via collision with impurity.
c
      common /b1/ kdene,kdeni,kdenz,ksvi,ksvz,ksve,krad,ngh,ngl
c
      ni = nfhx(i)
      id = 2
      if (ni .eq. 1)  id = 1
c
      ansq = (FLOAT (ni))**2
      sizf = (1.465e-11)*(ansq*z*z/er)*(1.0-EXP (-ansq*er/(31250.0*z)))
      return
c
      end

      subroutine slow1 (atw,atwf,ene,en,enn,ezero,ibcur,ibcx,ifirst,kj,
     .                  nj,nion,pzero,te,vz,zsq,zsqf,bke,bki,ecrit,
     .                  emzrat,encap,ge,gi,gth,taus,fionx,rtstcx)
c
      USE constnts,only : xmassp,xmasse,charg4
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  This subroutine evaluates the transfer functions N, Ge, Gi, Ke, and
c     Ki defined by Callen et al., IAEA Tokyo, Vol. I, 645 (1974) to
c     describe fast ion slowing down.  Ge and Gi are the fractions of
c     the fast ion energy transferred to the electrons and thermal ions,
c     respectively.  The remaining energy fraction is assumed lost due
c     to secondary charge exchange, i.e., charge exchange between fast
c     ions and thermal neutrals.  Similarly, Ke and Ki are the fractions
c     of the fast ion parallel momentum transferred to electrons and
c     thermal ions.  Moreover, N, Ge, and Ke are proportional to the
c     slowing-down times for particles, energy, and parallel momentum
c     and are used in subroutine slow2.  The fast ions may be beam ions
c     from neutral beam heating or alpha particles from fusion.
c  This subroutine has been modified to account for plasma rotation.
c     Expressions contained here are valid for Vrot/Vthi< = 1 (although
c     assumptions of Maxwellian target distributions may be invalid).
c     An additional factor, Gth, describes momentum and energy sources
c     to the target ions due to nonzero kinetic energy at thermalization
c     with the rotating target.  (gts: 14jul89)
c ----------------------------------------------------------------------
c
      dimension atw(*), ene(*), en(kj,*), enn(kj,*), pzero(*), te(*),
     .          vz(*), zsq(kj,*)
      dimension bke(*), bki(*), ecrit(*), emzrat(*), encap(*), ge(*),
     .          gi(*), gth(*), taus(*)
      data      pi, rot2pi / 3.14159, 2.50663 /
      data      xkeverg / 6.241e+08 /
 
c
c  calculate constants
c
      econst = 0.5 * (4.5 * pi * xmassp/xmasse) ** 0.333333  ! nomilaly 14.8
      tconst = atwf * xmassp / (zsqf * xmasse)
      tiny   = 0.0001
c
c  begin loop over mesh points
c
      do j=1,nj
        if (ifirst .eq. 0)  go to 20
c
c       calculate ecrit, the "critical energy" at which fast ions
c       transfer equal energy to electrons and thermal ions;
c       also calculate emzrat = mi*<Z>/(mf*[Z])
c
        sum1 = 0.0
        sum2 = 0.0
        do 10 k=1,nion
        sum1 = sum1 + en(j,k)*zsq(j,k)
   10   sum2 = sum2 + en(j,k)*zsq(j,k)/atw(k)
        ecrit(j)  = econst*te(j)*atwf*(sum2/ene(j))**0.666667        !kev
        emzrat(j) = sum1 / (atwf*sum2)
c
c       calculate taus, the Spitzer momentum exchange time
c       for electron-ion collisions
c
        tegt0   = MAX (te(j), tiny)  ! keV, avoid probs near plasma edge
        teerg   = 1.6e-9*tegt0
        xlam    = 24.0 - LOG (SQRT (ene(j))/(1.0e3*tegt0))
        ztaue   = SQRT (xmasse*teerg**3)
     .          / (1.33333*rot2pi*ene(j)*charg4*xlam)
        taus(j) = tconst*ztaue
c
c       calculate erel0, the fast-ion initial energy in the rotating
c       frame; this is  0.5 * mass_fast_ion * (v_fastion-v_rot)**2
c
   20   erot  = 0.0
        prot  = 0.0
        if (vz(1) .ne. 0.0) then
          erot  = (pzero(j)*vz(j) - 0.5 * atwf*xmassp*vz(j)**2)*xkeverg
          prot  = atwf*xmassp*vz(j)
        end if
        erel0 = ABS (ezero - erot)
        prel0 = pzero(j) - prot
c
c       calculate vrat = SQRT (ecrit/ezero) and taurat = taus/taucx,
c       where taucx is the charge exchange time
c
        vrat   = SQRT (ecrit(j)/erel0)
        taurat = 0.0
        if (ibcx .eq. 0  )  go to 50
        if (atwf .gt. 3.0)  go to 50
c
        taurat_part = taus(j) * (enn(j,1)+enn(j,2))
        ivbcxit     = 0
        vbeamcx     = SQRT (2.0*erel0 /(xkeverg*atwf*xmassp))
        if (rtstcx .gt. -5.0)
     .    vbeamcx = ABS (rtstcx) * vbeamcx   ! use rtstcx for testing
        erel0avg  = erel0
 450    if (rtstcx .gt. 0.0) then
          taurat  = taurat_part * cxr(erel0/atwf)
          taurat  = taurat * rtstcx ! rtstcx is input, defaulted to 1.0
        else ! use sigma(v)*v rather than Maxwellian average
c              assume neutrals have same bulk speed as thermal ions
c              i.e., use erel0 as relative energy
          taurat  = taurat_part * cxrv(erel0avg)*vbeamcx
        end if
c
c       calculate N( = encap), ge, and gi
c       assume taurat independent of v:
c
   50   encap(j) = encapf(vrat,taurat)
        ge(j)    = gef(vrat,taurat)*erel0/ezero
        if (ge(j) .lt. 0.0)  ge(j) = 0.0
        if(ge(j) .le. 0.0)ge(j) = 1.e-15 !HSJ 4/6/01
c
        if (rtstcx .lt. -5.0) then ! use average fast ion speed
c                                    valid only after beam has
c                                    established itself .. HSJ .. 6/5/96
          erel0avg  = erel0*ge(j)/(2.*encap(j))
          vfast_avg = SQRT (erel0avg/(xkeverg*atwf*xmassp))
          vcxerr    = ABS ((vfast_avg-vbeamcx)/vfast_avg)
          vbeamcx   = 0.5*(vfast_avg+vbeamcx)
          ivbcxit   = ivbcxit+1
          if (vcxerr .gt. 0.05 .and. ivbcxit .lt. 10)  go to 450
        end if
c
        gi(j)    = erel0/ezero - (1.0 + 0.5 * taurat)*ge(j)
        if (gi(j) .lt. 0.0)  gi(j) = 0.0
        gth(j)   = 1.0 - taurat*encap(j)
c
c       fionx allows testing of non-classical slowing down
c
        if (fionx .eq. 0.0)  go to 70
        ft    = fionx*gi(j)
        if (ft .gt. ge(j))  ft = ge(j)
        gi(j) = gi(j) + ft
        ge(j) = ge(j) - ft
c
c       calculate bke and bki
c
   70   prat   = 1.0
        if (pzero(j) .ne. 0.0) prat = prel0/pzero(j)
        !bke(j) = bkef(vrat,taurat,emzrat)*prat
        bke(j) = bkef(vrat,taurat,emzrat(j))*prat    ! corrected 6/7/2012 ! 8888899999
        bki(j) = prat - (1.0 + taurat)*bke(j)
c
c       end loop over mesh points
c
      end do
      return
c
      end

      subroutine slow2 (bke,dtt,encap,enfsav,ge,ibcur,nj,ppfsav,qfsav,
     .                  sfsav,spfsav,taus,wfsav,enf,ppf,qf,sf,spf,taupf,
     .                  tauppf,tauef,wf,iangrot,pprf,spbrf,spbrfsav,
     .                  pprfsav)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  This subroutine models fast ion slowing down.  The fast ions may
c     be beam ions from neutral beam heating or alpha particles from
c     fusion.  let the fast ion distribution function be given by
c         f(v,r,t) = g(v)*n(r,t).
c     This subroutine assumes that g(v) is fixed at all times at its
c     asymptotic value. Transients due to beam turn on and off are
c     accounted for only in n(r,t). That is, the fast ion density n(r,t)
c     at time t, is adjusted to produce the correct asymtotic source rate.
c     The change in g(v) over time to its asymtotic value is neglected.
c     for beam ions (not checked for fusion alphas) the quantities
c ----------------------------------------------------------------------
c
      dimension bke(*), encap(*), enfsav(*), ge(*), ppfsav(*), qfsav(*),
     .          sfsav(*), spfsav(*), taus(*), wfsav(*)
      dimension enf(*), ppf(*), qf(*), sf(*), spf(*), taupf(*),
     .          tauppf(*), tauef(*), wf(*),
     .          pprf(*),spbrf(*),spbrfsav(*),pprfsav(*)
c
c  begin loop over mesh points
c
      do 100 j=1,nj
c
c  calculate slowing down times
c
      taupf (j) = taus(j)*encap(j)
      tauppf(j) = taus(j) * ABS (bke(j))
      tauef (j) = 0.5 * taus(j)*ge(j)
c
c  calculate fast ion particle density and delayed particle source
c
      enf(j) = enfsav(j) * EXP (-dtt/taupf(j))
     .       + sfsav(j)*taupf(j)*(1.0-EXP (-dtt/taupf(j)))
      sf (j) = enf(j)/taupf(j)
c
c  calculate fast ion energy density and delayed energy source
c
      if(tauef(j) .le. 0.0)then 
         print *,"tauef .le. 0 at j =",j
         print *,"taus(j) =",taus(j)
         print *,"ge(j) =",ge(j)
      endif
      wf(j) = wfsav(j) * EXP (-dtt/tauef(j))
     .      + qfsav(j)*tauef(j)*(1.0-EXP (-dtt/tauef(j)))
      qf(j) = wf(j)/tauef(j)
c
c  calculate density and delayed source of fast ion parallel momentum
c
      if (ibcur .eq. 0)  go to 90
      ppf(j) = ppfsav(j) * EXP (-dtt/tauppf(j))
     .       + spfsav(j)*tauppf(j)*(1.0-EXP (-dtt/tauppf(j)))
      spf(j) = ppf(j)/tauppf(j)
   90 if (iangrot .eq. 0)  go to 100
c
c --- fast ion angular momentum density and delayed source
c
      pprf (j) = pprfsav(j) * EXP (-dtt/tauppf(j))
     .         + spbrfsav(j)*tauppf(j)*(1.0-EXP (-dtt/tauppf(j)))
      spbrf(j) = pprf(j)/tauppf(j)
c
  100 continue    ! end loop over mesh points
c
      return
c
      end

      subroutine smooth (x, y, n, ideg, nuse, nupdate)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension  x(*), y(*), coeff(8)
c
      do i=1,8
        coeff(i) = 0.0
      end do
c
      bigg = -1.0
c
      do i=1,n
        if (y(i) .gt. bigg)  bigg = y(i)
      end do
c
      smalll = bigg * 0.001
      i1 = 0
 2100 i1 = i1+1
      i2 = i1+nuse-1
      if (i2 .le. n)  go to 2120
c
c ----------------------------------------------------------------------
c use last set of coefficients to finish with
c ----------------------------------------------------------------------
c
      do i=i1,n
        y(i) = polval(x(i),coeff,ideg)
      end do
c
      go to 2150
c
 2120 if (y(i1) .le. smalll .and. y(i1+1) .le. smalll)  go to 2100
      call polfit(ideg,nuse,x(i1),y(i1),coeff,ier)
c
      do i=i1,i1+nupdate-1
        y(i) = polval(x(i),coeff,ideg)
      end do
c
      go to 2100
c
 2150 do i=1,n
        if (y(i) .lt. 0.0)  y(i) = 0.0
      end do
c
      return
c
      end

      subroutine solveq (ndim, n, a, b, ipvt)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c   solution of linear system, a*x = b .
c   do not use if decomp has detected singularity.
c
c   input..
c
c     ndim = declared row dimension of array containing a .
c
c     n = order of matrix.
c
c     a = triangularized matrix obtained from decomp .
c
c     b = right hand side vector.
c
c     ipvt = pivot vector obtained from decomp .
c
c   output..
c
c     b = solution vector, x .
c
      integer  ndim, n, ipvt(n), kb, km1, nm1, kp1, i, k, m
      real*8   a(ndim,n), b(n), t
c
c     forward elimination
c
      if (n .eq. 1)  go to 50
      nm1 = n-1
      do 20 k=1, nm1
         kp1 = k+1
         m = ipvt(k)
         t = b(m)
         b(m) = b(k)
         b(k) = t
         do 10 i=kp1, n
             b(i) = b(i) + a(i,k)*t
   10    continue
   20 continue
c
c     back substitution
c
      do 40 kb=1,nm1
         km1 = n-kb
         k = km1+1
         b(k) = b(k)/a(k,k)
         t = -b(k)
         do 30 i=1, km1
             b(i) = b(i) + a(i,k)*t
   30    continue
   40 continue
   50 b(1) = b(1)/a(1,1)
      return
c
      end

      subroutine sorspt (nbshape, bheigh, bwidth, bhfoc, bvfoc,
     .                   bhdiv, bvdiv, ib, ie, isourc, ke,
     .                   nsourc, sfrac1, vbeam,
     .                   x0, y0, z0, vx0, vy0, vz0)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- generates a particle at the injector surface with
c --- coordinates and velocities x0,y0,z0,vx0,vy0,vz0
c --- These coordinates are attached to the source center, with the
c --- x-direction perpendicular to the source along the source centerline
c
      external    RANDOM12                       ! random number generator
      character*8 nbshape(*)
      dimension   bheigh(*), bwidth(*), bhfoc(*), bvfoc(*),
     .            bhdiv(*), bvdiv(*), sfrac1(*), vbeam(ke,*)
      data        pio180 /0.017453293/
      data        rt2    /1.414213562/
      data        seed0  /0.0/
c
      x0     = 0.0
      isourc = 1
      if (nsourc .eq. 1)  go to 10
c
c     two sources
c     sfrac1 is fraction of source current coming from source 1
c     (upper or rightmost source, assuming positive offsets)
c
*     if (RANF   (   ) .gt. sfrac1(ib))  isourc = -1
      if (RANDOM12 (seed0) .gt. sfrac1(ib))  isourc = -1
   10 if ( nbshape(ib) .ne. 'circ'    )  go to 20
c
*  12 y0   = RANF   (     ) - 0.5
   12 y0   = RANDOM12 (seed0) - 0.5
*     z0   = RANF   (     ) - 0.5
      z0   = RANDOM12 (seed0) - 0.5
      zsqu = y0**2 + z0**2
      if (zsqu .gt. 0.25)  go to 12
      y0   = y0*bwidth(ib)
      z0   = z0*bwidth(ib)
      go to 30
c
*  20 y0 = bwidth(ib) * (RANF   (     ) - 0.5)
   20 y0 = bwidth(ib) * (RANDOM12 (seed0) - 0.5)
*     z0 = bheigh(ib) * (RANF   (     ) - 0.5)
      z0 = bheigh(ib) * (RANDOM12 (seed0) - 0.5)
c
c  Special coding for DIII-D Long Pulse Sources
c
   30 if (nbshape(ib) .ne. 'rect-lps')  go to 32
      if (   ABS (z0) .gt.  12.0     )  go to 31
c
c  particle on central 2 modules
c
      vdx = -1.0
      vdy =  0.0
      vdz =  0.0
      go to 33
c
c  particle on upper or lower modules
c
   31 vdx = -1.0
      vdy =  0.0
      vdz = -z0 / bvfoc(ib)
      go to 33
c
   32 vdx = -1.0
      vdy = -y0 / bhfoc(ib)
      vdz = -z0 / bvfoc(ib)
c
   33 vsqrt = 1.0 / SQRT (vdx**2+vdy**2+vdz**2)
      vx0   = vbeam(ie,ib)*vdx*vsqrt
      vy0   = vbeam(ie,ib)*vdy*vsqrt
      vz0   = vbeam(ie,ib)*vdz*vsqrt
c
      thz   = ranorm (1) * bvdiv(ib) / rt2
      thy   = ranorm (1) * bhdiv(ib) / rt2
      vz0   =  vz0 + thz * pio180 * vx0
      vy0   =  vy0 + thy * pio180 * vx0
      return
c
      end

      subroutine support (a, ni, nj, nk, nl, isupp, ifail)
c
      USE param
      USE mhdpar
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c determines crude 2d support of four-dimensional array
c a(i,j,k,l)   array searched for support approximation
c isupp(4,k,l) output.  with k,l fixed, contains (in order),
c   i1,i2,j1,j2, the indices of a corresponding to the
c   vertices of the smallest rectangle supporting a.
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'mhdpar.i'
c
      parameter (nxx2 = nw*2, nyx2 = nh*2)
c
      dimension a(nxx2,nyx2,ke,kb), isupp(4,ke,kb)
c
      do 10 i=1,ni
      do 10 j=1,nj
        if (a(i,j,nk,nl) .ne. 0)  go to 12
   10 continue
      go to 99
   12 continue
      isupp(1,nk,nl) = i
      do 20 i=ni,1,-1
      do 20 j=1,nj
      if (a(i,j,nk,nl) .ne. 0.0)  go to 22
   20 continue
      go to 99
   22 continue
      isupp(2,nk,nl) = i
      do 30 j=1,nj
      do 30 i=1,ni
      if (a(i,j,nk,nl) .ne. 0.0)  go to 32
   30 continue
      go to 99
   32 continue
      isupp(3,nk,nl) = j
      do 40 j=nj,1,-1
      do 40 i=1,ni
      if (a(i,j,nk,nl) .ne. 0.0)  go to 42
   40 continue
      go to 99
   42 continue
      isupp(4,nk,nl) = j
      ifail = 0
      return
c
   99 ifail = 1
      return
c
      end


      real *8 function tau0_cut_avg(taus,vcrit,kj,kb,nj,ib)
c------------------------------------------------------------------------
c     get an average value over the spatial grid, for the lower limit
c     in the fast ion integrals.
c
c-------------------------------------------------------------------------
      implicit none
      integer*4 j,nj,kj,kb,ib
      real*8    taus(kj),vcrit(kj,kb),term,t_avg,
     .          one_third,vthcut3,vc3,eval_vthcut3


      one_third =1.d0/3.d0
      t_avg =0.0d0
      do j=1,nj
         vc3 = vcrit(j,ib)**3
         vthcut3 = eval_vthcut3(j)
         term = taus(j)*one_third*LOG((vthcut3 +vc3)/vc3)
         t_avg= t_avg+term
      enddo
         tau0_cut_avg = t_avg/nj
      return
      end


      subroutine timtor (rin, rmax, x0, y0, z0, vx0, vy0, vz0,
     .                   zmin, zmax, tenter, texit)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  This subroutine calculates the times for a particle to enter and
c     exit a toroidal box surrounding the plasma starting from the
c     point (x0,y0,z0) and moving with velocity (vx0,vy0,vz0).
c ----------------------------------------------------------------------
c
      dimension  edge(4), timsol(6)
c
c  specify edges of toroidal box surrounding plasma
c
      edge(1) = rin
      edge(2) = rmax
      edge(3) = zmin
      edge(4) = zmax

c
c  find times for particle to intersect inside and outside of box
c
      isol = 0
      aa   = vx0**2+vy0**2
      if (aa .eq. 0.0)  go to 50
      bb   = 2.0 * (vx0*x0+vy0*y0)
      do 40 i=1,2
      cc   = x0**2+y0**2-edge(i)**2
      arg  = bb**2-4.0 * aa*cc
      if (arg) 40, 30, 10
   10 sqr  = SQRT (arg)
      tt   = (-bb-sqr)/(2.0*aa)
      zz   = z0+vz0*tt
      if (zz .lt. zmin)  go to 20
      if (zz .gt. zmax)  go to 20
      isol = isol+1
      timsol(isol) = tt
   20 tt   = (-bb+sqr)/(2.0*aa)
      zz   = z0+vz0*tt
      if (zz .lt. zmin)  go to 40
      if (zz .gt. zmax)  go to 40
      isol = isol+1
      timsol(isol) = tt
      go to 40
   30 tt   = -bb/(2.0*aa)
      zz   = z0+vz0*tt
      if (zz .lt. zmin)  go to 40
      if (zz .gt. zmax)  go to 40
      isol = isol+1
      timsol(isol) = tt
   40 continue
c
c  find times for particle to intersect top and bottom of box
c
   50 if (vz0 .eq. 0.0)  go to 70
      do 60 i=3,4
      tt = (edge(i)-z0)/vz0
      xx = x0+vx0*tt
      yy = y0+vy0*tt
      rr = SQRT (xx**2+yy**2)
      if (rr .lt. rin)  go to 60
      if (rr .gt. rmax)  go to 60
      isol = isol+1
      timsol(isol) = tt
   60 continue
   70 continue
c
c  return if particle misses box
c
      tenter = -1.0e10
      if (isol .eq. 0)  return
c
c  calculate times to enter and exit box
c
      tenter = timsol(1)
      iin    = 1
      do 80 i=2,isol
      if (timsol(i) .ge. tenter)  go to 80
      iin = i
      tenter = timsol(i)
   80 continue
      texit = 1.0e10
      do 90 i=1,isol
      if (i .eq. iin)  go to 90
      if (timsol(i) .lt. texit) texit = timsol(i)
   90 continue
      return
c
      end

      subroutine tozone (x, y, n, xz, yz, nz, nout, ncrt)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c convert from point values to zone values
c ----------------------------------------------------------------------
c
      dimension  x(*), y(*), xz(*), yz(*)
      data       tol/1.0e-20/
c
      ip     = 1
      iz     = 1
      yzlast = y(n)
c
c ----------------------------------------------------------------------
c interpolate to new grid that is compute yz(xz(iz))
c ----------------------------------------------------------------------
c
 2000 if (iz .gt. nz)  go to 2210
 2100 if (ip .gt. n )  go to 9200
      del    = reldif(x(ip),xz(iz))
      if (del .gt. tol)  go to 2110
      yz(iz) = y(ip)
      iz     = iz + 1
      go to 2000
 2110 if (x(ip) .gt. xz(iz))  go to 2120
      ip   = ip + 1
      go to 2100
 2120 ipm1   = ip - 1
      s      = (y(ip)-y(ipm1))/(x(ip)-x(ipm1))
      yz(iz) = y(ipm1)+(xz(iz)-x(ipm1))*s
      iz     = iz + 1
      go to 2000
c
c ----------------------------------------------------------------------
c convert to zone
c ----------------------------------------------------------------------
c
 2210 do iz=1,nz-1
        yz(iz) = 0.5 * (yz(iz)+yz(iz+1))
      end do
      yz(nz) = 0.5 * (yz(nz)+yzlast)
      return
c
c ----------------------------------------------------------------------
c fatal errors
c ----------------------------------------------------------------------
c
 9200 write (nout, 8000)
      write (ncrt, 8000)
 8000 format (' FATAL ERROR in subroutine TOZONE' //
     .   ' xz(iz) or xz(izp1) is not contained in interval x(1),x(n)' /)
      call STOP ('subroutine TOZONE: unspecified problem', 49)
c
      end
      

      real *8 function vsqfunc( taus_l,x)
      implicit none
      real *8  taus_l,x
        vsqfunc=(exp(3.D0*x/taus_l)-1.d0)**0.66666666666
      return 
      end

      subroutine wrap_xboley (atw,zzi,ebkev,ibion,mb,mfm1,nebin,ebfac,
     .                      debin,zne,zte,zni,nion,sgxn,sgxnmi,atw_beam)
c
c ----------------------------------------------------------------------
c
c  OUTPUT
c  debin      energy bin width if nebin .ne. 0
c
c introduce a wrapper that will allow inclusion of onetwo include
c files without affecting xboley subroutine
c ----------------------------------------------------- HSJ ------------
c
      USE param
      USE mhdpar
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c      include 'param.i'
c      include 'mhdpar.i'
      dimension sgxnd(kj,ke,kb)
      dimension ebkev(*),zne(*),zni(kz,*),zte(*),zzi(kz,*),
     .          sgxn(kcmp1,kz,kbe,*),sgxnmi(ke,*),atw(*),
     .          dnz(kz,4)              ! dnz is local to this subroutine
c
      if (nebin .eq. 0) then   ! no plasma rotation case
          call xboley (atw, zzi, ebkev, ibion, mb, mfm1, zne, zte, zni,
     .                 dnz, nion, sgxnd, sgxnmi, atw_beam, kz, nw, nh,
     .                 maxp, nap, kion, ke, kb)
          do ie=1,ke
             do ib=1,mb
                do j=1,mfm1
                   do i=1,kcmp1
                      if (i .eq. 4) then
                          sgxn(i,j,3*(ib-1)+ie,1)=sgxnd(j,ie,ib)
                      else
                          sgxn(i,j,3*(ib-1)+ie,1)=0.0
                      end if
                    end do
                end do
             end do
           end do
       else
c
c        include toroidal rotation; use existing method of energy bins
c
c            get the energy bin structure
c
             if (nebin .ne. 0) then
                ebmax = 0.0
                do ib=1,mb ! use largest beam energy to define the bins
                  ebmax = MAX (ebmax, ebkev(ib))
                end do
                ebmax = ebmax / atw_beam
                debin = ebmax * ebfac / FLOAT (nebin)
                jemax=nebin
              end if
c
c have to break out individual beams because of energy bin structure
c
             do ib=1,mb
                if (nebin .eq. 0) then ! use to eliminate nebin=0 branch
                  debin = ebkev(ib)
                  jemax = 1
                end if
                do je=1,jemax
                   ebin = FLOAT (je) * debin * atw_beam
                   mb_1 = 1    ! force xboley to work on 1 beam at atime
                   kbb  = 1
                   call xboley (atw, zzi, ebin, ibion, mb_1, mfm1, zne,
     .                          zte, zni, dnz, nion, sgxnd, sgxnmi,
     .                          atw_beam, kz, nw, nh, maxp, nap, kion,
     .                          ke, kbb)
                   do ie=1,ke
                     do j=1,mfm1
                       do i=1,kcmp1
                         if (i .eq. 4) then
                           sgxn(i,j,3*(ib-1)+ie,je)=sgxnd(j,ie,ib)
                         else
                           sgxn(i,j,3*(ib-1)+ie,1)=0.0
                         end if
                       end do
                     end do
                   end do
                end do
             end do
      end if
      return
c
      end

      subroutine xboley (atw, zzi, ebkev, ibion, mb, mfm1, zne, zte,
     .                   zni, dnz, nion, sgxn, sgxnmi, atw_beam,
     .                   kz, ki, kj, maxp, nap, kion, ke, kb)
c
c   This subroutine calculates the relevant cross section quantities
c   -sgxn- and -sgxnmi-, needed by the NFREYA code, using fits that are
c   based on the recent work by Janev, Boley and Post. See routine
c   sigfit for more details.
c   created  on 01-05-90 by J. Mandrekas
c   modified on 01-08-91 by J. Mandrekas to treat He3 as He4 for the
c   ARIES-III calculations
c   modified on 02-20-92 by J. Mandrekas for SUN/486/UNICOS platforms
c
****  include 'params.inc' ! original include file from Mandrekas
****                         equivalent info set up in WRAP_XBOLEY
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      dimension atw(kion), zzi(kz,kion), ebkev(kb), zne(kz), zte(kz),
     .          zni(kz,kion), sgxn(kz,ke,kb), sgxnmi(ke,kb), denz(4),
     .          dnz(kz,4)
      data      jhe3 /0/, jhe4 /0/, jcarb /0/, joxy /0/, jfe /0/,
     .          dhe3 /0.0/, dhe4 /0.0/
c
c     zero-out the impurity densities
c
      do   i=1,mfm1
        do k=1,4
          dnz(i,k) = 0.0
        end do
      end do
c
c --- recognize impurity species (including He)
c
      do i = 1, nion
        if (atw(i) .eq.  3.0 .and. zzi(1,i) .eq.  2.0)  jhe3  = i
        if (atw(i) .eq.  4.0 .and. zzi(1,i) .eq.  2.0)  jhe4  = i
        if (atw(i) .eq. 12.0 .and. zzi(1,i) .eq.  6.0)  jcarb = i
        if (atw(i) .eq. 16.0 .and. zzi(1,i) .eq.  8.0)  joxy  = i
        if (atw(i) .eq. 56.0 .and. zzi(1,i) .eq. 26.0)  jfe   = i
      end do
      do iz = 1, mfm1
        do k = 1, nion                  ! this loop is bogus HSJ
          if (jhe3  .ne. 0)  dhe3      = zni(iz,jhe3)
          if (jhe4  .ne. 0)  dhe4      = zni(iz,jhe4)
                             dnz(iz,1) = dhe3 + dhe4
          if (jcarb .ne. 0)  dnz(iz,2) = zni(iz,jcarb)
          if (joxy  .ne. 0)  dnz(iz,3) = zni(iz,joxy)
          if (jfe   .ne. 0)  dnz(iz,4) = zni(iz,jfe)
        end do
      end do
c
      do iz = 1, mfm1
        denz(1) = dnz(iz,1)
        denz(2) = dnz(iz,2)
        denz(3) = dnz(iz,3)
        denz(4) = dnz(iz,4)
        te   = zte(iz)
        dene = zne(iz)
        do   ib=1,mb
          do ie=1,3
            ebeam = ebkev(ib) / (ie*atw_beam)
            call sigfit (ebeam, dene, te, denz, sig)
            sgxn(iz,ie,ib) = dene * sig
          end do
        end do
      end do
c
      do ib = 1, mb
        if (ib .le. 1 .or. ebkev(ib) .ne. ebkev(1)) then
          do j = 1, 3
            sgxnm = sgxn(1,j,ib)
            do i = 2, mfm1
              if (sgxnm .lt. sgxn(i,j,ib))  sgxnm = sgxn(i,j,ib)
            end do
            sgxnmi(j,ib) = 1.0 / sgxnm
          end do
        else
          do j = 1, 3
            sgxnmi(j,ib) = sgxnmi(j,1)
          end do
        end if
      end do
      return
c
      end

      real*8 function yinter (dxi, nxm1, x0, x, ytab)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c this function performs a fast linear interpolation (or extrapolation)
c for y from x given ytab vs xtab, where xtab starts at x0 and has
c uniform spacing 1/dxi
c ----------------------------------------------------------------------
c
      dimension ytab(*)
c
      xx     = ABS (x-x0)*dxi + 1.0
      i      = xx
      i      = MIN0 (i, nxm1)
      wt     = xx - i
      yinter = (1.0-wt)*ytab(i) + wt*ytab(i+1)
      return
c
      end
