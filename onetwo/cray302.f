      subroutine baloo
c
      USE param
      USE io 
      USE solcon
      USE soln
      USE extra
      USE numbrs
      USE mesh
      USE sourc
      USE machin
      USE tfact
      USE geom
      USE tordlrot
      USE constnts
      USE tcoef
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     This subroutine constrains the plasma pressure gradient to be
c     no greater than that allowed for marginal stability to ideal
c     ballooning modes.  This is achieved by increasing the electron
c     and ion thermal conductivities in the appropriate regions of the
c     plasma. At the present time, only elliptical cross-section
c     plasmas are considered.
c

      dimension alpc(kj), alpha(kj), dpdr(kj), ensum(kj), p(kj), pp(kj),
     .          tep(kj), tip(kj), qqconde(kj) ,qqcondi(kj), tss(kj)
c
      data      ifirst  /1/
      data      relax   /1.0/, expon /2.0/, qm /0.99/, shearm /0.01/
      data      febal   /0.5/
c
      zero = 0.0
c
c     initialize xkebal and xkibal
c
      if (ifirst .eq. 0)  go to 15
      ifirst = 0
      do 10 j=1,nj-1
        xkebal(j) = 0.0
   10   xkibal(j) = 0.0
c
      if (xdebug(11) .ne. 0.0)  relax  = xdebug(11)
      if (xdebug(12) .ne. 0.0)  expon  = xdebug(12)
      if (wsaw       .eq. 0.0)  qm     = 0.0
      if (xdebug(13) .ne. 0.0)  qm     = xdebug(13)
      if (xdebug(14) .ne. 0.0)  shearm = xdebug(14)
      if (xdebug(15) .ne. 0.0)  febal  = xdebug(15)
      if (xdebug(10) .ne. 0.0)  write (nqik, 1010)
 1010 format (/ '   n     time     q(1)    te(1)',
     .        '     taue     beta   ar(11)   ar(26)   ar(41)',
     .           '  xkebal(11)  xkebal(26)  xkebal(41)')
c
c  calculate pressure and pressure gradient in CGS units
c
   15 do 30 j=1,nj
      ensum(j) = 0.0
      do 20 k=1,nion
   20 ensum(j) = ensum(j) + en(j,k)
   30 p(j) = 1.6e-9*(ene(j)*te(j)+ensum(j)*ti(j))
      do 40 j=1,nj-1
   40 dpdr(j) = (p(j+1)-p(j))/dr(j)
c
c  calculate alpha, alphac, and alphar
c
      alpha(1) = 0.0
      do 60 j=1,nj-1
        alpc(j)   = alphac(shearp(j))
        alpha(j)  = 0.0
        alphar(j) = 1.0
        if (alpc(j) .eq. 0.0)  go to 60
        qa = 0.5*(q(j)+q(j+1))
        alpha(j) = -8.0 * pi * SQRT (kappa)*rmajor*qa**2*dpdr(j)/btor**2
        alpha(j) = MAX (alpha(j), zero)
        if (       qa .le. qm    )  alphar(j) = 0.0
        if (shearp(j) .le. shearm)  alpha (j) = 0.0
        alphar(j) = alpha(j)/alpc(j)
   60 continue
c
c  return if ibaloo < 0
c
      if (ibaloo .lt. 0)  return
c
c  test whether stability constraint is approached or exceeded,
c  i.e., whether alpha > 0.99*alphac anywhere; if not, then return
c
c     tss indicates threshold for continuous stability to ballooning modes
c     ;Te for transition = 350eV / ddebug(47)
c
      do j=1,nj-1
        tss(j) = te(j)/0.35
        tss(j) = ddebug(47)*tss(j)
        if (alpha(j) .gt. 0.99*alpc(j) .and. tss(j) .lt. 1.0)  go to 120
      end do
c
      do 115 j=1,nj-1
      xkebal(j) = 0.0
  115 xkibal(j) = 0.0
      return
c
c  calculate new pressure profile consistent with stability constraint
c
  120 pp(nj) = p(nj)
      do 140 j=nj-1,1,-1
        if (alpha(j) .gt. alpc(j))  go to 130
        pp(j) = pp(j+1) + (p(j)-p(j+1))
        go to 140
  130   dpl    = -dr(j)*dpdr(j)*alpc(j)/alpha(j)
        pp(j) = pp(j+1) + dpl
  140 continue
c
c  adjust electron and ion temperatures (in the ratio febal:1-febal) to
c  obtain reduced pressure
c
      do j=1,nj
        pep    = 1.6e-9*ene(j)*te(j) - febal*(p(j)-pp(j))
        tep(j) = pep/(1.6e-9*ene(j))
        tip(j) = (pp(j)-pep)/(1.6e-9*ensum(j))
      end do
c
c  calculate electron and ion conduction terms from energy balance
c     equations; neglect convection
c
      if (febal .eq. 0.0)  go to 222
      k = nion+1
      do 220 j=1,nj
      dpedt_lcl = 0.0
      if (alpha(j) .ge. alpc(j) .and. tss(j) .lt. 1.0)  go to 215
      if(dtt .ne. 0.0)then
         dpedt_lcl   = 1.5*(ene(j)*u(k,j)-enesav(j)*usave(k,j))/dtt
      else
         dpedt_lcl =0.0
      endif
  215 qqconve = 0.0
  220 qqconde(j) = -qdelt(j)-qexch(j)+qohm(j)-qione(j)-qrad(j)+qbeame(j)
     .             +qrfe(j)+qfuse(j)+qe2d(j)-dpedt_lcl-qqconve
     .             -iangrot*kevperg*angrcple*sprbeame(j)*angrot(j)
  222 if (febal .eq. 1.0)  go to 242
      k = nion+2
      do 240 j=1,nj
      dpidt_lcl = 0.0
      if (alpha(j) .ge. alpc(j) .and. tss(j) .lt. 1.0)  go to 235
      do  i=1,nion
         IF(dtt .ne. 0.0)then 
             dpidt_lcl   = dpidt_lcl + 1.5*(u(i,j)*u(k,j)
     .                           -usave(i,j)*usave(k,j))/dtt
         ELSE
            dpidt_lcl =0.0
         ENDIF
      enddo
  235 qqconvi = 0.0
  240 qqcondi(j) = qdelt(j)+qexch(j)+qioni(j)-qcx(j)+qbeami(j)
     .           +qrfi(j)+qfusi(j)+qi2d(j)-dpidt_lcl-qqconvi
     .  +iangrot*kevperg*angrcple*sprbeame(j)*angrot(j)
c
c  calculate ballooning mode portion of electron and ion thermal
c     conductivities; if necessary, use underrelaxation to avoid
c     numerical problems
c
  242 if (febal .eq. 0.0)  go to 252
      k = nion+1
      call flxcal (qqconde,drr,hcap,nj,r,ra,p)
c
      do 250 j=1,nj-1
      if (alpha(j) .ge. alpc(j) .and. tss(j) .lt. 1.0)  go to 245
      xkebal(j) = alphar(j)**expon*xkebal(j)
      go to 250
  245 xkeb = 0.0
      grad = ((1.0-theta)*(te(j+1)-te(j))
     .        +theta*(tep(j+1)-tep(j)))/dr(j)
      if (grad .eq. 0.0)  go to 250
      xke       = -p(j)/grad
      xkeb      = xke - (d(k,k,j)-xkebal(j))
      xkeb      = MAX (xkeb, zero)
      xkebal(j) = (1.0-relax)*xkebal(j) + relax*xkeb
  250 continue
c
  252 if (febal .eq. 1.0)  go to 262
      k = nion+2
      call flxcal (qqcondi,drr,hcap,nj,r,ra,pp)
      do 260 j=1,nj-1
      if (alpha(j) .ge. alpc(j) .and. tss(j) .lt. 1.0)  go to 255
      xkibal(j) = alphar(j)**expon*xkibal(j)
      go to 260
  255 xkib = 0.0
      grad = ((1.0-theta)*(ti(j+1)-ti(j))
     .        +theta*(tip(j+1)-tip(j)))/dr(j)
      if (grad .eq. 0.0)  go to 260
      xki       = -pp(j)/grad
      xkib      = xki-(d(k,k,j)-xkibal(j))
      xkib      = MAX (xkib, zero)
      xkibal(j) = (1.0-relax)*xkibal(j) + relax*xkib
  260 continue
c
  262 if (xdebug(10) .ne. 0.0)
     .  write (nqik, 1020)  n, time, q(1), te(1), taue, beta,
     .                      alphar(11), alphar(26), alphar(41),
     .                      xkebal(11), xkebal(26), xkebal(41)
 1020 format (i4, 8f9.4, 1p3e12.3)
c
      return
c
      end

      subroutine cdeboo (coefa, coefb)
c
      USE param
      USE numbrs
      USE mesh
      USE tordlrot
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c  subroutine calculates the flow shear suppression vectors
c  rgce_mult and rgci_mult based on the phenomenological
c  model of Jim DeBoo.
c
c   INPUT (through include files)
c   dr(i)        i=1,2..nj-1 is the rho grid spacing
c   angrot(i)    i=1,2...nj  is toroidal rotation(rad/sec)
c   rgc(i)       i=1,2,..nj  is the time averaged gradient
c                in normalized rho on the full grid.
c
c   input and output
c   coefa(i) i=1,2...nj-1 electron multiplier on half grid
c   coefb(i) i=1,2...nj-1 ion           "      "  "    "
c ------------------------------------------------------------------ HSJ
c
c      include 'param.i'
c      include 'numbrs.i'                     ! pick up nj here
c      include 'mesh.i'                       ! pick up r, dr
c      include 'tordlrot.i'                   ! angrot, all rgc parms
      include 'storage.i'                    ! xdum
c
      dimension    coefa(*), coefb(*), grad_angrot(kj)
      equivalence (xdum(1), grad_angrot(1))
c
c     get the current gradient of angrot on the half grid
c
      call difydxhalfgrid (dr, angrot, grad_angrot, nj)
      do j=1,nj-1
        dlft   = rgc(j  ) / r(nj)
        drgt   = rgc(j+1) / r(nj)
        rgcabs = ABS (0.5*(dlft+drgt))  ! abs grad on half grid
        factr  = ABS (grad_angrot(j))   ! current abs grad on half grid
        factr  = MAX (factr, rgcabs)
        if (rgcabs .ne. 0.0) then
          factr = factr / rgcabs
        else
          factr = 1.0
        end if
        rgc_mult(j) = rgca / (rgcb+rgcc*(factr**rgcd))  ! multiplier
        coefa(j)    = coefa(j)*rgce_mult*rgc_mult(j)    ! for electrons
        coefb(j)    = coefb(j)*rgci_mult*rgc_mult(j)    ! for ions
      end do
      return
c
      end

      subroutine csperp
c
c     force cs(j) = 1.0 aka G. Staebler HSJ 12/15/98
c

c
c     -- S.J. Thompson -- General Atomics -- Core Physics Fusion Group
c
c     SUBPROGRAM DESCRIPTION:
c     ----------------------
c
c     Here we calculate the Staebler-Hinton H-mode correction model
c     perpendicular shear term
c
c     bpol   d Epsi   Ls
c     Sperp = rmajor * ---- * ------ * --
c     btor     dr     cs
c
c     where
c
c     te + ti
c     cs = SQRT [ ------- ]
c     mi
c
c     and
c
c     c      d (pi)                   btor * rmajor
c     Epsi =  ----- * ------ + angrot - Kpol * -------------
c     e  ni    dr                         r2capi
c
c     If input value of lsfctr (stored in fs(15)) = 0 then we take
c
c     r       dq
c     Ls = 1 / [ -------- * -- ]
c     R * q**2   dr
c
c     else if (lsfctr .gt. 0) we take the constant value
c
c     Ls = rmajor .
c
c     In these expressions - as well as those which follow - Bp is the poloidal
c     magnetic field (Gauss), Bm is the magnetic field amplitude (Gauss), w is
c     toroidal rotation frequency (rad/sec), r is a radial flux coordinate (cm),
c     R is the major radius at the magnetic axis (cm), c is the speed of light
c     (cm/sec), e is electron charge (statcoulomb), ni is ion density (1/cm**3),
c     pi is ion pressure (keV/cm**3), Ti and Te are the main ion and electron
c     temperatures (keV), respectively, and mi is 2 * proton rest mass (g). In
c     the following, we shall use the electron density, ne (1/cm**3), in place
c     of ni.
c
c     Notice that calculated quantites are in CGS units. Where appropriate,
c     quantities will be converted.
c
c
c     *Modifications:
c
c     12-Aug-1994  S.J. Thompson. Array FS was extended to include multipliers
c     of various terms in the Sperp term of the Staebler-Hinton
c     model. By setting these values appropriately in the input
c     namelist, terms may be turned on/off. These terms include:
c
c     fs(16) =  xrot       * Multiplier of toroidal frequency derivative
c     component in Sperp term. Set to 1.0 by default.
c     fs(17) =  xeden      * Multiplier of electron density   derivative
c     component in Sperp term. Set to 1.0 by default.
c     fs(18) =  xsecder    * Multiplier of second             derivative
c     component in Sperp term. Set to 1.0 by default.
c     3-Dec-1996   G. M. Staebler. Changed definition of Sperp to toroidal
c     version (T. S. Hahm and K. H. Burrell, Phys. Plasmas vol. 2
c     1995, pg. 1645)
c
c
c     *Reference:  G.M. Staebler and F.L. Hinton, Particle and Energy
c     Confinement Bifurcation in Tokamaks, Phys. Fluids B 5
c     (4), p. 1281-1288, April 1993
c
c     NOMENCLATURE:
c     ------------
c
c     Quantity:               Description:                      Units:
c     --------                -----------                       -----
c
c     r (1, ..., nj)          Radial flux coordinate            cm
c
c     dr (1, ..., nj-1)       Spatial increment                 cm
c
c     te (1, ..., nj)         Electron temperature              keV
c
c     ti (1, ..., nj)         Main ion temperature              keV
c
c     q (1, ..., nj)          Safety factor                     dimensionless
c
c     rmajor                  Major radius at magnetic axis     cm
c
c     angrot (1, ..., nj)     Toroidal rotation frequency       rad/sec
c
c     bpol (1, ..., nj)       Poloidal magnetic field           Gauss
c
c     btor                    Magnetic field magnitude          Gauss
c
c     ni (1, ..., nj)         thermal main ion density          1 / cm**3
c
c     sperp (1, ..., nj)      Perpendicular shear               dimensionless
c
c     nprim                   Number of primary ion species     dimensionless
c
c     atw (3)                 Mass number of primary ions      proton mass units
c
c     r2capi (1, ..., nj)     flux surface average <R**2>       M**2
c
c     Kpol_d (1, ..., nj)     main ion poloidal velocity/Bpol   M/sec/T
c
c     Quantity:               Grid:                             Where stored:
c     --------                ----                              ------------
c
c     r (1, ..., nj)          Full                              mesh.i
c
c     dr (1, ..., nj-1)       Full                              mesh.i
c
c     te (1, ..., nj)         Full                              soln.i
c
c     ti (1, ..., nj)         Full                              soln.i
c
c     q (1, ..., nj)          Full                              extra.i
c
c     rmajor                  Not applicable                    machin.i
c
c     angrot (1, ..., nj)     Full                              tordlrot.i
c
c     bpol (1, ..., nj)       Full                              extra.i
c
c     btor                    Not applicable                    machin.i
c
c     ene (1, ..., nj)        Full                              soln.i
c
c     sperp (1, ..., nj)      Half                              staebler.i
c
c     nprim                   Not applicable                    numbrs.i
c
c     atw (3)                 Not applicable                    ions.i
c
c     Kpol_d (1, ..., nj)     full                              cer.i
c
c     r2capi (1, ..., nj)     full                              geom.i
c
c     --- S.J. Thompson ---------------------------------- 30 Mar 94 ---
c
      USE param
      USE ions
      USE soln
      USE mhdpar 
      USE extra
      USE numbrs
      USE mesh
      USE machin
      USE tfact
      USE geom
      USE tordlrot
      USE constnts, only : e => charge,c => cee
      USE staebler
      USE cer
      implicit  integer (i-n), real*8 (a-h, o-z)
      include 'storage.i'
c
c     ----------------------------------------------------------------------
c     Declare local variables
c     ----------------------------------------------------------------------
c
      real*8 cs(kjm1), dndr(kjm1), dpdr(kjm1), dqdr(kjm1), dwdr(kjm1),
     .     ls(kjm1), pi(kj), sd(kj), xmi, dphidpsia(kjm1),dphidpsi(kj)
      real*8  pm, cke, cTMgcm, dpsidr(kjm1), rbp1, rbp2, rbpa,
     .     deni,dtedr(kjm1)
c
c     ----------------------------------------------------------------------
c     Utilize existing storage
c     ----------------------------------------------------------------------
c
      equivalence (cs  (1), sdum(1))
      equivalence (dndr(1), tdum(1))
      equivalence (dpdr(1), udum(1))
      equivalence (dqdr(1), vdum(1))
      equivalence (dwdr(1), wdum(1))
      equivalence (ls  (1), xdum(1))
      equivalence (pi  (1), ydum(1))
      equivalence (sd  (1), zdum(1))
c
c     ----------------------------------------------------------------------
c     Define useful constants
c     ----------------------------------------------------------------------
c
c      data c   / 2.99790e+10 /  ! velocity of light (cm/s)
      data pm  / 1.67260e-24 /  ! proton rest mass (g)
c      data e   / 4.80320e-10 /  ! electron charge (statcoulomb)
      data cke / 1.60219e-09 /  ! conversion from keV to erg
      data cTMgcm / 1.0e-2 /    ! conversion from M/Tesla to cm/Gauss
c
c     ------------------------------------------------------------------
c     Calculate mass of primary ions
c     ------------------------------------------------------------------
c
      xmi = 0.0
      do j=1,nprim
         xmi = atw(j) * pm + xmi ! (g)
      end do
c
c     ------------------------------------------------------------------
c     update the poloidal rotation if cer_ion is set
c     ------------------------------------------------------------------
c
      if (cer_ion .ne. ' ')  call nclboot (jhirsh, sdum)
c
c     ------------------------------------------------------------------
c     Calculate pressure, pi = ni * ti, on the full-grid
c     ------------------------------------------------------------------
c
      do j=1,nj
         pi(j) = en(j,1) * ti(j) * cke ! (erg/cm**3)
      end do
c
c     ------------------------------------------------------------------
c     Determine the derivative of the toroidal rotation frequency, the
c     pressure, and the safety factor, respectively, on the half-grid.
c     ------------------------------------------------------------------
c
      call difydxhalfgrid (dr, angrot, dwdr, nj) ! (rad/s-cm)
      call difydxhalfgrid (dr, pi, dpdr, nj) ! (erg/cm**4)
      call difydxhalfgrid (dr, q, dqdr, nj) ! (  1/cm)
      call difydxhalfgrid (dr, te, dtedr, nj) ! (kev/cm)
c
c     ------------------------------------------------------------------
c     Calculate length vector, ls, and velocity vector, cs, on half grid
c     ------------------------------------------------------------------
c
      do j=1, nj-1
         cs(j)=1.0
c
c     cs(j) = ((te(j) + te(j+1)) + (ti(j) + ti(j+1))) / 2.0
c     cs(j) = (SQRT (cs(j) * cke / xmi)/rmajor)*
c    .         SQRT (rmajor * ABS (dtedr(j))*2.0/(te(j) + te(j+1)))
c     cs(j) = cs(j)*(q(j) + q(j+1))/2.0
c     if (fs(15) .eq. 0.0) then
c       ls(j) = rmajor * (q(j) + q(j+1))**2
c       ls(j) = ls(j) / (2.0 * (r(j) + r(j+1)) * dqdr(j))    ! (cm)
c     else
c       ls(j) = rmajor
c     end if
c
c     ------------------------------------------------------------------
c     Calculate Epsi = dphidpsi on the half grid
c     ------------------------------------------------------------------
c
         term1(j)  = 0.5 * (angrot(j) + angrot(j+1))
         rbp1      = rbp(j)/(fcap(j)*gcap(j)*hcap(j))
         rbp2      = rbp(j+1)/(fcap(j+1)*gcap(j+1)*hcap(j+1))
         rbpa      = 0.5 * (rbp1 + rbp2)
         dpsidr(j) = (rbpa / ra(j)) * rmajor
         deni      = 0.5 * (en(j,1) + en(j+1,1))
c
****     write (6, *) 'j = ',j
****     write (6, *) 'term1 = ',term1(j)
c
         term2(j)  = c * dpdr(j)/(e * deni * dpsidr(j))
c
****     write (6, *) 'term2 = ',term2(j)
c
         term3(j) = -btor*rmajor*cTMgcm* (Kpol_d(j) + Kpol_d(j+1))/
     .        (r2capi(j)*fcap(j) + r2capi(j+1)*fcap(j+1))
c
****     write (6, *) 'term3 = ',term3(j),btor,rmajor,r2capi(j)
c
         dphidpsia(j) = fs(16) * term1(j) + fs(17) * term2(j) +
     .                  fs(18) * term3(j)
c
         term1(j) = cs(j)
      end do
c
c     interpolate dphidpsia onto full mesh
c
      call intrp (1, 1, ra, dphidpsia, nj-1, r, dphidpsi, nj)
c
c     smooth dphidpsia
c
      do i=1,nj-1
         dphidpsia(i) = (dphidpsi(i+1) + dphidpsi(i))/6.0 +
     .               2.0*dphidpsia(i)/3.0
      end do
c
c     interpolate smoothed dphidpsia onto full mesh
c
      call intrp (1, 1, ra, dphidpsia, nj-1, r, dphidpsi, nj)
c
c ----------------------------------------------------------------------
c     calculate Sperp on half grid
c ----------------------------------------------------------------------
c
      call difydxhalfgrid (dr, dphidpsi, sperp, nj)
      do i = 1, nj-1
         sperp(i) = ABS (dpsidr(i)*sperp(i)/btor)/
     .        (cs(i) + 1.0e-6)  ! 1/sec
      end do
      return
c
      end

      subroutine cstaebler
c
c --- S.J. Thompson --- General Atomics --- Core Physics Fusion Group --
c
c                           SUBPROGRAM DESCRIPTION:
c                           ----------------------
c
c  Here we calculate coefficients COEFA, COEFB, COEFC and COEFD as well
c  as the perpendicular shear term SPERP appropriate for the Staebler-
c  Hinton* model. These quantities - which will serve as multiplicative
c  factors for the anomalous electron and ion diffusivities, diffusion
c  coefficient, and toroidal momentum diffusivity, respectively - are
c  described in the following way:
c
c        COEFA = aeh +              ael
c                       ---------------------------
c                       1 + (alfae * sperp)**gammah
c
c        COEFB = aih +              ail
c                       ---------------------------
c                       1 + (alfai * sperp)**gammah
c
c        COEFC =  bh +              bl
c                       ---------------------------
c                       1 + (betah * sperp)**gammah
c
c        COEFD =  ch +              cl
c                       --------------------------- .
c                       1 + (sigma * sperp)**gammah
c
c  The various quantities either defined by the user or set by default
c  are stored in array FS in the following fashion:
c
c        FS(1)  = IFSFLAG
c        FS(2)  = AEH
c        FS(3)  = AEL
c        FS(4)  = AIH
c        FS(5)  = AIL
c        FS(6)  = BH
c        FS(7)  = BL
c        FS(8)  = CH
c        FS(9)  = CL
c        FS(10) = ALFAE
c        FS(11) = ALFAI
c        FS(12) = BETAH
c        FS(13) = SIGMA
c        FS(14) = GAMMAH
c
c  Should the flow shear suppression flag IFSFLAG be set to zero these
c  calculations will not take place. Otherwise, the quantities listed
c  above are calculated for later use, and printed to a user-designated
c  output file.
c
c                               NOMENCLATURE:
c                               ------------
c
c   Variable:               Description:
c   --------                -----------
c
c   r (1, ..., nj)          Rho grid position.
c
c   nj                      Number of grid points.
c
c   nout                    I/O unit number where data is to be printed.
c
c   tstaebler               Time at which coefficients were computed.
c
c   coefa (1, ..., nj-1)    Staebler-Hinton multiplicative factor for
c                           anomalous electron diffusivity term.
c
c   coefb (1, ..., nj-1)    Staebler-Hinton multiplicative factor for
c                           anomalous ion diffusivity term.
c
c   coefc (1, ..., nj-1)    Staebler-Hinton multiplicative factor for
c                           diffusion coefficient.
c
c   coefd (1, ..., nj-1)    Staebler-Hinton multiplicative factor for
c                           toroidal momentum diffusivity.
c
c   sperp (1, ..., nj-1)    Staebler-Hinton perpendicular shear term.
c
c  *Reference:  G.M. Staebler and F.L. Hinton, Particle and Energy
c               Confinement Bifurcation in Tokamaks, Phys. Fluids B 5
c               (4), p. 1281-1288, April 1993
c
c --- S.J. Thompson ---------------------------------------------- 30 Mar 94 ---
c

      USE param
      USE numbrs
c
      USE staebler
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'numbrs.i'
c      include 'staebler.i'
c
      call csperp
      do j=1,nj-1
        coefa(j) = 1.0 - ABS (fs(10)*sperp(j))
        if (coefa(j) .lt. 0.0)  coefa(j) = 0.0
        coefa(j) = fs(3)*(coefa(j))**fs(14) + fs(2)
        coefb(j) = 1.0 - ABS (fs(11)*sperp(j))
        if (coefb(j) .lt. 0.0)  coefb(j) = 0.0
        coefb(j) = fs(5)*(coefb(j))**fs(14) + fs(4)
        coefc(j) = 1.0 - ABS (fs(12)*sperp(j))
        if (coefc(j) .lt. 0.0)  coefc(j) = 0.0
        coefc(j) = fs(7)*(coefc(j))**fs(14) + fs(6)
        coefd(j) = 1.0 - ABS (fs(13)*sperp(j))
        if (coefd(j) .lt. 0.0)  coefd(j) = 0.0
        coefd(j) = fs(9)*(coefd(j))**fs(14) + fs(8)
      end do
      return
c
      end

      subroutine evtspl (x, y, nx, cs, ic, r, tspl, npts, t)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- evaluate the tension spline at npts points of r, return results in tspl
c ----------------------------------------------------------------------
c
      dimension  x(*), y(*), cs(ic,3), r(*), tspl(*)
c
      do j=1,npts
        z = r(j)
        do i=1,nx-1
          if (z .le. x(i+1))  go to 30
        end do
        if (ABS ((z-x(nx))/x(nx)) .lt. 1.0e-8)  go to 30
c
        print *,'(z-x(nx))/x(nx) =',(z-x(nx))/x(nx)
        print *,'x =',x(1:nx)
        print *,'j,=',j
        print *,'r =',r(1:npts) 
        call STOP ('subroutine EVTSPL: unspecified problem', 124)
c
c       evaluate the spline
c
   30   dxl = z-x(i)
        dxr = x(i+1)-z
        if (t .ne. 0.0) then
             tspl(j) = cs(i,1)*cs(i,2) * SINH (t*dxr)+
     .                 (y(i)*t*t-cs(i,1))*cs(i,3)*dxr
     .                +cs(i+1,1)*cs(i,2) * SINH (t*dxl)+
     .                 (y(i+1)*t*t-cs(i+1,1))*cs(i,3)*dxl
        else
             delx    = x(i+1)-x(i)
             tspl(j) = (dxr*(-cs(i  ,1) * dxl*(dxr+delx)+6.0 * y(i  ))
     .                + dxl*(-cs(i+1,1) * dxr*(dxl+delx)+6.0 * y(i+1)))
     .                 /(6.0 * delx)
        end if
      end do
      return
c
      end

      subroutine find1 (m, n, val, array, ia)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
      real *8 is
c
      dimension array (ia)
c
c  this routine finds m an n such that array(m) < val < array(n)
c  if val = array(i), m=n=i are returned
c  if val is not within the array limits, m = n=0 are returned
c  the routine works whether array is ordered largest to least
c  or least to largest
c
      is = array(ia)-array(1)
      if (is .eq. 0)  go to 12
      is = is / ABS (is)
      if (is*(val-array(1 )))  12, 15, 10
   10 if (is*(val-array(ia)))  30, 20, 12
   12 m = 0
      n = 0
      return
   15 m = 1
      n = 1
      return
   20 m = ia
      n = ia
      return
   30 if (is .lt. 0)  go to 60
      m = 1
      n = ia
   35 im = (m+n)/2
      if (im .eq. m)  return
      if (val-array(im))  40, 50, 45
   40 n = im
      go to 35
   45 m = im
      go to 35
   50 m = im
      n = im
      return
   60 m = ia
      n = 1
   65 in = (m+n)/2
      if (in .eq. n)  return
      if (val-array(in))  70, 80, 75
   70 n = in
      go to 65
   75 m = in
      go to 65
   80 m = in
      n = in
      return
c
      end

      subroutine get_orbloss (losseval, izone, rpos, zpos, vtroid,
     .                        vrad, vz0, vbeam, phi, elossb)
c


c
c ----------------------------------------------------------------------
c     subroutine calculates the table that stores the loss boundary
c     This table is used in subroutine FREYA to determine if an ion is lost
c     if the Stambaugh model for loss boundaries is selected (i.e., if
c     iborb .ge. 2). Fz, the charge number of the ion, is assumed equal to
c     unity, for consistency with FREYA.  Subroutine also does the table
c     interpolation and returns the critical energy at a given pitch angle
c     (which of the two options is performed is set by the switches
c     losseval and losscalc; see below).
c                   (REF. "Calculating Orbit Loss in a Separatrix
c                    Bounded Tokamak",  DIII-D Physics Memo No. 9303,
c                    16 March 1993).
c --- input:
c     losseval     used if the evaluation section of get_orbloss is
c                  called (the evaluation section is at the end of
c                  this subroutine, i.e., if losscalc = 0).
c                  Set losseval to 1 if the bicubic
c                  spline of psi has to be generated. If the bicubic
c                  spline allready exists then set losseval = 0 and the
c                  generation calculation will be skipped. losseval
c                  is returned as 0
c
c --- the following input from the argument list is required if
c --- subroutine is called with losscalc = 0. These values are defined
c --- in subroutine FREYA. If the routine is called with losscalc = 1
c --- then this input is never accesed.
c
c     izone
c     rpos
c     zpos
c     vtroid
c     vrad
c     vz0
c     vbeam
c
c -------------------- INPUT THROUGH INCLUDE FILES ---------------------
c
c --- INCLUDE file param.i
c     kf              max value of mf
c
c --- INCLUDE file numbrs.i
c     nj           size of bp
c
c --- INCLUDE file etc.i
c     bp(i)         i = 1,2..nj. flux surface average poloidal b field as
c                   determined in last MHD equilibrium calculation,
c                   as a function of psir(bp in kgauss)
c --- INCLUDE file ions.i
c     atw(i)    i = 1,2,..nion
c
c --- INCLUDE file nub.i
c     atw_beam  mass no. of beam
c     nbeams    number of beams
c     ebkev(i)  i = 1,2...nbeams,beams energy,keV
c     mf        number of radial values in vectors such as
c               rotsid,etc.
c
c --- INCLUDE file nub2.i
c     these vectors are set in prenub
c     rotsid(i) i = 1,2..mf major radius of flux zones(on outboard side)
c     fpsio(i)  corresponding values of f( psi)
c     potsid(i) corresponding values of psi
c
c --- INCLUDE file machin.i
c     btor            in GAUSS
c
c --- INCLUDE file mhdpar.i
c     nw              size of mhdgrid
c     nh
c
c --- INCLUDE file limiter.i
c
c     xlimiter(i)
c     ylimiter(i) i = 1,2,...nlimiter   ! THESE VALUES ARE IN METERS !
c
c --- INCLUDE file mhdgrid.i
c     rmhdgrid(i)     i = 1,2..nw, CM
c     zmhdgrid(j)     j = 1,2..nh, CM   the rectangular MHD grid
c
c --- INCLUDE file small.i
c     p(i,j)          i = 1,2..nw,j=1,2..nh,is the psi array defined
c                     over the rmhdgrid,zmhdgrid rectangle. here p has
c                     units of GAUSS*CM^2
c     xax(1)          major radius to magnetic axis
c
c --- INCLUDE file bicube.i
c     n2cspln         quantities needed to get bicubic spline representation
c     cspln           of psi
c     wnoperm
c     nh2
c     pds(6)
c
c --- INCLUDE file contour.i
c     rplasbdry(i)    i = 1,2..nplasbdry
c     zplasbdry(i)    r,z points describing plasma boundary
c     nplasbdry       # of points in above vectors
c
c --- INCLUDE file oloss.i
c     losscalc        =0 don't do the setup calculations.
c                     =1 do the setup calculations and return.
c                     losscalc is set to 1  whenever
c                     tport is called (which happens when a new
c                     equilibrium calculation has ben done)
c                     losscalc is set to 0 below in this subroutine
c                     so that subsequent calls this subroutine do not
c                     perform the setup calculations. Instead the
c                     evaluation section of the subroutine is executed.
c
c --- INCLUDE file constnts.i
c     pi
c
c --- INCLUDE file rhog.i
c     psir(i)                i = 1,2,..nj psi on rho grid,kgauss-cm**2
c
c --- INCLUDE file soln.i
c     rbp(i)                 i = 1,2..nj, rbp=fcap*gcap*hcap*r*bp
c
c --- INCLUDE file geom.i
c     fcap(i)
c     gcap(i)
c     hcap(i)           i = 1,2..nj geometric factors
c
c --- INCLUDE file mesh.i
c     r(i)               i = 1,2..nj the rho grid (in cm)
c
c ----------------------------------- OUTPUT ---------------------------
c
c --- OUTPUT THROUGH ARGUMENT LIST:
c     losseval        explained above
c     phi and elossb are output if this routine is called with losscalc = 0
c     phi             pitch angle in degrees
c     elossb          critical energy in keV
c
c --- OUTPUT TO INCLUDE FILE OLOSS.I
c     npitch          number of pitch angles (actually set in oloss.i)
c     pitchrmaj(j)    j = 1,2...mf  (currently a copy of rotsid)
c     pitcheng(j,i)   j = 1,..npitch,i=1,..mf. the energy, keV, above which
c                     and ion born at position pitchrmaj(i), with pitch angle
c                     pitchang(j) is lost.
c     bpeq(i)         flux surface average bpoloidal from last equilibrium
c     bpoltp(i)       current flux surface average value
c                     Note: the difference between bpeq and bpoltp is that
c                     bpoltp is a time evolved version of bpeq. (profided
c                     that the current evolution is turned on, i.e., itxj = 1)
c                     At the beginning of each transport/mhd equilibrium cycle
c                     bpeq and bpoltp will be the same.
c     losscalc        returned set to zero whenever this routine
c                     is called
c
c ------------------------------------------------------------------ HSJ
c
      USE param
      USE ions
      USE nub
      USE nub2
      USE soln
      USE contour
      USE limiter
      USE mhdpar
      USE mhdgrid
      USE numbrs 
      USE mesh
      USE machin
      USE geom
      USE constnts
      USE rhog
      USE bicube
      USE etc
      USE oloss_dat
      USE gpsi
      USE replace_imsl,            ONLY : my_ibcccu

      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'ions.i'
c      include 'nub.i'
c      include 'nub2.i'
c      include 'numbrs.i'
c      include 'contour.i'
c      include 'limiter.i'
c      include 'machin.i'
c      include 'mhdpar.i'
c      include 'mhdgrid.i'
c      include 'small.i'
c      include 'bicube.i'
c      include 'constnts.i'
c      include 'rhog.i'
c      include 'etc.i'
c      include 'soln.i'              ! get rbp
c      include 'geom.i'              ! get fcap, gcap, hcap
c      include 'mesh.i'              ! get r (i.e., rho)
c
      if (losscalc .eq. 1) then     ! only do the setup calculations
c
c       set up the inputs required by subroutine olossa
c
        cconst = 0.0
        call multpl1 (pitchlos, npitch*kf, cconst) ! zero pitchlos array
c
        fm = atw_beam            ! mass of ion in proton masses
        fz = 1.0                   ! charge of ion; FREYA only accepts 1
c
        if (nbeams .gt. 1) then    ! largest energy to be used in olossa
           tmax = 0.0
           do jb=1,nbeams
             tmax = MAX (tmax, ebkev(jb))
           end do
        else
           tmax = ebkev(1)
        end if
        tmax = 1050*tmax ! tmax is in eV, slightly larger than beam eng.
c
c --- set up bicubic spline representation of psi
c
        ier = 0
        call my_ibcccu (p,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .               wnoperm,ier)
        if (ier .ne. 0)
     .    call STOP ('subroutine GET_ORBLOSS: IBCCCU problem #1', 167)
c
c --- get location of x point. If there is no x point then
c --- assume that the point with the smallest value of bpoloidal
c --- can be used. Get rsepin and rsepout,the min and max major
c --- radius extension of the plasma.
c
        bsqmin = 5.0 * ABS (btor)           ! btor is in gauss
        rsepout = rmhdgrid(1)
        rsepin = rmhdgrid(nw)
        do j=1,nplasbdry
           call dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,
     .                  nw,rplasbdry(j),zplasbdry(j),
     .                  pds,ier,3)
           bsq = (pds(2)**2+pds(3)**2)/rplasbdry(j)**2
           bsqmin = MIN (bsqmin,bsq)
           if (bsqmin .eq. bsq) then
             rx   = rplasbdry(j)
             zx   = zplasbdry(j)
             psix = pds(1)       ! psix should (must) equal potsid(mf)
           end if
           rsepout = MAX (rsepout, rplasbdry(j))
           if (rsepout .eq. rplasbdry(j))  zsep0 = zplasbdry(j)
           rsepin = MIN (rsepin,rplasbdry(j))
        end do
c
c --- set up radial grid from just inside rx to rsepout at elevation of x point
c
        rzxpt(1) = rx-0.1*(rmhdgrid(2)-rmhdgrid(1))
        call dbcevl1 (rmhdgrid,nw,zmhdgrid,nh,cspln,
     .                nw,rzxpt(1),zx,pds,ier,1)
        psizxpt(1) = pds(1)
        drzx       = (rsepout-rzxpt(1))/(nrxpt-1)
        do j=2,nrxpt
           rzxpt(j) = rzxpt(j-1)+drzx
           if (j .eq. nrxpt)  rzxpt(j) = rsepout        ! avoid roundoff
           call dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,
     .                  nw,rzxpt(j),zx,pds,ier,1)
           psizxpt(j) = pds(1)
        end do
c
c --- get the inboard and outboard radius of limiter at elevation
c --- zsep0:
c
        zoutmin = zmhdgrid(nh)
        zinmin  = zmhdgrid(nh)
        do j=1,nlimiter
           if (100.0*xlimiter(j) .gt. xax(1)) then  ! get outboard value
              zout = ABS (100.*ylimiter(j)-zsep0)
              zoutmin = MIN (zoutmin,zout)
              if (zoutmin .eq. zout)  izout = j
           else
              zin = ABS (100.0*ylimiter(j)-zsep0)   ! get inboard value
              zinmin = MIN (zinmin,zin)
              if (zinmin .eq. zin)  izin = j
           end if
        end do
c
c       assume closest point is OK without interpolation
c
        rlimin = 100.0 * xlimiter(izin )
        rlimot = 100.0 * xlimiter(izout)
c
        call dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,
     .               nw,rlimin,zsep0,pds,ier,1)
        psilin = ABS (pds(1))*1.0e-8
        call dbcevl1 (rmhdgrid,nw,zmhdgrid,nh,cspln,
     .                nw,rlimot,zsep0,pds,ier,1)
c
        psilot  = ABS (pds(1)) * 1.0e-8
        rlimin  =  rlimin / 100.0
        rlimot  =  rlimot / 100.0
        rx      =      rx / 100.0
        zx      =      zx / 100.0
        psix    = ABS (psix) * 1.0e-8
        rsepin  = rsepin  / 100.0
        rsepout = rsepout / 100.0
        do j=1,nrxpt
          rzxpt(j)   = rzxpt(j) / 100.0
          psizxpt(j) = ABS (psizxpt(j)) * 1.0e-8
        end do
c
c --- now load the table pitchlos . Pitchloss(j,i) is the energy,in keV,
c --- corresponding to pitch angle pitchang(j) and major radius pitchrmaj(i).
c --- Note that this logic relies on the fact that
c --- subroutine OLOSSA returns the same vector of pitch angles, pitchang,
c --- each time it is called.
c
        do i=1,mf                  ! loop over possible values of rmajor
           r0           = rotsid(i)
           pitchrmaj(i) = r0 / 100.0  ! copy of rotsid for now in meters
           b0   = 1.0e-6 * ABS (fpsio(i))/pitchrmaj(i)  ! toroidal field
           psi0 = ABS (potsid(i))*1.0e-8                ! psi value
           call olossa(fm,fz,tmax,pitchrmaj(i),b0,psi0,rx,zx,psix,
     .                 rsepin,rsepout,rlimin,rlimot,psilin,psilot,
     .                 rzxpt,psizxpt,nrxpt,pitchang,pitcheng,npitch)
c
c          j => pitch angle, i => rmajor
c
           do j=1,npitch
             pitchlos(j,i) = pitcheng(j) / 1000.0    ! keV
           end do
        end do
c
c --- load bpeq
c
        cconst = 1.0e+3
        call multpl1 (psir, nj, cconst)     ! convert psir to  gauss-cm2
        call intrp   (0, 2, psir, bp, nj, potsid, bpeq, mf)  ! load bpeq
        call multpl1 (bpeq, mf, cconst)     ! convert bpeq to gauss
        cconst = 1.0e-3
        call multpl1 (psir, nj, cconst)     ! convert psir to kgauss-cm2
c
        losscalc = 0  ! don't do these calcs again until tport is called
        return        ! from RUNTWO (i.e., next equilibrium cycle)
      end if          ! losscalc = 1 branch
c
c --- evaluation section of subroutine (done when losscalc = 0)
c
      if (losseval .lt. 0) return
c
c --- determine the spline representation
c
      if (losseval .ne. 0) then
c
        call my_ibcccu (p,rmhdgrid,nw,zmhdgrid,nh,cspln,nw,
     .                   wnoperm,ier)
        if (ier .ne. 0)
     .    call STOP ('subroutine GET_ORBLOSS: IBCCCU problem #2', 168)
c
c --- get bp0 (in gauss):
c
        do j=2,nj
          bp0(j) = rbp(j)/(fcap(j)*gcap(j)*hcap(j)*r(j))
        end do
        bp0(1) = 0.0
c
c --- interpolate bp0 onto the psi grid FREYA zoned grid (stored in bpoltp):
c
        cconst = 1.0e3
        call multpl1(psir,nj,cconst)  ! temp convert psir to gauss-cm2
        call intrp(0,2,psir,bp0,nj,potsid,bpoltp,mf)       ! load bptolp
        cconst = 1.0e-3
        call multpl1(psir,nj,cconst) ! convert psir to kgauss-cm2
c
        losseval = 0                 ! don't do this on subsequent calls
c
      end if
c
c --- get radial and vertical components of poloidal b field:
      call dbcevl1(rmhdgrid,nw,zmhdgrid,nh,cspln,
     .                                      nw,rpos,zpos,pds,ier,3)
      bpz = pds(2)/rpos
      bpr = -pds(3)/rpos
c
c  compute accurate pitch angle                                  ! rs
c
      izp1    = izone+1                                          ! rsHSJ
      fpsia   = 0.5*(fpsio(izone)+fpsio(izp1))                   ! rsHSJ
      ravg    = 0.5*(pitchrmaj(izone)+pitchrmaj(izp1))           ! rsHSJ
      btrd    = fpsia / ravg                                     ! rsHSJ
      bpoltpa = 0.5*(bpoltp(izone)+bpoltp(izp1))                 ! rsHSJ
      bpeqa   = 0.5*(bpeq(izone)+bpeq(izp1))                     ! rsHSJ
c
c scale bp components from values stored at equilibrium call     ! rs
c bpoltp is the current value of flux-averaged bp                ! rs
c   in transport evolution                                       ! rs
c bpeq is the flux-averaged value of bp                          ! rs
c   stored at equilibrium call                                   ! rs
c
      bpz    = bpz*bpoltpa/bpeqa                                 ! rsHSJ
      bpr    = bpr*bpoltpa/bpeqa                                 ! rsHSJ
      cosphi = (btrd*vtroid+bpr*vrad+bpz*vz0) ! check this calc  ! rs
     .       / (SQRT (btrd**2+bpz**2+bpr**2)*vbeam)              ! rsHSJ
      phi    = (180.0 / pi) * ACOS (cosphi)                      ! rs
      if (phi .gt. 180 .or. phi .lt. 90)  return                 ! rsHSJ
      call olinterp (phi, pitchang, npitch, rpos, pitchrmaj, mf, ! rsHSJ
     .               pitchlos, 1, elossb, ier)                   ! rsHSJ
      return
c
      end

      subroutine getpow (time, pwrinpt)
c
      USE param
      USE nub
      USE rf
      USE nub4            !time_dep_beam
      USE transp,                       only    : beam_data,use_nubeam               
      USE constnts ,                    only    : pisq
      USE machin,                       ONLY    : rmajor
      USE numbrs,                       only    : nj
      USe mesh,                         only    : r
      USE sourc,                        only    : qrfe,qrfi,qohm
      USE geom,                         only    : hcap
      USE P_Nfreya_12_interface,        ONLY    : use_P_Nfreya       
      USE  neutral_beams,               ONLY    : nf_bptor => bptor

      implicit  integer (i-n), real*8 (a-h, o-z)

      real *8, allocatable,dimension(:)  :: prfle,prfli,prflo
      real *8 joupkev
c
c ----------------------------------------------------------------------
c this subroutine determines the external (beam and RF) input
c power at time t
c ------------------------------------------------------------------ HSJ
c

c
      allocate(prfle(nj),prfli(nj),prflo(nj))
      pwrinpt_rf = 0.0
      IF(use_P_Nfreya)bptor(:) = nf_bptor(:)

c
c --- check the RF models
c

c      do j=1,krf
c        factr   = turnonmod (time, rfon(j), turnonp(j), rfoff(j),
c     .              rftime(j),rframp_timeup(j),rframp_timedown(j))
c        pwrinpt_rf = factr * rfpow(j) + pwrinpt_rf
c        pwrinpt_rf = pwrinpt_rf + rf_ext_qetot(j) + rf_ext_qitot(j)
c      end do
c --- check external rf input

      joupkev = 1.60217733e-16
      const = 4.0*pisq*rmajor*joupkev
      call trap3(r,qrfe,nj,hcap,const,prfle)
      call trap3(r,qrfi,nj,hcap,const,prfli)
      call trap3(r,qohm,nj,hcap,const,prflo)
      pwrinpt_rf = prfle(nj)  + prfli(nj)


c------------------------------------------------------------------------------
c --- beam power assume full on or off,same for all beams
c------------------------------------------------------------------------------
       pwrinpt_beam =0.0
       beam_model :  IF(time_dep_beam .eq.1 .AND. 
     .                                      .NOT. use_P_Nfreya) THEN      
                     !pulsed beams any number of 
                     !pulses could currently be active
         do l=1,nbeams
            ll=l
            do i = 1,n21s
               ii=i
               do np =1,n_pulse
                  nnpp =np
                  ontime = pbeamOn(np,i,l) - time
                  if(ontime .gt. 0.0e0)go to 100 !these are future pulses
                  offtime = pbeamOff(np,i,l) - time
                  if(offtime .lt. 0.0e0)go to 90 !these are past  pulses
                  if(ontime .le. 0.0 .and. offtime .gt. 0.0)then !pulse is active
                     sfrac = sfrac1(l)
                     if(i .eq. 2) sfrac = 1.- sfrac !second source beam l
                     pwrinpt_beam =pwrinpt_beam + bptor(l)*sfrac
                     go to 100
                  endif
 90            enddo            ! end loop on pulse #
 100           continue
            enddo    ! end loop on source #
          enddo      ! end loop on beam #
      ELSEIF(time_dep_beam .eq.1 .AND.  use_P_Nfreya) THEN  beam_model
         ! parallel P_Nfreya model
         DO j =1,beam_data%nbeam
            IF(beam_data%beamlet_active(j) == 1) pwrinpt_beam =  
     .                                   pwrinpt_beam +bptor(j)
         ENDDO


      ELSEIF( .not. use_nubeam)then   beam_model ! old single pulse  beam
          factr=1.0
          if (time .le. beamon(1))  factr = 0.0  
          if (time .gt. beamon(1)+btime(1))  factr = 0.0  
          do j=1,nbeams
             pwrinpt_beam = factr*bptor(j)+pwrinpt_beam
          end do
      ELSE  beam_model  !use_nubeam
             pwrinpt_beam = beam_data%pwf_tot_intg
      ENDIF  beam_model




      !prf and poh are in MW
      pwrinpt = prfle(nj) + prfli(nj)   +pwrinpt_beam +prflo(nj)
c      print *,'pwrinpt in getpow',pwrinpt
c      print *,'prfle,prfli, ,poh,pbeam =',
c     .       prfle(nj),prfli(nj),prflo(nj),pwrinpt_beam

      deallocate(prfle,prfli,prflo)
      return
c
      end




      subroutine gmprd (a, b, r, n, m, l)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
****  dimension a(1),b(1),r(1)                                   ! HSJrs
      dimension a(*),b(*),r(*)                                   ! HSJrs
c
c     a = name of first input matrix
c     b = name of second input matrix
c     n = rows in a
c     r = name of output matrix
c     m = columns in a,  = rows in b
c     l = columns in b
c
      ir = 0
      ik = -m
c
      do 30 k=1,l
        ik = ik+m
        do 20 j=1,n
          ir = ir+1
          ji = j-n
          ib = ik
          r(ir) = 0
          do 10 i=1,m
            ji = ji+n
            ib = ib+1
            r(ir) = r(ir)+a(ji)*b(ib)
   10     continue
   20   continue
   30 continue
c
      return
c
      end

      subroutine lagrange (x, f, n, sd)
c
c --- S.J. Thompson --- General Atomics --- Core Physics Fusion Group --
c
c                           SUBPROGRAM DESCRIPTION:
c                           ----------------------
c
c    Given arrays x and f of length n containing a function f(i) = f(x(i))
c    (i = 1, ..., n), this routine returns an array of second derivatives
c    sd, also of length n. These derivatives are the finite difference ap-
c    proximation using quadratic Lagrange interpolation. The result is a
c    constant on the interval (x(i-1), x(i+1)):
c
c                      2 * f(i-1)            2 * f(i)
c    f"(x(i)) =  ---------------------- - ------------- +
c                h(i) * [h(i) + h(i+1)]   h(i) * h(i+1)
c
c                       2 * f(i+1)
c                ------------------------
c                h(i+1) * [h(i) + h(i+1)]
c
c    which for equal intervals becomes
c
c         f(i-1) - 2 * f(i) + f(i+1)
c    f" = -------------------------- .
c                   h * h
c
c    Here we have used the abbreviations
c
c                h(i) = x(i) - x(i-1)
c                h(i+1) = x(i+1) - x(i)
c
c    and
c
c                f(i-1) = f(x(i-1))
c                f(i) = f(x(i))
c                f(i+1) = f(x(i+1)).
c
c    Reference: "Numerical Methods for Engineering Applications," pp. 52,
c                Joel H. Ferziger, John Wiley & Sons, Inc., 1981
c
c --- S.J. Thompson ---------------------------------------------- 18 Apr 94 ---
c
      implicit none
c
      integer  i, n
      real*8   f(n), ha, hm, hp, ht, x(n), sd(n)
c
      do i=2, n-1
        hm    = x(i) - x(i-1)
        hp    = x(i+1) - x(i)
        ha    = hm + hp
        ht    = hm * hp
        sd(i) = 2.0 * ((f(i-1) / hm + f(i+1) / hp) / ha - f(i) / ht)
      end do
c
      sd(1) = sd(2)
      sd(n) = sd(n-1)
      return
c
      end

      subroutine mint1 (dr, jb, jmix, ntab, rb, r, rxtab, ub, u, xx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  This subroutine integrates u vs. rxtab before mixing.
c
      dimension  dr(*), r(*), rxtab(*), u(*), xx(*)
c
c  set limits of do loop based upon whether rxtab(i) values increase
c     or decrease with increasing i
c
      if (rxtab(ntab) .lt. rxtab(1))  go to 10
      i1 =  1
      i2 = ntab - 1
      i3 =  1
      ip =  1
      go to 20
c
   10 i1 = ntab - 1
      i2 =  1
      i3 = -1
      ip =  0
c
c  perform integration
c
   20 do 60 i=i1,i2,i3
        xx(i) = 0.0
        do j=jb,jmix
          ra = rb
          ua = ub
          if (r(j+1) .gt. rxtab(i+ip))  go to 30
          rb = r(j+1)
          ub = u(j+1)
          go to 40
c
   30     jb = j
          rb = rxtab(i+ip)
          ub = u(j) + (u(j+1)-u(j))*(rb-r(j))/dr(j)
c
   40     xx(i) = xx(i) + (ra*ua+rb*ub)*(rb-ra)
          if (r(j+1) .gt. rxtab(i+ip))  go to 60
        end do
c
   60 continue
c
      return
c
      end

      subroutine mint2 (dr, j0, jmix, r, rtab, u, uint)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  This subroutine integrates u*r over the mixing region.
c
      dimension dr(*), r(*), u(*)
c
      uint = 0.0
      do 10 j=j0+1,jmix
        if (j .eq. j0+1 .and. rtab .ne. 0.0)  go to 10
        uint = uint + (r(j-1)*u(j-1)+r(j)*u(j))*dr(j-1)
   10 continue
      if (rtab .eq. 0.0)  go to 20
      uint = uint + r(j0+1)*u(j0+1)*dr(j0)
   20 uint = uint + r(jmix)*u(jmix)*dr(jmix)
      return
c
      end

      subroutine minv (a, n, ieror)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
****  dimension a(1),l(10),m(10)                                 ! HSJrs
      dimension  a(*), l(10), m(10)                              ! jf
c
c     a= input matrix (replaced by inverse)
c     n= order of matrix a
c       ieror will be 1 on return if matrix was singular
c     d= resultant determinant
c     l,m= work vector of length n
c
c     search for largest element
c
      d = 1.0
      nk = -n
      do 80 k=1,n
      nk = nk+n
      l(k) = k
      m(k) = k
      kk = nk+k
      biga = a(kk)
      do 20 j=k,n
      iz = n*(j-1)
      do 20 i=k,n
      ij = iz+i
   10 if (ABS (biga) - ABS (a(ij)))  15, 20, 20
   15 biga = a(ij)
      l(k) = i
      m(k) = j
   20 continue
c
c     interchange rows
c
      j = l(k)
      if (j - k)  35, 35, 25
   25 ki = k-n
      do 30 i=1,n
      ki = ki+n
      hold = -a(ki)
      ji = ki-k+j
      a(ki) = a(ji)
   30 a(ji) = hold
c
c     interchange columns
c
   35 i = m(k)
      if (i - k)  45, 45, 38
   38 jp = n*(i-1)
      do 40 j=1,n
      jk = nk+j
      ji = jp+j
      hold = -a(jk)
      a(jk) = a(ji)
   40 a(ji) = hold
c
c     divide column by minus pivot (contained in biga)
c
   45 if (biga)  48, 46, 48
   46 d = 0.0
      ieror = 1
      return
   48 do 55 i=1,n
      if (i - k)  50, 55, 50
   50 ik = nk+i
      a(ik) = a(ik)/(-biga)
   55 continue
c
c     reduce matrix
c
      do 65 i=1,n
      ik = nk+i
      hold = a(ik)
      ij = i-n
      do 65 j=1,n
      ij = ij+n
      if (i - k)  60, 65, 60
   60 if (j - k)  62, 65, 62
   62 kj = ij-i+k
      a(ij) = hold*a(kj)+a(ij)
   65 continue
c
c     divide row by pivot
c
      kj = k-n
      do 75 j=1,n
      kj = kj+n
      if (j - k)  70, 75, 70
   70 a(kj) = a(kj)/biga
   75 continue
c
c     product of the pivots
c
      d = d * biga
c
c     replace pivot by reciprocal
c
      a(kk) = 1.0 / biga
   80 continue
c
c     final interchange of rows and columns
c
      k = n
  100 k = (k-1)
      if (k)  150, 150, 105
  105 i = l(k)
      if (i - k)  120, 120, 108
  108 jq = n*(k-1)
      jr = n*(i-1)
      do 110 j=1,n
      jk = jq+j
      hold = a(jk)
      ji = jr+j
      a(jk) = -a(ji)
  110 a(ji) = hold
  120 j = m(k)
      if (j - k)  100, 100, 125
  125 ki = k-n
      do 130 i=1,n
      ki = ki+n
      hold = a(ki)
      ji = ki-k+j
      a(ki) = -a(ji)
  130 a(ji) = hold
      go to 100
  150 return
c
      end

      complex*16 function mixdis (gamma)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  This function evaluates the dispersion relation for the growth rate
c     gamma of the m = 1 tearing mode.
c
      common /dispar/ gamt, wse, wsi
      complex*16      gamma, xi
      data            xi /(0.0, 1.0)/
c
      mixdis = gamma*(gamma+xi*wse)*(gamma-xi*wsi) - gamt**3
      return
c
      end

      subroutine mix (atw, btor, curden, dt, ene, en, enalp, enb,
     .                etor, fcap, gcap, hcap,
     .                nbeams, nion, nj, nmix,  ppb, psis, q,
     .                r, rbp, rmajor, te, ti, time, walp, wb, zeff,
     .                curp, enep, enp, enap, enbp, esaw, etorp, ppbp,
     .                psisp, qp, qmag, qsawe, qsawi, rbpp,
     .                rmixxx, rsxxx, wap, wbp)
c
      USE param
      USE mixcom
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- old prescription:
c     subroutine mix(atw, btor, curden, dt, ene, en, enalp, enb,
c    *         epste, epsti, etor, fcap, gcap, hcap, ipmix, mixpro,
c    *         nbeams, nion, nj, nmix, ntab, ppb, psis, q, qmix,
c    *         rsmixx, r, rbp, rmajor, te, ti, time, tsmix,tdmix,
c    *         walp, wb, w0mix, w1mix, w2mix, w3mix, w5mix, zeff,
c    *         curp, enep, enp, enap, enbp, esaw, etorp, ppbp,
c    *         psisp, qp, qmag, qsawe, qsawi, rbpp, rmixx, rsx,
c    *         tem, tep, tim, tip, timmix, wap, wbp)
c
c --- the old CRAY CIVIC compiler only allowed 60 arguments in a subroutine call
c --- consequently the parameters in INCLUDE file mixcom which were previously
c --- passed through the argument list are now included directly in
c --- INCLUDE file mixcom  HSJ 11/27/90
c --- (NOTE: old arguments rmixx and rsx were replaced with rmixxx and
c --- rsxxx in this subroutine so that they won't clash with rmixx and rsx
c --- in INCLUDE file mixcom. the assignment of rmixxx to rmixx and rsxxx to
c --- rsx is done outside this subroutine (in TPORT)).
c
c  in this subroutine we use posbtor = ABS (btor) for calculations
c  btor is normally negative on the eqdsk.   HSJ
c  This subroutine models the MHD mixing associated with internal
c     (sawtooth) disruptions.  The redistribution or mixing of current,
c     particles, and energy is controlled by the following flag:
c       mixpro    2, Generate mixed profiles consistent with helical
c                    flux conservation as suggested by Kadomtsev for
c                    a single-valued q profile and since extended by
c                    Parail and Pereverzev to the case of a double-
c                    valued q profile.  Mix temperatures instead of
c                    energy densities; then renormalize to conserve
c                    energy.
c                 1, Generate mixed profiles consisitent with helical
c                    flux conservation as above, but mix energy densi-
c                    ties.  This can lead to unphysical temperature
c                    profiles if the density profiles are peaked and
c                    not mixed.
c                 0, Generate mixed profiles for densities and temper-
c                    atures that are flat.
c                -1, Generate mixed profiles for densities and energies
c                    that are flat.
c     The various mechanisms for triggering the mixing are as follows:
c       qmix>0    Disruption is triggered when qmin< = qmix.
c       rsmixx>0  Disruption is triggered when rsxxx> = rsmixx.
c       tsmix>0   Disruption is triggered when ts> = tsmix.
c       tdmix>0   Two disruptions are triggered depending upon
c                 the extent of the mixing regions.
c       ipmix>0   Disruption is triggered when psis(rs2)>0
c                 as suggested by Parail and Pereverzev.
c       w0mix>0   Disruption is triggered when the width(s) of the
c                 m = 1 island(s) is(are) sufficiently large.
c                 As suggested by Jahns et al. the island width is
c                 calculated by integrating in time the island growth
c                 rate, which is obtained by solving the dispersion
c                 relation of Waddell et al.
c  References
c     1.  B. B. Kadomtsev, Sov. J. Plasma Phys. 1, 389 (1975).
c     2.  V. V. Parail and G. V. Pereverzev, Sov. J. Plasma Phys. 6,
c         14 (1980).
c     3.  G. L. Jahns, M. Soler, B. V. Waddell, J. D. Callen, and
c         H. R. Hicks, Nucl. Fusion 18, 609 (1978).
c     4.  B. V. Waddell, G. Laval, M. Rosenbluth, B. Coppi, and
c         B. Basu, ORNL/TM-5968, Oak Ridge National Laboratory (1977).
c
c      include 'param.i'
c
      parameter (numbrislnds = 3, numbrislndsp1 = numbrislnds+1)
c
      include 'mixx.i'        ! dr and bp are stored here (remove later)
c      include 'mixcom.i'
c
      common /dispar/ gamt, wse, wsi
c     dimension atw(kion), r(kj), tdmix(2), zeff(kj)
      dimension atw(kion), r(kj), zeff(kj)
      dimension curden(kj), ene(kj), en(kj,kion), enalp(kj),
     .          enb(kj,ke,kb), etor(kj), fcap(kj), gcap(kj), hcap(kj),
     .          psis(kj), rbp(kj), q(kj), te(kj), ti(kj),
     .          walp(kj), wb(kj,ke,kb), ppb(kj,ke,kb)
****  dimension curp(kj), enep(kj), enp(kj,kion), enap(kj),
**** .          enbp(kj,ke,kb), etorp(kj), psisp(kj), rbpp(kj),
**** .          qp(kj), tem(kj), tep(kj), tim(kj), tip(kj),
**** .          wap(kj), wbp(kj,ke,kb), ppbp(kj,ke,kb),
**** .          esaw(kj), qmag(kj), qsawe(kj), qsawi(kj)
      dimension curp(kj), enep(kj), enp(kj,kion), enap(kj),
     .          enbp(kj,ke,kb), etorp(kj), psisp(kj), rbpp(kj),
     .          qp(kj), wap(kj), wbp(kj,ke,kb), ppbp(kj,ke,kb),
     .          esaw(kj), qmag(kj), qsawe(kj), qsawi(kj)
      dimension ptab(ktab), rtab(ktab), rxtab(ktab,numbrislndsp1)
      dimension xx(ktab,numbrislndsp1), rtaba(ktab), utaba(ktab)
      dimension js(numbrislnds), psisx(numbrislnds), rs(numbrislnds)
      dimension infer(numbrislnds)
      real*8    ten_thousand
      complex*16 mixdis, gamma(numbrislnds)
      external  mixdis
      data      pi /3.141592654/, amu /1.25664e-6/
      data      ifirst /1/, im /2/
c
c     zero out miscellaneous parameters
c
      zero = 0.0
      ntab = ktab           ! ktab is set in param.i
      posbtor = ABS (btor)
      r0 = 0.0
      rmix = 0.0
      do 5 is=1,numbrislnds
        rs(is) = 0.0
    5 continue
c      print *,'qmix =',qmix
c
c  initialize timmix and width; write heading for printout
c
      wid = w0mix
      ip = 1
      if (ifirst .eq. 0)  go to 10
      ifirst = 0
      timmix = time-dt
      wid = w0mix
      write (nmix, 1105)
c
c  calculate bp, dr, and psis
c  note that psi = poloidal flux/twopi and,if we define ,
c  PHI = toroidal flux /twopi then PHI=0.5*r(j)**2*bt0
c  hence psis = (poloidal flux -toroidal flux)/twopi  
c  integrate Bp0 = (1./R0)(d psi/ drho) to get  psi
c   and rho = SQRT(toroidal flux/(pi*btor)) to get PHI  :          HSJ
c
   10 bp(1)   = 0.0
      psi     = 0.0
      psis(1) = 0.0
c
      do j=2,nj
        bp(j)   = rbp(j) / (fcap(j)*gcap(j)*hcap(j)*r(j))
        dr(j-1) = r(j) - r(j-1)
        psi     = psi + rmajor*0.5*(bp(j-1)+bp(j))*dr(j-1)
        psis(j) = psi - 0.5*r(j)**2*posbtor  !see note above this is not the
                                             !helical flux ?? HSJ
c        print *,'psis(j) =',psis(j)
      end do
c
c set the final profiles equal to the initial ones, except for te and ti
c
      do 18 j=1,nj
      bpp(j) = bp(j)
      curp(j) = curden(j)
      enep(j) = ene(j)
      enap(j) = enalp(j)
      etorp(j) = etor(j)
      psisp(j) = psis(j)
      rbpp(j) = rbp(j)
c     qp(j)   = q(j)
      qp(j)   = ABS (q(j))    ! for MHD version HSJ
      wap(j)  = walp(j)
      do 16 k=1,nion
   16 enp(j,k) = en(j,k)
      if (w5mix .eq. 0.0)  go to 18
      do 17 ib=1,nbeams
      do 17 ie=1,3
      enbp(j,ie,ib) = enb(j,ie,ib)
      wbp(j,ie,ib) = wb(j,ie,ib)
   17 ppbp(j,ie,ib) = ppb(j,ie,ib)
   18 continue
c
c  calculate minimum safety factor q
c
c     qmin = q(1)
      qmin = ABS (q(1))                ! for MHD version HSJ
      do 20 j=2,nj
** 20 qmin = MIN (q(j),qmin)
   20 qmin = MIN (ABS (q(j)), qmin)    ! for MHD version HSJ
c      print *,'qmin =',qmin            !HSJ  06/07/02
c
c  skip further calculations if q>0.998 everywhere; this reduces the
c     possibility of triggering a disruption due to numerical noise
c
      if (qmin .gt. 0.998)  go to 800
c
c  find each rs such that q(rs) = 1; determine number of islands
c
      is = 0
      do 40 j=1,nj-1
****  if (q(j) .ne. 1.0)  go to 30
      if (ABS (q(j)) .ne. 1.0)  go to 30    ! for MHD version HSJ
      is = is+1
      if (is .gt. numbrislnds)  go to 40
      js(is) = j
      rs(is) = r(j)
      go to 40
** 30 sign = (1.0-q(j))*(1.0-q(j+1))
   30 sign = (1.0 - ABS (q(j)))*(1.0 - ABS (q(j+1)))   ! for MHD version
      if (sign .ge. 0.0)  go to 40
      is = is+1
      if (is .gt. numbrislnds)  go to 40
      js(is) = j
****  rs(is) = r(j) + (r(j+1)-r(j))*(1.0-q(j))/(q(j+1)-q(j))
      qj = ABS (q(j))
      qjp1 = ABS (q(j+1))
      rs(is) = r(j) + (r(j+1)-r(j))*(1.0-qj)/(qjp1-qj)
   40 continue
      islnd = is
      if (islnd .gt.  numbrislnds) ip = 2
c      print *,'ip,is,islnd =',ip,is,islnd
c      print *,'rs(1,2.3) =',rs(1),rs(2),rs(3)
c      print *,'js(1-3) =',js(1),js(2),js(3)
      if (ip .eq. 2)  go to 800
c
c  skip further calculations if rs values are too close together
c
      if (islnd .eq. 2 .and. js(2)-js(1) .le. 1) ip = 3
      if (islnd .eq. 3 .and. js(3)-js(2) .le. 1) ip = 3
c      print *,'ip 3 =',ip
      if (ip .eq. 3)  go to 800
c
c  calculate psis(rs) using quadratic interpolation
c
      do 80 is=1,3
      if (rs(is) .eq. 0.0)  go to 80
      jx = js(is)
      r21 = r(jx+1)-r(jx)
      r31 = r(jx+2)-r(jx)
      r32 = r(jx+2)-r(jx+1)
      cx = (r32*psis(jx)-r31*psis(jx+1)+r21*psis(jx+2)) / (r21*r32*r31)
      bx = (psis(jx+1)-psis(jx)-cx*(r(jx+1)**2-r(jx)**2)) / r21
      ax =  psis(jx)-bx*r(jx)-cx*r(jx)**2
      dx = -cx
      rsxxx = bx/(2.0*dx)
      psisx(is) = ax + dx*rsxxx**2
   80 continue
c      print *,'psisx(1-3) =',psisx(1),psisx(2),psisx(3)
c      print *,'js(1-3)',js(1),js(2),js(3)
c      print *,'rsmixx =',rsmixx
c
c --- now have helical flux,psisx(i) on rs(i) grid,and  js(i),
c --- where r(js(i)) .le. rs(i)
c --- i = 1,..islnd ,where islnd .le. nmbrislnds
c
c  test trigger on safety factor
c
      if (qmin .le. qmix)  go to 201
c
c  test trigger on radius of q = 1 singular surface
c
      if (rsmixx .eq. 0.0)  go to 120
      if (rs(islnd)/r(nj) .ge. rsmixx)  go to 201
  120 continue
c
c  test trigger on sawtooth period
c
      if (tsmix .eq. 0.0)  go to 130
      if (tdmix(1) .ne. 0.0)  go to 124
      if (time-timmix .ge. tsmix-1.0e-6)  go to 201
      go to 130
  124 if (im .eq. 2)  go to 128
      if (islnd .gt. 1)  go to 126
      if (time-timmix .ge. tsmix-1.0e-6)  go to 201
      go to 130
  126 if (time-timmix .ge. tdmix(1)-1.0e-6)  go to 201
      go to 130
  128 if (time-timmix .ge. tdmix(2)-1.0e-6)  go to 201
  130 continue
c
c  test trigger on psis
c
      if (ipmix .eq. 0)  go to 135
      if (islnd .eq. 2 .and. psisx(2) .gt. 0.0)  go to 201
      if (islnd .eq. 3 .and. psisx(3) .gt. psisx(1))  go to 201
  135 continue
c
c  determine whether to calculate width of outermost m = 1 island
c
c      print *,'w0mix  =',w0mix
      if (w0mix .eq. 0.0)  go to 800
c
c  calculate various quantities at the outermost rs
c
      is = islnd
      jx = js(is)
      enes = 0.5*(ene(jx)+ene(jx+1))
      tes = 0.5*(te(jx)+te(jx+1))
      zeffs = 0.5*(zeff(jx)+zeff(jx+1))
      ensum = 0.0
      do 140 k=1,kion
      ens(k) = 0.5*(en(jx,k)+en(jx+1,k))
  140 ensum = ensum + ens(k)
c
c  calculate taur, taua, and s; use MKS formulas
c
      xlam = 24.0 - LOG (1.0e-3 * SQRT (enes)/tes)
      eta = 3.28e-9*zeffs*xlam*(0.29+0.46/(1.08+zeffs))/tes**1.5
      taur = amu*(0.01*rs(is))**2/eta
      rho = 0.0
      do 150 k=1,kion
  150 rho = rho + atw(k)*1.6726231e-27*1.0e6*ens(k)
      taua = 0.01*rmajor * SQRT (amu*rho)/(1.0e-4*posbtor)
      s = taur/taua
c
c  calculate alpha and gamt
c
      alpha = rs(is)*(q(jx+1)-q(jx))/dr(jx)
      alpha = ABS (alpha)
      gamt = (alpha*s)**0.666667/taur
c
c  calculate wse and wsi; use MKS formulas
c
      dtedr = 1000.0*(te(jx+1)-te(jx))/(0.01*dr(jx))
      dpedr = 1.0e9*(ene(jx+1)*te(jx+1)-ene(jx)*te(jx))/(0.01*dr(jx))
      dpidr = 0.0
      do 160 k=1,kion
  160 dpidr = dpidr + 1.0e9*(en(jx+1,k)*ti(jx+1)-en(jx,k)*ti(jx))
     .                /(0.01*dr(jx))
      wse = (dpedr+0.71*1.0e6*enes*dtedr)/(enes*rs(is)*posbtor)
      wsi = dpidr/(ensum*rs(is)*posbtor)
      wse = ABS (wse)
      wsi = ABS (wsi)
c
c  calculate gam by solving dispersion relation
c
      ten_thousand = 1.0e-4
      call zanlyt2 (mixdis, ten_thousand, 4, 0, 0, numbrislnds,
     .              gamma, 20, infer, ier)
      rt1 = REAL (gamma(1))
      rt2 = REAL (gamma(2))
      rt3 = REAL (gamma(3))
      gam =  MAX (rt1, rt2, rt3)
c
c  calculate width of outermost m = 1 island
c
      wid = wid * EXP (gam*dt)
c
c  test trigger on width of outermost m = 1 island
c
      if (islnd .eq. 1 .and. wid .ge. rs(1)      )  go to 201
      if (islnd .eq. 2 .and. wid .ge. rs(2)-rs(1))  go to 201
      if (islnd .eq. 3 .and. wid .ge. rs(3)-rs(2))  go to 201
      ip = 5
      go to 800
c
c  calculate r0 and rmix
c
  201 is = islnd
      p0 = psisx(is)
c      print *,'psisx(is) =',psisx(is)
c      print *,'islnd =',islnd
      if (islnd .eq. 1)  pmix = 0.0
      if (islnd .eq. 2)  pmix = psisx(1)
      if (islnd .eq. 3)  pmix = MIN (psisx(2), zero)
      if (islnd .eq. 3 .and. psisx(3) .lt. psisx(1)) pmix = psisx(2)
      j0 = 1
      r0 = 0.0
      j1 = 1
      if (islnd .eq. 3 .and. psisx(3) .lt. psisx(1)) j1 = js(1)
      if (islnd .eq. 1)  go to 205
      if (islnd .eq. 2 .and. psisx(2) .ge. 0.0)  go to 205
      if (islnd .eq. 3 .and. psisx(3) .ge. psisx(1))  go to 205
      do 202 j=j1,nj-1
        if (psis(j+1) .lt. p0)  go to 204
  202 continue
  204 j0 = j
      r0 = r(j0) + dr(j0)*(p0-psis(j0))
     .             /(psis(j0+1)-psis(j0))
c
  205 do 206 j=j1,nj-1
      jmix = j
      if (psis(j+1) .lt. pmix)  go to 208
  206 continue
**208 jmix = j
  208 rmix = r(jmix) + dr(jmix)*(pmix-psis(jmix)) /
     .                  (psis(jmix+1)-psis(jmix))
c
c  skip further calculations if rmix> = r(nj) or rmix<=rs(island)
c

      if (rmix .ge. r(nj) .or. rmix .le. rs(islnd)) ip = 4
c      print *,'ip,rmix, r(nj)= ',ip,rmix,r(nj)
c      print *,'rs(islnd) =',rs(islnd)
c      print *,'j1,jmix =',j1,jmix
c      print *,'r(jmix) =',r(jmix)
c      print *,'dr(jmix),pmix,psis(jmix)',dr(jmix),pmix,psis(jmix)
c      print *,'psis(jmix+1)=',psis(jmix+1)
      if (ip .eq. 4)  go to 800
c
c  test on r0 and set im when tdmix(1) .ne. 0
c
      if (tdmix(1) .eq. 0.0)  go to 210
c      print *,'islnd, r0',islnd,r0
      if (islnd .gt. 1 .and. r0 .eq. 0.0)  go to 800
      im = 1
      if (r0 .gt. 0.0) im = 2
  210 continue
c
c  set up ptab, a table of equally spaced psis values
c
      do 215 i=1,ntab
  215 ptab(i) = pmix + (ntab-i)*(p0-pmix)/(ntab-1)
c
c  interpolate to obtain r values (stored in rxtab) corresponding to
c     psis values (stored in ptab) before mixing;
c     rxtab(i,1) is inside rs(1),
c     rxtab(i,2) is between rs(1) and rs(2),
c     rxtab(i,3) is between rs(2) and rs(3), and
c     rxtab(i,4) is outside rs(3)
c
      if (islnd .eq. 1) ifor = -1
      if (islnd .gt. 1) ifor = 1
      if (islnd .eq. 3 .and. r0 .eq. 0.0) ifor = -1
      ib = 0
      if (r0 .ne. 0.0) ib = islnd-2
      pa = 0.0
      if (r0 .ne. 0.0) pa = p0
      call mterp(ifor, j0, js(ib+1), ntab, pa, psisx(ib+1), psis, ptab,
     .           r0, rs(ib+1), r, rxtab(1,ib+1))
      if (islnd .eq. 1)  go to 220
      ifor = -ifor
      ib = ib+1
      call mterp(ifor, js(ib), js(ib+1), ntab, psisx(ib), psisx(ib+1),
     .           psis, ptab, rs(ib), rs(ib+1), r, rxtab(1,ib+1))
      if (islnd .eq. 2 .or. r0 .ne. 0.0)  go to 220
      ifor = -ifor
      ib = ib+1
      call mterp(ifor, js(ib), js(ib+1), ntab, psisx(ib), psisx(ib+1),
     .           psis, ptab, rs(ib), rs(ib+1), r, rxtab(1,ib+1))
  220 ifor = -ifor
      ib = ib+1
      call mterp(ifor, js(ib), jmix+1, ntab, psisx(ib), pmix,
     .           psis, ptab, rs(ib), rmix, r, rxtab(1,ib+1))
c
c  calculate rtab, a table of r values after mixing
c
      rtab(1) = r0
      do 260 i=2,ntab-1
      rx2 = rxtab(i,is+1)**2 - rxtab(i,is)**2
      if (islnd .gt. 1) rx2 = rx2 + rxtab(i,is-1)**2
      if (islnd .eq. 3 .and. r0 .eq. 0.0) rx2 = rx2 - rxtab(i,1)**2
  260 rtab(i) = SQRT (rx2)
      rtab(ntab) = rmix
c
c  interpolate to obtain psis after mixing; use linear interpolation
c     except inside rtab(2), where quadratic interpolation is employed
c
      if (r0 .eq. 0.0) psisp(1) = ptab(1)
      do 285 j=j0+1,jmix
      if (r0 .ne. 0.0)  go to 270
      if (r(j) .gt. rtab(2))  go to 270
      psisp(j) = ptab(1) + (ptab(2)-ptab(1))*(r(j)/rtab(2))**2
      go to 285
  270 continue
      do 275 i=2,ntab-1
      if (r(j) .lt. rtab(i+1))  go to 280
  275 continue
      i = ntab-1
  280 psisp(j) = ptab(i) + (ptab(i+1)-ptab(i))*(r(j)-rtab(i))
     .                     /(rtab(i+1)-rtab(i))
  285 continue
c
c  calculate bp, rbp, and q after mixing; use quadratic extrapolation
c     on axis
c
      do 310 j=j0+1,jmix
      bpp(j) = ((psisp(j+1)-psisp(j-1))/(r(j+1)-r(j-1)) + r(j)*posbtor)
     .         /rmajor
      rbpp(j) = fcap(j)*gcap(j)*hcap(j)*r(j)*bpp(j)
  310 qp(j) = r(j)*posbtor/(rmajor*bpp(j))
      if (r0 .eq. 0.0)
     .  qp(1) = (qp(2)*r(3)**2-qp(3)*r(2)**2)/(r(3)**2-r(2)**2)
c
c  calculate current density after mixing; use quadratic extrapolation
c     on axis
c
      do 340 j=j0+1,jmix
      drbp = rbpp(j)/fcap(j)-rbpp(j-1)/fcap(j-1)
  340 curp(j) = (2.0e-6*drbp/(amu*dr(j-1))
     .          - hcap(j-1)*r(j-1)*curp(j-1)) / (hcap(j)*r(j))
      if (r0 .eq. 0.0)
     .  curp(1) = (curp(2)*r(3)**2-curp(3)*r(2)**2)/(r(3)**2-r(2)**2)
c
c  mix the thermal particle densities while conserving particles
c
      mx = mixpro
      jx = js(islnd)
      if (w1mix .eq. 0.0)  go to 450
      do 420 j=j0,jmix+1
      w(j) = 1.0
  420 wp(j) = 1.0
      call mixu(dr, hcap, is, j0, jmix, jx, ktab, mx, ntab, r,
     .          rtab, rxtab, ene, enep, w, wp, xx, rtaba, utaba)
      do 440 k=1,nion
  440 call mixu(dr, hcap, is, j0, jmix, jx, ktab, mx, ntab, r,
     .          rtab, rxtab,
     .          en(1,k), enp(1,k), w, wp, xx, rtaba, utaba)
  450 continue
c
c  mix the electron temperature while conserving electron energy
c
      do 510 j=1,nj
      tem(j) = te(j)
      tep(j) = te(j)
      qsawe(j) = 0.0
      w(j) = ene(j)
  510 wp(j) = enep(j)
      if (w2mix .eq. 0.0)  go to 550
      if (w2mix .gt. 0.0)  go to 530
      do 520 j=1,jmix
      tem(j) = te(j) + epste*(1.0-2.0*(r(j)/rmix)**2)
  520 tem(j) = MAX (tem(j),te(nj))
  530 call mixu(dr, hcap, is, j0, jmix, jx, ktab, mx, ntab, r,
     .          rtab, rxtab, tem, tep, w, wp, xx, rtaba, utaba)
      do 540 j=1,jmix
  540 qsawe(j) = 1.5*(enep(j)*tep(j)-ene(j)*tem(j))/(time-timmix)
  550 continue
c
c  mix the ion temperature while conserving ion energy
c
      do 610 j=1,nj
      tim(j) = ti(j)
      tip(j) = ti(j)
      qsawi(j) = 0.0
      w(j) = 0.0
      wp(j) = 0.0
      do 610 k=1,nion
      w(j) = w(j) + en(j,k)
  610 wp(j) = wp(j) + enp(j,k)
      if (w3mix .eq. 0.0)  go to 650
      if (w3mix .gt. 0.0)  go to 630
      do 620 j=1,jmix
      tim(j) = ti(j) + epsti*(1.0-2.0*(r(j)/rmix)**2)
  620 tim(j) = MAX (tim(j),ti(nj))
  630 call mixu(dr, hcap, is, j0, jmix, jx, ktab, mx, ntab, r,
     .          rtab, rxtab, tim, tip, w, wp, xx, rtaba, utaba)
      do 640 j=1,jmix
  640 qsawi(j) = 1.5*(wp(j)*tip(j)-w(j)*tim(j))/(time-timmix)
  650 continue
c
c  mix the fast ion density, energy, and parallel momentum
c
      if (w5mix .eq. 0.0)  go to 670
      do 660 j=j0,jmix+1
      w(j) = 1.0
  660 wp(j) = 1.0
      call mixu(dr, hcap, is, j0, jmix, jx, ktab, mx, ntab, r,
     .          rtab, rxtab, enalp, enap, w, wp, xx, rtaba, utaba)
      call mixu(dr, hcap, is, j0, jmix, jx, ktab, mx, ntab, r,
     .          rtab, rxtab, walp, wap, w, wp, xx, rtaba, utaba)
      do 665 ib=1,nbeams
      do 665 ie=1,3
      call mixu(dr, hcap, is, j0, jmix, jx, ktab, mx, ntab, r,
     .          rtab, rxtab,
     .          enb(1,ie,ib), enbp(1,ie,ib), w, wp, xx, rtaba, utaba)
      call mixu(dr, hcap, is, j0, jmix, jx, ktab, mx, ntab, r,
     .          rtab, rxtab,
     .          wb(1,ie,ib), wbp(1,ie,ib), w, wp, xx, rtaba, utaba)
      call mixu(dr, hcap, is, j0, jmix, jx, ktab, mx, ntab, r,
     .          rtab, rxtab,
     .          ppb(1,ie,ib), ppbp(1,ie,ib), w, wp, xx, rtaba, utaba)
  665 continue
  670 continue
c
c  calculate approximate electric field after mixing
c
      do 680 j=1,jmix
  680 etorp(j) = (te(j)/tep(j))**1.5*etor(j)
c
c  calculate the change in electric field asociated with mixing
c
      do 685 j=nj,1,-1
      if (j .le. jmix)  go to 682
      esaw(j) = 0.0
      go to 685
  682 esaw(j) = esaw(j+1) - 0.5e-8*(bpp(j)-bp(j)+bpp(j+1)-bp(j+1))*dr(j)
     .                      /(time-timmix)
  685 continue
c
c  calculate the change in magnetic energy associated with mixing;
c     add energy to electrons if w2mix>0
c
      xme = 0.0
      xmep = 0.0
      xint = 0.0
      do 690 j=2,jmix+1
      xme = xme + (hcap(j-1)*r(j-1)*bp(j-1)**2
     .            +hcap(j)*r(j)*bp(j)**2)*dr(j-1)
      xmep = xmep + (hcap(j-1)*r(j-1)*bpp(j-1)**2
     .              +hcap(j)*r(j)*bpp(j)**2)*dr(j-1)
  690 xint = xint + (hcap(j-1)*r(j-1)+hcap(j)*r(j))*dr(j-1)
      qmagav = (xme-xmep)/(xint*8.0 * pi*1.6e-9*(time-timmix))
      do 692 j=1,nj
  692 qmag(j) = 0.0
      do 695 j=j0,jmix
      if (j .eq. j0 .and. rtab(1) .ne. 0.0)  go to 695
      tep(j) = tep(j) + qmagav*(time-timmix)/(1.5*enep(j))
      qmag(j) = qmagav
  695 continue
      ip = 5
      go to 800
c
c  write out summary information
c
  800 rs1x  = rs(1)/r(nj)
      rs2x  = rs(2)/r(nj)
      rs3x  = rs(3)/r(nj)
      rsxxx   = MAX (rs1x,rs2x,rs3x)
      r0x   = r0/r(nj)
      rmixxx = rmix/r(nj)
      widx  = wid/r(nj)
      qaxis = ABS (q(1))                    ! for MHD version HSJ
      go to (810, 820, 830, 840, 850), ip
  810 write (nmix,1110) time, te(1), qaxis, rs1x, rs2x, rs3x
      return
  820 write (nmix,1120) time, te(1), qaxis, rs1x, rs2x, rs3x
      return
  830 write (nmix,1130) time, te(1), qaxis, rs1x, rs2x, rs3x
      return
  840 write (nmix,1140) time, te(1), qaxis, rs1x, rs2x, rs3x,
     .   r0x, rmixxx
      return
  850 write (nmix,1110) time, te(1), qaxis, rs1x, rs2x, rs3x,
     .   r0x, rmixxx, widx, taur, s, alpha, gamt, wse, gam
c
c  reset mixing time and island width
c
      if (rmix .eq. 0.0)  return
c
c     profile assignment is done in tport only if time = timmix
c
      timmix = time
      wid    = w0mix
      return
c
 1105 format (/ 5x,'time',5x,'teo',6x,'qo', '  rs1/a', '  rs2/a',
     .           '  rs3/a',3x,'r0/a',1x,'rmix/a', '  wid/a',6x,'taur',
     .          9x,'s',5x,'alpha',6x,'gamt',7x,'wse',7x,'gam' /
     .          6x,'(s)',3x,'(keV)',57x,'(s)',20x,3(5x,'(s-1)'))
 1110 format (f9.4,2f8.4,6f7.3,1p6e10.2)
 1120 format (f9.4,2f8.4,3f7.3,5x,'island>3')
 1130 format (f9.4,2f8.4,3f7.3,5x,'rs values are too close together')
 1140 format (f9.4,2f8.4,5f7.3,5x,'rmix> = rminor or rmix<=rs(island)')
c
      end

      subroutine mixu (dr, hcap, is, j0, jmix, jx, ktab, mx, ntab, r,
     .                 rtab, rxtab, u, up, w, wp, xx, rtaba, utaba)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  This subroutine mixes the quantity u while conserving its integral.
c
      dimension  dr(*), hcap(*), r(*), rtab(*), rxtab(ktab,*),
     .            u(*),   up(*), w(*),   wp(*)
      dimension  xx(ktab,*), rtaba(*), utaba(*)
c
c  multiply u by hcap; also multiply u by w if mx .ne. 2
c
      do 20 j=j0,jmix+1
        u(j) = hcap(j)*u(j)
        if (mx .eq. 2)  go to 20
        u(j) = u(j)*w(j)
   20 continue
c
c  test whether profile of u or u*w is to be made flat
c
      if (mx .le. 0)  go to 500
c
c  generate profile for u consistent with helical flux conservation
c     (when mx>0)
c
c  integrate u vs. rxtab(i,k) before mixing
c
      jb = j0
      rb = rtab(1)
      ub = u(j0) + (u(j0+1)-u(j0))*(rb-r(j0))/dr(j0)
      k1 = 1
      if (is .eq. 3 .and. rtab(1) .ne. 0.0) k1 = 2
      do 110 k=k1,is+1
      call mint1(dr, jb, jmix, ntab, rb, r, rxtab(1,k), ub, u, xx(1,k))
  110 continue
c
c  combine separate integrals to obtain u vs. rtaba after mixing
c
      do 320 i=1,ntab-1
      rtaba(i) = 0.5*(rtab(i)+rtab(i+1))
      xxsum = 0.0
      do 310 k=k1,is+1
      xxsum = xxsum + xx(i,k)
  310 continue
  320 utaba(i) = xxsum/(rtab(i+1)**2-rtab(i)**2)
      rtaba(ntab) = rb
      utaba(ntab) = ub
c
c  interpolate to obtain u vs. r after mixing; use linear interpolation
c     unless rtab(1) = 0 and r is inside of rtaba(1), in which case
c     quadratic interpolation is employed
c
      if (rtab(1) .eq. 0.0)
     .   up(1) = u(jx) + (u(jx+1)-u(jx))*(rxtab(1,is)-r(jx))/dr(jx)
      do 440 j=j0+1,jmix
      if (rtab(1) .ne. 0.0)  go to 410
      if (r(j) .gt. rtaba(1))  go to 410
      up(j) = up(1) + (utaba(1)-up(1))*(r(j)/rtaba(1))**2
      go to 440
  410 continue
      do 420 i=1,ntab-1
      if (r(j) .lt. rtaba(i+1))  go to 430
  420 continue
      i = ntab-1
  430 up(j) = utaba(i) + (utaba(i+1)-utaba(i))*(r(j)-rtaba(i))
     .                   /(rtaba(i+1)-rtaba(i))
  440 continue
c
c  multiply u by w if mx .eq. 2
c
      if (mx .ne. 2)  go to 460
      do 450 j=j0,jmix
      u(j) = u(j)*w(j)
      if (j .eq. j0 .and. rtab(1) .ne. 0.0)  go to 450
      up(j) = up(j)*wp(j)
  450 continue
      u(jmix+1) = u(jmix+1)*w(jmix+1)
  460 continue
c
c  renormalize u to conserve integral of u*w
c
      call mint2(dr, j0, jmix, r, rtab(1), u, uint)
      call mint2(dr, j0, jmix, r, rtab(1), up, upint)
      if (upint .eq. 0.0)  go to 610
      fact = uint/upint
      do 470 j=j0,jmix
      if (j .eq. j0 .and. rtab(1) .ne. 0.0)  go to 470
      up(j) = fact*up(j)
  470 continue
      go to 610
c
c  generate a flat profile for u (if mx = 0) or u*w (if mx<0)
c
  500 call mint2(dr, j0, jmix, r, rtab(1), u, uint)
      do 510 j=j0,jmix
      if (j .eq. j0 .and. rtab(1) .ne. 0.0)  go to 510
      up(j) = hcap(j)*wp(j)
      if (mx .eq. 0)  go to 510
      up(j) = hcap(j)
  510 continue
      call mint2(dr, j0, jmix, r, rtab(1), up, wint)
      uflat = uint/wint
      do 530 j=j0,jmix
      if (j .eq. j0 .and. rtab(1) .ne. 0.0)  go to 530
      up(j) = uflat
      if (mx .eq. 0)  go to 530
      up(j) = hcap(j)*uflat
  530 continue
c
c  divide u by hcap*w; if mx .ne. 0, divide up by hcap*wp
c
  610 continue
      do 620 j=j0,jmix+1
  620 u(j) = u(j)/(hcap(j)*w(j))
      if (mx .eq. 0)  return
      do 630 j=j0,jmix
      if (j .eq. j0 .and. rtab(1) .ne. 0.0)  go to 630
      up(j) = up(j)/(hcap(j)*wp(j))
  630 continue
      return
c
      end

      subroutine mterp (ifor, ja, jb, ntab, pa, pb, psis, ptab,
     .                  ra, rb, r, rxtab)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  This subroutine uses linear interpolation to obtain r values (stored
c     in rxtab) corresponding to psis values (stored in ptab) before
c     mixing.  The r values lie between ra and rb and increase or
c     decrease depending upon the sign of ifor.
c
      dimension psis(*), ptab(*), r(*), rxtab(*)
c
c  save values of r and psis arrays
c
      rja = r(ja)
      rjb = r(jb)
      pja = psis(ja)
      pjb = psis(jb)
c
c  load endpoints of interpolation interval into r and psis arrays
c
      r(ja) = ra
      r(jb) = rb
      psis(ja) = pa
      psis(jb) = pb
c
c  perform interpolation from ra to rb if ifor>0
c
      if (ifor .lt. 0)  go to 105
      do 30 i=1,ntab
      rxtab(i) = ra
      if (ptab(i) .ge. psis(ja))  go to 30
      rxtab(i) = rb
      if (ptab(i) .le. psis(jb))  go to 30
      do 10 j=ja,jb-1
      if (psis(j+1) .lt. ptab(i))  go to 20
   10 continue
   20 rxtab(i) = r(j) + (r(j+1)-r(j))*(ptab(i)-psis(j))
     .                  /(psis(j+1)-psis(j))
   30 continue
      go to 210
c
c  perform interpolation from rb to ra if ifor<0
c
  105 continue
      do 130 i=1,ntab
      rxtab(i) = rb
      if (ptab(i) .ge. psis(jb))  go to 130
      rxtab(i) = ra
      if (ptab(i) .le. psis(ja))  go to 130
      do 110 j=jb-1,ja,-1
      if (psis(j) .lt. ptab(i))  go to 120
  110 continue
  120 rxtab(i) = r(j) + (r(j+1)-r(j))*(ptab(i)-psis(j))
     .                  /(psis(j+1)-psis(j))
  130 continue
c
c  restore original values to r and psis arrays
c
  210 r(ja) = rja
      r(jb) = rjb
      psis(ja) = pja
      psis(jb) = pjb
      return
c
      end

      subroutine olinterp (v1, a1, i1, v2, a2, i2, f, id, val, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- identical to subroutine INTERP as used in original olossa routine.
c --- renamed olinterp for use in ONETWO                     HSJ
c
c this routine evaluates f(v1,v2) by interpolation in a table.
c v1 is found in array a1 and v2 is found in array a2
c i1 is the dimension of a1 and i2 is the dimension of a2
c if i2 is 1 then a one-dimensional function is assumed
c ier = 1 is returned if either v1 is not in a1 or v2 is not in a2
c otherwise ier = 0 is returned.
c interpolation scheme
c id = -2  log on arguments,  log on f
c id = -1  log on arguments,  lin on f
c id = +1  lin on arguments,  lin on f
c id = +2  lin on arguments,  log on f
c the log interpolations assume that the arrays a1,a2, or f are logs.
c
      dimension a1(i1), a2(i2), f(i1,i2)
c
      ier = 0
      if (id .gt. 0)  go to 10
      rv1 = LOG (v1)
      if (i2 .gt. 1)  rv2 = LOG (v2)
      go to 15
   10 rv1 = v1
      if (i2 .gt. 1)  rv2 = v2
   15 call find1(m1,n1,rv1,a1,i1)
      if (m1 .eq. 0 )  go to 35
      if (i2 .gt. 1 )  go to 25
      if (m1 .eq. n1)  go to 20
      val = f(m1,1)+(f(n1,1)-f(m1,1))*(rv1-a1(m1))/(a1(n1)-a1(m1))
      go to 30
   20 val = f(m1,1)
      go to 30
   25 call find1(m2,n2,rv2,a2,i2)
      if (m2 .eq. 0)  go to 35
      val = f(m1,m2)
      if (m1 .eq. n1)  go to 27
      val = val+(f(n1,m2)-f(m1,m2))*(rv1-a1(m1))/(a1(n1)-a1(m1))
   27 if (m2 .eq. n2)  go to 30
      val = val+(f(m1,n2)-f(m1,m2))*(rv2-a2(m2))/(a2(n2)-a2(m2))
   30 if (IABS (id) .eq. 1)  go to 40
      val = exp(val)
      return
   35 val = 0.0
      ier = 1
   40 return
c
      end

      subroutine olossa (fm, fz, tmax, r0, b0, psi0, rx, zx, psix, rmin,
     .                   rmax, rlimin, rlimot, psilin, psilot, rzx,
     .                   psizx, nzx, chi, tmina, nchi)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c routine returns the orbit loss boundary in polar coordinates.
c
c Input Values
c fm is the ion mass in units of proton mass
c fz is the ion charge in units of the proton charge.
c tmax is the largest ion energy considered, eV
c r0 is the birth radius, m.
c b0 is the toroidal field and r0
c psi0 is the psi value at r0
c rx is the radius of the x point, m
c zx is the z position of the x point, m
c psix is the psi value at the x point
c rmin is the smallest radius of the plasma surface,m
c rmax is the largest radius of the plasma surface,m
c rlimin is the inboard radius of the limiter at zsepo
c rlimot is the outboard radius of the limiter at zsepo
c psilin is the psi value at rlimin,zsepo
c psilot is the psi value at rlimot,zsepo
c rzx is a array of radii from just inside rx to rmax at zx,m
c psizx is an array of psi values at the rzx,zx positions
c nzx is the number of such psizx values.
c nchi is the dimension of chi, tmina
c
c Returned Values
c chi is an array of pitch angles from 180 to 90 degrees.
c tmina is an array of kinetic energies (ev) at the loss boundary
c for each chi value.
c
c converted to subroutine from olossa.prog 11-7-92 by Ron Stambaugh
c modified 10-9-94 to add a calculation of orbits that hit the outer limiter.
c Calculation valid only for Er = 0
c ----------------------------------------------------------------------
c
      dimension  ts(1000),ps(1000),pots(1000),rms(1000),psrx(1000)
      dimension  tminc(500), tminm(500), tminl(500)
      dimension  chil(100),tll(100)
      dimension  chiu(100),tlu(100)
      dimension  rzx(nzx),psizx(nzx),chi(nchi),tmina(nchi)
      real*8     mu
      integer    ncrt
      data       ncrt/6/     ! WRITE to unit ncrt instead of using PRINT
c
****  data (rzx(i),i = 1,14)/1.435,1.445,1.455,1.465,1.475,1.485,1.495,
**** .                     1.505,1.515,1.525,1.535,1.700,2.000,2.200/
****  data (psizx(i),i = 1,14)/.0569,.0569,.0569,.0569,.0568,.0567,
**** .           .0565,.0564,.0561,.0559,.0556,.0455,-0.00241,-0.574/
c ----------------------------------------------------------------------
c
      idbg = 0
****  if (ABS (r0-2.11) .lt. 0.02)  idbg = 1
      potx = pot(psix, psix)
      pot0 = pot(psi0, psix)
      ten  = 10.0
      eacc = tmax / 498.0
      eacc = MAX (eacc, ten)
c
c     Following line commented out 4-16-94 to do NBI faster
c
****  eacc  = 10.0
c
      q2m   = (fz*fz*1.6e-19)/(2.*fm*1.67e-27)
      sq2m  = SQRT (q2m)
      iskip = 0
****  if (idbg .eq. 1)  read 991, iskip
c
c test for case of psi0 = psix,
c in which case only the torn open bananas are relevant.
c
      if (psi0 .le. psix)  go to 290
c
c     loop over pitch angles
c
      do 500 i=1,nchi
      chi(i) = 180.0 - (90.0 / (nchi-1)) * (i-1)
      schi   = SIN (chi(i)*3.14159/180.0)
      cchi   = COS (chi(i)*3.14159/180.0)
      schi2  = schi**2
      cchi2  = cchi**2
      if (idbg  .eq. 1)
     .write  (ncrt, 900)  i, chi(i), schi, cchi
  900 format (1x, i3, 3e10.3)
      if (iskip .eq. i)
     .read 991,  iskip
  991 format (i3)
c
c look at the counter loss orbits
c
      tminc(i) = tmax
      if (schi2 .ge. rx/r0)  go to 80
      tmin = -pot0/(1.-(r0/rx)*schi2)
      if (tmin .lt. 0.0)  tmin = 0.0
      t = tmax
    5 mu = t*schi2/b0
      e = t+pot0
      fff = 1.-mu*b0*r0/(e*rx)
      if (fff .le. 0.0)  go to 80
      enew = q2m * (psix-psi0)**2 / (rx * SQRT (1.0-mu*b0*r0/(e*rx))
     .       -r0 * SQRT (1.0-mu*b0/e-pot0/e))**2
      tnew = enew-pot0
      de = enew-e
      if (iskip .eq. i)  write (ncrt, 992)  tmin, tmax, t, e, enew, de
  992 format (1x,6e10.3)
      if (ABS (de) .lt. eacc)  go to 60
      if (   tnew  .gt. tmin)  go to 50
      tnew = tmax
      go to 60
   50 t = tnew
      go to 5
   60 tminc(i) = tnew
c
c here find mirror orbits that hit the outer limiter
c coding as of 10-9-94 valid only for zero potential
c
   80 tminl(i) = tmax
      rm = r0*schi2
      if (rm .le. rlimin)  go to 100
      if (rm .gt. r0    )  go to 100
      tminl(i) = q2m*(psi0-psilot)**2/
     .          (rlimot * SQRT (1.0-rm/rlimot)+r0 * SQRT (1.0-rm/r0))**2
c
c look at mirror loss orbits through x-point
c
  100 tminm(i) = tmax
      if (iskip .eq. i)
     .write  (ncrt, 993) tnew
  993 format (' at label 100 waiting on iskip', e10.3)
      if (iskip .eq. i)  read 991, iskip
      tminm(i) = tmax
      tmx      = tmax
      if (chi(i) .lt. 91 .and. chi(i) .gt. 89)  go to 500
      if (cchi2  .le. 0.0)  go to 101
      tmx = q2m*(psi0-psix)**2/r0/r0/cchi2
      if (tmx .le. -pot0+1.0)  go to 200
  101 if (tmx .gt.  tmax    )  tmx = tmax
c
c  calculate the mirror radius for this largest tmx
c
      es     = tmx+pot0
      mu     = (tmx/b0)*schi2
      psitmx = psi0 - r0 * SQRT (tmx-mu*b0) / sq2m
      pottmx = pot(psitmx,psix)
      rmtmx  = tmx*r0*schi2/(tmx+pot0-pottmx)
c
c  use different logic depending on the sign of the potential
c
      if (pot0)  131, 131, 135
  131 tlo = MAX (-pot0+1.0, ten)
c
c  calculate the largest value of rmirror
c
      es    = tlo+pot0
      mu    = (tlo/b0)*schi2
      psim  = psi0 - r0 * SQRT (tlo-mu*b0) / sq2m
      potm  = pot(psim,psix)
      rmmax = tlo*r0*schi2/(tlo+pot0-potm)
      if (idbg .eq. 1)  write (ncrt, 962)  i, pot0, tlo, psim, potm,
     .            rmmax, tmx, psitmx, pottmx, rmtmx, rx, rmin, rmax
c
c  if rmtmx is greater than rx, then no loss orbits
c
      if (rmtmx .gt. rx)  go to 200
c
c  if rmmax is less than rmin, then no loss orbits here.
c
      if (rmmax .lt. rmin)  go to 200
      go to 139
  135 if (rmtmx .lt. rmin)  go to 200
      tlo   = 1.0
  136 tlo   = tlo+eacc
      if (tlo .gt. tmx-eacc)  go to 200
      es    = tlo+pot0
      mu    = (tlo/b0)*schi2
      psim  = psi0 - r0 * SQRT (tlo-mu*b0) / sq2m
      potm  = pot(psim,psix)
      rmtlo = tlo*r0*schi2/(tlo+pot0-potm)
      if (rmtlo .le. rmin)  go to 136
      if (rmtlo .ge. rx  )  go to 200
      if (idbg .eq. 1)  write (ncrt, 962)  i, pot0, tlo, psim, potm,
     .                            rmtlo, tmx, psitmx, pottmx, rmtmx
  962 format (' mirror', i3, (1x,5e12.5))
c
  139 nt = (tmx-tlo)/eacc
      nt = nt+1
      if ( idbg .eq. 1)  write (ncrt, 962)  i, pot0, tlo, psim, potm,
     .             rmmax, tmx, psitmx, pottmx, rmtmx, rx, rmin, rmax
      if (iskip .eq. i)  write (ncrt, 963)  i, tmx, tlo, nt
  963 format (' i,tmx,tlo,nt', i3, 2e10.3, i3)
c
      do 110 j=1,nt
      ts(j)   = tlo+eacc*(j-1)
      es      = ts(j)+pot0
      mu      = (ts(j)/b0)*schi2
      ps(j)   = psi0-r0* SQRT (ts(j)-mu*b0)/sq2m
      pots(j) = pot(ps(j),psix)
      rms(j)  = ts(j)*r0*schi2/(ts(j)+pot0-pots(j))
      psrx(j) = 0.0
      if (rms(j) .gt. rx .or. rms(j) .lt. rmin)  go to 109
      psim    = ps(j)
      rm      = rms(j)
      potm    = pots(j)
c
****  if (iskip .eq. i)  write (ncrt, 982) i, tmx, rm, pot0, tmx1, tmx2
* 982 format (' inside ', i3, 5e10.3)
****  if (iskip .eq. i)  write (ncrt, 995) tmx, tlo, dt, t, mu, psim, rm
* 995 format (' at 103 ', 7e10.3)
****  if (iskip .eq. i)  write (ncrt, 998) tmx, tlo, dt, t, mu, psim, rm
* 998 format (' above 110 ', 7e10.3)
c
c     get psi at rx
c
  115 psiacc = (psi0 - psix) / 100.0
c
c  try for a solution with potrx = 0
c
      potrx = 0.0
      vpar  = ts(j)+pot0-mu*b0*r0/rx-potrx
      if (vpar .le. 0.0 )  go to 120
      psi   = psim-rx* SQRT (vpar)/sq2m
      if ( psi .gt. psix)  go to 120
c
c  found a lost orbit
c
      psrx(j) = psi
      go to 109
c
c  look for confined orbit, psix < psi < psim
c  get intermediate potential point
c
  120 psimid = (psim + psix) / 2.0
      potmid = pot (psimid, psix)
      call quadft (psix,potx,psimid,potmid,psim,potm,cd,ce,cf,ieror)
      if (ieror .eq. 1)  go to 109
c
c  form the coefficients of a quadratic eqn for psi
c
      aq   = q2m/rx/rx+cf
      bq   = -2.0 * psim * q2m / rx / rx + ce
      cq   = psim*psim*q2m/rx/rx-es+mu*b0*r0/rx+cd
      disc = bq*bq-4.*aq*cq
      if (disc)  109, 121, 123
  121 psi  = -bq / 2.0 / aq
      if (psi .lt. psix .or. psi .gt. psim)  go to 109
      psrx(j) = psi
      go to 110
  123 psi1 = (-bq + SQRT (disc)) / 2.0 / aq
      psi2 = (-bq - SQRT (disc)) / 2.0 / aq
      if (psi1 .lt. psix .or. psi1 .gt. psim)  go to 127
      pot1 = pot(psi1,psix)
c
c  check vparallel
c
      vpar    = es-mu*b0*r0/rx-pot1
      if (vpar .le. 0.0)  go to 127
      psrx(j) = psi1
c
c  check if psi2 is also a solution
c
      if (psi2 .lt. psix .or. psi2 .gt. psim)  go to 109
c
c  psi2 is also a solution
c
      pot2 = pot(psi2,psix)
      vpar = es-mu*b0*r0/rx-pot2
      if (vpar .le. 0.0)  go to 109
      if (idbg .eq. 1  )  write (ncrt, 971)  psi1, psi2
  971 format (' two solutions ', 2e12.5)
c
c  if two solns, take one closest to separatrix?
c
      if (psi2 .lt. psrx(j))  psrx(j) = psi2
      go to 109
c
c  check if psi2 is the only solution
c
  127 if (psi2 .lt. psix .or. psi2 .gt. psim)  go to 109
      psrx(j) = psi2
  109 if (idbg .eq. 1)
     .write  (ncrt, 961) j, ts(j), ps(j), pots(j), rms(j), psrx(j), vpar
  961 format (i4, 6e10.3)
  110 continue
c
      do 150 j=1,nt
      jsv = j
  150 if (psrx(j) .gt. 0.0 .and. psrx(j) .le. psix)  go to 160
      go to 200
  160 tminm(i) = ts(jsv)
  200 if (iskip .eq. i)  write (ncrt, 996)  tmx, tlo,     tminm(i)
* 200 if (iskip .eq. i)  write (ncrt, 996)  tmx, tlo, dt, tminm(i) ! dt?
  996 format (' below 200', 7e10.3)
  500 if (iskip .eq. i)  iskip = iskip + 1
      go to 299
c
c for points outside the separatrix,
c check mirror orbits that hit the outer limiter
c
  290 do 295 i=1,nchi
      chi(i) = 180.0 - (90.0 / (nchi-1)) * (i-1)
      schi   = SIN (chi(i)*3.14159/180.0)
      schi2  = schi**2
c
c here find mirror orbits that hit the outer limiter
c coding as of 10-9-94 valid only for zero potential
c but the potential should be zero for points outside the separatrix
c
      tminl(i) = tmax
      rm = r0*schi2
      if (rm .le. rlimin)  go to 295
      if (rm .gt. r0    )  go to 295
      tminl(i) = q2m*(psi0-psilot)**2/
     .          (rlimot* SQRT (1.-rm/rlimot)+r0 * SQRT (1.-rm/r0))**2
  295 continue
c
c now do the torn open bananas
c
  299 jmin = 0
      jmax = nzx
      do 300 j=1,nzx
      if (rzx(j) .lt. r0)  go to 311
      go to 300
  311 jmax = j
      e = (q2m*(psi0-psizx(j))**2/r0/r0+pot0)/(1.-rzx(j)/r0)
      if (e .gt. 0 .and. rzx(j) .ge. rx)  go to 312
      jmin = j
      go to 300
  312 mu = e*rzx(j)/b0/r0
      t = e-pot0
      schi2 = b0*mu/t
      if (schi2 .le. 1.0)  go to 313
      jmin = j
      go to 300
  313 schi    = SQRT (schi2)
      chil(j) = 180.0 - 180.0 * ASIN (schi) / 3.14159
  300 tll(j)  = t
      jmin = jmin+1
      ntorn = jmax-jmin+1
c
c valid values should lie between jmin and jmax
c
      if (idbg .eq. 1)  write (ncrt, 314)  jmin, jmax, ntorn
  314 format (' torn orbit search, jmin,jmax,ntorn = ',3i5)
      if (jmax .ge. jmin)  go to 315
      go to 316
  315 do 302 j=jmin,jmax
      chiu(ntorn-j+jmin) = chil(j)
  302 tlu(ntorn-j+jmin) = tll(j)
      if (idbg .eq. 1)  write (ncrt, 999)  (chiu(i), tlu(i), i=1,ntorn)
  999 format (6e10.3)
c
c look for psi0 outside separatrix
c
  316 if (psi0 .gt. psix)  go to 360
      itob = nchi
      do 350 i=1,nchi
      chi(i) = 180.0 - (90.0 / (nchi-1)) * (i-1)
      if (ntorn .gt. 0)  go to 345
      tmina(i) = 0.0
      go to 350
  345 call olinterp (chi(i), chiu, ntorn, chi(i), chiu, 1, tlu,
     .                                                  1, val, ier)
c
      if (ier .ne. 1) then
        tmina(i) = val
        if (tmina(i) .gt. tmax)  tmina(i) = tmax
        itob = i
        if (idbg .eq. 1)  write (ncrt, 772) val
  772   format (' setting itob = i',e12.5)
      else if (i .lt. itob) then
        tmina(i) = 0.0
        if (idbg .eq. 1)  write (ncrt, 773)  i, chi(i)
  773   format ('took first ELSEIF', i5, e10.3)
      else if (i .gt. itob) then
        tmina(i) = tmax
      end if
c
      if (idbg .eq. 1)  write (ncrt, 349)  i, itob, ier, tmina(i), val
  349 format (' i,itob,ier,tmina,val = '3i5,2e12.5)
  350 continue
      do 355 i=1,nchi
  355 tmina(i) = MIN (tmina(i), tminl(i))
      return
c
c  fill in missing mirrors.
c
  360 do 388 i=1,nchi-2
      if (tminm(i) .eq. tmax)  go to 388
c
c  find next mirror
c
      do 384 j=i+1,nchi
        if (tminm(j) .eq. tmax)  go to 384
        go to 385
  384 continue
      go to 388
  385 if (j-i .eq. 1)  go to 388
      avg = (tminm(i)+tminm(j))/2.0
      do 386 k=i+1,j-1
  386 tminm(k) = avg
  388 continue
c
      imsv = 1
      itob = 1
      icsv = 1
c
      do i=1,nchi
        if (ntorn .gt. 0)  go to 371
        ier = 1
        go to 372
  371   call olinterp(chi(i),chiu,ntorn,chi(i),chiu,1,tlu,1,val,ier)
  372   if (ier .eq. 1)  val = tmax
c
c  this test saves the last index where a mirror lost orbit was.
c  these indices help fill a gap between the end of the mirror loss
c  angle range and the start of the torn open banana angle range.
c
        if (tminm(i) .lt. tmax)  imsv = i
        if (tminc(i) .lt. tmax)  icsv = i
        if (     val .eq. tmax)  go to 397
        if (    itob .gt. 1   )  go to 397
c
c  here is the first torn open banana orbit
c
        itob = i
c
c  now olinterpolate all values between imsv and itob
c  unless there were no mirror loss orbits, then use icsv.
c
        if (imsv .eq. 1)  go to 394
        j1 = imsv+1
        j2 = itob-1
        if (idbg .eq. 1)
     .  write  (ncrt, 391) j1, j2, icsv, imsv, itob, val, tminm(imsv)
  391   format (5i4,2e12.5)
        if (j1 .gt. j2)  go to 397
c
        do j=j1,j2
          vinter = tminm(imsv) +
     .      (val-tminm(imsv))*(chi(j)-chi(imsv))/(chi(itob)-chi(imsv))
          tmina(j) = MIN (tminc(j), vinter)
          if (idbg .eq. 1)
     .    write  (ncrt, 989)  tmina(j), vinter, tminc(j)
  989     format (1x, 6(1x, e10.3))
        end do
c
        go to 397
  394   if (icsv .eq. 1)  go to 397
        j1 = icsv+1
        j2 = itob-1
        if (idbg .eq. 1)
     .  write (ncrt, 391) j1, j2, icsv, imsv, itob, val, tminm(imsv)
        if (  j1 .gt. j2)  go to 397
c
        do j=j1,j2
          vinter = tminc(icsv) +
     .      (val-tminc(icsv))*(chi(j)-chi(icsv))/(chi(itob)-chi(icsv))
          tmina(j) = MIN (tminc(j), vinter)
          if (idbg .eq. 1)  write (ncrt, 989) tmina(j), vinter, tminc(j)
        end do
c
  397   tmina(i) = MIN (tminc(i), tminm(i), val)
        if (itob .gt. 1 .and. i .ge. itob)  tmina(i) = val
        tmina(i) = MIN (tmina(i), tminl(i))
        if (tmina(i) .gt. tmax)  tmina(i) = tmax
  399   if (idbg .eq. 1)  write (ncrt, 989) chi(i), tminc(i), tminm(i),
     .                                      val, tminl(i), tmina(i)
      end do
c
      if (idbg .ne. 1)  go to 480
      do 481 i=1,nchi
  481 write (ncrt, 489)  chi(i), tmina(i)
c
c  GLITCH BUSTER
c
  480 do 488 i=1,nchi-2
      if (tmina(i) .gt. tmax)  tmina(i) = tmax
      if (tmina(i) .eq. tmax)  go to 487
      do 484 j=i+1,nchi
      if (tmina(j) .ge. tmax)  go to 484
      go to 485
  484 continue
      go to 487
  485 if (j-i .eq. 1)  go to 487
      avg = (tmina(i)+tmina(j))/2.0
      do 486 k=i+1,j-1
  486 tmina(k) = avg
  487 if (idbg .eq. 1)  write (ncrt, 489)  chi(i), tmina(i)
  489 format (1x,2e10.3)
  488 continue
c
      return
c
      end

      real*8 function pot (psi, psix)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
      parameter (ner = 8)
c
      common /potential/ psipt(200), potpt(200), npots
c
      pot = 0.0
      if (npots .le. 1)  return
c
      if (psi .le. psix) then
        pot = 0.0
      else if (psi .ge. psipt(npots)) then
        pot = potpt(npots)
      else if (psi .le. psipt(1)    ) then
        pot = potpt(1)
      else
        call olinterp(psi,psipt,npots,psi,psipt,1,potpt,1,val,ier)
        pot = val
      end if
      return
c
      end

      subroutine printang (nunit, timet, t, jprt)
c
      USE param
      USE ions
      USE solcon
      USE soln
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE tordlrot
      USE constnts
      USE flx,  ONLY : flux
      USE bd_condtn,only : bctime,ub,fluxb,
     .    ub_save,ub_rho_edge,bctime_zone

      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- print out information related to toroidal rotation
c
c      include 'param.i'
c      include 'bcon.i'
c      include 'flx.i'
c      include 'geom.i'
c      include 'ions.i'
c      include 'machin.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'constnts.i'
c      include 'solcon.i'
c     include 'soln.i'
      include 'storage.i'
c      include 'tordlrot.i'
c
      data mtmfile /0/
c
c ----------------------------------------------------------------------
c --- create a file for momentum output
c ----------------------------------------------------------------------
c
      if (mtmfile .eq. 0 .and. momtm_file .eq. 1) then
        mtmfile =  1
        call getioun(ioc,41)
        open (unit = ioc, file = 'momtmout', status = 'UNKNOWN')
      end if
c
c ----------------------------------------------------------------------
c --- get the local angular momentum conf. time. (theta weighting of
c --- torque terms is neglected)
c --- (was not calculated if psourc was not called)
c --- calculate various other quantites related to tor. rotation
c ----------------------------------------------------------------------
c
      onemt  = 1.0 - theta
      xmassp = 1.6726231e-24
        do j=1,nj
          xdum(j) = onemt * usave(nk,j) + theta * u(nk,j)
          ydum(j) = 4.0 * pi**2 * rmajor * hcap(j) / volume
          avg0    = 0.0
          avg1    = 0.0
          do k=1,nion
            avg0  = avg0 + usave(k,j) * usave(nk,j) * atw(k)
            avg1  = avg1 + u    (k,j) * u    (nk,j) * atw(k)
          end do
          avg0         =  avg0 * r2capi0(j) * xmassp
          avg1         =  avg1 * r2capi (j) * xmassp
          avgangmt(j)  = onemt * avg0 + theta * avg1
          rotenergy(j) = 0.5 * (onemt * avg0 * usave(nk,j)
     .                 + theta * avg1 * u(nk,j))
          avgangdt(j)  = 0.0
          if (dtt .ne. 0.0)
     .    avgangdt(j)  = (avg1 - avg0) / dtt
          crud         = storqueb(j) - avgangdt(j) ! beam torque only ??? HSJ
          tauang(j)    = 0.0
          if (crud .ne. 0.0)  tauang(j) = avgangmt(j) / crud
          if (   j .lt.  nj)    zdum(j) = flux(nk,j)
        end do
      zdum(nj) = fluxb(nk)
c
      call trapv (r, xdum, ydum, nj, avgomega)
      call trapv (r, dangrot, ydum, nj, davgomdt)
      call trapv (r, avgangmt, ydum, nj, avgmtm)
      call trapv (r, avgangdt, ydum, nj, davgmtdt)
      call trapv (r, vionz, ydum, nj, avgvionz)
      call trapv (r, zdum, ydum, nj, avgflux)
      call trapv (r, xkangrot, ydum, nj, avgxkang)
      call trapv (r, tauang, ydum, nj, avgtauan)
      call trapv (r, amtinrta, ydum, nj, avgamtin)
      call trapv (r, rotenergy, ydum, nj, avgrote)
c
c --- volume average sources
c
      call trapv (r, sprbeame, ydum, nj, avgsprbe)
      call trapv (r, sprbeami, ydum, nj, avgsprbi)
      call trapv (r, ssprcxl, ydum, nj, avgsscx)
      call trapv (r, sprcxre, ydum, nj, avgsprcx)
      call trapv (r, spreimpt, ydum, nj, avgsprei)
      call trapv (r, sprcx, ydum, nj, avgspcx)
      call trapv (r, spr2d, ydum, nj, avgspr2d)
      call trapv (r, storque, ydum, nj, avgstorq)
c
      avrotkev = kevperg * avgrote
      avrotjou = avgrote * 1.0e-07 * volume
      angmtotl = volume * avgmtm * 1.0e-07
      crud     = beamtorq - davgmtdt * 1.0e-07
      tauangtl = 0.0
      if (crud .ne. 0.0)
     .tauangtl = angmtotl / crud
c
      call header (nunit, timet, t)
      write  (nunit, 1000)
 1000 format (2x,40(1h-), 'TOROIDAL ROTATION RESULTS',53(1h-) //
     .     88x,'total',9x,'local'                              /
     .     3x,'j',5x,'r',7x,'omega',4x,'d(omega)/dt',
     .     4x,'angmtm',3x,'d(angmtm)/dt',3x,'vionz',7x,'flux',
     .     6x,'ang. momtm', '  ang. momtm', '  momt inrta'     /
     .     14x,'(ang speed)',16x,'density',14x,'(ion speed)',
     .     13x,'diffusivity', '  conf. time',3x,'density'      /
     .     9x,'cm',6x,'1./sec',5x,'1./sec**2', '  g/cm*sec',3x,
     .     'g/cm*sec**2',3x,'cm/sec',6x,'g/sec**2', '  cm**2/sec',
     .     6x,'sec',8x,'g/cm')
c


      do 100 j=1,nj
        j1prt = ((j-1) / jprt) * jprt
        if ((j1prt .ne. j-1) .and. (j .ne. nj))  go to 100
        write  (nunit, 1010) j,r(j),xdum(j),dangrot(j),avgangmt(j),
     .                    avgangdt(j),vionz(j),zdum(j),xkangrot(j),
     .                                       tauang(j),amtinrta(j)
 1010   format (1x,i4,2x,f6.1,9(2x,1pe10.2))
  100 continue
c
      write  (nunit, 1011) avgomega,davgomdt,avgmtm,davgmtdt,
     .                     avgvionz,avgflux,avgxkang,avgtauan,avgamtin
 1011 format (/ '  volume avg.',9(2x,1pe10.2))
      write  (nunit, 1015)  nj-1, nj
 1015 format (// '  diffusivity and flux are on the half grid in rho' /
     .           '  for j = 1 to',i5,' extrapolated values are',
     .          '   given for j= ',i5)
      write  (nunit, 1020)  angmtotl, tauangtl, totinrta
 1020 format ('  stored ang. mtm., kg*m**2/sec:',1pe12.3 /
     .      '  global ang. mtm. conf. time,sec:',1pe12.3 /
     .       '  total momt. of inertia,kg*m**2:',1pe12.3)
c
c ----------------------------------------------------------------------
c --- output to file momtmout
c ----------------------------------------------------------------------
c
      if(momtm_file .eq. 1)then
      write  (ioc, 2000)  timet, nj, volume, rmajor
 2000 format (2x, 1pe14.6, 2x, i6, 2x, 1pe14.6, 2x, 1pe14.6)
 2020 format (6(1pe16.6, 2x))
      write  (ioc, 2020)  (r       (j), j=1,nj)
      write  (ioc, 2020)  (xdum    (j), j=1,nj)
      write  (ioc, 2020)  (dangrot (j), j=1,nj)
      write  (ioc, 2020)  (avgangmt(j), j=1,nj)
      write  (ioc, 2020)  (avgangdt(j), j=1,nj)
      write  (ioc, 2020)  (vionz   (j), j=1,nj)
      write  (ioc, 2020)  (zdum    (j), j=1,nj)
      write  (ioc, 2020)  (xkangrot(j), j=1,nj)
      write  (ioc, 2020)  (tauang  (j), j=1,nj)
      write  (ioc, 2020)  (amtinrta(j), j=1,nj)
      write  (ioc, 2020)  (fluxangv(j), j=1,nj)
      write  (ioc, 2020)  (fluxangc(j), j=1,nj)
      write  (ioc, 2020)  (dkapomeg(j), j=1,nj)
      endif
c
c --- second page of output
c
      call header (nunit,timet,t)
      write  (nunit, 1100)
 1100 format (2x,40(1h-),'TOROIDAL ROTATION SOURCES',42(1h-),
     .                   'rot energy' /
     .        2x,47(1h-),'g/(cm*sec**2)',47(1h-),'ergs/cm**3')
      write  (nunit, 1110)
 1110 format (3x,'j',6x,'r',6x,'sprbeame',4x,'sprbeami',
     .      4x,'ssprcxl',5x,'sprcxre',5x,'spreimpt',
     .       5x,'sprcx',7x,'spr2d',7x,'stotal')
      do 200 j=1,nj
        j1prt = ((j-1)/jprt)*jprt
        if ((j1prt .ne. j-1) .and. (j .ne. nj))  go to 200
        write (nunit, 1010) j,r(j),sprbeame(j),sprbeami(j),ssprcxl(j),
     .                      sprcxre(j),spreimpt(j),sprcx(j),spr2d(j),
     .                      storque(j),rotenergy(j)
  200 continue
      write (nunit, 1011) avgsprbe,avgsprbi,avgsscx,avgsprcx,
     .                    avgsprei,avgspcx,avgspr2d,avgstorq,avgrote
c
c --- print explanation of terms
c
      write  (nunit, 1500)
 1500 format (// 6x,'sprbeame and sprbeami represent input of angular' /
     .      '  momentum to ion distribution due to coulomb collisions' /
     .      '  of fast ions with thermal ions and electrons')
      write  (nunit, 1510)
 1510 format (6x,'ssprcxl is source due to secondary charge exchange' /
     .    '  (fast ion plus thermal neutral forms thermal ion '       /
     .    '  at neutral temp. plus fast neutral)')
      write  (nunit, 1520)
 1520 format (6x,'sprcxre is sink due to recombination of thermal'    /
     .    '  ions and charge exchange of thermal ions with fast'      /
     .    '  neutrals (impurities are neglected here)')
      write  (nunit, 1530)
 1530 format (6x,'spreimpt is source due to electron impact ' /
     .      '  ionization of thermal neutrals (ion impact is' /
     .      '  neglected since the neutral transport model'   /
     .      '  is valid only for energies where ion impact'   /
     .      '  is negligible)')
      write  (nunit, 1540)
 1540 format (6x, 'sprcx is drag due to charge exchange of' /
     .          '  thermal neutrals with thermal ions')
      write  (nunit, 1550)
 1550 format (6x, 'spr2d is source due to two-dimensional effects' /
     .          '  ie changing geometry and moment of inertia')
      write  (nunit, 1560)
 1560 format (6x, 'rot energy is the kinetic energy of rotation')
      if (nbeamtcx .eq. 0)  write (nunit, 1570)
 1570 format ('  beam and total torque neglect cx losses')
      if (nbeamtcx .eq. 1)  write (nunit, 1580)
 1580 format ('  beam and total torque account for cx losses')
c
c ----------------------------------------------------------------------
c --- output to file momtmout
c ----------------------------------------------------------------------
c
      if(momtm_file .eq. 1)then
      write  (ioc, 2020) (sprbeame (j), j=1,nj)
      write  (ioc, 2020) (sprbeami (j), j=1,nj)
      write  (ioc, 2020) (ssprcxl  (j), j=1,nj)
      write  (ioc, 2020) (sprcxre  (j), j=1,nj)
      write  (ioc, 2020) (spreimpt (j), j=1,nj)
      write  (ioc, 2020) (sprcx    (j), j=1,nj)
      write  (ioc, 2020) (spr2d    (j), j=1,nj)
      write  (ioc, 2020) (storque  (j), j=1,nj)
      write  (ioc, 2020) (rotenergy(j), j=1,nj)
      write  (ioc, 2020) (hcap     (j), j=1,nj)
      write  (ioc, 2020) (rcap     (j), j=1,nj)
c rcapi not written out
      write  (ioc, 2020) (r2capi   (j), j=1,nj)
      write  (ioc, '(///)')
      endif
      return
c
      end

      subroutine print_ifs (r, nj, jprt, nout, t,ifsflag)
c
c
c ----------------------------------------------------------------------
c --- print out Dorland-Kotchenreuther results,
c     calculated by subroutine IP_CHI2
c --- the values to be printed were calculated at the time stored in
c --- time_ifs
c ----------------------------------------------------------------------
c
      USE param
      USE mhdpar 
      USE machin
      USE tcoef
      USE rhog
      USE ifs
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'mhdpar.i'
c      include 'machin.i'
c      include 'ifs.i'
c      include 'rhog.i'
      include 'rebut.i'
      include 'storage.i'
c      include 'tcoef.i'
c
      dimension r(nj)
c
      do j=1,nj
          if (j .eq. 1) then
            call header (nout, time_ifs, t)
            write (nout, 10) time_ifs
          end if
c
          j1prt = ((j-1) / jprt) * jprt
c
          if (j1prt .ne. j-1 .and. j .ne. nj)  go to 50
          write (nout, 20) j, r(j), psir(j), rmajorvec(j),
     .                     rhod_ifs(j),
     .                     d_ifs(j),xdchitot(j),chi_e_ifs(j),
     .                     xchietot(j), chi_i_ifs(j), xchiitot(j)
c
   50     continue
        end do
c
      do j=1,nj
          if (j .eq. 1) then
            call header (nout, time_ifs, t)
            write (nout, 30) time_ifs
          end if
c
          j1prt = ((j-1) / jprt) * jprt
c
          if (j1prt .ne. j-1 .and. j .ne. nj)  go to 60
          write (nout, 20) j, r(j), psir(j), rmajorvec(j),
     .                     rhod_ifs(j),xke_ifs(j),
     .                     xketot(j), xki_ifs(j),
     .                     xkitot(j),xkang_ifs(j)
c
   60     continue
        end do
        write (nout, 25) rhod_max_ifs
c
      do j=1,nj
          if (j .eq. 1) then
            call header (nout, time_ifs, t)
            write (nout, 40) time_ifs
          end if
c
          j1prt = ((j-1) / jprt) * jprt
c
          if (j1prt .ne. j-1 .and. j .ne. nj)  go to 70
          write (nout, 20) j,rhod_ifs(j)*rhod_max_ifs,
     .                     rhod_ifs(j),rlt_ifs(j),
     .                     rln_ifs(j), rlne_ifs(j),
     .                     RLTcrit_ifs(j),RLTcritz_ifs(j),
     .                     g_ifs(j),shat_ifs(j)
c
   70     continue
        end do
        write (nout, 25) rhod_max_ifs
c
   10 format (1x, 48(1h-), 1x,
     .       'IFS Model: Diffusivities ',
     .        60(1h-)                                       /
     .        41x, 'of confinement at time ',1pe12.4,' sec' /
     .        3x,'j',6x,'rho',6x,'psi',10x,'rmaj',4x,
     .        'rhod_IFSn',4x,
     .       ' d  IFS  ', '  d  tot ',4x,'chie IFS',2x,'chie tot',
     .        2x,'chii IFS',4x,'chii tot ',3x,'        '    /
     .        7x,'cm',7x,'kgauss-cm**2',7x,'cm',3x,'       ', 3x,
     .        'cm**2/sec',
     .        3x, 'cm**2/sec', 2x, 'cm**2/sec', '  cm**2/sec',
     .     '  cm**2/sec', '  cm**2/ sec' , '  cm**2/ sec '  /
     .        2x, 126(1h=))
c
   20 format (2x, i3, 10(1pe11.3))
   25 format (2x,'maximum value of rhod = ',1pe12.6)
c
   30 format (1x, 38(1h-), 1x,
     .       'IFS Model: Conductivites ',
     .        30(1h-)                                       /
     .        31x, 'of confinement at time ',1pe12.4,' sec' /
     .        3x,'j',3x,'rho',9x,'psi',8x,'rmaj',7x,'rhod_IFSn',
     .        2x,'ke  IFS',6x,'ke tot',
     .        5x,'ki IFS',5x,'ki tot ',5x,'kang tot',       /
     .        7x,'cm',7x,'kgauss-cm**2',4x,'cm',6x,'cm',6x,
     .        '1/(cm sec)',
     .        x, '1/(cm sec)', 2x, '1/(cm sec', '  1/(cm sec)',
     .        'g??' /
     .        2x, 126(1h=))
c
   40 format (1x, 48(1h-), 1x,
     .       'IFS Model: PARAMETERS ',
     .        60(1h-)                         /
     .        48x, ' at time ',1pe12.4,' sec' /
     .        3x,'j',3x,'rhod_IFS',3x,'rhod_IFSn'
     .        8x,'rlt',6x,'rln',
     .        6x,'rlne',5x,'RLTcrit ',4x,'RLTcritz',2x,'crit ratio',
     .        2x,'shear'                      /
     .        9x, 'cm', 7x, '  ', 4x, '   ', 6x, '  ', 6x,
     .        '   ',
     .        6x, '  ', 4x, '   ', 4x, '    ' /
     .        2x, 126(1h=))
      return
c
      end

      subroutine print_rebut (r, nj, jprt, nout, t,ifsflag)
c
      USE param
      USE machin
      USE tcoef
      USE rhog
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- print out Rebut-Lallia-Watkins results, calculated by subroutine REBUTLAL.
c --- the values to be printed were calculated at the time stored in timerebut
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'machin.i'
      include 'rebut.i'
c      include 'rhog.i'
      include 'storage.i'
c      include 'tcoef.i'
c
      dimension    r(nj)
      character*40 timsg
c
      if (tirlw .eq. 0) then
        timsg =      'used REAL ti for RLW calculations'
      else
        timsg = 'used EFFECTIVE ti for RLW calculations'
      end if
c
      do j=1,nj
          if (j .eq. 1) then
            call header (nout,timerebut,t)
            if (ifsflag .eq. 0) then
              write (nout, 10)  timerebut, timsg
            else
              write (nout, 11)  timerebut, timsg
            end if
          end if
c
          j1prt = ((j-1) / jprt) * jprt
c
          if (j1prt .ne. j-1 .and. j .ne. nj)  go to 50
          write (nout, 20) j, r(j), psir(j), rmajorvec(j),
     .                     crit_grad(j), grad_te(j), xchierl(j),
     .                     xchiirl(j),xeffrl(j),xeffsf(j)
   50     continue
        end do
c
   10 format (1x, 38(1h-), 1x,
     .       'Rebut-Lallia critical temperature gradient model ',
     .        30(1h-)                                        /
     .        51x, 'of confinement at time ',1pe12.4,' sec'  /
     .        45x, '(gradients calculated in rho space) '    /
     .        45x, a                                         /
     .        3x,'j',3x,'rho',6x,'psi',10x,'rmaj',7x,
     .       'crit-grad', '  gradte',5x,'xchierl',4x,'xchiirl',
     .        4x,'xeffrl',5x,'xeffsf'                        /
     .        7x,'cm',7x,'kgauss-cm**2',1x,'cm',9x,'keV/cm',
     .        5x, 'keV/cm', 5x, 'cm**2/sec', '  cm**2/sec',
     .     '  cm**2/sec', '  cm**2/sec'                      /
     .        2x, 126(1h=))
c
   11 format (1x, 38(1h-), 1x,
     .       'Rebut-Lallia critical temperature gradient model ',
     .        30(1h-)                                        /
     .        51x, 'of confinement at time ',1pe12.4,' sec'  /
     .        45x,'(gradients calculated in rho space) '     /
     .        45x,'(FLOW SHEAR SUPPRESSION IS IN EFFECT) '   /
     .        45x, a                                         /
     .        3x,'j',3x,'rho',6x,'psi',10x,'rmaj',7x,
     .       'crit-grad', '  gradte',5x,'xchierl',4x,'xchiirl',
     .        4x,'xeffrl',5x,'xeffsf'                        /
     .        7x, 'cm', 7x, 'kgauss-cm**2', ' cm', 9x, 'keV/cm',
     .        5x, 'keV/cm', 5x, 'cm**2/sec', '  cm**2/sec',
     .     '  cm**2/sec', '  cm**2/sec'                      /
     .        2x,126(1h=))
c
   20 format (2x, i3, 10(1pe11.3))
c
      return
c
      end

      subroutine print_rgcmodel (nout, t)
c
c                           SUBPROGRAM DESCRIPTION:
c                           ----------------------
c
c   Print to logical unit nout results of phenomenological flow
c   suppression modification of conductivity. Inputs are not modified on output.
c
c                               NOMENCLATURE:
c                               ------------
c
c   Variable:               Description:
c   --------                -----------
c
c   r (1, ..., nj)          Rho grid position.
c
c   nj                      Number of grid points.
c
c   nout                    I/O unit number where data is to be printed.
c
c   tstaebler               Time at which coefficients were computed.
c
c   t                       Time point.
c
c ------------------------------------------------------------------ HSJ
c
c
      USE param
      USE numbrs
      USE mesh
      USE tordlrot
      USE staebler
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'staebler.i'                  ! pick up tstaebler
c      include 'tordlrot.i'
c
      call header (nout, tstaebler, t)
c
      write  (nout, 10)
   10 FORMAT (8x, 'PHENOMENOLOGICAL FLOW SHEAR SUPPRESSION CORRECTION' /
     .       20x, ' (half-grid quantities):')
c
      write  (nout, 20)
   20 format (4x, 'j', 3x, 'rho (cm)', 4x, 'rho/a', 6x, 'E-MULT',
     .        6x, 'I-MULT'3x, 'RGCMULT', 6x, 'RGC')
c
      do j=1, nj-1
        rjct = REAL (j) + 0.5
        rho  = (r(j) + r(j+1)) / 2.0
        droa = rho/r(nj)
        drgc = 0.5*(rgc(j)+rgc(j+1))
        re   = rgce_mult*rgc_mult(j)
        ri   = rgci_mult*rgc_mult(j)
        write  (nout, 30)  rjct, rho,droa, re, ri, rgc_mult(j), drgc
   30   format (1x, f4.1, 6(1pe11.3))
      end do
      write  (nout, 40)  rgca, rgcb, rgcc, rgcd
   40 format (/ 5x, 'MODEL PARAMETER VALUES:' /
     .            ' RGCA = ',1PE12.4,2X,'RGCB = ', 1pe12.4 /
     .            ' RGCC = ',1PE12.4,2X,'RGCD = ', 1pe12.4)
      write  (nout, 42)
   42 format (2x,'E-MULT is the multiplier for anomalous electron',
     .          ' conductivity' /
     .           'I-MULT is the multiplier for anomalous ion',
     .          ' conductivity')
      if (times_rgc .gt. -1.0e30 .and. timee_rgc .gt. -1.0e30) then
        write  (nout, 45)  times_rgc, timee_rgc
   45   format (' rgc is time average of angrot gradient',
     .          ' wrt normalized rho' /
     .          ' the time average was taken from ', 1pe12.6,
     .          ' to ', 1pe12.6)
      end if
      write  (nout, 50) rgc_string
   50 format (a)
c
      return
c
      end

      subroutine print_shay (nout, jprt, r, timeshay, t,
     .                       xkeshay, xkishay, xkeinv, xkiinv, xkeneo,
     .                       xkineo, psir, rmajorvec, grad_te, nj, te,
     .                       smult, skimult, snexp, sbpexp, srexp,
     .                       sbigrexp, stexp, sdtdrexp,
     .                       sdenscale, irinshay, iroutshay,ifsflag)
c
c ----------------------------------------------------------------------
c
c --- print out Hsieh conductivity results
c
c --- input
c
c  nout            fortran io unit number for output
c  r(j)            j1,2..nj,full grid in cm r refers to rho grid
c
c --- the following thermal conductivities are all in (1/cm-sec)
c
c  xkeshay(j)      normalized conductivites on full grid
c  xkishay(j)      (xkeshay and xkishay  INCLUDE neoclassical terms
c  xkeinv(j)       j = 1,2..nj,the power balance ke on full grid
c  xkiinv(j)                                      ki
c  xkeneo(j)
c  xkineo(j)       j = 1,2..nj,the neoclassical conductivities
c
c  timeshay           time,sec, at which Hsieh model was evaluated
c  t                  time index
c  psir(j)          psi values corresponding to rho(j)
c  rmajorvec(j)     major radius values corresponding to rho(j)
c  te(j)            j = 1,2..nj electron temp in keV on full grid
c
c --- the following are the exponents in Hsieh model
c  snexp
c  sbpexp
c  srexp
c  sbigrexp
c  stexp
c  sdtdrexp
c
c  sdenscale       electron density scale factor for Hsieh model
c  irinshay        the Hsieh model is assumed valid over the
c  iroutshay        range r(irinshay) to r(iroutshay)
c
c --- output (NONE)
c
c ------------------------------------------------------------------ HSJ
c
      implicit none
c
      integer  j,nj,nout,j1prt,jprt,irinshay,iroutshay,ifsflag
      real*8   r(*),xkeshay(*),xkishay(*),grad_te(*),
     .         xkeinv(*),xkiinv(*),xkeneo(*),xkineo(*),
     .         timeshay,psir(*),
     .         rmajorvec(*),te(*),smult,skimult,t,sdenscale,
     .         snexp,sbpexp,srexp,sbigrexp,stexp,sdtdrexp
c
c --- get grad te
c
      do j=1,nj
          if (j .gt. 1 .and. j .lt. nj) then
              grad_te(j) =  (te(j+1)-te(j-1))/(r(j+1)-r(j-1))
          else   if (j .eq. 1) then
              grad_te(j) =  0.0
          else   if (j .eq. nj) then
              grad_te(j) = (te(nj)-te(nj-1))/(r(nj)-r(nj-1))
          end if
      end do
c
c --- do the printing
c
        do j=1,nj
            if (j .eq. 1) then
                call header (nout,timeshay,t)
                if (ifsflag .eq. 0) then
                    write (nout, 10)timeshay
                else
                    write (nout, 11)timeshay
                end if
            end if
c
            j1prt = ((j-1) / jprt) * jprt
c
            if (j1prt .ne. j-1 .and. j .ne. nj)  go to 50
            if (j .eq. irinshay)  write (nout, 25)
            write (nout, 20) j, r(j), psir(j), rmajorvec(j),
     .                       grad_te(j), xkeshay(j),
     .                       xkeneo(j),xkeinv(j),xkishay(j),
     .                       xkineo(j),xkiinv(j)
            if (j .eq. iroutshay)  write (nout, 25)
   50     continue
        end do
c
   10 format (50(1h-),'Hsieh conductivity model ', 40(1h-)      /
     .        50x,' of confinement at time ',1pe12.4,' sec'     /
     .        45x,'(gradients calculated in rho space) '        /
     .         3x,'j',3x,'rho',13x,'psi',8x,'rmaj',3x,
     .          '  gradte',4x,'xkeshay',5x,'xkeneo',
     .         3x,'xkepwb',4x,'xkishay',4x,'xkineo',7x,'xkipwb' /
     .         7x,'cm',7x,'kgauss-cm**2',6x,'cm',6x,'keV/cm',
     .         4x,'(1/cm-s)', '  (1/cm-s)', '  (1/cm-s)',
     .         3x,'(1/cm-s)',3x,'(1/cm-s)',3x,'(1/cm-s)'        /
     .         2x,126(1h=))
   11 format (50(1h-),'Hsieh conductivity model ', 40(1h-)           /
     .        50x,' of confinement at time ',1pe12.4,' sec'          /
     .        45x,'(gradients calculated in rho space) '             /
     .        45x,'(FLOW SHEAR TURBULENCE SUPPRESION IS IN EFFECT) ' /
     .         3x,'j',3x,'rho',13x,'psi',8x,'rmaj',3x,
     .          '  gradte',4x,'xkeshay',5x,'xkeneo',
     .         3x,'xkepwb',4x,'xkishay',4x,'xkineo',7x,'xkipwb'      /
     .         7x,'cm',7x,'kgauss-cm**2',6x,'cm',6x,'keV/cm',
     .         4x,'(1/cm-s)', '  (1/cm-s)', '  (1/cm-s)',
     .         3x,'(1/cm-s)',3x,'(1/cm-s)',3x,'(1/cm-s)'             /
     .         2x,126(1h=))
c
   20 format (2x,i3,10(1pe11.3))
   25 format (10(1h-))
      write  (nout, 64)
 64   format (' IN SIMULATION MODE THE POWER BALANCE KE AND KI ' /
     .        ' WILL EQUAL THE HSIEH KE AND KI IF NO OTHER MODELS',
     .        ' ARE ACTIVE')
      write  (nout, 65) snexp,sbpexp,srexp,sbigrexp,stexp,sdtdrexp,
     .                  sdenscale,smult,skimult
   65 format ('  NOTE THAT THE ABOVE VALUES ARE VALID ONLY WITH',
     .         ' THE FOLLOWING SET OF EXPONENTS '              /
     .        10x,'snexp = ',f8.2, '  sbpexp   = ',f8.2        /
     .        10x,'srexp = ',f8.2, '  sbigrexp = ',f8.2        /
     .        10x,'stexp = ',f8.2, '  sdtdrexp = ',f8.2        /
     .        10x,'sdenscale = ',1pe12.4, '  smult = ',1pe12.4 /
     .        10x,'skimult = ',1pe12.4 ////)
      return
c
      end

      subroutine print_staebler (nout, t)
c
c --- S.J. Thompson --- General Atomics --- Core Physics Fusion Group --
c
c                           SUBPROGRAM DESCRIPTION:
c                           ----------------------
c
c   Print to output file results from calculations of the Staebler-
c   Hinton H-mode correction model. Inputs are not modified on output.
c
c                               NOMENCLATURE:
c                               ------------
c
c   Variable:               Description:
c   --------                -----------
c
c   r (1, ..., nj)          Rho grid position.
c
c   nj                      Number of grid points.
c
c   nout                    I/O unit number where data is to be printed.
c
c   tstaebler               Time at which coefficients were computed.
c
c   t                       Time point.
c
c   coefa (1, ..., nj-1)    Staebler-Hinton multiplicative factor for
c                           anomalous electron diffusivity term.
c
c   coefb (1, ..., nj-1)    Staebler-Hinton multiplicative factor for
c                           anomalous ion diffusivity term.
c
c   coefc (1, ..., nj-1)    Staebler-Hinton multiplicative factor for
c                           diffusion coefficient.
c
c   coefd (1, ..., nj-1)    Staebler-Hinton multiplicative factor for
c                           toroidal momentum diffusivity.
c
c   sperp (1, ..., nj-1)    Staebler-Hinton perpendicular shear term.
c
c --- S.J. Thompson ---------------------------------------------- 30 Mar 94 ---
c
      USE param
      USE numbrs
      USE mesh
      USE staebler
      implicit  integer (i-n), real*8 (a-h, o-z)
c     include 'param.i'
c      include 'mesh.i'
c      include 'numbrs.i'
c      include 'staebler.i'
c
      call header (nout, tstaebler, t)
c
      write  (nout, 10)
   10 FORMAT (7X,'STAEBLER-HINTON H-MODE CORRECTION',
     .           ' (half-grid dimensionless quantities):')
c
      write  (nout, 20)
** 20 format (4x, 'j', 3x, 'rho (cm)', 4x, 'coefa', 6x, 'coefb',
**** .        6x, 'coefc', 6x, 'coefd', 6x, 'sperp')
   20 format (4x, 'j', 3x, 'rho (cm)', 4x, 'CoefA', 6x, 'CoefB',
     .        6x, 'CoefC', 6x, 'CoefD', 6x, 'Sperp', 6x, 'Term1',
     .        6x, 'Term2', 6x, 'Term3')
c
      do j=1, nj-1
        rjct = REAL (j) + 0.5
        rho  = (r(j) + r(j+1)) / 2.0
****    write  (nout, 30)  rjct, rho, coefa(j), coefb(j), coefc(j),
**** .                     coefd(j), sperp(j)
        write  (nout, 30)  rjct, rho, coefa(j), coefb(j), coefc(j),
     .                     coefd(j), sperp(j),
     .                     term1(j), term2(j), term3(j)
** 30   format (1x, f4.1, 6(1pe11.3))
   30   format (1x, f4.1, 9(1pe11.3))
      end do
c
      write  (nout, *)
      do i=1,15
        write  (nout, 60)  i, staeblrmodl(i)
   60   format (3x,          'staeblrmodl(', i2, ') = ', a)
      end do
c
      write  (nout, 64)
   64 format (/ '           Bp   dw'     /
     .          ' Term1 = - -- * -- * R' /
     .          '           Bm   dr'      )
c
      write  (nout, 69)
   69 format (  '             1     d (ni)   d (pi)       c'    /
     .          ' Term2 = [ ----- * ------ * ------ ] * ------' /
     .          '           ni**2     dr       dr       e * Bm'  )
c
      write  (nout, 74)
   74 format (  '              1   d**2 (pi)       c'    /
     .          ' Term3 = - [ ---- --------- ] * ------' /
     .          '              ni    dr**2       e * Bm' /)
      return
c
      end

      subroutine print_weiland (r, nj, jprt, nout, t,ifsflag)
c
      USE param
      USE machin
      USE tcoef
      USE rhog
      USE weiland
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- print out Weiland Nordman results, calculated by subroutine WEILAND_12
c --- the values to be printed were calculated at the time stored in
c --- time_weiland
c ----------------------------------------------------------------------
c
c      include 'param.i'
c      include 'machin.i'
c      include 'weiland.i'
c      include 'rhog.i'
      include 'rebut.i'
      include 'storage.i'
c      include 'tcoef.i'
c
      dimension r(nj)
c
      do j=1,nj
          if (j .eq. 1) then
            call header (nout, time_weiland, t)
            write (nout, 10) time_weiland
          end if
c
          j1prt = ((j-1) / jprt) * jprt
c
          if (j1prt .ne. j-1 .and. j .ne. nj)  go to 50
          write (nout, 20) j, r(j), psir(j), rmajorvec(j),
     .                     d_weiland(j), xdchitot(j),xchie_weiland(j),
     .                     xchietot(j), xchii_weiland(j),
     .                     xchiitot(j),xeffwl(j)
c
   50     continue
        end do
c
      do j=1,nj
          if (j .eq. 1) then
            call header (nout, time_weiland, t)
            write (nout, 30) time_weiland
          end if
c
          j1prt = ((j-1) / jprt) * jprt
c
          if (j1prt .ne. j-1 .and. j .ne. nj)  go to 60
          write (nout, 20) j, r(j), psir(j), rmajorvec(j),
     .                     xke_weiland(j),
     .                     xketot(j), xki_weiland(j),
     .                     xkitot(j)
c
   60     continue
        end do
c
   10 format (1x, 38(1h-), 1x,
     .       'Weiland Model: Diffusivities ',
     .        30(1h-)                                        /
     .        31x, 'of confinement at time ',1pe12.4,' sec'  /
     .        3x,'j',6x,'rho',6x,'psi',10x,'rmaj',4x,
     .       ' d  wn  ', '  d  tot ',5x,'chie wn',6x,'chie tot',
     .        2x,'chii wn',4x,'chii tot ',3x,'chi eff'        /
     .        7x,'cm',7x,'kgauss-cm**2',7x,'cm',3x,'cm**2/sec',
     .        3x, 'cm**2/sec', 2x, 'cm**2/sec', '  cm**2/sec',
     .     '  cm**2/sec', '  cm**2/ sec' , '  cm**2/ sec '   /
     .        2x, 126(1h=))
c
   30 format (1x, 38(1h-), 1x,
     .       'Weiland Model: Conductivites ',
     .        30(1h-)                                        /
     .        31x, 'of confinement at time ',1pe12.4,' sec'  /
     .        3x,'j',3x,'rho',9x,'psi',8x,'rmaj',7x,
     .        4x,'ke  wn',6x,'ke tot',
     .        5x,'ki wn',5x,'ki tot ',5x,             /
     .        7x,'cm',7x,'kgauss-cm**2',4x,'cm',6x,'1/(cm sec)',
     .        x, '1/(cm sec)', 2x, '1/(cm sec', '  1/(cm sec)',  /
     .        2x, 126(1h=))
c
   20 format (2x, i3, 10(1pe11.3))
c
      return
c
      end

      subroutine ptorque (timet, t)
c
      USE param
      USE io 
      USE extra
      USE numbrs
      USE mesh
      USE tordlrot
      USE machin,          ONLY : volume 
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c --- print a table of torque densities to unit nout
c ------------------------------------------------------------------ HSJ
c
c      include 'param.i'                 ! get parameters
c      include 'tordlrot.i'              ! spbolt, storque, storqueb,smagtorque
c      include 'numbrs.i'                ! nj
c      include 'io.i'                    ! nout
c      include 'mesh.i'                  ! r(j)
c      include 'extra.i'                 ! fact(j)
  
c
      call header (nout, timet, t)
c
      write  (nout, 1)
    1 format ('  -------------------------- TORQUE DENSITIES, ',
     .        '  gram/(cm*sec2) ---------------------------' //
     .           6x,'j',10x,'r(j)',
     .           10x,'spbolt',10x,'storqueb',6x,'storque',8x,
     .           'smag_brake',5x,'storque_ntv')
    2 format (2x,i5,5x,f10.4,7x,5(e12.4,3x))
      do 100 j=1,nj
        j1prt = ((j-1)/jprt)*jprt
        if ((j1prt .ne. j-1) .and. (j .ne. nj))  go to 100
        write (nout, 2) j,r(j),spbolt(j),storqueb(j),storque(j),
     .                   smagtorque(j),sntvtorque(j)
  100 continue
      call trapv(r,spbolt,fact,nj,rtorque)
      call trapv(r,storque,fact,nj,ttorque)
      call trapv(r,storqueb,fact,nj,btorque)
      call trapv(r,smagtorque,fact,nj,tmagtorque)
!     when this routine is called fact is normalized to 1.0 hence:
      rtorque = rtorque*volume
      ttorque = ttorque*volume
      btorque = btorque*volume
      write  (nout, 3)  rtorque, btorque, ttorque, tmagtorque,ntvtorquet
    3 format ('  vol. intg (gram-cm2/sec2): ',5(e12.4,3x))
c
      write  (nout, 20)
   20 format (' sbpolt is the J x B torque density resulting',
     .        ' from (assumed radial) return current' /
     .        ' caused by PROMPT fast ion orbit loss')
      if (rtorque .eq. 0.0)  write (nout, 21)
   21 format (' spbolt is identically zero because either ' /
     .        ' iborb = 2 was turned off ' /
     .        ' or the beam is off at this time ')
      write  (nout, 22)
   22 format (' storqueb is the beam torque')
      write  (nout, 24)
 24   format (' smag_brake is torque due to magnetic brakeing ')
      write  (nout, 25)
 25   format (' storque_ntv is torque due to NTV ')
      write  (nout, 23)
   23 format (' storque is the total torque ' ///)
c
      return
c
      end

      subroutine quadft (x1, y1, x2, y2, x3, y3, cd, ce, cf, ieror)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c  this routine takes three points as input (x1,y1), (x2,y2), x3,y3)
c  and fits a parabola y = cd + ce*x + cf*x**2 through them.
c  ieror = 1 indicates the transformation matrix was singular, which
c  should happen only if two points given are identical
c
      dimension a(3,3), y(3), c(3)
c
      y(1) = y1
        y(2) = y2
        y(3) = y3
      a(1,1) = 1.0
        a(1,2) = x1
        a(1,3) = x1*x1
      a(2,1) = 1.0
        a(2,2) = x2
        a(2,3) = x2*x2
      a(3,1) = 1.0
        a(3,2) = x3
        a(3,3) = x3*x3
      call minv(a,3,ieror)                                       ! HSJrs
      if (ieror .eq. 1)  return
      call gmprd(a,y,c,3,3,1)
      cd = c(1)
      ce = c(2)
      cf = c(3)
      return
c
      end

      subroutine rebutlal (eta, curden, fcap, te, ti, q, zeff, fftrap,
     .                     ene, r, dr, atw, nj, nout, ncrt)
c
c
c ----------------------------------------------------------------------
c  written  by Holger St.John, General Atomics core fusion group
c  modified by S. J. Thompson 26 January 1993
c
c  in the critical electron gradient model of plasma transport developed by P.
c  Rebut, P. Lallia, and M. Watkins, anomalous transport is determined by tur-
c  bulence in the magnetic field topology which occurs above a certain thresh-
c  hold. in a tokamak the driving force is assumed to be the temperature grad-
c  ient. in this model the anomalous electron diffusivity is given by
c
c xan,e = c**2 * sqr (u0 * mp) * [1 - sqr (epsilon)] * sqr (1 + zeff) *  q**2  *
c         --------------------                                          ------
c          2 * btor * sqr (r)                                           grad q
c
c          sqr (te)  * | grad te + 2 * grad ne |
c          ---------   | -------       ------- |
c          sqr (ti)    |    te            ne   |
c
c  the ion term is defined by
c
c    xi = 2 * xan,e *       zi        * sqr (te) * factr * heavyte * heavyq
c                     ---------------   --------
c                     sqr (1 + zeff)    sqr (ti)
c
c  UPDATE 8/22/95 HSJ
c    Rebut-Lallia model was changed to reflect ion rho-star scaling.
c    The change is that xi above is replaced by:
c
c                          zi          2 * Te
c  xi = 2 * xan,e * --------------- * -------- * rescale * factr * heavyte * hea
c                   sqr (1 + zeff)     (Te+Ti)
c
c     where rescale = 0.3 * R0 * Bt0 / SQRT (Te + Ti)
c     (That is, SQRT (Te/Ti) gets replaced by (2 * Te / (Te+Ti)) and
c     a new factor, rescale, is multiplied in.)
c     This new model is now the default. To get back the old model,
c     use switch  rlw_model = 'old'  in the NAMELIS1 input namelist.
c
c  END OF 8/22/95 HSJ UPDATE
c
c   the electron contribution given by
c
c    xe = xan,e * factr * heavyte * heavyq
c
c  where we define
c
c                     factr = 1.0 - | (grad te)c |
c                                  --------------
c                                    | grad te |
c
c  and Heaviside functions
c
c                     heavyte = h (|grad te| - |(grad te)c|)
c
c                     heavyq = h (grad q).
c
c  the critical electron temperature gradient is given by
c
c   (grad te)c = 0.06 *     sqr (e**2)      * 1 * sqr [eta * j * (btor**3)]
c                       -------------------   -   ------------------------- .
c                       sqr [u0 * sqr (me)]   q      sqr [ne * sqr (te)]
c
c  for xi we assume a single hydrogenic ion species with mass number of primary
c  ion species atw and charge zi equal to 1
c
c  the ion heat and particle diffusivities are summed over all ion species,
c  eta is the neoclassical resistivity, and j
c  is the plasma current density. note that the critical electron temperature
c  gradient is positive definite, and thus has the opposite sign of
c  what grad te is normally.
c
c  references:  Rebut, P.H., Lallia, P.P., Watkins, M.L., in plasma physics and
c               controlled nuclear fusion research 1988 (proc. 12th int. conf.
c               nice, 1988), vol. 2, iaea, vienna (1989) 191
c
c               Rebut, P.H., Watkins, M.L, Gambier, D.J., Boucher, D., phys.
c               fluids b 3 (1991) 2209
c
c               Boucher, D., personal communication
c
c               for the HSJ 8/22/95 update, the reference is
c               Seville, IAEA, 1994, Rosenbluth et al.
c
c  input parameters:
c  ----------------
c
c  variable name:          description:                          units:
c  -------------           -----------                           -----
c
c    eta (nj)              neoclassical resistivity              seconds
c
c    curden (nj)           current density                       amp/cm**2
c                          note that curden= <jtoroidal*r0/r>
c                          where < > means the usual flux surface average
c
c    fcap (nj)             f(psilim) / f(psi)                    dimensionless
c
c    te (nj)               electron temperature                  keV
c
c    ti (nj)               ion temperature                       keV
c
c    q (nj)                safety factor                         dimensionless
c
c    eps (nj)              inverse aspect ratio                  dimensionless
c
c    zeff (nj)             z-effective                           dimensionless
c
c    fftrap (nj)           trapped particle fraction             dimensionless
c
c    ene (nj)              electron density                      number/cm**3
c
c    r (nj)                rho                                   cm
c    dr(nj-1)              rho(j+1)-rho(j)                       cm
c
c    rmajorvec (nj)        rmajor on outboard side of plasma     cm
c                          at elevation of magnetic axis
c
c  * psir (nj)             psi on outboard side of plasma        kgauss-cm**2
c
c                          at elevation of magnetic axis
c
c  * flim                  f(psilim) =  r0 * btor                  gauss*cm
c
c    atw                   mass number of primary ion species    dimensionless
c
c  qrebsmth          integer variable if ne 0 smooth q profile
c                    used below. otherwise don't smooth q.
c                    (qrebsmth is in rebut.i)
c  * psir and rmajorvec are in 'rhog.i'; flim in 'machin.i'
c
c   output
c   to INCLUDE file rebut.i:
c         crit_grad(j)           j = 1,2...nj the critical te gradient
c                                in units of keV/cm
c         grad_te(j)             the gradient of te (at rho position r(j))
c                                in units of keV/cm
c         xchierl(j)             Rebut-Lallia electron thermal diffusivity
c                                (without neoclassical term)
c                                in units of cm**2/sec
c         xchiirl(j)             Rebut-Lallia ion thermal diffusivity
c                                (without neoclassical term)
c                                in units of cm**2/sec
c         drebut(j)              Rebut-Lallia particle diffusion
c                                coefficient.
c                                in units of cm**2/sec
c
c     NOTE:
c         xkirebutj)             Rebut-Lallia ion thermal conductivity
c                                (without neoclassical term)
c                                in units of 1.0/cm*sec
c                                          and
c         xkerebutj)             Rebut-Lallia elec thermal conductivity
c                                (without neoclassical term)
c                                in units of 1.0/cm*sec
c                                are calculated outside of this routine
c  modified 1 june 1992:  variables crit_grad, grad_te, xchierl, and xchiirl
c                         added to plotting output.              s. j. thompson
c ----------------------------------------------------------------------
c
      USE param
      USE machin
      USE tcoef
      USE rhog
      USE neo2d
      USE replace_imsl,                     ONLY : my_icsmou,my_icsscv

      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
c      include 'machin.i'
c      include 'neo2d.i'                         ! pick up epsilon
      include 'rebut.i'
c      include 'rhog.i'
      include 'storage.i'
c      include 'tcoef.i'
c
      dimension r(nj), eta(nj), curden(nj), fcap(nj), te(nj), ti(nj),
     .          q(nj), zeff(nj), teprime(kj), eneprime(kj), qprime(kj),
     .          ene(nj), fftrap(nj),absq(kj), tiprime(kj),dr(*),
     .          qsmooth(kj),work(kj),cq(kj,3)
c
      equivalence (teprime(1),  xdum(1))      ! utilize existing storage
      equivalence (eneprime(1), ydum(1))
      equivalence (qprime(1),   zdum(1))
      equivalence (absq(1),     wdum(1))
      equivalence (tiprime(1),  vdum(1))
      equivalence (qsmooth(1),  sdum(1))
      equivalence (work(1),     udum(1))
      equivalence (cq(1,1),     tdum(1))
c
c ----------------------------------------------------------------------
c MKS constants
c ----------------------------------------------------------------------
c
      data u0     /1.25663706144e-06/  ! perm of free space (kg-m/cm**2)
      data c      /3.00e+08/           ! velocity of light (m/s)
      data emass  /9.11e-31/           ! electron rest mass (kg)
      data amu    /1.67e-27/           ! proton rest mass (kg)
      data echarge/1.60e-19/           ! electron charge (c)
c
c ----------------------------------------------------------------------
c check for valid storage
c ----------------------------------------------------------------------
c
      zero = 0.0
c
      if (kstore .lt. 3*kj) then
          write  (nout, 1)  kstore, 3*kj
          write  (ncrt, 1)  kstore, 3*kj
    1     format (' ERROR in subroutine REBUTLAL'             /
     .            ' matrix cq does not have sufficient space' /
     .            ' program must stop'                        /)
          call STOP ('subroutine REBUTLAL: problem #1', 41)
      end if
c
c ----------------------------------------------------------------------
c calculate constants in above expressions
c ----------------------------------------------------------------------
c
      amass      = atw * amu                        ! mp (kg)
      zprim      = 1.0                              ! zi (dimensionless)
      consta     = 0.5 * c * c * SQRT (u0 * amass)
      constc     = 0.06 * echarge / SQRT (u0 * SQRT (emass))
      rootrmajor = SQRT (0.01*rmajor)               ! meters**(1/2)
c
c ----------------------------------------------------------------------
c electron temperature gradient
c ----------------------------------------------------------------------
c
      call difydxhalfgrid (dr, te, teprime, nj)     ! keV per centimeter
c
c ----------------------------------------------------------------------
c ion temperature gradient
c ----------------------------------------------------------------------
c
      call difydxhalfgrid (dr, ti, tiprime, nj)     ! keV per centimeter
c
c ----------------------------------------------------------------------
c electron density gradient
c ----------------------------------------------------------------------
c
      call difydxhalfgrid (dr, ene, eneprime, nj)   ! number per cm
c
c ----------------------------------------------------------------------
c  gradient of safety factor
c  smooth q profile if called for (ie, qrebsmth = 1)
c ----------------------------------------------------------------------
c
      do i=1,nj
         absq   (i) = ABS (q(i))
         qsmooth(i) = absq(i)
      end do
      if (qrebsmth .ne. 0) then
            ier = 0
            if (qrebsmth .eq. 1) then
                ijob = 2    ! equally spaced r, more than 20 data points
                call my_icsscv(r,absq,nj,qsmooth,cq,kj,ijob,work,ier)
            else if (qrebsmth .eq. 2) then
                dis   = 1.0
                sc    = 0.0
                maxit = 8
                call my_icsmou (r,qsmooth,nj,dis,sc,maxit,work,ier)
            else                    ! qrebsmth > 2
                nrebs = qrebsmth/2
                do j=nrebs+1,nj-nrebs
                    qsum = 0.0
                    do k=j-nrebs,j+nrebs
                      qsum = qsum+qsmooth(k)
                    end do
                    qsmooth(j) = qsum/(2*nrebs+1)
                end do
            end if
            if (ier .ne. 0) then
              write  (nout, 2)  ier, qrebsmth
              write  (ncrt, 2)  ier, qrebsmth
    2         format (' ERROR in IMSL, called by REBUTLAL' /
     .                ' ier, qrebsmth =', i5, 2x, i5)
              call STOP ('subroutine REBUTLAL: problem #2', 42)
            end if
            do j=1,nj
              absq(j) = qsmooth(j)
            end do
      end if
c
      call difydxhalfgrid (dr, absq, qprime, nj)        ! per centimeter
c
c ----------------------------------------------------------------------
c  define electron and ion diffusivities at spatial locations
c  r(j+1/2),j = 1,2..nj-1. note that r(j) is the rho grid.
c ----------------------------------------------------------------------
c
      do j=1,nj-1
         tea = 0.5*(te(j+1)+te(j))
         tia = 0.5*(ti(j+1)+ti(j))
         enea = 0.5*(ene(j+1)+ene(j))
         curdena = 0.5*(curden(j+1)+curden(j))
         fcapa = 0.5*(fcap(j+1)+fcap(j))
         zeffa = 0.5*(zeff(j+1)+zeff(j))
         rmaj = 0.005 * (rmajorvec(j) +rmajorvec(j+1))
c
c ----------------------------------------------------------------------
c  spitzer resistivity (i.e., without trapped particle correction)
c ----------------------------------------------------------------------
c
        convrt = 8.98755179e+09     ! convert from seconds to ohm-meters
        res = eta(j) * convrt * (1.0 - fftrap(j))
c
        IF (res .lt. 0.0)then
           print *,'te =',te
           print *,'ftrap =',fftrap
           print *,'eta =',eta
           call STOP ('subroutine REBUTLAL: negative resistivity', 118)
        ENDIF
c
c ----------------------------------------------------------------------
c  current density
c  negative current density near the plasma edge is taken as positive
c  (this is in accord with what Boucher does, see "expressions used
c   in the code PRETOR by D.Boucher,(2/15/93)")
c ----------------------------------------------------------------------
c
        cur = ABS (curdena) * 1.0e+04  ! convert (amp/cm**2 to amp/m**2)
c
c ----------------------------------------------------------------------
c  toroidal magnetic field
c ----------------------------------------------------------------------
c
        bt = ABS (flim / rmaj / fcapa) * 1.0e-06                 ! tesla
c
c ----------------------------------------------------------------------
c  safety factor
c ----------------------------------------------------------------------
c
        qval = ABS (q(j))                                ! dimensionless
c
c ----------------------------------------------------------------------
c  critical electron termperature gradient
c ----------------------------------------------------------------------
c
        if (tea .le. 0.0) then
          call STOP ('subroutine REBUTLAL: non-positive elec.temp', 119)
        else if (enea .le. 0.0) then
          call STOP ('subroutine REBUTLAL: non-positive elec.dens', 120)
        else
          xte = tea  * 1.60219e-16  ! convert from keV to joules
          xne = enea * 1.00000e+06  ! convert from per cm**3 to per m**3
          critgrad = (constc / qval) * SQRT (res * cur * (bt**3))
          critgrad = critgrad / SQRT (xne * SQRT (xte))   ! joules/meter
c
c --- use exact form Boucher uses:
c --- bt in tesla, cur in ma/m**2, te in keV, ne in 10**19/m**3
c
             curma = cur * 1.0e-6    ! Ma/m**2
             xnea  = xne / 1.0e19
          critgrad = SQRT (bt*bt*bt*0.028*zeffa*curma/xnea)
     .                       *6.0/tea/qval    !    keV/m
          critgrad = critgrad * 1.60219e-16   ! joules/m
        end if
c
c ----------------------------------------------------------------------
c  electron temperature gradient
c ----------------------------------------------------------------------
c
        gradte = teprime(j) * 1.60219e-14   ! convert keV/cm to joules/m
c
c ----------------------------------------------------------------------
c  electron density gradient
c ----------------------------------------------------------------------
c
        gradne = eneprime(j) * 1.0e+08   ! convert per cm**4 to per m**4
c
c ----------------------------------------------------------------------
c  gradient of safety factor
c  if gradq is less than zero then Heaviside function, heavyq below,
c  sets the result to zero.
c ----------------------------------------------------------------------
c
        gradq = qprime(j) * 100.0             ! convert from 1/cm to 1/m
        gradq = MAX (gradq, zero)
c
c ----------------------------------------------------------------------
c  anomalous transport term
c ----------------------------------------------------------------------
c
        if (tia .le. 0.0) then
          call STOP ('subroutine REBUTLAL: non-positive ion temp', 121)
        else if (gradq .eq. 0.0) then
          xanom = 0.0
        else if (bt .eq. 0.0) then
          call STOP ('subroutine REBUTLAL: zero toroidal field', 122)
        else
            xte = tea  * 1.60219e-16     ! convert keV to joules
            xne = enea * 1.00000e+06     ! convert per cm**3 to per m**3
            qterm = 1.0 / (gradq/qval/qval+0.01)
            xanom = SQRT (tea / tia) * qterm
****        xanom = SQRT (tea / tia) * (qval**2) / gradq
            xanom = xanom * abs (gradte / xte + 2. * gradne / xne)
            epsb  = 0.5*(r(j+1)+r(j))/rmajor     ! Boucher value for EPS
            xanom = xanom * SQRT (1.0 + zeffa) * (1.0 - SQRT (epsb))
****        xanom = xanom * SQRT (1.0 + zeffa) * (1.0 - SQRT (eps(j)))
            xanom = xanom * consta / (bt * rootrmajor)   ! meters**2/sec
        end if
c
c ----------------------------------------------------------------------
c Heaviside function criterion
c ----------------------------------------------------------------------
c
        if (ABS (gradte) .le. ABS (critgrad)) then
            heavyte = 0.0
        else
            heavyte = 1.0
        end if
c
        if (qprime(j) .le. 0.0) then
            heavyq = 0.0
        else
            heavyq = 1.0
        end if
c
c ----------------------------------------------------------------------
c  electron diffusivity
c ----------------------------------------------------------------------
c
        if (gradte .eq. 0.0) then
c
            xe = 0.0
c
        else
c
            xe = xanom * (1.0 - (ABS (critgrad) / ABS (gradte)))
            xe = xe * heavyte * heavyq
c
        end if
c
c ----------------------------------------------------------------------
c ion diffusivity
c ----------------------------------------------------------------------
c
        if      (tia .le. 0.0) then
          xi = 0.0
        else if (tea .le. 0.0) then
          xi = 0.0
        else
          if (rlw_model .eq. 'old') then
            rescale = 1.0
            factr   = SQRT (tea / tia)
          else ! rlw_model .eq. 'new'
            rescale = 0.3*0.01*rmajor*bt / SQRT (tea + tia)! rmajor in m
            factr   = 2.0 * tea / (tea + tia)
          end if
          xi = xe * 2.0 * factr * rescale * zprim / SQRT (1.0 + zeffa)
        end if
c
c ----------------------------------------------------------------------
c  convert units
c ----------------------------------------------------------------------
c
        xi = xi * 1.0e+04           ! convert from m**2/sec to cm**2/sec
        xe = xe * 1.0e+04           ! convert from m**2/sec to cm**2/sec
        xchierl(j)   = xe
        xchiirl(j)   = xi
        drebut(j)    = 0.7*xchiirl(j)
        crit_grad(j) = critgrad/1.60219e-14 ! convert joules/m to keV/cm
          grad_te(j) = gradte  /1.60219e-14 ! convert joules/m to keV/cm
c
      end do                                ! end loop over half grid
c
      return
c
      end

      subroutine rfcur_simple (xnue, te, zeff, pehf, irfcur, model,
     .                         nj, kj, currfs)
c
c ----------------------------------------------------------------------
c --- subroutine uses simple formulae to roughly approximate
c --- the RF driven current given the ech power deposition profile.
c ----------------------------------------------------------------------
c
c --- input
c
c  xnue(j)      j = 1,2..nj ,reciprocal electron relaxation time,1/sec
c  te(j)                   electron temperature,keV
c  zeff(j)                 effective charge number
c  pehf(j)                 electron heating profile,keV/(cm**3-sec)
c  model                   index into RF arrays
c  irfcur(model)           1.0 for co current drive
c                          0.0     no
c                         -1.0    ctr
c
c  nj,kj                   used amd declared size of arrays
c
c --- output
c
c  currfs(j,model)         driven current, amps/cm**2
c
c ------------------------------------------------------------------ HSJ
c
      implicit none
c
      integer model, j, nj, kj
      real*8  xnue(*), te(*), zeff(*), pehf(*), currfs(kj,*),
     .        u1sq, charge, ve, irfcur(*)
c
      u1sq   = 5.0                 ! assume for now
      charge = 1.602e-19           ! electronic charge in coulombs
      do j=1,nj
        ve = 1.33e9 * SQRT (te(j)) ! elec velocity in cm/sec, te in keV
        currfs(j,model) = irfcur(model) * 1.5 * charge * ve * u1sq
     .                  * pehf(j) / ((5.0 + zeff(j)) * xnue(j) * te(j))
      end do
      return
c
      end

      subroutine shay_chie (te, ene, bp0, r, ra, dr, rmajor, rmbpav, nj,
     .                      smult, snexp, sbpexp, srexp, sbigrexp,
     .                      stexp, sdtdrexp, srin, srout, suserho,
     .                      schie, schii, skimult, irinshay, iroutshay,
     .                      nout, ncrt, sdenscale, ishayform, skimass)
c
c ----------------------------------------------------------------------
c  evaluate Hsieh ("shay") model of thermal diffusivity:
c  define the factor f by
c  f = smult*(ne**(snexp-1)*((-dte/drho)**sdtdrexp)*(rho**srexp)/(te**stexp)
c   then if suserho = true set
c       schie = f* <(R**(2-sbigrexp))*(bp**(2-sbpexp))>/(rmajor*bp0)**2
c   else if suserho = false use
c       schie = f/((rmajor**sbigrexp)*(bp0**sbpexp))
c  the constant smult must carry units, depending on the exponents
c  in the model, so that schie has units of cm**2/sec
c  for more information see TDN notebook, vol ii, pgs 55-58.
c ------------------------------------------------------------------ HSJ
c
c  UPGRADE:   August 1, 1995     Daniel Finkenthal
c
c  A dimensionally-correct form of the model can be invoked by setting
c  ishayform = 1. The old model will be multiplied by the factor sunits:
c
c       sunits = 8 * pi * kboltzman * clight * SQRT (Te/m_e)
c
c  where m_e = 511 keV and Te is the electron temperature.
c  Otherwise, sunits = 1.0 and the model is calculated as before.
c
c  The multiplier SKIMULT for the ion thermal diffusivity is now
c  scaled by a factor of
c                           1 / SQRT(A)
c
c  to make the model species-independent, where A is the atomic
c  mass number of the ions, as specified by SKIMASS. Default is
c  SKIMASS = 2.0 for Deuterium.
c
c  Otherwise, sunits is automatically set to unity if the user selects
c  the old version of the model (ISHAYFORM = 0).
c
c  Note that the dimensionless form of the model only applies for the
c  default exponents set in subroutine INIT as follows:
c            snexp     =  2.0     ! density exponent
c            sbpexp    =  2.0     ! bp exponent
c            srexp     =  3.0     ! rho exponent
c            sbigrexp  =  1.0     ! R exponent
c            stexp     =  1.0     ! te exponent
c            sdtdrexp  =  2.0     ! dte/dr exponent
c
c ------------------------------------------------------------------ DFF
c
c --- input:
c
c   te(j)        j = 1,2..nj electron temp on full grid (keV)
c   ene(j)                            density           (#/cm**3)
c   bp0(j)       flux surf. avg. poloidal b field (gauss)
c   r(j)         rho grid   (cm)
c   ra(j)        j = 1,2..nj-1 half grid in rho   (cm)
c   dr(j)        j = 1,2..nj-1, dr(j) = r(j+1)-r(j)   (cm)
c   rmajor       R0       (cm)
c   rmbpav(j)    required if suserho = .true.
c                rmbpav is <(R**(2-sbigrexp))*(bp**(2-sbpexp))>
c                on the (full) rho grid    (R in cm,bp in gauss)
c  snexp         exponent of ene(j)
c  sbpexp        exponent of bp0 (or bp if suserho = true,see above)
c  smult         multiplier of model
c  srexp         exponent of rho
c  sbigrexp      exponent of R
c  stexp         exponent of (1/te)
c  sdtdrexp      exponent of (-dt/dr)
c  srin          inner value of normalized rho
c  srout         outer value of normalized rho
c                (the model is active for srion<r/ra<srout )
c  suserho       explained above
c  nout          the outone file
c  ncrt          Fortran I/O unit number for error messages
c  sdenscale     density is in units of sdenscale
c  ishayform     explained above
c  skimult       Multiplier for the ion componenent of the model,
c                expressed in terms of the electron multiplier SMULT
c  skimass       A user input ion mass (in atomic units) used to make
c                the skimult parameter species-independent. This should
c                be set to the atomic mass number of the primary ion
c                species. (Default is 2.0 for Deuterium, but is automatically
c                set to 1.0 for ISHAYFORM = 0 in routine INIT if the user
c                requested the old model ( ISHAYFORM = 0 ).
c
c  INCLUDE file storage.i is included here for temporary storage of some
c                         intermediate results not required after this
c                         subroutine returns
c --- output:
c
c  schie(j)     j = 1,2..nj-1,electron thermal diffusivity on the half
c               grid (i.e., ra(j) grid) in units of cm**2/sec
c  schii(j)     same for ion thermal diffusivity
c  irinshay     index for the half grid vector ra(j).
c               ra(irinshay)/r(nj) is the smallest value of rho
c               such that ra(irinshay)/r(nj) .ge. srin
c  iroutshay    index for the half grid vector ra(j).
c               ra(iroutshay)/r(nj) is the largest value of rho
c               such that ra(iroutshay)/r(nj) .le. srout
c
c ----------------------------------------------------------------------
c
      implicit none
c
      logical  suserho
      integer  j, nj, irinshay, iroutshay, checklimits,
     .         nout, ncrt, kinc, ishayform
      real*8   ene(*),r(*),ra(*),te(*),bp0(*),schie(*),schii(*),
     .         dr(*),rmbpav(*),bpa,enea,sdenscale,
     .         f,rmajor,rmbpava,sbigrexp,sbpexp,sdtdrexp,skimult,
     .         smult,snexp,srexp,srin,srout,stexp,tea,
     .         rin,teprime(1000),
     .         pi, kbolt, clight, mekev, sunits, skimass
c
      include 'storage.i'
c
      equivalence (xdum(1), teprime(1))
      data         checklimits /0/
c
      pi     = ACOS (-1.0)
      clight =        2.99792e10
      kbolt  =        1.60219e-9
      mekev  =      511.003
c
c --- on first call to this routine check that srin and srout are set correctly
c
      if (checklimits .eq. 0) then
        irinshay  = nj
        iroutshay = 0
        kinc      = 0
        do j=1,nj-1
          rin   = ra(j)/r(nj)
          if (rin .ge. srin .and. rin .le. srout) then
            kinc      = kinc + 1
            irinshay  = MIN (irinshay , j)
            iroutshay = MAX (iroutshay, j)
          end if
        end do
        if (kinc .eq. 0) then
          write  (nout, 100)  srin, srout
          write  (ncrt, 100)  srin, srout
  100     format (' subroutine SHAY_CHIE reports an input error:' /
     .            '   srin and srout are set incorrectly'         /
     .            '   srin,    srout =', 2(2x, 1pe12.4))
          call STOP ('subroutine SHAY_CHIE: unspecified problem', 123)
        else
          checklimits = 1
        end if
      end if
c
c --- get d(te)/drho on half grid (units of keV/cm)
c
      call difydxhalfgrid (dr, te, teprime, nj)
c
c --- loop over the half grid
c
      do j=1,nj-1
        if ( teprime(j) .gt. 0.0)  teprime(j) = -teprime(j) ! HSJ8/30/98
        if ((teprime(j) .lt. 0.0)  .and. (j .ge. irinshay)
     .                             .and. (j .le. iroutshay)) then
          enea = 0.5*(ene(j+1)+ene(j))/sdenscale
          tea  = 0.5*(te(j+1)+te(j))
          bpa  = ABS (0.5*(bp0(j+1)+bp0(j)))
          f    = smult*(enea**(snexp-1))
          f    = f * (ra(j)**srexp)
          f    = f * ((-teprime(j))**sdtdrexp)
          f    = f / (tea**stexp)
          if (ishayform .eq. 1) then
            sunits = SQRT (tea / mekev)
            sunits = 8.0 * pi * kbolt * clight * rmajor * sunits
          else
            sunits = 1.0
          end if
          f = sunits * f
          if (suserho) then
            rmbpava  = 0.5 * (rmbpav(j+1)+rmbpav(j))
            schie(j) =   f * rmbpava/((rmajor*bpa)**2)
          else
            schie(j) = f / ((rmajor**sbigrexp)*(bpa**sbpexp))
          end if
        else
          schie(j) = 0.0
        end if
        schii(j) = schie(j) * skimult / SQRT (skimass)
      end do
      return
c
      end

      subroutine shay_multiplier (srincev,sroutcev,nout,r,timeshay,
     .             t,xkeshay,xkishay,xkeinv,xkiinv,xkeneo,xkineo,
     .             psir,rmajorvec,grad_te,ikeuse,ikiuse,xkepwb,xkipwb,
     .             nj,smult,skimult,smultstder,skimultstder,te,
     .             slim95, skilim95, irinshaycev, iroutshaycev, skimass,
     .             snexp,sbpexp,srexp,sbigrexp,stexp,sdtdrexp,
     .             sdenscale,irinshay,iroutshay,smulta,skimulta,
     .             maxtimes,ktimes,ifsflag)
c
c ----------------------------------------------------------------------
c
c --- given the power balance electron and ion cond.,compute
c --- least squares approximations to the multipliers smult
c --- and skimult needed to define the Hsieh model. Note that the
c --- least squares methods used here are a linearized form of the
c --- actual least squares required to find smult and skimult
c --- simultaneously. This linearized model should be adequate. Its
c --- drawback is that the covariance between smult and skimult is
c --- lost in doing the problem this way. Its advantage is that we
c --- don't get into the vagaries on non linear least squares fitting
c --- the data for the least squares problem is taken from those
c --- normalized rho points for which
c ---               srincev .le. rho .le. sroutcev
c --- if srincev = sroutcev then the single closest rho point
c --- is used for the evaluation of the constants. in this case
c --- a least squares analysis is not done.
c
c  UPGRADE:   August 1, 1995     Daniel Finkenthal
c
c  A dimensionally-correct form of the Hsieh model can be invoked by
c  setting ishayform = 1.  (See subroutine SHAY_CHIE.)
c  The multiplier SKIMULT for the ion thermal diffusivity is now
c  scaled by a factor of
c                             1 / SQRT (A)
c
c  to make the model species-independent, where A is the atomic
c  mass number.
c  Although not obvious, the original version of this routine does not
c  use the ion conductivity calculated by the shay_chie routine as it
c  does for the electron component in the determination of SMULT. Rather,
c  SKIMULT is calculated strictly in terms of the power-balance (PB) ion
c  conductivity and the Hseih electron conductivity, thereby bypassing the
c  ion cond. calculated in shay_chie. Therefore, in order to correcty obtain
c  the mass-dependent dimensionally correct ion conductivity factor SKIMULT
c  we must include the  1 / SQRT (A)  factor here, where A = SKIMASS.
c  for backwards compatability, SKIMASS will be automatically reset to 1.0 in
c  the initialization routine for the old version of the model (ISHAYFORM = 0)
c
c --- input
c
c  srincev         these are the values of normalized rho over
c  sroutcev        which the least squares approximation will be made
c  nout            fortran io unit number for output
c  r(j)            j1,2..nj,full grid in cm r refers to rho grid
c
c --- the following thermal conductivities are all in (1/cm-sec)
c
c  xkeshay(j)      non-normalized conductivites on full grid
c  xkishay(j)      (xkeshay and xkishay DO NOT include neoclassical terms
c                     when scsmult = true,which is the only time this
c                     routine is called).
c  xkeinv(j)       j = 1,2..nj,the power balance ke on full grid
c  xkiinv(j)                                      ki
c  xkeneo(j)
c  xkineo(j)       j = 1,2..nj,the neoclassical conductivities
c
c  timeshay           time,sec, at which Hsieh model was evaluated
c  t                  time index
c  psir(j)          psi values corresponding to rho(j)
c  rmajorvec(j)     major radius values corresponding to rho(j)
c  te(j)            j = 1,2..nj electron temp in keV on full grid
c
c --- the following are the exponents in Hsieh model
c
c  snexp
c  sbpexp
c  srexp
c  sbigrexp
c  stexp
c  sdtdrexp
c
c  sdenscale     electron density scale factor for Hsieh model
c  irinshay      the Hsieh model is assumed valid over the
c  irinshay      range r(irinshay) to r(iroutshay)
c
c  skimass       A user input ion mass (in atomic units) used to make
c                the skimult parameter species-independent. This should
c                be set to the atomic mass number of the primary ion
c                species. (Default is 2.0 for Deuterium, but is automatically
c                set to 1.0 for ISHAYFORM = 0 in routine INIT if the user
c                requested the old model ( ISHAYFORM = 0 ).
c
c --- temporary local storage needed while this subroutine is active
c
c  xkepwb(j)            j = 1,2..nj
c  xkipwb(j)
c  grad_te(j)
c  ikeuse(j)
c  ikiuse(j)
c
c  maxtimes          max storage available
c
c --- output
c
c  smulta(k)         least squares approximation of the multiplier
c                    smult. note that all points in the range
c                    irinshaycev to iroutshaycev were used to generate
c                    the least squares approximation.
c  skimulta(k)       same for ion multiplier
c  smultstder(k)     standard error of smult
c  skimultstder(k)   standard error of skimult
c  slim95(k)         95% confidence limit on smult
c  skilim95(k)       95% confidence limit on skimult
c  irinshaycev       starting index
c  iroutshaycev      ending index used for least squares calculations
c  xkeshay(j)         j = 1,2..nj,normalized xke including neoclassical
c  xkishay(j)         j = 1,2..nj,normalized xki including neoclassical
c
c ------------------------------------------------------------------ HSJ
c
      implicit none
c
      logical checklimits, leastsquares
      integer j,nj,irinshaycev,iroutshaycev,kinc,nout,ifsflag,
     .        ike,iki,ikefree,ikifree,ikeuse(*),ikiuse(*),
     .        j1prt,jprt,irinshay,iroutshay,ktimes,maxtimes
      real*8  r(*),xkeshay(*),xkishay(*),grad_te(*),
     .        xkeinv(*),xkiinv(*),xkeneo(*),xkineo(*),
     .        smult,skimult,d,dmin,rin,student(34),slim95(*),
     .        skilim95(*),smultstder(*),skimultstder(*),timeshay,
     .        xkepwb(*),xkipwb(*),te(*),diff,stud,rmajorvec(*),
     .        sresidsq,sumkesq,sumkisq,sumksq,psir(*),
     .        srincev,sroutcev,t,sdenscale,
     .        snexp,sbpexp,srexp,sbigrexp,stexp,sdtdrexp,
     .        smulta(*), skimulta(*), skimass
      data    checklimits /.true./
c
c --- student's t table for 95% confidence
c
      data student /12.706,4.303,3.182,2.776,2.571,2.447,2.365,2.306,
     .               2.262,2.228,2.201,2.179,2.160,2.145,2.131,2.120,
     .               2.110,2.101,2.093,2.086,2.080,2.074,2.069,2.064,
     .               2.060,2.056,2.052,2.048,2.045,2.042,2.021,2.000,
     .               1.980,1.960/
c
c --- on first call to this routine check that srincev and sroutcev
c --- are set correctly
c
      if (checklimits) then
          leastsquares = .true.    ! try to use least squares estimators
          irinshaycev  = nj
          iroutshaycev = 0
          kinc         = 0
          do j=1,nj
            rin = r(j) / r(nj)
            if (rin .ge. srincev .and. rin .le. sroutcev) then
              kinc         = kinc + 1
              irinshaycev  = MIN ( irinshaycev, j)
              iroutshaycev = MAX (iroutshaycev, j)
            end if
          end do
          if (kinc .le. 1) then
            dmin = 1.0
            do j=1,nj
              rin  = r(j) / r(nj)
              d    = ABS (rin - srincev)
              dmin = MIN (d, dmin)
              if (d .eq. dmin) then
                irinshaycev  = j
                iroutshaycev = j
              end if
            end do
            leastsquares = .false.        ! use single point estimator
          end if
          checklimits = .false.
      end if
c
      ktimes = ktimes + 1
      ktimes = MIN (ktimes, maxtimes)   ! if too many skip intermediates
      if (leastsquares) then
c
c --- least squares solution for the constants in Hsieh model
c --- define power balance minus neoclassical conductivites
c --- obtain the pointers ikeuse and ikiuse which will point
c --- to data that can be used in the fits
c
            ike = 0
            iki = 0
            do j=1,nj
c
c               neoclassical power balance
c
                xkepwb(j) = xkeinv(j) - xkeneo(j)
                xkipwb(j) = xkiinv(j) - xkineo(j)
                ikeuse(j) = 0
                ikiuse(j) = 0
                if (xkepwb(j) .gt. 0.0)  ikeuse(j) = 1
                if (xkipwb(j) .gt. 0.0)  ikiuse(j) = 1
                if (j .lt. irinshaycev) then
                  ikeuse(j) = 0
                  ikiuse(j) = 0
                else if (j .gt. iroutshaycev) then
                  ikeuse(j) = 0
                  ikiuse(j) = 0
                end if
                ike = ike + ikeuse(j)
                iki = iki + ikiuse(j)
            end do
c
c --- we have ike data values to determine the single constant
c --- smult, and similarly for the ion multiplier skimult
c
            ikefree = ike - 1
            ikifree = iki - 1
            if (ikefree .lt. 1 .or. ikifree .lt. 1) then
              leastsquares = .false.
              go to 1000
            end if
c
c --- get least squares solution for electron cond. multiplier
c
            sumkesq = 0.0
            sumksq  = 0.0
            do j=1,nj
              if (ikeuse(j) .gt. 0) then
                sumkesq = sumkesq+xkeshay(j)**2
                sumksq  = sumksq +xkeshay(j)*xkepwb(J)
              end if
            end do
            smult = sumksq / sumkesq  ! least squares estimate for smult
            sresidsq = 0.0            ! get sum of residuals squared
            do j=1,nj
              if (ikeuse(j) .gt. 0)
     .            sresidsq = sresidsq+(xkepwb(j)-smult*xkeshay(j))**2
            end do
c
c --- get standard error of smult
c
            smultstder(ktimes) = SQRT (sresidsq/(ikefree*sumkesq))
c
c --- get 95% confidence limits for smult,
c --- (there is 95% confidence that the "true" smult is in the interval
c --- [smult-slim95,smult+slim95] )
c
            if (ikefree .le. 30) then
              stud = student(ikefree)
            else if (ikefree .le. 40) then
              diff = student(31)-student(30)
              stud = student(30)+(ikefree-30)*diff
            else if (ikefree .le. 60) then
              diff = student(32)-student(31)
              stud = student(31)+(ikefree-40)*diff
            else if (ikefree .le. 120) then
              diff = student(33)-student(32)
              stud = student(32)+(ikefree-60)*diff
            else ! ikefree large enough for standard normal distribution
              stud = student(34)
            end if
c
            slim95(ktimes) = smultstder(ktimes)*stud
c
c --- finally get real xkeshay
c
            do j=1,nj
              xkeshay(j) = xkeshay(j) * smult
            end do
c
c --- get least squares solution for ion cond. multiplier
c
            sumkisq = 0.0
            sumksq  = 0.0
            do j=1,nj
              if (ikiuse(j) .gt. 0) then
                sumkisq = sumkisq + xkeshay(j)**2
                sumksq  = sumksq  + xkeshay(j) * xkipwb(j)
              end if
            end do
c
c           least squares estimate for skimult
c
            skimult  = sumksq / sumkisq / SQRT (skimass)
            sresidsq = 0.0                ! get sum of residuals squared
            do j=1,nj
              if (ikiuse(j) .gt. 0)
     .          sresidsq = sresidsq + (xkipwb(j) - skimult
     .                               * xkeshay(j))**2
            end do
c
c --- get standard error of skimult
c
            skimultstder(ktimes) = SQRT (sresidsq/(ikifree*sumkisq))
c
c --- get 95% confidence limits for skimult,
c --- (there is 95% confidence that the "true" skimult is in the interval
c --- [skimult-skilim95,skimult+skilim95] )
c
            if (ikifree .le. 30) then
              stud = student(ikifree)
            else if (ikifree .le. 40) then
              diff = student(31)-student(30)
              stud = student(30)+(ikifree-30)*diff
            else if (ikifree .le. 60) then
              diff = student(32)-student(31)
              stud = student(31)+(ikifree-40)*diff
            else if (ikifree .le. 120) then
              diff = student(33)-student(32)
              stud = student(32)+(ikifree-60)*diff
            else ! ikifree large enough for standard normal distribution
              stud = student(34)
            end if
            skilim95(ktimes) = skimultstder(ktimes)*stud
      end if
c
 1000 if (.not. leastsquares) then
        smult                = xkepwb (irinshaycev)
     .                       / xkeshay(irinshaycev)
        smultstder  (ktimes) = 0.0                 ! can't be determined
        slim95      (ktimes) = 0.0                 ! can't be determined
        skimult              = xkipwb(irinshaycev)
     .             / (smult * xkeshay(irinshaycev))
     .             /  SQRT (skimass)
        skimultstder(ktimes) = 0.0
        skilim95    (ktimes) = 0.0
      end if
c
c --- done with evaluation of multipliers. now do some printing of results
c
c --- get full conductivities, including neoclassical
c
      do j=1,nj
        xkeshay(j) =           xkeshay(j) + xkeneo(j)
        xkishay(j) = skimult * xkeshay(j) + xkineo(j)
        if      (j .gt. 1 .and. j .lt. nj) then
          grad_te(j) =  (te(j+1)-te(j-1))/(r(j+1)-r(j-1))
        else if (j .eq. 1) then
          grad_te(j) =  0.0
        else if (j .eq. nj) then
          grad_te(j) = (te(nj)-te(nj-1))/(r(nj)-r(nj-1))
        end if
      end do
c
c --- do the printing
c
        jprt = 1           ! note this value is LOCAL to this subroutine
        do j=1,nj
            if (j .eq. 1) then
                call header (nout,timeshay,t)
                if (ifsflag .eq. 0) then          ! flow shear model OFF
                     write (nout, 10) timeshay
                else                              ! flow shear model ON
                     write (nout, 11) timeshay
                end if
            end if
c
            j1prt = ((j-1) / jprt) * jprt
c
            if (j1prt .ne. j-1 .and. j .ne. nj)  go to 50
            if (    j .eq. irinshay           )  write (nout, 25)
            write (nout, 20) j, r(j), psir(j), rmajorvec(j),
     .                       grad_te(j), xkeshay(j),
     .                       xkeneo(j),xkeinv(j),xkishay(j),
     .                       xkineo(j),xkiinv(j),ikeuse(j),ikiuse(j)
            if (j .eq. iroutshay)  write (nout, 25)
   50     continue
        end do
c
   10 format (50(1h-),'Hsieh conductivity model ', 40(1h-) /
     .       50x,' of confinement at time ',1pe12.4,' sec' /
     .       45x,'(gradients calculated in rho space) '    /
     .       3x,'j',3x,'rho',13x,'psi',8x,'rmaj',3x,
     .        '  gradte',4x,'xkeshay',5x,'xkeneo',
     .       3x,'xkepwb',4x,'xkishay',4x,'xkineo',7x,'xkipwb',
     .       5x,'ike', '  iki'                             /
     .       7x,'cm',7x,'kgauss-cm**2',6x,'cm',6x,'keV/cm',
     .       4x,'(1/cm-s)', '  (1/cm-s)', '  (1/cm-s)',
     .       3x,'(1/cm-s)',3x,'(1/cm-s)',3x,'(1/cm-s)'     /
     .       2x,126(1h=))
   11 format (50(1h-),'Hsieh conductivity model ', 40(1h-)    /
     .       50x,' of confinement at time ',1pe12.4,' sec'    /
     .       45x,'(gradients calculated in rho space) '       /
     .       45X,'(FLOW TURBULENCE SUPPRESSION IS IN EFFECT)' /
     .       3x,'j',3x,'rho',13x,'psi',8x,'rmaj',3x,
     .        '  gradte',4x,'xkeshay',5x,'xkeneo',
     .       3x,'xkepwb',4x,'xkishay',4x,'xkineo',7x,'xkipwb',
     .       5x,'ike', '  iki'                                /
     .       7x,'cm',7x,'kgauss-cm**2',6x,'cm',6x,'keV/cm',
     .       4x,'(1/cm-s)', '  (1/cm-s)', '  (1/cm-s)',
     .       3x,'(1/cm-s)',3x,'(1/cm-s)',3x,'(1/cm-s)'        /
     .       2x,126(1h=))
c
   20 format (2x,i3,10(1pe11.3),4x,i2,2x,i2)
   25 format (10(1h-))
      if (leastsquares) then
        write (nout,30)ikefree,smult,smultstder(ktimes),
     .         smult-slim95(ktimes),     smult+slim95(ktimes)
        write (nout,40)ikifree,skimult,skimultstder(ktimes),
     .       skimult-skilim95(ktimes), skimult+skilim95(ktimes)
      else
          write (nout, 60)  smult, skimult
      end if
      write  (nout, 65) snexp, sbpexp, srexp, sbigrexp, stexp,
     .                  sdtdrexp, sdenscale
   30 format ('  least squares estimator results for electron',
     .        '  thermal conductivity multiplier (nfree = ',i5,')' /
     .       10x,'smult =',1pe12.4,2x,
     .        '  std. error =',1pe12.4                             /
     .       10x,'Assuming the model is correct there is a 95%',
     .        '  probability that the true value'                  /
     .       10x,'of smult is in the range (',1pe12.4,',',1pe12.4,')')
   40 format ('  least squares estimator results for ion',
     .        '  thermal conductivity multiplier (nfree = ',i5,')' /
     .       10x,'skimult =',1pe12.4,2x,
     .        '  std. error =',1pe12.4                             /
     .       10x,'Assuming the model is correct there is a 95%',
     .        '  probability that the true value'                  /
     .       10x,'of skimult is in the range (',1pe12.4,',',1pe12.4,')')
   60 format ('  single point estimate of Hsieh multipliers '       /
     .        '  smult =',1pe12.4, '  skimult = ',1pe12.4)
   65 format ('  NOTE THAT THE ABOVE VALUES ARE VALID ONLY WITH ',
     .        '  THE FOLLOWING SET OF EXPONENTS '       /
     .        10x,'snexp = ',f8.2, '  sbpexp   = ',f8.2 /
     .        10x,'srexp = ',f8.2, '  sbigrexp = ',f8.2 /
     .        10x,'stexp = ',f8.2, '  sdtdrexp = ',f8.2 /
     .        10x,'sdenscale = ',1pe12.4 ////)
c
c --- reset smult and skimult to 1.0,which is the value we want to
c --- use in subroutine SHAY_CHIE, when we are finding smult and skimult
c --- (i.e., when scsmult = true, which is the only time this routine is called)
c
      smulta(ktimes)   = smult
      skimulta(ktimes) = skimult
      smult            = 1.0
      skimult          = 1.0
      return
c
      end

      subroutine smthsj (x, y, n, nsmooth)
c
      implicit none
c
c     smoothing routines for the array of n data points, (x,y)
c     x is assumed to be equal-spacing, and not used in the routine.
c     five-point smoothing procedure with weight (0.1, 0.15, 0.5, 0.15, 0.1)
c     is used for the interior points.  Both end points are not modified.
c     The next boundary points are smoothed with the weight (0.25, 0.5, 0.25).
c     nsmooth is the number of paths of smoothing.
c     If nsmooth = 0, then no smoothing will be done. .. HSJ and YRLL (9 Jun 94)
c
c      integer     nmax,
c      parameter  (nmax = 101)
      real*8      x(*), y(*)
c      real *8 , dimension(:),allocatable :: yt      
      integer     i, n, npath, nsmooth,istat
       real *8 , dimension(n) :: yt  
c
c      if( .not. allocated(yt))then
c         allocate (yt(1:n),STAT = istat)
c         if(istat .ne. 0)
c     .          call allocate_error("yt in smthsj",0,istat) 
c      endif
      if (nsmooth .gt. 0) then
        if (n .lt. 5 )
     .    call STOP ('subroutine SMTHSJ: n < 5 ', 125)
        do npath=1,nsmooth
          do i=1,n
            yt(i) = y(i)
          end do
          y(  1) = yt(1)
          y(  2) = 0.25*yt(  1) + 0.50*yt(  2) + 0.25*yt(3)
          y(n-1) = 0.25*yt(n-2) + 0.50*yt(n-1) + 0.25*yt(n)
          y(n  ) = yt(n)
          do i=3,n-2
            y(i) = 0.10*yt(i-2) + 0.15*yt(i-1) + 0.50*yt(i) +
     .             0.15*yt(i+1) + 0.10*yt(i+2)
          end do
        end do
      end if
      return
c
      end

      subroutine sxrcal (codeid, kappa, ene, jsxr, nj, nx, ny,
     .                   p, pmax, psir, r, rmajor, rin, rmax, te,
     .                   rmhdgrid, zmhdgrid, zax, zmin, zmax,
     .                   idiode, narray, namar, ndiode, roamin, sxr)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c     This subroutine calculates the unnormalized SXR signals from the old
c     radial and vertical SXR diode arrays on Doublet III (if jsxr = 1)
c     or the new side and top diode arrays on Doublet III (if jsxr = 2).
c     A constant impurity enhancement factor is assumed, although the
c     resulting units are arbitrary.
c
      parameter   (kdi = 16, kar = 4, ksxr = 2)
      real*8       kappa
      character*8  codeid, namar, name
      dimension    namar(kar), ndiode(kar), idiode(kdi,kar)
      dimension    ene(*), psir(*), r(*), te(*)
      dimension    p(nx,*), rmhdgrid(*), zmhdgrid(*)
      dimension    roamin(kdi,kar), sxr(kdi,kar)
      dimension    nar(ksxr), name(kar,ksxr), ndiod(kar,ksxr),
     .             idiod(kdi,kar,ksxr)
      dimension    rw(kar,ksxr), zw(kar,ksxr), bw(kar,ksxr),
     .             alphaw(kar,ksxr), alpha(kdi,kar,ksxr)
      dimension    sen(kdi,kar), sex(kdi,kar)
c
c  parameters for old diode arrays (jsxr = 1)
c
      data nar(1) /3/
      data (name(k,1),k = 1,3) /'radial', 'vertup', 'vertdown'/
      data (ndiod(k,1),k = 1,3) /11, 10, 10/
      data (rw(k,1),k = 1,3) /253.0, 204.0, 204.0/
      data (zw(k,1),k = 1,3) /89.0, 89.0, 89.0/
      data (bw(k,1),k = 1,3) /21.6, 12.7, 12.7/
      data (alphaw(k,1),k = 1,3) /140.76, 180.0, 180.0/
      data (idiod(i,1,1),i = 1,11) /1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12/
      data (idiod(i,2,1),i = 1,10) /60, 61, 62, 63, 64, 65, 66, 67,
     .                            68, 69/
      data (idiod(i,3,1),i = 1,10) /60, 71, 72, 73, 74, 75, 76, 77,
     .                            78, 79/
      data (alpha(i,1,1),i = 1,11) /147.13, 146.29, 145.44, 144.59,
     .                            143.74, 142.89, 142.04, 141.19,
     .                            140.33, 139.48, 137.78/
      data (alpha(i,2,1),i = 1,10) /180.0, 176.0, 172.0, 168.0, 164.0,
     .                            160.0, 156.0, 152.0, 148.0, 144.0/
      data (alpha(i,3,1),i = 1,10) /180.0, 184.0, 188.0, 192.0, 196.0,
     .                            200.0, 204.0, 208.0, 212.0, 216.0/
c
c parameters for new diode arrays (jsxr = 2)
c
      data nar(2) /4/
      data (name(k,2),k = 1,4) /'sideup', 'sidedown', 'topout', 'topin'/
      data (ndiod(k,2),k = 1,4) /16, 16, 16, 16/
      data (rw(k,2),k = 1,4) /194.0, 194.0, 150.8, 146.8/
      data (zw(k,2),k = 1,4) / 91.0,  87.0, 150.7, 150.7/
      data (bw(k,2),k = 1,4) /120.0, 120.0, 120.0, 120.0/
      data (alphaw(k,2),k = 1,4) /140.0, 220.0, 310.0, 230.0/
      data (idiod(i,1,2),i = 1,16) / 1,  2,  3,  4,  5,  6,  7,  8,
     .                             9, 10, 11, 12, 13, 14, 15, 16/
      data (idiod(i,2,2),i = 1,16) /17, 18, 19, 20, 21, 22, 23, 24,
     .                            25, 26, 27, 28, 29, 30, 31, 32/
      data (idiod(i,3,2),i = 1,16) / 1,  2,  3,  4,  5,  6,  7,  8,
     .                             9, 10, 11, 12, 13, 14, 15, 16/
      data (idiod(i,4,2),i = 1,16) /17, 18, 19, 20, 21, 22, 23, 24,
     .                            25, 26, 27, 28, 29, 30, 31, 32/
      data (alpha(i,1,2),i = 1,16) /120.0, 124.0, 128.0, 132.0, 136.0,
     .                            140.0, 144.0, 148.0, 152.0, 156.0,
     .                            160.0, 164.0, 168.0, 172.0, 176.0,
     .                            180.0/
      data (alpha(i,2,2),i = 1,16) /180.0, 184.0, 188.0, 192.0, 196.0,
     .                            200.0, 204.0, 208.0, 212.0, 216.0,
     .                            220.0, 224.0, 228.0, 232.0, 236.0,
     .                            240.0/
      data (alpha(i,3,2),i = 1,16) /330.0, 326.0, 322.0, 318.0, 314.0,
     .                            310.0, 306.0, 302.0, 298.0, 294.0,
     .                            290.0, 286.0, 282.0, 278.0, 274.0,
     .                            270.0/
      data (alpha(i,4,2),i = 1,16) /270.0, 266.0, 262.0, 258.0, 254.0,
     .                            250.0, 246.0, 242.0, 238.0, 234.0,
     .                            230.0, 226.0, 222.0, 218.0, 214.0,
     .                            210.0/
c
      data nc     /51/
      data pio180 / 0.017453293/
      data plimo  /-1.0e10/
c
c initialize some parameters
c
      js = IABS (jsxr)
      if (js .ne. 1 .and. js .ne. 2)  return
      narray = nar(js)
      do 5 k=1,narray
        namar (k) = name (k,js)
        ndiode(k) = ndiod(k,js)
        do 5 i=1,ndiode(k)
    5     idiode(i,k) = idiod(i,k,js)
      if (jsxr .ne. 1 .and. jsxr .ne. 2)  return
c
c initialize plim
c
      if (codeid .eq. 'onedee')  plim = r(nj) / SQRT (kappa)
      if (codeid .ne. 'onedee')  plim = pmax
c
c  set coordinates to those of array window
c
      do 200 k=1,narray
      x0 = rw(k,jsxr)
      y0 = 0.0
      z0 = zw(k,jsxr)
c
c  calculate direction cosines for chord
c
      do 200 i=1,ndiode(k)
      if (namar(k) .ne. 'radial')  go to 10
      cx = COS (alpha(i,k,jsxr)*pio180)
      cy = SIN (alpha(i,k,jsxr)*pio180)
      cz = 0.0
      go to 20
   10 cx = COS (alpha(i,k,jsxr)*pio180)
      cy = 0.0
      cz = SIN (alpha(i,k,jsxr)*pio180)
   20 continue
c
c  skip calculation to determine chord endpoints if plim has not changed
c     since last call to sxrcal
c
      if (plim .eq. plimo)  go to 110
c
c  calculate distances along chord to enter and exit toroidal box
c     surrounding plasma
c
      call timtor(rin, rmax, x0, y0, z0, cx, cy, cz, zmin, zmax,
     .            senter, sexit)
c
c  skip remaining calculations if chord misses box
c
      sen(i,k) = -1.0e10
      if (senter .le. -1.0e10)  go to 110
c
c  calculate distances along chord to enter and exit plasma
c
      s1 = senter
      x1 = x0 + cx*s1
      y1 = y0 + cy*s1
      z1 = z0 + cz*s1
      r1 = SQRT (x1**2+y1**2)
      if (codeid .eq. 'onedee')
     .p1 = SQRT ((r1-rmajor)**2+(z1-zax)**2/kappa**2)
      if (codeid .ne. 'onedee')
     .call bilin(nx,ny,rmhdgrid,zmhdgrid,p,r1,z1,p1)
      ii = 1
      ds = (sexit-senter)/(nc-1.0)
      do 60 l=2,nc
      s2 = s1 + ds
      x2 = x1 + cx*ds
      y2 = y1 + cy*ds
      z2 = z1 + cz*ds
      r2 = SQRT (x2**2+y2**2)
      if (codeid .eq. 'onedee')
     .p2 = SQRT ((r2-rmajor)**2+(z2-zax)**2/kappa**2)
      if (codeid .ne. 'onedee')
     .call bilin(nx,ny,rmhdgrid,zmhdgrid,p,r2,z2,p2)
      if (ii .eq. 2)  go to 30
      if (p2 .gt. plim)  go to 50
      sen(i,k) = s1 + ds*(p1-plim)/(p1-p2)
      ii = 2
      go to 50
   30 if (p2 .lt. plim)  go to 50
      sex(i,k) = s1 + ds*(plim-p1)/(p2-p1)
      go to 110
   50 s1 = s2
      x1 = x2
      y1 = y2
      z1 = z2
      p1 = p2
   60 continue
      sex(i,k) = sexit
c
c  skip integral calculation if chord misses plasma
c
  110 roamin(i,k) = 1.0
      sxr(i,k) = 0.0
      if (sen(i,k) .le. -1.0e10)  go to 200
c
c  calculate window thickness along chord
c
      b = bw(k,jsxr)/COS ((alphaw(k,jsxr)-alpha(i,k,jsxr))*pio180)
c
c  evaluate line integral; omit check and correction for a
c     reentrant chord
c
      ds = (sex(i,k)-sen(i,k))/(nc-1.0)
      x1 = x0 + cx*sen(i,k)
      y1 = y0 + cy*sen(i,k)
      z1 = z0 + cz*sen(i,k)
      do 180 l=1,nc
      wt = 0.5
      rho1 = r(nj)
      if (l .eq. 1 .or. l .eq. nc)  go to 160
      wt = 1.0
      x1 = x1 + cx*ds
      y1 = y1 + cy*ds
      z1 = z1 + cz*ds
      r1 = SQRT (x1**2+y1**2)
      if (codeid .eq. 'onedee')
     .rho1 = SQRT (kappa*(r1-rmajor)**2+(z1-zax)**2/kappa)
      if (codeid .ne. 'onedee')
     .call bilin(nx,ny,rmhdgrid,zmhdgrid,p,r1,z1,p1)
      if (codeid .ne. 'onedee')
     .call interp(p1,psir,nj,r,rho1)
  160 call interp(rho1,r,nj,ene,ene1)
      call interp(rho1,r,nj,te,te1)
      roa1 = rho1/r(nj)
      roamin(i,k) = MIN (roa1, roamin(i,k))
  180 sxr(i,k) = sxr(i,k) + wt*ds*(1.0e-14*ene1)**2*sxremi(te1,b)
  200 continue
c
c  store plimo
c
      plimo = plim
      return
c
      end

      real*8 function sxremi (te, wt)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c          SXREMI         B. STOCKDALE, G. JAHNS       FEB. 1984
c
c          SXREMI RETURNS DIODE RESPONSE TO SXR EMISSION AS A FUNCTION
c          OF ELECTRON TEMPERATURE (TE, keV) AND EFFECTIVE
c          WINDOW THICKNESS (WT, 10**-6 M).
c
c          SXR EMISSION AT VARIOUS PLASMA TEMPERATURES WAS
c          CALCULATED BY G. JAHNS.  THIS CALCULATION PERFORMED
c          AN INTEGRAL OVER PHOTON ENERGY, AND INCLUDED THE
c          VARIATION OF THE DETECTOR RESPONSE AND GAUNT FACTOR.
c          A RATIONAL POLYNOMIAL WAS THEN FIT TO THE SXR EMISSION
c          IN LOG-LOG SPACE.  SINCE THE DETECTOR RESPONSE ALSO
c          DEPENDS UPON THE EFFECTIVE WINDOW THICKNESS (WHICH
c          VARIES WITH DIODE VIEWING ANGLE), A LINEAR
c          INTERPOLATION IS DONE IN THE LOG-LOG SPACE TO
c          OBTAIN THE RESPONSE FOR THE GIVEN EFFECTIVE THICKNESS.
c
c          THE MAXIMUM ABSOLUTE ERROR DUE TO THE INTERPOLATION IS
c          0.1% AT TE NEAR 10 keV, WHILE THE MAXIMUM RELATIVE
c          ERROR IS ABOUT 5% AT TE NEAR 0.1 keV.
c
      sxremi = 0.0
      if (TE .lt.  0.105)  return
      if (TE .gt. 10.0  )  return
      x = LOG10 (te)
c
c          SECOND GENERATION TOP AND SIDE DIODE ARRAYS
c          ON D-III (WHICH HAVE THE SAME WINDOW THICKNESS).
c
      if (WT .lt. 120. .or. WT .gt. 157.0)  go to 120
      WTN    = 120.0
      SXRN   = (-1.2969531 + X*(1.7384246 + X*8.4242456E-2)) /
     .  (1.0 + X*(0.95616248 + X*(0.41304208 + X*(0.16731969 +
     .       X*(9.0831717E-2 + X*3.0497805E-2)))))
      WTA    = 156.65
      SXRA   = (-1.3901854 + X*(1.4565813 + X*0.51460988)) /
     .  (1.0 + X*(1.21409800 + X*(0.65562881 + X*(0.27660306 +
     .       X*(1.3965734E-1 + X*4.5739938E-2)))))
      FRAC   = (WT-WTN)/(WTA-WTN)
      SXREMI = SXRN + FRAC*(SXRA-SXRN)
      SXREMI = 10.0**SXREMI
      return
c
c          FIRST GENERATION VERTICAL DIODE ARRAY ON D-III
c
  120 if (WT .lt. 12.7 .or. WT .gt. 15.8)  go to 140
      WTN = 12.70
      SXRN = (-0.74848697 + X*1.4349280) /
     .  (1.0 + X*(0.70337790 + X*(0.18035533 + X*(0.14465028 +
     .       X*8.5894026E-2))))
      WTA = 15.70
      SXRA = (-0.70523787 + X*1.3995740) /
     .  (1.0 + X*(0.68744098 + X*(0.16320590 + X*(0.13966804 +
     .       X*8.6690828E-2))))
      FRAC   = (WT-WTN)/(WTA-WTN)
      SXREMI = SXRN + FRAC*(SXRA-SXRN)
      SXREMI = 10.0**SXREMI
      return
c
c          FIRST GENERATION TANGENTIAL DIODE ARRAY ON D-III
c
  140 if (WT .lt. 21.6 .or. WT .gt. 21.8)  go to 160
      WTN    = 21.60
      SXRN   = (-0.81810804 + X*1.4897495) /
     .  (1.0 + X*(0.73077340 + X*(0.20949095 + X*(0.14974936 +
     .       X*8.2169405E-2))))
      WTA    = 21.74
      SXRA   = (-0.82044917 + X*1.4918113) /
     .  (1.0 + X*(0.73083423 + X*(0.21007007 + X*(0.15106074 +
     .       X*8.2875695E-2))))
      FRAC   = (WT-WTN)/(WTA-WTN)
      SXREMI = SXRN + FRAC*(SXRA-SXRN)
      SXREMI = 10.0**SXREMI
      return
c
c          INVALID WINDOW THICKNESS
c
  160 return
c
      end

      subroutine tspline (x, y, nx, bpar, cs, ic, ier, t, a, b, c,
     .                    fpp, r, dx, tmax, rgrid, tspl, npts)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c --- tspline calculates the second derivatives of the spline at the knots.
c --- these coefficients define the spline and are stored
c --- in array cs(i,j),i = 1,#knots,j=1,3
c --- subroutine EVTSPLN is used to evaluate the spline,
c
c --- input
c
c      x      vector,length nx,of knot locations (must be in ascending order)
c      y          function values at the knots
c      n          #points in x and y
c    rgrid(i)     i = 1,2...npts values of x at which evaluation is desired
c      bpar       array used to specifiy boundary conditions
c      ic         exact row dimension of matrix c
c      t          tension parameter
c
c      a,b,c,     temporary work vectors of length nx
c      fppr,dx
c
c --- output
c
c      cs         array of spline coefficients (see above)
c --- the value of the spline for x in the interval (x(i),x(i+1)) is
c ---           (the range on i is 1 to nx-1 )
c          f(x) = cs(i,1)*cs(i,2) * SINH (t*(x(i+1)-x))+
c               (y(i)*t*t-cs(i,1))*cs(i,3)*(x(i+1)-x)
c               +cs(i+1,1)*cs(i,2) * SINH (t*(x-x(i)))+
c               (y(i+1)*t*t-cs(i+1,1))*cs(j,3)*(x-x(i))
c       tmax   max allowed value of tension parameter
c       tspl(i)      i = 1,2..npts values of y at points  rgrid(i)
c
c ------------------------------------------------------------------ HSJ
c
      logical    reverseset
      dimension  a(*),b(*),c(*),fpp(*),r(*),dx(*),rgrid(*),tspl(*)
      dimension  x(*),y(*),bpar(*),cs(ic,3)
c
c --- check input
c
      reverseset = .false.
      dxmin = ABS (x(nx)-x(1))
      nm1 = nx-1
      ier = 0
      if (ic .lt. nx)  ier = 129
      if (nx .lt. 2)  ier = 130
      if (ier .ne. 0)  go to 1000
      do 15 j=1,nm1
        dxmin = MIN (dxmin, ABS (x(j+1)-x(j)))
        if (x(j) .lt. x(j+1))  go to 15
        ier = 131
        go to 1000
   15 continue
c
c --- some initialization
c
   16 bp1   = bpar(1)
      bp2   = bpar(2)
      bp3   = bpar(3)
      bp4   = bpar(4)
      tmax  = 30.0 / dxmin    ! max tension, avoids overflow
      t     = MIN (t, tmax)
      tsq   = t*t
      nm1   = nx-1
      dx(1) = x(2) - x(1)
c
c --- in case smtdx is incorrectly set
c
      if (t .eq. 0.0)  smtdx = 0.0
c
c --- set up vectors a,b,c
c --- vector a is lower diagonal
c --- vector b is diagonal
c --- vector c is upper diagonal
c
      if (t .eq. 0.0)  go to 50     ! zero tension case must be separate
      do 20 j=2,nm1
        dx(j) = x(j+1)-x(j)
        a(j) = (1.0/dx(j-1)-t / SINH (t*dx(j-1)))/tsq
        b(j) = t * COSH (t*dx(j-1)) / SINH (t*dx(j-1))
        b(j) = b(j)-1.0/dx(j-1)-1.0/dx(j)
        b(j) = (b(j)+t * COSH (t*dx(j)) / SINH (t*dx(j)))/tsq
        c(j) = (1.0/dx(j)-t / SINH (t*dx(j)))/tsq
   20   r(j) = (y(j+1)-y(j))/dx(j)-(y(j)-y(j-1))/dx(j-1)
c
c --- boundary conditions on spline
c
        th2 = t*dx(1)
        thn = t*dx(nm1)
        b1 = (1.0/t)*(1.0 / SINH (th2)-1.0/th2)
        a1 = -1.0/tanh(th2)+1.0/th2
        a1 = a1/(1.0 / SINH (th2)-1.0/th2)
        an = (1.0/t)*(-1.0 / SINH (thn)+1.0/thn)
        bn = 1.0 / tanh(thn)-1.0/thn
        bn = bn/(-1.0 / SINH (thn)+1.0/thn)
        if (bpar(1) .ne. -1.0e30)  go to 30
        bp1 = 1.0
****    bp2 = (-1.0/(dx(1)*b1))*(dx(1)*bpar(2)+y(1)-y(2))
        bp2 = (+1.0/(dx(1)*b1))*(dx(1)*bpar(2)+y(1)-y(2))
   30   if (bpar(3) .ne. -1.0e30)  go to 40
        bp3 = 1.0
        bp4 = (1.0/(dx(nm1)*an))*(dx(nm1)*bpar(4)+y(nm1)-y(nx))
   40   b(1) = a1
        c(1) = bp1
        a(nx) = bp3
        b(nx) = bn
        r(1) = bp2
        r(nx) = bp4
        go to 120
c
c --- for tension sufficiently small go over into cubic spline
c
   50 do 110 j=2,nm1
      dx(j) = x(j+1)-x(j)
      a(j) = dx(j-1)/6.0
      b(j) = (dx(j)+dx(j-1))/3.0
      c(j) = dx(j)/6.0
      if (t .eq. 0.0)  go to 110
      a(j) = a(j)-1.94444444e-02*tsq*dx(j-1)**3
      b(j) = b(j)-tsq*(dx(j-1)**3+dx(j)**3)/45.0
      c(j) = c(j)-1.94444444e-02*tsq*dx(j)**3
  110 r(j) = (y(j+1)-y(j))/dx(j)-(y(j)-y(j-1))/dx(j-1)
c
c --- boundary conditions
c
      if (bpar(1) .ne. -1.0e30)  go to 60
      bp1 = 1.0
      bcl = dx(1)*(-0.166666667+1.94444444e-02*dx(1)**2*tsq)
      bp2 = (1.0/(dx(1)*bcl))*(y(1)-y(2)+dx(1)*bpar(2))
   60 if (bpar(3) .ne. -1.0e30)  go to 70
      bp3 = 1.0
      ac = dx(nm1)*(0.166666667-1.94444444*dx(nm1)**2*tsq)
      bp4 = (1.0/(dx(nm1)*ac))*(dx(nm1)*bpar(4)+y(nm1)-y(nx))
   70 a(nx) = bp3
      b( 1) = 2.0 + dx(1  )**2*tsq / 10.0
      c( 1) = bp1
      b(nx) = 2.0 + dx(nm1)**2*tsq / 10.0
      r( 1) = bp2
      r(nx) = bp4
c
c --- note a(1) and c(nx) are not used
c --- r is vector of rhs
c
c --- solve the tridiagonal system
c --- forward elimination
c
  120 do j=2,nx
        pv = a(j)/b(j-1)
        b(j) = b(j)-c(j-1)*pv
        r(j) = r(j)-r(j-1)*pv
      end do
c
c --- back substitution
c
      fpp(nx) = r(nx)/b(nx)
      do j=1,nm1
        i = nx-j
        fpp(i) = (r(i)-c(i)*fpp(i+1))/b(i)
      end do
c
c --- now have vector of second derivatives fpp. set up the vector cs(i,j)
c --- to be used in evaluating the tension spline.
c --- convention used for natural spline will not work here.
c --- however we still pack all the necessary information into array cs.
c
      do 150 j=1,nx
        cs(j,1) = fpp(j)
        if (j .eq. nx )  go to 150
        if (t .eq. 0.0)  go to 150
        cs(j,2) = 1.0 / (SINH (t*dx(j))*t*t)
        cs(j,3) = 1.0 / (dx(j)*t*t)
  150 continue
c
      call evtspl (x, y, nx, cs, ic, rgrid, tspl, npts, t)
      if (.not. reverseset)  return
c
 1000 do j=1,nx/2
        xhold     = x(j)
        x(j)      = x(nx-j+1)
        x(nx-j+1) = xhold
        yhold     = y(j)
        y(j)      = y(nx-j+1)
        y(nx-j+1) = yhold
      end do
c
      reverseset  = .not. reverseset
      if (reverseset)  go to 16
      return
c
      end

      subroutine zanlyt2 (f, eps, nsig, kn, nguess, n, x, itmax,
     .                    infer, ier)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c   modified IMSL routine name - ZANLYT2
c
c ----------------------------------------------------------------------
c
c   computer            - dec10/single
c
c   latest revision     - january 1, 1978
c
c   purpose             - zeros of an analytic complex function
c                           using the muller method with deflation
c
c   usage               - call zanlyt2 (f,eps,nsig,kn,nguess,n,x,itmax,
c                                       infer,ier)
c
c   arguments    f      - a complex function subprogram, f(z), written
c                           by the user specifying the equation whose
c                           roots are to be found.  f must appear in
c                           an external statement in the calling pro-
c                           gram.
c                eps    - first stopping criterion.  let fp(z) = f(z)/p
c                           where p = (z-z(1))*(z-z(2))*,,,*(z-z(k-1))
c                           and z(1),...,z(k-1) are previously found
c                           roots.  if ((ABS (f(z)) .le. eps)  .and.
c                           (ABS (fp(z)) .le. eps)), then z is accepted
c                           as a root. (input)
c                nsig   - 2nd stopping criterion.  a root is accepted
c                           if two successive approximations to a given
c                           root agree in the first nsig digits. (input)
c                             note. if either or both of the stopping
c                             criteria are fulfilled, the root is
c                             accepted.
c                kn     - the number of known roots which must be stored
c                           in x(1),...,x(kn), prior to entry to ZANLYT2
c                nguess - the number of initial guesses provided. these
c                           guesses must be stored in x(kn+1),...,
c                           x(kn+nguess).  nguess must be set equal
c                           to zero if no guesses are provided. (input)
c                n      - the number of new roots to be found by
c                           ZANLYT2 (input)
c                x      - a complex vector of length kn+n.  x(1),...,
c                           x(kn) on input must contain any known
c                           roots.  x(kn+1),..., x(kn+n) on input may,
c                           on user option, contain initial guesses for
c                           the n new roots which are to be computed.
c                           if the user does not provide an initial
c                           guess, zero is used.
c                           on output, x(kn+1),...,x(kn+n) contain the
c                           approximate roots found by ZANLYT2.
c                itmax  - the maximum allowable number of iterations
c                           per root (input)
c                infer  - an integer vector of length kn+n.  on
c                           output infer(j) contains the number of
c                           iterations used in finding the j-th root
c                           when convergence was achieved.  if
c                           convergence was not obtained in itmax
c                           iterations, infer(j) will be greater than
c                           itmax (output).
c                ier    - error parameter (output)
c                         warning error
c                           ier = 33 indicates failure to converge with-
c                             in itmax iterations for at least one of
c                             the (n) new roots.
c
c   precision/hardware  - single and double/h32
c                       - single/h36,h48,h60
c
c   reqd. IMSL routines - uertst1,ugetio
c
c   notation            - information on special notation and
c                           conventions is available in the manual
c                           introduction or through IMSL routine uhelp
c
c   remarks      ZANLYT2 always returns the last approximation for root j
c                in x(j). if the convergence criterion is satisfied,
c                then infer(j) is less than or equal to itmax. if the
c                convergence criterion is not satisified, then infer(j)
c                is set to either itmax+1 or itmax+k, with k greater
c                than 1. infer(j) = itmax+1 indicates that ZANLYT2 did
c                not obtain convergence in the allowed number of iter-
c                ations. in this case, the user may wish to set itmax
c                to a larger value. infer(j) = itmax+k means that con-
c                vergence was obtained (on iteration k) for the defla-
c                ted function
c                              fp(z) = f(z)/((z-z(1)...(z-z(j-1)))
c
c                but failed for f(z). in this case, better initial
c                guesses might help or, it might be necessary to relax
c                the convergence criterion.
c
c   copyright           - 1978 by imsl, inc. all rights reserved.
c
c   warranty            - IMSL warrants only that IMSL testing has been
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c
c ----------------------------------------------------------------------
c
      external F
c
      dimension           x(*), infer(*)
      real*8              rzero,rten,rhun,rp01,ax,eps1,qz,eps,tpq
      complex*16          x,d,dd,den,fprt,frt,h,rt,t1,t2,t3,
     .                    tem,z0,z1,z2,bi,f,xx,xl,y0,y1,y2,x0,
     .                    zero,p1,one,four,p5
      data                zero/(0.0,0.0)/,p1/(0.1,0.0)/,
     .                    one/(1.0,0.0)/,four/(4.0,0.0)/,
     .                    p5/(0.5,0.0)/,
     .                    rzero/0.0/,rten/10.0/,rhun/100.0/,
     .                    ax/0.1/,ickmax/3/,rp01/0.01/
c
c                                  first executable statement
c
      ier = 0
      if (n .lt. 1)  go to 90
      eps1  = rten**(-nsig)
      eps1  = MIN (eps1,rp01)
c                                  set number of iterations
      knp1  = kn+1
      knpn  = kn+n
      knpng = kn+nguess
      do 5 i=1,knpn
         infer(i) = 0
         if (i .gt. knpng)  x(i) = zero
    5 continue
      l   = knp1
   10 jk  = 0
      ick = 0
      xl  = x(l)
   15 ic  = 0
      h   = ax
      h   = p1*h
      if (ABS (xl) .gt. ax)  h = p1*xl
c                                  first three points are
c                                    xl+h,  xl-h,  xl
      rt = xl+h
      assign 20 to nn
      go to 50
   20 z0 = fprt
      y0 = frt
      x0 = rt
      rt = xl-h
      assign 25 to nn
      go to 50
   25 z1 = fprt
      y1 = frt
      h = xl-rt
      d = h/(rt-x0)
      rt = xl
      assign 30 to nn
      go to 50
   30 z2 = fprt
      y2 = frt
c                                  begin main algorithm
   35 dd = one + d
      t1 = z0*d*d
      t2 = z1*dd*dd
      xx = z2*dd
      t3 = z2*d
      bi = t1-t2+xx+t3
      den = bi*bi-four*(xx*t1-t3*(t2-xx))
c                                  use denominator of maximum amplitude
      t1 = SQRT (den)
      qz = rhun * MAX (ABS (bi),ABS (t1))
      t2 = bi + t1
      tpq = ABS (t2)+qz
      if (tpq .eq. qz) t2 = zero
      t3 = bi - t1
      tpq = ABS (t3) + qz
      if (tpq .eq. qz) t3 = zero
      den = t2
      qz = ABS (t3)-ABS (t2)
      if (qz .gt. rzero) den = t3
c                                  test for zero denominator
      assign 30 to nn
      if (ABS (den) .eq. rzero)  go to 65
      d = -xx/den
      d = d+d
      h = d*h
      rt = rt + h
c                                  check convergence of the first kind
c
      if (ABS (h) .le. eps1 * MAX (ABS (rt),ax))  go to 70
      if (ic .ne. 0)  go to 15
      assign 40 to nn
      go to 50
   40 qz = ABS (fprt)-ABS (z2)*rten
      if (qz .ge. rzero)  go to 45
      z0 = z1
      z1 = z2
      z2 = fprt
      y0 = y1
      y1 = y2
      y2 = frt
      go to 35
c                  take remedial action to induce convergence
   45 continue
      d = d*p5
      h = h*p5
      rt = rt-h
   50 jk = jk+1
      if (jk .gt. itmax)  go to 75
      frt = f(rt)
      fprt = frt
c                                  test to see if first root is being
c                                     determined
      if (l .eq. 1)  go to 60
c                                  compute deflated function
      lm1 = l-1
      do 55 i=1,lm1
         tem = rt - x(i)
         if (ABS (tem) .eq. rzero)  go to 65
   55 fprt = fprt/tem
   60 continue
c                                  check convergence of the second kind
c
      if (ABS (fprt) .le. eps .and. ABS (frt) .le. eps)  go to 80
      go to nn,(20,25,30,40)
   65 continue
      if (ic .ne. 0)  go to 15
      tem = rten * eps1
      if (ABS (rt) .gt. ax)  tem = tem*rt
      rt  = rt+tem
      d   = (h + tem) * d / h
      h   =  h + tem
      go to 50
c                                  check solution
   70 continue
      if (ic .ne. 0)  go to 80
      ic = 1
      z0 = y1
      z1 = y2
      z2 = f(rt)
      xl = rt
      ick = ick+1
      if (ick .le. ickmax)  go to 35
c                                  warning error, itmax = maximum
      jk  = itmax + jk
   75 ier = 33
c                                  a root has been found
   80 x(l)     = rt
      infer(l) = jk
      l        = l + 1
      if (  l .le. knpn)  go to 10
      if (ier .eq. 0   )  go to 90
      call uertst1 (ier, 'zanlyt2')
c
   90 return
c
      end
