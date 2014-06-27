
      subroutine getbsmodel (jhirsh_out)
c
      USE param
      USE tfact
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c ----------------------------------------------------------------------
c     GETBSMODEL simply returns the value of jhirsh chosen by the user
c ----------------------------------------------------------------------
c
c
      jhirsh_out = jhirsh
      return
c
      end

      subroutine inicon
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** INICON loads the common block containing physical constants and  *
c ***  conversion constants.                                           *
c *** Last Revision: 1/95 W.A.Houlberg and S.E.Attenberger, ORNL.      *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
c      include 'param.i'
      include 'houlberg.i'
c
c ---------------------------------------------------------------------
c Physical constants
c ---------------------------------------------------------------------
c Velocity of light - (meter/second):
c
      clight = 2.9979e+08
c
c Elementary charge - (coulomb):
c
      coulom = 1.6022e-19
c
c Electron mass - (kilogram):
c
      elemas = 9.1095e-31
c
c Permittivity of free space - (Farad/meter):
c
      epsilo = 8.8419e-12
c
c Proton mass - (kilogram):
c
      promas = 1.6726e-27
c
c Permeability of free space - (Henry/meter):
c
      xmuo = 1.2566e-06
c
c Pi, our old favorite:
c
      pi = ACOS (-1.0)
c
c The probability of achieving equality under capitalism:
c
      zero = 0.0
c
c ---------------------------------------------------------------------
c Unit conversions.
c ---------------------------------------------------------------------
c
c Joules per keV:
c
      xj7kv = 1.6022e-16
      return
c
      end

      integer function ismax (n, sx, incx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** ISMAX is the index of the largest element of sx.                 *
c *** Last Revision: 1/95 W.A.Houlberg and S.E.Attenberger, ORNL.      *
c *** Input:                                                           *
c *** n-number of elements to be checked.                              *
c *** sx-vector to be checked.                                         *
c *** incx-increment in sx index.                                      *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
      real*8  sx(*)
c
      smax  = sx(1)
      ismax = 1
      ix    = 1
      do i=1,n
        if (sx(ix).gt.smax) then
          smax  = sx(ix)
          ismax = ix
        end if
        ix = ix + incx
      end do
      return
c
      end

      integer function ismin (n, sx, incx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** ISMIN is the index of the smallest element of sx.                *
c *** Last Revision: 1/95 W.A.Houlberg and S.E.Attenberger, ORNL.      *
c *** Input:                                                           *
c *** n-number of elements to be checked.                              *
c *** sx-vector to be checked.                                         *
c *** incx-increment in sx index.                                      *
c***********************************************************************
c
      real*8  sx(*)
c
      smin  = sx(1)
      ismin = 1
      ix    = 1
      do i=1,n
        if (sx(ix) .lt. smin) then
          smin  = sx(ix)
          ismin = ix
        end if
        ix = ix + incx
      end do
      return
c
      end
      subroutine nclboot (jhirsh, curbs,nclboot_fail)
c
c following nclboot was copied from G. Staebler 12/15/98 HSJ
c

c
c -------------------------------------------------------------------- *
c  NCLBOOT returns the bootstrap current calculated using the          *
c   Houlberg model (see reference below), which is comprised of a      *
c   set of necoclasical routines that are called. These routines       *
c   calculate the neoclassical flows using the reduced charge          *
c   state formalism of Hirshman, the friction coefficients of          *
c   Hirshman, and the velocity dependent viscosities of Shaing.        *
c                                                                      *
c                                                                      *
c   References:                                                        *
c      Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079.                   *
c      Kessel, Nucl Fusion 34 (1994) 1221.                             *
c      Shaing, Yokoyama, Wakatani, Hsu, to be published.               *
c      Houlberg, Shaing, Hirshman, to be published.                    *
c                                                                      *
c   Created:  6/95 by Daniel Finkenthal, GA                            *
c                                                                      *
c   Based on the code provided on 2/95 by:                             *
c             2/95 by W.A. Houlberg and S.E. Attenberger, ORNL         *
c                                                                      *
c   Comments:                                                          *
c       A set of metric quantities are passed in through the           *
c    houlberg.i common block/include file. These metrics have          *
c    been calculated in subroutine  NCLMETRICS.                        *
c       For now we assume that only a single (fully stripped)          *
c    charge state exists for each isotope, xi = 1.0. For better        *
c    coupling with the edge some sort of coronal equilibrium model     *
c    should be used to get a better representation of charge state     *
c    populations.                                                      *
c                                                                      *
c   Inputs:                                                            *
c       jhirsh - specifies which model options to use                  *
c              = 95  No fast ion contributions                         *
c              = 96  Include fast ion contributions                    *
c       Various metric quantities and plasma parameters passed in      *
c       through common blocks.                                         *
c                                                                      *
c   Internal Variables, passed through to Houlberg's routines:         *
c       mi-number of isotopes in test case.                            *
c       mz-highest charge state.                                       *
c       dencut-cutoff density (ignore densities below this)-/m**3.     *
c       amuai(a)-atomic mass number of isotope a, e first.             *
c       tai(a)-temperature of isotope a-keV.                           *
c       vtai(a)-thermal velocity of isotope a-m/s.                     *
c       denai(a,i)-density of isotope a, charge state i-/m**3.         *
c       xzi(a,i)-nZ**2/sum_i n*Z**2 for isotope a for xi factors.      *
c       xgrt(a,i)-d T(a,i)/ d rho-keV/m.                               *
c       xgrp(a,i)-d n(a,i)T(a,i)/d rho-keV/m**4.                       *
c       xb2-<B**2>-T**2.                                               *
c       xbm2-<B**-2>-/T**2.                                            *
c       xft-trapped fraction.                                          *
c       xngrth-n.grad(theta) (flux function)-radians/m.                *
c       xfhat-xmuo*F/psi'-1/m.                                         *
c       xgrbm2-<|grad rho|**2/B**2>-/T**2.                             *
c       fm(j)-poloidal expansion of geometric factor, j=1,3.           *
c                                                                      *
c -------------------------------------------------------------------- *
c                                                                      *
c   Changes (08-21-1995):                                              *
c       Added the fast ion pressure term for the beam and alpha        *
c    particles. This is done by including the beam and alpha           *
c    densities as seperate isotopes in the mi-long arrays (mi is       *
c    the number of isotopes.) The beam and alpha pressures are         *
c    derived from the  averaged  stored energy density (keV/cm3) of    *
c    the beam and alphas, ie:                                          *
c                    press = nkT = 0.666667*W                          *
c    (The 0.67 factor arises from Kinetic Theory.) The fast ion        *
c    densities enbeam and enalp are passed in through include files    *
c    and used to calculate effective temperatures from the stored      *
c    energy densities.                                                 *
c                              Daniel Finkenthal       08-21-1995      *
c                                                                      *
c -------------------------------------------------------------------- *
c                                                                      *
c   Changes (08-30-1995):                                              *
c       After consulting with W. Houlberg, it has been decided         *
c    that the above implementaion of the fast ion bootstrap is         *
c    incorrect since the viscosities and friction coefficients         *
c    will be markedly different for a non-maxwellian distribution.     *
c    Shaing et al believe the fast ions contribute little. However,    *
c    the thermal electrons due to the beam neutrals are important!     *
c    For now, we will neglect the beam bootstrap current entirely!     *
c                                                                      *
c                              Daniel Finkenthal       08-30-1995      *
c                                                                      *
c -------------------------------------------------------------------- *
c                                                                      *
c   Changes (09-08-1995):                                              *
c       Reintroduced fast ion bootstrap as an option (JHIRSH = 96).    *
c    Shaing now says that using the fast ion pressure is ok if energy  *
c    on same order as thermal ions. Beam probably ok, but really       *
c    pushing it for alphas!                                            *
c                                     Daniel Finkenthal    09-08-1995  *
c                                                                      *
c -------------------------------------------------------------------- *
c                                                                      *
c   Changes (8-1-1996):                                                *
c       Replaced the NCLASS subroutines with a new version from        *
c   Houlberg which includes a third moment and ion orbit squeezing.    *
c   An option to compute the radial electric field was introduced.     *
c   If cer_ion is not ' ' the Electric field and the poloidal and      *
c   toroidal rotation of the cer and main ions will be computed and    *
c   stored in cer.i for output to the netCDF file "iterdb.nc". The     *
c   electric field shear is also used to iterate the ion orbit         *
c   sqeezing factors. The default (ftcalc = 'analytic') trapped        *
c   particle fraction was set to the Lin-Liu-Miller formula. Setting   *
c   ftcalc = 'exact' obtains the integral result. The differences      *
c   in ft are small.                                                   *
c   The electron density for jhirsch=95 now includes the beam and      *
c   alpha particle contributions since these are thermal electrons.    *
c                        Gary Staebler                                 *
c ---------------------------------------------------------------------*
c
c ONETWO include files used to pass variables through common blocks
c

      USE param
      USE fusion
      USE ions
      USE soln
      USE mhdpar 
      USE nub2
      USE numbrs
      USE mesh
      USE machin
      USE geom
      USE tordlrot
      USE metrics
      USE neo2d
      USE cer
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c More terms for the beam and alpha components, Added 8-21-95 by DFF
c
c      include 'nub2.i'
c      include 'fusion.i'
c
c  CER velocities and magnetic components
c
c      include 'cer.i'
c      include 'tordlrot.i'
c
c  Metrics calculated in ONETWO to be used in calculation
c
c      include 'metrics.i'
c
c  physics and machine constants to be passed to boostrap routines
c
      include 'houlberg.i'
c
c  Variables used locally
c
      integer ima, k, j, itersqz
      real*8  curbs(kj), xmks, dvedr
      real*8  fcapa, rbp1, rbp2, rbpa, dpsidr
      real*8  rbta, udia_da(kj), Kpol_da(kj)
      real*8  udia_ca(kj), Kpol_ca(kj), Kpol_expa
      real*8  cer_bpa, cer_btdra, sqzmax, sqzmax_old, sqz_change
      real*8  cdia, Epsi_expa(kj), Epsia(kj)
      real*8  angrot_ca(kj), angrot_da(kj), ave_vpar_da(kj)
      real*8  sqzmin, sqz_da(kj), ugrta(kj)
      save    init
      data    init /0/
c
      mimz4     = 4 * mxmi * mxmz
      xmks      = 1.0e6
      inclprint = 0
      kncord    = 2
      sqzmin    = 1.0e-6
      nclboot_fail = 0
c
      if (init .eq. 0) then
c
c ***** Initialize fixed parameters.
c
c  Set constants
c
        call inicon
c
c  Set cutoff density
c
        dencut = 1.0e10
c
c  Set maximum charge state mz = 10 (up to Neon)
c
        mz = 10 ! should be set to mxmz HSJ
        mz = mxmz ! HSj 12/16/98
c
c  Set masses and charges. mith is number of THERMAL species!
c
        mith     = 1 + nion
        amuai(1) = elemas / promas
        do k = 2, mith
          amuai(k) = atw(k-1)
        end do
        init = 1
c
      end if
c
      if (inclprint .eq. 1) then
        write (6, *) '================================================'
        write (6, '(1x, 14(a10, 2x))')
     .               'Rho', 'xngrth', 'xfmsum', 'xb2', 'xbm2', 'xfhat',
     .               'xft', 'xgrbm2', 'denai(2,1)', 'tai(2)',
     .               'xgrp(2,1)','xgrt(2)', 'curbs'
      end if
c
c Determine if fast ions are to be included.
c Set IFASTION = 1 if JHIRSH = 96
c
      if (jhirsh .eq. 96)  ifastion = 1
c
c ----------------------------------------------------------------------
c Start big loop over centered-mesh
c ----------------------------------------------------------------------
c
      itersqz    = 0
      sqzmax     = 1.0
      sqzmax_old = 1.0
      call sscal (mi  , zero, tai      , 1)
      mimz       = mxmi * mz
      call sscal (mimz, zero, denai    , 1)
      call sscal (mimz, zero, xzi      , 1)
      call sscal (mi  , zero, xgrt     , 1)
      call sscal (mimz, zero, xgrp     , 1)
      call sscal (kj  , zero, Epsi     , 1)
      call sscal (kj  , zero, Epsia    , 1)
      call sscal (kj  , zero, Epsi_expa, 1)
      call sscal (kj  , zero, Kpol_da  , 1)
      call sscal (kj  , zero, Kpol_ca  , 1)
      call sscal (kj  , zero, angrot_da, 1)
      call sscal (kj  , zero, angrot_ca, 1)
      call sscal (kj  , zero, udia_ca  , 1)
      call sscal (kj  , zero, udia_da  , 1)
      call sscal (kj  , zero, ugrta    , 1)
  100 do j=1,nj-1
        drho = dr(j) / 100.0
c
c       Set metrics:
c
        xft    = 0.5 * (ftncl(j+1)+ftncl(j) )
        xb2    = 0.5 * (bsq(j+1) +bsq(j) )
        xbm2   = 0.5 * (bmsq(j+1) +bmsq(j) )
        xgrbm2 = 0.5 * (grbmsq(j+1) +grbmsq(j) )
        xngrth = 0.5 * (grth(j+1) + grth(j) )
        do m=1,3
          fm(m) = gfm(m,j)
        end do
c
c       Calculate xfhat - NOTE that units cancel:
c
        fcapa  = 0.5 * (fcap(j+1) +fcap(j))
        rbp1   = rbp(j)/(fcap(j)*gcap(j)*hcap(j))
        rbp2   = rbp(j+1)/(fcap(j+1)*gcap(j+1)*hcap(j+1))
        rbpa   = 0.5 * (rbp1 + rbp2)
        dpsidr = (rbpa / ra(j)) * rmajor
        xfhat  = btor * rmajor / (fcapa * dpsidr)
c
c      Calculate normalized Electric Field derivative at mesh center
c
        rbta  = btor*rmajor/(fcapa*xmks)
        dvedr = (promas*rbta*xfhat/(coulom*(btor*1.0e-4)**2))*
     .          (Epsi(j+1) - Epsi(j))/drho
c
c **** Calculate some dependent variables at mesh centers
c
c   First Electrons:
c      (Determine from ion densities)
c
        ene0       = 0.0
        ene1       = 0.0
        do k=1,nion
          ene0     = ene0 + en(j,k)  * z(j,k)
          ene1     = ene1 + en(j+1,k)* z(j+1,k)
        end do
        ene0       = ene0 + 1.0 * enbeam(j) + 2.0 * enalp(j)
        ene1       = ene1 + 1.0 * enbeam(j) + 2.0 * enalp(j)
        denai(1,1) = xmks*(ene1 + ene0) / 2.0
        tai  (1  ) = 0.5 * (te(j+1) + te (j))
        xgrp (1,1) = xmks*(ene1 * te(j+1) - ene0 * te(j))/drho
        xgrt (1  ) = (te(j+1) - te(j)) / drho
        xzi  (1,1) = 1.0
        sqz  (1,1) = 1.0 + amuai(1)*dvedr
        if (sqz(1,1) .lt. sqzmin)
     .  sqz  (1,1) = sqzmin
        sqz  (1,1) = (sqz(1,1))**1.5
        sqzmax     =  MAX (sqzmax, sqz(1,1))
c
c   Ions and Impurities (k=1,..NION):
c     For now we assume that only a single (fully stripped)
c     charge state exists for each isotope, xzi=1.0
c
         
        do k=1,nion
          ima           = k+1
          ia            = ANINT (z(j,k))
          denai(ima,ia) = 0.5 * xmks*(en(j+1,k)+en(j,k))
          tai(ima)      = 0.5 * (ti(j+1) + ti(j))
          xgrt(ima)     = (ti(j+1)-ti(j))/drho
          xgrp(ima,ia)  = xmks*(en(j+1,k)*ti(j+1)-en(j,k)*ti(j))/drho
          xzi(ima,ia)   = 1.0
          sqz(ima,ia)   = 1.0 - amuai(ima)*dvedr/z(j,k)
          if (sqz(ima,ia) .lt. sqzmin)
     .    sqz(ima,ia)   = sqzmin
          if (k .eq. 1)  sqz_da(j) = sqz(ima,ia)
          sqz(ima,ia)   = (sqz(ima,ia))**1.5
          sqzmax        =  MAX (sqzmax, sqz(ima,ia))
        end do


c
c   Set the number of ion species to the number of thermal species
c       (May be incremented if beam and alpha ions are present)
c
        mi = mith
c
c ------------ FAST ION CONTRIBUTIONS ----------------------------------
c   If fastion = 1 then Add Beam and Alpha particles
c   For now assume an effective temperature derived from the  averaged
c   stored energy density of the beam and alphas, ie:
c                    press = nkT = 0.666667*W
c   The (2/3) factor arises from Kinetic Theory.
c   We have to be careful since enbeam or enalp may be zero at some points
c   on the grid. Therefore, only include these ion densities in the mi-long
c   arrays if greater than dencut (Houlberg's routines will crash for
c   zero denisity species).
c
        if (ifastion .eq. 1) then
c
c         Beam Particles:
c
          enbeama = 0.5 * xmks*(enbeam(j+1)+enbeam(j))
          if (enbeama .ge. dencut) then
            if (enbeam(j) .ne. 0) then
              tbeamj0     = 0.666667 * wbeam(j)  /enbeam(j)
              tbsave      = MIN (tbeamj0, tbsave)
            else
              tbeamj0     = tbsave
            end if
            if (enbeam(j+1) .ne. 0) then
              tbeamj1     = 0.666667 * wbeam(j+1)/enbeam(j+1)
            else
              tbeamj1     = tbsave
            end if
            mi  = mi + 1
            ima = mi
            ia  = 1
            amuai(ima)    = 2.0
            denai(ima,ia) = enbeama
            tai  (ima)    = 0.5 * (tbeamj0 + tbeamj1)
            xgrt (ima)    = (tbeamj1 - tbeamj0)/drho
            xgrp (ima,ia) = xmks*0.666667*(wbeam(j+1)-wbeam(j))/drho
            xzi  (ima,ia) = 1.0
            sqz  (ima,ia) = 1.0 - amuai(ima)*dvedr
            if (sqz(ima,ia) .lt. sqzmin)
     .      sqz  (ima,ia) =  sqzmin
            sqz  (ima,ia) = (sqz(ima,ia))**1.5
            sqzmax        = MAX (sqzmax, sqz(ima,ia))
          end if
c
c         Alpha Particles:
c
          enalpa = 0.5 * xmks*(enalp(j+1)+enalp(j))
          if (enalpa .ge. dencut) then
            if (enalp(j) .ne. 0) then
              talpj0      = 0.666667 * walp(j)  /enalp(j)
              talpsave    = MIN (talpj0, talpsave)
            else
              talpj0      = talpsave
            end if
            if (enalp(j+1) .ne. 0) then
              talpj1      = 0.666667 * walp(j+1)/enalp(j+1)
            else
              talpj1      = talpsave
            end if
            mi            = mi + 1
            ima           = mi
            ia            = 2
            amuai(ima)    = 4.0
            denai(ima,ia) = enalpa
            tai  (ima)    = 0.5 * (talpj0 + talpj1)
            xgrt (ima)    = (talpj1 - talpj0)/drho
            xgrp (ima,ia) = xmks*0.666667*(walp(j+1)-walp(j))/drho
            xzi  (ima,ia) = 1.0
            sqz  (ima,ia) = 1.0 - amuai(ima)*dvedr/2.0
            if (sqz(ima,ia) .lt. sqzmin)
     .      sqz  (ima,ia) =  sqzmin
            sqz  (ima,ia) = (sqz(ima,ia))**1.5
            sqzmax        = MAX (sqzmax, sqz(ima,ia))
          end if
c
        end if

c
c -------------- end of fast ion contributions -------------------------
c
c  Thermal velocities-m/s:
c
        do ima=1,mi
          vtai(ima) = SQRT (2.0*xj7kv*tai(ima)/(amuai(ima)*promas))
        end do
c
c       Set unity <E.B>:
c
        xedotb = 1.0
c
c ----------------------------------------------------------------------
c Now call the Houlberg neoclassical routines in sequence
c to get bootstrap current at this point on the half-grid
c ----------------------------------------------------------------------
c
c  Calculate M's and N's for each isotope:
c

        call nclmn
c
c  Evaluate tauab and averaged m*n/tauab:
c

        call ncltau
c
c  Calculate calm and caln:
c

        call nclcmn
c
c  Calculate viscosities:
c

        call nclmu
c
c  Calculate flows:
c

        call nclfl(nclboot_fail)
        if(nclboot_fail .eq. 1)then
           write(6,FMT ='(" Unknown failure in NCLBOOT")')
           return
        endif

c
c -------------------------------------------------------------------
c
c       Transform the Houlberg bootstrap <J_bs dot B>     (A-T/m2)
c       into the form used by ONETWO:
c       Convert  A-T/m2 -> A-T/cm2:
c
        xjdotb = xjdotb * 1.0e-4
c
c       Divide by Bt0 to get bootstrap current in form used by ONETWO:
c       <J_bs dot B/Bt0> (A/cm2). Note that Btor is in Gauss!
c
        curbs(j) = xjdotb / (btor*1.0e-4)
        fmsum    = fm(1) + fm(2) + fm(3)
c
        if (inclprint .eq. 1) then
          write (6, '(1x, 14(e10.3, 2x))')
     .                ra(j), xngrth, fmsum, xb2, xbm2, xfhat, xft,
     .                xgrbm2, denai(2,1), tai(2), xgrp(2,1), xgrt(2),
     .                curbs(j)
        end if
c
        if (j_cer .ne. 0 ) then
c
c         calculate CER velocities and flux-averaged velocities
c         neglecting the <E*B> contribution to <U*B>
c
          cdia          = (xfhat*xj7kv) / (rbta*coulom)
          cer_bpa       = 0.5*(cer_bp(j+1)   + cer_bp(j))
          cer_btdra     = 0.5*(cer_btdr(j+1) + cer_btdr(j))
          angrot_ca(j)  = 0.5*(angrot(j+1)   + angrot(j))
          Kpol_expa     = 0.5*(Kpol_exp(j+1) + Kpol_exp(j))
          ugrta(j)      = cdia*xgrt(2)
c
c         calculate Epsi from CER data and ave_upol_c
c
          k            = j_cer
          ima          = k + 1
          ia           = AINT (z(j,k))
          udia_ca(j)   = 0.0
          if (denai(ima,ia) .gt. dencut)
     .    udia_ca(j)   = cdia*xgrp(ima,ia)/(z(j,k)*denai(ima,ia))
          Kpol_ca(j)   = (uai(1,1,ima,ia) + rbta*udia_ca(j))/xb2
          Epsia(j)     = angrot_ca(j) +
     .                   udia_ca(j) - cer_btdra*Kpol_ca(j)
          Epsi_expa(j) = angrot_ca(j) +
     .                   udia_ca(j) - cer_btdra*Kpol_expa
c
c         calculate main ion velocities
c
          k              = 1
          ima            = k + 1
          ia             = AINT (z(j,k))
          udia_da(j)     = 0.0
          if (denai(ima,ia) .gt. dencut)
     .    udia_da(j)     = cdia*xgrp(ima,ia)/(z(j,k)*denai(ima,ia))
          Kpol_da(j)     = (uai(1,1,ima,ia) + rbta*udia_da(j))/xb2
          angrot_da(j)   = Epsia(j) - udia_da(j) + cer_btdra*Kpol_da(j)
          ave_vpar_da(j) = Kpol_da(j)*xb2 -
     .                     rbta*(udia_da(j) - Epsia(j))
          ave_vpar_da(j) = ave_vpar_da(j) / SQRT (xb2)
        end if
      end do
c
c *** end of loop over mesh ***
c
      if (j_cer .ne. 0) then
c
c  interpolate Epsia onto full mesh
c
         call intrp (1, 1, ra, Epsia, nj-1, r, Epsi, nj)
c
c  check for convergence of squeeze factor
c
         if (squeeze) then
           itersqz = itersqz + 1
           if (itersqz .eq. 1)  go to 100
           sqz_change = ABS (sqzmax - sqzmax_old) / sqzmax_old
           sqzmax_old = sqzmax
           sqzmax     = 1.0
           if (itersqz .lt. 20 .and. sqz_change .gt. 0.01)  go to 100
           write (6, *) 'ION ORBIT SQUEEZING USED IN NCLBOOT'
           write (6, *) 'itersqz = ',itersqz, '  sqzmax = ',
     .                   sqzmax_old, '  sqz_change = ', sqz_change
         end if
         call intrp (1,1,ra,Kpol_ca,nj-1,r,Kpol_c,nj)
         call intrp (1,1,ra,Kpol_da,nj-1,r,Kpol_d,nj)
         call intrp (1,1,ra,angrot_ca,nj-1,r,angrot_c,nj)
         call intrp (1,1,ra,angrot_da,nj-1,r,angrot_d,nj)
         call intrp (1,1,ra,ave_vpar_da,nj-1,r,ave_vpar_d,nj)
         call intrp (1,1,ra,sqz_da,nj-1,r,sqz_d,nj)
         call intrp (1,1,ra,Epsi_expa,nj-1,r,Epsi_exp,nj)
         call intrp (1,1,ra,udia_da,nj-1,r,udia_d,nj)
c
c        the following is in error I think HSJ 1/31/99:
c
****     call intrp (1,1,ra,udai_ca,nj-1,r,udia_c,nj)
c
c        replaced with:
c
         call intrp (1,1,ra,udia_ca,nj-1,r,udia_c,nj)
         call intrp (1,1,ra,ugrta,nj-1,r,ugrt,nj)
      end if

      return
c
      end

      subroutine nclcmn
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** NCLCMN calculates the neoclassical friction coefficients.        *
c *** References:                                                      *
c *** Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079.                    *
c *** Houlberg, Shaing, Hirshman, to be published.                     *
c *** Last Revision: 5/96 W.A.Houlberg, ORNL.                          *
c***********************************************************************
c
c      include 'param.i'
      include 'houlberg.i'
c
c *** Zero out calm.
c
      mi9 = 9 * mxmi
      call sscal(mi9,zero,calm,1)
c
c *** Calculate calm and caln.
c
      do     ima=1,mi
        do   imb=1,mi
          do   j=1,3
            do k=1,3
c
c ***         Sum over isotopes b for test particle component
c
              calm(j,k,ima)=calm(j,k,ima)
     .                      +amntau(ima,imb)*capm(j,k,ima,imb)
c
c ***         Field particle component
c
              caln(j,k,ima,imb)=amntau(ima,imb)*capn(j,k,ima,imb)
            end do
          end do
        end do
      end do
      return
c
      end

      subroutine nclfl(nclboot_fail)
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** NCLFL calculates the neoclassical flows u_ai and q_ai/p_ai.      *
c *** References:                                                      *
c *** Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079.                    *
c *** Houlberg, Shaing, Hirshman, to be published.                     *
c *** Last Revision: 5/96 W.A.Houlberg, ORNL.                          *
c***********************************************************************
c
      parameter      (mxsp=16)
c
c      include 'param.i'
      include 'houlberg.i'
c
      dimension      aa(3,3),         rhat(3,6,mxsp),  crhat(3,6,mxmi),
     .               srcth(3,mxsp),   xl(3),
     .               jzi(mxsp),       jmi(mxsp),       indx(3)
      dimension      ab(3*mxmi,3*mxmi)
      dimension      xab(3*mxmi,3),   indxb(3*mxmi)
c
c *** Initialize isotopic responses.
c
      mi18 = 3 * 6 * mi
      call sscal (mi18, zero, crhat, 1)
c
c *** Calculate responses for each species.
c
      cc=(xfhat/coulom)*xj7kv
      india=0
      do ima=1,mi
        do ia=1,mz
          if (denai(ima,ia).gt.dencut) then
            india=india+1
c
c ***       Set isotope number and charge state for this species.
c
            jmi(india)=ima
            jzi(india)=ia
c
c ***       Set up response matrix for each charge state.
c
            do i=1,kncord
                do j=1,kncord
                aa(i,j)=xzi(ima,ia)*calm(i,j,ima)-xmu(i,j,ima,ia)
                end do
              end do
c
c ***       Get lu decomposition of response matrix.
c
            iflag=0
            call uldcmp(aa,kncord,3,indx,d,iflag)
c
c ***       Get sources and evaluate responses from back substitution.
c ***       Lambda terms involving isotopic flows.
c
            do j=1,kncord
                do i=1,kncord
                  if (i.eq.j) then
                    rhat(i,j,india)=xzi(ima,ia)
                  else
                    rhat(i,j,india)=0.0
                  end if
                end do
              call ulbksb(aa,kncord,3,indx,rhat(1,j,india))
                do i=1,kncord
                crhat(i,j,ima)=crhat(i,j,ima)+
     .             xzi(ima,ia)*rhat(i,j,india)
                end do
              end do
c
c ***       Poloidal source (p' and T') terms.
c
            if (amuai(ima) .lt. 0.5) then
              zchga=-ia
            else
              zchga=ia
            end if
            srcth(1,india)=(cc/zchga)*xgrp(ima,ia)/denai(ima,ia)
              srcth(2,india)=(cc/zchga)*xgrt(ima)
              srcth(3,india)=0.0
              do i=1,kncord
                rhat(i,4,india)=0.0
                do j=1,kncord
                  rhat(i,4,india)=rhat(i,4,india)
     .                          +srcth(j,india)*xmu(i,j,ima,ia)
                end do
              end do
            call ulbksb(aa,kncord,3,indx,rhat(1,4,india))
            do i=1,kncord
              crhat(i,4,ima)=crhat(i,4,ima)+xzi(ima,ia)*rhat(i,4,india)
            end do
c
c ***       Parallel electric field terms for resistivity.
c
            rhat(1,5,india)=-zchga*coulom*denai(ima,ia)
            rhat(2,5,india)=0.0
            rhat(3,5,india)=0.0
            call ulbksb(aa,kncord,3,indx,rhat(1,5,india))
              do i=1,kncord
              crhat(i,5,ima)=crhat(i,5,ima)
     .                       +xzi(ima,ia)*xedotb*rhat(i,5,india)
              end do
          end if
        end do
      end do
c
c *** Evaluate parallel electrical resistivity.
c
      sumcon=0.0
      do inda=1,india
        jma=jmi(inda)
        jza=jzi(inda)
        if (amuai(jma).lt.0.5) then
          zchga=-jza
        else
          zchga=jza
        end if
c
c ***   Parallel electrical conductivity.
c
        sumcon=sumcon+zchga*coulom*denai(jma,jza)*rhat(1,5,inda)
      end do
c
c *** Parallel electrical resistivity.
c
      etap=1.0/sumcon
c
c *** Initialize coefficient matrix for isotopic flows.
c
      mxmi3  = 3 * mxmi
      mi3    = 3 * mi
      mxmi33 = mxmi3 * mi3
      call sscal(mxmi33,zero,ab,1)
c
c *** Load coefficient matrix and source terms for isotopic flows.
c
      do ima=1,mi
          do ia=1,kncord
            ia1=ima+(ia-1)*mi
c
c ***     Diagonal coefficients.
c
          ab(ia1,ia1)=1.0
c
c ***     Source terms.
c ***     Poloidal source (p' and T') terms.
c
          xab(ia1,1)=crhat(ia,4,ima)
c
c ***     <E.B> terms.
c
          xab(ia1,2)=crhat(ia,5,ima)
c
c ***     Field particle friction terms.
c
          do imb=1,mi
              do ib=1,kncord
                ib1=imb+(ib-1)*mi
                do k=1,kncord
                ab(ia1,ib1)=ab(ia1,ib1)
     .                      +caln(k,ib,ima,imb)*crhat(ia,k,ima)
                end do
              end do
            end do
        end do
      end do
c
c *** Get lu decomposition of coefficient matrix.
c
      iflag=0
      mma=kncord*mi
      mmad=3*mxmi
      call uldcmp(ab,mma,mmad,indxb,d,iflag)
      if(iflag .eq. 1) then
         nclboot_fail =1 
       return
      endif
c
c *** Evaluate isotopic flows from back substitution for sources,
c *** then for <E.B> terms.
c *** xab(1,k) to xab(mi,k) are the isotopic velocities.
c *** xab(mi+1,k) to xab(2*mi,k) are the isotopic heat flows.
c *** xab(2*mi+1,k) to xab(3*mi,k) are the u2 flows.
c *** Initialize species flows.
c
      mimz6 = 3 * 2 * mxmi * mz
      call sscal(mimz6,zero,uai,1)
c
c *** Evaluate species flows.
c
      do k=1,2
        call ulbksb(ab,mma,mmad,indxb,xab(1,k))
      end do
      do inda=1,india
        jma=jmi(inda)
        jza=jzi(inda)
        if (amuai(jma).lt.0.5) then
          zchga=-jza
        else
          zchga=jza
        end if
c
c ***     Source contributions.
c
        do k=1,2
            do ia=1,kncord
            if (k.eq.1) then
              uai(ia,k,jma,jza)=rhat(ia,4,inda)
            else
              uai(ia,k,jma,jza)=xedotb*rhat(ia,5,inda)
            end if
            end do
c
c ***     Response contributions.
c
          do imb=1,mi
              call sscal (kncord, zero, xl, 1)
              do ib=1,kncord
                ib1=imb+(ib-1)*mi
                do ia=1,kncord
                xl(ia)=xl(ia)-caln(ia,ib,jma,imb)*xab(ib1,k)
                end do
              end do
              do ib=1,kncord
                do ia=1,kncord
                uai(ia,k,jma,jza)=uai(ia,k,jma,jza)
     .                            +xl(ib)*rhat(ia,ib,inda)
                end do
              end do
          end do
        end do
      end do
c
c *** Calculate bootstrap current and fluxes.
c
      xjdotb=0.0
      mimz4=mxmi*mxmz*4
      call sscal(mimz4,zero,gfl,1)
      call sscal(mimz4,zero,qfl,1)
      cbp=xfhat/xb2/coulom
      cps=(xfhat/coulom)*(1.0/xb2-xbm2)
      ccl=(xgrbm2/coulom)/xfhat
      ceb=xfhat/xb2*xedotb
      do inda=1,india
        jma=jmi(inda)
        jza=jzi(inda)
        if (amuai(jma).lt.0.5) then
          zchga=-jza
        else
          zchga=jza
        end if
c
c ***   <J.B>.
c
        xjdotb=xjdotb+denai(jma,jza)*zchga*coulom*uai(1,1,jma,jza)
c
c ***   Banana-Plateau.
c
        cbpa=cbp/zchga
          cbpaq=cbpa*(xj7kv*tai(jma))
          do ib=1,kncord
c
c ***     Pressure and temperature grasient contributions.
c
          upapt=uai(ib,1,jma,jza)+srcth(ib,inda)
          gfl(1,jma,jza)=gfl(1,jma,jza)
     .                   -cbpa*xmu(1,ib,jma,jza)*upapt
          qfl(1,jma,jza)=qfl(1,jma,jza)
     .                   -cbpaq*xmu(2,ib,jma,jza)*upapt
c
c ***     <E.B> contributions.
c
          gfl(4,jma,jza)=gfl(4,jma,jza)
     .                   -cbpa*xmu(1,ib,jma,jza)*uai(ib,2,jma,jza)
          qfl(4,jma,jza)=qfl(4,jma,jza)
     .                   -cbpaq*xmu(2,ib,jma,jza)*uai(ib,2,jma,jza)
          end do
c
c ***   Test particle contributions.
c
        cpsa=cps*(xzi(jma,jza)/zchga)
          cpsaq=cpsa*(xj7kv*tai(jma))
        ccla=ccl*(xzi(jma,jza)/zchga)
          cclaq=ccla*(xj7kv*tai(jma))
          do ib=1,kncord
c
c ***     Pfirsch-Schluter.
c
          gfl(2,jma,jza)=gfl(2,jma,jza)
     .                   -cpsa*calm(1,ib,jma)*srcth(ib,inda)
          qfl(2,jma,jza)=qfl(2,jma,jza)
     .                   -cpsaq*calm(2,ib,jma)*srcth(ib,inda)
c
c ***     Classical.
c
          gfl(3,jma,jza)=gfl(3,jma,jza)
     .                   +ccla*calm(1,ib,jma)*srcth(ib,inda)
          qfl(3,jma,jza)=qfl(3,jma,jza)
     .                   +cclaq*calm(2,ib,jma)*srcth(ib,inda)
          end do
c
c ***   Field particle contributions.
c
        do indb=1,india
          jmb=jmi(indb)
          jzb=jzi(indb)
          if (amuai(jmb).lt.0.5) then
            zchgb=-jzb
          else
            zchgb=jzb
          end if
          cpsb=cpsa*xzi(jmb,jzb)
            cpsbq=cpsb*(xj7kv*tai(jma))
          cclb=ccla*xzi(jmb,jzb)
            cclbq=cclb*(xj7kv*tai(jma))
            do ib=1,kncord
c
c ***       Pfirsch-Schluter.
c
            gfl(2,jma,jza)=gfl(2,jma,jza)
     .                     -cpsb*caln(1,ib,jma,jmb)*srcth(ib,indb)
            qfl(2,jma,jza)=qfl(2,jma,jza)
     .                     -cpsbq*caln(2,ib,jma,jmb)*srcth(ib,indb)
c
c ***       Classical.
c
            gfl(3,jma,jza)=gfl(3,jma,jza)
     .                     +cclb*caln(1,ib,jma,jmb)*srcth(ib,indb)
            qfl(3,jma,jza)=qfl(3,jma,jza)
     .                     +cclbq*caln(2,ib,jma,jmb)*srcth(ib,indb)
            end do
        end do
      end do
      return
c
      end

      subroutine nclk (xx)
c
c***********************************************************************
c *** NCLK calculates the velocity-dependent neoclassical viscosity    *
c ***  coefficients, K, for the banana and Pfirsch-Schluter regimes for*
c ***  charge state i of isotope a.                                    *
c *** Input:                                                           *
c *** xx-v/vth, normalized velocity.                                   *
c *** References:                                                      *
c *** Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079.                    *
c *** Kessel, Nucl Fusion 34 (1994) 1221.                              *
c *** Shaing, Yokoyama, Wakatani, Hsu, Phys Plasmas 3 (1996) 965.      *
c *** Houlberg, Shaing, Hirshman, to be published.                     *
c *** Last Revision: 5/96 W.A.Houlberg, ORNL.                          *
c *** Comments:                                                        *
c *** The collisional modifications to the banana contributions are    *
c ***  lumped in with the banana contribution for the Hirshman, Sigmar *
c ***  and Kessel models.                                              *
c *** The array fm(k) and xngrth require different approximations to   *
c ***  the magnetic geometry in the formulations used here:            *
c *** Shaing, et al., kboot.lt.1, .gt.4:                               *
c *** fm(k)=2<sin(m.theta) n.gr(B)><sin(m.theta) n.gr(B) B.gr(theta)>  *
c ***        /<B**2><B.gr(theta)>, for system with n.gr(theta)=f(psi). *
c *** xngrth=n.gr(theta)=f(psi).                                       *
c *** Hirshman, Sigmar, kboot=1,3 (3 includes 1/fc correction):        *
c *** fm(k)=L_c^*(k), k=1,3, Eqn 4.66.                                 *
c *** xngrth-<(n.gr(B))**2> / <B**2> for Eqn 4.74.                     *
c *** Kessel, kboot=2:                                                 *
c *** fm(1)=Rq/eps**3/2 for Eqn 49.                                    *
c *** fm(2)=eps**3/2.                                                  *
c *** xngrth-not used.                                                 *
c *** Banana regime only, kboot=4:                                     *
c *** fm-not used.                                                     *
c *** xngrth-not used.                                                 *
c***********************************************************************
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
      include 'houlberg.i'
c
      mimz=mxmi*mz
      call sscal(mimz,zero,xkps,1)
      call sscal(mimz,zero,xkban,1)

      call nclnu(xx)

      do ima=1,mi
        if (kboot.eq.1.or.kboot.eq.3) then
c
c ***     H&S formulation.
c
          elcstar=0.0
          do mp=1,3
            elcstar=elcstar+fm(mp)
          end do
          omegat=vtai(ima)/elcstar
          c1=8.0/(3.0*pi)*xft*omegat/vtai(ima)**2/xx/xngrth
          c2=5.0*pi/8.0/xx/omegat
        else if (kboot .eq. 2) then
c
c ***     Kessel formulation.
c
          c1=2.48*fm(1)/vtai(ima)/xx
          c2=1.96*fm(1)*fm(2)/vtai(ima)/xx
        else
c
c ***     Shaing, et al., formulation.
c
          c1=1.5*vtai(ima)**2*xx**2
        end if
        do ia=1,mz
          if (denai(ima,ia).gt.dencut) then
            if (kboot.gt.0.and.kboot.le.3) then
c
c ***         H&S, Kessel formulations.
c
              xkban(ima,ia)=xft*xnud(ima,ia)/(1.0+c1*xnud(ima,ia))
     .                                      /(1.0+c2*xnut(ima,ia))
              if (kboot.eq.2.or.kboot.eq.3) then
c
c ***           H&S, Kessel formulations with circulating fraction.
c
                xkban(ima,ia)=xkban(ima,ia)/(1.0-xft)
              end if
            else
c
c ***         Shaing, et al., formulaton.
c
              xkban(ima,ia)=xft/(1.0-xft)*xnud(ima,ia)
              c2=c1/xnut(ima,ia)
              do mp=1,3
                xkps(ima,ia)=xkps(ima,ia)+c2*fm(mp)*xnuti(mp,ima,ia)
              end do
            end if
          end if
        end do
      end do

      return
c
      end

      subroutine nclmetrics (xmagax, ymagax, r0, mp, xp0, yp0, bp0,
     .                       xeps, fpsi, xqs, ffprime, x, nx, y, ny,
     .                       cspln, xngrth, xb2, xbm2, xgrbm2, fm)
c
c
c ----------------------------------------------------------------------
c  NCLMETRICS sets up the metrics for the Houlberg model of the bootstrap
c  current (subroutine NCLBOOT). The metrics consist of a number of geometric
c  and flux-surface-averaged quantities.
c  Reference: Houlberg, W.A. et al (1995). Nuclear Fusion (to be published)
c
c  All output variables are passed through the argument list back to FLUXAV
c  (the module that called this routine), which stores the metrics in an INCLUDE
c  file named "metrics.i" for use in NCLBOOT
c
c  The metrics are calculated here on the psi grid and then later interpolated
c  on to the rho grid in module RHOSET as with other flux-surface averaged
c  quantities calculated in FLUXAV.
c
c  This file was written using GNU EMACS, an excellent and FREE editor.
c  Support free software!
c
c     Created:   MAY-19-1995   Daniel Finkenthal
c     Modified:  JUN-20-1995   Daniel Finkenthal
c
c ----------------------------------------------------------------------
c
c --- Input through argument list:
c
c  xmagax            Radial location of magnetic axis
c
c  ymagax            Z location of magnetic axis
c
c  mp                Number of points on contour
c
c  xp0(1..mp)        locations of points around contour
c
c  yp0(1..mp)        locations of points around contour
c
c  bp0(1..mp)        poloidal B field at contour locations
c
c  fpsi              poloidal flux function from EQDSK file
c
c  ffprime           another  flux function from EQDSK file
c
c  qpsi              safety factor defined on psival grid
c
c  x(1..nx)          mesh locations of contour-spline
c
c  y(1..ny)          mesh locations of contour-spline
c
c  cspln(nx,ny)      values of contour-spline
c
c --- Input through argument INCLUDE files:
c
c   mhdpar:
c
c   storage:
c
c --- Output through argument list:
c     The Houlberg Metrics!
c
c  xngrth            n.grad(theta) (flux function)-radians/m
c
c  xb2               <B**2>                       -T**2
c
c  xbm2              <B**-2>                      -1/T**2
c
c  xgrbm2            <|grad psi|**2/B**2>         -1/T**2
c
c  fm(j)             poloidal expansion of geometric factor, j=1,3
c
c --- There is no output to include files!
c
c ----------------------------------------------------------------------
c

      USE mhdpar          
      USE replace_imsl,                   ONLY : my_dbcevl1
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'mhdpar.i'
      include 'storage.i'
c
c  Local Variables:
c
      integer    mp, nconmax, nwork
      parameter (nconmax = kstore)
      parameter (nwork   = kstore)
c
      integer    nx,ny
      real*8     x(*),y(*),cspln(2,nx,*),pds(6)
      real*8     xmagax, ymagax
      real*8     xp(nconmax),yp(nconmax),arclen(nconmax)
      real*8     xp0(nconmax),yp0(nconmax)
      real*8     ctheta(nconmax)
      real*8     fpsi,ffprime,fprime
      real*8     bt(nwork),btot(nwork)
      real*8     bp(kstore), bpinv(kstore),bp0(kstore)
      real*8     nthinv(kstore)
      real*8     pi
      real*8     ndgradb(nwork),fm(6)
      real*8     xngrth,gamma
      real*8     xb2,avbtot2
      real*8     xbm2, avbinv2, xgrbm2
      real*8     rmaj,rmaj2
      real*8     dbrdr,dbzdr,dbtdr,dbdr
      real*8     dbrdz,dbzdz,dbtdz,dbdz
      real*8     eqn41,eqn42
      real*8     capbr,capbz,sum2,bpinteg,avbp
      real*8     ydif,mindif
      integer    inew, iold, istart
      integer    i, m, ierr
c
c diagnostic variables:
c
      character  bootname*10
      integer    idebug
      real*8     rdebug(10,kstore)
      save       jpoint
      data       jpoint /0/
c
c ----------------------------------------------------------------------
c Set some internal variables, used for debugging, etc
c ----------------------------------------------------------------------
c
      pi      = ACOS (-1.0)
      kboot   = 0
      idebug  = 0
      icircle = 0
      ncrt    = 6
      jpoint  = jpoint + 1
c
c ----------------------------------------------------------------------
c Calculate fprime:
c (ffprime and fpsi are known from subroutine FOFPSI in cray208.f)
c ----------------------------------------------------------------------
c
      fprime = ffprime / fpsi
c
c ----------------------------------------------------------------------
c Re-order contour arrays so that contour begins on the outer midplane.
c ----------------------------------------------------------------------
c
c  Find origin:
c
      istart = 0
      mindif = 1000000.0
      do i=1,mp
        ydif = ABS (ymagax-yp0(i))
        if (ydif .le. mindif .and. xp0(i).ge. xmagax) then
          istart = i
        end if
        mindif = MIN (ydif, mindif)
      end do
c
c  Reorder arrays to start at istart (origin):
c
      do inew=1,mp-1
        iold = inew + istart - 1
        if (iold .gt. mp)
     .  iold = iold - mp     + 1
        xp(inew) = xp0(iold)
        yp(inew) = yp0(iold)
        bp(inew) = bp0(iold)
      end do
c
      xp(mp) = xp(1)
      yp(mp) = yp(1)
      bp(mp) = bp(1)
c
c  Compute arclength around the contour (as done in FLUXAV)
c
      arclen(1) = 0.0
      do i=1,mp-1
        arclen(i+1) = arclen(i) + SQRT ((xp(i+1)-xp(i))**2
     .                                + (yp(i+1)-yp(i))**2)
      end do
c
c ----------------------------------------------------------------------
c  Begin loop around the contour
c ----------------------------------------------------------------------
c
      do i = 1, mp
        bpinv(i) = 1.0/ bp(i)
        bt(i)    = fpsi/xp(i)
c
c  Get psi and derivatives:
c
        call my_dbcevl1 (x,nx,y,ny,cspln,nx,xp(i),yp(i),pds,ierr,6)
c
c  Real psi is opposite sign to what is returned:
c
        pds(1)  = -pds(1)
        pds(2)  = -pds(2)
        pds(3)  = -pds(3)
        pds(4)  = -pds(4)
        pds(5)  = -pds(5)
        pds(6)  = -pds(6)
c
        rmaj    = xp(i)
        rmaj2   = xp(i)*xp(i)
        capbr   = -pds(3)/xp(i)
        capbz   =  pds(2)/xp(i)
        btot(i) = SQRT (capbr*capbr+capbz*capbz+bt(i)*bt(i))
c
        dbrdr   =  pds(3)/(xp(i)*xp(i))-pds(4)/xp(i)
        dbzdr   = -pds(2)/rmaj2 + pds(5)/rmaj
        dbtdr   = -fpsi/rmaj2 + fprime*pds(2)/xp(i)
c
        dbdr    = (capbr*dbrdr+capbz*dbzdr+bt(i)*dbtdr)/btot(i)
c
        dbrdz   = -pds(6)/rmaj
        dbzdz   = pds(4)/rmaj
        dbtdz   = fprime*pds(3)/rmaj
c
        dbdz    = (capbr*dbrdz+capbz*dbzdz+bt(i)*dbtdz)/btot(i)
c
        ndgradb(i) = -(capbr*dbdr+capbz*dbdz)/btot(i)
c
c  Make use of the identity Lp/2pi*B dot grad(theta) = bp
c  Note: turns out that direction of bp does not matter (bp same as -bp)
c
        nthinv(i) = ABS (btot(i)) / (bp(i))
c
      end do
c
c ----------------------------------------------------------------------
c  END LOOP AROUND CONTOUR
c ----------------------------------------------------------------------
c
c ----------------------------------------------------------------------
c  Calculate CapTheta array and n dot grad(CapTheta)=gamma
c ----------------------------------------------------------------------
c
      call integ (arclen, nthinv, mp, sum2)
      gamma = 2.0 * pi / sum2
c
      ctheta(1) = 0.0
      do 210 i=2,mp
      call integ (arclen, nthinv, i, sum2)
 210  ctheta(i) = sum2 * gamma
c
c ----------------------------------------------------------------------
c  Calculate flux surface averages
c ----------------------------------------------------------------------
c
c  Must first get bpinteg:
c
      call integ (arclen,bpinv,mp,bpinteg)
c
c  Special quantity <n dot grad(small theta)>:
c
      do 212 i=1,mp
 212  ydum(i) = bpinv(i)*btot(i)/ nthinv(i)
      call integ (arclen,ydum,mp,sum2)
      avbp = 2.0 * pi * sum2 / arclen(mp) / bpinteg
c
      do i=1,mp
        rdebug(1,i) = nthinv(i)
        rdebug(2,i) = ydum(i)
      end do
c
c  Calculate <B**2>
c
      do 150 i=1,mp
 150  ydum(i) = bpinv(i)*btot(i)*btot(i)
      call integ (arclen,ydum,mp,sum2)
      avbtot2 = sum2 / bpinteg
c
c  Calculate <1/B**2>
c
      do 170 i=1,mp
 170  ydum(i) = bpinv(i)/(btot(i)*btot(i))
      call integ (arclen,ydum,mp,sum2)
      avbinv2 = sum2/bpinteg
c
c  Calculate <(grad psi/B)**2>:
c
      do 256 i=1,mp
 256  ydum(i) = xp(i)*xp(i)*bp(i) / ( btot(i)*btot(i) )
      call integ(arclen,ydum,mp,sum2)
      xgrbm2  = sum2 / bpinteg                       ! <(grad psi/B)**2>
c
      do i=1,mp
        rdebug(5,i) = ydum(i)
      end do
c
c  Rename the variables to be passed out:
c
      xngrth = gamma
      xb2    = avbtot2
      xbm2   = avbinv2
c
c ----------------------------------------------------------------------
c  Calculate the poloidal expansion factors Fm(1...3)
c  using equation 35 with eqns 41 and 42:
c ----------------------------------------------------------------------
c
      do m=1,3
        xm = FLOAT (m)
c
c  Calculate Eqn 41:
c
       do 220 i=1,mp
 220    ydum(i) = bpinv(i) * SIN (xm*ctheta(i)) * ndgradb(i)
        call integ (arclen,ydum,mp,sum2)
        eqn41=sum2/bpinteg
c
      do i=1,mp
        rdebug(3,i) = SIN (xm*ctheta(i))*ndgradb(i)
      end do
c
c  Calculate Eqn 42:
c
      do 230 i=1,mp
 230  ydum(i) = ydum(i) * gamma * ABS (btot(i))
      call integ (arclen, ydum, mp, sum2)
      eqn42 = sum2 / bpinteg
c
      do i=1,mp
        rdebug(4,i) = rdebug(3,i)*gamma * ABS (btot(i))
      end do
c
c  Put it all together to get Fm's:
c
        fm(m) = 2.0 / (avbtot2*avbp) * eqn41 * eqn42
c
      end do
c
c ----------------------------------------------------------------------
c This stuff is used ONLY for a check: It is valid only for circular,
c low-Beta case ONLY. (icircle = 1)
c ----------------------------------------------------------------------
c
      if (icircle .gt. 0) then
        if (jpoint .eq. 1)
     .  write (ncrt, '(a)') ' getting circular metrics'
        write (ncrt, '( )')
        write (ncrt, '(a6, 10(2x, e12.3))') ' Hlbg', pds(1), fm(1),
     .                fm(2), fm(3), fm(1)+fm(2)+fm(3), xngrth, xqs, xeps
c
c *** kboot-internal option for metrics.
c *** kboot=1 for Hirshman, Sigmar model: fm(k)=L_ck^*,
c ***                                     xngrth=<(n.gr(B))^2>/<B^2>.
c ***      =2 for Kessel model:           fm(1)=Rq/eps**3/2.
c ***                                     fm(2)=eps**3/2.
c ***                                     xngrth not used.
c ***      =else for Shaing, et al model: fm(k)=F_k,
c ***                                     xngrth=<n.gr(theta)>
c
        kboot = 0
        if (kboot .eq. 1) then          ! Hirshman
          xngrth = (xeps/(xqs*r0))**2/2.0
          fm(1)  =  xqs*r0
          fm(2)  =  0.0
          fm(3)  =  0.0
        else if (kboot .eq. 2) then     ! Kessel
          xngrth = 0.0
          fm(1)  = xqs * r0 / xeps**1.5
          fm(2)  =            xeps**1.5
        else                            ! Shaing formulation
c
c         xngrth-local n.grad(theta) (flux function)-radians/m.
c
          xngrth = 1.0 / (xqs * r0)
c
c         fm(j)-local poloidal expansion of geometric factor.
c
          call sscal (3, zero, fm, 1)
          if (xeps .gt. 0.0) then
            eps2 = xeps**2
            b    = SQRT (1.0 - eps2)
            a    = (1.0 - b) / xeps
            c    = (b**3.0) * (xqs*r0)**2
            do m=1,3
              xm    = FLOAT (m)
              fm(m) = xm * a**(2.0*xm) * (1.0+xm*b) / c
            end do
          end if
        end if
c
        write (ncrt, '(a6, 2x, i12, 7(2x,e12.3))') 'Shaing', j,
     .                 fm(1), fm(2), fm(3), fm(1)+fm(2)+fm(3), xngrth
c
      end if    ! CIRCLE
c
c ----------------------------------------------------------------------
c This stuff is ONLY for diagnostic purposes (idebug = 1)
c ----------------------------------------------------------------------
c
      if (idebug .eq. 1) then
c
c Open the "metrics" file if this is the first call
c
        if (jpoint .eq. 1) then
          call getioun(n71,71)
          open   (unit = n71, file = 'metrics', status = 'UNKNOWN')
          write  (n71, 1050)
 1050     format (4x,'Psi',6x,'xfm(1)',6x,'xfm(2)',5x,'xfm(3)',7x,
     .           'sum_xfm',7x,'Check',7x,'Diff(%)',7x,'xngrth',
     .            7x,'xgrbm2')
        end if
c
c       Check Fm(1..3)
c
        do i=1,mp
          ydum(i) = ndgradb(i)*ndgradb(i)*bpinv(i)
        end do
        call integ (arclen, ydum, mp, sum2)
        checkfm = sum2 / bpinteg / xb2
        sumfm   = fm(1) + fm(2) + fm(3)
c
        write  (n71, 1056) pds(1), fm(1), fm(2), fm(3), sumfm, checkfm,
     .                   (sumfm-checkfm)/checkfm*100.0, xngrth, xgrbm2
 1056   format (1x, 1p12e12.3)
c
        do i=1,mp
          rdebug(6,i) = ndgradb(i)
        end do
c
        write  (bootname, '(a7, i3.3)')  'bootout', jpoint
        call getioun(n72,72)
        open   (unit = n72, file = bootname , status = 'UNKNOWN')
        write  (n72, 1058)
 1058   format (1x,'Arclen',6x,'CapTheta',4x,'Bpinv',7x,'nthinv',
     .          6x,'B/nthin/Bp',2x,'Eqn_41',6x,'Eqn_42',6x,'grpsi/B',
     .          4x,'ndgradB')
        do i=1,mp
          write (n72, 1056)  arclen(i), ctheta(i), bpinv(i),
     .                     (rdebug(k,i),k=1,6)
        end do
        close (unit = n72)
      end if     ! idebug
      return
c
      end

      subroutine nclmn
c
c
c***********************************************************************
c *** NCLMN calculates the test particle (M) and field particle (N)    *
c ***  matrix elements of the collision operator using the Laguerre    *
c ***  polynomials of order 3/2 as basis functions for each isotopic   *
c ***  species combination.                                            *
c *** References:                                                      *
c *** Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079 (HS81).             *
c *** Hirshman, Phys Fluids 20 (1977) 589.                             *
c *** Houlberg, Shaing, Hirshman, to be published.                     *
c *** Last Revision: 5/96 W.A.Houlberg, ORNL.                          *
c *** Comments:                                                        *
c *** The indices on the matrices are one greater than the notation in *
c ***  the review article so as to avoid 0 as an index.                *
c***********************************************************************
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'param.i'
      include 'houlberg.i'
c
c *** Loop over isotope a.
c
      do ima=1,mi
c
c ***   Loop over isotope b.
c
        do imb=1,mi
c
c ***     Ratio of masses.
c
          xmab=amuai(ima)/amuai(imb)
c
c ***     Ratio of temperatures.
c
          xtab=tai(ima)/tai(imb)
c
c ***     Ratio of thermal velocities, vtb/vta.
c
          xab = SQRT (xmab/xtab)
c
c ***     Elements of M.
c
          xab2=xab**2
          yab32=(1.0+xab2) * SQRT (1.0+xab2)
          yab52=(1.0+xab2) * yab32
          yab72=(1.0+xab2) * yab52
          yab92=(1.0+xab2) * yab72
c
c ***     Eqn 4.11 for M00 (HS81).
c
          capm(1,1,ima,imb)=-(1.0+xmab)/yab32
c
c ***     Eqn 4.12 for M01 (HS81).
c
          capm(1,2,ima,imb)=3.0/2.0*(1.0+xmab)/yab52
c
c ***     Eqn 4.15 for M02 (HS81).
c
          capm(1,3,ima,imb)=-15.0/8.0*(1.0+xmab)/yab72
c
c ***     Eqn 4.8 for M10 (HS81).
c
          capm(2,1,ima,imb)=capm(1,2,ima,imb)
c
c ***     Eqn 4.13 for M11 (HS81).
c
          capm(2,2,ima,imb)=-(13.0/4.0+xab2*(4.0+xab2*15.0/2.0))/yab52
c
c ***     Eqn 4.16 for M12 (HS81).
c
          capm(2,3,ima,imb)=(69.0/16.0+xab2*(6.0+xab2*63.0/4.0))/yab72
c
c ***     Eqn 4.8 for M20 (HS81).
c
          capm(3,1,ima,imb)=capm(1,3,ima,imb)
c
c ***     Eqn 4.8 for M21 (HS81).
c
          capm(3,2,ima,imb)=capm(2,3,ima,imb)
c
c ***     Eqn 5.21 for M22 (HS81).
c
          capm(3,3,ima,imb)=-(433.0/64.0+xab2*(17.0+xab2*(459.0/8.0
     .                      +xab2*(28.0+xab2*175.0/8.0))))/yab92
c
c ***     Elements of N.
c ***     Momentum conservation, Eqn 4.11 for N00 (HS81).
c
          capn(1,1,ima,imb)=-capm(1,1,ima,imb)
c
c ***     Eqn 4.9 and 4.12 for N01 (HS81).
c
          capn(1,2,ima,imb)=-xab2*capm(1,2,ima,imb)
c
c ***     Eqn 4.15 for N02 (HS81) - corrected rhs by Ta/Tb.
c
          capn(1,3,ima,imb)=-xab2**2*capm(1,3,ima,imb)
c
c ***     Momentum conservation, Eqn 4.12 for N10 (HS81).
c
          capn(2,1,ima,imb)=-capm(2,1,ima,imb)
c
c ***     Eqn 4.14 for N11 (HS81).
c
          capn(2,2,ima,imb)=(27.0/4.0)*xtab*xab2/yab52
c
c ***     Eqn 4.17 for N12 (HS81).
c
          capn(2,3,ima,imb)=-225.0/16.0*xtab*xab2**2/yab72
c
c ***     Momentum conservation for N20 (HS81).
c
          capn(3,1,ima,imb)=-capm(3,1,ima,imb)
c
c ***     Eqn 4.9 and 4.17 for N21 (HS81).
c
          capn(3,2,ima,imb)=-225.0/16.0*xab2**2/yab72
c
c ***     Eqn 5.22 for N22 (HS81).
c
          capn(3,3,ima,imb)=2625.0/64.0*xtab*xab2**2/yab92
        end do
      end do
      return
c
      end

      subroutine nclmu
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** NCLMU calculates the neoclassical viscosities by integrating the *
c ***  velcoity-dependent banana and Pfirsch-Schluter contributions    *
c ***  over velocity.                                                  *
c *** References:                                                      *
c *** Shaing, Yokoyama, Wakatani, Hsu, Phys Plasmas 3 (1996) 965.      *
c *** Houlberg, Shaing, Hirshman, to be published.                     *
c *** Last Revision: 5/96 W.A.Houlberg, ORNL.                          *
c***********************************************************************
c
      parameter      (mpnts=13)
c
c      include 'param.i'
      include 'houlberg.i'
c
      save           c3,              x,               w
      save           init
      dimension      x(mpnts),        w(mpnts,5)
      data bmax/     3.2/
      data init/     0/
c
      iflag=0
      if (init.eq.0) then
c
c ***   Initialize integration points and weights for integration.
c
        h=bmax/(mpnts-1)
        x(1)=0.0
        w(1,1)=0.0
        w(1,2)=0.0
        w(1,3)=0.0
        w(1,4)=0.0
        w(1,5)=0.0
        do i=2,mpnts
          x(i)=h*(i-1)
          x2=x(i)*x(i)
          expmx2=exp(-x2)
          x4=x2*x2
          w(i,1)=x4*expmx2
          x6=x4*x2
          w(i,2)=x6*expmx2
          x8=x4*x4
          w(i,3)=x8*expmx2
          x10=x4*x6
          w(i,4)=x10*expmx2
          x12=x6*x6
          w(i,5)=x12*expmx2
        end do
        c3 = 8.0 / 3.0 / SQRT (pi) * h
        init=1
      end if
c
c *** Initialize viscosities and constants.
c
      mimz9 = 9 * mxmi * mz
      call sscal (mimz9, zero, xmu, 1)
c
c *** Loop over grid (first node has null value).
c
      do m=2,mpnts
          if (m.eq.mpnts) then
c
c ***     Use half weight for end point.
c
          ewt=0.5
          else
c
c ***     Use full weight.
c
          ewt=1.0
        end if
        xx=x(m)
c
c ***   Get velocity-dependent k values.
c

        call nclk(xx)

c
c ***   Loop over isotopes.
c
        do ima=1,mi
c
c ***     Loop over charge states.
c
          do ia=1,mz
            if (denai(ima,ia).gt.dencut) then
              c4=c3*denai(ima,ia)*amuai(ima)*promas
              if (kboot.gt.0.and.kboot.le.4) then
c
c ***           H&S, Kessel, banana only formulations.
c
                xk=xkban(ima,ia)
              else
c
c ***           Shaing, et al., formulation.
c
                xk=xkban(ima,ia)/sqz(ima,ia)*xkps(ima,ia)
     .             /(xkban(ima,ia)/sqz(ima,ia)+xkps(ima,ia))
              end if
              xmu(1,1,ima,ia)=xmu(1,1,ima,ia)+ewt*c4*xk*w(m,1)
              xmu(1,2,ima,ia)=xmu(1,2,ima,ia)+ewt*c4*xk
     .                        *(w(m,2)-5.0/2.0*w(m,1))
              xmu(2,2,ima,ia)=xmu(2,2,ima,ia)+ewt*c4*xk
     .                        *(w(m,3)-5.0*w(m,2)+(25.0/4.0)*w(m,1))
             xmu(1,3,ima,ia)=xmu(1,3,ima,ia)+ewt*c4*xk
     .                        *((1.0/2.0)*w(m,3)-(7.0/2.0)*w(m,2)
     .                        +(35.0/8.0)*w(m,1))
             xmu(2,3,ima,ia)=xmu(2,3,ima,ia)+ewt*c4*xk
     .                        *((1.0/2.0)*w(m,4)-(19.0/4.0)*w(m,3)
     .                        +(105.0/8.0)*w(m,2)-(175.0/16.0)*w(m,1))
             xmu(3,3,ima,ia)=xmu(3,3,ima,ia)+ewt*c4*xk
     .                        *((1.0/4.0)*w(m,5)-(7.0/2.0)*w(m,4)
     .                        +(133.0/8.0)*w(m,3)-(245.0/8.0)*w(m,2)
     .                        +(1225.0/64.0)*w(m,1))
            end if
          end do
        end do
      end do

c
c *** Fill out matrix.
c
        do        l=1,2
          do      k=l+1,3
            do  ima=1,mi
              do ia=1,mz
                xmu(k,l,ima,ia)=xmu(l,k,ima,ia)
              end do
            end do
          end do
        end do
      return
c
      end

      subroutine nclnu (xx)
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** NCLNU calculates the velocity dependent nu_D, nu_T, and          *
c *** nu_T*I_Rm for a test particle in a multiple species plasma.      *
c *** References:                                                      *
c *** Hirshman, Sigmar, Phys Fluids 19 (1976) 1532.                    *
c *** Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079.                    *
c *** Shaing, Yokoyama, Wakatani, Hsu, Phys Plasmas 3 (1996) 965.      *
c *** Houlberg, Shaing, Hirshman, to be published.                     *
c *** Last Revision: 5/96 W.A.Houlberg, ORNL.                          *
c *** Input:                                                           *
c *** xx-velocity normalized to thermal velocity-dimensionless.        *
c***********************************************************************
c
c      include 'param.i'
      include 'houlberg.i'
c
c *** Zero out collisionalities.
c
      mimz = mxmi * mz
      call sscal(mimz,zero,xnud,1)
      call sscal(mimz,zero,xnut,1)
      mimzmp = mimz*3
      call sscal(mimzmp,zero,xnuti,1)
c
c *** Calculate nu_D and nu_T.
c

      do ima=1,mi
        do imb=1,mi
          xab = vtai(imb) / vtai(ima)
          xb = xx / xab
          call uchand (xb, g, phi)
          c1 = (3.0 * SQRT (pi) /4.0) * ((phi-3.0*g)/xx**3
     .       +  4.0 * (tai(ima)/tai(imb) + 1.0/xab**2)*g/xx)
          c2 = (3.0 * SQRT (pi) / 4.0) * (phi-g)/xx**3
          do ia=1,mz
            if (denai(ima,ia).gt.dencut) then
              do jb=1,mz
                if (denai(imb,jb).gt.dencut) then
                  xnud(ima,ia)=xnud(ima,ia)+c2/tauab(ima,imb,ia,jb)
                  xnut(ima,ia)=xnut(ima,ia)+c1/tauab(ima,imb,ia,jb)
                end if
              end do
            end if
          end do
        end do
      end do


      if (kboot.lt.1.or.kboot.gt.4) then
c
c ***   Shaing, et al., formulation.
c ***   Calculate nu_T*I_m.
c

        do ima=1,mi
          do mp=1,3
            wm=xx*vtai(ima)*mp*xngrth
            do ia=1,mz
              if (denai(ima,ia).gt.dencut) then
                xnw=xnut(ima,ia)/wm
                xnw2=xnw**2
                if (xnw.gt.3.0) then
                  xnuti(mp,ima,ia)=0.4
                else
                  xtani=xnw*atan(1.0/xnw)
                  xnuti(mp,ima,ia)=0.5*xtani+xnw2*(3.0*(xtani-0.5)
     .                                      +xnw2*4.5*(xtani-1.0))
                end if
              end if
            end do
          end do
        end do
      end if


      return
c
      end

      subroutine nclout (nunncl, time, t, k, r, nj, jprt, jhirsh,
     .                   curden, curohm, curboot, curbi,
     .                   curbe, curbet, currf, curdri,
     .                   totcur, totohm, totboot, totbi, totbe,
     .                   totbet, totrf, totdri)
c
c
c ----------------------------------------------------------------------
c     subroutine NCLOUT creats a special output file named "nclout"
c     which contains the current profiles
c ----------------------------------------------------------------------
c
      USE mhdpar   
      implicit  integer (i-n), real*8 (a-h, o-z)
c      include 'mhdpar.i'
      include 'storage.i'
c
      real*8  r(kstore), curden(kstore), curohm(kstore),
     .        curboot(kstore), curbi(kstore), curbe(kstore),
     .        curbet (kstore), currf(kstore), curdri(kstore)
      real*8  totcur, totohm, totboot, totbi, totbe,
     .        totbet, totrf, totdri
      logical opened
      save    opened
      data    opened /.false./
c
      if (.not. opened) then
        call getioun(nunncl,nunncl)
        open (unit = nunncl, file = 'nclout', status = 'UNKNOWN')
        opened = .true.
      end if
c
      call header (nunncl, time, t)
      write (nunncl, 1052) jhirsh
      write (nunncl, 1053)
c
      do 640 j=1,nj
        j1prt = (j-1) / jprt * jprt
        if (j1prt .ne. j-1 .and. j .ne. nj)  go to 640
        write (nunncl, 1056) j, r(j),
     .                     curden(j), curohm(j), curboot(j), curbi (j),
     .                     curbe (j), curbet(j), currf  (j), curdri(j)
  640 continue
c
      write (nunncl, 1054) totcur, totohm, totboot, totbi, totbe,
     .                     totbet, totrf, totdri
      write (nunncl, 1055)
c
 1052 format (5x, 'KBOOT = ', i2.2)
 1053 format (5x, 'current densities (A/cm**2) and total currents (a)'/
     .        5x, '(curbi, curbe,curbet, currf and curohm at time '   /
     .        5x, 'point n-1/2, curden at time point n)'              //
     .        4x, 'j',5x,'r(cm)',6x,'curden',6x,'curohm',5x,'curboot',
     .        7x, 'curbi',7x,'curbe',7x,'curbet',
     .        7x, 'currf',5x,'curdrive')
 1054 format (/ 10x, 'total', 1p8e12.3)
 1055 format ('----------------------------------------------------')
 1056 format (1x, i4, f10.2, 1p8e12.3)
c
      return
c
      end

      subroutine ncltau
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** NCLTAU calculates the collision time for 90 degree scattering of *
c ***  each charge state of isotope a on each charge state of isotope b*
c ***  in a multiple species plasma assuming a common value for the    *
c ***  Coulomb logarithm necessary for the reduced charge state model. *
c *** References:                                                      *
c *** Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079.                    *
c *** Houlberg, Shaing, Hirshman, to be published.                     *
c *** Last Revision: 5/96 W.A.Houlberg, ORNL.                          *
c***********************************************************************
c
c      include 'param.i'
      include 'houlberg.i'
c
      dimension  xlnab(mxmi,mxmi), xz(mxmi), xz2n(mxmi)
c
c *** Initialize collision times.
c *** Calculate Coulomb logarithm.
c *** Get weight factors and set electron lambda.
c
      mi2    = mxmi**2
      call sscal (mi2   , zero, amntau, 1)
      mi2mz2 = mi2*mxmz*mz
      call sscal (mi2mz2, zero, tauab , 1)
      do ima=1,mi
        if (amuai(ima).lt.0.5) then
          xlnabe = 37.8 - LOG (SQRT (denai(ima,1)) / tai(ima))
        else
          sumn=0.0
          sumnz=0.0
          xz2n(ima)=0.0
          do ia=1,mz
            if (denai(ima,ia).gt.dencut) then
              sumn=sumn+denai(ima,ia)
              sumnz = sumnz + FLOAT (ia) * denai(ima,ia)
              xz2n(ima) = xz2n(ima) + FLOAT (ia)**2 * denai(ima,ia)
            end if
          end do
          xz(ima)=sumnz/sumn
        end if
      end do
      do ima=1,mi
        do imb=1,mi
          if (amuai(ima).lt.0.5.or.amuai(imb).lt.0.5) then
            xlnab(ima,imb)=xlnabe
          else
            t1 = xz(ima)*xz(imb)*(amuai(ima)+amuai(imb))
     .         / (amuai(ima)*tai(imb)+amuai(imb)*tai(ima))
            t2 = SQRT (xz2n(ima) / tai(ima) + xz2n(imb) / tai(imb))
            xlnab(ima,imb) = 40.3 - LOG (t1*t2)
          end if
        end do
      end do
      c1 = 4.0 / 3.0 / SQRT (pi) * 4.0 * pi
     .     * (coulom / (4.0 * pi * epsilo))**2 * (coulom/promas)**2
c
c *** Calculate collision times.
c
      do ima=1,mi
        c2=(vtai(ima)**3)*amuai(ima)**2/c1
        do imb=1,mi
          c3=c2/xlnab(ima,imb)
          do ia=1,mz
            do jb=1,mz
              if ((denai(ima,ia).gt.dencut).and.
     .            (denai(imb,jb).gt.dencut)) then
                tauab(ima,imb,ia,jb) = c3 / (FLOAT(ia)**2)
     .                               /(denai(imb,jb) * FLOAT (jb)**2)
                amntau(ima,imb)=amntau(ima,imb)
     .                          +amuai(ima)*promas*denai(ima,ia)
     .                          /tauab(ima,imb,ia,jb)
              end if
            end do
          end do
        end do
      end do
      return
c
      end

      real*8 function ssum (n, sx, incx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** SSUM is the sum of the elements of the array sx.                 *
c *** Last Revision: 1/95 W.A.Houlberg and S.E.Attenberger, ORNL.      *
c *** Input:                                                           *
c *** n-number of elements to be summed.                               *
c *** sx-vector to be summed.                                          *
c *** incx-increment in sx index.                                      *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
      dimension sx(*)
c
      ssum = 0.0
      ix   = 1
      do i=1,n
        ssum = ssum + sx(ix)
        ix   = ix   + incx
      end do
      return
c
      end

      subroutine uchand (x, g, phi)
c
      USE param
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** UCHAND evaluates the Chandresekar and error functions.           *
c *** Last Revision: 1/95 W.A.Houlberg and S.E.Attenberger, ORNL.      *
c *** Input:                                                           *
c *** x-argument of Chandresekar and error functions.                  *
c *** Output:                                                          *
c *** g-value of Chandresekar function.                                *
c *** phi-value of error function.                                     *
c *** Other Comments:                                                  *
c *** G(x)=[phi(x)-x*phi'(x)]/2x**2.                                   *
c *** phi'(x)=(2/pi**1/2)exp(-x**2).                                   *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
      real*8 g, uerf, x
c
c      include 'param.i'
      include 'houlberg.i'
c

      phi  = uerf(x)
      xsq = x*x
      if(xsq .gt. 400.)then  !asymtotic matching not done HSJ
         phip =0.0
      else
         phip = (2.0 / SQRT (pi)) * EXP (-xsq)
      endif
      g    = (phi-x*phip)/(2.0*x**2)
      return
c
      end

      real*8 function uerf (x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** UERF evaluates the error function erf(x).                        *
c *** Last Revision: 1/95 W.A.Houlberg and S.E.Attenberger, ORNL.      *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
      real*8  ugammp, x
c
c      if (x .lt. 0.0) then
c        uerf = -ugammp (0.5, x**2)
c       else
c        uerf =  ugammp (0.5, x**2)
c      endif

      xsq = x*x
      if(xsq  .lt. 600.)then             !underflow trap on some machines HSJ
         uerf =SIGN(ugammp(0.5,xsq),x)
      else
         uerf =SIGN(1.D0,x) !jmp.ibm
      endif
      return
      

c
      end

      real*8 function ugamln (xx)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** UGAMLN evaluates the logarithm of the gamma function.            *
c *** Last Revision: 1/95 W.A.Houlberg and S.E.Attenberger, ORNL.      *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
      integer        j
      real*8         xx
      real*8         cof(6),          ser,             stp,
     .               tmp,             x,               y
      save           cof,             stp
      data cof/       76.18009172947146,
     .               -86.50532032941677,
     .                24.01409824083091,
     .                -1.231739572450155,
     .                 1.208650973866179e-3,
     .                -5.395239384953e-6/
      data stp/        2.5066282746310005/
c
      x   =  xx
      y   =  x
      tmp =  x + 5.5
      tmp = (x + 0.5) * LOG (tmp) - tmp
      ser = 1.000000000190015
      do j=1,6
        y   = y   + 1.0
        ser = ser + cof(j) / y
      end do
      ugamln = tmp + LOG (stp * ser / x)
      return
c
      end

      real*8 function ugammp (a, x)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** UGAMMP evaluates the gamma function.                             *
c *** Last Revision: 1/95 W.A.Houlberg and S.E.Attenberger, ORNL.      *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
      real*8  a, gammcf, gamser, gln, x
c
      if (x .lt. 0.0 .or. a .le. 0.0)
     .  call STOP ('function UGAMMP: bad arguments to function', 177)
      if (x .lt. a+1.0) then
        call ugser (gamser, a, x, gln)
        ugammp = gamser
      else
        call ugcf (gammcf, a, x, gln)
        ugammp = 1.0 - gammcf
      end if
      return
c
      end

      subroutine ugcf (gammcf, a, x, gln)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** UGCF evaluates a continued fraction used by UGAMMP.              *
c *** Last Revision: 1/95 W.A.Houlberg and S.E.Attenberger, ORNL.      *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
      integer    i, itmax
      real*8     a, an, b, c, d, del, eps, fpmin, gammcf, gln, ugamln, x
c
      parameter (itmax = 100, eps = 3.0e-7, fpmin = 1.0e-30)
c
      gln = ugamln(a)
      b   = x+1.0-a
      c   = 1.0 / fpmin
      d   = 1.0 / b
      h   = d
      do i=1,itmax
        an  = -i*(i-a)
        b   = b+2.0
        d   = an*d+b
        if (ABS (d) .lt. fpmin)  d = fpmin
        c   = b+an/c
        if (ABS (c) .lt. fpmin)  c = fpmin
        d   = 1.0 / d
        del = d * c
        h   = h * del
        if (ABS (del-1.0) .lt. eps)  go to 10
      end do
      call STOP ('subroutine UGCF: A too large, ITMAX too small', 178)
   10 gammcf = EXP (-x + a * LOG (x) - gln) * h
      return
c
      end

      subroutine ugser (gamser, a, x, gln)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** UGSER evaluates a series used by UGAMMP.                         *
c *** Last Revision: 1/95 W.A.Houlberg and S.E.Attenberger, ORNL.      *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
      integer        itmax,           n
      real*8         a,               ap,              del,
     .               eps,             gamser,          gln,
     .               sum,             ugamln,          x
c
      parameter     (itmax = 100,     eps = 3.0e-7)
c
      gln = ugamln(a)
      if   (x .le. 0.0) then
        if (x .lt. 0.0)
     .    call STOP ('subroutine UGSER: x is less than zero', 179)
        gamser = 0.0
        return
      end if
      ap  = a
      sum = 1.0 / a
      del = sum
      do n=1,itmax
        ap  = ap+1.0
        del = del*x/ap
        sum = sum+del
        if (ABS (del) .lt. ABS (sum)*eps)  go to 10
      end do
      call STOP ('subroutine UGSER: A too large, ITMAX too small', 180)
   10 gamser = sum * EXP (-x + a * LOG (x) - gln)
      return
c
      end

      subroutine ulbksb (a, n, np, indx, b)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** ULBKSB solves the matrix equation a x = b, where a is in lu form *
c ***  as output by subroutine ULDCMP.                                 *
c *** References:                                                      *
c *** Flannery, Teukolsky, Vetterling, Numerical Recipes.              *
c *** Last Revision: 1/95 W.A. Houlberg and S.E.Attenberger, ORNL.     *
c *** Input:                                                           *
c *** a-matrix in lu decomposed form.                                  *
c *** n-dimension of matrix.                                           *
c *** np-first dimension of a array.                                   *
c *** indx-vector showing row permutations due to partial pivoting.    *
c *** b-right hand side of equation (overwritten on return).           *
c *** Output:                                                          *
c *** b-solution vector.                                               *
c *** Comments:                                                        *
c *** The complete procedure is: given a x = b                         *
c ***                            call uldcmp(a,n,np,indx,d,ierr)       *
c ***                            if (ierr.eq.0)                        *
c ***                              call ulbksb(a,n,np,indx,b)          *
c ***  and b is now the solution vector.                               *
c *** To solve with a different right hand side, just reload b as      *
c ***  desired and use the same lu decomp, call ulbksb(a,n,np,indx,b). *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
      dimension      a(np,np),        indx(n),         b(n)
c
c --- Find the index of the first nonzero element of b.
c --- Unscramble the permutation of rows as we go.
c
      ii = 0
      do i=1,n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if (ii .ne. 0) then
          do j=ii,i-1
            sum = sum-a(i,j)*b(j)
          end do
        else if (sum.ne.0.0) then
          ii = i
        end if
        b(i) = sum
      end do
c
c --- Begin back substitution.
c
      do i=n,1,-1
        sum = b(i)
        if (i .lt. n) then
          do j=i+1,n
            sum = sum - a(i,j) * b(j)
          end do
        end if
        b(i) = sum / a(i,i)
      end do
      return
c
      end

      subroutine uldcmp (a, n, np, indx, d, ierr)
c
      implicit  integer (i-n), real*8 (a-h, o-z)
c
c***********************************************************************
c *** ULDCMP does an lu decomposition of the matrix a, is used in      *
c ***  conjuction with ULBKSB to solve linear equations or to invert   *
c ***  a matrix.                                                       *
c *** References:                                                      *
c *** Flannery, Teukolsky, Vetterling, Numerical Recipes.              *
c *** Last Revision: 1/95 W.A. Houlberg and S.E.Attenberger, ORNL.     *
c *** Input:                                                           *
c *** a-square matrix, overwritten on return.                          *
c *** n-dimension of matrix.                                           *
c *** np-first dimension of a array.                                   *
c *** Output:                                                          *
c *** indx-vector showing row permutations due to partial pivoting.    *
c *** d-flag for number of row exchanges.                              *
c ***  = 1 even number.                                                *
c ***  =-1 odd  number.                                                *
c *** ierr-error flag.                                                 *
c ***     =0 normal exit.                                              *
c ***     =1 singular matrix.                                          *
c ***                                                                  *
c *** Adapted for ONETWO by D. Finkenthal   6/95                       *
c***********************************************************************
c
      parameter (nmax = 100)
      dimension  a(np,np), indx(n), vv(nmax)
c
      ierr = 0
      d    = 1.0
c
c --- Loop over rows to get the implicit scaling information
c
      do i=1,n
        aamax = 0.0
        do j=1,n
          if (ABS (a(i,j)) .gt. aamax)  aamax = ABS (a(i,j))
        end do
        if (aamax .eq. 0.0) then
          ierr = 1
          return
        end if
        vv(i) = 1.0 / aamax
      end do
c
c --- Loop over columns of Crout's method
c
      do   j=1,n
        do i=1,j-1
          sum = a(i,j)
          do k=1,i-1
            sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
        end do
c
c ----- Search for largest pivot element (dum is a figure of merit)
c
        aamax = 0.0
        do i=j,n
          sum = a(i,j)
          do k=1,j-1
            sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
          dum    = vv(i) * ABS (sum)
          if (dum .ge. aamax) then
            imax  = i
            aamax = dum
          end if
        end do
        if (j .ne. imax) then
c
c ------- Interchange rows
c
          do k=1,n
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
          end do
          d = -d
          vv(imax) = vv(j)
        end if
        indx(j) = imax
        if (a(j,j).eq.0.0) then
          ierr = 1
          return
        end if
        if (j .ne. n) then
c
c ------- Divide by pivot element
c
          dum = 1.0 / a(j,j)
          do i=j+1,n
            a(i,j) = a(i,j)*dum
          end do
        end if
      end do
      return
c
      end
