
MODULE colrate
      IMPLICIT NONE
      SAVE
! --- used to pass information to various subroutines
! --- for determination of collision rate integrals
!
      INTEGER, PARAMETER :: n1v       = 100  ! thermal
      INTEGER, PARAMETER :: n2v       = 100  ! fast, used for beam
      INTEGER, PARAMETER :: nzeta     =  21  ! max points in zeta quadrature
      INTEGER, PARAMETER :: nlegendre =  2  ! # Legendre polynomials = nlegendre+1
      INTEGER, PARAMETER :: n2vtble   = (n2v*(n2v+1))/2        ! symmetric storage
!
!     nzeta controls gaussian integration over zeta.
!     this calculation is done once, at beginning of run ONLY.
!
      INTEGER, PARAMETER ::  nvtherml = 50 ! max points in thermal speed quadrature
      INTEGER, PARAMETER ::  nvfast   = 50 ! max points in fast ion      quadrature
!
!     nvtherml and nvfast is number of points in gaussian quadrature of
!     thermal and fast ion distributions. has to be done each time
!     neutron rate is evaluated. hence make these as small as possible
!     to gain speed.
!
      REAL*8           mass_deut, mass_trit, mass_beam,mass_he3
!
      REAL*8          vth, vf, zetath, zetaf, phif, phith,            &
                      beta, wx, wy, wz, wflowv, vfsq, vthsq, betasq,  &
                      vb, vcrit, vc3a, rhe3max,rddnmax,               &
                      vc3i, vcrit3, z2, expd, expt, expdt,exphe3,     &
                      vb3, zetab, vc4, ecritlocal, den_neut,          &
                      tauslocal, vcutnterm, acore, taufctr, vbfrac,   &
                      sqpd, svthf, vthf, ztf, vfast_low_lim,          &
                      vfast_up_lim, sigvrddn(n2v,n1v), v1(n1v),       &
                      xzeta(nzeta), wzeta(nzeta), xvtherml(nvtherml), &
                      wxvtherml(nvtherml), xvfast(nvfast),            &
                      beam_beam_ddnhe3 (1:n2vtble,0:nlegendre),       &
                      beam_beam_ddpt   (1:n2vtble,0:nlegendre),       &
                      beam_beam_dtnhe4 (1:n2vtble,0:nlegendre),       &
                      beam_beam_tt2nhe4(1:n2vtble,0:nlegendre),       &
                      vbeam_beam(n2v), sigvrdtn(n2v,n1v),             &
                      sigvrddp (n2v,n1v), sigvrdhe3p(n2v,n1v),        &
                      sigvrtt2n(n2v,n1v), v2(n2v), wxvfast(nvfast)
      REAL*8 ::  twopi = 6.283185308 
      REAL*8 ::  umdd  = 5.2175e-16       ! (keV/(cm/sec)**2), for d-d
      INTEGER &
                      icxcalc, nvfact, ibeam_chg_exchg, nzeta_set,&
                      ib_d, ib_dt, ib_t, nterms
!
! --- vth is thermal speed, cm/sec
! --- vf  is fast ion speed
! --- vb  is speed of beam ions
! --- ecritlocal critical energy, keV, at which equal amounts of energy
! --- are transferred to electrons and ions
! --- mass_deut is mass of deuteron in grams
! --- mass_trit is mass of tritium  in grams
! --- wflowv magnitude of bulk motion speed,cm/sec
! --- zetath is pitch angle of thermal ion relative to B field
! --- zetaf  is pitch angle of fast    ion relative to B field
! --- expd, expt = m/(2kT) for m=d and m=t respectively
! --- phif  is azimuth of fast    ion speed
! --- phith is azimuth of thermal ion speed
! --- beta = m/(2*kt)
! --- wx,wy,wz bulk thermal ion speed
! --- vthf = 2*vth*vf
! --- beta2 = 2*beta*sqzth
! --- ts3tcx taus divided by 3*tcx
! --- vb3 = vb**3
! --- nterms number of terms used in sum of fast ion distribution
! --- zetab = (cos)pitch angle of beam
! --- den_neut  is neutral density, number/cm**3
! --- tauslocal is fast ion slowing-down time, sec
! --- vcutnterm is speed below which series summation of fast ion
! --- distribution is feasible.
! --- nvfact is multiplier on vcutnterm, don't change it
! --- acore is alpha defined in core's paper
! --- vbfrac is fraction of vbeam above which it is assumed that
! --- the fast ion pitch angle distribution is a delta FUNCTION at the
! --- initial beam pitch angle.
! --- umdd is 0.5*(mf*mth)/(mf+mth)*(conversion factor to keV)
! --- (i.e., um is defined so that um*vrel**2 is energy in center
! --- of mass system in keV, mf, mth mass of fast and thermal ions)
! --- sqpd  = SQRT ((1-zetaf**2)*(1-zetath**2))
! --- svthf = vth**2+vf**2
! --- ibeam_chg_exchg IF 1 INCLUDE charge exchange effects
! --- vcrit
! --- vcrit3
! --- z2
! --- betasq 2 * beta * SQRT (1-zetath*zetath)
! --- vfast_low_lim, vfast_up_lim, lower and upper limits of fast ion
! --- distribtution FUNCTION (accounts for transients due to beam turn
! --- on and off)
! --- vc3i = 1/(vf3+vcrit3)
! --- vthf = 2*vf*vth
! --- ztf  = zetath*zetaf
! --- i_rate  is reaction rate indicator
! ---            (see SUBROUTINE BEAM_THERMAL_RATE)
! --- ib_rate is reaction rate indicator for beam-beam reactions
! ---            (see SUBROUTINE BEAM_BEAM_RATE)
! --- nzeta_set =0 IF zeta (COS(polar angle) quadrature rule is not set
! ---           =1 IF quadrature is set up
! --- sigvrddn, sigvrdtn are tables used to look up the azimuthal
! - integrals for the d(d,n)he4 and d(t,n)he4 reactions.

END MODULE colrate
