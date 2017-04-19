!-------------------------------------------------------------------------
! ifs_pppl.f90
!
! PURPOSE:
!  Evaluate electron and ion heat flux, as return critical gradients,
!  using the (infamous) IFS-PPPL model.
!
! NOTES:
!  - Converted to Fortran 90 and simplified (Nov 2007).
!  - Original documentation also cleaned up (see below).
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
! ORIGINAL DOCUMENTATION:
!
! The formulas embodied in this subroutine are documented in the Physics
! of Plasmas article entitled Quantitative Predictions of Tokamak Energy 
! Confinement from First-Principles Simulations with Kinetic Effects, 
! by M. Kotschenreuther, W. Dorland, M.A. Beer, and G.W. Hammett, 
! Vol. 2, p. 2381, (1995).  Extensions to non-circular cross-sections are
! described below. 
!
! There is a significant typographical error in that paper.  R/Ln* is 
! defined to be max(6,R/Ln); it should be min(6,R/Ln) as defined in 
! this subroutine.
!
! Also, note that in deriving these formulas, we assumed that the 
! density gradient scale lengths for the different species were equal. 
! This is an approximation that needs to be relaxed.
!
! As emphasized in the paper, these formulas were derived numerically and 
! are therefore not trustworthy outside a particular region of parameter 
! space.  For example, we did not parameterize the heat flux in the weak 
! magnetic shear limit; thus, one should not use the model in this limit.
! I have attempted to reduce related strange numerical behaviors by 
! limiting some inputs to be roughly within their range of validity.  
!
! Stifness:
! ---------
! For many cases that we have simulated, the transport equations
! tend to be very stiff.  That is, the plasma temperature gradient scale 
! length tends to adjust itself to be close to the critical gradient scale
! length over some region of the plasma, because chi becomes very large
! very fast for shorter temperature gradient scale lengths.  Typically, 
! we have had to be very careful with the numerical algorithm used in the 
! transport equation solver with this experience in mind.  The details 
! of our implementation are available to anyone that is interested.
!
! Geometry:
! --------
! The nonlinear simulations that were done to obtain these formulas were 
! mostly done in a simplified geometry, using a shifted circle, low beta, 
! high aspect ratio expansion.  Some modifications due to more sophisticated 
! geometrical models have been calculated and have been included here, 
! but should be considered preliminary.  There are two important issues 
! that must be noted.  First, we derived our formulas using a different 
! radial coordinate.  Second, since we are actually calculating the 
! transport coefficients in general geometry, we require less assumptions 
! for the form of the transport equation to be solved.  
!
! Let me describe the <|grad rho|> issue first:
!
!     The database standard modeling assumptions that were agreed upon for
!     this exercise include the assumption that the anomalous fluxes for 
!     non-circular cross-sections are simply related to the anomalous 
!     fluxes for related circular cross-section plasmas.  That is, in order 
!     to get the factor of < |grad rho|**2 > that appears as a coefficient 
!     of chi in the energy transport equations, one assumes that 
!
!     chi_anom_general = chi_anom_circular * (grad rho).
!
!     One need not make this assumption; one can just calculate the quantity 
!     chi_anom_general directly.  One would then have a transport equation
!     of the form
!
!     (3/2) d(n T)/dt = 
!            (1/V') d/drho V' <|grad rho|> n chi d/drho(T)] + ...
!
!     in which (grad rho) appears to the first power, rather than the second,
!     and chi is the thermal diffusivity from a general geometry theory.
!
!     This is arguably the better way to proceed, since 
!
!         Vprime <|grad rho|> = A
!
!     where A is the surface area.  In this form, the quantity 
!
!          -n chi dT/drho 
!
!     can be identified as the heat flux per unit area, a natural 
!     quantity from a theoretical/simulation point of view.  This chi is 
!     the quantity returned by this subroutine.
!
!     If you are solving the transport equations in the form that
!     the ITER Expert Group agreed upon, e.g., 
!
!     (3/2) d(n T)/dt = 
!        (1/V') d/drho V' <|grad rho|**2> n chi d/drho(T)] + ...
!
!     then you need to multiply the chi_i and chi_e reported by this 
!     subroutine by the factor <|grad rho|>/<|grad rho|**2>.  This should
!     result in only small corrections to the predicted profiles.
!
! The choice of radial coordinate is more difficult to resolve:
!
!     We did not use the sqrt(toroidal flux) radial coordinate in our 
!     non-circular cross-section simulations.  Instead, we used "rho_d", 
!     where rho_d is defined to be the average horizontal minor radius 
!     at the elevation of the magnetic axis, normalized to the value of 
!     this quantity at the LCFS.
!
!     In other words, denote by "d" the horizontal diameter of a given flux 
!     surface measured at the elevation of the magnetic axis.  Denote by
!     "D" the horizontal diameter of the last closed flux surface at the 
!     elevation of the magnetic axis.  Then rho_d = d/D.  I believe this 
!     is variable number 67 (RMINOR) in the ITER Profile Database Standard 
!     List.
!
!     It is not difficult to allow for an arbitrary radial coordinate
!     in a transport code.  One must obtain all of the radial
!     quantities as functions of rho_d rather than rho via interpolation.
!
!     However, I do not expect everyone to go to this length to test our 
!     model, since you agreed to use the sqrt(toroidal flux) definition 
!     of the radial coordinate.  Thus, I suggest the following alternative:
!     Simply use the rho_d coordinate to define the scale lengths that 
!     appear in the formulas below.  For most quantities (such as R/LT), 
!     this simply amounts to including an additional factor d rho/d rho_d 
!     in the expressions passed to this subroutine.  While not completely 
!     correct, this workaround captures the dominant effect, related to 
!     evaluating the flux near the critical gradient.
!
! Summary of comments:
! -------------------
! (1) The general geometry extensions to the IFS/PPPL model were derived 
!     using rho = d/D = RMINOR as the radial coordinate.  To be 
!     most accurate, the transport equation should be solved using d/D 
!     as the radial coordinate.
!
!     If you use rho proportional to sqrt(toroidal flux) instead of rho=d/D
!     as your radial coordinate, you should at least carefully define the 
!     scale lengths as indicated below (using rho=d/D).
!     
! (2) This routine should be used to return the thermal transport
!     coefficients (chi_i, chi_e) for energy transport equations of the form
!
! (3/2) d(n T)/dt = (1/V') d/drho V' <|grad rho|> n chi d/drho(T)] + ...
!                                           
!     Note that <|grad rho|> only appears to the first power according
!     to this definition of chi.  If your code is hardwired to solve an
!     equation of the form 
!
! (3/2) d(n T)/dt = (1/V') d/drho V' <|grad rho|**2> n chi d/drho(T)] + ...
!
!     then multiply the chi_i and chi_e obtained from this routine by 
!     the factor <|grad rho|>/<|grad rho|**2>.
!     
!     This parameterization of chi is not complete.  There are significant
!     neglected physical processes that are known to be important in 
!     many operational regimes.  
!
!     The most significant problems are: 
!
!     (1) Trapped ion/long wavelength ITG modes.  These modes are known 
!     to be unstable for typical edge tokamak parameters.  However, until 
!     we have nonlinear estimates of the associated thermal diffusivity, 
!     these modes are ignored, leading to overly optimistic predictions of
!     edge thermal confinement.  
!     (2) Trapped electron modes, which can alter the stability boundary 
!     significantly for low collisionality.  At high collisionality these 
!     modes are generally stable and thus largely irrelevant.  When present, 
!     they are associated most strongly with particle transport, although 
!     there is also an associated heat transport.  
!     (3) Minority ion density gradients, which can strongly change chi 
!     and LT_crit. 
!     (4) Sheared flows, which are stabilizing.  This includes diamagnetic 
!     and ExB shear flows.
!     (5) Finite beta effects, generally stabilizing.
!-------------------------------------------------------------------------
! INPUT:
!
!   rlt:   R/L_Ti, where R = major radius and 1/L_Ti = -1/T_i dT_i/drho_d
!   rln:   R/L_ni, where 1/L_ni = -1/n_i dn_i/drho_d
!  rlne:   R/L_ne, where 1/L_ne = -1/n_e dn_e/drho_d
!     q:   The safety factor.
! kappa:   The elongation, defined here to be the kappa = maximum 
!           height/maximum diameter
!  shat:   Magnetic shear = rho_d/q dq/drho_d
!   zth:   Thermal Z_eff.  The simulations that were carried out to 
!           generate the formulae in this subroutine assumed the plasma 
!           was composed of a thermal hydrogenic species, thermal carbon, 
!           a hydrogenic beam species, and electrons.  We found that low-Z
!           impurities primarily act to dilute the main ion concentration, 
!           and can accounted for to first order by modifying the 
!           definition of Z_eff.  Some of the more important effects of 
!           the fast ions in the plasma are also partially accounted 
!           for by this parameter, which is:
!           zth == (n_i + 36 n_C)/(n_e - n_beam)
! nbeam:   Local fast ion (beam) density normalized to the electron
!          density.  
!   tau:   T_i/T_e.  Note that this is opposite to a widely used convention. 
!   eps:   rho_d/R, the local minor radius normalized to the major radius. 
!   gnu:   Dimensionless collisionality parameter.
!           gnu == rmajor*2.5e-7*n_e/(T_e**1.5*T_i**0.5), where n_e is in 
!            units of cm**-3, T_e and T_i are in eV, and rmajor is in units 
!            of m.  For an R = 2.4 m, 100 eV, 1e13 plasma, gnu=600.
! rmajor:  Major radius of the plasma
!  rho_i:  Local thermal gyroradius of thermal hydrogenic species.  
!   v_ti:  sqrt(T_i/m_i) where T_i and m_i are the local thermal hydrogenic 
!           temperature and average thermal hydrogenic mass.
!
! Units - The only dimensional parameters in the inputs are the major
!         radius, rho_i, and v_t.  Their units should be consistent; 
!         the chis that are returned will be in units of L*L*T, where 
!         L and T are the respective length and time units for the inputs.
!
! OUTPUT:
!
!  rltcrit:  R/L_Tcrit for ITG mode.
! rltcritz:  R/L_Tcrit for carbon branch.
!    chi_0:  Normalized chi (ignore).
!        g:  L_Tc/L_T, where L_Tc is the critical temperature gradient 
!             scale length for the deuterium branch of the ITG mode.
!    chi_i:  Anomalous ion thermal diffusivity from toroidal ITG mode.
!    chi_e:  Anomalous electron thermal diffusivity from toroidal ITG mode.
!-------------------------------------------------------------------------

subroutine ifs_pppl(rlt,&
     rln,&
     rlne,&
     q,&
     kappa,&
     shat,&
     zth,&
     nbeam,&
     tau,&
     eps,&
     gnu,&
     rmajor,&
     rho_i,&
     v_ti,&
     rltcrit,&
     rltcritz,&
     chi_0,&
     g,&
     chi_i,&
     chi_e) 

  implicit none

  ! INPUT
  real, intent(in) :: rlt,rln,rlne,q,kappa
  real, intent(in) :: shat,zth,nbeam,tau,eps,gnu
  real, intent(in) :: rmajor,rho_i,v_ti

  ! OUTPUT
  real, intent(inout) :: rltcrit,rltcritz,chi_0
  real, intent(inout) :: g,chi_i,chi_e

  ! INTERNAL
  real :: taub,nu,chi0
  real :: f_0,f_z,chiz,g_facz
  real :: c1,trln,trlne,tshat
  real :: chie1,chie2,g_fac1

10 format(t2,a,t10,1pe12.5)
  if (0==1) then
     print 10,'rlt',rlt
     print 10,'rln',rln
     print 10,'rlne',rlne
     print 10,'q',q
     print 10,'kappa',kappa
     print 10,'shat',shat
     print 10,'zth',zth
     print 10,'nbeam',nbeam
     print 10,'tau',tau
     print 10,'eps',eps
     print 10,'gnu',gnu
     print 10,'rmajor',rmajor
     print 10,'rho_i',rho_i
     print 10,'v_ti',v_ti
  endif

  taub = tau/(1.0-nbeam)
  nu   = gnu*0.84

  trln  = min(abs(rln),6.0)*sign(1.0,rln)
  trlne = min(abs(rlne),6.0)

  ! Correct low-shear limit (s<0.5) by enforing the 
  ! well-known (approximate) ITG symmetry about s=0.5.
  ! (J. Candy, 4 November 2007)

  if (shat > 0.5) then
     tshat = shat
  else
     tshat = 1.0-shat
  endif

  ! Critical ion temperature gradients:

  rltcrit = 2.46*(1.0+2.78/q**2)**0.26* &
       (zth/2.0)**0.7*taub**0.52 &
       *( (0.671+0.570*tshat-0.189*trln)**2 &
       +0.335*trln+0.392-0.779*tshat+0.210*tshat**2) &
       *(1.0-0.942*(2.95*eps**1.257/nu**0.235-0.2126) &
       *zth**0.516/tshat**0.671)

  rltcritz = 0.75*(1.0+taub)*(1.0+tshat)* &
       max(1.0,3.0-2.0*trlne/3.0)* &
       (1.0+6.0*max(0.0,2.9-zth))

  if (zth > 3.0) then
     c1 = (3.0/zth)**1.8
  else
     c1 = 1.0
  endif

  f_0 = 11.8*c1*q**1.13/(1.0+tshat**0.84)/taub**1.07 &
       *(1.0+6.72*eps/nu**0.26/q**0.96) &
       /(1.0+((kappa-1.0)*q/3.6)**2)

  f_z = 7.88/(1.0+tshat)*max(0.25,zth-3.0)/taub**0.8 &
       /(1.0+((kappa-1.0)*q/3.6)**2)

  chi0 = f_0*rho_i**2*v_ti/rmajor
  chiz = f_z*rho_i**2*v_ti/rmajor

  if (rlt-rltcrit > 0.0) then
     g_fac1 = min(sqrt(rlt-rltcrit),(rlt-rltcrit))
  else
     g_fac1 = 0.0
  endif

  if (rlt-rltcritz > 0.0) then
     g_facz = min(sqrt(rlt-rltcritz),(rlt-rltcritz))
  else
     g_facz = 0.0
  endif

  chi_i = max(chi0*g_fac1,chiz*g_facz)

  g     = rlt/rltcrit
  chi_0 = chi0*sqrt(abs(rltcrit))

  chie1 = chi0*g_fac1*1.44*tau**0.4*(q/tshat)**0.3*nu**0.14 &
       *max(0.16667,eps)
  chie1 = 0.5*chie1*(1.0+trlne/3.0)

  chie2 = 0.5*max(2.0,(1.0+rlne/3.0))
  chie2 = chie2*0.526*tau*nu**0.22
  chie2 = chie2*chiz*g_facz

  ! Correction for n_i/n_e and ratio of heat fluxes rather than chi_s:

  chi_e = max(chie1,chie2)*(7.0-zth)/6.0

  chi_i = max(0.0,chi_i)
  chi_e = max(0.0,chi_e)


  if (0==1) then
     print 10,'chi_i',chi_i
     print 10,'chi_e',chi_e
  endif

  return

end subroutine ifs_pppl
