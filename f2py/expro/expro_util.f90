subroutine expro_compute_derived

  use expro
  use geo

  implicit none

  integer :: n
  integer :: i
  integer :: is
  integer :: nx

  double precision, parameter :: k  = 1.6022d-12 ! erg/eV
  double precision, parameter :: e  = 4.8032d-10 ! statcoul
  double precision, parameter :: c  = 2.9979d10  ! cm/s
  double precision, parameter :: pi = 3.1415926535897932d0
  double precision, parameter :: mu0 = 1.2566e-6 ! permeability [SI: H/m]

  double precision :: mp ! mass_deuterium/2.0 (g)

  double precision, dimension(:), allocatable :: torflux
  double precision, dimension(:), allocatable :: temp
  double precision, dimension(:), allocatable :: cc
  double precision, dimension(:), allocatable :: loglam

  double precision :: r_min
  double precision :: fa,fb
  double precision :: theta(1)
  double precision :: p_ave,ne_ave
  double precision :: eps0,m0,ipma
  double precision :: bt2_ave
  double precision :: bp2_ave

  mp = expro_mass_deuterium/2d0  ! mass_deuterium/2.0 (g)

  if (expro_ctrl_n_ion == -1) expro_ctrl_n_ion = expro_n_ion

  !---------------------------------------------------------------------
  ! Sanity checks
  !
  if (expro_ctrl_numeq_flag == 1 .and. expro_nfourier == -1) then
     print '(a)','ERROR: (expro) input.gacode.geo missing'
     stop
  endif
  if (expro_ctrl_quasineutral_flag == -1) then
     print '(a)','ERROR: (expro) expro_ctrl_quasineutral_flag not set.'
     stop
  endif
  if (expro_ctrl_numeq_flag == -1) then
     print '(a)','ERROR: (expro) expro_ctrl_numeq_flag not set.'
     stop
  endif
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Infer orientation
  ! 
  expro_signb = nint(expro_torfluxa/abs(expro_torfluxa))
  expro_signq = nint(expro_q(1)/abs(expro_q(1)))
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Derived quantities:
  !
  allocate(torflux(expro_n_exp))
  allocate(temp(expro_n_exp))

  torflux(:) = expro_torfluxa*expro_rho(:)**2

  ! b_unit
  call bound_deriv(expro_bunit,torflux,0.5d0*expro_rmin**2,expro_n_exp)

  ! s
  call bound_deriv(temp,log(abs(expro_q)),expro_rmin,expro_n_exp)
  expro_s(:) = expro_rmin(:)*temp(:)

  !         d(rmaj)
  ! drmaj = -------
  !           dr
  call bound_deriv(expro_drmaj,expro_rmaj,expro_rmin,expro_n_exp)

  !         d(zmag)
  ! dzmag = -------
  !           dr
  call bound_deriv(expro_dzmag,expro_zmag,expro_rmin,expro_n_exp)

  !             r   d(kappa)
  ! s_kappa = ----- -------- 
  !           kappa    dr
  call bound_deriv(temp,expro_kappa,expro_rmin,expro_n_exp)
  expro_skappa(:) = (expro_rmin(:)/expro_kappa(:))*temp(:)

  !             d(delta)
  ! s_delta = r -------- 
  !                dr

  call bound_deriv(temp,expro_delta,expro_rmin,expro_n_exp)
  expro_sdelta(:) = expro_rmin(:)*temp(:)

  !            d(zeta)
  ! s_zeta = r -------- 
  !              dr
  call bound_deriv(temp,expro_zeta,expro_rmin,expro_n_exp)
  expro_szeta(:) = expro_rmin(:)*temp(:) 

  !                d(cos/sinx)
  ! s_cos/sinx = r -------- 
  !                dr
  call bound_deriv(temp,expro_shape_cos0,expro_rmin,expro_n_exp)
  expro_shape_scos0(:) = expro_rmin(:)*temp(:)
  call bound_deriv(temp,expro_shape_cos1,expro_rmin,expro_n_exp)
  expro_shape_scos1(:) = expro_rmin(:)*temp(:)
  call bound_deriv(temp,expro_shape_cos2,expro_rmin,expro_n_exp)
  expro_shape_scos2(:) = expro_rmin(:)*temp(:)
  call bound_deriv(temp,expro_shape_cos3,expro_rmin,expro_n_exp)
  expro_shape_scos3(:) = expro_rmin(:)*temp(:)
  call bound_deriv(temp,expro_shape_cos4,expro_rmin,expro_n_exp)
  expro_shape_scos4(:) = expro_rmin(:)*temp(:)
  call bound_deriv(temp,expro_shape_cos5,expro_rmin,expro_n_exp)
  expro_shape_scos5(:) = expro_rmin(:)*temp(:)
  call bound_deriv(temp,expro_shape_cos6,expro_rmin,expro_n_exp)
  expro_shape_scos6(:) = expro_rmin(:)*temp(:)
  call bound_deriv(temp,expro_shape_sin3,expro_rmin,expro_n_exp)
  expro_shape_ssin3(:) = expro_rmin(:)*temp(:)
  call bound_deriv(temp,expro_shape_sin4,expro_rmin,expro_n_exp)
  expro_shape_ssin4(:) = expro_rmin(:)*temp(:)
  call bound_deriv(temp,expro_shape_sin5,expro_rmin,expro_n_exp)
  expro_shape_ssin5(:) = expro_rmin(:)*temp(:)
  call bound_deriv(temp,expro_shape_sin6,expro_rmin,expro_n_exp)
  expro_shape_ssin6(:) = expro_rmin(:)*temp(:)

  ! 1/L_ne = -dln(ne)/dr (1/m)
  call bound_deriv(expro_dlnnedr,-log(expro_ne),expro_rmin,expro_n_exp)

  ! 1/L_Te = -dln(Te)/dr (1/m)
  call bound_deriv(expro_dlntedr,-log(expro_te),expro_rmin,expro_n_exp)

  ! NOTE: expro_sdln* will be renormalized after calculation of rhos later

  ! sne = -ne''/ne (1/m^2) [not fully normalized yet]
  call bound_deriv(expro_sdlnnedr,expro_ne*expro_dlnnedr,expro_rmin,expro_n_exp)
  expro_sdlnnedr = expro_sdlnnedr/expro_ne

  ! sTe = -Te''/Te (1/m^2) [not fully normalized yet]
  call bound_deriv(expro_sdlntedr,expro_te*expro_dlntedr,expro_rmin,expro_n_exp)
  expro_sdlntedr = expro_sdlntedr/expro_te

  expro_dlnnidr = 0d0
  expro_dlntidr = 0d0
  expro_sdlnnidr = 0d0
  expro_sdlntidr = 0d0

  do is=1,expro_n_ion
     if (minval(expro_ni(is,:)) > 0d0 .and. minval(expro_ti(is,:)) > 0d0) then
        ! 1/L_ni = -dln(ni)/dr (1/m)
        call bound_deriv(expro_dlnnidr(is,:),-log(expro_ni(is,:)),expro_rmin,expro_n_exp)

        ! 1/L_Ti = -dln(Ti)/dr (1/m)
        call bound_deriv(expro_dlntidr(is,:),-log(expro_ti(is,:)),expro_rmin,expro_n_exp)

        ! sni = -ni''/ni (1/m^2) [not fully normalized yet]
        call bound_deriv(expro_sdlnnidr(is,:),expro_ni(is,:)*expro_dlnnidr(is,:),expro_rmin,expro_n_exp)
        expro_sdlnnidr(is,:) = expro_sdlnnidr(is,:)/expro_ni(is,:)

        ! sTi = -Ti''/Ti (1/m^2) [not fully normalized yet]
        call bound_deriv(expro_sdlntidr(is,:),expro_ti(is,:)*expro_dlntidr(is,:),expro_rmin,expro_n_exp)
        expro_sdlntidr(is,:) = expro_sdlntidr(is,:)/expro_ti(is,:)
     endif
  enddo

  ! 1/L_Ptot = -dln(Ptot)/dr (1/m)
  if (minval(expro_ptot) > 0d0) then
     call bound_deriv(expro_dlnptotdr,-log(expro_ptot),expro_rmin,expro_n_exp)
  else
     expro_dlnptotdr = 0d0
  endif
  !--------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Fourier coefficients for geometry (if they exist)
  !
  if (expro_nfourier > 0) then
     do n=0,expro_nfourier
        do i=1,4  

           ! aR_n = expro_geo(1,n,:)
           ! bR_n = expro_geo(2,n,:)
           ! aZ_n = expro_geo(3,n,:)
           ! bZ_n = expro_geo(4,n,:)

           ! d(aR_n)/dr, d(bR_n)/dr, d(aZ_n)/dr, d(bZ_n)/dr

           call bound_deriv(expro_dgeo(i,n,:),expro_geo(i,n,:),&
                expro_rmin,expro_n_exp)
        enddo
     enddo
  endif

  !-------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Geometry factors: 
  !
  ! - w0, w0p, vol, volp
  !
  geo_nfourier_in = expro_nfourier
  geo_signb_in    = expro_signb

  ! r_min = a [m]
  r_min = expro_rmin(expro_n_exp)
  
  do i=2,expro_n_exp

     ! Parameters to be passed to geo library   
     !
     ! NOTE: dp/dr set to zero without loss of generality.
     ! 
     geo_rmin_in      = expro_rmin(i)/r_min
     geo_rmaj_in      = expro_rmaj(i)/r_min
     geo_drmaj_in     = expro_drmaj(i)
     geo_zmag_in      = expro_zmag(i)/r_min
     geo_dzmag_in     = expro_dzmag(i)
     geo_q_in         = expro_q(i)
     geo_s_in         = expro_s(i)
     geo_kappa_in     = expro_kappa(i)
     geo_s_kappa_in   = expro_skappa(i)
     geo_delta_in     = expro_delta(i)
     geo_s_delta_in   = expro_sdelta(i)
     geo_zeta_in      = expro_zeta(i)
     geo_s_zeta_in    = expro_szeta(i)
     geo_shape_cos0_in   = expro_shape_cos0(i)
     geo_shape_s_cos0_in = expro_shape_scos0(i)
     geo_shape_cos1_in   = expro_shape_cos1(i)
     geo_shape_s_cos1_in = expro_shape_scos1(i)
     geo_shape_cos2_in   = expro_shape_cos2(i)
     geo_shape_s_cos2_in = expro_shape_scos2(i)
     geo_shape_cos3_in   = expro_shape_cos3(i)
     geo_shape_s_cos3_in = expro_shape_scos3(i)
     geo_shape_cos4_in   = expro_shape_cos4(i)
     geo_shape_s_cos4_in = expro_shape_scos4(i)
     geo_shape_cos5_in   = expro_shape_cos5(i)
     geo_shape_s_cos5_in = expro_shape_scos5(i)
     geo_shape_cos6_in   = expro_shape_cos6(i)
     geo_shape_s_cos6_in = expro_shape_scos6(i)
     geo_shape_sin3_in   = expro_shape_sin3(i)
     geo_shape_s_sin3_in = expro_shape_ssin3(i)
     geo_shape_sin4_in   = expro_shape_sin4(i)
     geo_shape_s_sin4_in = expro_shape_ssin4(i)
     geo_shape_sin5_in   = expro_shape_sin5(i)
     geo_shape_s_sin5_in = expro_shape_ssin5(i)
     geo_shape_sin6_in   = expro_shape_sin6(i)
     geo_shape_s_sin6_in = expro_shape_ssin6(i)
     geo_beta_star_in = 0d0
     !
     theta(1) = 0d0
     if (expro_ctrl_numeq_flag == 0) then
        ! Call geo with extended harmonic shape
        geo_model_in = 0
     else
        ! Call geo with Fourier shape
        geo_model_in = 1
        geo_fourier_in(1:4,0:geo_nfourier_in) = expro_geo(:,:,i)/r_min
        geo_fourier_in(5:8,0:geo_nfourier_in) = expro_dgeo(:,:,i)
     endif

     call geo_interp(1,theta,.true.)
     if (minval(geov_jac_r) <= 0d0) then
        print '(a,i3,a,f5.2)','WARNING: (expro_util) J < 0 for i=',i,' in input.gacode: ',minval(geov_jac_r)
        expro_jerr = 1
     endif

     ! V, dV/dr and S (note that S=dV/dr only in a circle)
     expro_volp(i) = geo_volume_prime*r_min**2
     expro_vol(i)  = geo_volume*r_min**3
     expro_surf(i) = geo_surf*r_min**2

     ! |grad r| at theta=0
     expro_grad_r0(i) = geo_grad_r0

     ! <|grad r|> 
     expro_ave_grad_r(i) = geo_fluxsurfave_grad_r

     ! B_poloidal and B_toroidal [T] at theta=0
     expro_bp0(i) = geo_bp(1)*expro_bunit(i)
     expro_bt0(i) = geo_bt(1)*expro_bunit(i)

     ! <B_p^2> and <Bt^2> for betap and betat calculations
     expro_bp2(i) = geo_fluxsurfave_bp2*expro_bunit(i)**2
     expro_bt2(i) = geo_fluxsurfave_bt2*expro_bunit(i)**2

     expro_thetascale(i) = geo_thetascale
  enddo

  !--------------------------------------------------------------
  ! Extrapolate some quantities to axis:
  !
  call bound_extrap(fa,fb,expro_grad_r0,expro_rmin,expro_n_exp)
  expro_grad_r0(1) = fa

  call bound_extrap(fa,fb,expro_ave_grad_r,expro_rmin,expro_n_exp)
  expro_ave_grad_r(1) = fa

  call bound_extrap(fa,fb,expro_bp0,expro_rmin,expro_n_exp)
  expro_bp0(1) = fa

  call bound_extrap(fa,fb,expro_bt0,expro_rmin,expro_n_exp)
  expro_bt0(1) = fa

  call bound_extrap(fa,fb,expro_fpol,expro_rmin,expro_n_exp)
  expro_fpol(1) = fa
  !
  ! Both V and dV/dr are zero on axis.
  !
  expro_vol(1)  = 0d0
  expro_volp(1) = 0d0  
  expro_thetascale(1) = expro_thetascale(2)

  !--------------------------------------------------------------

  !-----------------------------------------------------------------
  ! CGS calculation of deuterium sound speed (cm/s) and 
  ! deuterium gyroradius (cm)
  !
  expro_cs(:)   = sqrt( k*(1d3*expro_te(:))/(2.0*mp) )   
  expro_rhos(:) = expro_cs(:)/(e*(1d4*expro_bunit(:))/(2.0*mp*c))
  ! 
  ! Convert to m/s and m:
  !
  expro_cs(:)   = expro_cs(:)/1d2
  expro_rhos(:) = expro_rhos(:)/1d2
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Renormalize shearing parameters
  !
  ! sn = -n''/n*rhos (1/m) 
  ! sT = -T''/T*rhos (1/m)
  expro_sdlnnedr = expro_sdlnnedr*expro_rhos(:)
  expro_sdlntedr = expro_sdlntedr*expro_rhos(:)
  do is=1,expro_n_ion
     expro_sdlnnidr(is,:) = expro_sdlnnidr(is,:)*expro_rhos(:)
     expro_sdlntidr(is,:) = expro_sdlntidr(is,:)*expro_rhos(:)
  enddo
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Compute the electron-electron collision frequency (1/s)
  allocate(cc(expro_n_exp))
  allocate(loglam(expro_n_exp))
  cc(:) = sqrt(2.0) * pi * 1.6022**4 * 1.0 / (4.0 * pi * 8.8542)**2 &
       * 1.0 / (sqrt(3.3452) * 1602.2**1.5) * 1e9
  loglam(:) = 24.0 - log(sqrt(expro_ne(:)*1e13)/(expro_te(:)*1e3))
  expro_nuee(:) = cc(:) * loglam(:) * expro_ne(:) &
       / (sqrt(expro_masse/2.0) * expro_te(:)**1.5)
  deallocate(cc)
  deallocate(loglam)
  !-----------------------------------------------------------------

  !--------------------------------------------------------------
  ! Compute w0p, gamma_e, gamma_p and mach:
  !
  call bound_deriv(expro_w0p,expro_w0,expro_rmin,expro_n_exp)
  !  
  expro_gamma_e(:) = -expro_rmin(:)/expro_q(:)*expro_w0p(:)
  expro_gamma_p(:) = -expro_rmaj(:)*expro_w0p(:)
  expro_mach(:)    = expro_rmaj(:)*expro_w0(:)/expro_cs(:)
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Transport particle, momentum and energy sources
  !
  ! Ohmic electron power  
  temp = expro_qohme
  call volint(temp,expro_pow_e_ohmic,expro_n_exp)
  ! Total auxiliary electron power  
  temp = expro_qbeame+expro_qrfe+expro_qione
  call volint(temp,expro_pow_e_aux,expro_n_exp)
  ! Total electron power 
  temp = temp+expro_qohme-(expro_qbrem-expro_qsync-expro_qline)-expro_qei+expro_qfuse
  call volint(temp,expro_pow_e,expro_n_exp)
  ! Total auxiliary ion power 
  temp = expro_qbeami+expro_qrfi+expro_qioni+expro_qcxi
  call volint(temp,expro_pow_i_aux,expro_n_exp)
  ! Total ion power 
  temp = temp+expro_qei+expro_qfusi
  call volint(temp,expro_pow_i,expro_n_exp)

  ! Exchange power
  call volint(expro_qei,expro_pow_ei,expro_n_exp)

  ! Fusion power
  call volint(expro_qfuse,expro_pow_e_fus,expro_n_exp)
  call volint(expro_qfusi,expro_pow_i_fus,expro_n_exp)

  ! Radiated power (sink/negative)
  call volint(expro_qbrem,expro_pow_e_brem,expro_n_exp)
  call volint(expro_qsync,expro_pow_e_sync,expro_n_exp)
  call volint(expro_qline,expro_pow_e_line,expro_n_exp)

  ! Particle/momentum
  call volint(expro_qpar_beam,expro_flow_beam,expro_n_exp)
  call volint(expro_qpar_wall,expro_flow_wall,expro_n_exp)
  call volint(expro_qmom,expro_flow_mom,expro_n_exp)
  !--------------------------------------------------------------

  ! Clean up
  deallocate(torflux)

  ! Density profile control

  if (expro_ctrl_quasineutral_flag == 1) then

     ! Accumulate density offsets (dummy variables below)
     expro_ni_new(:) = 0d0
     expro_dlnnidr_new(:) = 0d0
     do is=2,expro_ctrl_n_ion
        expro_ni_new(:) = expro_ni_new(:)+expro_z(is)*expro_ni(is,:)
        expro_dlnnidr_new(:) = expro_dlnnidr_new(:)+expro_z(is)*expro_ni(is,:)*expro_dlnnidr(is,:)
     enddo

     ! New quasineutral ion 1 profiles:

     ! density  
     expro_ni_new(:) = (expro_ne(:)-expro_ni_new(:))/expro_z(1)
     ! gradient scale length 1/L_ni = -dln(ni)/dr (1/m)
     expro_dlnnidr_new(:) = (expro_ne(:)*expro_dlnnedr(:)-expro_dlnnidr_new(:))/(expro_ni_new(:)*expro_z(1))

     ! sni = -ni''/ni (1/m^2)
     call bound_deriv(expro_sdlnnidr_new(:),expro_ni_new(:)*expro_dlnnidr_new(:),&
          expro_rmin,expro_n_exp)
     expro_sdlnnidr_new(:) = expro_sdlnnidr_new(:)/expro_ni_new(:)*expro_rhos(:)

     if (minval(expro_ni_new(:)) <= 0d0) expro_error = 1

  else

     expro_ni_new(:) = expro_ni(1,:)
     expro_dlnnidr_new(:) = expro_dlnnidr(1,:)
     expro_sdlnnidr_new(:) = expro_sdlnnidr(1,:)

  endif

  deallocate(temp)

  do is=1,expro_ctrl_n_ion
     if (minval(expro_ni(is,:)) <= 0d0) expro_error=1
  enddo

  !-------------------------------------------------------------------------------
  ! Compute expro scalars

  nx = expro_n_exp

  ! Aspect ratio and mass (JC: need to generalize)
  eps0 = expro_rmin(nx)/expro_rmaj(nx)
  m0   = 2d0
  ipma = abs(expro_current)

  !  p_ave = <p>  in Pa 
  ! ne_ave = <ne> in 10^19/m^3
  p_ave = 0d0
  ! bt2_ave/bp2_ave = <B^2> in T^2 
  bt2_ave = 0d0 ; bp2_ave = 0d0 ; ne_ave = 0d0
  do i=2,nx
     ne_ave  = ne_ave+(expro_rmin(i)-expro_rmin(i-1))*(expro_ne(i)+expro_ne(i-1))/2
     p_ave   = p_ave+(expro_vol(i)-expro_vol(i-1))*(expro_ptot(i)+expro_ptot(i-1))/2
     bt2_ave = bt2_ave+(expro_vol(i)-expro_vol(i-1))*(expro_bt2(i)+expro_bt2(i-1))/2
     bp2_ave = bp2_ave+(expro_vol(i)-expro_vol(i-1))*(expro_bp2(i)+expro_bp2(i-1))/2     
  enddo
  ne_ave  = ne_ave/expro_rmin(nx)
  p_ave   = p_ave/expro_vol(nx)
  bt2_ave = bt2_ave/expro_vol(nx)
  bp2_ave = bp2_ave/expro_vol(nx)

  ! SI unit calculation (B [T] and p [Pa])
  ! NOTE: 1/beta = 1/betap + 1/betat 
  expro_betap = 2*mu0*p_ave/bp2_ave
  expro_betat = 2*mu0*p_ave/bt2_ave
  ! For beta_n, use EFIT normalization factor a[m]*Bt[T]/Ip[MA]=a*BCENTR/CURRENT
  expro_betan = 1/(1/expro_betap+1/expro_betat)*(r_min*expro_bcentr/ipma)
  ! Greenwald density (current [MA])
  expro_greenwald = ipma/(pi*r_min**2)
  ! Transport power (MW)
  expro_ptransp = expro_pow_e(nx)+expro_pow_i(nx)

  if (expro_ptransp > 1e-6) then
     ! tau = W[MJ]/P[MW]
     expro_tau = (1.5*p_ave*expro_vol(nx)*1d-6)/expro_ptransp

     expro_tau98y2 = 0.05621     * &
          ipma**0.93             * &
          sqrt(bt2_ave)**0.15    * &
          expro_ptransp**(-0.69) * &
          ne_ave**0.41           * &
          m0**0.19               * &
          expro_rmaj(nx)**1.97   * &
          eps0**0.59             * &
          expro_kappa(nx)**0.78
  endif

  ! Reset geo variables to standard values to prevent expro from leaving
  ! geo in a "mysterious" state

  geo_rmin_in      = 1d0
  geo_rmaj_in      = 0.5d0
  geo_drmaj_in     = 3d0
  geo_zmag_in      = 0d0
  geo_dzmag_in     = 0d0
  geo_q_in         = 2d0
  geo_s_in         = 1d0
  geo_kappa_in     = 1d0
  geo_s_kappa_in   = 0d0
  geo_delta_in     = 0d0
  geo_s_delta_in   = 0d0
  geo_zeta_in      = 0d0
  geo_s_zeta_in    = 0d0
  geo_shape_cos0_in   = 0d0
  geo_shape_s_cos0_in = 0d0
  geo_shape_cos1_in   = 0d0
  geo_shape_s_cos1_in = 0d0
  geo_shape_cos2_in   = 0d0
  geo_shape_s_cos2_in = 0d0
  geo_shape_cos3_in   = 0d0
  geo_shape_s_cos3_in = 0d0
  geo_shape_cos4_in   = 0d0
  geo_shape_s_cos4_in = 0d0
  geo_shape_cos5_in   = 0d0
  geo_shape_s_cos5_in = 0d0
  geo_shape_cos6_in   = 0d0
  geo_shape_s_cos6_in = 0d0
  geo_shape_sin3_in   = 0d0
  geo_shape_s_sin3_in = 0d0
  geo_shape_sin4_in   = 0d0
  geo_shape_s_sin4_in = 0d0
  geo_shape_sin5_in   = 0d0
  geo_shape_s_sin5_in = 0d0
  geo_shape_sin6_in   = 0d0
  geo_shape_s_sin6_in = 0d0
  geo_beta_star_in = 0d0

end subroutine expro_compute_derived

!-------------------------------------------------------
! bound_deriv.f90
!
! PURPOSE:
!  Compute the finite-difference derivative of a 
!  function (on a grid which may be unequally-spaced)
!  using the Lagrange interpolating polynomial.  
!-------------------------------------------------------

subroutine bound_deriv(df,f,r,n)

  implicit none

  integer, intent(in) :: n

  double precision, intent(inout), dimension(n) :: df
  double precision, intent(in), dimension(n) :: f
  double precision, intent(in),dimension(n) :: r

  double precision :: r1,r2,r3,ra
  double precision :: f1,f2,f3

  integer :: i

  do i=1,n

     if (i == 1) then

        ! Left boundary

        ra = r(1)

        r1 = r(1)
        r2 = r(2)
        r3 = r(3)
        f1 = f(1)
        f2 = f(2)
        f3 = f(3)

     else if (i == n) then

        ! Right boundary

        ra = r(n)

        r1 = r(n-2)
        r2 = r(n-1)
        r3 = r(n)
        f1 = f(n-2)
        f2 = f(n-1)
        f3 = f(n)

     else

        ! Interior

        ra = r(i)

        r1 = r(i-1)
        r2 = r(i)
        r3 = r(i+1)
        f1 = f(i-1)
        f2 = f(i)
        f3 = f(i+1)

     endif

     ! Derivative of Lagrange interpolating polynomial:

     df(i) = ((ra-r1)+(ra-r2))/(r3-r1)/(r3-r2)*f3 &
          + ((ra-r1)+(ra-r3))/(r2-r1)/(r2-r3)*f2 &
          + ((ra-r2)+(ra-r3))/(r1-r2)/(r1-r3)*f1

  enddo

end subroutine bound_deriv

!-------------------------------------------------------
! bound_extrap.f90
!
! PURPOSE:
!  Extrapolate to left and right boundaries. 
!-------------------------------------------------------

subroutine bound_extrap(fa,fb,f,r,n)

  implicit none

  integer, intent(in) :: n

  double precision, intent(inout) :: fa
  double precision, intent(inout) :: fb
  double precision, intent(in), dimension(n) :: f
  double precision, intent(in), dimension(n) :: r

  double precision :: r1,r2,r3
  double precision :: ra,rb
  double precision :: f1,f2,f3


  ! Left boundary

  ra = r(1)

  r2 = r(2)
  r3 = r(3)

  f2 = f(2)
  f3 = f(3)

  fa = (ra-r2)/(r3-r2)*f3+(r3-ra)/(r3-r2)*f2

  ! Right boundary

  rb = r(n)

  r1 = r(n-2)
  r2 = r(n-1)

  f1 = f(n-2)
  f2 = f(n-1)

  fb = (rb-r1)/(r2-r1)*f2+(r2-rb)/(r2-r1)*f1

end subroutine bound_extrap

subroutine bound_interp(xj,yj,nj,x,y,n)

  implicit none

  integer :: i,j
  integer, intent(in) :: nj,n

  double precision :: u,r1,r2,r3
  double precision, intent(in), dimension(nj) :: xj,yj
  double precision, intent(in), dimension(n) :: x
  double precision, intent(inout), dimension(n) :: y
  integer, dimension(n) :: ja,jb,jc

  j=1
  do i=1,n
10   if (x(i) >= xj(j) .and. x(i) <= xj(j+1)) then
        ja(i) = j
        jb(i) = j+1
        if (j-1 > 0) then
           jc(i) = j-1
        else
           jc(i) = j+2
        endif
     else
        j = j+1
        goto 10
     endif
  enddo

  do i=1,n
     u = x(i)
     r1 = xj(ja(i))
     r2 = xj(jb(i))
     r3 = xj(jc(i))
     
     y(i) = (u-r1)*(u-r2)/(r3-r1)/(r3-r2)*yj(jc(i)) &
          + (u-r1)*(u-r3)/(r2-r1)/(r2-r3)*yj(jb(i)) &
          + (u-r2)*(u-r3)/(r1-r2)/(r1-r3)*yj(ja(i))
     
  enddo
  
end subroutine bound_interp

!------------------------------------------------------------------------------
! write routines:
!  writes - scalar
!  writev - vector
!  writea - 2D array
!------------------------------------------------------------------------------

subroutine expro_writes(x,xs1,xs2)

  use expro, only : ident

  implicit none

  double precision, intent(in) :: x
  character(len=*), intent(in) :: xs1,xs2

  if (abs(x) > 1d-16) then
     write(1,'(4a)') ident//trim(xs1)//' | '//trim(xs2)
     write(1,10) x
  endif

10 format(1pe14.7)

end subroutine expro_writes

subroutine expro_writei(i,xs1)

  use expro, only : ident

  implicit none

  integer, intent(in) :: i
  character(len=*), intent(in) :: xs1

  if (i > 0) then
     write(1,'(a)') ident//trim(xs1)
     write(1,10) i
  endif

10 format(i0)

end subroutine expro_writei

subroutine expro_writev(x,n,xs1,xs2)

  use expro, only : ident

  implicit none

  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: x
  character(len=*), intent(in) :: xs1,xs2
  integer :: i

  if (sum(abs(x)) > 1d-16) then
     write(1,'(4a)') ident//trim(xs1)//' | '//trim(xs2)
     do i=1,n
        write(1,10) i,x(i)
     enddo
  endif

10 format(i3,1x,1pe14.7)

end subroutine expro_writev

subroutine expro_writea(x,m,n,xs1,xs2)

  use expro, only : ident

  implicit none

  integer, intent(in) :: m,n
  double precision, intent(in), dimension(m,n) :: x
  character(len=*), intent(in) :: xs1,xs2
  integer :: i

  if (sum(abs(x)) > 1d-16) then
     write(1,'(4a)') ident//trim(xs1)//' | '//trim(xs2)
     do i=1,n
        write(1,10) i,x(:,i)
     enddo
  endif

10 format(i3,1x,10(1pe14.7,1x))

end subroutine expro_writea

subroutine expro_skip_header(io)

  implicit none

  integer, intent(in) :: io
  character(len=1) :: ctemp

  do
     read(io,'(a)') ctemp     
     if (ctemp /= '#') exit
  enddo
  backspace io 

end subroutine expro_skip_header

subroutine volint(f,fdv,n)

  use expro, only : expro_vol
  
  implicit none

  integer :: i
  integer, intent(in) :: n
  double precision, intent(in) :: f(n)
  double precision, intent(out) :: fdv(n)

  fdv(1) = 0d0

  ! Integration is exact for constant f (density)
  do i=2,n
     fdv(i) = fdv(i-1)+0.5d0*(f(i)+f(i-1))*(expro_vol(i)-expro_vol(i-1))
  enddo

end subroutine volint
