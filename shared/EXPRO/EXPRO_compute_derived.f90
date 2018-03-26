!--------------------------------------------------------------
! EXPRO_compute_derived.f90
!
! PURPOSE:
!  Compute numerous derived quantities, including primary 
!  ion density *EXPRO_ni_new* that satisfies quasineutrality.
!--------------------------------------------------------------

subroutine EXPRO_compute_derived

  use EXPRO_interface
  use GEO_interface

  implicit none

  integer :: n
  integer :: i
  integer :: is

  real, parameter :: k  = 1.6022e-12 ! erg/eV
  real, parameter :: mp = 1.6726e-24 ! g
  real, parameter :: me = 0.00027230 ! m_ele/m_deuterium (dimensionless)
  real, parameter :: e  = 4.8032e-10 ! statcoul
  real, parameter :: c  = 2.9979e10  ! cm/s
  real, parameter :: pi = 3.1415926535897932

  real, dimension(:), allocatable :: rho
  real, dimension(:), allocatable :: dummy
  real, dimension(:), allocatable :: cc
  real, dimension(:), allocatable :: loglam

  real :: r_min
  real :: fa,fb

  !---------------------------------------------------------------------
  ! Sanity checks
  !
  if (EXPRO_ctrl_numeq_flag == 1 .and. EXPRO_nfourier == -1) then
     print '(a)','ERROR: (EXPRO) input.profiles.geo missing'
     stop
  endif
  if (EXPRO_ctrl_quasineutral_flag == -1) then
     print '(a)','ERROR: (EXPRO) EXPRO_ctrl_quasineutral_flag not set.'
     stop
  endif
  if (EXPRO_ctrl_numeq_flag == -1) then
     print '(a)','ERROR: (EXPRO) EXPRO_ctrl_numeq_flag not set.'
     stop
  endif
  if (EXPRO_ctrl_n_ion == -1) then
     print '(a)','ERROR: (EXPRO) EXPRO_ctrl_n_ion not set.'
     stop
  endif
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  ! Infer orientation
  ! 
  EXPRO_signb = nint(EXPRO_b_ref/abs(EXPRO_b_ref))
  EXPRO_signq = nint(EXPRO_q(1)/abs(EXPRO_q(1)))
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Derived quantities:
  !
  allocate(rho(EXPRO_n_exp))
  allocate(dummy(EXPRO_n_exp))
  
  rho(:) = EXPRO_arho*EXPRO_rho(:)

  ! b_unit
  call bound_deriv(dummy,rho**2,EXPRO_rmin**2,EXPRO_n_exp)
  EXPRO_bunit(:) = EXPRO_b_ref*dummy(:)

  ! s
  call bound_deriv(dummy,EXPRO_q,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_s(:) = (EXPRO_rmin(:)/EXPRO_q(:))*dummy(:)

  !         d(rmaj)
  ! drmaj = -------
  !           dr
  call bound_deriv(EXPRO_drmaj,EXPRO_rmaj,EXPRO_rmin,EXPRO_n_exp)

  !         d(zmag)
  ! dzmag = -------
  !           dr
  call bound_deriv(EXPRO_dzmag,EXPRO_zmag,EXPRO_rmin,EXPRO_n_exp)

  !             r   d(kappa)
  ! s_kappa = ----- -------- 
  !           kappa    dr
  call bound_deriv(dummy,EXPRO_kappa,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_skappa(:) = (EXPRO_rmin(:)/EXPRO_kappa(:))*dummy(:)

  !             d(delta)
  ! s_delta = r -------- 
  !                dr
  call bound_deriv(dummy,EXPRO_delta,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_sdelta(:) = EXPRO_rmin(:)*dummy(:) 

  !            d(zeta)
  ! s_zeta = r -------- 
  !              dr
  call bound_deriv(dummy,EXPRO_zeta,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_szeta(:) = EXPRO_rmin(:)*dummy(:) 

  ! 1/L_ne = -dln(ne)/dr (1/m)
  call bound_deriv(EXPRO_dlnnedr,-log(EXPRO_ne),EXPRO_rmin,EXPRO_n_exp)

  ! 1/L_Te = -dln(Te)/dr (1/m)
  call bound_deriv(EXPRO_dlntedr,-log(EXPRO_te),EXPRO_rmin,EXPRO_n_exp)

  ! sne = -ne''/ne (1/m^2)
  call bound_deriv(EXPRO_sdlnnedr,EXPRO_ne*EXPRO_dlnnedr,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_sdlnnedr = EXPRO_sdlnnedr/EXPRO_ne

  ! sTe = -Te''/Te (1/m^2)
  call bound_deriv(EXPRO_sdlntedr,EXPRO_te*EXPRO_dlntedr,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_sdlntedr = EXPRO_sdlntedr/EXPRO_te

  EXPRO_dlnnidr = 0.0
  EXPRO_dlntidr = 0.0
  EXPRO_sdlnnidr = 0.0
  EXPRO_sdlntidr = 0.0

  do is=1,EXPRO_n_ion
     if (minval(EXPRO_ni(is,:)) > 0.0) then
        ! 1/L_ni = -dln(ni)/dr (1/m)
        call bound_deriv(EXPRO_dlnnidr(is,:),-log(EXPRO_ni(is,:)),EXPRO_rmin,EXPRO_n_exp)

        ! 1/L_Ti = -dln(Ti)/dr (1/m)
        call bound_deriv(EXPRO_dlntidr(is,:),-log(EXPRO_ti(is,:)),EXPRO_rmin,EXPRO_n_exp)

        ! sni = -ni''/ni (1/m^2)
        call bound_deriv(EXPRO_sdlnnidr(is,:),EXPRO_ni(is,:)*EXPRO_dlnnidr(is,:),EXPRO_rmin,EXPRO_n_exp)
        EXPRO_sdlnnidr(is,:) = EXPRO_sdlnnidr(is,:)/EXPRO_ni(is,:)

        ! sTi = -Ti''/Ti (1/m^2)
        call bound_deriv(EXPRO_sdlntidr(is,:),EXPRO_ti(is,:)*EXPRO_dlntidr(is,:),EXPRO_rmin,EXPRO_n_exp)
        EXPRO_sdlntidr(is,:) = EXPRO_sdlntidr(is,:)/EXPRO_ti(is,:)
     endif
  enddo

  ! 1/L_Ptot = -dln(Ptot)/dr (1/m)
  if (minval(EXPRO_ptot) > 0.0) then
     call bound_deriv(EXPRO_dlnptotdr,-log(EXPRO_ptot),EXPRO_rmin,EXPRO_n_exp)
  else
     EXPRO_dlnptotdr = 0.0
  endif

  ! dr/d(rho)
  call bound_deriv(EXPRO_drdrho,EXPRO_rmin,rho,EXPRO_n_exp)
  !--------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Fourier coefficients for geometry (if they exist)
  !
  if (EXPRO_nfourier > 0) then
     do n=0,EXPRO_nfourier
        do i=1,4  

           ! aR_n = EXPRO_geo(1,n,:)
           ! bR_n = EXPRO_geo(2,n,:)
           ! aZ_n = EXPRO_geo(3,n,:)
           ! bZ_n = EXPRO_geo(4,n,:)

           ! d(aR_n)/dr, d(bR_n)/dr, d(aZ_n)/dr, d(bZ_n)/dr

           call bound_deriv(EXPRO_dgeo(i,n,:),EXPRO_geo(i,n,:),&
                EXPRO_rmin,EXPRO_n_exp)
        enddo
     enddo
  endif

  !-------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Geometry factors: 
  !
  ! - w0, w0p, vol, volp
  !
  GEO_nfourier_in = EXPRO_nfourier
  GEO_signb_in    = EXPRO_signb
  call GEO_alloc(1)

  r_min = EXPRO_rmin(EXPRO_n_exp)

  do i=2,EXPRO_n_exp

     ! Parameters to be passed to GEO library   
     !
     ! NOTE: dp/dr set to zero without loss of generality.
     ! 
     GEO_rmin_in      = EXPRO_rmin(i)/r_min
     GEO_rmaj_in      = EXPRO_rmaj(i)/r_min
     GEO_drmaj_in     = EXPRO_drmaj(i)
     GEO_zmag_in      = EXPRO_zmag(i)/r_min
     GEO_dzmag_in     = EXPRO_dzmag(i)
     GEO_q_in         = EXPRO_q(i)
     GEO_s_in         = EXPRO_s(i)
     GEO_kappa_in     = EXPRO_kappa(i)
     GEO_s_kappa_in   = EXPRO_skappa(i)
     GEO_delta_in     = EXPRO_delta(i)
     GEO_s_delta_in   = EXPRO_sdelta(i)
     GEO_zeta_in      = EXPRO_zeta(i)
     GEO_s_zeta_in    = EXPRO_szeta(i)
     GEO_beta_star_in = 0.0
     !
     if (EXPRO_ctrl_numeq_flag == 0) then
        ! Call GEO with model shape
        GEO_model_in = 0
        call GEO_do()
     else
        ! Call GEO with general (numerical) shape
        GEO_model_in = 1
        GEO_fourier_in(1:4,:) = EXPRO_geo(:,:,i)/r_min
        GEO_fourier_in(5:8,:) = EXPRO_dgeo(:,:,i)
        call GEO_do()
        if (minval(GEOV_jac_r) <= 0.0) then
           print '(a,i3,a)','WARNING: (EXPRO) Negative GEO Jacobian for i =',i,' in input.profiles'
        endif
     endif

     ! V and dV/dr
     EXPRO_volp(i) = GEO_volume_prime*r_min**2
     EXPRO_vol(i)  = GEO_volume*r_min**3

     ! |grad r| at theta=0
     EXPRO_grad_r0(i) = GEO_grad_r0

     ! <|grad r|> 
     EXPRO_ave_grad_r(i) = GEO_fluxsurfave_grad_r

     call GEO_interp(0.0)

     ! B_poloidal and B_toroidal [T] at theta=0
     EXPRO_bp0(i) = GEO_bp*EXPRO_bunit(i)
     EXPRO_bt0(i) = GEO_bt*EXPRO_bunit(i)

     EXPRO_thetascale(i) = GEO_thetascale

     ! Plasma current [A] I = (1/mu0) Int[Bp dl] 
     EXPRO_ip(i) = 7.958e5*(GEO_bl*r_min*EXPRO_bunit(i))
     
  enddo

  !--------------------------------------------------------------
  ! Extrapolate some quantities to axis:
  !
  call bound_extrap(fa,fb,EXPRO_grad_r0,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_grad_r0(1) = fa

  call bound_extrap(fa,fb,EXPRO_ave_grad_r,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_ave_grad_r(1) = fa

  call bound_extrap(fa,fb,EXPRO_bp0,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_bp0(1) = fa

  call bound_extrap(fa,fb,EXPRO_bt0,EXPRO_rmin,EXPRO_n_exp)
  EXPRO_bt0(1) = fa
  !
  ! Both V and dV/dr are zero on axis.
  !
  EXPRO_vol(1)  = 0.0
  EXPRO_volp(1) = 0.0  
  EXPRO_ip(1)   = 0.0
  EXPRO_thetascale(1) = EXPRO_thetascale(2)

  !--------------------------------------------------------------

  !-----------------------------------------------------------------
  ! CGS calculation of deuterium sound speed (cm/s) and 
  ! deuterium gyroradius (cm)
  !
  EXPRO_cs(:)   = sqrt( k*(1e3*EXPRO_te(:))/(2.0*mp) )   
  EXPRO_rhos(:) = EXPRO_cs(:)/(e*(1e4*EXPRO_bunit(:))/(2.0*mp*c))
  ! 
  ! Convert to m/s and m:
  !
  EXPRO_cs(:)   = EXPRO_cs(:)/100.0
  EXPRO_rhos(:) = EXPRO_rhos(:)/100.0
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Compute the electron-electron collision frequency (1/s)
  allocate(cc(EXPRO_n_exp))
  allocate(loglam(EXPRO_n_exp))
  cc(:) = sqrt(2.0) * pi * 1.6022**4 * 1.0 / (4.0 * pi * 8.8542)**2 &
          * 1.0 / (sqrt(3.3452) * 1602.2**1.5) * 1e9
  loglam(:) = 24.0 - log(sqrt(EXPRO_ne(:)*1e13)/(EXPRO_te(:)*1e3))
  EXPRO_nuee(:) = cc(:) * loglam(:) * EXPRO_ne(:) &
       / (sqrt(me) * EXPRO_te(:)**1.5)
  deallocate(cc)
  deallocate(loglam)
  !-----------------------------------------------------------------
  
  !--------------------------------------------------------------
  ! Compute w0p, gamma_e, gamma_p and mach:
  !
  call bound_deriv(EXPRO_w0p,EXPRO_w0,EXPRO_rmin,EXPRO_n_exp)
  !  
  EXPRO_gamma_e(:) = -EXPRO_rmin(:)/EXPRO_q(:)*EXPRO_w0p(:)
  EXPRO_gamma_p(:) = -EXPRO_rmaj(:)*EXPRO_w0p(:)
  EXPRO_mach(:)    = EXPRO_rmaj(:)*EXPRO_w0(:)/EXPRO_cs(:)
  !--------------------------------------------------------------

  ! Clean up
  call GEO_alloc(0)
  deallocate(rho)

  ! Density profile control

  if (EXPRO_ctrl_quasineutral_flag == 1) then

     EXPRO_ni_new(:) = 0.0
     do is=2,EXPRO_ctrl_n_ion
        EXPRO_ni_new(:) = EXPRO_ni_new(:)+EXPRO_ctrl_z(is)*EXPRO_ni(is,:)
     enddo
     EXPRO_ni_new(:) = (EXPRO_ne(:)-EXPRO_ni_new(:))/EXPRO_ctrl_z(1)

     ! 1/L_ni = -dln(ni)/dr (1/m)
     call bound_deriv(EXPRO_dlnnidr_new(:),-log(EXPRO_ni_new(:)),&
          EXPRO_rmin,EXPRO_n_exp)

     ! sni = -ni''/ni (1/m^2)
     call bound_deriv(EXPRO_sdlnnidr_new(:),EXPRO_ni_new(:)*EXPRO_dlnnidr_new(:),EXPRO_rmin,EXPRO_n_exp)
     EXPRO_sdlnnidr_new(:) = EXPRO_sdlnnidr_new(:)/EXPRO_ni_new(:)
     
     if (minval(EXPRO_ni_new(:)) <= 0.0) then
        EXPRO_error = 1
     endif

  else

     EXPRO_ni_new(:) = EXPRO_ni(1,:)
     EXPRO_dlnnidr_new(:) = EXPRO_dlnnidr(1,:)
     EXPRO_sdlnnidr_new(:) = EXPRO_sdlnnidr(1,:)

  endif

  deallocate(dummy)

  do is=1,EXPRO_ctrl_n_ion
     if (minval(EXPRO_ni(is,:)) <= 0.0) then
        EXPRO_error=1
     endif
  enddo
 
end subroutine EXPRO_compute_derived
 
