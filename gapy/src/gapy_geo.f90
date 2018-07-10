module geo
  
  double precision :: GEO_signb_in = 1.0 
  double precision :: GEO_rmin_in = 0.5
  double precision :: GEO_rmaj_in = 3.0
  double precision :: GEO_drmaj_in = 0.0
  double precision :: GEO_zmag_in = 0.0
  double precision :: GEO_dzmag_in = 0.0
  double precision :: GEO_q_in = 1.0
  double precision :: GEO_s_in = 2.0
  double precision :: GEO_kappa_in = 1.0
  double precision :: GEO_s_kappa_in = 0.0
  double precision :: GEO_delta_in = 0.0
  double precision :: GEO_s_delta_in = 0.0
  double precision :: GEO_zeta_in = 0.0
  double precision :: GEO_s_zeta_in = 0.0
  double precision :: GEO_beta_star_in = 0.0
  double precision :: GEO_beta_star_1_in = 0.0
  double precision :: GEO_beta_star_2_in = 0.0

  integer :: GEO_ntheta_in=1001
  integer :: GEO_nfourier_in=0
  integer :: GEO_model_in

  double precision, dimension(8,0:32) :: GEO_fourier_in

  ! Values interpolated at input vector locations
  
  double precision, dimension(:), allocatable :: GEO_b 
  double precision, dimension(:), allocatable :: GEO_dbdt
  double precision, dimension(:), allocatable :: GEO_dbdt2
  double precision, dimension(:), allocatable :: GEO_bp
  double precision, dimension(:), allocatable :: GEO_bt
  double precision, dimension(:), allocatable :: GEO_gsin
  double precision, dimension(:), allocatable :: GEO_gcos1
  double precision, dimension(:), allocatable :: GEO_gcos2
  double precision, dimension(:), allocatable :: GEO_g_theta
  double precision, dimension(:), allocatable :: GEO_grad_r
  double precision, dimension(:), allocatable :: GEO_gq
  double precision, dimension(:), allocatable :: GEO_captheta
  double precision, dimension(:), allocatable :: GEO_nu
  double precision, dimension(:), allocatable :: GEO_l_r
  double precision, dimension(:), allocatable :: GEO_l_t
  double precision, dimension(:), allocatable :: GEO_nsin
  double precision, dimension(:), allocatable :: GEO_usin
  double precision, dimension(:), allocatable :: GEO_ucos
  double precision, dimension(:), allocatable :: GEO_bigr
  double precision, dimension(:), allocatable :: GEO_bigr_r
  double precision, dimension(:), allocatable :: GEO_bigr_t
  double precision, dimension(:), allocatable :: GEO_theta_nc
  double precision, dimension(:), allocatable :: GEO_theta_s
  double precision, dimension(:), allocatable :: GEO_chi2

  ! Scalar-valued functions (i.e., flux-surface averages, etc)

  double precision :: GEO_f
  double precision :: GEO_ffprime
  double precision :: GEO_volume_prime
  double precision :: GEO_volume
  double precision :: GEO_fluxsurfave_grad_r2
  double precision :: GEO_fluxsurfave_grad_r
  double precision :: GEO_grad_r0
  double precision :: GEO_thetascale
  double precision :: GEO_bl

  ! INTERNAL vector-valued functions used for interpolation
  !
  ! Defined over interval theta=(-pi,pi)

  double precision, dimension(:), allocatable :: GEOV_theta
  double precision, dimension(:), allocatable :: GEOV_b 
  double precision, dimension(:), allocatable :: GEOV_dbdt
  double precision, dimension(:), allocatable :: GEOV_dbdt2
  double precision, dimension(:), allocatable :: GEOV_bp
  double precision, dimension(:), allocatable :: GEOV_bt
  double precision, dimension(:), allocatable :: GEOV_gsin
  double precision, dimension(:), allocatable :: GEOV_gcos1
  double precision, dimension(:), allocatable :: GEOV_gcos2
  double precision, dimension(:), allocatable :: GEOV_g_theta
  double precision, dimension(:), allocatable :: GEOV_jac_r
  double precision, dimension(:), allocatable :: GEOV_grad_r
  double precision, dimension(:), allocatable :: GEOV_gq
  double precision, dimension(:), allocatable :: GEOV_captheta
  double precision, dimension(:), allocatable :: GEOV_nu
  double precision, dimension(:), allocatable :: GEOV_l_r
  double precision, dimension(:), allocatable :: GEOV_l_t
  double precision, dimension(:), allocatable :: GEOV_nsin
  double precision, dimension(:), allocatable :: GEOV_usin
  double precision, dimension(:), allocatable :: GEOV_ucos
  double precision, dimension(:), allocatable :: GEOV_bigr
  double precision, dimension(:), allocatable :: GEOV_bigr_r
  double precision, dimension(:), allocatable :: GEOV_bigr_t
  double precision, dimension(:), allocatable :: GEOV_theta_nc
  double precision, dimension(:), allocatable :: GEOV_theta_s
  double precision, dimension(:), allocatable :: GEOV_chi2


contains

  ! ** Main callable routine **
  
  subroutine geo_interp(n,theta_in,new_flag)

    !-------------------------------------
    implicit none
    !
    integer, intent(in) :: n
    double precision, intent(in), dimension(n) :: theta_in
    logical, intent(in) :: new_flag
    double precision :: theta_0
    !
    integer :: n_theta
    integer :: i1
    integer :: i2
    integer :: itheta
    !
    double precision :: x0
    double precision :: x1
    double precision :: dx
    double precision :: z
    double precision, parameter :: pi=3.141592653589793
    double precision, parameter :: tol=1e-6
    !-------------------------------------

    if (allocated(geo_b)) call geo_salloc(n,0)

    if (new_flag) then
       if (allocated(geov_b)) call geo_alloc(0)
       call geo_alloc(1)
       call geo_do
    endif

    call geo_salloc(n,1)

    !----------------------------------------------------------
    ! If we are only using s-alpha, set functions now and exit
    !
    if (GEO_model_in == -1) then

       ! Theta-independent functions

       GEO_fluxsurfave_grad_r  = 1.0
       GEO_fluxsurfave_grad_r2 = 1.0
       GEO_grad_r0 = 1.0

       GEO_ffprime   = 0.0
       GEO_f         = GEO_rmaj_in

       GEO_volume       = 2*pi**2*GEO_rmin_in**2*GEO_rmaj_in
       GEO_volume_prime = 4*pi**2*GEO_rmin_in*GEO_rmaj_in

       ! Theta-dependent functions (some are set to zero for now)

       do itheta=1,n
          
          theta_0 = theta_in(itheta)

          if (abs(theta_0) > pi+tol) then
             print *,'ERROR in GEO: theta_0 out of bounds in GEO_interp.'
          endif

          GEO_b(itheta)     = 1.0/(1.0+GEO_rmin_in/GEO_rmaj_in*cos(theta_0))
          GEO_dbdt(itheta)  = (GEO_rmin_in/GEO_rmaj_in)*sin(theta_0)
          GEO_dbdt2(itheta) = 0.0 ! check with NEO usage
          GEO_bp(itheta) = GEO_b(itheta)*GEO_rmin_in/GEO_rmaj_in/GEO_q_in
          GEO_bt(itheta) = GEO_b(itheta)

          ! Added extra B here to make proper connection
          ! to s-alpha without having to artificially remove 
          ! the B in the denominator of the drift.

          GEO_gsin(itheta)   = sin(theta_0)*GEO_b(itheta) 
          GEO_gcos1(itheta)  = cos(theta_0)*GEO_b(itheta)
          GEO_gcos2(itheta)  = -GEO_rmaj_in*GEO_beta_star_in
          GEO_usin(itheta)   = sin(theta_0)*GEO_b(itheta)
          GEO_ucos(itheta)   = cos(theta_0)*GEO_b(itheta)

          GEO_g_theta(itheta) = 1.0
          GEO_grad_r(itheta)  = 1.0
          GEO_gq(itheta)      = 1.0
          GEO_captheta(itheta)  = GEO_s_in*theta_0-&
               GEO_q_in**2*GEO_rmaj_in*GEO_beta_star_in*sin(theta_0) 
          GEO_nu(itheta)     = -GEO_q_in*theta_0
          GEO_l_r(itheta) = 0.0
          GEO_l_t(itheta) = 0.0
          GEO_nsin(itheta) = 0.0 ! check with NEO usage
          GEO_bigr(itheta) = GEO_rmaj_in/GEO_b(itheta)
          GEO_bigr_r(itheta) = cos(theta_0)
          GEO_bigr_t(itheta) = -GEO_rmin_in*sin(theta_0)
          GEO_theta_nc(itheta) = theta_0

       enddo

    else

       !----------------------------------------------------------
       ! General case:
       !
       ! To illustrate what's happening, let assume:
       !  - n_theta = 3 
       !  - theta_0 = pi+eps
       !
       !                *              
       !  x------x------x
       ! -pi     0      pi

       n_theta = size(GEOV_theta)

       do itheta=1,n
          
          theta_0 = theta_in(itheta)

          ! n_theta = 3

          x0 = theta_0-GEOV_theta(1)

          ! x0 = 2*pi in this case

          dx = GEOV_theta(2)-GEOV_theta(1)

          ! dx = pi

          i1 = int(x0/dx)+1

          ! i1 = 3

          i2 = i1+1

          ! i2 = 4 (out of bounds)

          x1 = (i1-1)*dx 

          ! x1 = 2*pi

          z  = (x0-x1)/dx

          ! z = 0.0

          !---------------------------------------------
          ! Catch the error associated with the special 
          ! case documented above.
          !
          if (i2 > n_theta) then  
             i2 = n_theta
          endif
          !---------------------------------------------

          GEO_b(itheta)     = GEOV_b(i1)+(GEOV_b(i2)-GEOV_b(i1))*z
          GEO_dbdt(itheta)  = GEOV_dbdt(i1)+(GEOV_dbdt(i2)-GEOV_dbdt(i1))*z
          GEO_dbdt2(itheta) = GEOV_dbdt2(i1)+(GEOV_dbdt2(i2)-GEOV_dbdt2(i1))*z
          GEO_bp(itheta)    = GEOV_bp(i1)+(GEOV_bp(i2)-GEOV_bp(i1))*z
          GEO_bt(itheta)    = GEOV_bt(i1)+(GEOV_bt(i2)-GEOV_bt(i1))*z
          GEO_gsin(itheta)     = GEOV_gsin(i1)+(GEOV_gsin(i2)-GEOV_gsin(i1))*z
          GEO_gcos1(itheta)    = GEOV_gcos1(i1)+(GEOV_gcos1(i2)-GEOV_gcos1(i1))*z
          GEO_gcos2(itheta)    = GEOV_gcos2(i1)+(GEOV_gcos2(i2)-GEOV_gcos2(i1))*z
          GEO_g_theta(itheta)  = GEOV_g_theta(i1)+(GEOV_g_theta(i2)-GEOV_g_theta(i1))*z
          GEO_grad_r(itheta)   = GEOV_grad_r(i1)+(GEOV_grad_r(i2)-GEOV_grad_r(i1))*z
          GEO_gq(itheta)       = GEOV_gq(i1)+(GEOV_gq(i2)-GEOV_gq(i1))*z
          GEO_captheta(itheta) = GEOV_captheta(i1)+(GEOV_captheta(i2)-GEOV_captheta(i1))*z
          GEO_nu(itheta)       = GEOV_nu(i1)+(GEOV_nu(i2)-GEOV_nu(i1))*z
          GEO_bigr(itheta)     = GEOV_bigr(i1)+(GEOV_bigr(i2)-GEOV_bigr(i1))*z
          GEO_l_r(itheta)      = GEOV_l_r(i1)+(GEOV_l_r(i2)-GEOV_l_r(i1))*z
          GEO_l_t(itheta)      = GEOV_l_t(i1)+(GEOV_l_t(i2)-GEOV_l_t(i1))*z
          GEO_nsin(itheta)     = GEOV_nsin(i1)+(GEOV_nsin(i2)-GEOV_nsin(i1))*z
          GEO_usin(itheta)     = GEOV_usin(i1)+(GEOV_usin(i2)-GEOV_usin(i1))*z
          GEO_ucos(itheta)     = GEOV_ucos(i1)+(GEOV_ucos(i2)-GEOV_ucos(i1))*z
          GEO_bigr_r(itheta)   = GEOV_bigr_r(i1)+(GEOV_bigr_r(i2)-GEOV_bigr_r(i1))*z
          GEO_bigr_t(itheta)   = GEOV_bigr_t(i1)+(GEOV_bigr_t(i2)-GEOV_bigr_t(i1))*z
          GEO_theta_nc(itheta) = GEOV_theta_nc(i1)+(GEOV_theta_nc(i2)-GEOV_theta_nc(i1))*z
          GEO_theta_s(itheta)  = GEOV_theta_s(i1)+(GEOV_theta_s(i2)-GEOV_theta_s(i1))*z
          GEO_chi2(itheta)     = GEOV_chi2(i1)+(GEOV_chi2(i2)-GEOV_chi2(i1))*z

       enddo

    endif

  end subroutine geo_interp

  subroutine geo_do

    !-----------------------------------------------------------
    implicit none
    !
    integer :: n_theta
    integer :: ny
    integer :: i
    integer :: n
    !
    integer, dimension(:), allocatable :: ic
    !
    double precision :: theta
    double precision :: d_theta
    double precision :: a
    double precision :: a_t
    double precision :: a_tt
    double precision :: x
    double precision :: bigr_tt
    double precision :: bigz_tt
    double precision :: g_tt
    double precision :: f
    double precision :: f_prime
    double precision :: c 
    double precision :: pi_2
    double precision :: denom
    !
    double precision :: b1
    double precision :: b2
    double precision :: b3
    double precision :: b4
    double precision :: b5
    !
    double precision, dimension(:), allocatable :: bigz
    double precision, dimension(:), allocatable :: bigz_r
    double precision, dimension(:), allocatable :: bigz_t
    double precision, dimension(:), allocatable :: bigz_l
    double precision, dimension(:), allocatable :: bigr_l
    double precision, dimension(:), allocatable :: r_c
    double precision, dimension(:), allocatable :: r_sc
    double precision, dimension(:), allocatable :: dbdl
    double precision, dimension(:,:), allocatable :: e
    double precision, dimension(:,:), allocatable :: ei
    double precision, dimension(:), allocatable :: loop
    double precision, dimension(:), allocatable :: beta_star
    double precision, dimension(:), allocatable :: a_R,b_R,a_Z,b_Z
    double precision, dimension(:), allocatable :: a_Rp,b_Rp,a_Zp,b_Zp
    !
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! Check for missing value
    !
    if (abs(GEO_signb_in) < 1e-10) then
       print '(a)','ERROR: (GEO_do) Bad value for GEO_signb_in.'
       stop
    endif
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! If we are using the s-alpha model, just compute stuff 
    ! directly and exit:
    !
    if (GEO_model_in == -1) then
       return
    endif
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! Setup for case of general geomtry
    !
    ny = GEO_nfourier_in
    !
    allocate(a_R(0:ny))
    allocate(b_R(0:ny))
    allocate(a_Z(0:ny))
    allocate(b_Z(0:ny))
    allocate(a_Rp(0:ny))
    allocate(b_Rp(0:ny))
    allocate(a_Zp(0:ny))
    allocate(b_Zp(0:ny))
    !
    a_R(:)  = GEO_fourier_in(1,0:ny)
    b_R(:)  = GEO_fourier_in(2,0:ny)
    a_Z(:)  = GEO_fourier_in(3,0:ny)
    b_Z(:)  = GEO_fourier_in(4,0:ny)
    a_Rp(:) = GEO_fourier_in(5,0:ny)
    b_Rp(:) = GEO_fourier_in(6,0:ny)
    a_Zp(:) = GEO_fourier_in(7,0:ny)
    b_Zp(:) = GEO_fourier_in(8,0:ny)
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! Allocate internal variables
    !
    n_theta = GEO_ntheta_in
    !
    allocate(ic(2-n_theta:2*n_theta-2))
    !
    allocate(bigz(n_theta))
    allocate(bigz_r(n_theta))
    allocate(bigz_t(n_theta))
    allocate(bigz_l(n_theta))
    allocate(bigr_l(n_theta))
    allocate(r_c(n_theta))
    allocate(r_sc(n_theta))
    allocate(dbdl(n_theta))
    allocate(e(n_theta,4))
    allocate(ei(n_theta,4))
    allocate(loop(4))
    allocate(beta_star(n_theta))
    !-----------------------------------------------------------

    pi_2 = 8.0*atan(1.0)

    do i=1,n_theta-1
       ic(i) = i
       ic(i-(n_theta-1)) = i
       ic(i+(n_theta-1)) = i
    enddo

    d_theta = pi_2/(n_theta-1)

    if (abs(GEO_delta_in) > 1.0) then
       print '(a)','ERROR: (GEO_do) |GEO_delta_in| > 1.'
       stop
    endif

    !------------------------------------------------------------------
    ! Compute fundamental geometric quantities (basic derivatives
    ! like dl/dt) and metric elements (g_tt).
    !
    do i=1,n_theta

       theta = -0.5*pi_2+(i-1)*d_theta

       GEOV_theta(i) = theta

       if (GEO_model_in == 0) then

          !-----------------------------------------
          ! Generalized Miller-type parameterization
          !-----------------------------------------

          x = asin(GEO_delta_in)

          ! A
          ! dA/dtheta
          ! d^2A/dtheta^2
          a    = theta+x*sin(theta)
          a_t  = 1.0+x*cos(theta)
          a_tt = -x*sin(theta)

          ! R(theta)
          ! dR/dr
          ! dR/dtheta
          ! d^2R/dtheta^2
          GEOV_bigr(i) = GEO_rmaj_in+GEO_rmin_in*cos(a)
          GEOV_bigr_r(i) = GEO_drmaj_in+cos(a)-GEO_s_delta_in/cos(x)*sin(theta)*sin(a)
          GEOV_bigr_t(i) = -GEO_rmin_in*a_t*sin(a)
          bigr_tt = -GEO_rmin_in*a_t**2*cos(a)-GEO_rmin_in*a_tt*sin(a)

          !-----------------------------------------------------------

          ! A
          ! dA/dtheta
          ! d^2A/dtheta^2
          a    = theta+GEO_zeta_in*sin(2*theta)
          a_t  = 1.0+2*GEO_zeta_in*cos(2*theta)
          a_tt = -4*GEO_zeta_in*sin(2*theta)

          ! Z(theta)
          ! dZ/dr
          ! dZ/dtheta
          ! d^2Z/dtheta^2
          bigz(i)   = GEO_zmag_in+GEO_kappa_in*GEO_rmin_in*sin(a)
          bigz_r(i) = GEO_dzmag_in+GEO_kappa_in*(1.0+GEO_s_kappa_in)*sin(a)+&
               GEO_kappa_in*GEO_s_zeta_in*cos(a)*sin(2*theta)
          bigz_t(i) = GEO_kappa_in*GEO_rmin_in*cos(a)*a_t
          bigz_tt   = -GEO_kappa_in*GEO_rmin_in*sin(a)*a_t**2+&
               GEO_kappa_in*GEO_rmin_in*cos(a)*a_tt

       else

          !-----------------------------------------
          ! Fourier-expansion (completely general)
          !-----------------------------------------

          GEOV_bigr(i)   = 0.5*a_R(0)
          GEOV_bigr_r(i) = 0.5*a_Rp(0)
          GEOV_bigr_t(i) = 0.0
          bigr_tt   = 0.0
          do n=1,ny
             GEOV_bigr(i) = GEOV_bigr(i)+a_R(n)*cos(n*theta)+b_R(n)*sin(n*theta)        
             GEOV_bigr_r(i)  = GEOV_bigr_r(i)+a_Rp(n)*cos(n*theta)+b_Rp(n)*sin(n*theta)        
             GEOV_bigr_t(i)  = GEOV_bigr_t(i)-n*a_R(n)*sin(n*theta)+n*b_R(n)*cos(n*theta) 
             bigr_tt = bigr_tt-n*n*(a_R(n)*cos(n*theta)+b_R(n)*sin(n*theta)) 
          enddo

          bigz(i)  = 0.5*a_Z(0)
          bigz_r(i)= 0.5*a_Zp(0)
          bigz_t(i) = 0.0
          bigz_tt   = 0.0
          do n=1,ny
             bigz(i)   = bigz(i)+a_Z(n)*cos(n*theta)+b_Z(n)*sin(n*theta)        
             bigz_r(i) = bigz_r(i)+a_Zp(n)*cos(n*theta)+b_Zp(n)*sin(n*theta)        
             bigz_t(i) = bigz_t(i)-n*a_Z(n)*sin(n*theta)+n*b_Z(n)*cos(n*theta) 
             bigz_tt   = bigz_tt-n*n*(a_Z(n)*cos(n*theta)+b_Z(n)*sin(n*theta)) 
          enddo

       endif

       g_tt = GEOV_bigr_t(i)**2+bigz_t(i)**2

       GEOV_jac_r(i) = GEOV_bigr(i)*(GEOV_bigr_r(i)*bigz_t(i)-GEOV_bigr_t(i)*bigz_r(i))

       GEOV_grad_r(i) = GEOV_bigr(i)*sqrt(g_tt)/GEOV_jac_r(i)

       GEOV_l_t(i) = sqrt(g_tt)

       ! 1/(du/dl)
       r_c(i) = GEOV_l_t(i)**3/(GEOV_bigr_t(i)*bigz_tt-bigz_t(i)*bigr_tt)

       ! cos(u)
       bigz_l(i) = bigz_t(i)/GEOV_l_t(i)

       ! -sin(u)
       bigr_l(i) = GEOV_bigr_t(i)/GEOV_l_t(i)

       ! l_r = cos(u) dZ/dr - sin(u) dR/dr
       GEOV_l_r(i) = bigz_l(i)*bigz_r(i)+bigr_l(i)*GEOV_bigr_r(i)

       GEOV_nsin(i) = (GEOV_bigr_r(i)*GEOV_bigr_t(i)+bigz_r(i)*bigz_t(i))/GEOV_l_t(i)

       ! beta_star(theta)
       beta_star(i) = GEO_beta_star_in + &
            GEO_beta_star_1_in*(1.0-cos(theta)) + &
            GEO_beta_star_2_in*(1.0-cos(2*theta)) 

    enddo
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Loop integral (1 to n_theta-1) to compute f
    !
    c = 0.0
    do i=1,n_theta-1
       c = c+GEOV_l_t(i)/(GEOV_bigr(i)*GEOV_grad_r(i))
    enddo
    f = GEO_rmin_in/(c*d_theta/pi_2)
    !
    ! Loop integral to compute V'
    !
    c = 0.0
    do i=1,n_theta-1
       c = c+GEOV_l_t(i)*GEOV_bigr(i)/GEOV_grad_r(i)
    enddo
    GEO_volume_prime = pi_2*c*d_theta
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! bt (toroidal field, Bt) 
    ! bp (poloidal field, Bp) 
    ! b  (total field, B)
    !
    do i=1,n_theta
       GEOV_bt(i) = f/GEOV_bigr(i)
       GEOV_bp(i) = (GEO_rmin_in/GEO_q_in)*GEOV_grad_r(i)/GEOV_bigr(i)
       GEOV_b(i)  = GEO_signb_in*sqrt(GEOV_bt(i)**2+GEOV_bp(i)**2)
    enddo
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! dbdl  (db/dl)
    ! dbdt  (db/d(theta))
    ! dbdt2 (d^2b/dt^2)
    ! gsin  (generalized sine)
    ! gcos1 (generalized cosine)
    ! gcos2 
    ! g_theta (G_theta)
    ! gq      (G_q)
    !
    !
    ! Use 5-point stencils for derivatives.  We get poor 
    ! accuracy for db/dt without this.
    ! 
    do i=1,n_theta

       b5 = GEOV_b(ic(i+2))
       b4 = GEOV_b(ic(i+1))
       b3 = GEOV_b(ic(i))
       b2 = GEOV_b(ic(i-1))
       b1 = GEOV_b(ic(i-2))

       GEOV_dbdt(i)  = (-b5+8.0*b4-8.0*b2+b1)/(12.0*d_theta)
       dbdl(i)  = GEOV_dbdt(i)/GEOV_l_t(i)
       GEOV_dbdt2(i) = (-b5+16.0*b4-30.0*b3+16.0*b2-b1)/(12.0*d_theta**2)
       GEOV_gsin(i)  = GEOV_bt(i)*GEO_rmaj_in*dbdl(i)/GEOV_b(i)**2
       GEOV_gcos1(i) = (GEOV_bt(i)**2/GEOV_bigr(i)*bigz_l(i)+GEOV_bp(i)**2/r_c(i))*GEO_rmaj_in/GEOV_b(i)**2
       GEOV_gcos2(i) = 0.5*(GEO_rmaj_in/GEOV_b(i)**2)*GEOV_grad_r(i)*(-beta_star(i))

       GEOV_g_theta(i) = GEOV_bigr(i)*GEOV_b(i)*GEOV_l_t(i)/(GEO_rmin_in*GEO_rmaj_in*GEOV_grad_r(i))
       GEOV_gq(i)    = GEO_rmin_in*GEOV_b(i)/(GEO_q_in*GEOV_bigr(i)*GEOV_bp(i))

       GEOV_usin(i)  = -GEOV_bigr_t(i)/GEOV_l_t(i)
       GEOV_ucos(i)  = (GEOV_bt(i)/GEOV_b(i))*bigz_l(i)

    enddo
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Compute integrands for E1,E2,E3 and E4=nu
    !
    ! NOTE: E3 now contains beta_star(theta)
    !
    do i=1,n_theta
       c = d_theta*GEOV_l_t(i)/(GEOV_bigr(i)*GEOV_grad_r(i))
       ei(i,1) = c*2.0*GEOV_bt(i)/GEOV_bp(i)*(GEO_rmin_in/r_c(i)-GEO_rmin_in*bigz_l(i)/GEOV_bigr(i))
       ei(i,2) = c*GEOV_b(i)**2/GEOV_bp(i)**2
       ei(i,3) = c*GEOV_grad_r(i)*0.5/GEOV_bp(i)**2*(GEOV_bt(i)/GEOV_bp(i))*beta_star(i)     
       ei(i,4) = -c*GEOV_grad_r(i)*(GEOV_bt(i)/GEOV_bp(i))     
    enddo
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Compute integrals E1,E2,E3,E4=nu from integrands
    !
    e(n_theta/2+1,:) = 0.0
    do i=n_theta/2+2,n_theta
       e(i,:) = e(i-1,:)+0.5*(ei(i-1,:)+ei(i,:))
    enddo
    do i=n_theta/2,1,-1
       e(i,:) = e(i+1,:)-0.5*(ei(i+1,:)+ei(i,:))
    enddo
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Compute f_prime = df/d(psi) (units 1/length). 
    !
    ! (conceptually, ff_prime is determined from q and s).
    !
    loop(:) = e(n_theta,:)-e(1,:)
    !
    f_prime = (pi_2*GEO_q_in*GEO_s_in/GEO_rmin_in-loop(1)/GEO_rmin_in+loop(3))/loop(2)

    do i=1,n_theta
       GEOV_nu(i) = e(i,4)
       GEOV_captheta(i) = GEOV_bp(i)/GEOV_b(i)*GEOV_grad_r(i)*GEOV_bigr(i)* & 
            (e(i,1)/GEO_rmin_in+e(i,2)*f_prime-e(i,3))
    enddo
    !------------------------------------------------------------------

    !-----------------------------------------------------------
    ! Scalar variables contained in GEO_interface:
    !
    ! NOTE: Flux-surface average:
    !
    !        /
    !        | d(theta) G_theta/B f 
    !        /
    ! <f> = -----------------------
    !        /
    !        | d(theta) G_theta/B
    !        /

    !
    ! (1) Loop integrals (1 to n_theta-1) to compute 
    !     flux-surface averages:
    !
    ! Denominator:

    denom = sum(GEOV_g_theta(1:n_theta-1)/GEOV_b(1:n_theta-1))

    GEO_fluxsurfave_grad_r = sum(GEOV_grad_r(1:n_theta-1)*GEOV_g_theta(1:n_theta-1)/ &
         GEOV_b(1:n_theta-1))/denom

    GEO_fluxsurfave_grad_r2 = sum(GEOV_grad_r(1:n_theta-1)**2*GEOV_g_theta(1:n_theta-1)/ &
         GEOV_b(1:n_theta-1))/denom

    ! theta(i) = 0 for i = n_theta/2+1

    GEO_grad_r0 = GEOV_grad_r(n_theta/2+1)
    GEO_ffprime = f*f_prime
    GEO_f       = f
    !
    ! pre-factor of 0.5 comes from triangular element in phi-direction:
    ! dV = (0.5*R*dphi)*(R*dZ) 
    !
    GEO_volume = 0.5*pi_2*sum(bigz_t(1:n_theta-1)*GEOV_bigr(1:n_theta-1)**2)*d_theta
    GEO_bl     = sum(GEOV_l_t(:)*GEOV_bp(:))*d_theta
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! GS2/NCLASS angle
    !
    GEOV_theta_nc(1) = GEOV_theta(1)
    do i=2,n_theta
       GEOV_theta_nc(i) = GEOV_theta_nc(i-1)+0.5*(GEOV_g_theta(i)+GEOV_g_theta(i-1))*d_theta
    enddo
    GEOV_theta_nc(:) = -0.5*pi_2+pi_2*(0.5*pi_2+GEOV_theta_nc(:))/(0.5*pi_2+GEOV_theta_nc(n_theta))
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! Straight fieldline angle
    !
    GEOV_theta_s(1) = GEOV_theta(1)
    r_sc(:) = GEOV_jac_r(:)/GEOV_bigr(:)**2
    do i=2,n_theta
       GEOV_theta_s(i) = GEOV_theta_s(i-1)+0.5*(r_sc(i)+r_sc(i-1))*d_theta
    enddo
    GEOV_theta_s(:) = -0.5*pi_2+pi_2*(0.5*pi_2+GEOV_theta_s(:))/(0.5*pi_2+GEOV_theta_s(n_theta))
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! Poloidal scale length (useful for code resolution choice).
    !
    r_sc = 0.0
    do i=2,n_theta-1
       r_sc(i) = (GEOV_gsin(i+1)-GEOV_gsin(i-1))/(2*d_theta)
    enddo
    GEO_thetascale = maxval(abs(r_sc(:)))
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! chi2 for comparison with le3:
    ! 
    ! GEOV_chi2 -> 2*chi2/chi1^2 = (2/q) psi2/psi1^2 + s/r^2
    !
    GEOV_chi2(:) = (1.0/geo_q_in)*(GEOV_bp(:)*bigz_l(:)-GEOV_bp(:)*GEOV_bigr(:)/r_c(:)+ &
         0.5*GEO_q_in/GEO_rmin_in*GEO_beta_star_in*GEOV_bigr(:)**2 &
         -GEO_ffprime)/(GEOV_bigr(:)*GEOV_bp(:))**2 &
         +GEO_s_in/GEO_rmin_in**2
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! Deallocate internal variables
    !
    deallocate(ic)
    !
    deallocate(bigz)
    deallocate(bigz_r)
    deallocate(bigz_t)
    deallocate(bigz_l)
    deallocate(bigr_l)
    deallocate(r_c)
    deallocate(r_sc)
    deallocate(dbdl)
    deallocate(e)
    deallocate(ei)
    deallocate(loop)
    deallocate(beta_star)
    !
    deallocate(a_R)
    deallocate(b_R)
    deallocate(a_Z)
    deallocate(b_Z)
    deallocate(a_Rp)
    deallocate(b_Rp)
    deallocate(a_Zp)
    deallocate(b_Zp)
    !-----------------------------------------------------------

  end subroutine geo_do

  subroutine geo_alloc(flag)

    !-------------------------------------------
    implicit none
    !
    integer, intent(in) :: flag
    !-------------------------------------------

    if (flag == 1) then
       if (GEO_ntheta_in < 9) then 
          print *,'Need more points in GEO_alloc.' 
          stop
       endif
       if (modulo(GEO_ntheta_in,2) == 0) then 
          print *,'GEO_ntheta_in must be odd in GEO_alloc.' 
          stop
       endif
       allocate(GEOV_b(GEO_ntheta_in))
       allocate(GEOV_bt(GEO_ntheta_in))
       allocate(GEOV_bp(GEO_ntheta_in))
       allocate(GEOV_dbdt(GEO_ntheta_in))
       allocate(GEOV_dbdt2(GEO_ntheta_in))
       allocate(GEOV_gsin(GEO_ntheta_in))
       allocate(GEOV_gcos1(GEO_ntheta_in))
       allocate(GEOV_gcos2(GEO_ntheta_in))
       allocate(GEOV_grad_r(GEO_ntheta_in))
       allocate(GEOV_jac_r(GEO_ntheta_in))
       allocate(GEOV_g_theta(GEO_ntheta_in))
       allocate(GEOV_gq(GEO_ntheta_in))
       allocate(GEOV_captheta(GEO_ntheta_in))
       allocate(GEOV_nu(GEO_ntheta_in))
       allocate(GEOV_theta(GEO_ntheta_in))
       allocate(GEOV_l_r(GEO_ntheta_in))
       allocate(GEOV_l_t(GEO_ntheta_in))
       allocate(GEOV_nsin(GEO_ntheta_in))
       allocate(GEOV_usin(GEO_ntheta_in))
       allocate(GEOV_ucos(GEO_ntheta_in))
       allocate(GEOV_bigr(GEO_ntheta_in))
       allocate(GEOV_bigr_r(GEO_ntheta_in))
       allocate(GEOV_bigr_t(GEO_ntheta_in))
       allocate(GEOV_theta_nc(GEO_ntheta_in))
       allocate(GEOV_theta_s(GEO_ntheta_in))
       allocate(GEOV_chi2(GEO_ntheta_in))

    else
       
       deallocate(GEOV_b)
       deallocate(GEOV_bt)
       deallocate(GEOV_bp)
       deallocate(GEOV_dbdt)
       deallocate(GEOV_dbdt2)
       deallocate(GEOV_gsin)
       deallocate(GEOV_gcos1)
       deallocate(GEOV_gcos2)
       deallocate(GEOV_jac_r)
       deallocate(GEOV_grad_r)
       deallocate(GEOV_g_theta)
       deallocate(GEOV_gq)
       deallocate(GEOV_captheta)
       deallocate(GEOV_nu)
       deallocate(GEOV_theta)
       deallocate(GEOV_l_r)
       deallocate(GEOV_l_t)
       deallocate(GEOV_nsin)
       deallocate(GEOV_usin)
       deallocate(GEOV_ucos)
       deallocate(GEOV_bigr)
       deallocate(GEOV_bigr_r)
       deallocate(GEOV_bigr_t)
       deallocate(GEOV_theta_nc)
       deallocate(GEOV_theta_s)
       deallocate(GEOV_chi2)

    endif

  end subroutine geo_alloc

  subroutine geo_salloc(n,flag)

    !-------------------------------------------
    implicit none
    !
    integer, intent(in) :: n
    integer, intent(in) :: flag
    !-------------------------------------------

    if (flag == 1) then
       allocate(geo_b(n))
       allocate(geo_dbdt(n))
       allocate(geo_dbdt2(n))
       allocate(geo_bp(n))
       allocate(geo_bt(n))
       allocate(GEO_gsin(n))
       allocate(GEO_gcos1(n))
       allocate(GEO_gcos2(n))
       allocate(GEO_g_theta(n))
       allocate(GEO_grad_r(n))
       allocate(GEO_gq(n))
       allocate(GEO_captheta(n))
       allocate(GEO_nu(n))
       allocate(GEO_l_r(n))
       allocate(GEO_l_t(n))
       allocate(GEO_nsin(n))
       allocate(GEO_usin(n))
       allocate(GEO_ucos(n))
       allocate(GEO_bigr(n))
       allocate(GEO_bigr_r(n))
       allocate(GEO_bigr_t(n))
       allocate(GEO_theta_nc(n))
       allocate(GEO_theta_s(n))
       allocate(GEO_chi2(n))
    else
       deallocate(geo_b)
       deallocate(geo_dbdt)
       deallocate(geo_dbdt2)
       deallocate(geo_bp)
       deallocate(geo_bt)
       deallocate(GEO_gsin)
       deallocate(GEO_gcos1)
       deallocate(GEO_gcos2)
       deallocate(GEO_g_theta)
       deallocate(GEO_grad_r)
       deallocate(GEO_gq)
       deallocate(GEO_captheta)
       deallocate(GEO_nu)
       deallocate(GEO_l_r)
       deallocate(GEO_l_t)
       deallocate(GEO_nsin)
       deallocate(GEO_usin)
       deallocate(GEO_ucos)
       deallocate(GEO_bigr)
       deallocate(GEO_bigr_r)
       deallocate(GEO_bigr_t)
       deallocate(GEO_theta_nc)
       deallocate(GEO_theta_s)
       deallocate(GEO_chi2)
    endif

  end subroutine geo_salloc

end module geo
