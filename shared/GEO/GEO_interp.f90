subroutine GEO_interp(theta_0)

  use GEO_interface

  !-------------------------------------
  implicit none
  !
  real, intent(in) :: theta_0
  !
  integer :: n_theta
  integer :: i1
  integer :: i2
  !
  real :: x0
  real :: x1
  real :: dx
  real :: z
  real, parameter :: pi=3.141592653589793
  real, parameter :: tol=1e-6
  !-------------------------------------

  if (abs(theta_0) > pi+tol) then
     print *,'ERROR in GEO: theta_0 out of bounds in GEO_interp.'
  endif

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

     GEO_b      = 1.0/(1.0+GEO_rmin_in/GEO_rmaj_in*cos(theta_0))
     GEO_dbdt   = (GEO_rmin_in/GEO_rmaj_in)*sin(theta_0)
     GEO_dbdt2  = 0.0 ! check with NEO usage
     GEO_bp     = GEO_b*GEO_rmin_in/GEO_rmaj_in/GEO_q_in
     GEO_bt     = GEO_b

     ! Added extra B here to make proper connection
     ! to s-alpha without having to artificially remove 
     ! the B in the denominator of the drift.

     GEO_gsin   = sin(theta_0)*GEO_b 
     GEO_gcos1  = cos(theta_0)*GEO_b
     GEO_gcos2  = -GEO_rmaj_in*GEO_beta_star_in
     GEO_usin   = sin(theta_0)*GEO_b
     GEO_ucos   = cos(theta_0)*GEO_b

     GEO_g_theta = 1.0
     GEO_grad_r  = 1.0
     GEO_gq      = 1.0
     GEO_captheta  = GEO_s_in*theta_0-&
          GEO_q_in**2*GEO_rmaj_in*GEO_beta_star_in*sin(theta_0) 
     GEO_nu     = -GEO_q_in*theta_0
     GEO_l_r = 0.0
     GEO_l_t = 0.0
     GEO_nsin = 0.0 ! check with NEO usage
     GEO_bigr = GEO_rmaj_in/GEO_b
     GEO_bigr_r = cos(theta_0)
     GEO_bigr_t = -GEO_rmin_in*sin(theta_0)
     GEO_theta_nc = theta_0
     return
  endif
  !----------------------------------------------------------

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

  GEO_b        = GEOV_b(i1)+(GEOV_b(i2)-GEOV_b(i1))*z
  GEO_dbdt     = GEOV_dbdt(i1)+(GEOV_dbdt(i2)-GEOV_dbdt(i1))*z
  GEO_dbdt2    = GEOV_dbdt2(i1)+(GEOV_dbdt2(i2)-GEOV_dbdt2(i1))*z
  GEO_bp       = GEOV_bp(i1)+(GEOV_bp(i2)-GEOV_bp(i1))*z
  GEO_bt       = GEOV_bt(i1)+(GEOV_bt(i2)-GEOV_bt(i1))*z
  GEO_gsin     = GEOV_gsin(i1)+(GEOV_gsin(i2)-GEOV_gsin(i1))*z
  GEO_gcos1    = GEOV_gcos1(i1)+(GEOV_gcos1(i2)-GEOV_gcos1(i1))*z
  GEO_gcos2    = GEOV_gcos2(i1)+(GEOV_gcos2(i2)-GEOV_gcos2(i1))*z
  GEO_g_theta  = GEOV_g_theta(i1)+(GEOV_g_theta(i2)-GEOV_g_theta(i1))*z
  GEO_grad_r   = GEOV_grad_r(i1)+(GEOV_grad_r(i2)-GEOV_grad_r(i1))*z
  GEO_gq       = GEOV_gq(i1)+(GEOV_gq(i2)-GEOV_gq(i1))*z
  GEO_captheta = GEOV_captheta(i1)+(GEOV_captheta(i2)-GEOV_captheta(i1))*z
  GEO_nu       = GEOV_nu(i1)+(GEOV_nu(i2)-GEOV_nu(i1))*z
  GEO_bigr     = GEOV_bigr(i1)+(GEOV_bigr(i2)-GEOV_bigr(i1))*z
  GEO_l_r      = GEOV_l_r(i1)+(GEOV_l_r(i2)-GEOV_l_r(i1))*z
  GEO_l_t      = GEOV_l_t(i1)+(GEOV_l_t(i2)-GEOV_l_t(i1))*z
  GEO_nsin     = GEOV_nsin(i1)+(GEOV_nsin(i2)-GEOV_nsin(i1))*z
  GEO_usin     = GEOV_usin(i1)+(GEOV_usin(i2)-GEOV_usin(i1))*z
  GEO_ucos     = GEOV_ucos(i1)+(GEOV_ucos(i2)-GEOV_ucos(i1))*z
  GEO_bigr_r   = GEOV_bigr_r(i1)+(GEOV_bigr_r(i2)-GEOV_bigr_r(i1))*z
  GEO_bigr_t   = GEOV_bigr_t(i1)+(GEOV_bigr_t(i2)-GEOV_bigr_t(i1))*z
  GEO_theta_nc = GEOV_theta_nc(i1)+(GEOV_theta_nc(i2)-GEOV_theta_nc(i1))*z
  GEO_theta_s  = GEOV_theta_s(i1)+(GEOV_theta_s(i2)-GEOV_theta_s(i1))*z
  GEO_chi2     = GEOV_chi2(i1)+(GEOV_chi2(i2)-GEOV_chi2(i1))*z

end subroutine GEO_interp
