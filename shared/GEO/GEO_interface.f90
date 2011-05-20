!-----------------------------------------------------
! GEO_interface.f90
!
! PURPOSE:
!  Interface to GEO shape functions
!
! USAGE:
!  
!  (1) Model shape 
!
!  GEO_model_in = 0
!  call GEO_alloc(1)
!  call GEO_do()
!  call GEO_interp(theta)
!  call GEO_alloc(0)
!
!  (2) General shape 
!
!  GEO_model_in = 1
!  GEO_nfourier_in = <dim>
!  call GEO_alloc(1)
!  call GEO_do()
!  call GEO_interp(theta)
!  call GEO_alloc(0)
!
! INPUT:
!
!  (1) Dimensions and control:
!
!  GEO_ntheta_in   : number of theta interpolation nodes. 
!  GEO_nfourier_in : number of gen. geo. expansion coeffs.
!  GEO_model_in    : geometry method (-1=s-alpha, 0=model, 1=general).
!  GEO_signb_in    : magnetic field orientation.
!
!  (2) Shape parameters
!
!  GEO_signb_in
!  GEO_rmin_in 
!  GEO_rmaj_in
!  GEO_drmaj_in
!  GEO_zmag_in
!  GEO_dzmag_in
!  GEO_q_in (can be negative)
!  GEO_s_in
!  GEO_kappa_in
!  GEO_s_kappa_in
!  GEO_delta_in
!  GEO_s_delta_in
!  GEO_zeta_in
!  GEO_s_zeta_in
!  GEO_beta_star_in -> beta_unit*dlnpdr 
!  GEO_fourier_in(8,GEO_nfourier_in)
!
! OUTPUT:
!
!  GEO_b        : B/B_unit (carries the sign of B_unit)
!  GEO_bp       : B_p/B_unit
!  GEO_bt       : B_t/B_unit
!  GEO_dbdt2    : d(GEO_b)/dtheta
!  GEO_dbdt2    : d^2(GEO_b)/dtheta^2
!  GEO_gsin     : generalized sine
!  GEO_gcos1    : generalized cosine
!  GEO_gcos2    : pressure-dependent part of gcos
!  GEO_g_theta  : G_theta
!  GEO_grad_r   : |Grad r|
!  GEO_gq       : G_q
!  GEO_captheta : Theta
!  GEO_nu       : eikonal function nu(r,theta)
!  GEO_l_t      : dl/dt
!  GEO_nsin     : neoclassical sine (Rr Rt - Zr Zt)/(dl/dt)
!  GEO_usin     : Coriolis drift sine
!  GEO_ucos     : Coriolis drift cosine
!  GEO_bigr     : R
!  GEO_bigr_r   : dR/dr
!  GEO_bigr_t   : dR/dtheta
!  GEO_theta_nc : theta_nc (constant fieldline angle)
!-----------------------------------------------------

module GEO_interface

  real :: GEO_signb_in=1.0 
  real :: GEO_rmin_in 
  real :: GEO_rmaj_in
  real :: GEO_drmaj_in
  real :: GEO_zmag_in
  real :: GEO_dzmag_in
  real :: GEO_q_in
  real :: GEO_s_in
  real :: GEO_kappa_in
  real :: GEO_s_kappa_in
  real :: GEO_delta_in
  real :: GEO_s_delta_in
  real :: GEO_zeta_in
  real :: GEO_s_zeta_in
  real :: GEO_beta_star_in

  integer :: GEO_ntheta_in=1001
  integer :: GEO_nfourier_in=0
  integer :: GEO_model_in

  real, dimension(:,:), allocatable :: GEO_fourier_in

  ! Vector-valued functions:
  !
  ! Defined over interval theta=(-pi,pi)

  real, dimension(:), allocatable :: GEOV_theta

  real, dimension(:), allocatable :: GEOV_b 
  real, dimension(:), allocatable :: GEOV_dbdt
  real, dimension(:), allocatable :: GEOV_dbdt2
  real, dimension(:), allocatable :: GEOV_bp
  real, dimension(:), allocatable :: GEOV_bt
  real, dimension(:), allocatable :: GEOV_gsin
  real, dimension(:), allocatable :: GEOV_gcos1
  real, dimension(:), allocatable :: GEOV_gcos2
  real, dimension(:), allocatable :: GEOV_g_theta
  real, dimension(:), allocatable :: GEOV_grad_r
  real, dimension(:), allocatable :: GEOV_gq
  real, dimension(:), allocatable :: GEOV_captheta
  real, dimension(:), allocatable :: GEOV_nu
  real, dimension(:), allocatable :: GEOV_l_t
  real, dimension(:), allocatable :: GEOV_nsin
  real, dimension(:), allocatable :: GEOV_usin
  real, dimension(:), allocatable :: GEOV_ucos
  real, dimension(:), allocatable :: GEOV_bigr
  real, dimension(:), allocatable :: GEOV_bigr_r
  real, dimension(:), allocatable :: GEOV_bigr_t
  real, dimension(:), allocatable :: GEOV_theta_nc

  ! Scalar-valued functions:
  !
  ! Evaluated at theta=theta_0

  real :: GEO_b 
  real :: GEO_dbdt
  real :: GEO_dbdt2
  real :: GEO_bp
  real :: GEO_bt
  real :: GEO_gsin
  real :: GEO_gcos1
  real :: GEO_gcos2
  real :: GEO_g_theta
  real :: GEO_grad_r
  real :: GEO_gq
  real :: GEO_captheta
  real :: GEO_nu
  real :: GEO_l_t
  real :: GEO_nsin
  real :: GEO_usin
  real :: GEO_ucos
  real :: GEO_bigr
  real :: GEO_bigr_r
  real :: GEO_bigr_t
  real :: GEO_theta_nc

  ! Flux-surface averages:

  real :: GEO_f
  real :: GEO_ffprime
  real :: GEO_beta_star
  real :: GEO_volume_prime
  real :: GEO_volume
  real :: GEO_fluxsurfave_grad_r2
  real :: GEO_fluxsurfave_grad_r
  real :: GEO_grad_r0

end module GEO_interface
