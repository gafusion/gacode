module geo
  
  double precision :: geo_signb_in = 1.0 
  double precision :: geo_rmin_in = 0.5
  double precision :: geo_rmaj_in = 3.0
  double precision :: geo_drmaj_in = 0.0
  double precision :: geo_zmag_in = 0.0
  double precision :: geo_dzmag_in = 0.0
  double precision :: geo_q_in = 2.0
  double precision :: geo_s_in = 1.0
  double precision :: geo_kappa_in = 1.0
  double precision :: geo_s_kappa_in = 0.0
  double precision :: geo_delta_in = 0.0
  double precision :: geo_s_delta_in = 0.0
  double precision :: geo_zeta_in = 0.0
  double precision :: geo_s_zeta_in = 0.0
  double precision :: geo_shape_cos0_in = 0.0
  double precision :: geo_shape_s_cos0_in = 0.0
  double precision :: geo_shape_cos1_in = 0.0
  double precision :: geo_shape_s_cos1_in = 0.0
  double precision :: geo_shape_cos2_in = 0.0
  double precision :: geo_shape_s_cos2_in = 0.0
  double precision :: geo_shape_cos3_in = 0.0
  double precision :: geo_shape_s_cos3_in = 0.0
  double precision :: geo_shape_sin3_in = 0.0
  double precision :: geo_shape_s_sin3_in = 0.0
  double precision :: geo_beta_star_in = 0.0
  double precision :: geo_beta_star_1_in = 0.0
  double precision :: geo_beta_star_2_in = 0.0

  integer :: geo_ntheta_in=1001
  integer :: geo_nfourier_in=0
  integer :: geo_model_in

  double precision, dimension(8,0:32) :: geo_fourier_in

  ! Values interpolated at input vector locations
  
  double precision, dimension(:), allocatable :: geo_b 
  double precision, dimension(:), allocatable :: geo_dbdt
  double precision, dimension(:), allocatable :: geo_dbdt2
  double precision, dimension(:), allocatable :: geo_bp
  double precision, dimension(:), allocatable :: geo_bt
  double precision, dimension(:), allocatable :: geo_gsin
  double precision, dimension(:), allocatable :: geo_gcos1
  double precision, dimension(:), allocatable :: geo_gcos2
  double precision, dimension(:), allocatable :: geo_g_theta
  double precision, dimension(:), allocatable :: geo_grad_r
  double precision, dimension(:), allocatable :: geo_gq
  double precision, dimension(:), allocatable :: geo_captheta
  double precision, dimension(:), allocatable :: geo_nu
  double precision, dimension(:), allocatable :: geo_l_r
  double precision, dimension(:), allocatable :: geo_l_t
  double precision, dimension(:), allocatable :: geo_nsin
  double precision, dimension(:), allocatable :: geo_usin
  double precision, dimension(:), allocatable :: geo_ucos
  double precision, dimension(:), allocatable :: geo_bigr
  double precision, dimension(:), allocatable :: geo_bigr_r
  double precision, dimension(:), allocatable :: geo_bigr_t
  double precision, dimension(:), allocatable :: geo_theta_nc
  double precision, dimension(:), allocatable :: geo_theta_s
  double precision, dimension(:), allocatable :: geo_chi2

  ! Scalar-valued functions (i.e., flux-surface averages, etc)

  double precision :: geo_f
  double precision :: geo_ffprime
  double precision :: geo_volume_prime
  double precision :: geo_volume
  double precision :: geo_surf
  double precision :: geo_fluxsurfave_grad_r2
  double precision :: geo_fluxsurfave_grad_r
  double precision :: geo_fluxsurfave_bt2
  double precision :: geo_fluxsurfave_bp2
  double precision :: geo_grad_r0
  double precision :: geo_thetascale
  double precision :: geo_bl

  ! INTERNAL vector-valued functions used for interpolation
  !
  ! Defined over interval theta=(-pi,pi)

  double precision, dimension(:), allocatable :: geov_theta
  double precision, dimension(:), allocatable :: geov_b 
  double precision, dimension(:), allocatable :: geov_dbdt
  double precision, dimension(:), allocatable :: geov_dbdt2
  double precision, dimension(:), allocatable :: geov_bp
  double precision, dimension(:), allocatable :: geov_bt
  double precision, dimension(:), allocatable :: geov_gsin
  double precision, dimension(:), allocatable :: geov_gcos1
  double precision, dimension(:), allocatable :: geov_gcos2
  double precision, dimension(:), allocatable :: geov_g_theta
  double precision, dimension(:), allocatable :: geov_jac_r
  double precision, dimension(:), allocatable :: geov_grad_r
  double precision, dimension(:), allocatable :: geov_gq
  double precision, dimension(:), allocatable :: geov_captheta
  double precision, dimension(:), allocatable :: geov_nu
  double precision, dimension(:), allocatable :: geov_l_r
  double precision, dimension(:), allocatable :: geov_l_t
  double precision, dimension(:), allocatable :: geov_nsin
  double precision, dimension(:), allocatable :: geov_usin
  double precision, dimension(:), allocatable :: geov_ucos
  double precision, dimension(:), allocatable :: geov_bigr
  double precision, dimension(:), allocatable :: geov_bigr_r
  double precision, dimension(:), allocatable :: geov_bigr_t
  double precision, dimension(:), allocatable :: geov_theta_nc
  double precision, dimension(:), allocatable :: geov_theta_s
  double precision, dimension(:), allocatable :: geov_chi2


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
    if (geo_model_in == -1) then

       ! Theta-independent functions

       geo_fluxsurfave_grad_r  = 1.0
       geo_fluxsurfave_grad_r2 = 1.0
       geo_grad_r0 = 1.0

       geo_ffprime   = 0.0
       geo_f         = geo_rmaj_in

       geo_volume       = 2*pi**2*geo_rmin_in**2*geo_rmaj_in
       geo_volume_prime = 4*pi**2*geo_rmin_in*geo_rmaj_in
       geo_surf         = geo_volume_prime

       ! Theta-dependent functions (some are set to zero for now)

       do itheta=1,n
          
          theta_0 = theta_in(itheta)

          if (abs(theta_0) > pi+tol) then
             print *,'ERROR in geo: theta_0 out of bounds in geo_interp.'
          endif

          geo_b(itheta)     = 1.0/(1.0+geo_rmin_in/geo_rmaj_in*cos(theta_0))
          geo_dbdt(itheta)  = (geo_rmin_in/geo_rmaj_in)*sin(theta_0)
          geo_dbdt2(itheta) = 0.0 ! check with NEO usage
          geo_bp(itheta) = geo_b(itheta)*geo_rmin_in/geo_rmaj_in/geo_q_in
          geo_bt(itheta) = geo_b(itheta)

          ! Added extra B here to make proper connection
          ! to s-alpha without having to artificially remove 
          ! the B in the denominator of the drift.

          geo_gsin(itheta)   = sin(theta_0)*geo_b(itheta) 
          geo_gcos1(itheta)  = cos(theta_0)*geo_b(itheta)
          geo_gcos2(itheta)  = -geo_rmaj_in*geo_beta_star_in
          geo_usin(itheta)   = sin(theta_0)*geo_b(itheta)
          geo_ucos(itheta)   = cos(theta_0)*geo_b(itheta)

          geo_g_theta(itheta) = 1.0
          geo_grad_r(itheta)  = 1.0
          geo_gq(itheta)      = 1.0
          geo_captheta(itheta)  = geo_s_in*theta_0-&
               geo_q_in**2*geo_rmaj_in*geo_beta_star_in*sin(theta_0) 
          geo_nu(itheta)     = -geo_q_in*theta_0
          geo_l_r(itheta) = 0.0
          geo_l_t(itheta) = 0.0
          geo_nsin(itheta) = 0.0 ! check with NEO usage
          geo_bigr(itheta) = geo_rmaj_in/geo_b(itheta)
          geo_bigr_r(itheta) = cos(theta_0)
          geo_bigr_t(itheta) = -geo_rmin_in*sin(theta_0)
          geo_theta_nc(itheta) = theta_0

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

       n_theta = size(geov_theta)

       do itheta=1,n
          
          theta_0 = theta_in(itheta)

          ! n_theta = 3

          x0 = theta_0-geov_theta(1)

          ! x0 = 2*pi in this case

          dx = geov_theta(2)-geov_theta(1)

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

          geo_b(itheta)     = geov_b(i1)+(geov_b(i2)-geov_b(i1))*z
          geo_dbdt(itheta)  = geov_dbdt(i1)+(geov_dbdt(i2)-geov_dbdt(i1))*z
          geo_dbdt2(itheta) = geov_dbdt2(i1)+(geov_dbdt2(i2)-geov_dbdt2(i1))*z
          geo_bp(itheta)    = geov_bp(i1)+(geov_bp(i2)-geov_bp(i1))*z
          geo_bt(itheta)    = geov_bt(i1)+(geov_bt(i2)-geov_bt(i1))*z
          geo_gsin(itheta)     = geov_gsin(i1)+(geov_gsin(i2)-geov_gsin(i1))*z
          geo_gcos1(itheta)    = geov_gcos1(i1)+(geov_gcos1(i2)-geov_gcos1(i1))*z
          geo_gcos2(itheta)    = geov_gcos2(i1)+(geov_gcos2(i2)-geov_gcos2(i1))*z
          geo_g_theta(itheta)  = geov_g_theta(i1)+(geov_g_theta(i2)-geov_g_theta(i1))*z
          geo_grad_r(itheta)   = geov_grad_r(i1)+(geov_grad_r(i2)-geov_grad_r(i1))*z
          geo_gq(itheta)       = geov_gq(i1)+(geov_gq(i2)-geov_gq(i1))*z
          geo_captheta(itheta) = geov_captheta(i1)+(geov_captheta(i2)-geov_captheta(i1))*z
          geo_nu(itheta)       = geov_nu(i1)+(geov_nu(i2)-geov_nu(i1))*z
          geo_bigr(itheta)     = geov_bigr(i1)+(geov_bigr(i2)-geov_bigr(i1))*z
          geo_l_r(itheta)      = geov_l_r(i1)+(geov_l_r(i2)-geov_l_r(i1))*z
          geo_l_t(itheta)      = geov_l_t(i1)+(geov_l_t(i2)-geov_l_t(i1))*z
          geo_nsin(itheta)     = geov_nsin(i1)+(geov_nsin(i2)-geov_nsin(i1))*z
          geo_usin(itheta)     = geov_usin(i1)+(geov_usin(i2)-geov_usin(i1))*z
          geo_ucos(itheta)     = geov_ucos(i1)+(geov_ucos(i2)-geov_ucos(i1))*z
          geo_bigr_r(itheta)   = geov_bigr_r(i1)+(geov_bigr_r(i2)-geov_bigr_r(i1))*z
          geo_bigr_t(itheta)   = geov_bigr_t(i1)+(geov_bigr_t(i2)-geov_bigr_t(i1))*z
          geo_theta_nc(itheta) = geov_theta_nc(i1)+(geov_theta_nc(i2)-geov_theta_nc(i1))*z
          geo_theta_s(itheta)  = geov_theta_s(i1)+(geov_theta_s(i2)-geov_theta_s(i1))*z
          geo_chi2(itheta)     = geov_chi2(i1)+(geov_chi2(i2)-geov_chi2(i1))*z

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
    if (abs(geo_signb_in) < 1e-10) then
       print '(a)','ERROR: (geo_do) Bad value for geo_signb_in.'
       stop
    endif
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! If we are using the s-alpha model, just compute stuff 
    ! directly and exit:
    !
    if (geo_model_in == -1) then
       return
    endif
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! Setup for case of general geomtry
    !
    ny = geo_nfourier_in
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
    a_R(:)  = geo_fourier_in(1,0:ny)
    b_R(:)  = geo_fourier_in(2,0:ny)
    a_Z(:)  = geo_fourier_in(3,0:ny)
    b_Z(:)  = geo_fourier_in(4,0:ny)
    a_Rp(:) = geo_fourier_in(5,0:ny)
    b_Rp(:) = geo_fourier_in(6,0:ny)
    a_Zp(:) = geo_fourier_in(7,0:ny)
    b_Zp(:) = geo_fourier_in(8,0:ny)
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! Allocate internal variables
    !
    n_theta = geo_ntheta_in
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

    if (abs(geo_delta_in) > 1.0) then
       print '(a)','ERROR: (geo_do) |geo_delta_in| > 1.'
       stop
    endif

    !------------------------------------------------------------------
    ! Compute fundamental geometric quantities (basic derivatives
    ! like dl/dt) and metric elements (g_tt).
    !
    do i=1,n_theta

       theta = -0.5*pi_2+(i-1)*d_theta

       geov_theta(i) = theta

       if (geo_model_in == 0) then

          !-----------------------------------------
          ! Generalized Miller-type parameterization
          !-----------------------------------------

          x = asin(geo_delta_in)

          ! A
          ! dA/dtheta
          ! d^2A/dtheta^2
          a    = theta &
               + geo_shape_cos0_in &
               + geo_shape_cos1_in*cos(theta) &
               + geo_shape_cos2_in*cos(2*theta) &
               + geo_shape_cos3_in*cos(3*theta) &
               + x*sin(theta) &
               - geo_zeta_in*sin(2*theta) &
               + geo_shape_sin3_in*sin(3*theta) 
          a_t  = 1.0 &
               - geo_shape_cos1_in*sin(theta) &
               - 2*geo_shape_cos2_in*sin(2*theta) &
               - 3*geo_shape_cos3_in*sin(3*theta) &
               + x*cos(theta) &
               - 2*geo_zeta_in*cos(2*theta) &
               + 3*geo_shape_sin3_in*cos(3*theta)
          a_tt = -geo_shape_cos1_in*cos(theta) &
               - 4*geo_shape_cos2_in*cos(2*theta) &
               - 9*geo_shape_cos3_in*cos(3*theta) &
               - x*sin(theta) &
               + 4*geo_zeta_in*sin(2*theta) &
               - 9*geo_shape_sin3_in*sin(3*theta)

          ! R(theta)
          ! dR/dr
          ! dR/dtheta
          ! d^2R/dtheta^2
          geov_bigr(i) = geo_rmaj_in + geo_rmin_in*cos(a)
          geov_bigr_r(i) = geo_drmaj_in + cos(a) &
               - sin(a) * &
               (geo_shape_s_cos0_in &
               + geo_shape_s_cos1_in*cos(theta) &
               + geo_shape_s_cos2_in*cos(2*theta) &
               + geo_shape_s_cos3_in*cos(3*theta) &
               + geo_s_delta_in/cos(x)*sin(theta) &
               - geo_s_zeta_in*sin(2*theta) &
               + geo_shape_s_sin3_in*sin(3*theta))
          geov_bigr_t(i) = -geo_rmin_in * a_t * sin(a)
          bigr_tt = &
               -geo_rmin_in*a_t**2*cos(a) &
               -geo_rmin_in*a_tt*sin(a)

          !-----------------------------------------------------------

          ! A
          ! dA/dtheta
          ! d^2A/dtheta^2
          a    = theta
          a_t  = 1.0
          a_tt = 0.0

          ! Z(theta)
          ! dZ/dr
          ! dZ/dtheta
          ! d^2Z/dtheta^2
          bigz(i)   = geo_zmag_in+geo_kappa_in*geo_rmin_in*sin(a)
          bigz_r(i) = geo_dzmag_in + geo_kappa_in*(1.0+geo_s_kappa_in)*sin(a) 
          bigz_t(i) = geo_kappa_in*geo_rmin_in*cos(a)*a_t
          bigz_tt   = -geo_kappa_in*geo_rmin_in*sin(a)*a_t**2+&
               geo_kappa_in*geo_rmin_in*cos(a)*a_tt

       else

          !-----------------------------------------
          ! Fourier-expansion (completely general)
          !-----------------------------------------

          geov_bigr(i)   = 0.5*a_R(0)
          geov_bigr_r(i) = 0.5*a_Rp(0)
          geov_bigr_t(i) = 0.0
          bigr_tt   = 0.0
          do n=1,ny
             geov_bigr(i) = geov_bigr(i)+a_R(n)*cos(n*theta)+b_R(n)*sin(n*theta)        
             geov_bigr_r(i)  = geov_bigr_r(i)+a_Rp(n)*cos(n*theta)+b_Rp(n)*sin(n*theta)        
             geov_bigr_t(i)  = geov_bigr_t(i)-n*a_R(n)*sin(n*theta)+n*b_R(n)*cos(n*theta) 
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

       g_tt = geov_bigr_t(i)**2+bigz_t(i)**2

       geov_jac_r(i) = geov_bigr(i)*(geov_bigr_r(i)*bigz_t(i)-geov_bigr_t(i)*bigz_r(i))

       geov_grad_r(i) = geov_bigr(i)*sqrt(g_tt)/geov_jac_r(i)

       geov_l_t(i) = sqrt(g_tt)

       ! 1/(du/dl)
       r_c(i) = geov_l_t(i)**3/(geov_bigr_t(i)*bigz_tt-bigz_t(i)*bigr_tt)

       ! cos(u)
       bigz_l(i) = bigz_t(i)/geov_l_t(i)

       ! -sin(u)
       bigr_l(i) = geov_bigr_t(i)/geov_l_t(i)

       ! l_r = cos(u) dZ/dr - sin(u) dR/dr
       geov_l_r(i) = bigz_l(i)*bigz_r(i)+bigr_l(i)*geov_bigr_r(i)

       geov_nsin(i) = (geov_bigr_r(i)*geov_bigr_t(i)+bigz_r(i)*bigz_t(i))/geov_l_t(i)

       ! beta_star(theta)
       beta_star(i) = geo_beta_star_in + &
            geo_beta_star_1_in*(1.0-cos(theta)) + &
            geo_beta_star_2_in*(1.0-cos(2*theta)) 

    enddo
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Loop integral (1 to n_theta-1) to compute f
    !
    c = 0.0
    do i=1,n_theta-1
       c = c+geov_l_t(i)/(geov_bigr(i)*geov_grad_r(i))
    enddo
    f = geo_rmin_in/(c*d_theta/pi_2)
    !
    ! Loop integral to compute V'
    !
    c = 0.0
    do i=1,n_theta-1
       c = c+geov_l_t(i)*geov_bigr(i)/geov_grad_r(i)
    enddo
    geo_volume_prime = pi_2*c*d_theta
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! bt (toroidal field, Bt) 
    ! bp (poloidal field, Bp) 
    ! b  (total field, B)
    !
    do i=1,n_theta
       geov_bt(i) = f/geov_bigr(i)
       geov_bp(i) = (geo_rmin_in/geo_q_in)*geov_grad_r(i)/geov_bigr(i)
       geov_b(i)  = geo_signb_in*sqrt(geov_bt(i)**2+geov_bp(i)**2)
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

       b5 = geov_b(ic(i+2))
       b4 = geov_b(ic(i+1))
       b3 = geov_b(ic(i))
       b2 = geov_b(ic(i-1))
       b1 = geov_b(ic(i-2))

       geov_dbdt(i)  = (-b5+8.0*b4-8.0*b2+b1)/(12.0*d_theta)
       dbdl(i)  = geov_dbdt(i)/geov_l_t(i)
       geov_dbdt2(i) = (-b5+16.0*b4-30.0*b3+16.0*b2-b1)/(12.0*d_theta**2)
       geov_gsin(i)  = geov_bt(i)*geo_rmaj_in*dbdl(i)/geov_b(i)**2
       geov_gcos1(i) = (geov_bt(i)**2/geov_bigr(i)*bigz_l(i)+geov_bp(i)**2/r_c(i))*geo_rmaj_in/geov_b(i)**2
       geov_gcos2(i) = 0.5*(geo_rmaj_in/geov_b(i)**2)*geov_grad_r(i)*(-beta_star(i))

       geov_g_theta(i) = geov_bigr(i)*geov_b(i)*geov_l_t(i)/(geo_rmin_in*geo_rmaj_in*geov_grad_r(i))
       geov_gq(i)    = geo_rmin_in*geov_b(i)/(geo_q_in*geov_bigr(i)*geov_bp(i))

       geov_usin(i)  = -geov_bigr_t(i)/geov_l_t(i)
       geov_ucos(i)  = (geov_bt(i)/geov_b(i))*bigz_l(i)

    enddo
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Compute integrands for E1,E2,E3 and E4=nu
    !
    ! NOTE: E3 now contains beta_star(theta)
    !
    do i=1,n_theta
       c = d_theta*geov_l_t(i)/(geov_bigr(i)*geov_grad_r(i))
       ei(i,1) = c*2.0*geov_bt(i)/geov_bp(i)*(geo_rmin_in/r_c(i)-geo_rmin_in*bigz_l(i)/geov_bigr(i))
       ei(i,2) = c*geov_b(i)**2/geov_bp(i)**2
       ei(i,3) = c*geov_grad_r(i)*0.5/geov_bp(i)**2*(geov_bt(i)/geov_bp(i))*beta_star(i)     
       ei(i,4) = -c*geov_grad_r(i)*(geov_bt(i)/geov_bp(i))     
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
    f_prime = (pi_2*geo_q_in*geo_s_in/geo_rmin_in-loop(1)/geo_rmin_in+loop(3))/loop(2)

    do i=1,n_theta
       geov_nu(i) = e(i,4)
       geov_captheta(i) = geov_bp(i)/geov_b(i)*geov_grad_r(i)*geov_bigr(i)* & 
            (e(i,1)/geo_rmin_in+e(i,2)*f_prime-e(i,3))
    enddo
    !------------------------------------------------------------------

    !-----------------------------------------------------------
    ! Scalar variables contained in geo_interface:
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

    denom = sum(geov_g_theta(1:n_theta-1)/geov_b(1:n_theta-1))

    geo_fluxsurfave_grad_r = sum(geov_grad_r(1:n_theta-1)*geov_g_theta(1:n_theta-1)/ &
         geov_b(1:n_theta-1))/denom

    geo_fluxsurfave_grad_r2 = sum(geov_grad_r(1:n_theta-1)**2*geov_g_theta(1:n_theta-1)/ &
         geov_b(1:n_theta-1))/denom

    geo_fluxsurfave_bp2 = sum(geov_bp(1:n_theta-1)**2*geov_g_theta(1:n_theta-1)/ &
         geov_b(1:n_theta-1))/denom

    geo_fluxsurfave_bt2 = sum(geov_bt(1:n_theta-1)**2*geov_g_theta(1:n_theta-1)/ &
         geov_b(1:n_theta-1))/denom

    ! theta(i) = 0 for i = n_theta/2+1

    geo_grad_r0 = geov_grad_r(n_theta/2+1)
    geo_ffprime = f*f_prime
    geo_f       = f
    !
    ! pre-factor of 0.5 in dV comes from triangular element in phi-direction:
    ! dV = (0.5*R*dphi)*(R*dZ) 
    ! dS = (R*dphi)*(dl)
    !
    geo_volume = 0.5*pi_2*sum(bigz_t(1:n_theta-1)*geov_bigr(1:n_theta-1)**2)*d_theta
    geo_surf   = pi_2*sum(geov_l_t(1:n_theta-1)*geov_bigr(1:n_theta-1))*d_theta
    geo_bl     = sum(geov_l_t(:)*geov_bp(:))*d_theta
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! GS2/NCLASS angle
    !
    geov_theta_nc(1) = geov_theta(1)
    do i=2,n_theta
       geov_theta_nc(i) = geov_theta_nc(i-1)+0.5*(geov_g_theta(i)+geov_g_theta(i-1))*d_theta
    enddo
    geov_theta_nc(:) = -0.5*pi_2+pi_2*(0.5*pi_2+geov_theta_nc(:))/(0.5*pi_2+geov_theta_nc(n_theta))
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! Straight fieldline angle
    !
    geov_theta_s(1) = geov_theta(1)
    r_sc(:) = geov_jac_r(:)/geov_bigr(:)**2
    do i=2,n_theta
       geov_theta_s(i) = geov_theta_s(i-1)+0.5*(r_sc(i)+r_sc(i-1))*d_theta
    enddo
    geov_theta_s(:) = -0.5*pi_2+pi_2*(0.5*pi_2+geov_theta_s(:))/(0.5*pi_2+geov_theta_s(n_theta))
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! Poloidal scale length (useful for code resolution choice).
    !
    r_sc = 0.0
    do i=2,n_theta-1
       r_sc(i) = (geov_gsin(i+1)-geov_gsin(i-1))/(2*d_theta)
    enddo
    geo_thetascale = maxval(abs(r_sc(:)))
    !-----------------------------------------------------------

    !-----------------------------------------------------------
    ! chi2 for comparison with le3:
    ! 
    ! geov_chi2 -> 2*chi2/chi1^2 = (2/q) psi2/psi1^2 + s/r^2
    !
    geov_chi2(:) = (1.0/geo_q_in)*(geov_bp(:)*bigz_l(:)-geov_bp(:)*geov_bigr(:)/r_c(:)+ &
         0.5*geo_q_in/geo_rmin_in*geo_beta_star_in*geov_bigr(:)**2 &
         -geo_ffprime)/(geov_bigr(:)*geov_bp(:))**2 &
         +geo_s_in/geo_rmin_in**2
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
       if (geo_ntheta_in < 9) then 
          print *,'Need more points in geo_alloc.' 
          stop
       endif
       if (modulo(geo_ntheta_in,2) == 0) then 
          print *,'geo_ntheta_in must be odd in geo_alloc.' 
          stop
       endif
       allocate(geov_b(geo_ntheta_in))
       allocate(geov_bt(geo_ntheta_in))
       allocate(geov_bp(geo_ntheta_in))
       allocate(geov_dbdt(geo_ntheta_in))
       allocate(geov_dbdt2(geo_ntheta_in))
       allocate(geov_gsin(geo_ntheta_in))
       allocate(geov_gcos1(geo_ntheta_in))
       allocate(geov_gcos2(geo_ntheta_in))
       allocate(geov_grad_r(geo_ntheta_in))
       allocate(geov_jac_r(geo_ntheta_in))
       allocate(geov_g_theta(geo_ntheta_in))
       allocate(geov_gq(geo_ntheta_in))
       allocate(geov_captheta(geo_ntheta_in))
       allocate(geov_nu(geo_ntheta_in))
       allocate(geov_theta(geo_ntheta_in))
       allocate(geov_l_r(geo_ntheta_in))
       allocate(geov_l_t(geo_ntheta_in))
       allocate(geov_nsin(geo_ntheta_in))
       allocate(geov_usin(geo_ntheta_in))
       allocate(geov_ucos(geo_ntheta_in))
       allocate(geov_bigr(geo_ntheta_in))
       allocate(geov_bigr_r(geo_ntheta_in))
       allocate(geov_bigr_t(geo_ntheta_in))
       allocate(geov_theta_nc(geo_ntheta_in))
       allocate(geov_theta_s(geo_ntheta_in))
       allocate(geov_chi2(geo_ntheta_in))

    else
       
       deallocate(geov_b)
       deallocate(geov_bt)
       deallocate(geov_bp)
       deallocate(geov_dbdt)
       deallocate(geov_dbdt2)
       deallocate(geov_gsin)
       deallocate(geov_gcos1)
       deallocate(geov_gcos2)
       deallocate(geov_jac_r)
       deallocate(geov_grad_r)
       deallocate(geov_g_theta)
       deallocate(geov_gq)
       deallocate(geov_captheta)
       deallocate(geov_nu)
       deallocate(geov_theta)
       deallocate(geov_l_r)
       deallocate(geov_l_t)
       deallocate(geov_nsin)
       deallocate(geov_usin)
       deallocate(geov_ucos)
       deallocate(geov_bigr)
       deallocate(geov_bigr_r)
       deallocate(geov_bigr_t)
       deallocate(geov_theta_nc)
       deallocate(geov_theta_s)
       deallocate(geov_chi2)

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
       allocate(geo_gsin(n))
       allocate(geo_gcos1(n))
       allocate(geo_gcos2(n))
       allocate(geo_g_theta(n))
       allocate(geo_grad_r(n))
       allocate(geo_gq(n))
       allocate(geo_captheta(n))
       allocate(geo_nu(n))
       allocate(geo_l_r(n))
       allocate(geo_l_t(n))
       allocate(geo_nsin(n))
       allocate(geo_usin(n))
       allocate(geo_ucos(n))
       allocate(geo_bigr(n))
       allocate(geo_bigr_r(n))
       allocate(geo_bigr_t(n))
       allocate(geo_theta_nc(n))
       allocate(geo_theta_s(n))
       allocate(geo_chi2(n))
    else
       deallocate(geo_b)
       deallocate(geo_dbdt)
       deallocate(geo_dbdt2)
       deallocate(geo_bp)
       deallocate(geo_bt)
       deallocate(geo_gsin)
       deallocate(geo_gcos1)
       deallocate(geo_gcos2)
       deallocate(geo_g_theta)
       deallocate(geo_grad_r)
       deallocate(geo_gq)
       deallocate(geo_captheta)
       deallocate(geo_nu)
       deallocate(geo_l_r)
       deallocate(geo_l_t)
       deallocate(geo_nsin)
       deallocate(geo_usin)
       deallocate(geo_ucos)
       deallocate(geo_bigr)
       deallocate(geo_bigr_r)
       deallocate(geo_bigr_t)
       deallocate(geo_theta_nc)
       deallocate(geo_theta_s)
       deallocate(geo_chi2)
    endif

  end subroutine geo_salloc

end module geo
