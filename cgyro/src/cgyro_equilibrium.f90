subroutine cgyro_equilibrium

  use cgyro_globals
  use GEO_interface

  implicit none

  integer :: it,ir,is,r
  real :: gtheta_ave,gtheta0,err
  real, dimension(n_theta+1) :: x,y

  real, parameter :: tol=1e-14

  ! Parameters needed for equilibrium
  ! geo_numeq_flag, geo_ny, and geo_yin already set 

  GEO_signb_in     = -btccw
  GEO_rmin_in      = rmin
  GEO_rmaj_in      = rmaj
  GEO_drmaj_in     = shift
  GEO_zmag_in      = zmag
  GEO_dzmag_in     = dzmag
  GEO_q_in         = q
  GEO_s_in         = s
  GEO_kappa_in     = kappa
  GEO_s_kappa_in   = s_kappa
  GEO_delta_in     = delta
  GEO_s_delta_in   = s_delta
  GEO_zeta_in      = zeta
  GEO_s_zeta_in    = s_zeta
  GEO_beta_star_in = beta_star(0)
  GEO_fourier_in(:,:) = geo_yin(:,:)

  call GEO_do()  

  !-----------------------------------------------------------------
  ! Generate theta-grid (equally-spaced or constant-wind-speed)
  !
  !
  d_theta = (2*pi/n_theta)
  do it=1,n_theta+1
     y(it) = -pi+(it-1)*d_theta
     x(it) = y(it)
  enddo

  if (constant_stream_flag == 1) then

     ! At the end of this process, theta will NOT be equally-spaced, 
     ! and the parallel motion is 
     !
     !       1       d         1          d
     ! ----------- ------ = -------- ---------- 
     ! GEO_g_theta dtheta   g_theta   dtheta_eq
     !
     ! However, theta_eq never appears explicitly EXCEPT for this derivative,
     ! so stencils are for the equally-spaced theta_eq grid.

     gtheta_ave = GEO_g_theta  ! at theta=0
     err = 1e4

     do while (err > tol)
        do it=1,n_theta
           call GEO_interp(0.5*(y(it)+y(it+1)))
           x(it+1) = x(it)+d_theta/GEO_g_theta
        enddo
        gtheta0 = gtheta_ave
        gtheta_ave  = (2*pi)/(x(n_theta+1)-x(1))
        y   = (x-x(1))*gtheta_ave-pi
        err = abs((gtheta_ave-gtheta0)/gtheta_ave)
     enddo

     ! Now, d_theta is really a constant (dtheta_eq).  This is only ever 
     ! used for finite-difference stencils.

  endif
  
  ! This theta grid is:
  !
  ! 1. NOT EQUALLY SPACED if constant_stream_flag == 1
  ! 2. Actually the real theta.
  !
  theta(:) = y(1:n_theta)

  !--------------------------------------------------------
  ! Manage subset of theta-values for plotting output
  !
  if (zf_test_mode == 0) then
     ! Location of theta=0
     it0 = n_theta/2+1
  else
     it0 = n_theta/3+1
  endif

  r = n_theta/theta_plot

  itp(:) = 0
  if (theta_plot == 1) then
     itp(it0) = 1
  else
     do it=1,n_theta
        if (modulo(it,r) == 1) itp(it) = it/r+1
     enddo
  endif
  !-----------------------------------------------------------------

  do ir=1,n_radial/box_size
     do it=1,n_theta
        if (ipccw*btccw > 0) then
           thetab(ir,it) = theta(it)+2*pi*(ir-1-n_radial/2/box_size)
        else
           thetab(ir,it) = -theta(it)+2*pi*(-ir+n_radial/2/box_size)
        endif
     enddo
  enddo

  call GEO_interp(0.0)
  bigr_th0   = GEO_bigr
  bigr_r_th0 = GEO_bigr_r
  
  do it=1,n_theta

     call GEO_interp(theta(it))     

     bigr(it)   = GEO_bigr
     bigr_r(it) = GEO_bigr_r
     bmag(it)   = GEO_b
     btor(it)   = GEO_bt
     bpol(it)   = GEO_bp

     ! Define modified G_theta
     if (constant_stream_flag == 0) then
        ! Theta-dependent
        g_theta(it) = GEO_g_theta
     else
        ! Constant by construction
        g_theta(it) = gtheta_ave
     endif
     g_theta_geo(it) = GEO_g_theta

     ! flux-surface average weights
     w_theta(it) = g_theta(it)/GEO_b

  enddo

  w_theta(:) = w_theta(:)/sum(w_theta) 

  mach_one_fac = 1.0
  
  ! 1. Compute rotation (M^2) terms IF required
  !
  ! NOTE: this changes beta_star and thus GEO (which affects gcos2 
  ! and captheta)
  call cgyro_init_rotation


  ! 2. Compute terms required for O(M) rotation, and no rotation.
  ! 
  do it=1,n_theta

     call GEO_interp(theta(it))

     do is=1,n_species
  
        omega_stream(it,is) = sqrt(2.0)*vth(is)/(q*rmaj*g_theta(it))
        
        omega_trap(it,is) = -0.5*sqrt(2.0)*vth(is) &
             *(GEO_dbdt/GEO_b)/(q*rmaj*GEO_g_theta) 
        
        omega_rdrift(it,is) = -rho*vth(is)**2*mass(is)/(Z(is)*GEO_b) &
             *GEO_grad_r/rmaj*GEO_gsin 

        omega_adrift(it,is) = -rho*vth(is)**2*mass(is)/(Z(is)*GEO_b) &
             *GEO_gq/rmaj*(GEO_gcos1+GEO_gcos2+GEO_captheta*GEO_gsin)

        omega_aprdrift(it,is) = 2.0*rho*vth(is)**2 &
             *mass(is)/(Z(is)*GEO_b)*GEO_gq/rmaj*GEO_gcos2
        
        omega_cdrift(it,is) = -2.0*sqrt(2.0)*rho*vth(is) &
             * mass(is)/(Z(is)*GEO_b)*GEO_gq/rmaj &
             *(GEO_ucos+GEO_captheta*GEO_usin)*mach*mach_one_fac
        
        omega_cdrift_r(it,is) = -2.0*sqrt(2.0)*rho*vth(is) &
             * mass(is)/(Z(is)*GEO_b)*GEO_grad_r/rmaj*GEO_usin &
             * mach*mach_one_fac
 
     enddo
     
     omega_gammap(it) = GEO_bt/GEO_b*GEO_bigr/rmaj*gamma_p*mach_one_fac

     do ir=1,n_radial
        k_perp(ic_c(ir,it)) = sqrt((2.0*pi*px(ir)*GEO_grad_r/length &
             + k_theta*GEO_gq*GEO_captheta)**2 &
             + (k_theta*GEO_gq)**2)
        k_x(ic_c(ir,it)) = 2.0*pi*px(ir)*GEO_grad_r/length &
             + k_theta*GEO_gq*GEO_captheta
     enddo
     
  enddo
!$acc enter data copyin(xi,vel,omega_stream)
  
end subroutine cgyro_equilibrium


