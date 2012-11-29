!------------------------------------------------
! GEO_do.f90
!
! PURPOSE:
!  Calculation of geometry coefficients for 
!  Miller local equilibrium model.
!
!  See usage description in GEO_interface.f90
!------------------------------------------------
 
subroutine GEO_do()

  use GEO_interface,&
       b => GEOV_b, &
       dbdt2 => GEOV_dbdt2, &
       bp => GEOV_bp, &
       bt => GEOV_bt, &
       gcos1 => GEOV_gcos1, &
       gcos2 => GEOV_gcos2, &
       g_theta => GEOV_g_theta, &
       grad_r => GEOV_grad_r, &
       gq => GEOV_gq, &
       dbdt => GEOV_dbdt, &
       gsin => GEOV_gsin, &
       captheta => GEOV_captheta, &
       nu => GEOV_nu, &
       thetav => GEOV_theta, &
       l_t => GEOV_l_t, &
       nsin => GEOV_nsin, &
       usin => GEOV_usin, &
       ucos => GEOV_ucos, &
       bigr => GEOV_bigr, &
       bigr_r => GEOV_bigr_r, &
       bigr_t => GEOV_bigr_t

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
  real :: theta
  real :: d_theta
  real :: a
  real :: a_t
  real :: a_tt
  real :: x
  real :: bigr_tt
  real :: bigz_tt
  real :: bigz_r
  real :: g_tt
  real :: jac_r
  real :: f
  real :: f_prime
  real :: c 
  real :: pi_2
  real :: denom
  !
  real :: b1
  real :: b2
  real :: b3
  real :: b4
  real :: b5
  !
  real, dimension(:), allocatable :: bigz
  real, dimension(:), allocatable :: bigz_t
  real, dimension(:), allocatable :: r_c
  real, dimension(:), allocatable :: bigz_l
  real, dimension(:), allocatable :: dbdl
  real, dimension(:,:), allocatable :: e
  real, dimension(:,:), allocatable :: ei
  real, dimension(:), allocatable :: loop
  real, dimension(:), allocatable :: a_R,b_R,a_Z,b_Z
  real, dimension(:), allocatable :: a_Rp,b_Rp,a_Zp,b_Zp
  !
  !-----------------------------------------------------------
    
  !-----------------------------------------------------------
  ! Check for missing value
  !
  if (abs(GEO_signb_in) < 1e-10) then
     print *,'bad value for GEO_signb_in'
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
  a_R(:)  = GEO_fourier_in(1,:)
  b_R(:)  = GEO_fourier_in(2,:)
  a_Z(:)  = GEO_fourier_in(3,:)
  b_Z(:)  = GEO_fourier_in(4,:)
  a_Rp(:) = GEO_fourier_in(5,:)
  b_Rp(:) = GEO_fourier_in(6,:)
  a_Zp(:) = GEO_fourier_in(7,:)
  b_Zp(:) = GEO_fourier_in(8,:)
  !-----------------------------------------------------------
 
  !-----------------------------------------------------------
  ! Allocate internal variables
  !
  n_theta = GEO_ntheta_in
  !
  allocate(bigz(n_theta))
  allocate(bigz_t(n_theta))
  allocate(r_c(n_theta))
  allocate(dbdl(n_theta))
  allocate(bigz_l(n_theta))
  allocate(e(n_theta,4))
  allocate(ic(2-n_theta:2*n_theta-2))
  allocate(ei(n_theta,4))
  allocate(loop(4))
  !-----------------------------------------------------------

  pi_2 = 8.0*atan(1.0)

  do i=1,n_theta-1
     ic(i) = i
     ic(i-(n_theta-1)) = i
     ic(i+(n_theta-1)) = i
  enddo

  d_theta = pi_2/(n_theta-1)

  !------------------------------------------------------------------
  ! Compute fundamental geometric quantities (basic derivatives
  ! like dl/dt) and metric elements (g_tt).
  !
  do i=1,n_theta

     theta = -0.5*pi_2+(i-1)*d_theta

     thetav(i) = theta

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
        bigr(i) = GEO_rmaj_in+GEO_rmin_in*cos(a)
        bigr_r(i) = GEO_drmaj_in+cos(a)-GEO_s_delta_in/cos(x)*sin(theta)*sin(a)
        bigr_t(i) = -GEO_rmin_in*a_t*sin(a)
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
        bigz(i) = GEO_zmag_in+GEO_kappa_in*GEO_rmin_in*sin(a)
        bigz_r  = GEO_dzmag_in+GEO_kappa_in*(1.0+GEO_s_kappa_in)*sin(a)+&
             GEO_kappa_in*GEO_s_zeta_in*cos(a)*sin(2*theta)
        bigz_t(i) = GEO_kappa_in*GEO_rmin_in*cos(a)*a_t
        bigz_tt = -GEO_kappa_in*GEO_rmin_in*sin(a)*a_t**2+&
             GEO_kappa_in*GEO_rmin_in*cos(a)*a_tt

     else

        !-----------------------------------------
        ! Fourier-expansion (completely general)
        !-----------------------------------------

        bigr(i)   = 0.5*a_R(0)
        bigr_r(i) = 0.5*a_Rp(0)
        bigr_t(i) = 0.0
        bigr_tt   = 0.0
        do n=1,ny
           bigr(i) = bigr(i)+a_R(n)*cos(n*theta)+b_R(n)*sin(n*theta)        
           bigr_r(i)  = bigr_r(i)+a_Rp(n)*cos(n*theta)+b_Rp(n)*sin(n*theta)        
           bigr_t(i)  = bigr_t(i)-n*a_R(n)*sin(n*theta)+n*b_R(n)*cos(n*theta) 
           bigr_tt = bigr_tt-n*n*(a_R(n)*cos(n*theta)+b_R(n)*sin(n*theta)) 
        enddo

        bigz(i)  = 0.5*a_Z(0)
        bigz_r   = 0.5*a_Zp(0)
        bigz_t(i) = 0.0
        bigz_tt   = 0.0
        do n=1,ny
           bigz(i) = bigz(i)+a_Z(n)*cos(n*theta)+b_Z(n)*sin(n*theta)        
           bigz_r  = bigz_r+a_Zp(n)*cos(n*theta)+b_Zp(n)*sin(n*theta)        
           bigz_t(i)  = bigz_t(i)-n*a_Z(n)*sin(n*theta)+n*b_Z(n)*cos(n*theta) 
           bigz_tt = bigz_tt-n*n*(a_Z(n)*cos(n*theta)+b_Z(n)*sin(n*theta)) 
        enddo

     endif

     g_tt = bigr_t(i)**2+bigz_t(i)**2

     jac_r = bigr(i)*(bigr_r(i)*bigz_t(i)-bigr_t(i)*bigz_r)

     grad_r(i) = bigr(i)*sqrt(g_tt)/jac_r

     l_t(i) = sqrt(g_tt)

     ! 1/(du/dl)
     r_c(i) = l_t(i)**3/(bigr_t(i)*bigz_tt-bigz_t(i)*bigr_tt)

     ! cos(u)
     bigz_l(i) = bigz_t(i)/l_t(i)

     nsin(i) = (bigr_r(i)*bigr_t(i)+bigz_r*bigz_t(i))/l_t(i)

  enddo
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Loop integral (1 to n_theta-1) to compute f
  !
  c = 0.0
  do i=1,n_theta-1
     c = c+l_t(i)/(bigr(i)*grad_r(i))
  enddo
  f = GEO_rmin_in/(c*d_theta/pi_2)
  !
  ! Loop integral to compute V'
  !
  c = 0.0
  do i=1,n_theta-1
     c = c+l_t(i)*bigr(i)/grad_r(i)
  enddo
  GEO_volume_prime = pi_2*c*d_theta
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! bt (toroidal field, Bt) 
  ! bp (poloidal field, Bp) 
  ! b  (total field, B)
  !
  do i=1,n_theta
     bt(i) = f/bigr(i)
     bp(i) = (GEO_rmin_in/GEO_q_in)*grad_r(i)/bigr(i)
     b(i)  = GEO_signb_in*sqrt(bt(i)**2+bp(i)**2)
  enddo
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! dbdl  (db/dl)
  ! dbdt  (db/d(theta))
  ! dbdt2 (d^2b/dt^2)
  ! gsin  (generalized sine)
  ! gcos1 (generalized cosine)
  ! gcos2 
  ! g_theta    (G_theta)
  ! gq    (G_q)
  !
  !
  ! Use 5-point stencils for derivatives.  We get poor 
  ! accuracy for db/dt without this.
  ! 
  do i=1,n_theta

     b5 = b(ic(i+2))
     b4 = b(ic(i+1))
     b3 = b(ic(i))
     b2 = b(ic(i-1))
     b1 = b(ic(i-2))

     dbdt(i)  = (-b5+8.0*b4-8.0*b2+b1)/(12.0*d_theta)
     dbdl(i)  = dbdt(i)/l_t(i)
     dbdt2(i) = (-b5+16.0*b4-30.0*b3+16.0*b2-b1)/(12.0*d_theta**2)
     gsin(i)  = bt(i)*GEO_rmaj_in*dbdl(i)/b(i)**2
     gcos1(i) = (bt(i)**2/bigr(i)*bigz_l(i)+bp(i)**2/r_c(i))*GEO_rmaj_in/b(i)**2
     gcos2(i) = 0.5*(GEO_rmaj_in/b(i)**2)*grad_r(i)*(-GEO_beta_star_in)
     g_theta(i)    = bigr(i)*b(i)*l_t(i)/(GEO_rmin_in*GEO_rmaj_in*grad_r(i))
     gq(i)    = GEO_rmin_in*b(i)/(GEO_q_in*bigr(i)*bp(i))

     usin(i)  = -bigr_t(i)/l_t(i)
     ucos(i)  = (bt(i)/b(i))*bigz_l(i)

  enddo
  !------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Compute integrands for E1,E2,E3 and E4=nu
  !
  do i=1,n_theta
     c = d_theta*l_t(i)/(bigr(i)*grad_r(i))
     ei(i,1) = c*2.0*bt(i)/bp(i)*(GEO_rmin_in/r_c(i)-GEO_rmin_in*bigz_l(i)/bigr(i))
     ei(i,2) = c*b(i)**2/bp(i)**2
     ei(i,3) = c*grad_r(i)*0.5/bp(i)**2*(bt(i)/bp(i))     
     ei(i,4) = -c*grad_r(i)*(bt(i)/bp(i))     
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
  f_prime = (pi_2*GEO_q_in*GEO_s_in/GEO_rmin_in-loop(1)/GEO_rmin_in-&
       loop(3)*(-GEO_beta_star_in))/loop(2)

  do i=1,n_theta
     nu(i)    = e(i,4)
     captheta(i) = bp(i)/b(i)*grad_r(i)*bigr(i)* & 
          (e(i,1)/GEO_rmin_in+e(i,2)*f_prime+e(i,3)*(-GEO_beta_star_in))
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

  denom = sum(g_theta(1:n_theta-1)/b(1:n_theta-1))

  GEO_fluxsurfave_grad_r = sum(grad_r(1:n_theta-1)*g_theta(1:n_theta-1)/ &
       b(1:n_theta-1))/denom

  GEO_fluxsurfave_grad_r2 = sum(grad_r(1:n_theta-1)**2*g_theta(1:n_theta-1)/ &
       b(1:n_theta-1))/denom

  ! theta(i) = 0 for i = n_theta/2+1

  GEO_grad_r0   = grad_r(n_theta/2+1)
  GEO_ffprime   = f*f_prime
  GEO_beta_star = GEO_beta_star_in
  GEO_f         = f
  !
  ! pre-factor of 0.5 comes from triangular element in phi-direction:
  ! dV = (0.5*R*dphi)*(R*dZ) 
  !
  GEO_volume    = 0.5*pi_2*sum(bigz_t(1:n_theta-1)*bigr(1:n_theta-1)**2)*d_theta
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  ! Straight field-line angle
  !
  GEOV_theta_nc(1) = thetav(1)
  do i=2,n_theta
     GEOV_theta_nc(i) = GEOV_theta_nc(i-1)+0.5*(g_theta(i)+g_theta(i-1))*d_theta
  enddo
  GEOV_theta_nc(:) = -0.5*pi_2+pi_2*(0.5*pi_2+GEOV_theta_nc(:))/(0.5*pi_2+GEOV_theta_nc(n_theta))
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  ! Deallocate internal variables
  !
  deallocate(bigz)
  deallocate(bigz_t)
  deallocate(r_c)
  deallocate(dbdl)
  deallocate(bigz_l)
  deallocate(e)
  deallocate(ic)
  deallocate(ei)
  deallocate(loop)
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

end subroutine GEO_do

