subroutine cgyro_init_rotation

  use timer_lib
  use cgyro_globals
  use cgyro_io
  use geo
  
  implicit none

  real, dimension(:), allocatable :: phi_rot, phi_rot_tderiv, phi_rot_rderiv
  real, dimension(:), allocatable :: pr_r
  real, dimension(:), allocatable :: bstar
  real :: phi_rot_avg
  integer, dimension(:), allocatable :: thcyc
  real, dimension(:), allocatable :: thcderiv
  real :: x, x0, fac, sum_zn, dsum_zn
  integer :: is, it, j, id, jt
  integer, parameter :: jmax = 200
  real, dimension(:), allocatable :: dens_deriv, jacob_r
  real, dimension(:,:), allocatable :: jacob
  real :: sum1, sum2, sum3, sum4, sum_j, sum_jr
  real :: gt_ave,gt0,err
  real, dimension(n_theta+1) :: xtmp,ytmp
  real, dimension(n_theta) :: ttmp
  real, parameter :: tol=1e-14
  
  if (rotation_model == 1) then
     ! O(M) terms only
     dens_rot(:,:)          = 1.0
     dens_ele_rot(:)        = 1.0
     dens_avg_rot(:)        = 1.0
     do is=1,n_species
        dlnndr_avg_rot(is)  = dlnndr(is)
     enddo
     lambda_rot(:,:)        = 0.0
     dlambda_rot(:,:)       = 0.0
     omega_rot_trap(:,:)    = 0.0
     omega_rot_u(:,:)       = 0.0
     omega_rot_drift(:,:)   = 0.0
     omega_rot_drift_r(:,:) = 0.0
     omega_rot_star(:,:)    = 0.0
     omega_rot_edrift(:)    = 0.0
     omega_rot_edrift_r(:)  = 0.0
     return
  endif
  
  allocate(phi_rot(n_theta))
  allocate(phi_rot_tderiv(n_theta))
  allocate(phi_rot_rderiv(n_theta))
  allocate(pr_r(n_theta))
  allocate(bstar(0:n_beta_star))
  
  ! solve the quasi-neutrality relation for poloidal part of phi

  ! partial component of n/n(theta0) 
  ! -- will add the phi_rot component after the QN solve
  do is=1, n_species
     do it=1,n_theta
        dens_rot(it,is) = exp(0.5 * (mach/vth(is))**2 &
             * (bigr(it)**2 - bigr_th0**2) / rmaj**2) 
     enddo
  enddo

  if (ae_flag == 1) then
     dens_ele_rot(:) = 1.0
  else
     dens_ele_rot(:) = exp(0.5 * (mach/vth(is_ele))**2 &
          * (bigr(:)**2 - bigr_th0**2) / rmaj**2)
  endif
     
  phi_rot_avg = 0.0
  x = 0.05          ! initial guess for phi_rot(1)
  do it=1,n_theta
     
     j=1
     
     do
        ! use Newton's method to solve the quasi-neutrality relation
        
        sum_zn  = 0.0
        dsum_zn = 0.0
        do is=1, n_species
           fac = z(is) * dens(is) * dens_rot(it,is) &
                * exp(-x * z(is)/temp(is))
           sum_zn  = sum_zn  + fac
           dsum_zn = dsum_zn - fac * z(is)/temp(is) 
        enddo
        
        if (ae_flag == 1) then
           fac = -dens_ele * dens_ele_rot(it) * exp(x/temp_ele)
           sum_zn  = sum_zn  + fac
           dsum_zn = dsum_zn + fac/temp_ele 
        endif

        x0 = x
        x  = x0 - sum_zn / dsum_zn

        if (abs(x-x0) < 1.0e-12) exit
        
        j = j + 1
        if(j > jmax) exit
        
     enddo

     if (j > jmax) then
        call cgyro_error('Rotation poloidal density computation failed to converge')
        return
     endif
     
     ! phi_rot = [phi(theta) - phi(theta0)]
     phi_rot(it) = x

     ! <phi_rot>
     phi_rot_avg = phi_rot_avg + w_theta(it) * phi_rot(it)

  enddo
  
  ! cyclic index (for theta-periodicity)
  allocate(thcyc(1-n_theta:2*n_theta))
  do it=1,n_theta
     thcyc(it-n_theta) = it
     thcyc(it) = it
     thcyc(it+n_theta) = it
  enddo

  ! 4th-order centered derivative
  allocate(thcderiv(-2:2))
  thcderiv(-2) =  1.0 / (12.0 * d_theta)
  thcderiv(-1) = -8.0 / (12.0 * d_theta)
  thcderiv(0)  =  0.0 / (12.0 * d_theta)
  thcderiv(1)  =  8.0 / (12.0 * d_theta)
  thcderiv(2)  = -1.0 / (12.0 * d_theta)

  ! d phi_rot / d theta
  phi_rot_tderiv(:) = 0.0
  do it=1,n_theta
     do id=-2,2
        jt = thcyc(it+id)
        phi_rot_tderiv(it) = phi_rot_tderiv(it) &
             + phi_rot(jt) * thcderiv(id) 
     enddo
  enddo
  deallocate(thcyc)
  deallocate(thcderiv)
  
  ! n(theta)/n(theta0)
  do is=1,n_species
     do it=1,n_theta
        dens_rot(it,is) = dens_rot(it,is) * exp(-z(is)/temp(is)*phi_rot(it))
     enddo
  enddo
  dens_ele_rot(:) = dens_ele_rot(:) * exp(1.0/temp_ele*phi_rot(:))

  ! Solve d/dr of QN to get d phi_rot / dr

  do it=1,n_theta
     
     ! sum_a z^2/T n(theta)
     sum_zn = 0.0
     do is=1,n_species
        sum_zn = sum_zn + z(is)*z(is)/temp(is)*dens(is)*dens_rot(it,is)
     enddo
     if (ae_flag == 1) then
        sum_zn = sum_zn + 1.0/temp_ele*dens_ele*dens_ele_rot(it)
     endif

     phi_rot_rderiv(it) = 0.0
     pr_r(it) = 0.0
     do is=1,n_species

        phi_rot_rderiv(it) = phi_rot_rderiv(it) &
             + z(is)*dens(is)*dens_rot(it,is)*(-dlnndr(is) &
             - dlntdr(is)*(z(is)/temp(is)*phi_rot(it) &
             - 0.5*(mach/rmaj/vth(is))**2 * (bigr(it)**2 - bigr_th0**2)) &
             + (mach/rmaj/vth(is))**2 * (bigr(it) * bigr_r(it) &
             - bigr_th0 * bigr_r_th0) + (mach/rmaj)/vth(is)**2 &
             * (-gamma_p/rmaj) * (bigr(it)**2 - bigr_th0**2))
        
        ! effective pressure gradient
        pr_r(it) = pr_r(it) &
             + temp(is)*dens(is)*dens_rot(it,is)*(-dlnndr(is) &
             - dlntdr(is)*(1.0 + z(is)/temp(is)*phi_rot(it) &
             - 0.5*(mach/rmaj/vth(is))**2 * (bigr(it)**2 - bigr_th0**2)) &
             - (mach/rmaj/vth(is))**2 * (bigr_th0 * bigr_r_th0) &
             + (mach/rmaj)/vth(is)**2 &
             * (-gamma_p/rmaj) * (bigr(it)**2 - bigr_th0**2))

     enddo
     if (ae_flag == 1) then
      
        phi_rot_rderiv(it) = phi_rot_rderiv(it) &
             - dens_ele*dens_ele_rot(it)*(-dlnndr_ele &
             + dlntdr_ele/temp_ele*phi_rot(it))
        
        pr_r(it) = pr_r(it) &
             + temp_ele*dens_ele*dens_ele_rot(it)*(-dlnndr_ele - dlntdr_ele)
             
     endif

     phi_rot_rderiv(it) = phi_rot_rderiv(it) / sum_zn
     
  enddo
  
  ! beta_star drift
  ! dp/dr -> dp/dr - sum_a n_a m_a omega^2 R dR/dr
  ! Compute projections of beta_star
  bstar(:) = 0.0
  do j=0,n_beta_star
     do it=1,n_theta
        bstar(j) = bstar(j) + cos(j*theta(it)) * pr_r(it) &
             * g_theta(it)/g_theta_geo(it)
     enddo
     if (j == 0) then
        bstar(j) = bstar(j)/n_theta
     else
        bstar(j) = bstar(j)*2.0/n_theta
     endif
  enddo
  ! beta_star(theta) = beta_star_0 + beta_star_1 * (1-cos(th))
  !                    + beta_star_1 * (1-cos(2*th))
  beta_star(0) = bstar(0)
  do j=1,n_beta_star
     beta_star(j) = -bstar(j)
     beta_star(0) = beta_star(0) + bstar(j)
  enddo
  beta_star(:) = beta_star(:) * beta_star_fac

  !do it=1,n_theta
  !   print *, theta(it), beta_star(0) &
  !        + beta_star(1)*(1-cos(theta(it))) &
  !        + beta_star(2)*(1-cos(2.0*theta(it)))
  !enddo
  
  ! Put rotation in equilibrium only if Mach^2 effects are included
  if (rotation_model == 2 .or. rotation_model == 3) then
     GEO_beta_star_in   = beta_star(0)
     GEO_beta_star_1_in = beta_star(1)
     GEO_beta_star_2_in = beta_star(2)
  endif

  ! Need to recompute interpolation since GEO_beta_star now set
  call geo_interp(n_theta,theta,.true.)
  
  do it=1,n_theta
 
     do is=1,n_species
     
        lambda_rot(it,is) = z(is)/temp(is) * (phi_rot(it) - phi_rot_avg) &
             - 0.5 * (mach * bigr(it) / rmaj / vth(is))**2

        ! bhat dot grad lambda
        dlambda_rot(it,is) = z(is)/temp(is) * phi_rot_tderiv(it) &
             / (q*rmaj*g_theta(it)) &
             - (mach / rmaj /vth(is))**2 * GEO_bigr(it) &
             * GEO_bigr_t(it) / (q*rmaj*GEO_g_theta(it))
           
        omega_rot_trap(it,is) = -0.5*sqrt(2.0)*vth(is) *dlambda_rot(it,is)

        omega_rot_u(it,is) = -vth(is)/sqrt(2.0)*dlambda_rot(it,is)
        
        omega_rot_star(it,is) = dlntdr(is) * (z(is)/temp(is)*phi_rot(it) &
             - 0.5*(mach/vth(is))**2 * (GEO_bigr(it)**2 - bigr_th0**2)/rmaj**2)  &
             + mach*gamma_p/vth(is)**2 * (GEO_bigr(it)**2 - bigr_th0**2)/rmaj**2 &
             + (mach/vth(is))**2 * bigr_th0/rmaj**2 * bigr_r_th0 

        omega_rot_drift(it,is) = -(mach/rmaj)**2 * GEO_bigr(it) * rho &
             * mass(is)/(z(is)*GEO_b(it)) * GEO_gq(it) &
             * (GEO_captheta(it)*GEO_usin(it)*GEO_bt(it)/GEO_b(it) + GEO_ucos(it)*GEO_b(it)/GEO_bt(it))

        omega_rot_drift_r(it,is) = -(mach/rmaj)**2 * GEO_bigr(it) * rho &
             * mass(is)/(z(is)*GEO_b(it)) * GEO_usin(it) * GEO_grad_r(it) &
             * GEO_bt(it) / GEO_b(it) 
        
     enddo

     omega_rot_edrift(it) =  rho*phi_rot_rderiv(it) &
          - rho * phi_rot_tderiv(it) / (q*rmaj*g_theta(it)) &
          * (GEO_bt(it)/GEO_bp(it) * GEO_captheta(it) / GEO_grad_r(it) &
          + GEO_l_r(it) * GEO_b(it)/GEO_bp(it))
     
     omega_rot_edrift_r(it) = -rho * GEO_bt(it)/GEO_bp(it)/GEO_b(it) * GEO_grad_r(it) &
          * phi_rot_tderiv(it) / (q*rmaj*g_theta(it))
     
  enddo  
  
  if (rotation_model == 3) then
     ! O(mach^2) terms only 
     ! no O(mach) terms
     mach_one_fac = 0.0
  endif

  if (rotation_model == 4) then
     ! O(mach^2) terms from "GKW CF TRAP" only
     ! no O(mach) terms and no "GKW CF DRIFT" term
     mach_one_fac = 0.0
     
     omega_rot_drift(:,:)   = 0.0
     omega_rot_drift_r(:,:) = 0.0

  endif

  if (rotation_model == 5) then
     ! O(mach^2) terms from "GKW CF DRIFT" only
     ! no O(mach) terms and no "GKW CF TRAP" term
     mach_one_fac = 0.0
     
     dens_rot(:,:)         = 1.0
     dens_ele_rot(:)       = 1.0
     phi_rot(:)            = 0.0
     phi_rot_tderiv(:)     = 0.0
     phi_rot_rderiv(:)     = 0.0
     lambda_rot(:,:)       = 0.0
     dlambda_rot(:,:)      = 0.0
     omega_rot_trap(:,:)   = 0.0
     omega_rot_u(:,:)      = 0.0
     omega_rot_star(:,:)   = 0.0
     omega_rot_edrift(:)   = 0.0
     omega_rot_edrift_r(:) = 0.0
     
  endif

  ! print out some diagnostics
  if (silent_flag == 0 .and. i_proc == 0) then
     open(unit=io,file=trim(path)//'out.cgyro.rotation',status='replace')
     do it=1,n_theta
        write(io,'(es16.9,1x)') phi_rot(it)
     enddo
     do is=1,n_species
        do it=1,n_theta
           write(io,'(es16.9,1x)') dens(is)*dens_rot(it,is)
        enddo
     enddo
     close(io)
  endif
  allocate(jacob(2,n_theta))
  allocate(jacob_r(n_theta))
  ! jacob(1) = jacob at rmin
  ! jacob(2) = jacob at rmin+0.001
  do id=2,1,-1
     if(id == 1) then
        GEO_rmin_in = rmin
     else
        GEO_rmin_in = rmin+0.001
     endif
     ttmp(1) = 0.0
     call geo_interp(1,ttmp,.true.)
     if (constant_stream_flag == 1) then
        do it=1,n_theta+1
           ytmp(it) = -pi+(it-1)*d_theta
           xtmp(it) = ytmp(it)
        enddo
        err = 1e4
        gt_ave = GEO_g_theta(1)
        do while (err > tol)
           do it=1,n_theta
              ttmp(it) = 0.5*(ytmp(it)+ytmp(it+1))
           enddo
           call geo_interp(n_theta,ttmp,.false.)
           do it=1,n_theta
              xtmp(it+1) = xtmp(it)+d_theta/GEO_g_theta(it)
           enddo
           gt0 = gt_ave
           gt_ave  = (2*pi)/(xtmp(n_theta+1)-xtmp(1))
           ytmp   = (xtmp-xtmp(1))*gt_ave-pi
           err = abs((gt_ave-gt0)/gt_ave)
        enddo
     endif
     call geo_interp(n_theta,theta,.false.)     
     do it=1,n_theta
        if (constant_stream_flag == 0) then
           jacob(id,it) = geo_g_theta(it)/geo_b(it)
        else
           jacob(id,it) = gt_ave/geo_b(it)
        endif
     enddo
  enddo
  ! Compute dJ/dr
  do it=1,n_theta
     jacob_r(it) = (jacob(2,it)-jacob(1,it))/0.001
  enddo
  sum_j = 0.0
  sum_jr = 0.0
  do it=1,n_theta
     sum_j  = sum_j  + jacob(1,it)   * d_theta  ! int j
     sum_jr = sum_jr + jacob_r(it) * d_theta  ! int jr
  enddo
  allocate(dens_deriv(n_theta))
  ! sum1 = <n_a>/n_a(theta0)
  do is=1,n_species
     sum1=0.0
     do it=1,n_theta
        sum1 = sum1 + dens_rot(it,is)*w_theta(it)
     enddo
     dens_avg_rot(is) = sum1
     ! A: dn/dr = n d (ln n(theta0))/dr + A
     do it=1,n_theta
        dens_deriv(it) = dens(is)*dens_rot(it,is) &
             * (- z(is)/temp(is)*phi_rot_rderiv(it) &
             - dlntdr(is)*(z(is)/temp(is)*phi_rot(it) &
             - 0.5*(mach/rmaj/vth(is))**2 * (bigr(it)**2 - bigr_th0**2)) &
             + (mach/rmaj/vth(is))**2 * (bigr(it) * bigr_r(it) &
             - bigr_th0 * bigr_r_th0) + (mach/rmaj)/vth(is)**2 &
             * (-gamma_p/rmaj) * (bigr(it)**2 - bigr_th0**2))
     enddo
     ! Compute -d (ln n(theta0))/dr such that 1/<n_a> d<n_a>/dr = const = dlnndr_a
     sum2 = 0.0
     sum3 = 0.0
     sum4 = 0.0
     do it=1,n_theta
        sum2 = sum2 + jacob_r(it) * dens(is)*dens_rot(it,is) * d_theta ! int jr n
        sum3 = sum3 + jacob(1,it) * dens_deriv(it) * d_theta           ! int j nr
        sum4 = sum4 + jacob(1,it) * dens(is)*dens_rot(it,is) * d_theta ! int j n
     enddo
     dlnndr_avg_rot(is) = (dlnndr(is)*dens(is)*sum1*sum_j+sum2+sum3)/sum4-sum_jr/sum_j
  enddo
  deallocate(dens_deriv)
  deallocate(jacob)
  deallocate(jacob_r)
  
  deallocate(phi_rot)
  deallocate(phi_rot_tderiv)
  deallocate(phi_rot_rderiv)
  deallocate(pr_r)
  deallocate(bstar)
  
end subroutine cgyro_init_rotation
