subroutine cgyro_equilibrium

  use cgyro_globals
  use geo
  use cgyro_io

  implicit none

  integer :: m
  integer :: it,ir,is,itor
  real :: gtheta_ave,gtheta0,err
  real, dimension(n_theta+1) :: x,y
  real, dimension(n_theta) :: ttmp
  real, parameter :: tol=1e-14

  ! Compute equilibrium quantities (even in test mode)
  geo_model_in    = geo_numeq_flag
  geo_ntheta_in   = geo_ntheta
  geo_nfourier_in = geo_ny

  if (zf_test_mode == 0) then
     ! Location of theta=0
     it0 = n_theta/2+1
  else
     it0 = n_theta/3+1
  endif

  ! Parameters needed for equilibrium
  ! geo_numeq_flag, geo_ny, and geo_yin already set 

  geo_signb_in     = -btccw
  geo_rmin_in      = rmin
  geo_rmaj_in      = rmaj
  geo_drmaj_in     = shift
  geo_zmag_in      = zmag
  geo_dzmag_in     = dzmag
  geo_q_in         = q
  geo_s_in         = s
  geo_kappa_in     = kappa
  geo_s_kappa_in   = s_kappa
  geo_delta_in     = delta
  geo_s_delta_in   = s_delta
  geo_zeta_in      = zeta
  geo_s_zeta_in    = s_zeta
  geo_beta_star_in = beta_star(0)

  ! Extended Miller
  geo_shape_cos0_in = shape_cos(0)
  geo_shape_cos1_in = shape_cos(1)
  geo_shape_cos2_in = shape_cos(2)
  geo_shape_cos3_in = shape_cos(3)
  geo_shape_cos4_in = shape_cos(4)
  geo_shape_cos5_in = shape_cos(5)
  geo_shape_cos6_in = shape_cos(6)
  geo_shape_sin3_in = shape_sin(3)
  geo_shape_sin4_in = shape_sin(4)
  geo_shape_sin5_in = shape_sin(5)
  geo_shape_sin6_in = shape_sin(6)
  geo_shape_s_cos0_in = shape_s_cos(0)
  geo_shape_s_cos1_in = shape_s_cos(1)
  geo_shape_s_cos2_in = shape_s_cos(2)
  geo_shape_s_cos3_in = shape_s_cos(3)
  geo_shape_s_cos4_in = shape_s_cos(4)
  geo_shape_s_cos5_in = shape_s_cos(5)
  geo_shape_s_cos6_in = shape_s_cos(6)
  geo_shape_s_sin3_in = shape_s_sin(3)
  geo_shape_s_sin4_in = shape_s_sin(4)
  geo_shape_s_sin5_in = shape_s_sin(5)
  geo_shape_s_sin6_in = shape_s_sin(6)

  geo_fourier_in(:,0:geo_nfourier_in) = geo_yin(:,:)

  ! Get initial geo solution, then set R,dR/dr at theta=0 
  ttmp(1) = 0.0
  call geo_interp(1,ttmp,.true.)
  bigr_th0   = geo_bigr(1)
  bigr_r_th0 = geo_bigr_r(1)

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
     ! geo_g_theta dtheta   g_theta   dtheta_eq
     !
     ! However, theta_eq never appears explicitly EXCEPT for this derivative,
     ! so stencils are for the equally-spaced theta_eq grid.

     err = 1e4
     gtheta_ave = geo_g_theta(1)

     do while (err > tol)
        do it=1,n_theta
           ttmp(it) = 0.5*(y(it)+y(it+1))
        enddo
        call geo_interp(n_theta,ttmp,.false.)
        do it=1,n_theta
           x(it+1) = x(it)+d_theta/geo_g_theta(it)
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
  m = n_theta/theta_plot

  itp(:) = 0
  if (theta_plot == 1) then
     itp(it0) = 1
  else
     is = 1
     do it=1,n_theta,m  
        itp(it) = is
        is = is+1
     enddo
  endif
  !-----------------------------------------------------------------

  !---------------------------------------------------------------
  ! Construct balloon-mode extended angle
  ! (see also extended_ang in cgyro_write_timedata)
  !
  do ir=1,n_radial/box_size
     if (sign_qs > 0) then
        thetab(:,ir) = theta(:)+2*pi*(px(ir)+px0)
     else
        ! Reverse output direction (see extended_ang)
        thetab(:,n_radial/box_size-ir+1) = theta(:)-2*pi*(px(ir)+px0)
     endif
  enddo
  !-----------------------------------------------------------------

  call geo_interp(n_theta,theta,.false.)     

  do it=1,n_theta

     bigr(it)   = geo_bigr(it)
     bigr_r(it) = geo_bigr_r(it)
     bmag(it)   = geo_b(it)
     btor(it)   = geo_bt(it)
     bpol(it)   = geo_bp(it)

     ! Define modified G_theta
     if (constant_stream_flag == 0) then
        ! Theta-dependent
        g_theta(it) = geo_g_theta(it)
     else
        ! Constant by construction
        g_theta(it) = gtheta_ave
     endif
     g_theta_geo(it) = geo_g_theta(it)

     ! flux-surface average weights
     w_theta(it) = g_theta(it)/geo_b(it)

  enddo

  w_theta(:) = w_theta(:)/sum(w_theta) 

  ! GS2 normalizing magnetic field (B_GS2/B_unit)
  b_gs2 = geo_f/rmaj

  ! Geometry coefficient from plotting only
  captheta = geo_captheta
  
  mach_one_fac = 1.0

  ! 1. Compute rotation (M^2) terms IF required
  !
  ! NOTE: this changes beta_star and thus GEO (which affects gcos2 
  ! and captheta)
  call cgyro_init_rotation


  ! 2. Compute terms required for O(M) rotation, and no rotation.
  !

  ! EAB LE3 test
  open(unit=1,file='out.cgyro.geogk',status='replace')
  write (1,'(i2)') n_theta
  write (1,'(e16.8)') theta(:)
  write (1,'(e16.8)') geo_b(:)
  write (1,'(e16.8)') 1.0/(q*rmaj*geo_g_theta(:))
  write (1,'(e16.8)') (geo_dbdt(:)/geo_b(:))/(q*rmaj*geo_g_theta(:))
  write (1,'(e16.8)') -1.0/geo_b(:)*geo_grad_r(:)/rmaj*geo_gsin(:) 
  write (1,'(e16.8)') -1.0/geo_b(:)*geo_gq(:)/rmaj&
       *(geo_gcos1(:)+geo_gcos2(:)+geo_captheta(:)*geo_gsin(:))
  write (1,'(e16.8)') 1.0/geo_b(:)*geo_gq(:)/rmaj*geo_gcos2(:)
  write (1,'(e16.8)') geo_grad_r(:)**2
  write (1,'(e16.8)') geo_gq(:)**2 * (1.0+geo_captheta(:)**2)
  write (1,'(e16.8)') 2.0*geo_grad_r(:)*geo_gq(:)*geo_captheta(:)
  close(1)

  do it=1,n_theta

     do is=1,n_species

        omega_stream(it,is,:) = sqrt(2.0)*vth(is)/(q*rmaj*g_theta(it))

        omega_trap(it,is,:) = -0.5*sqrt(2.0)*vth(is) &
             *(geo_dbdt(it)/geo_b(it))/(q*rmaj*geo_g_theta(it)) 

        omega_rdrift(it,is) = -rho*vth(is)**2*mass(is)/(Z(is)*geo_b(it)) &
             *geo_grad_r(it)/rmaj*geo_gsin(it) 

        omega_adrift(it,is) = -rho*vth(is)**2*mass(is)/(Z(is)*geo_b(it)) &
             *geo_gq(it)/rmaj*(geo_gcos1(it)+geo_gcos2(it)+geo_captheta(it)*geo_gsin(it))

        omega_aprdrift(it,is) = 2.0*rho*vth(is)**2 &
             *mass(is)/(Z(is)*geo_b(it))*geo_gq(it)/rmaj*geo_gcos2(it)

        omega_cdrift(it,is) = -2.0*sqrt(2.0)*rho*vth(is) &
             * mass(is)/(Z(is)*geo_b(it))*geo_gq(it)/rmaj &
             *(geo_ucos(it)+geo_captheta(it)*geo_usin(it))*mach*mach_one_fac

        omega_cdrift_r(it,is) = -2.0*sqrt(2.0)*rho*vth(is) &
             * mass(is)/(Z(is)*geo_b(it))*geo_grad_r(it)/rmaj*geo_usin(it) &
             * mach*mach_one_fac

     enddo

     omega_gammap(it) = geo_bt(it)/geo_b(it)*geo_bigr(it)/rmaj*gamma_p*mach_one_fac

     do itor=nt1,nt2
       do ir=1,n_radial
        k_x(ic_c(ir,it),itor) = &
               2.0*pi*(px(ir)+px0)*geo_grad_r(it)/length &
             + k_theta_base*itor*geo_gq(it)*geo_captheta(it)
        k_perp(ic_c(ir,it),itor) = &
               sqrt(k_x(ic_c(ir,it),itor)**2+(k_theta_base*itor*geo_gq(it))**2)
       enddo
     enddo

  enddo

  do itor=nt1,nt2
   select case (stream_term)
   case (1)
     if (itor == 0) then
        omega_stream(:,1,itor) = stream_factor*omega_stream(:,1,itor)
        omega_trap(:,1,itor) = stream_factor*omega_trap(:,1,itor)
     endif
   case (2)
     if (itor > 0) then
        omega_stream(:,1,itor) = stream_factor*omega_stream(:,1,itor)
        omega_trap(:,1,itor) = stream_factor*omega_trap(:,1,itor)
     endif
   case (12)
     omega_stream(:,1,itor) = stream_factor*omega_stream(:,1,itor)
     omega_trap(:,1,itor) = stream_factor*omega_trap(:,1,itor)
   case (3)
     if (itor == 0 .and. n_species > 1) then
        omega_stream(:,2,itor) = stream_factor*omega_stream(:,2,itor)
        omega_trap(:,2,itor) = stream_factor*omega_trap(:,2,itor)
     endif
   case (4)
     if (itor > 0 .and. n_species > 1) then
        omega_stream(:,2,itor) = stream_factor*omega_stream(:,2,itor)
        omega_trap(:,2,itor) = stream_factor*omega_trap(:,2,itor)
     endif
   case (34)
     omega_stream(:,2,itor) = stream_factor*omega_stream(:,2,itor)
     omega_trap(:,2,itor) = stream_factor*omega_trap(:,2,itor)
   end select
  enddo

#if defined(OMPGPU)
!$omp target enter data map(to:xi,omega_stream,w_theta,dens2_rot)
#elif defined(_OPENACC)
!$acc enter data copyin(xi,omega_stream,w_theta,dens2_rot)
#endif

#if defined(OMPGPU) || defined(_OPENACC)
  if (explicit_trap_flag == 1) then
     call cgyro_error("explicit_trap_flag=1 not supported in GPU code.")
     return
  endif
#endif

end subroutine cgyro_equilibrium


