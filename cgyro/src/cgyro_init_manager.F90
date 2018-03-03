!-----------------------------------------------------------------
! cgyro_init_manager.f90
!
! PURPOSE:
!  Manage initialization of arrays and other setup issues.
!  In particular,
!
!  1. collisionless streaming arrays (timed setup)
!  2. collision arrays (timed)
!  3. write initial data
!  4. Construct initial distributions
!-----------------------------------------------------------------

subroutine cgyro_init_manager

  use mpi
  use timer_lib
  use cgyro_globals
  use GEO_interface

#ifdef _OPENACC
  use cgyro_io
  use precision_m, only : singlePrecision
  use cufft_m, only : cufftPlanMany, &
       CUFFT_C2R,CUFFT_Z2D,CUFFT_R2C,CUFFT_D2Z
#endif

  implicit none

  include 'fftw3.f03'

#ifdef _OPENACC
  integer :: howmany,istatus
  integer, parameter :: irank = 2
  integer, dimension(irank) :: ndim,inembed,onembed
  integer :: idist,odist,istride,ostride
#endif

  if (hiprec_flag == 1) then
     fmtstr  = '(es16.9)'
     fmtstr_len = 17
     fmtstrn = '(10(es16.9,1x))'
  endif

  !------------------------------------------------------
  ! Initialize startup timers 
  !  NOTE: "Runtime" timers are initialized in cgyro_write_timedata,
  !  and input timer is initialized by read_input.
  !------------------------------------------------------
  call timer_lib_init('str_init')
  call timer_lib_init('coll_init')
  call timer_lib_init('io_init')

  !----------------------------------------------------
  ! Initialize GLOBAL arrays
  !----------------------------------------------------

  allocate(energy(n_energy))
  allocate(vel(n_energy))
  allocate(w_e(n_energy))
  allocate(e_deriv1_mat(n_energy,n_energy))

  ! Construct energy nodes and weights
  call pseudo_maxwell_new(n_energy,&
       e_max,&
       energy,&
       w_e,&
       e_deriv1_mat,&
       trim(path)//'out.cgyro.egrid')

  vel(:) = sqrt(energy(:))

  allocate(xi(n_xi))
  allocate(w_xi(n_xi))
  allocate(xi_lor_mat(n_xi,n_xi))
  allocate(xi_deriv_mat(n_xi,n_xi))
  ! Construct xi (pitch-angle) nodes and weights
  call pseudo_legendre(n_xi,xi,w_xi,xi_deriv_mat,xi_lor_mat)
  w_xi = 0.5*w_xi

  allocate(theta(n_theta))
  allocate(thetab(n_radial/box_size,n_theta))
  allocate(w_theta(n_theta))
  allocate(g_theta(n_theta))
  allocate(g_theta_geo(n_theta))
  allocate(bmag(n_theta))
  allocate(btor(n_theta))
  allocate(bpol(n_theta))
  allocate(k_perp(nc))
  allocate(k_x(nc))
  allocate(bigr(n_theta))
  allocate(bigr_r(n_theta))
  allocate(itp(n_theta))
  allocate(omega_stream(n_theta,n_species))
  allocate(omega_trap(n_theta,n_species))
  allocate(omega_rdrift(n_theta,n_species))
  allocate(omega_adrift(n_theta,n_species))
  allocate(omega_aprdrift(n_theta,n_species))
  allocate(omega_cdrift(n_theta,n_species))
  allocate(omega_cdrift_r(n_theta,n_species))
  allocate(omega_gammap(n_theta))

  allocate(lambda_rot(n_theta,n_species))
  allocate(dlambda_rot(n_theta,n_species))
  allocate(dens_rot(n_theta,n_species))
  allocate(dens_ele_rot(n_theta))
  allocate(omega_rot_trap(n_theta,n_species))
  allocate(omega_rot_u(n_theta,n_species))
  allocate(omega_rot_drift(n_theta,n_species))
  allocate(omega_rot_drift_r(n_theta,n_species))
  allocate(omega_rot_edrift(n_theta))
  allocate(omega_rot_edrift_r(n_theta))
  allocate(omega_rot_star(n_theta,n_species))

  if (test_flag == 0) then

     !----------------------------------------------------
     ! Initialize DISTRIBUTED arrays
     !----------------------------------------------------

     call timer_lib_in('str_init')

     ! Global (undistributed) arrays
     allocate(fcoef(n_field,nc))
     if (n_field < 3) then
        allocate(gcoef(n_field,nc))
     else
        allocate(gcoef(5,nc))
     endif
     allocate(field(n_field,nc))
     allocate(field_loc(n_field,nc))
     allocate(field_old(n_field,nc))
     allocate(field_old2(n_field,nc))
     allocate(field_old3(n_field,nc))
     allocate(    moment(n_radial,theta_plot,n_species,2))
     allocate(moment_loc(n_radial,theta_plot,n_species,2))
     allocate(     flux(n_radial,n_species))
     allocate( flux_loc(n_radial,n_species))
     allocate(    fflux(n_species,3,n_field))
     allocate(fflux_loc(n_species,3,n_field))
     allocate(    gflux(0:n_global,n_species,3))
     allocate(gflux_loc(0:n_global,n_species,3))
     allocate(recv_status(MPI_STATUS_SIZE))

     allocate(icd_c(-nup_theta:nup_theta, nc))
     allocate(dtheta(-nup_theta:nup_theta, nc))
     allocate(dtheta_up(-nup_theta:nup_theta, nc))

     ! Velocity-distributed arrays
     allocate(rhs(nc,nv_loc,4))
     allocate(h_x(nc,nv_loc))
     allocate(g_x(nc,nv_loc))
     allocate(psi(nc,nv_loc))
     allocate(chi(nc,nv_loc))
     allocate(h0_x(nc,nv_loc))
     allocate(cap_h_c(nc,nv_loc))
     allocate(cap_h_ct(nv_loc,nc))
     allocate(omega_cap_h(nc,nv_loc))
     allocate(omega_h(nc,nv_loc))
     allocate(omega_s(n_field,nc,nv_loc))
     allocate(omega_ss(n_field,nc,nv_loc))
     allocate(jvec_c(n_field,nc,nv_loc))
     allocate(jvec_v(n_field,nc_loc,nv))
     allocate(dvjvec_c(n_field,nc,nv_loc))
     allocate(dvjvec_v(n_field,nc_loc,nv))
     allocate(jxvec_c(n_field,nc,nv_loc))
     allocate(upfac1(nc,nv_loc))
     allocate(upfac2(nc,nv_loc))
     ! Real-space distributed arrays
     allocate(cap_h_v(nc_loc,nv))
     allocate(cap_h_v_prime(nc_loc,nv))

     ! Nonlinear arrays
     if (nonlinear_flag == 1) then
        if (nonlinear_method == 1) then
           allocate(f_nl(nc,nsplit,n_toroidal))
           allocate(g_nl(nc,nsplit,n_toroidal))
           allocate(fpack(nc,nsplit*n_toroidal))
           allocate(gpack(nc,nsplit*n_toroidal))
        else
           allocate(f_nl(n_radial,nsplit,n_toroidal))
           allocate(g_nl(n_radial,nsplit,n_toroidal))
           allocate(fpack(n_radial,nsplit*n_toroidal))
           allocate(gpack(n_radial,nsplit*n_toroidal))
        endif
     endif

     if (collision_model == 5) then
        allocate(cmat_simple(n_xi,n_xi,n_energy,n_species,n_theta))
     else
        allocate(cmat(nv,nv,nc_loc))
     endif

  endif

  ! Compute equilibrium quantities (even in test mode)
  GEO_model_in    = geo_numeq_flag
  GEO_ntheta_in   = geo_ntheta
  GEO_nfourier_in = geo_ny
  call GEO_alloc(1)
  call cgyro_equilibrium

  if (test_flag == 0) then

     call cgyro_init_arrays
     call timer_lib_out('str_init')

     call timer_lib_in('coll_init')
     call cgyro_init_collision
     call timer_lib_out('coll_init')

  endif

  call cgyro_check_memory(trim(path)//runfile_memory)

  ! Write initial data

  call timer_lib_in('io_init')
  call cgyro_write_initdata
  call timer_lib_out('io_init')

  if (test_flag == 1) return

  ! Initialize h (via restart or analytic IC)
  call timer_lib_in('str_init')
  call cgyro_init_h
  call timer_lib_out('str_init')

  ! Initialize nonlinear dimensions and arrays 
  if (nonlinear_method == 1) then

     ! Direct convolution

     ny0 = n_toroidal-1
     nx0 = n_radial/2
     ny = int(1.5*ny0)+1
     nx = int(1.5*nx0)+1

  else

     ! 2D FFT lengths 
     nx0 = n_radial
     ny0 = 2*n_toroidal-1

     ! 3/2-rule for dealiasing the nonlinear product
     nx = (3*nx0)/2
     ny = (3*ny0)/2
     
#ifndef _OPENACC
     allocate(fx(0:ny/2,0:nx-1,n_omp))
     allocate(gx(0:ny/2,0:nx-1,n_omp))
     allocate(fy(0:ny/2,0:nx-1,n_omp))
     allocate(gy(0:ny/2,0:nx-1,n_omp))

     allocate(uxmany(0:ny-1,0:nx-1,nsplit))
     allocate(uymany(0:ny-1,0:nx-1,nsplit))
     allocate(vx(0:ny-1,0:nx-1,n_omp))
     allocate(vy(0:ny-1,0:nx-1,n_omp))
     allocate(uv(0:ny-1,0:nx-1,n_omp))

     ! Create plans once and for all, with global arrays fx,ux
     plan_c2r = fftw_plan_dft_c2r_2d(nx,ny,gx(:,:,1),vx(:,:,1),FFTW_PATIENT)
     plan_r2c = fftw_plan_dft_r2c_2d(nx,ny,uv(:,:,1),fx(:,:,1),FFTW_PATIENT)
#endif

#ifdef _OPENACC
     call cgyro_info('GPU-aware code triggered.')

     allocate( fxmany(0:ny/2,0:nx-1,nsplit) )
     allocate( fymany(0:ny/2,0:nx-1,nsplit) )
     allocate( gxmany(0:ny/2,0:nx-1,nsplit) )
     allocate( gymany(0:ny/2,0:nx-1,nsplit) )

     allocate( uxmany(0:ny-1,0:nx-1,nsplit) )
     allocate( uymany(0:ny-1,0:nx-1,nsplit) )
     allocate( vxmany(0:ny-1,0:nx-1,nsplit) )
     allocate( vymany(0:ny-1,0:nx-1,nsplit) )
     allocate( uvmany(0:ny-1,0:nx-1,nsplit) )

     !-------------------------------------------------------------------
     ! 2D
     !   input[ b*idist + (x * inembed[1] + y)*istride ]
     !  output[ b*odist + (x * onembed[1] + y)*ostride ]
     !  isign is the sign of the exponent in the formula that defines
     !  Fourier transform  -1 == FFTW_FORWARD
     !                      1 == FFTW_BACKWARD
     !-------------------------------------------------------------------

     ndim(1) = nx
     ndim(2) = ny
     idist = size(fxmany,1)*size(fxmany,2)
     odist = size(uxmany,1)*size(uxmany,2)
     istride = 1
     ostride = 1
     inembed = size(fxmany,1)
     onembed = size(uxmany,1)

     istatus = cufftPlanMany(&
          cu_plan_c2r_many, &
          irank, &
          ndim, &
          inembed, &
          istride, &
          idist, &
          onembed, &
          ostride, &
          odist, &
          merge(CUFFT_C2R,CUFFT_Z2D,kind(uxmany) == singlePrecision), &
          nsplit)

     idist = size(uxmany,1)*size(uxmany,2)
     odist = size(fxmany,1)*size(fxmany,2)
     inembed = size(uxmany,1)
     onembed = size(fxmany,1) 
     istride = 1
     ostride = 1
     istatus = cufftPlanMany(&
          cu_plan_r2c_many, &
          irank, &
          ndim, &
          inembed, &
          istride, &
          idist, &
          onembed, &
          ostride, &
          odist, &
          merge(CUFFT_R2C,CUFFT_D2Z,kind(uxmany) == singlePrecision), &
          nsplit)
#endif

  endif

end subroutine cgyro_init_manager

!---------------------------------------------------
! Calculate weight correction to integrate function
! over interval [0,inf] instead of [0,b]
!
! Do this by calculating offsets to match
!
! Int[1/u^2], Int[1/u], Int[1], Int[u], ...
!
! We start from 1/u^2 not 1 since these are the
! total polynomials of the original scheme without
! u^2 weighting.
!----------------------------------------------------

subroutine domain_renorm(u,w,n)

  implicit none

  integer, intent(in) :: n
  real, intent(in), dimension(n) :: u
  real, intent(inout), dimension(n) :: w
  integer, dimension(:), allocatable :: i_piv
  real, dimension(:,:), allocatable :: a
  real, dimension(:), allocatable :: b
  real :: pi
  integer :: info,m,i

  pi = 4*atan(1.0)

  allocate(i_piv(n))
  allocate(b(n))
  allocate(a(n,n))
  do m=0,n-1
     b(m+1) = 2*gamma((1+m)/2.0)/sqrt(pi)-sum(w*u**(m-2))    
     do i=1,n
        a(m+1,i) = u(i)**(m-2)
     enddo
  enddo

  call DGESV(n,1,a,n,i_piv,b,n,info)

  w = w+b

end subroutine domain_renorm
