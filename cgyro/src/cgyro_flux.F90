!---------------------------------------------------------------------------
! cgyro_flux.f90
!
! PURPOSE:
!  Compute flux and mean-square fluctuation amplitudes as
!  functions of (kx,ky).
!
! NOTES:
!
!  Flux definitions (these are all broken down by ky) 
!
!   gflux(l,ky) = Radial Fourier expansion of flux (l is the Fourier harmonic)
!   gflux(0,ky) = Domain-averaged flux [this is the standard flux-tube flux]
!     cflux(ky) = Interior-averaged flux [this is the "central" flux]
!
!  Global flux logic:
!
!   Retain radial harmonics l = [0,...,N_GLOBAL], where N_GLOBAL=4 by default
!
!   To trigger output of the global flux (viewable by -plot xflux),
!    set GFLUX_PRINT_FLAG=1
!
! NOTE:
!  cflux is only of interest when profile or ExB shear is active
!  cflux(ky) is the average flux in the "positive-shear region"
!  gflux(ky,0)-cflux is the flux in the "negative-shear region"
!---------------------------------------------------------------------------

module cgyro_flux_mod

  implicit none

  complex, private, dimension(:,:,:,:,:), allocatable :: moment_loc
  complex, dimension(:,:,:,:,:), allocatable :: moment

  ! Nonlinear fluxes (f=standard,c=central,g=global)
  real, private, dimension(:,:,:,:), allocatable :: cflux_loc
  real, dimension(:,:,:,:), allocatable :: cflux
  complex, private, dimension(:,:,:,:,:), allocatable :: gflux_loc
  complex, dimension(:,:,:,:,:), allocatable :: gflux
  real, dimension(:,:), allocatable :: cflux_tave, gflux_tave

  real, private :: tave_min, tave_max
  integer, private :: tave_step

  ! Not accessible on GPU
  complex, dimension(:,:,:), allocatable :: cap_h_c_cur
  ! internal work buffer
  complex, private, dimension(:,:,:), allocatable :: cap_h_c_dot
  ! cap_h_c history
  complex, private, dimension(:,:,:), allocatable :: cap_h_c_old
  complex, private, dimension(:,:,:), allocatable :: cap_h_c_old2
  complex, private, dimension(:,:,:), allocatable :: cap_h_c_old3
  logical, private :: cap_h_c_old_valid

contains

subroutine cgyro_flux_tave_reset(t_current)

  implicit none

  real, intent(in) :: t_current

  tave_step  = 0
  tave_min   = t_current
  tave_max   = t_current
  cflux_tave = 0.0
  gflux_tave = (0.0,0.0)

end subroutine cgyro_flux_tave_reset


!-----------------------------------------------------------------
! Initialization and cleanup, to be called once
!-----------------------------------------------------------------

subroutine cgyro_flux_init(n_radial,theta_plot,n_species,n_field,&
                           n_global,nc,nv_loc,nt1,nt2)

  implicit none

  integer, intent(in) :: n_radial,theta_plot,n_species,n_field,&
                         n_global,nc,nv_loc,nt1,nt2

  allocate(cap_h_c_dot(nc,nv_loc,nt1:nt2))
  allocate(cap_h_c_old(nc,nv_loc,nt1:nt2))
  allocate(cap_h_c_old2(nc,nv_loc,nt1:nt2))
  allocate(cap_h_c_old3(nc,nv_loc,nt1:nt2))
#if defined(OMPGPU)
!$omp target enter data map(alloc:cap_h_c_dot,cap_h_c_old,cap_h_c_old2,cap_h_c_old3)
#elif defined(_OPENACC)
!$acc enter data create(cap_h_c_dot,cap_h_c_old,cap_h_c_old2,cap_h_c_old3)
#endif
  cap_h_c_old_valid = .FALSE.

  allocate(cap_h_c_cur(nc,nv_loc,nt1:nt2))

  allocate(    moment(n_radial,theta_plot,n_species,nt1:nt2,3))
  allocate(moment_loc(n_radial,theta_plot,n_species,nt1:nt2,3))
  allocate(    cflux(n_species,4,n_field,nt1:nt2))
  allocate(cflux_loc(n_species,4,n_field,nt1:nt2))
  allocate(    gflux(0:n_global,n_species,4,n_field,nt1:nt2))
  allocate(gflux_loc(0:n_global,n_species,4,n_field,nt1:nt2))
  allocate(cflux_tave(n_species,4))
  allocate(gflux_tave(n_species,4))

  call cgyro_flux_tave_reset(0.0)

end subroutine cgyro_flux_init

subroutine cgyro_flux_cleanup

  implicit none

  if(allocated(moment))              deallocate(moment)
  if(allocated(moment_loc))          deallocate(moment_loc)
  if(allocated(cflux))               deallocate(cflux)
  if(allocated(cflux_loc))           deallocate(cflux_loc)
  if(allocated(gflux))               deallocate(gflux)
  if(allocated(gflux_loc))           deallocate(gflux_loc)
  if(allocated(cflux_tave))          deallocate(cflux_tave)
  if(allocated(gflux_tave))          deallocate(gflux_tave)

  if(allocated(cap_h_c_dot)) then ! one check enough
     deallocate(cap_h_c_cur)
#if defined(OMPGPU)
!$omp target exit data map(release:cap_h_c_dot,cap_h_c_old,cap_h_c_old2,cap_h_c_old3)
#elif defined(_OPENACC)
!$acc exit data delete(cap_h_c_dot,cap_h_c_old,cap_h_c_old2,cap_h_c_old3)
#endif
     deallocate(cap_h_c_old3)
     deallocate(cap_h_c_old2)
     deallocate(cap_h_c_old)
     deallocate(cap_h_c_dot)
  endif

end subroutine cgyro_flux_cleanup

subroutine cgyro_flux_save_cap_h_c

  use cgyro_globals, only: nt1,nt2,nv1,nv2,nc,cap_h_c

  implicit none

  integer :: itor,iv,ic,iv_loc

  if (cap_h_c_old_valid) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   private(iv_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector private(iv_loc) &
!$acc&         present(cap_h_c,cap_h_c_old,cap_h_c_old2,cap_h_c_old3) &
!$acc&         present(nt1,nt2,nv1,nv2,nc) default(none)
#else
!$omp parallel do collapse(3) private(iv_loc)
#endif
    do itor=nt1,nt2
      do iv=nv1,nv2
        do ic=1,nc
           iv_loc = iv-nv1+1
           cap_h_c_old3(ic,iv_loc,itor) = cap_h_c_old2(ic,iv_loc,itor)
           cap_h_c_old2(ic,iv_loc,itor) = cap_h_c_old(ic,iv_loc,itor)
           cap_h_c_old(ic,iv_loc,itor) = cap_h_c(ic,iv_loc,itor)
        enddo
      enddo
    enddo
  else
    ! first time around, just put cap_h_c in all
    cap_h_c_old_valid = .TRUE.
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   private(iv_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector private(iv_loc) &
!$acc&         present(cap_h_c,cap_h_c_old,cap_h_c_old2,cap_h_c_old3) &
!$acc&         present(nt1,nt2,nv1,nv2,nc) default(none)
#else
!$omp parallel do collapse(3) private(iv_loc)
#endif
    do itor=nt1,nt2
      do iv=nv1,nv2
        do ic=1,nc
           iv_loc = iv-nv1+1
           cap_h_c_old3(ic,iv_loc,itor) = cap_h_c(ic,iv_loc,itor)
           cap_h_c_old2(ic,iv_loc,itor) = cap_h_c(ic,iv_loc,itor)
           cap_h_c_old(ic,iv_loc,itor) = cap_h_c(ic,iv_loc,itor)
        enddo
      enddo
    enddo
  endif

end subroutine cgyro_flux_save_cap_h_c

! make a copy from cap_h_c to cap_h_c_cur
subroutine cgyro_flux_sync_cap_h_c_cur

  use cgyro_globals, only : cap_h_c

   implicit none

#if defined(OMPGPU)
!$omp target update from(cap_h_c)
#elif defined(_OPENACC)
!$acc update host(cap_h_c)
#endif
  cap_h_c_cur = cap_h_c

end subroutine cgyro_flux_sync_cap_h_c_cur

!-----------------------------------------------------------------

subroutine cgyro_flux

  use mpi
  use cgyro_globals, only : bmag, bigr, btor, delta_t, dens2_rot, &
       energy, i_c, i_err, ic_c, ie_v, ir_c, is_v, it_c, itp, &
       ix_v, jvec_c, jxvec_c, k_theta_base, lambda_rot, mach, &
       mass, nc, NEW_COMM_1, n_field, n_global, nonlinear_flag, &
       n_radial, nt1, nt2, nv1, nv2, pi, rho, rmaj, t_current, &
       temp, vel2, vth, w_exi, w_theta, xi, z
  use cgyro_field_mod, only : field_cur, field_dot

  implicit none

  integer :: ic,ie,ix,is,it,ir,i_field,itor,iv,iv_loc
  integer :: l,icl
  real :: dv,cn
  real :: vpar
  complex, dimension(0:n_global,n_field) :: prod1,prod2,prod3
  real :: dvr
  real :: erot
  real :: flux_norm
  complex :: cprod
  real, parameter :: x_fraction=0.2
  real :: u

  ! cap_h_c_old* are in GPU memory only, so compute cap_h_c_dot there
  ! note: We use explicitly a pre-allocated buffer in GPU memory + update from,
  !       to avoid dynamic memory allocation
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   private(iv_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector private(iv_loc) &
!$acc&         present(cap_h_c_dot,cap_h_c_old,cap_h_c_old2,cap_h_c_old3) &
!$acc&         present(nt1,nt2,nv1,nv2,nc) copyin(delta_t) default(none)
#else
!$omp parallel do collapse(3) private(iv_loc)
#endif
  do itor=nt1,nt2
     do iv=nv1,nv2
        do ic=1,nc
           iv_loc = iv-nv1+1
           cap_h_c_dot(ic,iv_loc,itor) = (3*cap_h_c_old(ic,iv_loc,itor) - &
                4*cap_h_c_old2(ic,iv_loc,itor) + &
                cap_h_c_old3(ic,iv_loc,itor) )/(2*delta_t)
        enddo
     enddo
  enddo
#if defined(OMPGPU)
!$omp target update from(cap_h_c_dot)
#elif defined(_OPENACC)
!$acc update host(cap_h_c_dot)
#endif

!$omp parallel do private(iv_loc,iv,is,ix,ie,dv,vpar,ic,ir,it,erot,cprod,cn) &
!$omp&            private(prod1,prod2,prod3,l,icl,dvr,u,flux_norm) &
!$omp&            shared(moment_loc,gflux_loc,cflux_loc)
  do itor=nt1,nt2

     !-----------------------------------------------------
     ! 1. Compute kx-ky moments (n,E,v)
     !-----------------------------------------------------

     moment_loc(:,:,:,itor,:) = 0.0

     iv_loc = 0
     do iv=nv1,nv2

        iv_loc = iv_loc+1

        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        ! Integration weight
        dv = w_exi(ie,ix)

        ! Parallel velocity
        vpar = vth(is)*vel2(ie)*xi(ix)

        do ic=1,nc

           ir = ir_c(ic)
           it = it_c(ic)

           erot  = (energy(ie)+lambda_rot(it,is))*temp(is)

           if (itp(it) > 0) then
              ! dvjvec_c(:,ic,iv_loc,itor) = dens2_rot(it,is)*w_exi(ie,ix)*z(is)*jvec_c(:,ic,iv_loc,itor)
              ! dvjvec_c(1,ic,iv_loc,itor)/z(is) = dens2_rot(it,is)*w_exi(ie,ix)*jvec_c(1,ic,iv_loc,itor)
              cprod = cap_h_c_cur(ic,iv_loc,itor)*dens2_rot(it,is)*w_exi(ie,ix)*jvec_c(1,ic,iv_loc,itor)
              cn    = dv*z(is)*dens2_rot(it,is)/temp(is)

              ! Note addition division by rho below
              
              ! Density moment: (delta n_a)/(n_norm)
              moment_loc(ir,itp(it),is,itor,1) = moment_loc(ir,itp(it),is,itor,1)-(cn*field_cur(1,ic,itor)-cprod)

              ! Energy moment : (delta E_a)/(n_norm T_norm)
              moment_loc(ir,itp(it),is,itor,2) = moment_loc(ir,itp(it),is,itor,2)-(cn*field_cur(1,ic,itor)-cprod)*erot

              ! Velocity moment : (delta v_a)/(n_norm v_norm)
              moment_loc(ir,itp(it),is,itor,3) = moment_loc(ir,itp(it),is,itor,3)-(cn*field_cur(1,ic,itor)-cprod)*vpar
           endif

        enddo
     enddo

     ! Unlike moments, fields are *not* divided by rho
     ! (cgyro_plot will adjust so both have same normalization)
     moment_loc(:,:,:,itor,:) = moment_loc(:,:,:,itor,:)/rho

     !-------------------------------------------------------------
     ! 2. Compute global ky-dependent fluxes (with field breakdown)
     !-------------------------------------------------------------

     gflux_loc(:,:,:,:,itor) = 0.0

     iv_loc = 0
     do iv=nv1,nv2

        iv_loc = iv_loc+1

        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        ! Integration weight
        dv = w_exi(ie,ix)

        ! Parallel velocity
        vpar = vth(is)*vel2(ie)*xi(ix)

        do ic=1,nc

           ir = ir_c(ic)
           it = it_c(ic)

           ! prod* are local to this loop
           prod1 = 0.0 
           prod2 = 0.0
           prod3 = 0.0

           ! Global flux coefficients (complex coefficients required to compute radial profile)
           do l=0,n_global

              ! H w^* + H^* w

              if (ir-l > 0) then
                 icl = ic_c(ir-l,it)
                 prod1(l,:) = prod1(l,:) &
                      +i_c*cap_h_c_cur(ic,iv_loc,itor)*conjg(jvec_c(:,icl,iv_loc,itor)*field_cur(:,icl,itor))
                 prod2(l,:) = prod2(l,:) &
                      +i_c*cap_h_c_cur(ic,iv_loc,itor)*conjg(i_c*jxvec_c(:,icl,iv_loc,itor)*field_cur(:,icl,itor))
                 prod3(l,:) = prod3(l,:) &
                      -cap_h_c_dot(ic,iv_loc,itor)*conjg(jvec_c(:,icl,iv_loc,itor)*field_cur(:,icl,itor)) &
                      +cap_h_c_cur(ic,iv_loc,itor)*conjg(jvec_c(:,icl,iv_loc,itor)*field_dot(:,icl,itor))
              endif
              if (ir+l <= n_radial) then
                 icl = ic_c(ir+l,it)
                 prod1(l,:) = prod1(l,:) &
                      -i_c*conjg(cap_h_c_cur(ic,iv_loc,itor))*jvec_c(:,icl,iv_loc,itor)*field_cur(:,icl,itor)
                 prod2(l,:) = prod2(l,:) &
                      -i_c*conjg(cap_h_c_cur(ic,iv_loc,itor))*i_c*jxvec_c(:,icl,iv_loc,itor)*field_cur(:,icl,itor)
                 prod3(l,:) = prod3(l,:) &
                      -conjg(cap_h_c_dot(ic,iv_loc,itor))*jvec_c(:,icl,iv_loc,itor)*field_cur(:,icl,itor) &
                      +conjg(cap_h_c_cur(ic,iv_loc,itor))*jvec_c(:,icl,iv_loc,itor)*field_dot(:,icl,itor)
              endif

           enddo

           dvr  = w_theta(it)*dens2_rot(it,is)*dv
           erot = (energy(ie)+lambda_rot(it,is))*temp(is)

           ! 1. Density flux: Gamma_a
           gflux_loc(:,is,1,:,itor) = gflux_loc(:,is,1,:,itor)+prod1(:,:)*dvr

           ! 2. Energy flux : Q_a
           gflux_loc(:,is,2,:,itor) = gflux_loc(:,is,2,:,itor)+prod1(:,:)*dvr*erot

           prod1(:,:) = prod1(:,:)*(mach*bigr(it)/rmaj+btor(it)/bmag(it)*vpar)+prod2(:,:)

           ! 3. Momentum flux: Pi_a
           gflux_loc(:,is,3,:,itor) = gflux_loc(:,is,3,:,itor)+prod1(:,:)*dvr*bigr(it)*mass(is)

           ! 4. Exchange
           gflux_loc(:,is,4,:,itor) = gflux_loc(:,is,4,:,itor)+0.5*prod3(:,:)*dvr*z(is)

        enddo

     enddo

     ! Construct "positive/interior" flux (real quantities)
     cflux_loc(:,:,:,itor) = real(gflux_loc(0,:,:,:,itor))
     do l=1,n_global
        u = 2*pi*l*x_fraction
        cflux_loc(:,:,:,itor) = cflux_loc(:,:,:,itor)+2*sin(u)*real(gflux_loc(l,:,:,:,itor))/u
     enddo

     !-----------------------------------------------------
     ! 3. Renormalize fluxes to GB or quasilinear forms
     !~----------------------------------------------------

     if (nonlinear_flag == 0 .and. itor > 0) then

        ! Quasilinear normalization (divide by |phi|^2)
        ! Note: We assume we compute flux_norm once per itor
        flux_norm = 0.0
        do ir=1,n_radial
           flux_norm = flux_norm+sum(abs(field_cur(1,ic_c(ir,:),itor))**2*w_theta(:))
        enddo

        ! Sign correction fixed 2025.11.19 (JC)
        ! need 2 for regression compatibility
        flux_norm = 2*flux_norm*sign(1.0,k_theta_base*rho)
        gflux_loc(:,:,:,:,itor) = gflux_loc(:,:,:,:,itor)/flux_norm 
        cflux_loc (:,:,:,itor)  = cflux_loc(:,:,:,itor)/flux_norm 

     else

        ! Complete definition of fluxes (not exchange)
        gflux_loc(:,:,1:3,:,itor) = gflux_loc(:,:,1:3,:,itor)*(k_theta_base*itor*rho)
        cflux_loc(:,1:3,:,itor)   = cflux_loc(:,1:3,:,itor)*(k_theta_base*itor*rho)

        ! GyroBohm normalizations
        gflux_loc(:,:,:,:,itor) = gflux_loc(:,:,:,:,itor)/rho**2
        cflux_loc(:,:,:,itor) = cflux_loc(:,:,:,itor)/rho**2

     endif
  enddo

  ! Reduced complex moment(kx,ky), below, is still distributed over n 

  call MPI_ALLREDUCE(moment_loc(:,:,:,:,:), &
       moment(:,:,:,:,:), &
       size(moment), &
       MPI_DOUBLE_COMPLEX, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  ! Global complex gflux(l,ky), below, is still distributed over n 

  call MPI_ALLREDUCE(gflux_loc, &
       gflux, &
       size(gflux), &
       MPI_DOUBLE_COMPLEX, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  ! Reduced real cflux(ky), below, is still distributed over n 

  call MPI_ALLREDUCE(cflux_loc, &
       cflux, &
       size(cflux), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)

  tave_step = tave_step + 1
  tave_max  = t_current
  do itor=nt1,nt2
     do i_field=1,n_field
        cflux_tave(:,:) = cflux_tave(:,:) + cflux(:,:,i_field,itor)
     enddo
     do i_field=1,n_field
        gflux_tave(:,:) = gflux_tave(:,:) + real(gflux(0,:,:,i_field,itor))
     enddo
  enddo
     
end subroutine cgyro_flux

subroutine cgyro_flux_sum(tave_min_out, tave_max_out, flux_tave_out)
  use mpi
  use cgyro_globals, only: n_species,gamma_e,NEW_COMM_2

  implicit none

  real, intent(out)     :: tave_min_out, tave_max_out
  real, intent(out)     :: flux_tave_out(11,3)

  integer :: i_err
  real, dimension(n_species,3) :: sum_out

  ! Return time-averaged flux data (need to reduce across n first)
  flux_tave_out(:,:) = 0.0
  if (abs(gamma_e) > 1e-10) then
     call MPI_ALLREDUCE(cflux_tave(:,:), &
       sum_out(:,:), &
       size(sum_out), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)
  else
     call MPI_ALLREDUCE(gflux_tave(:,:), &
       sum_out(:,:), &
       size(sum_out), &
       MPI_DOUBLE_PRECISION, &
       MPI_SUM, &
       NEW_COMM_2, &
       i_err)
  endif
  
  tave_min_out = tave_min
  tave_max_out = tave_max
  flux_tave_out(1:n_species,:) = sum_out(1:n_species,:)/tave_step

end subroutine cgyro_flux_sum

end module cgyro_flux_mod
