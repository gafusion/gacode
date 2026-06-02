!-----------------------------------------------------------------
! cgyro_field.f90
!
! PURPOSE:
!  Perform field solves for data distributed in both the velocity 
!  and configuration indices:
!  
!  - cgyro_field_v_notae() (configuration distributed)
!  - cgyro_field_c() (velocity distributed)
!-----------------------------------------------------------------

module cgyro_field_mod

  implicit none

  !
  ! Field_c

  ! Note: GPU only on GPU systems
  complex, dimension(:,:,:), allocatable :: field
  !
  real, private, dimension(:,:,:), allocatable :: fcoef
  real, private, dimension(:,:,:), allocatable :: gcoef
  complex, private, dimension(:,:,:), allocatable :: field_loc
  !
  ! Field_v
  complex, private, dimension(:,:,:,:), allocatable :: field_v
  complex, private, dimension(:,:,:,:), allocatable :: field_loc_v
  !
  ! Field_e - used for error estimation, fluxes and timedata

  ! Note: Cannot be used on GPU
  complex, dimension(:,:,:), allocatable :: field_dot
  complex, dimension(:,:,:), allocatable :: field_cur
  !
  complex, private, dimension(:,:,:), allocatable :: field_old
  complex, private, dimension(:,:,:), allocatable :: field_old2
  complex, private, dimension(:,:,:), allocatable :: field_old3

contains

!-----------------------------------------------------------------
! Initialization and cleanup, to be called once
!-----------------------------------------------------------------

subroutine cgyro_field_c_init(n_field,nc,nt1,nt2)

  implicit none

  integer, intent(in) :: n_field,nc,nt1,nt2

  allocate(field(n_field,nc,nt1:nt2))
  ! private buffers
  allocate(fcoef(n_field,nc,nt1:nt2))
  if (n_field < 3) then
     allocate(gcoef(n_field,nc,nt1:nt2))
  else
     allocate(gcoef(5,nc,nt1:nt2))
  endif
  allocate(field_loc(n_field,nc,nt1:nt2))
#if defined(OMPGPU)
!$omp target enter data map(alloc:fcoef,gcoef,field,field_loc)
#elif defined(_OPENACC)
!$acc enter data create(fcoef,gcoef,field,field_loc)
#endif
end subroutine cgyro_field_c_init

subroutine cgyro_field_e_init(n_field,nc,nt1,nt2)

  implicit none

  integer, intent(in) :: n_field,nc,nt1,nt2

  ! These are in CPU-memory only
  allocate(field_dot(n_field,nc,nt1:nt2))
  allocate(field_cur(n_field,nc,nt1:nt2))
  allocate(field_old(n_field,nc,nt1:nt2))
  allocate(field_old2(n_field,nc,nt1:nt2))
  allocate(field_old3(n_field,nc,nt1:nt2))

  ! initialize field_cur with values of field
  ! But field only have values on GPU, so fetch it from there
#if defined(OMPGPU)
!$omp target update from(field)
#elif defined(_OPENACC)
!$acc update host(field)
#endif
  field_cur = field

  ! Initialize time-history of fields (-3,-2,-1) to initial field.
  field_old  = field_cur
  field_old2 = field_cur
  field_old3 = field_cur
end subroutine cgyro_field_e_init

subroutine cgyro_field_v_init(n_field,nc,nt1,nt2,n_sim,nc_cl1,nc_cl2)

  implicit none

  integer, intent(in) :: n_field,nc,nt1,nt2,n_sim,nc_cl1,nc_cl2

  ! nc and nc_loc_coll must be last, since it will be collated     
  allocate(field_v(n_field,nt1:nt2,n_sim,nc))
  allocate(field_loc_v(n_field,nt1:nt2,n_sim,nc_cl1:nc_cl2))
#if defined(OMPGPU)
!$omp target enter data map(alloc:field_v,field_loc_v)
#elif defined(_OPENACC)
!$acc enter data create(field_v,field_loc_v)
#endif

end subroutine cgyro_field_v_init


subroutine cgyro_field_c_cleanup

  implicit none

  if(allocated(fcoef))  then ! one test enough
#if defined(OMPGPU)
!$omp target exit data map(release:fcoef,gcoef,field,field_loc)
#elif defined(_OPENACC)
!$acc exit data delete(fcoef,gcoef,field,field_loc)
#endif
    deallocate(field_loc)
    deallocate(field)
    deallocate(gcoef)
    deallocate(fcoef)
  endif
end subroutine cgyro_field_c_cleanup

subroutine cgyro_field_e_cleanup

  implicit none

  if(allocated(field_dot))  then ! one test enough
    deallocate(field_old3)
    deallocate(field_old2)
    deallocate(field_old)
    deallocate(field_cur)
    deallocate(field_dot)
  endif
end subroutine cgyro_field_e_cleanup

subroutine cgyro_field_v_cleanup

  implicit none

  if(allocated(field_v))  then ! one test enough
#if defined(OMPGPU)
!$omp target exit data map(release:field_v,field_loc_v)
#elif defined(_OPENACC)
!$acc exit data delete(field_v,field_loc_v)
#endif
    deallocate(field_loc_v)
    deallocate(field_v)
  endif
end subroutine cgyro_field_v_cleanup


!-----------------------------------------------------------------
! Velocity (configuration-distributed) field solve
!-----------------------------------------------------------------

! like cgyro_field_v_notae, but with parametrized start_t

subroutine cgyro_field_v_notae_s(start_t)
  use parallel_lib, only : fsendf, parallel_lib_collect_field, &
       parallel_lib_nj_loc
  use timer_lib, only : timer_lib_in, timer_lib_out
  use cgyro_globals, only : dvjvec_v, i_sim, nc, nc_cl1, nc_cl2, &
       n_field, n_sim, nt2, nv, nv_loc

  implicit none
  ! ------------------ 
  integer, intent(in) :: start_t
  !
  integer :: ic,ic_loc,iv,i_f,itor,j,k,nj_loc,ism
  integer :: vcount
  complex :: field_loc_l 

  call timer_lib_in('field')
  vcount = nv/nv_loc

  ! iv = j+(k-1)*nj_loc
  ! cap_h_v(ic_loc,itor,iv) = fsendf(j,itor,ic_loc,k)
  ! Use aggregate comm, since it is used in step_collision
  call parallel_lib_nj_loc(nj_loc)

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc data present(fsendf)
!$acc data present(field,field_v,field_loc_v)
#endif

  ! Poisson and Ampere RHS integrals of H

#if defined(OMPGPU)
!$omp target teams distribute collapse(4) firstprivate(start_t,nj_loc,vcount) &
!$omp&       private(ic_loc,field_loc_l) 
#elif defined(_OPENACC)
!$acc parallel loop collapse(4) gang private(ic_loc,field_loc_l) &
!$acc&         present(dvjvec_v,fsendf) firstprivate(start_t,nj_loc,vcount) &
!$acc&         firstprivate(nt2,nc_cl1,nc_cl2,n_field,nv,n_sim) default(none)
#else
!$omp parallel do collapse(4) firstprivate(start_t,nj_loc,vcount) &
!$omp&       private(ic_loc,field_loc_l,j,k,iv) 
#endif
  do ic=nc_cl1,nc_cl2
   do ism=1,n_sim
    do itor=start_t,nt2
     do i_f=1,n_field
      ic_loc = ic-nc_cl1+1
      field_loc_l = (0.0,0.0)
#if defined(OMPGPU)
!$omp parallel do simd collapse(2) reduction(+:field_loc_l) &
!$omp&       private(iv)
#elif defined(_OPENACC)
!$acc loop vector collapse(2) private(iv) reduction(+:field_loc_l)
#else
      ! just leave the compiler to vectorize, if it can
#endif
      do k=1,vcount
       do j=1,nj_loc
        iv = j+(k-1)*nj_loc
        field_loc_l = field_loc_l+dvjvec_v(i_f,iv,itor,ic_loc)*fsendf(j,itor,ic_loc,k+(ism-1)*vcount)
      enddo
      enddo
      field_loc_v(i_f,itor,ism,ic) = field_loc_l
     enddo
    enddo
   enddo
  enddo

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  ! Use aggregate comm, since it is used in step_collision
  call parallel_lib_collect_field(field_loc_v, field_v)

  call timer_lib_out('field_com')

  call timer_lib_in('field')
  ! Poisson LHS factors
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&    firstprivate(start_t,i_sim)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector &
!$acc&         independent present(fcoef) firstprivate(start_t,i_sim) &
!$acc&         present(nt2,nc,n_field) default(none)
#else
!$omp parallel do collapse(3) &
!$omp&    firstprivate(start_t,i_sim)
#endif
  do itor=start_t,nt2
     ! assuming  (.not.(itor == 0 .and. ae_flag == 1))
     do ic=1,nc
       do i_f=1,n_field
        field(i_f,ic,itor) = fcoef(i_f,ic,itor)*field_v(i_f,itor,i_sim,ic)
       enddo
     enddo
  enddo

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc end data
!$acc end data
#endif

  call timer_lib_out('field')

end subroutine cgyro_field_v_notae_s

! like cgyro_field_v, but skip (itor == 0 .and. ae_flag == 1)
subroutine cgyro_field_v_notae

  use cgyro_globals, only : ae_flag, nt1, nt2

  implicit none

  if (nt1 == 0 .and. ae_flag == 1) then
     if (nt2>0) then
        call cgyro_field_v_notae_s(1)
     endif
     ! else no-op
  else
     ! don't have to worry about ae_flag, just use all the elements
     call cgyro_field_v_notae_s(nt1)
  endif

end subroutine cgyro_field_v_notae

!-----------------------------------------------------------------
! Configuration (velocity-distributed) field solve
!-----------------------------------------------------------------

subroutine cgyro_field_c(update_cap)
  use parallel_lib, only : parallel_flib_sum_field
  use timer_lib, only : timer_lib_in, timer_lib_out
  use cgyro_globals, only : ae_flag, cap_h_c, dvjvec_c, h_x, is_v, &
       jvec_c, nc, n_field, nt1, nt2, nv1, nv2, temp, z
  implicit none

  logical, intent(in) :: update_cap

  integer :: ic,is,iv,iv_loc,i_f,itor
  integer :: itor1,itor2
  complex :: tmp,field_loc_l
  complex :: my_psi

  call timer_lib_in('field')

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc data present(h_x,cap_h_c)
!$acc data present(field,field_loc)
#endif

  ! Poisson and Ampere RHS integrals of h

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   private(field_loc_l,iv_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector &
!$acc&         independent private(field_loc_l) &
!$acc&         present(dvjvec_c) present(nt1,nt2,nc,n_field,nv1,nv2) default(none)
#else
!$omp parallel do collapse(3) &
!$omp&   private(field_loc_l,iv_loc)
#endif
  do itor=nt1,nt2
   do ic=1,nc
    do i_f=1,n_field
      field_loc_l = (0.0,0.0)    
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq private(iv_loc)
#endif
      do iv=nv1,nv2
         iv_loc = iv-nv1+1
         field_loc_l = field_loc_l+dvjvec_c(i_f,ic,iv_loc,itor)*h_x(ic,iv_loc,itor)
      enddo
      field_loc(i_f,ic,itor) = field_loc_l
    enddo
   enddo
  enddo
  call timer_lib_out('field')
  call timer_lib_in('field_com')

  call parallel_flib_sum_field(field_loc,field)

  call timer_lib_out('field_com')
  call timer_lib_in('field')
  if (n_field > 2) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(2)
#elif defined(_OPENACC)
!$acc parallel loop collapse(2) gang vector &
!$acc&         independent present(fcoef) &
!$acc&         present(nt1,nt2,nc) default(none)
#else
!$omp parallel do collapse(2)
#endif
    do itor=nt1,nt2
      do ic=1,nc
       field(3,ic,itor) = field(3,ic,itor)*fcoef(3,ic,itor)
      enddo
     enddo
  endif

  ! Poisson LHS factors
  itor1=nt1
  itor2=nt2

  if (nt1 == 0 .and. ae_flag == 1) then
    call cgyro_field_ae_c
    ! mark we already processed ==0, nothing else to do there
    itor1=1
  endif

  if (itor1<=itor2) then
     if (n_field > 2) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(2) &
!$omp&   private(tmp) map(to:itor1,itor2)
#elif defined(_OPENACC)
!$acc parallel loop collapse(2) gang vector &
!$acc&         independent private(tmp) present(gcoef) &
!$acc&         copyin(itor1,itor2) present(nc) default(none)
#else
!$omp parallel do collapse(2) &
!$omp&   private(tmp)
#endif
        do itor=itor1,itor2
         do ic=1,nc
          tmp = field(1,ic,itor)
          field(1,ic,itor) = gcoef(1,ic,itor)*field(1,ic,itor)+ &
                  gcoef(4,ic,itor)*field(3,ic,itor)
          field(2,ic,itor) = gcoef(2,ic,itor)*field(2,ic,itor)
          field(3,ic,itor) = gcoef(3,ic,itor)*field(3,ic,itor)+ &
                  gcoef(5,ic,itor)*tmp
         enddo
        enddo
     else
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   map(to:itor1,itor2)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector &
!$acc&         independent present(gcoef) &
!$acc&         copyin(itor1,itor2) present(nc,n_field) default(none)
#else
!$omp parallel do collapse(3)
#endif
        do itor=itor1,itor2
         do ic=1,nc
          do i_f=1,n_field
            field(i_f,ic,itor) = gcoef(i_f,ic,itor)*field(i_f,ic,itor)
          enddo
         enddo
        enddo
     endif
  endif

  if (update_cap) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   private(iv_loc,is,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector private(iv_loc,is,my_psi) &
!$acc&         present(jvec_c,z,temp,is_v) present(nt1,nt2,nv1,nv2,nc) default(none)
#else
!$omp parallel do collapse(3) &
!$omp&   private(iv_loc,is,my_psi)
#endif
   do itor=nt1,nt2
    do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
        cap_h_c(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+my_psi*z(is)/temp(is)
     enddo
    enddo
   enddo
  endif

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc end data
!$acc end data
#endif

  call timer_lib_out('field')
end subroutine cgyro_field_c

! like cgyro_field_c, but assume (my_toroidal == 0 .and. ae_flag == 1)
subroutine cgyro_field_c_ae
  use parallel_lib, only : parallel_flib_sum_field
  use timer_lib, only : timer_lib_in, timer_lib_out
  use cgyro_globals, only : cap_h_c, dvjvec_c, h_x, is_v, jvec_c, &
       nc, n_field, nv1, nv2, temp, z
  implicit none
  integer :: ic,is,iv,iv_loc,i_f,itor
  complex :: tmp,field_loc_l
  complex :: my_psi

  call timer_lib_in('field')

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc data present(h_x,cap_h_c)
!$acc data present(field,field_loc)
#endif

  ! Poisson and Ampere RHS integrals of h

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   private(field_loc_l,iv_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector &
!$acc&         independent private(field_loc_l) &
!$acc&         present(dvjvec_c) present(nc,n_field,nv1,nv2) default(none)
#else
!$omp parallel do collapse(3) &
!$omp&   private(field_loc_l,iv_loc)
#endif
  do itor=0,0
   do ic=1,nc
    do i_f=1,n_field
      field_loc_l = (0.0,0.0)    
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc loop seq private(iv_loc)
#endif
      do iv=nv1,nv2
         iv_loc = iv-nv1+1
         field_loc_l = field_loc_l+dvjvec_c(i_f,ic,iv_loc,itor)*h_x(ic,iv_loc,itor)
      enddo
      field_loc(i_f,ic,itor) = field_loc_l
    enddo
   enddo
  enddo
  call timer_lib_out('field')
  call timer_lib_in('field_com')

  call parallel_flib_sum_field(field_loc(:,:,0:0), &
                              field(:,:,0:0))

  call timer_lib_out('field_com')
  call timer_lib_in('field')
  if (n_field > 2) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(2)
#elif defined(_OPENACC)
!$acc parallel loop collapse(2) gang vector &
!$acc&         independent present(fcoef) present(nc) default(none)
#else
!$omp parallel do simd collapse(2)
#endif
    do itor=0,0
      do ic=1,nc
       field(3,ic,itor) = field(3,ic,itor)*fcoef(3,ic,itor)
      enddo
     enddo
  endif

  ! Poisson LHS factors
  call cgyro_field_ae_c

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   private(iv_loc,is,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector private(iv_loc,is,my_psi) &
!$acc&         present(jvec_c,z,temp,is_v) present(nv1,nv2,nc) default(none)
#else
!$omp parallel do collapse(3) &
!$omp&   private(iv_loc,is,my_psi)
#endif
  do itor=0,0
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
        cap_h_c(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+my_psi*z(is)/temp(is)
     enddo
   enddo
  enddo

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc end data
!$acc end data
#endif

  call timer_lib_out('field')
end subroutine cgyro_field_c_ae

!-----------------------------------------------------------------
! Adiabatic electron field solves for n=0
! Can only be called if itor==0
!-----------------------------------------------------------------
subroutine cgyro_field_ae_c

  use cgyro_globals, only : n_radial, n_theta, px, xzf, &
       zf_test_mode

  implicit none

  integer :: ir,i,j
  complex, dimension(n_theta) :: pvec_in
  complex :: pvec_out

  ! ic_c(ir,it) = (ir-1)*n_theta + it
#if defined(OMPGPU)
!$omp target teams distribute private(pvec_in,pvec_out,i,j)
#elif defined(_OPENACC)
!$acc parallel loop gang private(pvec_in,pvec_out,i,j) present(field,xzf)
#else
!$omp parallel do private(pvec_in,pvec_out,i,j)
#endif
  do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_mode == 0) then
#if defined(OMPGPU)
!$omp parallel do
#elif defined(_OPENACC)
!$acc loop vector
#endif
           do i=1,n_theta
             field(1,(ir-1)*n_theta + i,0) = 0.0
           enddo
        else
#if defined(OMPGPU)
!$omp parallel do
#elif defined(_OPENACC)
!$acc loop vector
#endif
           do i=1,n_theta
              pvec_in(i)  = field(1,(ir-1)*n_theta + i,0)
           enddo
#if defined(OMPGPU)
!$omp parallel do private(pvec_out,j)
#elif defined(_OPENACC)
!$acc loop vector private(pvec_out,j)
#endif
           do i=1,n_theta
              pvec_out = 0.0
              do j=1,n_theta
                 pvec_out = pvec_out+xzf(ir,i,j)*pvec_in(j)
              enddo
              field(1,(ir-1)*n_theta + i,0) = pvec_out
           enddo
        endif
  enddo

end subroutine cgyro_field_ae_c


! make a copy from field to field_cur
subroutine cgyro_field_e_sync_cur

   implicit none

#if defined(OMPGPU)
!$omp target update from(field)
#elif defined(_OPENACC)
!$acc update host(field)
#endif
   field_cur = field

end subroutine cgyro_field_e_sync_cur

! compute filed error and save the old values for next round
subroutine cgyro_field_e_compute(delta_t, norm_loc_s,error_loc_s)

  use cgyro_globals, only : nt1, nt2, nc, n_field

  implicit none

  real, intent(in) :: delta_t
  real, intent(out) :: norm_loc_s,error_loc_s

  integer :: itor,ic,i_f

  norm_loc_s = 0.0
  error_loc_s = 0.0

  ! field_olds are always only in system memory... too expensive to keep in GPU memory
  ! assuming field was already synched to system memory
!$omp parallel do collapse(3) reduction(+:norm_loc_s,error_loc_s)
  do itor=nt1,nt2
   do ic=1,nc
     do i_f=1,n_field

        ! 1. Estimate of total (field) error via quadratic interpolation

        field_loc(i_f,ic,itor) = 3*field_old(i_f,ic,itor) - &
                3*field_old2(i_f,ic,itor) + &
                field_old3(i_f,ic,itor)
        field_dot(i_f,ic,itor) = (3*field_cur(i_f,ic,itor) - &
                4*field_old(i_f,ic,itor) + &
                field_old2(i_f,ic,itor) )/(2*delta_t)

        ! Define norm and error for each mode number n
        norm_loc_s  = norm_loc_s  + abs(field_cur(i_f,ic,itor))
        error_loc_s = error_loc_s + abs(field_cur(i_f,ic,itor)-field_loc(i_f,ic,itor))

        ! save old values for next iteration
        field_old3(i_f,ic,itor) = field_old2(i_f,ic,itor)
        field_old2(i_f,ic,itor) = field_old(i_f,ic,itor)
        field_old(i_f,ic,itor)  = field_cur(i_f,ic,itor)
     enddo
   enddo
  enddo

end subroutine cgyro_field_e_compute

! get the top two old elements for field 1, for a single ic and itor
subroutine cgyro_field_e_get_diff1(ic,itor, fo1,fo2)

  implicit none
  integer, intent(in) :: ic,itor
  complex, intent(out) :: fo1,fo2
    
  fo1 = field_old(1,ic,itor)
  fo2 = field_old2(1,ic,itor)
end subroutine cgyro_field_e_get_diff1

subroutine cgyro_field_coefficients

  use mpi
  use cgyro_globals, only : ae_flag, betae_unit, collision_field_model, &
       collision_model, dens2_rot, dens_ele, dens_ele_rot, dens_rot, &
       dvjvec_c, dvjvec_v, i_err, ie_v, ir_c, is_v, it_c, ix_v, &
       jvec_c, jvec_v, k_perp, lambda_debye, nc, nc_cl1, nc_cl2, &
       NEW_COMM_1, n_field, nt1, nt2, nv, nv1, nv2, px, rho, &
       sum_cur_x, sum_den_h, sum_den_x, temp, temp_ele, vfac, w_exi, &
       z, zf_test_mode

  implicit none

  integer :: ic,ic_loc,ir,it,is,ie,ix,iv,iv_loc,itor
  real :: sum_one
  real, dimension(:,:), allocatable :: sum_loc
  real, dimension(:,:), allocatable :: pb11,pb12,pb21,pb22
 
  !-------------------------------------------------------------------------
  ! Field equation prefactors, sums.
  !
  allocate(sum_loc(nc,nt1:nt2))
!$omp parallel do collapse(2) private(iv,iv_loc,is,it,sum_one) shared(sum_loc)
  do itor=nt1,nt2
    do ic=1,nc
      it = it_c(ic)
      sum_one = 0.0
      do iv=nv1,nv2
        iv_loc = iv-nv1+1
        is = is_v(iv)
        sum_one = sum_one + vfac(iv_loc)*dens_rot(it,is) &
             *(1.0-jvec_c(1,ic,iv_loc,itor)**2) 
      enddo
      sum_loc(ic,itor) = sum_one
    enddo
  enddo

  call MPI_ALLREDUCE(sum_loc,&
       sum_den_x,&
       size(sum_den_x),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  if (ae_flag == 1) then
    do itor=nt1,nt2
     do ic=1,nc
        it = it_c(ic)
        sum_den_x(ic,itor) = sum_den_x(ic,itor)+dens_ele*dens_ele_rot(it)/temp_ele
     enddo
    enddo
  endif
  !----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! Field-solve coefficients (i.e., final numerical factors).
  !
  do itor=nt1,nt2
   do ic=1,nc
     ir = ir_c(ic) 
     it = it_c(ic)
     if (itor == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_mode == 0) then
        fcoef(:,ic,itor) = 0.0
     else
        fcoef(1,ic,itor) = 1.0/(k_perp(ic,itor)**2*lambda_debye**2*dens_ele/temp_ele &
             + sum_den_h(it))
        if (n_field > 1) fcoef(2,ic,itor) = 1.0/(-2.0*k_perp(ic,itor)**2* &
             rho**2/betae_unit*dens_ele*temp_ele)
        if (n_field > 2) fcoef(3,ic,itor) = -betae_unit/(2.0*dens_ele*temp_ele)
     endif
   enddo
  enddo

  if (n_field > 1) then

!$omp parallel do collapse(2) private(iv,iv_loc,is,it,sum_one) shared(sum_loc)
     do itor=nt1,nt2
       do ic=1,nc
         it = it_c(ic)
         sum_one = 0.0
         do iv=nv1,nv2
           iv_loc = iv-nv1+1
           is = is_v(iv)
           sum_one = sum_one + vfac(iv_loc)*dens_rot(it,is) &
                *jvec_c(2,ic,iv_loc,itor)**2 
         enddo
         sum_loc(ic,itor) = sum_one
       enddo
     enddo

     call MPI_ALLREDUCE(sum_loc,&
          sum_cur_x,&
          size(sum_cur_x),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

  endif

  if (n_field == 1 .or. n_field == 2) then
    do itor=nt1,nt2
     do ic=1,nc
        if (k_perp(ic,itor) > 0.0) then
           gcoef(1,ic,itor) = 1.0/(k_perp(ic,itor)**2*lambda_debye**2*&
                dens_ele/temp_ele+sum_den_x(ic,itor))
        endif
     enddo
    enddo
  endif

  if (n_field > 1) then
    do itor=nt1,nt2
     do ic=1,nc
        if (k_perp(ic,itor) > 0.0) then
           gcoef(2,ic,itor) = 1.0/(-2.0*k_perp(ic,itor)**2*&
                rho**2/betae_unit*dens_ele*temp_ele-sum_cur_x(ic,itor))
        endif
     enddo
    enddo
  endif

  if (n_field > 2) then
     allocate(pb11(nc,nt1:nt2))
     allocate(pb12(nc,nt1:nt2))
     allocate(pb21(nc,nt1:nt2))
     allocate(pb22(nc,nt1:nt2))

     do itor=nt1,nt2
      do ic=1,nc
        pb11(ic,itor) = k_perp(ic,itor)**2*lambda_debye**2* &
             dens_ele/temp_ele+sum_den_x(ic,itor)
      enddo
     enddo

!$omp parallel do collapse(2) shared(sum_loc) &
!$omp&         private(iv,iv_loc,ir,it,sum_one,is,ix,ie)
     do itor=nt1,nt2
       do ic=1,nc
         ir = ir_c(ic) 
         it = it_c(ic)
         sum_one = 0.0
         do iv=nv1,nv2
           iv_loc = iv-nv1+1
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)
           sum_one = sum_one - w_exi(ie,ix)*dens2_rot(it,is) &
                *z(is)*jvec_c(1,ic,iv_loc,itor)*jvec_c(3,ic,iv_loc,itor) &
                *z(is)/temp(is)
         enddo
         sum_loc(ic,itor) = sum_one
       enddo
     enddo

     call MPI_ALLREDUCE(sum_loc,&
          pb12,&
          size(pb12),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

!$omp parallel do collapse(2) shared(sum_loc) &
!$omp&         private(iv,iv_loc,ir,it,sum_one,is,ix,ie)
     do itor=nt1,nt2
       do ic=1,nc
         ir = ir_c(ic) 
         it = it_c(ic)
         sum_one = 0.0
         do iv=nv1,nv2
           iv_loc = iv-nv1+1
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)
           sum_one = sum_one + w_exi(ie,ix)*dens2_rot(it,is) &
                *temp(is)*jvec_c(3,ic,iv_loc,itor)**2 &
                *(z(is)/temp(is))**2
         enddo
         sum_loc(ic,itor) = sum_one
       enddo
     enddo

     call MPI_ALLREDUCE(sum_loc,&
          pb22,&
          size(pb22),&
          MPI_DOUBLE_PRECISION,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     ! Determinant
!$omp parallel do collapse(2) shared(sum_loc) private(sum_one)
     do itor=nt1,nt2
       do ic=1,nc
         pb21(ic,itor) = pb12(ic,itor)*betae_unit/(-2*dens_ele*temp_ele)
         pb22(ic,itor) = 1.0-pb22(ic,itor)*betae_unit/(-2*dens_ele*temp_ele) 

         if (k_perp(ic,itor) > 0.0) then
           sum_one = pb11(ic,itor)*pb22(ic,itor)-pb12(ic,itor)*pb21(ic,itor)
         else
           sum_one = 1.0
         endif

         gcoef(3,ic,itor) = pb11(ic,itor)/sum_one
         gcoef(1,ic,itor) = pb22(ic,itor)/sum_one
         gcoef(4,ic,itor) = -pb12(ic,itor)/sum_one
         gcoef(5,ic,itor) = -pb21(ic,itor)/sum_one
       enddo
     enddo

     deallocate(pb11)
     deallocate(pb12)
     deallocate(pb21)
     deallocate(pb22)
  endif

  deallocate(sum_loc)

!$omp parallel do private(iv,iv_loc,is,ie,ix,ic,it,ic_loc,ir) shared(gcoef,dvjvec_c,dvjvec_v)
  do itor=nt1,nt2
   ! Set selected zeros
   do ic=1,nc
     ir = ir_c(ic) 
     if (itor == 0 .and. (px(ir) == 0 .or. ir == 1) .and. zf_test_mode == 0) then
        gcoef(:,ic,itor) = 0.0
     endif
   enddo

   ! Arrays to speed up velocity integrals in Maxwell equations
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     ie = ie_v(iv)
     ix = ix_v(iv)
     do ic=1,nc
        it = it_c(ic)
        dvjvec_c(:,ic,iv_loc,itor) = dens2_rot(it,is)*w_exi(ie,ix)*z(is)* &
             jvec_c(:,ic,iv_loc,itor)
     enddo
   enddo
   if ((collision_model /= 5) .AND. (collision_field_model == 1)) then
    do ic=nc_cl1,nc_cl2
     ic_loc = ic-nc_cl1+1
     it = it_c(ic)
     do iv=1,nv
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)
        dvjvec_v(:,iv,itor,ic_loc) = dens2_rot(it,is)*w_exi(ie,ix)*z(is)* &
             jvec_v(:,ic_loc,itor,iv,1) ! all nsm are the same, use the 1st one
     enddo
    enddo
   endif
  enddo
  !-------------------------------------------------------------------------

#if defined(OMPGPU)
!$omp target update to(fcoef,gcoef,dvjvec_c)
#elif defined(_OPENACC)
!$acc update device(fcoef,gcoef,dvjvec_c)
#endif
  if ((collision_model /= 5) .AND. (collision_field_model == 1)) then
#if defined(OMPGPU)
!$omp target update to(dvjvec_v)
#elif defined(_OPENACC)
!$acc update device(dvjvec_v)
#endif
  endif

end subroutine cgyro_field_coefficients

end module cgyro_field_mod
