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

!-----------------------------------------------------------------
! Velocity (configuration-distributed) field solve
!-----------------------------------------------------------------

! like cgyro_field_v_notae, but with parametrized start_t
subroutine cgyro_field_v_notae_s(start_t)

  use parallel_lib
  use timer_lib
  use cgyro_globals

  implicit none
  ! ------------------ 
  integer, intent(in) :: start_t
  !
  integer :: itor,j,k,nj_loc

  call timer_lib_in('field')

  field_loc(:,:,start_t:nt2) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of H

  ! iv = j+(k-1)*nj_loc
  ! cap_h_v(ic_loc,itor,iv) = fsendf(j,itor,ic_loc,k)
  call parallel_lib_nj_loc(nj_loc)

!$omp parallel do collapse(2) private(ic_loc,iv,ic,k,j)
  do itor=start_t,nt2
   do ic=nc1,nc2
     ic_loc = ic-nc1+1
     do k=1,nproc
      do j=1,nj_loc
        iv = j+(k-1)*nj_loc
        field_loc(:,ic,itor) = field_loc(:,ic,itor)+dvjvec_v(:,ic_loc,itor,iv)*fsendf(j,itor,ic_loc,k)
      enddo
     enddo
   enddo
  enddo

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  call parallel_lib_sum_field(field_loc(:,:,start_t:nt2), &
                              field(:,:,start_t:nt2))

  call timer_lib_out('field_com')
  
  call timer_lib_in('field')

  ! Poisson LHS factors
!$omp parallel do
  do itor=start_t,nt2
     ! assuming  (.not.(itor == 0 .and. ae_flag == 1))
     field(:,:,itor) = fcoef(:,:,itor)*field(:,:,itor)
  enddo

  call timer_lib_out('field')
end subroutine cgyro_field_v_notae_s

! like cgyro_field_v, but skip (itor == 0 .and. ae_flag == 1)
subroutine cgyro_field_v_notae

  use parallel_lib
  use timer_lib
  use cgyro_globals

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

!
! GPU versions
!

! Note: Not supporting the ae-version of cgyro_field_v_gpu

subroutine cgyro_field_v_notae_s_gpu(start_t)
  use parallel_lib
  use timer_lib
  use cgyro_globals

  implicit none
  ! ------------------ 
  integer, intent(in) :: start_t
  !
  integer :: i_f,itor,j,k,nj_loc
  complex :: field_loc_l 

  call timer_lib_in('field')

  ! iv = j+(k-1)*nj_loc
  ! cap_h_v(ic_loc,itor,iv) = fsendf(j,itor,ic_loc,k)
  call parallel_lib_nj_loc(nj_loc)

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc data present(fsendf)
!$acc data present(field,field_loc)
#endif

  ! Poisson and Ampere RHS integrals of H

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   map(to:start_t)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector &
!$acc&         independent copyin(start_t) &
!$acc&         present(nt2,nc,n_field) default(none)
#endif
  do itor=start_t,nt2
   do ic=1,nc
       do i_f=1,n_field
        field_loc(i_f,ic,itor) = (0.0,0.0)
       enddo
   enddo
  enddo

#if defined(OMPGPU)
!$omp target teams distribute collapse(3) &
!$omp&       private(ic_loc,field_loc_l) map(to:start_t,nproc,nj_loc)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang private(ic_loc,field_loc_l) &
!$acc&         present(dvjvec_v,fsendf,field_loc) copyin(start_t,nproc,nj_loc) &
!$acc&         present(nt2,nc1,nc2,n_field,nv) default(none)
#endif
  do itor=start_t,nt2
   do ic=nc1,nc2
    do i_f=1,n_field
      ic_loc = ic-nc1+1
      field_loc_l = (0.0,0.0)
#if defined(OMPGPU)
!$omp parallel do simd collapse(2) reduction(+:field_loc_l) &
!$omp&       private(iv)
#elif defined(_OPENACC)
!$acc loop vector collapse(2) private(iv) reduction(+:field_loc_l)
#endif
      do k=1,nproc
       do j=1,nj_loc
        iv = j+(k-1)*nj_loc
        field_loc_l = field_loc_l+dvjvec_v(i_f,ic_loc,itor,iv)*fsendf(j,itor,ic_loc,k)
      enddo
     enddo
     field_loc(i_f,ic,itor) = field_loc_l
    enddo
   enddo
  enddo

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  call parallel_lib_sum_field_gpu(field_loc(:,:,start_t:nt2), &
                                  field(:,:,start_t:nt2))

  call timer_lib_out('field_com')

  call timer_lib_in('field')
  ! Poisson LHS factors
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&    map(to:start_t)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector &
!$acc&         independent present(fcoef) copyin(start_t) &
!$acc&         present(nt2,nc,n_field) default(none)
#endif
  do itor=start_t,nt2
     ! assuming  (.not.(itor == 0 .and. ae_flag == 1))
     do ic=1,nc
       do i_f=1,n_field
        field(i_f,ic,itor) = fcoef(i_f,ic,itor)*field(i_f,ic,itor)
       enddo
     enddo
  enddo

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc end data
!$acc end data
#endif

  call timer_lib_out('field')

end subroutine cgyro_field_v_notae_s_gpu

! like cgyro_field_v, but skip (itor == 0 .and. ae_flag == 1)
subroutine cgyro_field_v_notae_gpu

  use parallel_lib
  use timer_lib
  use cgyro_globals

  implicit none

  if (nt1 == 0 .and. ae_flag == 1) then
     if (nt2>0) then
        call cgyro_field_v_notae_s_gpu(1)
     endif
     ! else no-op
  else
     ! don't have to worry about ae_flag, just use all the elements
     call cgyro_field_v_notae_s_gpu(nt1)
  endif

end subroutine cgyro_field_v_notae_gpu

!-----------------------------------------------------------------
! Configuration (velocity-distributed) field solve
!-----------------------------------------------------------------
subroutine cgyro_field_c_cpu

  use parallel_lib
  use timer_lib
  use cgyro_globals

  implicit none

  integer :: is,itor
  complex :: my_psi
  
  complex, dimension(nc) :: tmp
  
  call timer_lib_in('field')

  field_loc(:,:,:) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of h

!$omp parallel private(iv_loc,ic)
!$omp do collapse(2) reduction(+:field_loc)
  do itor=nt1,nt2
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do ic=1,nc
        field_loc(:,ic,itor) = field_loc(:,ic,itor)+dvjvec_c(:,ic,iv_loc,itor)*h_x(ic,iv_loc,itor)
     enddo
   enddo
  enddo
!$omp end do
!$omp end parallel

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  call parallel_lib_sum_field(field_loc,field)

  call timer_lib_out('field_com')

  call timer_lib_in('field')

!$omp parallel do private(tmp) shared(field)
  do itor=nt1,nt2
   if (n_field > 2) then
     field(3,:,itor) = field(3,:,itor)*fcoef(3,:,itor)
   endif

   ! Poisson LHS factors
   if (itor == 0 .and. ae_flag == 1) then
    call cgyro_field_ae('c')
   else
     if (n_field > 2) then
        tmp(:) = field(1,:,itor)
        field(1,:,itor) = gcoef(1,:,itor)*field(1,:,itor) + &
                gcoef(4,:,itor)*field(3,:,itor)
        field(2,:,itor) = gcoef(2,:,itor)*field(2,:,itor)
        field(3,:,itor) = gcoef(3,:,itor)*field(3,:,itor) + &
                gcoef(5,:,itor)*tmp(:)
     else
        field(:,:,itor) = gcoef(:,:,itor)*field(:,:,itor)
     endif
   endif
  enddo

!$omp parallel do collapse(2) private(iv_loc,is,ic,my_psi)
  do itor=nt1,nt2
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
        cap_h_c(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+my_psi*z(is)/temp(is)
     enddo
   enddo
  enddo

  call timer_lib_out('field')

end subroutine cgyro_field_c_cpu

! like cgyro_field_c, but assume (my_toroidal == 0 .and. ae_flag == 1)
subroutine cgyro_field_c_ae_cpu

  use parallel_lib
  use timer_lib
  use cgyro_globals

  implicit none

  integer :: is,itor
  complex :: my_psi
  
  call timer_lib_in('field')

  field_loc(:,:,0:0) = (0.0,0.0)

  ! Poisson and Ampere RHS integrals of h

!$omp parallel private(iv_loc,ic)
!$omp do collapse(2) reduction(+:field_loc)
  do itor=0,0
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     do ic=1,nc
        field_loc(:,ic,itor) = field_loc(:,ic,itor)+dvjvec_c(:,ic,iv_loc,itor)*h_x(ic,iv_loc,itor)
     enddo
   enddo
  enddo
!$omp end do
!$omp end parallel

  call timer_lib_out('field')

  call timer_lib_in('field_com')

  call parallel_lib_sum_field(field_loc(:,:,0:0), &
                              field(:,:,0:0))

  call timer_lib_out('field_com')

  call timer_lib_in('field')

  do itor=0,0
   if (n_field > 2) then
     field(3,:,itor) = field(3,:,itor)*fcoef(3,:,itor)
   endif

   ! Poisson LHS factors
   call cgyro_field_ae('c')
  enddo

!$omp parallel do collapse(2) private(iv_loc,is,ic,my_psi)
  do itor=0,0
   do iv=nv1,nv2
     iv_loc = iv-nv1+1
     is = is_v(iv)
     do ic=1,nc
        my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
        cap_h_c(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+my_psi*z(is)/temp(is)
     enddo
   enddo
  enddo

  call timer_lib_out('field')

end subroutine cgyro_field_c_ae_cpu

#if defined(OMPGPU) || defined(_OPENACC)
subroutine cgyro_field_c_gpu
  use parallel_lib
  use timer_lib
  use cgyro_globals
  implicit none
  integer :: is,i_f,itor
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

  call parallel_lib_sum_field_gpu(field_loc,field)

  call timer_lib_out('field_com')
  call timer_lib_in('field')
  if (n_field > 2) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(2)
#elif defined(_OPENACC)
!$acc parallel loop collapse(2) gang vector &
!$acc&         independent present(fcoef) &
!$acc&         present(nt1,nt2,nc) default(none)
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
    ! Note: Called rarely, use the CPu version
#if defined(OMPGPU)
!$omp target update from(field)
#elif defined(_OPENACC)
!$acc update host(field)
#endif
    call cgyro_field_ae('c')
#if defined(OMPGPU)
!$omp target update to(field)
#elif defined(_OPENACC)
!$acc update device(field)
#endif
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

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   private(iv_loc,is,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector private(iv_loc,is,my_psi) &
!$acc&         present(jvec_c,z,temp,is_v) present(nt1,nt2,nv1,nv2,nc) default(none)
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

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc end data
!$acc end data
#endif

  call timer_lib_out('field')
end subroutine cgyro_field_c_gpu

! like cgyro_field_c, but assume (my_toroidal == 0 .and. ae_flag == 1)
subroutine cgyro_field_c_ae_gpu
  use parallel_lib
  use timer_lib
  use cgyro_globals
  implicit none
  integer :: is,i_f,itor
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

  call parallel_lib_sum_field_gpu(field_loc(:,:,0:0), &
                                  field(:,:,0:0))

  call timer_lib_out('field_com')
  call timer_lib_in('field')
  if (n_field > 2) then
#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(2)
#elif defined(_OPENACC)
!$acc parallel loop collapse(2) gang vector &
!$acc&         independent present(fcoef) present(nc) default(none)
#endif
    do itor=0,0
      do ic=1,nc
       field(3,ic,itor) = field(3,ic,itor)*fcoef(3,ic,itor)
      enddo
     enddo
  endif
  ! Poisson LHS factors
  ! Note: Called rarely, use the CPU version
#if defined(OMPGPU)
!$omp target update from(field(:,:,0:0))
#elif defined(_OPENACC)
!$acc update host(field(:,:,0:0))
#endif
    call cgyro_field_ae('c')
#if defined(OMPGPU)
!$omp target update to(field(:,:,0:0))
#elif defined(_OPENACC)
!$acc update device(field(:,:,0:0))
#endif

#if defined(OMPGPU)
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&   private(iv_loc,is,my_psi)
#elif defined(_OPENACC)
!$acc parallel loop collapse(3) gang vector private(iv_loc,is,my_psi) &
!$acc&         present(jvec_c,z,temp,is_v) present(nv1,nv2,nc) default(none)
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
end subroutine cgyro_field_c_ae_gpu

#endif


subroutine cgyro_field_c
  implicit none
#if defined(OMPGPU) || defined(_OPENACC)
   call cgyro_field_c_gpu
#else
   call cgyro_field_c_cpu
#endif
end subroutine cgyro_field_c

! like cgyro_field_c, but only for (itor == 0 .and. ae_flag == 1)
subroutine cgyro_field_c_ae
  implicit none
#if defined(OMPGPU) || defined(_OPENACC)
   call cgyro_field_c_ae_gpu
#else
   call cgyro_field_c_ae_cpu
#endif
end subroutine cgyro_field_c_ae

!-----------------------------------------------------------------
! Adiabatic electron field solves for n=0
! Can only be called if itor==0
!-----------------------------------------------------------------
subroutine cgyro_field_ae(space)

  use cgyro_globals
  use timer_lib  

  implicit none

  character(len=1), intent(in) :: space
  integer :: ir,i,j
  complex, dimension(n_theta) :: pvec_in,pvec_out

  if (space == 'c') then
     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_mode == 0) then
           field(1,ic_c(ir,:),0) = 0.0
        else
           do i=1,n_theta
              pvec_out(i) = 0.0
              pvec_in(i)  = field(1,ic_c(ir,i),0)
           enddo
           do j=1,n_theta
              do i=1,n_theta
                 pvec_out(i) = pvec_out(i)+xzf(ir,i,j)*pvec_in(j)
              enddo
           enddo
           do i=1,n_theta
              field(1,ic_c(ir,i),0) = pvec_out(i)
           enddo
        endif
     enddo
  else
     do ir=1,n_radial
        if ((px(ir) == 0 .or. ir == 1) .and. zf_test_mode == 0) then
           field(1,ic_c(ir,:),0) = 0.0
        else
           do i=1,n_theta
              pvec_out(i) = 0.0
              pvec_in(i)  =  field(1,ic_c(ir,i),0)
           enddo
           do j=1,n_theta
              do i=1,n_theta
                 pvec_out(i) = pvec_out(i)+hzf(ir,i,j)*pvec_in(j)
              enddo
           enddo
           do i=1,n_theta
              field(1,ic_c(ir,i),0) = pvec_out(i)
           enddo
        endif
     enddo
  endif

end subroutine cgyro_field_ae
