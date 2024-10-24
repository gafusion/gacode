subroutine cgyro_rhs_comm_async_hx
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  if (nonlinear_flag == 1) then
       ! prepare and transfer h_x (first half)
       call cgyro_nl_fftw_comm1_async
  endif

end subroutine cgyro_rhs_comm_async_hx

subroutine cgyro_rhs_comm_async_fd
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  if (nonlinear_flag == 1) then
       ! transfer fields
       call cgyro_nl_fftw_comm2_async
       ! and the other half of h_x
       call cgyro_nl_fftw_comm3_async
  endif

end subroutine cgyro_rhs_comm_async_fd

subroutine cgyro_rhs_r_comm_async(ij)
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij

  if (nonlinear_flag == 1) then
       call cgyro_nl_fftw_comm1_r(ij)
  endif

end subroutine cgyro_rhs_r_comm_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_rhs_comm_test
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  if (nonlinear_flag == 1) then
     call cgyro_nl_fftw_comm_test
  endif

end subroutine cgyro_rhs_comm_test

subroutine cgyro_rhs_comp1(ij,update_cap)

  use timer_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  logical, intent(in) :: update_cap
  !--------------------------------
  integer :: is,itor
  complex :: h_el,cap_el,my_psi,rhs_el

  call timer_lib_in('str')

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc data present(h_x,rhs,field,cap_h_c)

!$acc data  &
!$acc& present(rhs) &
!$acc& present(cap_h_c,z,temp,jvec_c) &
!$acc& present(is_v,ix_v,ie_v,it_c) &
!$acc& present(omega_cap_h,omega_h,omega_s)

#endif

  ! get the initial rhs initialization
  ! that depends on h_x,cap_h_c and field only
#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  private(iv,ic,iv_loc,rhs_el,h_el,is,cap_el,my_psi)
#elif defined(_OPENACC)
!$acc  parallel loop gang vector collapse(3) & 
!$acc& private(iv,ic,iv_loc,rhs_el,h_el,is,cap_el,my_psi) async(1)
#else
!$omp parallel do collapse(2) &
!$omp& private(iv,iv_loc,itor,is,ic,rhs_el,h_el,cap_el,my_psi) 
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        h_el = h_x(ic,iv_loc,itor)
        if (update_cap) then
           is = is_v(iv)
           my_psi = sum( jvec_c(:,ic,iv_loc,itor)*field(:,ic,itor))
           cap_el = h_el+my_psi*z(is)/temp(is)
           cap_h_c(ic,iv_loc,itor) = cap_el
        else
           cap_el = cap_h_c(ic,iv_loc,itor)
        endif
        ! Diagonal terms
        rhs_el = &
             omega_cap_h(ic,iv_loc,itor)*cap_el+&
             omega_h(ic,iv_loc,itor)*h_el

        rhs(ic,iv_loc,itor,ij) = rhs_el + &
             sum(omega_s(:,ic,iv_loc,itor)*field(:,ic,itor))
     enddo
   enddo
  enddo

  call cgyro_rhs_comm_test

#if (!defined(OMPGPU)) && defined(_OPENACC)
  !no async for OMPGPU for now
!$acc wait(1)

!$acc end data    

! h_x
!$acc end data 
#endif

  call timer_lib_out('str')

end subroutine cgyro_rhs_comp1

subroutine cgyro_rhs_comp2(ij)

  use timer_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  !--------------------------------
  integer :: is,itor
  integer :: id,jc
  real :: rval,rval2
  complex :: rhs_stream


  call timer_lib_in('str')

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc data present(g_x,rhs,field,cap_h_c)

!$acc data  &
!$acc& present(rhs) &
!$acc& present(cap_h_c) &
!$acc& present(is_v,ix_v,ie_v,it_c) &
!$acc& present(omega_stream,xi,vel) &
!$acc& present(dtheta,dtheta_up,icd_c)

#endif
  ! add stream to rhs
#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  private(iv,ic,iv_loc,is,rval,rval2,rhs_stream,id,jc)
#elif defined(_OPENACC)
!$acc  parallel loop gang vector collapse(3) & 
!$acc& private(iv,ic,iv_loc,is,rval,rval2,rhs_stream,id,jc) async(1)
#else
!$omp parallel do collapse(2) &
!$omp& private(itor,iv,iv_loc,is,ic,rval,rval2,rhs_stream,id,jc) 
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        is = is_v(iv)
        ! Parallel streaming with upwind dissipation 
        rval  = omega_stream(it_c(ic),is,itor)*vel(ie_v(iv))*xi(ix_v(iv))
        rval2 = abs(omega_stream(it_c(ic),is,itor))

        rhs_stream = 0.0
        do id=-nup_theta,nup_theta
           jc = icd_c(ic, id, itor)
           rhs_stream = rhs_stream &
                -rval*dtheta(ic,id,itor)*cap_h_c(jc,iv_loc,itor)  &
                -rval2*dtheta_up(ic,id,itor)*g_x(jc,iv_loc,itor)
        enddo

        rhs(ic,iv_loc,itor,ij) = rhs(ic,iv_loc,itor,ij) + rhs_stream
     enddo
   enddo
  enddo

  call cgyro_rhs_comm_test
#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc wait(1)

!$acc end data    

! g_x
!$acc end data 
#endif

  call timer_lib_out('str')

end subroutine cgyro_rhs_comp2

#if defined(OMPGPU) || defined(_OPENACC)
! gpu code

subroutine cgyro_rhs_trap(ij)

  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  !--------------------------------

  if (explicit_trap_flag == 1) then
     ! we should never get in here... should have failed during init
     ! hard abort if we somehow end up here
     call abort
  endif

end subroutine cgyro_rhs_trap

#else
! cpu code

subroutine cgyro_rhs_trap(ij)

  use timer_lib
  use parallel_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  !--------------------------------
  integer :: is, ix, ie, js, jx, je, jv, it, j, k, itor
  complex, dimension(:,:), allocatable :: rhs_trap
  complex, dimension(:), allocatable   :: bvec_trap
  integer :: nj_loc

  ! Explicit trapping term
  if (explicit_trap_flag == 1) then
     call timer_lib_in('str')

     allocate(rhs_trap(nc,nv_loc))
     allocate(bvec_trap(nv))
     call parallel_lib_rtrans_pack(cap_h_c)
     call parallel_lib_r_do(cap_h_v)
     call parallel_lib_nj_loc(nj_loc)
     
     do itor=nt1,nt2
      do ic=nc1,nc2
        ic_loc = ic-nc1+1
        it = it_c(ic)
        
        bvec_trap(:) = (0.0,0.0)
        do iv=1,nv
           is = is_v(iv)
           ix = ix_v(iv)
           ie = ie_v(iv)
           do jv=1,nv
              js = is_v(jv)
              jx = ix_v(jv)
              je = ie_v(jv)
              if (is == js) then
                 if (ie == je) then
                    bvec_trap(iv) = bvec_trap(iv) - (omega_trap(it,is,itor) * vel(ie) &
                         + omega_rot_trap(it,is) / vel(ie)) &
                         * (1.0 - xi(ix)**2) * xi_deriv_mat(ix,jx) * cap_h_v(ic_loc,itor,jv)
                 endif
                 if (ix == jx) then
                    bvec_trap(iv) = bvec_trap(iv) - omega_rot_u(it,is) * xi(ix) &
                         * e_deriv1_rot_mat(ie,je)/sqrt(1.0*e_max) * cap_h_v(ic_loc,itor,jv)
                 endif
              endif
           enddo
        enddo

        do k=1,nproc
           do j=1,nj_loc
              fsendf(j,itor,ic_loc,k) = bvec_trap(j+(k-1)*nj_loc)
           enddo
        enddo
      enddo
     enddo
     
     call parallel_lib_f_i_do(cap_h_ct)
     do itor=nt1,nt2
      do iv=nv1,nv2
        iv_loc = iv-nv1+1
        do ic=1,nc
           rhs_trap(ic,iv_loc) = cap_h_ct(iv_loc,itor,ic)
        enddo
      enddo

      rhs(:,:,itor,ij) = rhs(:,:,itor,ij) +  rhs_trap(:,:)
     enddo
     
     deallocate(rhs_trap)
     deallocate(bvec_trap)
     
     call timer_lib_out('str')
  endif
  

end subroutine cgyro_rhs_trap

#endif

subroutine cgyro_rhs(ij,update_cap)

  use timer_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  logical, intent(in) :: update_cap
  !--------------------------------

  ! fields is ready by now
  call cgyro_rhs_comm_async_fd

  call cgyro_upwind_prepare_async

  call cgyro_rhs_comp1(ij,update_cap)

  call cgyro_rhs_comm_test

  ! Wavenumber advection shear terms
  ! depends on h_x and field only, updates rhs
  call cgyro_advect_wavenumber(ij)
  !call cgyro_globalshear(ij)

  call cgyro_rhs_comm_test

  ! Nonlinear evaluation [f,g]
  if (nonlinear_flag == 1) then
     ! assumes someone already started the input comm
     ! and will finish the output comm
     call cgyro_nl_fftw()
  endif

  call cgyro_upwind_complete

  call cgyro_rhs_comm_test

  call cgyro_rhs_comp2(ij)

  call cgyro_rhs_trap(ij)

  ! updates rhs
  call cgyro_rhs_r_comm_async(ij)

end subroutine cgyro_rhs

