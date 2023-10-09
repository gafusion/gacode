subroutine cgyro_rhs_comm_async(which)
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  integer, intent(in) :: which

  if (nonlinear_flag == 1) then
     if (which == 1) then
       call cgyro_nl_fftw_comm1_async
     else
       call cgyro_nl_fftw_comm2_async
     endif
  endif

end subroutine cgyro_rhs_comm_async

! Note: Calling test propagates the async operations in some MPI implementations
subroutine cgyro_rhs_comm_test(which)
  use cgyro_nl_comm
  use cgyro_globals

  implicit none

  integer, intent(in) :: which

  if (nonlinear_flag == 1) then
     if (which == 1) then
       call cgyro_nl_fftw_comm1_test
     else
       call cgyro_nl_fftw_comm2_test
     endif
  endif

end subroutine cgyro_rhs_comm_test


subroutine cgyro_rhs(ij)

  use timer_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: is,itor
  integer :: id,jc
  real :: rval,rval2
  complex :: rhs_stream,rhs_el

  ! h_x is not modified after this and before nl_fftw
  call cgyro_rhs_comm_async(1)

  call timer_lib_in('str_mem')

#if (!defined(OMPGPU)) && defined(_OPENACC)
!$acc data present(h_x,g_x,rhs,field)
#endif

  call timer_lib_out('str_mem')


  ! Prepare suitable distribution (g, not h) for conservative upwind method
  if (n_field > 1) then
     call timer_lib_in('str')

#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  private(iv_loc,is)
#elif defined(_OPENACC)
!$acc parallel loop  collapse(3) independent gang vector &
!$acc&         private(iv_loc,is) &
!$acc&         present(is_v,z,temp,jvec_c) &
!$acc&         present(nt1,nt2,nv1,nv2,nc) &
!$acc&         default(none) async(1)
#endif
     do itor=nt1,nt2
      do iv=nv1,nv2
        do ic=1,nc
           iv_loc = iv-nv1+1
           is = is_v(iv)

           g_x(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)+ & 
                (z(is)/temp(is))*jvec_c(2,ic,iv_loc,itor)*field(2,ic,itor)
        enddo
      enddo
     enddo

#if (!defined(OMPGPU)) && defined(_OPENACC)
     call cgyro_rhs_comm_test(1)
!$acc wait(1)
#endif
     call cgyro_rhs_comm_test(1)

     call timer_lib_out('str')
  else
     call timer_lib_in('str_mem')

#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  private(iv_loc)
#elif defined(_OPENACC)
!$acc parallel loop  collapse(3) gang vector &
!$acc&         independent private(iv_loc) &
!$acc&         present(nt1,nt2,nv1,nv2,nc)  &
!$acc&         default(none) async(1)
#endif
     do itor=nt1,nt2
      do iv=nv1,nv2
        do ic=1,nc
           iv_loc = iv-nv1+1
           g_x(ic,iv_loc,itor) = h_x(ic,iv_loc,itor)
        enddo
      enddo
     enddo

#if (!defined(OMPGPU)) && defined(_OPENACC)
     call cgyro_rhs_comm_test(1)
!$acc wait(1)
#endif
     call cgyro_rhs_comm_test(1)

     call timer_lib_out('str_mem')
  endif

  call cgyro_upwind

  call cgyro_rhs_comm_test(1)
  call cgyro_rhs_comm_async(2)

#if (!defined(OMPGPU)) && defined(_OPENACC)
  call timer_lib_in('str_mem')

!$acc data  &
!$acc& present(rhs) &
!$acc& pcopyin(cap_h_c) &
!$acc& present(is_v,ix_v,ie_v,it_c) &
!$acc& present(omega_cap_h,omega_h,omega_s) &
!$acc& present(omega_stream,xi,vel) &
!$acc& present(dtheta,dtheta_up,icd_c)

  call timer_lib_out('str_mem')
#endif
  call timer_lib_in('str')

#if defined(OMPGPU)
  ! no async for OMPGPU for now
!$omp target teams distribute parallel do simd collapse(3) &
!$omp&  private(iv,ic,iv_loc,is,rval,rval2,rhs_stream,id,jc,rhs_el) &
!$omp&  map(to:cap_h_c)
#elif defined(_OPENACC)
!$acc  parallel loop gang vector collapse(3) & 
!$acc& private(iv,ic,iv_loc,is,rval,rval2,rhs_stream,id,jc,rhs_el) async(1)
#endif
  do itor=nt1,nt2
   do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        ! Diagonal terms
        rhs_el = &
             omega_cap_h(ic,iv_loc,itor)*cap_h_c(ic,iv_loc,itor)+&
             omega_h(ic,iv_loc,itor)*h_x(ic,iv_loc,itor)

        is = is_v(iv)
        ! Parallel streaming with upwind dissipation 
        rval  = omega_stream(it_c(ic),is,itor)*vel(ie_v(iv))*xi(ix_v(iv))
        rval2 = abs(omega_stream(it_c(ic),is,itor))

        rhs_stream = 0.0
        do id=-nup_theta,nup_theta
           jc = icd_c(id, ic, itor)
           rhs_stream = rhs_stream &
                -rval*dtheta(id,ic,itor)*cap_h_c(jc,iv_loc,itor)  &
                -rval2*dtheta_up(id,ic,itor)*g_x(jc,iv_loc,itor)
        enddo

        rhs(ic,iv_loc,itor,ij) = rhs_el + rhs_stream +&
             sum(omega_s(:,ic,iv_loc,itor)*field(:,ic,itor))
     enddo
   enddo
  enddo

#if (!defined(OMPGPU)) && defined(_OPENACC)
  call cgyro_rhs_comm_test(1)
!$acc wait(1)
#endif
  call cgyro_rhs_comm_test(1)
  call cgyro_rhs_comm_test(2)

  call timer_lib_out('str')

  if (explicit_trap_flag == 1) then
     ! we should never get in here... should have failed during init
     ! hard abort if we somehow end up here
     call abort
  endif

  ! Wavenumber advection shear terms
  call cgyro_advect_wavenumber(ij)

  ! Nonlinear evaluation [f,g]
  if (nonlinear_flag == 1) then     
     call cgyro_nl_fftw(ij)
  endif

#if (!defined(OMPGPU)) && defined(_OPENACC)
 call timer_lib_in('str_mem')

!$acc end data    

! g_x and h_x
!$acc end data 
 
  call timer_lib_out('str_mem')
#endif

end subroutine cgyro_rhs
