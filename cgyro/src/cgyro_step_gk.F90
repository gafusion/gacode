subroutine cgyro_step_gk

  use cgyro_globals

  implicit none

  ! RK4 time-advance for the distribution 
  !
  !           z e             vpar            z e  vperp^2
  !  h = H - ----- G0 ( phi - ----- Apar ) + ----- ---------- Gperp Bpar
  !            T               c               T   omega_a c
  !
  ! After time advance, we will have 
  !
  ! h    -> h_x
  ! H    -> cap_h_c
  ! phi  -> field(1)
  ! Apar -> field(2)
  ! Bpar -> field(3)

  h0_x = h_x
  
  ! Stage 1
  call cgyro_rhs(1)
  h_x = h0_x + 0.5 * delta_t * rhs(:,:,1)
  call cgyro_field_c

  ! Stage 2
  call cgyro_rhs(2)
  h_x = h0_x + 0.5 * delta_t * rhs(:,:,2)
  call cgyro_field_c

  ! Stage 3
  call cgyro_rhs(3)
  h_x = h0_x + delta_t * rhs(:,:,3)
  call cgyro_field_c

  ! Stage 4
  call cgyro_rhs(4)
  h_x = h0_x+delta_t*(rhs(:,:,1)+2*rhs(:,:,2)+2*rhs(:,:,3)+rhs(:,:,4))/6  
  call cgyro_field_c

  ! rhs(1) = 3rd-order error estimate
  rhs(:,:,1) = h0_x+delta_t*(rhs(:,:,2)+2*rhs(:,:,3))/3-h_x
  
  ! Filter special spectral components
  call cgyro_filter
  
end subroutine cgyro_step_gk
  
!==========================================================================

subroutine cgyro_rhs(ij)

  use timer_lib

  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: is,ir
  integer :: id,jc
  real :: rval,rval2
  complex :: rhs_stream
  complex :: rhs_ij(nc,nv_loc)

  ! Prepare suitable distribution (g, not h) for conservative upwind method
  g_x(:,:) = h_x(:,:)

  if (n_field > 1) then
!$omp parallel do private(iv_loc,is,ic)
     do iv=nv1,nv2
        iv_loc = iv-nv1+1
        is = is_v(iv)
        do ic=1,nc
           g_x(ic,iv_loc) = g_x(ic,iv_loc)+ & 
                (z(is)/temp(is))*jvec_c(2,ic,iv_loc)*field(2,ic)
        enddo
     enddo
  endif

  call timer_lib_in('str_comm')
  call cgyro_upwind
  call timer_lib_out('str_comm')

  call timer_lib_in('str')

#ifdef _OPENACC
!$acc data  &
!$acc& pcopyout(rhs_ij) &
!$acc& pcopyin(g_x,h_x,field,cap_h_c) &
!$acc& present(is_v,ix_v,ie_v,it_c) &
!$acc& present(omega_cap_h,omega_h,omega_s) &
!$acc& present(omega_stream,xi,vel) &
!$acc& present(dtheta,dtheta_up,icd_c)
  
!$acc  parallel loop gang vector collapse(2) & 
!$acc& private(iv,ic,iv_loc,is,rval,rval2,rhs_stream,id,jc)
#else
!$omp parallel do collapse(2) &
!$omp& private(iv,ic,iv_loc,is,rval,rval2,rhs_stream,id,jc) 
#endif
     do iv=nv1,nv2
        do ic=1,nc
           iv_loc = iv-nv1+1
           is = is_v(iv)
           ! Diagonal terms
           rhs_ij(ic,iv_loc) = &
                omega_cap_h(ic,iv_loc)*cap_h_c(ic,iv_loc)+&
                omega_h(ic,iv_loc)*h_x(ic,iv_loc)+&
                sum(omega_s(:,ic,iv_loc)*field(:,ic))

           ! Parallel streaming with upwind dissipation 
           rval  = omega_stream(it_c(ic),is)*vel(ie_v(iv))*xi(ix_v(iv))
           rval2 = abs(omega_stream(it_c(ic),is))

           rhs_stream = 0.0
           do id=-nup_theta,nup_theta
              jc = icd_c(ic,id)
              rhs_stream = rhs_stream &
                   -rval*dtheta(ic,id)*cap_h_c(jc,iv_loc)  &
                   -rval2*dtheta_up(ic,id)*g_x(jc,iv_loc) 
           enddo

           rhs_ij(ic,iv_loc) = rhs_ij(ic,iv_loc)+rhs_stream

        enddo
     enddo
!$acc end data	  

  rhs(:,:,ij) = rhs_ij(:,:)

  call timer_lib_out('str')

  ! Wavenumber advection shear terms
  call cgyro_advect_wavenumber(ij)

  ! Nonlinear evaluation [f,g]
  if (nonlinear_flag == 1) then     
     if (nonlinear_method == 1) then
        call cgyro_nl_direct(ij)
     else
        call cgyro_nl_fftw(ij)
     endif
  endif

  ! Remove p=-M
  if (psym_flag == 1) then
     do ic=1,nc
        ir = ir_c(ic) 
        if (ir == 1) then
           rhs(ic,:,ij) = 0.0
        endif
     enddo
  endif

end subroutine cgyro_rhs

!==========================================================================

subroutine cgyro_filter

  use cgyro_globals

  implicit none

  integer :: ir
  
  if (n == 0 .and. zf_test_flag == 0) then
     do ic=1,nc
        ir = ir_c(ic) 
        if (ir == 1 .or. px(ir) == 0) then
           h_x(ic,:)     = 0.0
           cap_h_c(ic,:) = 0.0
        endif
     enddo
  endif

end subroutine cgyro_filter

!==========================================================================
! NOTE: The routine below is NOT used.

subroutine cgyro_rhs_trap(ij)

  use parallel_lib

  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: is,ir,it,ie,ix,jx
  complex :: val

  call parallel_lib_r(transpose(cap_h_c),cap_h_v)
  cap_h_v_prime(:,:) = (0.0,0.0)
  ic_loc = 0
  do ic=nc1,nc2
     ic_loc = ic_loc+1
     it = it_c(ic)
     ir = ir_c(ic)
     do iv=1,nv
        is = is_v(iv)
        ix = ix_v(iv)
        ie = ie_v(iv)

        do jx=1,n_xi
           cap_h_v_prime(ic_loc,iv) = cap_h_v_prime(ic_loc,iv) &
                +xi_deriv_mat(ix,jx)*cap_h_v(ic_loc,iv_v(ie,jx,is))
        enddo
     enddo
  enddo

  ! Now have cap_h_v(ic_loc,iv)   

  call parallel_lib_f(cap_h_v_prime,cap_h_ct)
  cap_h_c = transpose(cap_h_ct)

  iv_loc = 0
  do iv=nv1,nv2

     iv_loc = iv_loc+1

     is = is_v(iv)
     ix = ix_v(iv)
     ie = ie_v(iv)

     do ic=1,nc

        ir = ir_c(ic) 
        it = it_c(ic)

        val = omega_trap(it,is)*sqrt(energy(ie))*(1.0-xi(ix)**2) 

        rhs(ic,iv_loc,ij) = rhs(ic,iv_loc,ij)-val*cap_h_c(ic,iv_loc)
     enddo
  enddo

end subroutine cgyro_rhs_trap
