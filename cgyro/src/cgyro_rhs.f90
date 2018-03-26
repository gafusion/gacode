subroutine cgyro_rhs(ij)

  use timer_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: is
  integer :: id,jc
  real :: rval,rval2
  complex :: rhs_stream
  complex :: rhs_ij(nc,nv_loc)

  call timer_lib_in('str')

  ! Prepare suitable distribution (g, not h) for conservative upwind method
  g_x(:,:) = h_x(:,:)

  if (n_field > 1) then
!$omp parallel do private(iv_loc,is,ic) schedule(dynamic, 1)
     do iv=nv1,nv2
        iv_loc = iv-nv1+1
        is = is_v(iv)
        do ic=1,nc
           g_x(ic,iv_loc) = g_x(ic,iv_loc)+ & 
                (z(is)/temp(is))*jvec_c(2,ic,iv_loc)*field(2,ic)
        enddo
     enddo
  endif

  call timer_lib_out('str')

  ! Correct g_x for number conservation
  call cgyro_upwind

  if ( (nonlinear_flag == 1) .and. (nonlinear_method /= 1) .and. is_staggered_comm_2) then 
     ! stagger comm1, to load ballance network traffic
     call timer_lib_in('nl_comm')
     call cgyro_nl_fftw_comm1
     call timer_lib_out('nl_comm')
  endif

  call timer_lib_in('str')

!$omp parallel do schedule(dynamic,1)&
!$omp& private(iv,ic,iv_loc,is,rval,rval2,rhs_stream,id,jc) 
  do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        ! Diagonal terms
        rhs_ij(ic,iv_loc) = &
             omega_cap_h(ic,iv_loc)*cap_h_c(ic,iv_loc)+&
             omega_h(ic,iv_loc)*h_x(ic,iv_loc)
     enddo 

     ! CPUs vectorize better if separate loop
     do ic=1,nc
        ! Don't need to recompute iv_loc, still valid from previous loop
        is = is_v(iv)
        ! Parallel streaming with upwind dissipation 
        rval  = omega_stream(it_c(ic),is)*vel(ie_v(iv))*xi(ix_v(iv))
        rval2 = abs(omega_stream(it_c(ic),is))

        rhs_stream = 0.0
        do id=-nup_theta,nup_theta
           jc = icd_c(id, ic)
           rhs_stream = rhs_stream &
                -rval*dtheta(id, ic)*cap_h_c(jc,iv_loc)  &
                -rval2*dtheta_up(id, ic)*g_x(jc,iv_loc) 
        enddo

        rhs_ij(ic,iv_loc) = rhs_ij(ic,iv_loc) + &
             sum(omega_s(:,ic,iv_loc)*field(:,ic))

        rhs(ic,iv_loc,ij) = rhs_ij(ic,iv_loc)+rhs_stream

     enddo
  enddo

  call timer_lib_out('str')

  ! Wavenumber advection shear terms
  call cgyro_advect_wavenumber(ij)

  if ( (nonlinear_flag == 1) .and. (nonlinear_method /= 1) .and. (.not. is_staggered_comm_2)) then 
     ! stagger comm1, to load ballance network traffic
     call timer_lib_in('nl_comm')
     call cgyro_nl_fftw_comm1
     call timer_lib_out('nl_comm')
  endif

  ! Nonlinear evaluation [f,g]
  if (nonlinear_flag == 1) then     
     if (nonlinear_method == 1) then
        call cgyro_nl_direct(ij)
     else
        call cgyro_nl_fftw(ij)
     endif
  endif

end subroutine cgyro_rhs
