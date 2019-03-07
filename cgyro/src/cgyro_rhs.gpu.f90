subroutine cgyro_rhs(ij)

  use timer_lib
  use cgyro_globals

  implicit none

  integer, intent(in) :: ij
  integer :: is
  integer :: id,jc
  real :: rval,rval2
  complex :: rhs_stream

  call timer_lib_in('str_mem')

!$acc data present(h_x,g_x,psi,rhs,field)

  call timer_lib_out('str_mem')


  ! Prepare suitable distribution (g, not h) for conservative upwind method
  if (n_field > 1) then
     call timer_lib_in('str')

!$acc parallel loop  collapse(2) independent private(iv_loc,is) &
!$acc&         present(is_v,z,temp,jvec_c) default(none)
     do iv=nv1,nv2
        do ic=1,nc
           iv_loc = iv-nv1+1
           is = is_v(iv)

           g_x(ic,iv_loc) = h_x(ic,iv_loc)+ & 
                (z(is)/temp(is))*jvec_c(2,ic,iv_loc)*field(2,ic)
        enddo
     enddo

     call timer_lib_out('str')
  else
     call timer_lib_in('str_mem')

!$acc parallel loop  collapse(2) independent private(iv_loc) default(none)
     do iv=nv1,nv2
        do ic=1,nc
           iv_loc = iv-nv1+1
           g_x(ic,iv_loc) = h_x(ic,iv_loc)
        enddo
     enddo

     call timer_lib_out('str_mem')
  endif

  call cgyro_upwind

  if ( (nonlinear_flag == 1) .and. (nonlinear_method /= 1) .and. is_staggered_comm_2) then 
     ! stagger comm1, to load ballance network traffic
     call cgyro_nl_fftw_comm1
  endif

  call timer_lib_in('str_mem')

!$acc data  &
!$acc& present(rhs) &
!$acc& pcopyin(cap_h_c) &
!$acc& present(is_v,ix_v,ie_v,it_c) &
!$acc& present(omega_cap_h,omega_h,omega_s) &
!$acc& present(omega_stream,xi,vel) &
!$acc& present(dtheta,dtheta_up,icd_c)

  call timer_lib_out('str_mem')
  call timer_lib_in('str')

!$acc  parallel loop gang vector collapse(2) & 
!$acc& private(iv,ic,iv_loc,is,rval,rval2,rhs_stream,id,jc)
  do iv=nv1,nv2
     do ic=1,nc
        iv_loc = iv-nv1+1
        ! Diagonal terms
        rhs(ic,iv_loc,ij) = &
             omega_cap_h(ic,iv_loc)*cap_h_c(ic,iv_loc)+&
             omega_h(ic,iv_loc)*h_x(ic,iv_loc)

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

        rhs(ic,iv_loc,ij) = rhs(ic,iv_loc,ij) + rhs_stream +&
             sum(omega_s(:,ic,iv_loc)*field(:,ic))
     enddo
  enddo

  call timer_lib_out('str')

  ! Wavenumber advection shear terms
  call cgyro_advect_wavenumber(ij)

  if ( (nonlinear_flag == 1) .and. (nonlinear_method /= 1) .and. (.not. is_staggered_comm_2)) then 
     ! stagger comm1, to load ballance network traffic
     call cgyro_nl_fftw_comm1
  endif

  ! Nonlinear evaluation [f,g]
  if (nonlinear_flag == 1) then     
     if (nonlinear_method == 1) then
        call cgyro_nl_direct(ij)
     else
        call cgyro_nl_fftw(ij)
     endif
  endif

 call timer_lib_in('str_mem')

!$acc end data    

  ! g_x and h_x
!$acc end data 
  call timer_lib_out('str_mem')

end subroutine cgyro_rhs
