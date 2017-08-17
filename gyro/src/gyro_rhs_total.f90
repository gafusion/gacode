
!-----------------------------------------------------------
! gyro_rhs_total.f90
!
! PURPOSE:
!  Compute the RHS of the electron and ion GKEs for both 
!  periodic and nonperiodic boundary conditions.  This is 
!  a dramatically simplified version of the original code 
!  to compute the righthand side.
!-----------------------------------------------------------

subroutine gyro_rhs_total

  use gyro_globals
  use gyro_pointers
  use math_constants
  use ompdata

  !-----------------------------------------------------------------------------
  implicit none
  !
  complex :: z_der   
  complex :: z_dis
  complex, dimension(:,:,:,:), allocatable :: cap_h
  complex, dimension(:,:,:,:), allocatable :: lit_h
  !-----------------------------------------------------------------------------

  allocate(cap_h(n_stack,i1_buffer:i2_buffer,n_nek_loc_1,n_kinetic))
  allocate(lit_h(n_stack,i1_buffer:i2_buffer,n_nek_loc_1,n_kinetic))

!$omp parallel
  rhs(:,ibeg:iend,:,:)    = (0.0,0.0)
  rhs_dr(:,ibeg:iend,:,:) = 0.0
!$omp end parallel

  !---------------------------------------------
  if (n_substep == 0) then
     call gyro_rhs_nonlinear
  endif
  !---------------------------------------------

  call gyro_timer_in('RHS-total')

  !----------------------------------------------------------------------
  ! Compute orbit-time derivative for ions ONLY, because this term is 
  ! explicit for ions.
  !
  if (kill_i_parallel_flag == 0) call gyro_tau_derivative
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Derivative-friendly functions which have array elements outside
  ! the region 1:n_x.                       
  !
  cap_h(:,i1_buffer:0,:,:) = (0.0,0.0)
  lit_h(:,i1_buffer:0,:,:) = (0.0,0.0)
  cap_h(:,n_x+1:i2_buffer,:,:) = (0.0,0.0)
  lit_h(:,n_x+1:i2_buffer,:,:) = (0.0,0.0)

!$omp parallel private(z_der,z_dis,is,p_nek_loc,i,m,i_diff)
  do is=1,n_kinetic
     do p_nek_loc=1,n_nek_loc_1
        do i = ibeg, iend
           cap_h(:,i,p_nek_loc,is) = (0.0,0.0)
           lit_h(:,i,p_nek_loc,is) = (0.0,0.0)
           do m=1,n_stack
              lit_h(m,i,p_nek_loc,is) = h(m,i,p_nek_loc,is)
              cap_h(m,i,p_nek_loc,is) = h(m,i,p_nek_loc,is)+&
                   z(is)*alpha_s(is,i)*gyro_u(m,i,p_nek_loc,is)
           enddo
        enddo
!$omp barrier  ! ensure all cap_h, lit_h values are available
  !----------------------------------------------------------------------

  !----------------------------------------------------------------------
  ! Drift (finite-kx) plus upwind dissipation
  !
        do i = ibeg, iend
           do m=1,n_stack
              z_der = 0.0
              z_dis = 0.0

              do i_diff=-m_dx,m_dx-i_dx
                 z_der = z_der+omega_dr(m,i,p_nek_loc,is)*&
                      w_d1(i_diff)*cap_h(m,i_loop(i+i_diff),p_nek_loc,is)

                 z_dis = z_dis+abs(omega_dr(m,i,p_nek_loc,is))*&
                      s_d1(i_diff)*lit_h(m,i_loop(i+i_diff),p_nek_loc,is)
              enddo

              ! Derivative plus dissipation
              rhs(m,i,p_nek_loc,is) = rhs(m,i,p_nek_loc,is)+&
                   z_der+z_dis

              ! Entropy tracking
              rhs_dr(m,i,p_nek_loc,is) = rhs_dr(m,i,p_nek_loc,is)+&
                   real(z_dis*conjg(cap_h(m,i,p_nek_loc,is)))

              ! Diagonal (ky) drift
              rhs(m,i,p_nek_loc,is) = rhs(m,i,p_nek_loc,is)+&
                   i_c*omega_d1(m,i,p_nek_loc,is)*cap_h(m,i,p_nek_loc,is)

              ! Diamagnetic drift
              rhs(m,i,p_nek_loc,is) = rhs(m,i,p_nek_loc,is)-&
                   i_c*omega_star(m,i,p_nek_loc,is)*gyro_u(m,i,p_nek_loc,is)
           enddo
        enddo !i

  !----------------------------------------------------------------------
  ! Buffer damping
  !
  if (boundary_method == 2) then
     do i = ibeg, iend
        rhs(:,i,p_nek_loc,is) = rhs(:,i,p_nek_loc,is)-explicit_damp_vec(is,i)*cap_h(:,i,p_nek_loc,is)
     enddo
  endif
  !----------------------------------------------------------------------
     enddo !p_nek_loc
  enddo !is
!$omp end parallel

  !----------------------------------------------------------------------
  ! Adaptive source
  !
  if (source_flag == 1) then

     call gyro_adaptive_source

!$omp parallel private(p_nek_loc,ie,is,p_nek,i)
     do is=1,n_kinetic

        p_nek_loc = 0

        do p_nek=1+i_proc_1,n_nek_1,n_proc_1

           p_nek_loc = p_nek_loc+1

           ie = nek_e(p_nek)  

           do i = ibeg, iend

              ! In this expression, nu_source = 0 if n > 0.
              ! (see gyro_radial_operators).
              rhs(:,i,p_nek_loc,is) = rhs(:,i,p_nek_loc,is) &
                   -nu_source*h0_eq(is,ie,i)

           enddo ! i

        enddo ! p_nek

     enddo ! is
!$omp end parallel

  endif
  !----------------------------------------------------------------------

  !---------------------------------------------------------
  ! Er shear
  !
!$omp parallel private(i)
   do i = ibeg, iend
      rhs(:,i,:,:) = rhs(:,i,:,:)-i_c*omega_eb_s(i)*h(:,i,:,:)
   enddo 
!$omp end parallel
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Krook ion-ion collision operator:
  !
  if (krook_flag == 1) then
     call gyro_collision_krook
!$omp parallel private(i)
     do i = ibeg, iend
        rhs(:,i,:,1) = rhs(:,i,:,1)+rhs_krook(:,i,:)
     end do
!$omp end parallel
  endif
  !---------------------------------------------------------

  deallocate(cap_h)
  deallocate(lit_h)

  call gyro_timer_out('RHS-total')

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'*[get_rhs_total done]'
  endif

end subroutine gyro_rhs_total
