subroutine cgyro_step_gk

  use timer_lib
  use cgyro_globals
  use cgyro_globals_math

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

  call cgyro_rhs_comm_async_hx
  call timer_lib_in('str_mem')
  call cgyro_vel_copy(h0_x, h_x)
  call timer_lib_out('str_mem')

  
  ! Stage 1
  call cgyro_rhs(1,.FALSE.)
  call timer_lib_in('str')
  call cgyro_vel_fma2(h_x, h0_x, 0.5 * delta_t, rhs(:,:,:,1))
  call timer_lib_out('str')

  ! Stage 2
  call cgyro_rhs_comm_async_hx
  call cgyro_field_c(.FALSE.)
  call cgyro_rhs(2,.TRUE.)
  call timer_lib_in('str')
  call cgyro_vel_fma2(h_x, h0_x, 0.5 * delta_t, rhs(:,:,:,2))
  call timer_lib_out('str')

  ! Stage 3
  call cgyro_rhs_comm_async_hx
  call cgyro_field_c(.FALSE.)
  call cgyro_rhs(3,.TRUE.)
  call timer_lib_in('str')
  call cgyro_vel_fma2(h_x, h0_x, delta_t, rhs(:,:,:,3))
  call timer_lib_out('str')

  ! Stage 4
  call cgyro_rhs_comm_async_hx
  call cgyro_field_c(.FALSE.)
  call cgyro_rhs(4,.TRUE.)

  ! rhs(1) = 3rd-order error estimate
  ! h_x = h0_x + c1*r(1) + c2*r(2) + c3*r(3) + c4*r(4)
  ! r1 = h0_x + ec1*r(2) + ec2*r(3) - h_x
  !    = h0_x + ec1*r(2) + ec2*r(3) - h0_x - c1*r(1) - c2*r(2) - c3*r(3) - c4*r(4)
  !    = -c1*r(1) + (ec1-c2)*r(2) + (ec2-c3)*r(3) - c4*r(4)
  ! delta_t/6 = delta_t/3 - delta_t/6
  call timer_lib_in('str')
  call cgyro_vel_solution_werror(3,h_x, &
            h0_x, &
            delta_t/6, rhs(:,:,:,1), &
            (/ 2*delta_t/6, 2*delta_t/6, delta_t/6 /), &
            rhs(:,:,:,2:4), &
            -delta_t/6, &
            (/ 0.d0, 2*delta_t/6, -delta_t/6 /))

  call timer_lib_out('str')

  call cgyro_field_c(.TRUE.)

end subroutine cgyro_step_gk

