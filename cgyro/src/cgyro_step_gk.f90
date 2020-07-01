! RK4 time-advance for the distribution 

subroutine cgyro_step_gk

  use timer_lib
  use cgyro_globals

  implicit none

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

  call timer_lib_in('str_mem')
!$omp parallel do collapse(2)
  do iv_loc=1,nv_loc
     do ic=1,nc
       h0_x(ic,iv_loc) = h_x(ic,iv_loc)
     enddo
  enddo
  call timer_lib_out('str_mem')

  
  ! Stage 1
  call cgyro_rhs(1)
  call timer_lib_in('str')
!$omp parallel do collapse(2)
  do iv_loc=1,nv_loc
     do ic=1,nc
       h_x(ic,iv_loc) = h0_x(ic,iv_loc) + 0.5 * delta_t * rhs(ic,iv_loc,1)
     enddo
  enddo
  call timer_lib_out('str')
  call cgyro_field_c

  ! Stage 2
  call cgyro_rhs(2)
  call timer_lib_in('str')
!$omp parallel do collapse(2)
  do iv_loc=1,nv_loc
     do ic=1,nc
       h_x(ic,iv_loc) = h0_x(ic,iv_loc) + 0.5 * delta_t * rhs(ic,iv_loc,2)
     enddo
  enddo
  call timer_lib_out('str')
  call cgyro_field_c

  ! Stage 3
  call cgyro_rhs(3)
  call timer_lib_in('str')
!$omp parallel do collapse(2)
  do iv_loc=1,nv_loc
     do ic=1,nc
        h_x(ic,iv_loc) = h0_x(ic,iv_loc) + delta_t * rhs(ic,iv_loc,3)
     enddo
  enddo
  call timer_lib_out('str')
  call cgyro_field_c

  ! Stage 4
  call cgyro_rhs(4)
  call timer_lib_in('str')
!$omp parallel do collapse(2)
  do iv_loc=1,nv_loc
     do ic=1,nc
        h_x(ic,iv_loc) = h0_x(ic,iv_loc)+delta_t*(&
                rhs(ic,iv_loc,1) &
             +2*rhs(ic,iv_loc,2) &
             +2*rhs(ic,iv_loc,3) &
               +rhs(ic,iv_loc,4))/6  
     enddo
  enddo
  call timer_lib_out('str')
  call cgyro_field_c

  ! rhs(1) = 3rd-order error estimate
  call timer_lib_in('str')
!$omp parallel do collapse(2)
  do iv_loc=1,nv_loc
     do ic=1,nc
        rhs(ic,iv_loc,1) = h0_x(ic,iv_loc) +delta_t*( &
             rhs(ic,iv_loc,2)&
            +2*rhs(ic,iv_loc,3))/3 &
                            - h_x(ic,iv_loc)
     enddo
  enddo
  call timer_lib_out('str')
    
end subroutine cgyro_step_gk
