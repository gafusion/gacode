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
  
subroutine cgyro_filter

  use cgyro_globals

  implicit none

  integer :: ir
  
  if (zf_test_mode == 0 .and. n == 0) then
     do ic=1,nc
        ir = ir_c(ic) 
        if (ir == 1 .or. px(ir) == 0) then
           h_x(ic,:)     = 0.0
           cap_h_c(ic,:) = 0.0
        endif
     enddo
  endif

  ! Remove p=-M (is this ever useful?)
  if (psym_flag == 1) then
     do ic=1,nc
        ir = ir_c(ic) 
        if (ir == 1) then
           h_x(ic,:)     = 0.0
           cap_h_c(ic,:) = 0.0
        endif
     enddo
  endif

end subroutine cgyro_filter

