subroutine cgyro_step_gk

  use timer_lib
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

  call timer_lib_in('str_mem')

!$acc parallel loop collapse(2) independent present(h0_x,h_x)
  do iv_loc=1,nv_loc
     do ic_loc=1,nc
       h0_x(ic_loc,iv_loc) = h_x(ic_loc,iv_loc)
     enddo
  enddo

  call timer_lib_out('str_mem')

  
  ! Stage 1
  call cgyro_rhs(1)
  call timer_lib_in('str')
!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs(:,:,1))
  do iv_loc=1,nv_loc
     do ic_loc=1,nc
       h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) + 0.5 * delta_t * rhs(ic_loc,iv_loc,1)
     enddo
  enddo
  call timer_lib_out('str')
  call cgyro_field_c_gpu

  ! Stage 2
  call cgyro_rhs(2)
  call timer_lib_in('str')
!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs(:,:,2))
  do iv_loc=1,nv_loc
     do ic_loc=1,nc
       h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) + 0.5 * delta_t * rhs(ic_loc,iv_loc,2)
     enddo
  enddo
  call timer_lib_out('str')
  call cgyro_field_c_gpu

  ! Stage 3
  call cgyro_rhs(3)
  call timer_lib_in('str')
!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs(:,:,3))
  do iv_loc=1,nv_loc
     do ic_loc=1,nc
        h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) + delta_t * rhs(ic_loc,iv_loc,3)
     enddo
  enddo
  call timer_lib_out('str')
  call cgyro_field_c_gpu

  ! Stage 4
  call cgyro_rhs(4)
  call timer_lib_in('str')
!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
  do iv_loc=1,nv_loc
     do ic_loc=1,nc
       h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                          + delta_t*( rhs(ic_loc,iv_loc,1)+2*rhs(ic_loc,iv_loc,2)+ &
                                      2*rhs(ic_loc,iv_loc,3)+rhs(ic_loc,iv_loc,4) )/6  
     enddo
  enddo
  call timer_lib_out('str')
  call cgyro_field_c_gpu

  ! rhs(1) = 3rd-order error estimate
  call timer_lib_in('str')
!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
  do iv_loc=1,nv_loc
     do ic_loc=1,nc
       rhs(ic_loc,iv_loc,1) = h0_x(ic_loc,iv_loc) &
                            + delta_t*(rhs(ic_loc,iv_loc,2)+2*rhs(ic_loc,iv_loc,3))/3 &
                            - h_x(ic_loc,iv_loc)
     enddo
  enddo

  ! Filter special spectral components
  call cgyro_filter_gpu
  
  call timer_lib_out('str')

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

subroutine cgyro_filter_gpu

  use cgyro_globals

  implicit none

  integer :: ir

  if (zf_test_mode == 0 .and. n == 0) then
!$acc parallel loop gang vector private(ir) &
!$acc          present(ir_c,px,h_x,cap_h_c) default(none)
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
!$acc parallel loop gang vector private(ir) &
!$acc          present(ir_c,h_x,cap_h_c) default(none)
     do ic=1,nc
        ir = ir_c(ic)
        if (ir == 1) then
           h_x(ic,:)     = 0.0
           cap_h_c(ic,:) = 0.0
        endif
     enddo
  endif

end subroutine cgyro_filter_gpu

