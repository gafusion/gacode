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

