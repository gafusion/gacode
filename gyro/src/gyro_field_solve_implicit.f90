subroutine gyro_field_solve_implicit

  use gyro_globals
  use gyro_maxwell_private

  !---------------------------------------------------
  implicit none
  !---------------------------------------------------

  call gyro_operators_on_h

  ! Generate RHS for field solve
  do ix=1,n_field
     call gyro_velocity_sum(ix)
  enddo

  call gyro_timer_in('Field-implicit')
  call gyro_sparse_solve_umfpack(n_maxwell,n_maxwell_row,3,1)
  call gyro_timer_out('Field-implicit')

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_field_solve_implicit done]'
  endif

end subroutine gyro_field_solve_implicit
