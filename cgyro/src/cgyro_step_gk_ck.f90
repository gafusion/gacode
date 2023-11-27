! Cash-Karp 6:5(4) adaptive integrator

subroutine cgyro_step_gk_ck

  use timer_lib
  use mpi
  use cgyro_globals
  use cgyro_globals_math
  use cgyro_io
  use cgyro_step

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

  logical :: is_first

  tol = error_tol

  ! total iteration counter
  itrk = 0

  ! good step counter
  istep = 0

  conv = 1

  scale_x = 0.0
  deltah2 = delta_t_gk
  delta_t_last_step = 0.0

  delta_x_min = delta_t*1e-10
  delta_x_max = delta_t

  delta_t_tot = 0.0
  total_local_error = 0.0
  local_max_error = 0.0

  deltah2_min = 1.0
  deltah2_max = 0.0

  delta_t_last = deltah2

  call timer_lib_in('str_mem')
  call cgyro_vel_copy(h0_old, h_x)
  call timer_lib_out('str_mem')

  is_first = .TRUE.

  do while (delta_t_tot < delta_t)

     call timer_lib_in('str')
     if (delta_t_tot + deltah2 > delta_t) then
        deltah2 = delta_t-delta_t_tot
        delta_t_last_step = deltah2
     else
        delta_t_last = deltah2
        deltah2_min = min(deltah2, deltah2_min)
        deltah2_max = max(deltah2, deltah2_max)
     endif

     if (conv == 1) then
        call cgyro_vel_copy(h0_x, h_x)
     else
        call cgyro_vel_copy2(h0_x, h_x, h0_old)
     endif
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     if (is_first) then
        ! fields already in good shape in the beginning
        call cgyro_rhs(1,.FALSE.)
        is_first = .FALSE.
     else
        call cgyro_field_c(.FALSE.)
        call cgyro_rhs(1,.TRUE.)
     endif

     call timer_lib_in('str')
     call cgyro_vel_fma2(h_x, h0_x, 0.2d0*deltah2, rhs(:,:,:,1))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(2,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fma3(h_x, &
            h0_x, &
            1.0/40.0*deltah2*3.0, rhs(:,:,:,1), &
            1.0/40.0*deltah2*9.0, rhs(:,:,:,2))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(3,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(3, h_x, &
            h0_x, &
            (/ deltah2*( 3.d0/10.d0), &
               deltah2*(-9.d0/10.d0), &
               deltah2*( 6.d0/ 5.d0) /), &
             rhs(:,:,:,1:3))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(4,.TRUE.) 

     call timer_lib_in('str')
     call cgyro_vel_fmaN(4, h_x, &
            h0_x, &
            (/ deltah2*(-11.d0/54.d0), &
               deltah2*(  5.d0/ 2.d0), &
               deltah2*(-70.d0/27.d0), &
               deltah2*( 35.d0/27.d0) /), &
            rhs(:,:,:,1:4))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(5,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(5, h_x, &
            h0_x, &
            (/ deltah2*( 1631.d0/ 55296.d0), &
               deltah2*(  175.d0/   512.d0), &
               deltah2*(  575.d0/ 13824.d0), &
               deltah2*(44275.d0/110592.d0), &
               deltah2*(  253.d0/  4096.d0) /), &
            rhs(:,:,:,1:5))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(6,.TRUE.)

     !-------------------
     ! SOLUTION and ERROR
     !-------------------

     call timer_lib_in('str')
     ! using a multiplication by 0 in one element is still efffienct, since the matrix element was read for the 2nd equation
     call cgyro_vel_solution_werror(4, h_x, &
            h0_x, &
            deltah2*(    37.d0/ 378.d0), rhs(:,:,:,1), &
            (/ deltah2*(250.d0/ 621.d0), &
               deltah2*(125.d0/ 594.d0), 0.d0, &
               deltah2*(512.d0/1771.d0) /), &
            rhs(:,:,:,3:6), &
            deltah2*(     37.d0/378.d0- 2825.d0/27648.d0), &
            (/ deltah2*( 250.d0/621.d0-18575.d0/48384.d0), &
               deltah2*( 125.d0/594.d0-13525.d0/55296.d0), &
               deltah2*(-277.d0/14336.d0), &
               deltah2*( 512.d0/1771.d0-1.d0/4.d0) /), &
            error_hx, error_rhs)
     call timer_lib_out('str')

     call timer_lib_in('str_comm')
     error_x(1) = error_rhs
     error_x(2) = error_hx

     call MPI_ALLREDUCE(error_x,error_sum,2,MPI_DOUBLE_PRECISION,&
          MPI_SUM,CGYRO_COMM_WORLD,i_err)

     error_x = error_sum
     delta_x = error_x(1)+eps
     rel_error = error_x(1)/(error_x(2)+eps)
     var_error = sqrt(total_local_error+rel_error*rel_error)
     call timer_lib_out('str_comm')

     if (var_error < tol) then

        istep = istep+1
        deltah2_vec(istep) = deltah2

        delta_t_tot = delta_t_tot + deltah2
        total_local_error = total_local_error + rel_error*rel_error

        scale_x = max((tol/delta_x*1.0/delta_t)**0.2, &
             (tol/delta_x*1.0/delta_t)**0.25)

        deltah2 = deltah2*max(1.0,min(6.0,scale_x))

        local_max_error = max(local_max_error,rel_error)
        conv = 1

        call timer_lib_in('str_mem')
        call cgyro_vel_copy(h0_old, h0_x)
        call timer_lib_out('str_mem')

     else
        deltah2 = 0.5*deltah2
        conv = 0
     endif

     deltah2 = min(deltah2,delta_x_max)
     deltah2 = max(delta_x_min,deltah2)

     itrk = itrk+1

     if (itrk > itrk_max) then
        call cgyro_error('Cash-Carp step exceeded max iteration count')
        return
     endif

  enddo

  ! caller expects the fields to be updated on exit
  call cgyro_field_c(.TRUE.)

  delta_t_gk = max(delta_t_last,deltah2)

  if (delta_t_last_step <= small) delta_t_last_step = delta_t_last

  if (delta_t_last_step < 0.1*delta_t_gk) then
     delta_t_gk = delta_t_last+delta_t_last_step
  else
     if (delta_t_last_step/itrk < 0.1*delta_t_gk) then
        delta_t_gk = delta_t_gk+delta_t_last_step/itrk
     endif
  endif
  
  delta_t_gk = min(delta_t,delta_t_gk)
  
  total_local_error = var_error

end subroutine cgyro_step_gk_ck
