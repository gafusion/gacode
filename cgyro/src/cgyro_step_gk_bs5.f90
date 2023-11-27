! Bogacki-Shampine (1996) 7:5(4) adaptive method

subroutine cgyro_step_gk_bs5

  use mpi
  use timer_lib
  use cgyro_globals
  use cgyro_globals_math
  use cgyro_io
  use cgyro_step

  implicit none

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

  ! Butcher table

  real, parameter :: a21  =  1.d0/6.d0
  
  real, parameter :: a31  = 2.d0/27.d0
  real, parameter :: a32  = 4.d0/27.d0
  
  real, parameter :: a41  = 183.d0/1372.d0
  real, parameter :: a42  = -162.d0/343.d0
  real, parameter :: a43  = 1053.d0/1372.d0
  
  real, parameter :: a51  =  68.0d0/297.d0
  real, parameter :: a52  = -4.d0/11.d0
  real, parameter :: a53  = 42.d0/143.d0
  real, parameter :: a54  = 1960.d0/3861.d0
  
  real, parameter :: a61  = 597.d0/22528.d0
  real, parameter :: a62  = 81.d0/352.d0
  real, parameter :: a63  = 63099.d0/585728.d0
  real, parameter :: a64  = 58653.d0/366080.d0
  real, parameter :: a65  = 4617.d0/20480.d0
  
  real, parameter :: a71  = 174197.d0/959244.d0
  real, parameter :: a72  = -30942.d0/79937.d0
  real, parameter :: a73  = 8152137.d0/19744439.d0
  real, parameter :: a74  = 666106.d0/1039181.d0
  real, parameter :: a75  = -29421.d0/29068.d0
  real, parameter :: a76  = 482048.d0/414219.d0
  
  real, parameter :: b1 =587.d0/8064.d0
  real, parameter :: b3 =4440339.d0/15491840.d0
  real, parameter :: b4 =24353.d0/124800.d0 
  real, parameter :: b5 =387.d0/44800.d0
  real, parameter :: b6 =2152.d0/5985.d0
  real, parameter :: b7 =7267.d0/94080.d0

  ! b1 - b1h
  
  real, parameter :: e1 = -3.d0/1280.d0
  real, parameter :: e3 = 6561.d0/632320.d0
  real, parameter :: e4 = -343.d0/20800.d0
  real, parameter :: e5 = 243.d0/12800.d0
  real, parameter :: e6 = -1.d0/95.d0

  logical :: is_first

  tol = error_tol

  itrk = 0
  conv = 0

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
 
  delta_t_last = 0.0
 
  call timer_lib_in('str_mem')
  call cgyro_vel_copy(h0_old, h_x)
  call timer_lib_out('str_mem')

  is_first = .TRUE.

  do while (delta_t_tot < delta_t .and. itrk <= itrk_max)
    
     call timer_lib_in('str')
     if (delta_t_tot + deltah2 > delta_t) then
        deltah2 = delta_t-delta_t_tot
        delta_t_last_step = deltah2
     else
        delta_t_last = deltah2
        deltah2_min = min(deltah2,deltah2_min)
        deltah2_max = max(deltah2,deltah2_max)
     endif

     if ((conv == 0) .and. (itrk >= 1)) then
        ! not converged so backing up
        call cgyro_vel_copy2(h0_x, h_x, h0_old)
     else
        call cgyro_vel_copy(h0_x, h_x)
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
     call cgyro_vel_fma2(h_x, h0_x, a21*deltah2, rhs(:,:,:,1))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(2,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(2,h_x, &
            h0_x, &
            (/ deltah2*a31, deltah2*a32 /), &
            rhs(:,:,:,1:2))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(3,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(3,h_x, &
            h0_x, &
            (/ deltah2*a41, deltah2*a42, deltah2*a43 /), &
            rhs(:,:,:,1:3))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(4,.TRUE.)
     
     call timer_lib_in('str')
     call cgyro_vel_fmaN(4,h_x, &
            h0_x, &
            (/ deltah2*a51, deltah2*a52, deltah2*a53, deltah2*a54 /), &
            rhs(:,:,:,1:4))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(5,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(5,h_x, &
            h0_x, &
            (/ deltah2*a61, deltah2*a62, deltah2*a63, deltah2*a64, deltah2*a65 /), &
            rhs(:,:,:,1:5))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(6,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(6,h_x, &
            h0_x, &
            (/ deltah2*a71, deltah2*a72, deltah2*a73, &
               deltah2*a74, deltah2*a75, deltah2*a76 /), &
            rhs(:,:,:,1:6))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(7,.TRUE.)

     !-------------------
     ! SOLUTION and ERROR
     !------------------

     call timer_lib_in('str')
     ! using a multiplication by 0 in one element is still efffienct, since the matrix element was read for the 1st equation
     call cgyro_vel_solution_werror(5, h_x, &
            h0_x, &
            deltah2*b1, rhs(:,:,:,1), &
            (/ deltah2*b3, deltah2*b4, deltah2*b5, deltah2*b6, deltah2*b7 /), &
            rhs(:,:,:,3:7), &
            deltah2*e1, &
            (/ deltah2*e3, deltah2*e4, deltah2*e5, deltah2*e6, 0.d0 /), &
            error_hx, error_rhs)
     call timer_lib_out('str')
   
     call timer_lib_in('str_comm')
     error_x(1) = error_rhs
     error_x(2) = error_hx

     call MPI_ALLREDUCE(error_x,error_sum,2,MPI_DOUBLE_PRECISION,&
          MPI_SUM,CGYRO_COMM_WORLD,i_err)
     call timer_lib_out('str_comm')
     
     error_x = error_sum
     delta_x = error_x(1)+eps
     rel_error = error_x(1)/(error_x(2)+eps)
     var_error = sqrt(total_local_error+rel_error*rel_error)
    
     if (var_error < tol) then
        conv = 1
        
        delta_t_tot = delta_t_tot+deltah2
        total_local_error = total_local_error + rel_error*rel_error

        scale_x = max((tol/delta_x*1.0/delta_t)**0.2, &
             (tol/delta_x*1.0/delta_t)**0.25)

        deltah2 = deltah2*max(1.0,min(6.0,scale_x))
        local_max_error = max(local_max_error,rel_error)

        call timer_lib_in('str_mem')
        call cgyro_vel_copy(h0_old, h0_x)
        call timer_lib_out('str_mem')
     else
        conv = 0
        deltah2 = 0.5*deltah2
     endif
     
     deltah2 = min(deltah2,delta_x_max)
     deltah2 = max(delta_x_min,deltah2)

     itrk = itrk+1

     if (itrk > itrk_max) then
        call cgyro_error('Bogacki-Shampine step exceeded max iteration count')
        return
     endif
     
  enddo

  ! caller expects the fields to be updated on exit
  call cgyro_field_c(.TRUE.)
  
  delta_t_gk = max(delta_t_last,4.0/5.0*deltah2)
  
  if (delta_t_last_step <= small) delta_t_last_step = delta_t_last

  if (delta_t_last_step < 0.1*delta_t_gk) then
     delta_t_gk = delta_t_last+delta_t_last_step
  else
     if (delta_t_last_step/itrk < 0.1*delta_t_gk) then
        delta_t_gk = delta_t_gk + delta_t_last_step/itrk
     endif
  endif

  delta_t_gk = min(delta_t,delta_t_gk)
  total_local_error = var_error

end subroutine cgyro_step_gk_bs5
