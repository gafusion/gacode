! Verner 10:7(6) adaptive integrator

subroutine cgyro_step_gk_v76

  use timer_lib
  use mpi
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

  real, parameter :: a21  = .5d-2
  
  real, parameter :: a31  = -1.07679012345679012d0
  real, parameter :: a32  = 1.185679012345679012d0
  
  real, parameter :: a41  = 0.4083333333333333333d-1
  real, parameter :: a43 =  0.1225d0
  
  real, parameter :: a51 = 0.638913923625572678d0
  real, parameter :: a53 = -2.4556726382236568097d0
  real, parameter :: a54 = 2.27225871459808413161d0
  
  real, parameter :: a61   = -2.661577375018757131d0
  real, parameter :: a63   = 10.804513886456137696d0
  real, parameter :: a64  = -8.35391465739619941197d0
  real, parameter :: a65  = 0.8204875949566569791420d0
  
  real, parameter :: a71  = 6.067741434696770992718d0
  real, parameter :: a73  = -24.7112736359110857973d0
  real, parameter :: a74  =  20.427517930788893940467d0
  real, parameter :: a75  = -1.9061579788166471506241d0
  real, parameter :: a76  = 1.00617224924206801479004d0
  
  real, parameter :: a81  = 12.0546700762532029950911d0
  real, parameter :: a83  = -49.754784950468989328073d0
  real, parameter :: a84  = 41.14288863860467663259698d0
  real, parameter :: a85  = -4.46176014997400418564191d0
  real, parameter :: a86  = 2.0423348222391749598217172d0
  real, parameter :: a87  = -0.983484366540610737953080d-1

  real, parameter :: a91  = 10.138146522881807876418451d0
  real, parameter :: a93  = -42.64113603171750214622846d0
  real, parameter :: a94  = 35.7638400399225700713502118d0
  real, parameter :: a95  = -4.3480228403929076533403703d0
  real, parameter :: a96  = 2.00986226837703589544194359d0
  real, parameter :: a97  = 0.348749046033827240595382285d0
  real, parameter :: a98  = -0.27143900510483128423715871d0

  real, parameter :: a101  = -45.03007203429867712435322d0
  real, parameter :: a103  = 187.32724376545888407524182d0
  real, parameter :: a104  = -154.02882369350186905967286d0
  real, parameter :: a105  = 18.56465306347536233859492333d0
  real, parameter :: a106  = -7.141809679295078854925420497d0
  real, parameter :: a107  = 1.30880857816137862511476270601d0

  real, parameter :: b1 = 0.4715561848627222170431765108d-1
  real, parameter :: b4 = 0.2575056429843415189596436101d0
  real, parameter :: b5 = 0.26216653977412620477138630958d0
  real, parameter :: b6 = 0.15216092656738557403231331992d0
  real, parameter :: b7 = 0.49399691700324842469071758932d0
  real, parameter :: b8 = -0.29430311714032504415572447441d0
  real, parameter :: b9 = 0.8131747232495109999734599440137d-1
  real, parameter :: b10 = 0.d0

  real, parameter :: b1h = 0.446086066063411762873181759748d-1
  real, parameter :: b4h = 0.267164037857137268050910226094d0
  real, parameter :: b5h = 0.220101830017729301997971577665d0
  real, parameter :: b6h = 0.2188431703143156830983120833513d0
  real, parameter :: b7h = 0.2289871705411202883378173889764d0
  real, parameter :: b8h = 0.d0
  real, parameter :: b9h = 0.d0
  real, parameter :: b10h = 0.20295184663356282227670547938d-1

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
        call cgyro_vel_copy(h0_old, h_x)
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

     ! rhs(*,2) is not used again, so shift down all the subsequent indexes
     call cgyro_rhs(3-1,.TRUE.)
     
     call timer_lib_in('str')
     call cgyro_vel_fmaN(2,h_x, &
            h0_x, &
            (/ deltah2*a41, deltah2*a43 /), &
            rhs(:,:,:,1:(3-1)))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(4-1,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(3,h_x, &
            h0_x, &
            (/ deltah2*a51, deltah2*a53, deltah2*a54 /), &
            rhs(:,:,:,1:(4-1)))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(5-1,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(4,h_x, &
            h0_x, &
            (/ deltah2*a61, deltah2*a63, deltah2*a64, deltah2*a65 /), &
            rhs(:,:,:,1:(5-1)))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(6-1,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(5,h_x, &
            h0_x, &
            (/ deltah2*a71, deltah2*a73, deltah2*a74, deltah2*a75, deltah2*a76 /), &
            rhs(:,:,:,1:(6-1)))
     call timer_lib_out('str')
     
     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(7-1,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(6,h_x, &
            h0_x, &
            (/ deltah2*a81, deltah2*a83, deltah2*a84, deltah2*a85, &
               deltah2*a86, deltah2*a87 /), &
            rhs(:,:,:,1:(7-1)))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(8-1,.TRUE.)

     call timer_lib_in('str')
     call cgyro_vel_fmaN(7,h_x, &
            h0_x, &
            (/ deltah2*a91, deltah2*a93, deltah2*a94, deltah2*a95, &
               deltah2*a96, deltah2*a97, deltah2*a98 /), &
            rhs(:,:,:,1:(8-1)))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(9-1,.TRUE.)
     
     call timer_lib_in('str')
     call cgyro_vel_fmaN(6,h_x, &
            h0_x, &
            (/ deltah2*a101, deltah2*a103, deltah2*a104, deltah2*a105, &
               deltah2*a106, deltah2*a107 /), &
            rhs(:,:,:,1:(7-1)))
     call timer_lib_out('str')

     call cgyro_rhs_comm_async_hx
     call cgyro_field_c(.FALSE.)
     call cgyro_rhs(10-1,.TRUE.)

     !-------------------
     ! SOLUTION and ERROR
     !-------------------
     
     call timer_lib_in('str')
     ! using a multiplication by 0 in one element is still efffienct, since the matrix element was read for the 2nd equation
     call cgyro_vel_solution_werror(7, h_x, &
            h0_x, &
            deltah2*b1, rhs(:,:,:,1), &
            (/ deltah2*b4, deltah2*b5, deltah2*b6, deltah2*b7, &
               deltah2*b8, deltah2*b9, deltah2*b10 /), &
            rhs(:,:,:,(4-1):(10-1)), &
            deltah2*(b1-b1h), &
            (/ deltah2*(b4-b4h), deltah2*(b5-b5h), deltah2*(b6-b6h), deltah2*(b7-b7h), &
               deltah2*(b8-b8h), deltah2*(b9-b9h), deltah2*(b10-b10h) /), &
            error_hx, error_rhs)
     call timer_lib_out('str')

     call timer_lib_in('str_comm')
     error_x(1) = error_rhs
     error_x(2) = error_hx

     call MPI_ALLREDUCE(error_x,error_sum,2,MPI_DOUBLE_PRECISION, &
          MPI_SUM,CGYRO_COMM_WORLD,i_err)
     call timer_lib_out('str_comm')

     error_x = error_sum
     delta_x = error_x(1)+eps
     rel_error = error_x(1)/(error_x(2)+eps)
     var_error = sqrt(total_local_error+rel_error*rel_error)
    
     if (var_error < tol) then
        conv = 1
        delta_t_tot = delta_t_tot + deltah2
        total_local_error = total_local_error + rel_error*rel_error

        scale_x = max((tol/delta_x*1.0/delta_t)**(1.0/6.0), &
             (tol/delta_x*1.0/delta_t)**(1.0/7.0))
        
        deltah2 = max(min(scale_x,8.0),1.0)*deltah2
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
        call cgyro_error('Verner step exceeded max iteration count')
        return
     endif

  enddo

  ! caller expects the fields to be updated on exit
  call cgyro_field_c(.TRUE.)

  delta_t_gk = max(delta_t_last,6.0/7.0*deltah2)
  
  if (delta_t_last_step <= small) delta_t_last_step = delta_t_last

  if (delta_t_last_step < 0.1*delta_t_last) then
     delta_t_gk = delta_t_last+delta_t_last_step
  else
     if (delta_t_last_step/itrk < 0.1*delta_t_last) then
        delta_t_gk = delta_t_last+delta_t_last_step/itrk
     endif
  endif

  delta_t_gk = min(delta_t,delta_t_gk)
  total_local_error = var_error

end subroutine cgyro_step_gk_v76
