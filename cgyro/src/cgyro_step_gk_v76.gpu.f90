subroutine cgyro_step_gk_v76

  use timer_lib
  use mpi
  use cgyro_globals

  implicit none

  ! V7(6), Vernier, ?2010 paper
  ! efficient
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

  integer converged, conv, rk_MAX , iiter
  double precision orig_delta_x_t, total_delta_x_step, delta_x_min, delta_x_max
  double precision error_sum(2), error_x(2)
  double precision delta_x, tau, deltah2, local_max_error
  double precision rel_error, delta_t_last
  double precision delta_t_last_step
  double precision deltah2_max, deltah2_min
  double precision var_error, scale_x, tol
  double precision gt_error(2)
  double precision error_rhs, error_hx
  double precision local_error_rhs, local_error_hx


  !! butcher table

  double precision, parameter :: a21  = .5e-2
  
  double precision, parameter :: a31  = -1.07679012345679012
  double precision, parameter :: a32  = 1.185679012345679012
  
  double precision, parameter :: a41  = .4083333333333333333e-1
  double precision, parameter :: a42  = 0.
  double precision, parameter :: a43 =  .1225
  
  double precision, parameter :: a51 =  .638913923625572678
  double precision, parameter :: a52 = 0.
  double precision, parameter :: a53 = -2.4556726382236568097
  double precision, parameter :: a54 = 2.27225871459808413161
  
  double precision, parameter :: a61   = -2.661577375018757131
  double precision, parameter :: a62   = 0.
  double precision, parameter :: a63   = 10.804513886456137696
  double precision, parameter :: a64  = -8.35391465739619941197
  double precision, parameter :: a65  = .8204875949566569791420
  
  double precision, parameter :: a71  = 6.067741434696770992718
  double precision, parameter :: a72  = 0.
  double precision, parameter :: a73  = -24.7112736359110857973
  double precision, parameter :: a74  =  20.427517930788893940467
  double precision, parameter :: a75  = -1.9061579788166471506241
  double precision, parameter :: a76  = 1.00617224924206801479004
  
  double precision, parameter :: a81  = 12.0546700762532029950911
  double precision, parameter :: a82  = 0.
  double precision, parameter :: a83  = -49.754784950468989328073
  double precision, parameter :: a84  = 41.14288863860467663259698
  double precision, parameter :: a85  = -4.46176014997400418564191
  double precision, parameter :: a86  = 2.0423348222391749598217172
  double precision, parameter :: a87  = -0.983484366540610737953080e-1

  double precision, parameter :: a91  = 10.138146522881807876418451
  double precision, parameter :: a92  = 0.
  double precision, parameter :: a93  = -42.64113603171750214622846
  double precision, parameter :: a94  = 35.7638400399225700713502118
  double precision, parameter :: a95  = -4.3480228403929076533403703
  double precision, parameter :: a96  = 2.00986226837703589544194359
  double precision, parameter :: a97  = .348749046033827240595382285
  double precision, parameter :: a98  = -.27143900510483128423715871

  double precision, parameter :: a101  = -45.03007203429867712435322
  double precision, parameter :: a102  = 0.
  double precision, parameter :: a103  = 187.32724376545888407524182
  double precision, parameter :: a104  = -154.02882369350186905967286
  double precision, parameter :: a105  = 18.56465306347536233859492333
  double precision, parameter :: a106  = -7.141809679295078854925420497
  double precision, parameter :: a107  = 1.30880857816137862511476270601
  double precision, parameter :: a108 = 0.
  double precision, parameter :: a109 = 0.

  double precision, parameter :: c1   = 0.
  double precision, parameter :: c2   = 5.e-2
  double precision, parameter :: c3   =.108888888888888888888888889
  double precision, parameter :: c4   =.163333333333333333333333333
  double precision, parameter :: c5   =.4555
  double precision, parameter :: c6   =.609509448997838131708700442
  double precision, parameter :: c7   =.884
  double precision, parameter :: c8   =.925
  double precision, parameter :: c9   = 1.
  double precision, parameter :: c10   =1.

  double precision, parameter :: b1 = .4715561848627222170431765108e-1
  double precision, parameter :: b2 = 0.
  double precision, parameter :: b3 = 0.
  double precision, parameter :: b4 = .2575056429843415189596436101
  double precision, parameter :: b5 = .26216653977412620477138630958
  double precision, parameter :: b6 = .15216092656738557403231331992
  double precision, parameter :: b7 = .49399691700324842469071758932
  double precision, parameter :: b8 = -.29430311714032504415572447441
  double precision, parameter :: b9 =.8131747232495109999734599440137e-1
  double precision, parameter :: b10 = 0.

  double precision, parameter :: b1h = .446086066063411762873181759748e-1
  double precision, parameter :: b2h = 0.
  double precision, parameter :: b3h = 0.
  double precision, parameter :: b4h = 0.267164037857137268050910226094
  double precision, parameter :: b5h = 0.220101830017729301997971577665
  double precision, parameter :: b6h = 0.2188431703143156830983120833513
  double precision, parameter :: b7h = 0.2289871705411202883378173889764
  double precision, parameter :: b8h = 0.
  double precision, parameter :: b9h = 0.
  double precision, parameter :: b10h = 0.20295184663356282227670547938e-1

  double precision, parameter :: e1 = b1-b1h
  double precision, parameter :: e2 = b2-b2h
  double precision, parameter :: e3 = b3-b3h
  double precision, parameter :: e4 = b4-b4h
  double precision, parameter :: e5 = b5-b5h
  double precision, parameter :: e6 = b6-b6h
  double precision, parameter :: e7 = b7-b7h
  double precision, parameter :: e8 = b8-b8h
  double precision, parameter :: e9 = b9-b9h
  double precision, parameter :: e10 = b10-b10h


  double precision, parameter :: EPS  = 2.2d-12

!$acc parallel loop collapse(2) independent present(h0_old)
    do iv_loc=1,nv_loc
       do ic_loc=1,nc
          h0_old(ic_loc,iv_loc) = 0.
       enddo
    enddo

  call timer_lib_in('str')

  local_max_error = 0.
  delta_t_last = 0.
  delta_t_last_step = 0.
  orig_delta_x_t = delta_t

  iiter = 0
  total_delta_x_step = 0.
  total_local_error = 0.
  var_error = 0.

  tol = delta_t_tol
  deltah2 = delta_t_gk

  if ( delta_t_gk .lt. 1.d-10) then
     deltah2 = orig_delta_x_t
     delta_t_gk = deltah2
  endif
  
  delta_x_min = 1.e-10*orig_delta_x_t
  delta_x_max = orig_delta_x_t

  converged = 0
  conv = 0
  total_local_error = 0.
  
  rk_MAX = 10000

!$acc parallel loop collapse(2) independent present(h0_old,h_x)
  do iv_loc=1,nv_loc
     do ic_loc=1,nc
        h0_old(ic_loc,iv_loc) = h_x(ic_loc,iv_loc)
     enddo
  enddo

  conv = 0
  delta_t_gk = 0.
  deltah2_min = 1.d10
  deltah2_max = 0.

  do while (total_delta_x_step .lt. orig_delta_x_t .and. iiter .le. rk_MAX )
    
     if ( total_delta_x_step + deltah2 .gt. orig_delta_x_t ) then
        deltah2 = orig_delta_x_t - total_delta_x_step
        delta_t_last_step = deltah2
     else
        delta_t_last = deltah2
        deltah2_min = min(deltah2, deltah2_min)
        deltah2_max = max(deltah2, deltah2_max)
     endif
     
     call timer_lib_in('str_mem')             
     if (( conv .eq. 0 ) .and. (iiter .ge. 1)) then
        
        call timer_lib_in('str_mem')

!$acc parallel loop collapse(2) independent present(h0_x,h0_old,h_x)
        do iv_loc=1,nv_loc
           do ic_loc=1,nc
              h0_x(ic_loc,iv_loc) = h0_old(ic_loc,iv_loc)
              h_x(ic_loc,iv_loc) = h0_old(ic_loc,iv_loc)
           enddo
        enddo
  
        call timer_lib_out('str_mem')        
     else
        call timer_lib_in('str_mem')
        
!$acc parallel loop collapse(2) independent present(h0_x,h_x)
        do iv_loc=1,nv_loc
           do ic_loc=1,nc
              h0_x(ic_loc,iv_loc) = h_x(ic_loc,iv_loc)
           enddo
        enddo
     endif

     call timer_lib_out('str_mem')        
     
     call cgyro_field_c_gpu
     
     if ( i_proc .eq. 0 ) write(*,*) " paper v76_effi deltah2 ", iiter, deltah2

     ! Stage 1
     !
     call cgyro_rhs(1)

     call timer_lib_in('str')

!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs(:,:,1))
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + a21*deltah2*rhs(ic_loc,iv_loc,1)
        enddo
     enddo

     call timer_lib_out('str')
     call cgyro_field_c_gpu

     ! Stage 2 ! k2

     call cgyro_rhs(2)
     call timer_lib_in('str')
     
!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*(a31*rhs(ic_loc,iv_loc,1) &
                + a32*rhs(ic_loc,iv_loc,2))
        enddo
     enddo
     call timer_lib_out('str')
     call cgyro_field_c_gpu

     ! Stage 3
     
     call cgyro_rhs(3)
     call timer_lib_in('str')     
     ! k4
!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc, iv_loc) &
                + deltah2*(a41*rhs(ic_loc, iv_loc, 1) &
                + a43*rhs(ic_loc,iv_loc,3))
        enddo
     enddo
     call timer_lib_out('str')

     call cgyro_field_c_gpu
     
     ! Stage 4
     call cgyro_rhs(4)
     !  k5

!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc)  &
                + deltah2*(a51*rhs(ic_loc, iv_loc, 1) &
                + a53*rhs(ic_loc,iv_loc,3) &
                + a54*rhs(ic_loc,iv_loc,4))
        enddo
     enddo

     call cgyro_field_c_gpu
     
     !  stage 5
     
     call cgyro_rhs(5)

!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*(a61*rhs(ic_loc, iv_loc, 1) &
                + a63*rhs(ic_loc,iv_loc,3) &
                + a64*rhs(ic_loc,iv_loc,4) &
                + a65*rhs(ic_loc,iv_loc,5))
        enddo
     enddo

     call cgyro_field_c_gpu

     ! stage 6
     call cgyro_rhs(6)

!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) + &
                deltah2*(a71*rhs(ic_loc, iv_loc, 1) &
                + a73*rhs(ic_loc,iv_loc,3) &
                + a74*rhs(ic_loc,iv_loc,4) &
                + a75*rhs(ic_loc,iv_loc,5) &
                + a76*rhs(ic_loc,iv_loc,6))
        enddo
     enddo

     
     call cgyro_field_c_gpu
     call cgyro_rhs(7)

     !! soln = h_x of order 4
     !! error_x(2) = sum(abs(h0_x(ic_loc,iv_loc)))

!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc)  &
                + deltah2*(a81*rhs(ic_loc, iv_loc, 1) &
                + a83*rhs(ic_loc,iv_loc,3) &
                + a84*rhs(ic_loc,iv_loc,4) &
                + a85*rhs(ic_loc,iv_loc,5) &
                + a86*rhs(ic_loc,iv_loc,6) &
                + a87*rhs(ic_loc,iv_loc,7))
        enddo
     enddo

     call cgyro_field_c_gpu
     call cgyro_rhs(8)

!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc)  &
                + deltah2*(a91*rhs(ic_loc, iv_loc, 1) &
                + a93*rhs(ic_loc,iv_loc,3) &
                + a94*rhs(ic_loc,iv_loc,4) &
                + a95*rhs(ic_loc,iv_loc,5) &
                + a96*rhs(ic_loc,iv_loc,6) &
                + a97*rhs(ic_loc,iv_loc,7) &
                + a98*rhs(ic_loc,iv_loc,8))
        enddo
     enddo

     call cgyro_field_c_gpu
     call cgyro_rhs(9)
     
!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc)  &
                + deltah2*(a101*rhs(ic_loc, iv_loc, 1) &
                + a103*rhs(ic_loc,iv_loc,3) &
                + a104*rhs(ic_loc,iv_loc,4) &
                + a105*rhs(ic_loc,iv_loc,5) &
                + a106*rhs(ic_loc,iv_loc,6) &
                + a107*rhs(ic_loc,iv_loc,7))
        enddo
     enddo

     call cgyro_field_c_gpu
     call cgyro_rhs(10)

     ! order 7 solution
     
!$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*(b1*rhs(ic_loc, iv_loc, 1) &
                + b4*rhs(ic_loc,iv_loc,4) &
                + b5*rhs(ic_loc,iv_loc,5) &
                + b6*rhs(ic_loc,iv_loc,6) &
                + b7*rhs(ic_loc,iv_loc,7) &
                + b8*rhs(ic_loc,iv_loc,8) &
                + b9*rhs(ic_loc,iv_loc,9))
        enddo
     enddo
     call cgyro_field_c_gpu

     error_rhs = 0.


!! !$acc parallel loop collapse(2) independent present(h0_x,h_x,rhs) reduction(+:error_rhs)

!$acc parallel loop collapse(2) gang present(rhs) reduction(+:error_rhs)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           rhs(ic_loc, iv_loc, 1) = deltah2*( &
                e1*rhs(ic_loc, iv_loc, 1) &
                + e4*rhs(ic_loc,iv_loc,4) &
                + e5*rhs(ic_loc,iv_loc,5) &
                + e6*rhs(ic_loc,iv_loc,6) &
                + e7*rhs(ic_loc,iv_loc,7) &
                + e8*rhs(ic_loc,iv_loc,8) &
                + e9*rhs(ic_loc,iv_loc,9) &
                + e10*rhs( ic_loc, iv_loc, 10))
           error_rhs = error_rhs + abs(rhs(ic_loc, iv_loc, 1))
        enddo
     enddo

     error_hx = 0.

!$acc parallel loop collapse(2) independent present(h_x) reduction(+:error_hx)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           error_hx = error_hx + abs(h_x(ic_loc,iv_loc))
        enddo
     enddo

     error_sum = 0.
     error_x = 0.

     !! error_x(1) = sum(abs(rhs(:,:,1)))
     !! error_x(2) = sum(abs(h_x))

     error_x(1) = error_rhs
     error_x(2) = error_hx

     call timer_lib_in('str_comm')
     
     call MPI_ALLREDUCE(error_x, error_sum, 2, MPI_DOUBLE_PRECISION, &
          MPI_SUM, MPI_COMM_WORLD, i_err)

     error_x = error_sum

     call timer_lib_out('str_comm')     

     if ( i_proc == 0 ) &
          write(*,*) " paper V76effic **** rhs_error ", &
          error_x(1), " hx_error ", error_x(2)

     delta_x = error_x(1)
     tau = tol*max(error_x(2), 1.)

     rel_error = error_x(1)/(error_x(2)+EPS)
     var_error = sqrt(total_local_error + rel_error*rel_error)

     !! method 1 local error
     
     !!if ( error_x(1) .lt. tau ) then
        
     !! method 2 "variance" error
     
     if ( var_error .lt. tol ) then
        
!!        if ( i_proc == 0 ) &
!!             write(*,*) " paper V76effic **** local error mode ", &
!!             rel_error, " variance error", var_error


!!paper        if ( i_proc == 0 ) &
        !! write(*,*) " dt ", deltah2, " V76 **** var error mode ", rel_error, " variance error", var_error

!$acc parallel loop collapse(2) independent present(h0_x,h0_old)
        do iv_loc=1,nv_loc
           do ic_loc=1,nc
              h0_old(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc)
           enddo
        enddo

        call cgyro_field_c_gpu

        converged = converged + 1
        conv = 1
        total_delta_x_step = total_delta_x_step + deltah2
        total_local_error = total_local_error + rel_error*rel_error

        scale_x = .95*max((tol/(delta_x + EPS )*1./delta_t)**(1./6.), &
             (tol/(delta_x + EPS )*1./delta_t)**(1./7.))

        scale_x = max(min(scale_x, 8.), 1.)
        
        deltah2 = scale_x*deltah2
        
        local_max_error = max(local_max_error, rel_error)

     else
        conv = 0
        deltah2 = .5*deltah2
        if (i_proc .eq. 0 ) then
           write(*,*) " v76 ***  error backing up, not converged ", &
                " total delta_x step ", total_delta_x_step, " new deltah2 ", deltah2
        endif
        flush(6)
     endif

     deltah2 = min(deltah2, delta_x_max)
     deltah2 = max(delta_x_min, deltah2)

     iiter = iiter + 1

     if ( iiter .gt. rk_MAX ) then
        if ( i_proc .eq. 0 ) then
           write(*,*) " v76 max count  exceeded, rel_error, var_error ", &
                iiter, rel_error, var_error
           flush(6)           
        endif
        stop
     endif
  enddo
  
  call timer_lib_out('str')

1000 continue

  delta_t_gk = delta_t_last
  if ( delta_t_last_step .lt.  .5*delta_t_last )  & 
       delta_t_gk = delta_t_last + delta_t_last_step
  
  total_local_error = var_error

  !!
  !! if ( i_proc == 0 ) &
  !! write(*,*) i_proc , " paper v76 deltah2_min, max converged ", deltah2_min, deltah2_max
  !!  
  ! Filter special spectral components
  
  call cgyro_filter
!!

end subroutine cgyro_step_gk_v76
