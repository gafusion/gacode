! Bogacki-Shampine (1996) 7:5(4) adaptive method  |  multithreaded version

subroutine cgyro_step_gk_bs5

  use mpi
  use timer_lib
  use cgyro_globals
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
  real, parameter :: a43 =  1053.d0/1372.d0
  
  real, parameter :: a51 =  68.0d0/297.d0
  real, parameter :: a52 = -4.d0/11.d0
  real, parameter :: a53 = 42.d0/143.d0
  real, parameter :: a54 = 1960.d0/3861.d0
  
  real, parameter :: a61   = 597.d0/22528.d0
  real, parameter :: a62   = 81.d0/352.d0
  real, parameter :: a63   = 63099.d0/585728.d0
  real, parameter :: a64  = 58653.d0/366080.d0
  real, parameter :: a65  = 4617.d0/20480.d0
  
  real, parameter :: a71  = 174197.d0/959244.d0
  real, parameter :: a72  = -30942.d0/79937.d0
  real, parameter :: a73  = 8152137.d0/19744439.d0
  real, parameter :: a74  = 666106.d0/1039181.d0
  real, parameter :: a75  = -29421.d0/29068.d0
  real, parameter :: a76  = 482048.d0/414219.d0
  
  real, parameter :: a81  = 587.d0/8064.d0
  real, parameter :: a82  = 0.d0
  real, parameter :: a83  = 4440339.d0/15491840.d0
  real, parameter :: a84  = 24353.d0/124800.d0
  real, parameter :: a85  = 387.d0/44800.d0
  real, parameter :: a86  = 2152.d0/5985.d0
  real, parameter :: a87  = 7267.d0/94080.d0

  real, parameter :: c2   = 1.d0/6.d0
  real, parameter :: c3   = 2.d0/9.d0
  real, parameter :: c4   = 3.d0/7.d0
  real, parameter :: c5   = 2.d0/3.d0
  real, parameter :: c6   = 3.d0/4.d0
  real, parameter :: c7   = 1.d0
  real, parameter :: c8   = 1.d0

  real, parameter :: b1 =587.d0/8064.d0
  real, parameter :: b2 =0.d0
  real, parameter :: b3 =4440339.d0/15491840.d0
  real, parameter :: b4 =24353.d0/124800.d0 
  real, parameter :: b5 =387.d0/44800.d0
  real, parameter :: b6 =2152.d0/5985.d0
  real, parameter :: b7 =7267.d0/94080.d0
  real, parameter :: b8 = 0.d0

  real, parameter :: b1h =6059.d0/80640.d0
  real, parameter :: b2h =0.d0
  real, parameter :: b3h =8559189.d0/30983680.d0
  real, parameter :: b4h =26411.d0/124800.d0
  real, parameter :: b5h =-927.d0/89600.d0
  real, parameter :: b6h =443.d0/1197.d0
  real, parameter :: b7h =7267.d0/94080.d0
  real, parameter :: b8h =0.d0

  real, parameter :: b1p = 2479.d0/34992.d0
  real, parameter :: b2p = 0.d0
  real, parameter :: b3p = 123.d0/416.d0
  real, parameter :: b4p = 612941.d0/3411720.d0
  real, parameter :: b5p = 43.d0/1440.d0
  real, parameter :: b6p = 2272.d0/6561.d0
  real, parameter :: b7p = 79937.d0/1113912.d0
  real, parameter :: b8p = 3293.d0/556956.d0

  ! b1 - b1h
  
  real, parameter :: e1 = -3.d0/1280.d0
  real, parameter :: e2 = 0.d0
  real, parameter :: e3 = 6561.d0/632320.d0
  real, parameter :: e4 = -343.d0/20800.d0
  real, parameter :: e5 = 243.d0/12800.d0
  real, parameter :: e6 = -1.d0/95.d0
  real, parameter :: e7 = 0.d0

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
!$omp parallel workshare
  h0_old(:,:) = h_x(:,:)
!$omp end parallel workshare
  call timer_lib_out('str_mem')

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
        
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic=1,nc
              h0_x(ic,iv_loc) = h0_old(ic,iv_loc)
              h_x(ic,iv_loc) = h0_old(ic,iv_loc)
           enddo
        enddo
     else
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic=1,nc
              h0_x(ic,iv_loc) = h_x(ic,iv_loc)
           enddo
        enddo        
     endif
     call timer_lib_out('str')

     call cgyro_field_c   
     call cgyro_rhs(1)

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic,iv_loc) = h0_x(ic,iv_loc) &
                + a21*deltah2*rhs(ic,iv_loc,1)
        enddo
     enddo
     call timer_lib_out('str')

     call cgyro_field_c
     call cgyro_rhs(2)

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic,iv_loc) = h0_x(ic,iv_loc) &
                + deltah2*(a31*rhs(ic,iv_loc,1) &
                + a32*rhs(ic,iv_loc,2))
        enddo
     enddo
     call timer_lib_out('str')

     call cgyro_field_c
     call cgyro_rhs(3)
     
     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic,iv_loc) = h0_x(ic,iv_loc) &
                + deltah2*(a41*rhs(ic,iv_loc,1) &
                + a42*rhs(ic,iv_loc,2) &
                + a43*rhs(ic,iv_loc,3))
        enddo
     enddo
     call timer_lib_out('str')

     call cgyro_field_c     
     call cgyro_rhs(4)
     
     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic,iv_loc) = h0_x(ic,iv_loc)  &
                + deltah2*(a51*rhs(ic,iv_loc,1) &
                + a52*rhs(ic,iv_loc,2) &
                + a53*rhs(ic,iv_loc,3) &
                + a54*rhs(ic,iv_loc,4))
        enddo
     enddo
     call timer_lib_out('str')

     call cgyro_field_c
     call cgyro_rhs(5)

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic,iv_loc) = h0_x(ic,iv_loc) &
                + deltah2*(a61*rhs(ic,iv_loc,1) &
                + a62*rhs(ic,iv_loc,2) &
                + a63*rhs(ic,iv_loc,3) &
                + a64*rhs(ic,iv_loc,4) &
                + a65*rhs(ic,iv_loc,5))
        enddo
     enddo
     call timer_lib_out('str')

     call cgyro_field_c
     call cgyro_rhs(6)

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic,iv_loc) = h0_x(ic,iv_loc) + &
                deltah2*(a71*rhs(ic, iv_loc,1) &
                + a72*rhs(ic,iv_loc,2) &
                + a73*rhs(ic,iv_loc,3) &
                + a74*rhs(ic,iv_loc,4) &
                + a75*rhs(ic,iv_loc,5) &
                + a76*rhs(ic,iv_loc,6))
        enddo
     enddo
     call timer_lib_out('str')

     call cgyro_field_c
     call cgyro_rhs(7)

     !---------
     ! SOLUTION
     !---------

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic,iv_loc) = h0_x(ic,iv_loc) &
                + deltah2*(b1*rhs(ic,iv_loc,1) &
                + b3*rhs(ic,iv_loc,3) &
                + b4*rhs(ic,iv_loc,4) &
                + b5*rhs(ic,iv_loc,5) &
                + b6*rhs(ic,iv_loc,6) &
                + b7*rhs(ic,iv_loc,7))
        enddo
     enddo
     call timer_lib_out('str')
    
     !---------
     ! ERROR
     !---------
              
     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           rhs(ic,iv_loc,1) = deltah2*(e1*rhs(ic,iv_loc,1) &
                + e3*rhs(ic,iv_loc,3) &
                + e4*rhs(ic,iv_loc,4) &
                + e5*rhs(ic,iv_loc,5) &
                + e6*rhs(ic,iv_loc,6))
        enddo
     enddo

     error_x(1) = sum(abs(rhs(:,:,1)))
     error_x(2) = sum(abs(h_x))     
     call timer_lib_out('str')

     call timer_lib_in('str_comm')
     call MPI_ALLREDUCE(error_x,error_sum,2,MPI_DOUBLE_PRECISION,&
          MPI_SUM,CGYRO_COMM_WORLD, i_err)
     call timer_lib_out('str_comm')

     error_x = error_sum
     delta_x = error_x(1)+eps
     rel_error = error_x(1)/(error_x(2)+eps)
     var_error = sqrt(total_local_error+rel_error*rel_error)
     
     if (var_error < tol) then

        call cgyro_field_c

        conv = 1
        delta_t_tot = delta_t_tot+deltah2
        total_local_error = total_local_error + rel_error*rel_error
       
        scale_x = max((tol/delta_x*1.0/delta_t)**0.2, &
             (tol/delta_x*1.0/delta_t)**0.25)

        deltah2 = deltah2*max(1.0,min(6.0,scale_x))
        local_max_error = max(local_max_error,rel_error)

        call timer_lib_in('str_mem')
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic=1,nc
              h0_old(ic,iv_loc) = h0_x(ic, iv_loc)
           enddo
        enddo
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
 
  delta_t_gk = max(delta_t_last,4.0/5.0*deltah2)
  
  if (delta_t_last_step == 0.0) delta_t_last_step = delta_t_last

  if ( delta_t_last_step < 0.1*delta_t_gk ) then
     delta_t_gk = delta_t_last+delta_t_last_step
  else
     if (delta_t_last_step/itrk < 0.1*delta_t_gk) then
        delta_t_gk = delta_t_gk + delta_t_last_step/itrk
     endif
  endif

  delta_t_gk = min(delta_t,delta_t_gk)
  total_local_error = var_error

end subroutine cgyro_step_gk_bs5
