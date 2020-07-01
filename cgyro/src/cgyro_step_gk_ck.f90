! Cash-Karp 6:5(4) adaptive integrator  |  multithreaded version

subroutine cgyro_step_gk_ck

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
  
  tol = error_tol

  itrk = 0
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
!$omp parallel workshare
  h0_old(:,:) = h_x(:,:)
!$omp end parallel workshare
 call timer_lib_out('str_mem')
        
  do while (delta_t_tot < delta_t)
     
     call timer_lib_in('str')
     if (delta_t_tot + deltah2 > delta_t) then
        deltah2 = delta_t-delta_t_tot
        delta_t_last_step = deltah2
     else
        delta_t_last = deltah2
        deltah2_min = min(deltah2,deltah2_min)
        deltah2_max = max(deltah2,deltah2_max)
     endif

     if (conv == 1) then
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic=1,nc
              h0_x(ic,iv_loc) = h_x(ic,iv_loc)
           enddo
        enddo
     else
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic=1,nc
              h0_x(ic,iv_loc) = h0_old(ic,iv_loc)
              h_x(ic,iv_loc) = h0_old(ic,iv_loc)
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
           h_x(ic, iv_loc) = h0_x(ic, iv_loc) &
                + 0.2d0*deltah2*rhs(ic, iv_loc, 1)
        enddo
     enddo
     call timer_lib_out('str')
     
     call cgyro_field_c
     call cgyro_rhs(2)

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic, iv_loc) = h0_x(ic,iv_loc) &
                + 1.0/40.0*deltah2*(3.0*rhs(ic, iv_loc, 1) &
                + 9.0*rhs(ic, iv_loc, 2))
        enddo
     enddo
     call timer_lib_out('str')

     call cgyro_field_c
     call cgyro_rhs(3)

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic, iv_loc) = h0_x(ic,iv_loc) &
                + deltah2*( 3.d0/10.d0*rhs(ic, iv_loc, 1) &
                - 9.d0/10.d0*rhs(ic, iv_loc, 2) &
                + 6.d0/5.d0*rhs(ic, iv_loc, 3))
        enddo
     enddo
     call timer_lib_out('str')
     
     call cgyro_field_c
     call cgyro_rhs(4) 

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic, iv_loc) = h0_x(ic,iv_loc) &
                + deltah2*(-11.d0/54.d0 *rhs(ic, iv_loc, 1) &
                + 5.d0/2.d0*rhs(ic, iv_loc, 2) &
                - 70.d0/27.d0*rhs(ic, iv_loc, 3) &
                + 35.d0/27.d0*rhs(ic, iv_loc, 4))
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
                + deltah2*(1631.d0/55296.d0*rhs(ic, iv_loc, 1) &
                + 175.d0/512.d0*rhs(ic, iv_loc, 2) &
                + 575.d0/13824.d0*rhs(ic, iv_loc, 3) &
                + 44275.d0/110592.d0*rhs(ic, iv_loc, 4) &
                + 253.d0/4096.d0*rhs(ic, iv_loc, 5))
        enddo
     enddo
     call timer_lib_out('str')
     
     call cgyro_field_c
     call cgyro_rhs(6)

     !---------
     ! SOLUTION
     !---------

     call timer_lib_in('str')
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic=1,nc
           h_x(ic, iv_loc) = h0_x(ic,iv_loc) &
                + deltah2*( 37.d0/378.d0*rhs(ic, iv_loc, 1) &
                + 250.d0/621.d0 * rhs(ic, iv_loc, 3) &
                + 125.d0/594.d0*rhs(ic, iv_loc, 4) &
                + 512.d0/1771.d0*rhs(ic, iv_loc, 6))
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
           rhs(ic,iv_loc,1) = deltah2*(& 
                (37.d0/378.d0-2825.d0/27648.d0)*rhs(ic, iv_loc, 1) &
                + (250.d0/621.d0-18575.d0/48384.d0)*rhs(ic, iv_loc, 3) &
                + (125.d0/594.d0-13525.d0/55296.d0)*rhs(ic, iv_loc, 4) &
                + (-277.d0/14336.d0)*rhs(ic, iv_loc, 5) &
                + (512.d0/1771.d0-1.d0/4.d0)*rhs(ic, iv_loc, 6))
        enddo
     enddo

     error_x(1) = sum(abs(rhs(:,:,1)))
     error_x(2) = sum(abs(h_x))
     call timer_lib_out('str')

     call timer_lib_in('str_comm')
     call MPI_ALLREDUCE(error_x,error_sum,2,MPI_DOUBLE_PRECISION,&
          MPI_SUM,CGYRO_COMM_WORLD,i_err)
     call timer_lib_out('str_comm')

     error_x = error_sum
     delta_x = error_x(1)+eps

     rel_error = error_x(1)/(error_x(2)+eps)
     var_error = sqrt(total_local_error+rel_error*rel_error)
     
     if (var_error < tol) then
        call cgyro_field_c

        delta_t_tot = delta_t_tot+deltah2
        total_local_error = total_local_error+rel_error*rel_error

        scale_x = max((tol/delta_x*1.0/delta_t)**0.2,&
             (tol/delta_x*1.0/delta_t)**0.25)

        deltah2 = deltah2*max(1.0,min(6.0,scale_x))
        
        local_max_error = max(local_max_error,rel_error)
        conv = 1

        call timer_lib_in('str_mem')
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic=1,nc
              h0_old(ic,iv_loc) = h0_x(ic,iv_loc)
           enddo
        enddo
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

  delta_t_gk = max(delta_t_last,deltah2)

  if (delta_t_last_step == 0.0) delta_t_last_step = delta_t_last

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
