subroutine cgyro_step_gk_ck

  use timer_lib
  use mpi
  use cgyro_globals

  implicit none

  ! RK4(5) time-advance for the distribution ; cash-karp method of embedde
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

  !! cash_karp rk45 ... with changes in discontinuity or oscillatory behavior
  !! not not best for stiff
  !!

  real deltah2, orig_delta_t
  real delta_x_min, delta_x_max, local_max_error
  real tol_ck, total_delta_step
  real error_x(5), error_sum(5)
  real tau_ck, delta_t_old, delta_t_gk_old, delta_t_last, delta_t_last_step
  real last_total_error, rel_error, var_error
  real deltah2_min, deltah2_max, scale_x
  
  integer conv, iiter, rk_count, rkMAXITER, converged, err_x
  

  complex, dimension(:,:), allocatable :: rk_error
  complex, dimension(:,:), allocatable :: h0_old

  real, parameter :: EPS = 2.e-12

  call timer_lib_in('str_mem')
  allocate(h0_old(nc,nv_loc))
  call timer_lib_out('str_mem')

  call timer_lib_in('str')

  !

  orig_delta_t = delta_t       !! keep track of delta_t for exiting of subroutine  
  tol_ck = delta_t_tol
  delta_t_gk_old = delta_t_gk
  deltah2 = delta_t_gk

  if ( delta_t_gk .lt. 1.0d-10) then
     deltah2 = .1*orig_delta_t
     delta_t_gk = delta_t
  endif

  !! taken from last step

  delta_x_min = orig_delta_t*1.e-10
  delta_x_max = orig_delta_t
  delta_t_old = delta_t_gk

  total_delta_step = 0.
  total_local_error = 0.
  last_total_error = 0.
  var_error = 0.
  conv = 0
  iiter = 0
  rk_count = 0
  rkMAXITER = 1000
  converged = 0
  deltah2_min = 1.e10
  deltah2_max = -1.
  
!$omp parallel workshare
  h0_old = h_x
!$omp end parallel workshare

  
  delta_t_gk = 0.
  delta_t_last = deltah2
  local_max_error=0.

  do while (total_delta_step .lt. orig_delta_t )

     if ( total_delta_step + deltah2 .gt. orig_delta_t ) then
        deltah2 = orig_delta_t - total_delta_step
        delta_t_last_step = deltah2
     else
        delta_t_last = deltah2
        deltah2_min = min(deltah2, deltah2_min)
        deltah2_max = max(deltah2, deltah2_max)
        !! if ( deltah2 .lt. 1.d-8 ) goto 1111   !! abandon keep?
     endif

     if ( i_proc .eq. 0 ) then
        write(*,*) " paper ck4 current iiter, deltah2 ", iiter, deltah2
     endif

     if (( conv .eq. 0 ) .and. (iiter .ge. 1)) then

!$omp parallel workshare
        h0_x = h0_old
        h_x = h0_old
!$omp end parallel workshare
     else
!$omp parallel workshare        
        h0_x = h_x
!$omp end parallel workshare

     endif

     call cgyro_field_c

     ! Stage 1

     call cgyro_rhs(1)            !k1

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc, iv_loc) &
                + 0.2*deltah2*rhs(ic_loc, iv_loc, 1)
        enddo
     enddo
     
     call cgyro_field_c
     
     ! Stage 2 ! k2
     
     call cgyro_rhs(2)

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + 1./40.*deltah2*(3.*rhs(ic_loc, iv_loc, 1) &
                + 9.*rhs(ic_loc, iv_loc, 2))
        enddo
     enddo

     call cgyro_field_c
     
     ! Stage 3

     call cgyro_rhs(3) ! k3

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*( 3./10.*rhs(ic_loc, iv_loc, 1) &
                - 9./10.*rhs(ic_loc, iv_loc, 2) &
                + 6./5.*rhs(ic_loc, iv_loc, 3))
        enddo
     enddo
     
     call cgyro_field_c
     
     call cgyro_rhs(4) 

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*(-11./54. *rhs(ic_loc, iv_loc, 1) &
                + 5./2.*rhs(ic_loc, iv_loc, 2) &
                - 70./27.*rhs(ic_loc, iv_loc, 3) &
                + 35./27.*rhs(ic_loc, iv_loc, 4))
        enddo
     enddo
     

     call cgyro_field_c
     call cgyro_rhs(5)

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*(1631./55296.*rhs(ic_loc, iv_loc, 1) &
                + 175./512.*rhs(ic_loc, iv_loc, 2) &
                + 575./13824.*rhs(ic_loc, iv_loc, 3) &
                + 44275./110592.*rhs(ic_loc, iv_loc, 4) &
                + 253./4096.*rhs(ic_loc, iv_loc, 5))
        enddo
     enddo
     
     call cgyro_field_c
     call cgyro_rhs(6)
     
     !! soln of order 5

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*( 37./378.*rhs(ic_loc, iv_loc, 1) &
                + 250./621. * rhs(ic_loc, iv_loc, 3) &
                + 125./594.*rhs(ic_loc, iv_loc, 4) &
                + 512./1771.*rhs(ic_loc, iv_loc, 6))
        enddo
     enddo

     call cgyro_field_c

     !

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           rhs(ic_loc, iv_loc, 1) = deltah2*(& 
                (37./378.-2825./27648.)*rhs(ic_loc, iv_loc, 1) &
                + (250./621.-18575./48384.)*rhs(ic_loc, iv_loc, 3) &
                + (125./594.-13525./55296.)*rhs(ic_loc, iv_loc, 4) &
                + (-277./14336.)*rhs(ic_loc, iv_loc, 5) &
                + (512./1771.-1./4.)*rhs(ic_loc, iv_loc, 6))
        enddo
     enddo

     error_sum = 0.
     error_x(1) = sum(abs(rhs(:, :, 1)))
     error_x(2) = sum(abs(h_x))

     call MPI_ALLREDUCE(error_x, error_sum, 2, &
          MPI_DOUBLE_PRECISION,&
          MPI_SUM, MPI_COMM_WORLD, err_x)

     error_x = error_sum

     tau_ck = tol_ck*max(error_x(2), 1.0)

     rel_error = error_x(1)/(error_x(2) +1.d-12)

     !!     var_error = sqrt(total_local_error + rel_error*rel_error)/max(sqrt(abs(iiter-1.)), 1.)
     var_error = sqrt(total_local_error + rel_error*rel_error)

     !! error_mode = 0 ; local error
     !! fro local error if ( error_x(1) .lt. tau_ck ) then
     !! for variance tight control if ( var_error .lt. tol_ck ) then
     
     if ( var_error .lt. tol_ck ) then
        
!$omp parallel workshare        
        h0_old = h_x
!$omp end parallel workshare

        if ( i_proc .eq. 0 ) then
           write(*,*) " paper ck4 current rel_error, var_error ", rel_error, var_error
        endif
        
        total_delta_step = total_delta_step + deltah2
        total_local_error = total_local_error + rel_error*rel_error
        
        scale_x = .95*max((tol_ck/(error_x(1) + EPS)*1./delta_t)**(.2), &
             (tol_ck/(error_x(1) + EPS)*1./delta_t)**(.25))

        scale_x = min(5., scale_x)
        deltah2 = deltah2*max(1., scale_x)
        
        local_max_error = max(local_max_error, rel_error)
        conv = 1
        converged=converged+1
     else
        deltah2 = .5*deltah2
        conv = 0
        if ( i_proc .eq. 0 ) then
           write(*,*) " *** backing up not converged step = ", rk_count
           write(*,*) " new deltah2 ", deltah2
           flush(6)
        endif
     endif

     deltah2 = min(deltah2, delta_x_max)
     deltah2 = max(delta_x_min, deltah2)

     iiter = iiter + 1
    
     if ( iiter .gt. rkMAXITER) then
        write(*,*) " **gk ck4(5) step exceeded max iteration count ", rk_count
        write(*,*) " **gk ck4(5) deltah2 ", deltah2
        stop
     endif

  end do

  call timer_lib_out('str')
  call cgyro_filter

  delta_t_gk = delta_t_last
  if ( delta_t_last_step .lt.  1.e-4*delta_t_last )  & 
       delta_t_gk = delta_t_last + delta_t_last_step

  total_local_error = var_error

  !! paper
  !! if ( i_proc .eq. 0 ) write(*,*) " total local error ", total_local_error,  &
  !! " variance local error ", var_error, &
  !! " delta_t_gk ", delta_t_gk

  !! if ( i_proc .eq. 0 ) write(*,*) " paper Out ck4 deltah2_min, max ", &
  !! iiter, deltah2_min, deltah2_max
  !!

  call timer_lib_out('str')

  if(allocated(h0_old)) deallocate(h0_old)
  
end subroutine cgyro_step_gk_ck
