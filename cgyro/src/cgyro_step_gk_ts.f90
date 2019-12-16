subroutine cgyro_step_gk_ts

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none


  ! RK5(4) TS.x
  !
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
  double precision orig_delta_x_t
  double precision total_delta_step
  double precision error_sum(2), delta_x_min, delta_x_max
  double precision delta_x, tau, deltah2, delta_t_last, rel_error, var_error
  double precision local_max_error, delta_t_last_step
  double precision deltah2_min, deltah2_max
  double precision error_x(2), tol
  double precision scale_x
  
  complex, dimension(:,:), allocatable :: h0_old

  ! Butcher table

  double precision, parameter :: a21  = 25./119.
  double precision, parameter :: a31  = 3388./46225. 
  double precision, parameter :: a32  = 11662./46225.
  double precision, parameter :: a41  = 11108483733./100797037600.
  double precision, parameter :: a42  = 58559679./1028541200.
  double precision, parameter :: a43 =  1224328293./4031881504.
  double precision, parameter :: a51 = 3647841587291./23037540439200.
  double precision, parameter :: a52 = -6696949./44138800.
  double precision, parameter :: a53 = 6056670084509./12222397423712.
  double precision, parameter :: a54 = 398862271280./7038791373477.
  double precision, parameter :: a61 = -433777768769./734415679200.
  double precision, parameter :: a62 = 22837171./21402800.
  double precision, parameter :: a63 = 985520681359./241137642592.
  double precision, parameter :: a64  = -710571347975./66787651428.
  double precision, parameter :: a65  = 402295920045./56854183892.
  double precision, parameter :: a71  = 72623./1029420.
  double precision, parameter :: a72  = 0.
  double precision, parameter :: a73  = 7688883449./7106433180.
  double precision, parameter :: a74  = -204671984741./98413084125.
  double precision, parameter :: a75  = 542134811./298843875.
  double precision, parameter :: a76  = 107014./946125.
  
  double precision, parameter :: c2   = 25./119.
  double precision, parameter :: c3   = 14./43.
  double precision, parameter :: c4   = 129./274.
  double precision, parameter :: c5   = 19./34.
  double precision, parameter :: c6   = 1.0
  double precision, parameter :: c7   = 1.0

  double precision, parameter :: b1 = 72623./1029420.
  double precision, parameter :: b2 = 0.
  double precision, parameter :: b3 = 7688883449./7106433180.
  double precision, parameter :: b4 = -204671984741./98413084125.
  double precision, parameter :: b5 = 542134811./298843875.
  double precision, parameter :: b6 = 107014./946125.
  double precision, parameter :: b7 = 0.0d0

  double precision, parameter :: b1p = 75385109./749189000.
  double precision, parameter :: b2p  = 0.
  double precision, parameter :: b3p  = 34239333279133./46547137329000.
  double precision, parameter :: b4p  = -273056276943569./214868567006250.
  double precision, parameter :: b5p  = 283855408983./217491931250.
  double precision, parameter :: b6p  = 82079738./1032853125.
  double precision, parameter :: b7p  = 1./20.

  double precision, parameter :: EPS  = 2.2d-12  !! 2.e-16 triggers ieee error


  allocate(h0_old(nc,nv_loc))
  
  call timer_lib_in('str')
  
  orig_delta_x_t = delta_t

  iiter = 0
  total_delta_step = 0.
  total_local_error = 0.
  local_max_error = 0.

  tol = delta_t_tol

  deltah2 = delta_t_gk
  delta_t_last = delta_t
  
  if ( delta_t_gk .lt. 1.d-10) then
     deltah2 = orig_delta_x_t
     delta_t_gk = deltah2
  endif
  
  delta_x_min = 1.e-10*orig_delta_x_t
  delta_x_max = orig_delta_x_t

  converged = 0
  conv = 0
  
  rk_MAX = 10000

!$omp parallel workshare
  h0_old = h_x
!$omp end parallel workshare

  
  conv = 0
  delta_t_gk = 0.
  deltah2_min = 1.d10
  deltah2_max = 0.

  do while (total_delta_step .lt. orig_delta_x_t )
    
     if ( total_delta_step + deltah2 .gt. orig_delta_x_t ) then
        deltah2 = orig_delta_x_t - total_delta_step
        delta_t_last_step = deltah2
        !! if ( deltah2 .lt. 1.d-9 ) goto 1111   !! abandon
     else
        delta_t_gk = deltah2+delta_t_gk
        delta_t_last = deltah2
        deltah2_min = min(deltah2, deltah2_min)
        deltah2_max = max(deltah2, deltah2_max)
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

     !! paper line
     !! if ( i_proc == 0 ) write(*,*) i_proc, " step ", iiter, " ts deltah2 ", deltah2

     !
     ! Stage 1
     !
     
     call cgyro_rhs(1)
!$omp parallel workshare
     h_x = h0_x + deltah2 * a21 * rhs(:,:,1)
!$omp end parallel workshare

     call cgyro_field_c
     
     ! Stage 2 ! k2

     call cgyro_rhs(2)

!$omp parallel workshare
     h_x = h0_x + deltah2*(a31*rhs(:,:,1) + a32*rhs(:,:,2))
!$omp end parallel workshare

     call cgyro_field_c
     call cgyro_rhs(3)

!$omp parallel workshare
     h_x = h0_x + deltah2*(a41*rhs(:,:,1) &
          + a42*rhs(:,:,2) + a43*rhs(:,:,3))
!$omp end parallel workshare
     !
     call cgyro_field_c
     
     ! Stage 4
     call cgyro_rhs(4)
     !  k5

!$omp parallel workshare     
     h_x = h0_x+deltah2*(a51*rhs(:,:,1) + a52*rhs(:,:,2) &
          + a53*rhs(:,:,3) + a54*rhs(:,:,4))
!$omp end parallel workshare

     call cgyro_field_c
     
     call cgyro_rhs(5)

!$omp parallel workshare
     h_x = h0_x + deltah2*(a61*rhs(:,:,1) + a62*rhs(:,:,2) &
          + a63*rhs(:,:,3) + a64*rhs(:,:,4) + a65*rhs(:,:,5))
!$omp end parallel workshare

     call cgyro_field_c
     call cgyro_rhs(6)

!$omp parallel workshare
     h_x = h0_x + deltah2*(a71*rhs(:,:,1) + a72*rhs(:,:,2) &
          + a73*rhs(:,:,3) + a74*rhs(:,:,4) &
          + a75*rhs(:,:,5) + a76*rhs(:,:,6))
!$omp end parallel workshare

     
     call cgyro_field_c
     call cgyro_rhs(7)
     
     !! soln = h_x of order 5
     !! error_x(2) = sum(abs(h0_x))

!$omp parallel workshare
     h_x = h0_x + deltah2*(b1*rhs(:,:,1) &
          +b3*rhs(:,:,3) + b4*rhs(:,:,4) &
          + b5*rhs(:,:,5) + b6*rhs(:,:,6))
!$omp end parallel workshare

     call cgyro_field_c

!$omp parallel workshare
     rhs(:,:,1) = deltah2*((b1-b1p)*rhs(:,:,1) &
          + (b3-b3p)*rhs(:,:,3) &
          + (b4-b4p)*rhs(:,:,4) &
          + (b5-b5p)*rhs(:,:,5) + (b6-b6p)*rhs(:,:,6) &
          + (b7-b7p)*rhs(:,:,7))
!$omp end parallel workshare

     error_x = 0.
     error_sum = 0.     
     error_x(1) = sum(abs(rhs(:,:,1)))
     error_x(2) = sum(abs(h_x))     

     call MPI_ALLREDUCE(error_x, error_sum, 2, MPI_DOUBLE_PRECISION,&
          MPI_SUM, MPI_COMM_WORLD, i_err)
     
     error_x(1) = error_sum(1)
     error_x(2) = error_sum(2)

     delta_x = error_x(1)
     tau = tol * max(error_x(2), 1.)

     rel_error = error_x(1)/(error_x(2)+1.d-12)
     var_error = sqrt(total_local_error + rel_error*rel_error)

     !!
     !! need to local mode v.s. global mode of convergence check
     !!
     !! if ( error_mode .eq. 0 ) then local error
     !! if ( delta_x .lt. tau ) then
       !! if ( error_mode .eq. 1 ) then variance error
     
     if ( var_error .lt. tol ) then
        !! if ( i_proc == 0 ) write(*,*) " variance error mode ", var_error, total_local_error
        !! if ( error_x(1) .lt. tau ) then

        if ( i_proc == 0 ) &
             write(*,*) iiter, " delt variance error mode ", deltah2, var_error, rel_error

!$omp parallel workshare
        h0_old = h0_x
!$omp end parallel workshare

        converged = converged + 1
        conv = 1
        total_delta_step = total_delta_step + deltah2
        total_local_error = total_local_error + rel_error*rel_error

        !! scale_x = max(0.95*(tol/(delta_x + 1.e-12))**(.2), &
        !!     0.95d0*(tol/(delta_x + 1.e-12))**(.25))

        scale_x = max((tol/(delta_x + 1.e-12)*1./delta_t)**(.2), &
             (tol/(delta_x + 1.e-12)*1./delta_t)**(.25))

        scale_x = min(5., scale_x)
        deltah2 = deltah2*max(1., scale_x)
        
        local_max_error = max(local_max_error, rel_error)
     else
        conv = 0
        deltah2 = .5*deltah2          !! interpolate?
        if (i_proc .eq. 0 ) then
           write(*,*) " ts ***  backing up *** not converged ", &
                " new deltah2 ", deltah2, " rel error ", rel_error
        endif
     endif
     
     deltah2 = min(deltah2, delta_x_max)
     deltah2 = max(delta_x_min, deltah2)

     iiter = iiter + 1

     if ( iiter .gt. rk_MAX ) then
        if ( i_proc == 0 ) then
           write(*,*) " rk ts max count  exceeded ", iiter
        endif
        stop
     endif
  enddo
  call timer_lib_out('str')  
  call cgyro_field_c

  delta_t_gk = delta_t_last
  if ( delta_t_last_step .lt.  1.e-4*delta_t_last )  & 
       delta_t_gk = delta_t_last + delta_t_last_step


  ! used for paper output
!!  if ( i_proc .eq. 0 ) then
!!     write(*,*) " local error ", total_local_error
!!     write(*,*) " delta_t_gk ", delta_t_gk
!!     write(*,*) " variance local error sqrt of local error ", var_error
!!  endif
  
  total_local_error = var_error
  
  if ( i_proc == 0 ) &
       write(*,*) i_proc , " ts deltah2_min, max ", deltah2_min, deltah2_max, converged

  ! Filter special spectral components
  call cgyro_filter
    
  if(allocated(h0_old)) deallocate(h0_old)  
  
end subroutine cgyro_step_gk_ts
