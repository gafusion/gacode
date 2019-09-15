subroutine cgyro_step_gk_ss

  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

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

  ! RK5(4) Sharp-Smart

  integer converged, conv, rk_MAX , iiter
  
  real orig_delta_x_t
  real total_delta_step
  real error_sum(2), delta_x_min, delta_x_max
  real delta_x, tau, deltah2, delta_t_last, rel_error, var_error
  real local_max_error, delta_t_last_step
  real deltah2_min, deltah2_max
  real error_x(2), tol
  real scale_x, delta2_h_min, delta2_h_max
  
  complex, dimension(:,:), allocatable :: h0_old

  real, parameter :: c2 = 16./105.
  real, parameter :: c3 = 8./35.
  real, parameter :: c4 = 9./20.
  real, parameter :: c5 = 2./3.
  real, parameter :: c6 = 7./9.
  real, parameter :: c7 = 1.0

  real, parameter :: a21  = 16./105.
  real, parameter :: a31  = 2./35.
  real, parameter :: a32  = 6./35.
  real, parameter :: a41  = 8793./40960.
  real, parameter :: a42  = -5103./8192.
  real, parameter :: a43 = 17577./20480.
  real, parameter :: a51 = 347./1458.
  real, parameter :: a52 = -7./20.
  real, parameter :: a53 = 3395./10044.
  real, parameter :: a54 = 49792./112995.
  
  real, parameter :: a61 = -1223224109959./9199771214400.
  real, parameter :: a62 = 1234787701./2523942720.
  real, parameter :: a63 = 568994101921./3168810084960.
  real, parameter :: a64  = -105209683888./891227836395.
  real, parameter :: a65  = 9./25.
  
  real, parameter :: a71  = 2462504862877./8306031988800.
  real, parameter :: a72  = -123991./287040.
  real, parameter :: a73  = 106522578491./408709510560.
  real, parameter :: a74  = 590616498832./804646848915.
  real, parameter :: a75  = -319138726./534081275.
  real, parameter :: a76  = 52758./71449.
  
  real, parameter :: b1 = 1093./15120.
  real, parameter :: b2 = 0.
  real, parameter :: b3 = 60025./190992.
  real, parameter :: b4 = 3200./20709.
  real, parameter :: b5 = 1611./11960.
  real, parameter :: b6 = 712233./2857960.
  real, parameter :: b7 = 3./40.

  real, parameter :: b1p = 84018211./991368000.
  real, parameter :: b2p  = 0.0
  real, parameter :: b3p  = 92098979./357791680.
  real, parameter :: b4p  = 17606944./67891005.
  real, parameter :: b5p  = 3142101./235253200.
  real, parameter :: b6p  = 22004596809./70270091500.
  real, parameter :: b7p  = 9./125.

  real, parameter :: EPS  = 2.2e-12
  
  allocate(h0_old(nc,nv_loc))
  
  call timer_lib_in('str')
  
  !! call cgyro_step_gkssp

  orig_delta_x_t = delta_t

  iiter = 0
  total_delta_step = 0.
  total_local_error = 0.
  local_max_error = 0.
  var_error = 0.
  
  delta_t_last = delta_t
  delta2_h_min = 1.0
  delta2_h_max = -1.0


  !! tol = .001d0
  !! if ( orig_delta_x_t .gt. 1 ) tol = tol/orig_delta_x_t
  !! if ( orig_delta_x_t .gt. 1 ) tol = tol*orig_delta_x_t

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
  
  rk_MAX = 1000

!$omp parallel workshare
  h0_old = h_x
!$omp end parallel workshare  
  conv = 0
  delta_t_gk = 0.0
  deltah2_min = 1.e10
  deltah2_max = -1.

  do while (total_delta_step .lt. orig_delta_x_t )
    
     if ( total_delta_step + deltah2 .gt. orig_delta_x_t ) then
        deltah2 = orig_delta_x_t - total_delta_step
        delta_t_last_step = deltah2
     else
        delta_t_gk = deltah2+delta_t_gk
        delta_t_last = deltah2
        !! deltah2_min = min(deltah2, deltah2_min)
        !! deltah2_max = max(deltah2, deltah2_max)
     endif

     if (( conv .eq. 0 ) .and. (iiter .ge. 1)) then
        !!
        !! last converged state
        !!
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

     !! for paper
     if ( i_proc == 0 ) write(*,*) i_proc, " step ", iiter, " ss deltah2 ", deltah2
     !! 

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

     !!     
     !!     h_x = h0_x + 3./32.*delta_x_t*(rhs(:,:,1) + 3.*rhs(:,:,2))
     !!

!$omp parallel workshare
     h_x = h0_x + deltah2*(a31*rhs(:,:,1) + a32*rhs(:,:,2))
!$omp end parallel workshare

     call cgyro_field_c
     
     ! Stage 3
     call cgyro_rhs(3)
     ! k4

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
     
     !  stage 5

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

!$omp parallel workshare
     h_x = h0_x + deltah2*(b1*rhs(:,:,1) &
          +b3*rhs(:,:,3) + b4*rhs(:,:,4) &
          + b5*rhs(:,:,5) + b6*rhs(:,:,6) + b7*rhs(:,:,7))
!$omp end parallel workshare

     call cgyro_field_c
     
!$omp parallel workshare
     rhs(:,:,1) = deltah2*((b1-b1p)*rhs(:,:,1) &
          + (b3-b3p)*rhs(:,:,3) &
          + (b4-b4p)*rhs(:,:,4) &
          + (b5-b5p)*rhs(:,:,5) &
          + (b6-b6p)*rhs(:,:,6) &
          + (b7-b7p)*rhs(:,:,7))
!$omp end parallel workshare

     error_sum = 0.
     error_x(1) = sum(abs(rhs(:,:,1)))
     error_x(2) = sum(abs(h_x))     

     call MPI_ALLREDUCE(error_x, error_sum, 2, MPI_DOUBLE_PRECISION,&
          MPI_SUM, MPI_COMM_WORLD, i_err)
     
     error_x(1) = error_sum(1)! local truncation error
     error_x(2) = error_sum(2)! local h0

     delta_x = error_x(1)
     tau = tol * error_x(2)
     
     !!     tau = tol * max(error_x(2), 1.d0)
     
     rel_error = error_x(1)/(error_x(2)+EPS)
     var_error = sqrt(total_local_error + rel_error*rel_error)

     !! if ( error_mode .eq. 0 ) then local error
     !! if ( error_mode .eq. 1 ) then variance error
     
     if ( var_error .lt. tol ) then
        
        if ( i_proc == 0 ) write(*,*) deltah2, " variance error mode ", var_error, total_local_error

!!     if ( error_x(1) .lt. tau ) then
        
!!         if ( i_proc == 0 ) &
!!             write(*,*) " local error mode ", rel_error, " variance error", var_error

!$omp parallel workshare        
        h0_old = h0_x
!$omp end parallel workshare
        
        converged = converged + 1
        conv = 1
        total_delta_step = total_delta_step + deltah2
        total_local_error = total_local_error + rel_error*rel_error

        scale_x = .95*max((tol/(delta_x + EPS)*1./delta_t)**(.2), &
             (tol/(delta_x + EPS)*1./delta_t)**(.25))
        
        scale_x = min(5., scale_x)
        deltah2 = deltah2*max(1., scale_x)
        local_max_error = max(local_max_error, rel_error)
     else
        conv = 0
        deltah2 = .5*deltah2          !! interpolate?
        if (i_proc .eq. 0 ) then
           write(*,*) " ss ***  error backing up *** not converged ", &
                " total delta_x step ", total_delta_step, " delta_x2 ", deltah2
           flush(6)
        endif
     endif
     
     deltah2 = min(deltah2, delta_x_max)
     deltah2 = max(delta_x_min, deltah2)

     iiter = iiter + 1

!!     if (deltah2 .lt. 1.d-8) then
!!        if ( i_proc .eq. 0 ) &
!!             write(*,*) " ******* Stopping due to small substep size ", deltah2
!!        flush(6)
!!        stop
!!     endif
     
     if ( iiter .gt. rk_MAX ) then
        if ( i_proc == 0 ) &
             write(*,*) " rk ss max count  exceeded ", iiter
        flush(6)
        stop
     endif
     
  enddo

  call timer_lib_out('str')  
  call cgyro_field_c

  delta_t_gk = delta_t_last
  if ( delta_t_last_step .lt.  1.e-4*delta_t_last )  & 
       delta_t_gk = delta_t_last + delta_t_last_step


  !! paper texts
  if ( i_proc .eq. 0 ) then
     write(*,*) " local error ", total_local_error
     write(*,*) " delta_t_gk ", delta_t_gk
     write(*,*) " variance local error sqrt of local error ", var_error
  endif
  
  total_local_error = var_error

  !! paper?
  !! if ( i_proc == 0 ) &
  !!       write(*,*) i_proc , " ss converged deltah2_min, max ", deltah2_min, deltah2_max
  
  call cgyro_filter
  
  if(allocated(h0_old)) deallocate(h0_old)  
  
end subroutine cgyro_step_gk_ss
