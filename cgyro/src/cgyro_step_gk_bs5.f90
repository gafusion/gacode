subroutine cgyro_step_gk_bs5
  use mpi
  use timer_lib
  use cgyro_globals

  implicit none

  ! BS5(4), from Bogacki-Shampine, 1996
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
  double precision orig_delta_x_t, total_delta_x_step, delta_x_min, delta_x_max
  double precision error_sum(2), error_x(2)
  double precision delta_x, tau, deltah2, local_max_error
  double precision rel_error, delta_t_last
  double precision delta_t_last_step
  double precision deltah2_max, deltah2_min
  double precision var_error, scale_x, tol

  complex, dimension(:,:), allocatable :: h0_old

  !! butcher table

  double precision, parameter :: a21  =  1./6.
  double precision, parameter :: a31  = 2./27.
  double precision, parameter :: a32  = 4./27.
  
  double precision, parameter :: a41  = 183./1372.
  double precision, parameter :: a42  = -162./343.
  double precision, parameter :: a43 =  1053./1372.
  
  double precision, parameter :: a51 =  68./297.
  double precision, parameter :: a52 = -4./11.
  double precision, parameter :: a53 = 42./143.
  double precision, parameter :: a54 = 1960./3861.
  
  double precision, parameter :: a61   = 597./22528.
  double precision, parameter :: a62   = 81./352.
  double precision, parameter :: a63   = 63099./585728.
  double precision, parameter :: a64  = 58653./366080.
  double precision, parameter :: a65  = 4617./20480.
  
  double precision, parameter :: a71  = 174197./959244.
  double precision, parameter :: a72  = -30942./79937.
  double precision, parameter :: a73  = 8152137./19744439.
  double precision, parameter :: a74  = 666106./1039181.
  double precision, parameter :: a75  = -29421./29068.
  double precision, parameter :: a76  = 482048./414219.
  
  double precision, parameter :: a81  = 587./8064.
  double precision, parameter :: a82  = 0.
  double precision, parameter :: a83  = 4440339./15491840.
  double precision, parameter :: a84  = 24353./124800.
  double precision, parameter :: a85  = 387./44800.
  double precision, parameter :: a86  = 2152./5985.
  double precision, parameter :: a87  = 7267./94080.

  double precision, parameter :: c2   = 1./6.
  double precision, parameter :: c3   = 2./9.
  double precision, parameter :: c4   = 3./7.
  double precision, parameter :: c5   = 2./3.
  double precision, parameter :: c6   = 3./4.
  double precision, parameter :: c7   = 1.0
  double precision, parameter :: c8   = 1.0

  double precision, parameter :: b1 =587./8064. 
  double precision, parameter :: b2 =0.
  double precision, parameter :: b3 =4440339./15491840. 
  double precision, parameter :: b4 =24353./124800. 
  double precision, parameter :: b5 =387./44800. 
  double precision, parameter :: b6 =2152./5985. 
  double precision, parameter :: b7 =7267./94080. 
  double precision, parameter :: b8 = 0.

  double precision, parameter :: b1h =6059./80640. 
  double precision, parameter :: b2h =0.
  double precision, parameter :: b3h =8559189./30983680. 
  double precision, parameter :: b4h =26411./124800. 
  double precision, parameter :: b5h =-927./89600.
  double precision, parameter :: b6h =443./1197. 
  double precision, parameter :: b7h =7267./94080.
  double precision, parameter :: b8h =0.

  double precision, parameter :: b1p = 2479./34992.
  double precision, parameter :: b2p = 0.
  double precision, parameter :: b3p = 123./416.
  double precision, parameter :: b4p = 612941./3411720.
  double precision, parameter :: b5p = 43./1440.
  double precision, parameter :: b6p = 2272./6561.
  double precision, parameter :: b7p = 79937./1113912.
  double precision, parameter :: b8p = 3293./556956.

  ! b1 - b1h
  
  double precision, parameter :: e1 = -3./1280.
  double precision, parameter :: e2 = 0.
  double precision, parameter :: e3 = 6561./632320.
  double precision, parameter :: e4 = -343./20800.
  double precision, parameter :: e5 = 243./12800.
  double precision, parameter :: e6 = -1./95.
  double precision, parameter :: e7 = 0.

  double precision, parameter :: EPS  = 2.2d-12
  
  allocate(h0_old(nc,nv_loc))

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
  
  delta_x_min = 1.d-10*orig_delta_x_t
  delta_x_max = orig_delta_x_t

  converged = 0
  conv = 0
  total_local_error = 0.
  
  rk_MAX = 1000

!$omp parallel workshare
  h0_old = h_x
!$omp end parallel workshare
  conv = 0
  delta_t_gk = 0.
  deltah2_min = 1.d10
  deltah2_max = 0.

  do while (total_delta_x_step .lt. orig_delta_x_t .and. iiter .le. rk_MAX )
    
     if ( total_delta_x_step + deltah2 .gt. orig_delta_x_t ) then
        deltah2 = orig_delta_x_t - total_delta_x_step
        delta_t_last_step = deltah2
     else
        !! delta_t_gk = deltah2+delta_t_gk
        delta_t_last = deltah2
        deltah2_min = min(deltah2, deltah2_min)
        deltah2_max = max(deltah2, deltah2_max)
     endif
     
     if (( conv .eq. 0 ) .and. (iiter .ge. 1)) then
        ! not converged so backing up
        
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic_loc=1,nc
              h0_x(ic_loc,iv_loc) = h0_old(ic_loc, iv_loc)
              h_x(ic_loc,iv_loc) = h0_old(ic_loc, iv_loc)
           enddo
        enddo
     else
        ! h0_x = h_x
!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic_loc=1,nc
              h0_x(ic_loc,iv_loc) = h_x(ic_loc, iv_loc)
           enddo
        enddo
        
     endif

     call cgyro_field_c

     ! Stage 1
     !
     
     call cgyro_rhs(1)

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc,iv_loc) = h0_x(ic_loc,iv_loc) &
                + a21*deltah2*rhs(ic_loc,iv_loc,1)
        enddo
     enddo

     call cgyro_field_c

     ! Stage 2 ! k2

     call cgyro_rhs(2)

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*(a31*rhs(ic_loc,iv_loc,1) &
                + a32*rhs(ic_loc,iv_loc,2))
        enddo
     enddo

     call cgyro_field_c

     ! Stage 3
     
     call cgyro_rhs(3)
     
     ! k4
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc, iv_loc) &
                + deltah2*(a41*rhs(ic_loc, iv_loc, 1) &
                + a42*rhs(ic_loc,iv_loc,2) &
                + a43*rhs(ic_loc,iv_loc,3))
        enddo
     enddo


     call cgyro_field_c
     
     ! Stage 4
     
     call cgyro_rhs(4)
     
     !  k5

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc)  &
                + deltah2*(a51*rhs(ic_loc, iv_loc, 1) &
                + a52*rhs(ic_loc,iv_loc,2) &
                + a53*rhs(ic_loc,iv_loc,3) &
                + a54*rhs(ic_loc,iv_loc,4))
        enddo
     enddo

     
     call cgyro_field_c
     
     !  stage 5
     
     call cgyro_rhs(5)

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*(a61*rhs(ic_loc, iv_loc, 1) &
                + a62*rhs(ic_loc,iv_loc,2) &
                + a63*rhs(ic_loc,iv_loc,3) &
                + a64*rhs(ic_loc,iv_loc,4) &
                + a65*rhs(ic_loc,iv_loc,5))
        enddo
     enddo

     call cgyro_field_c

     ! stage 6
     call cgyro_rhs(6)

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) + &
                deltah2*(a71*rhs(ic_loc, iv_loc, 1) &
                + a72*rhs(ic_loc,iv_loc,2) &
                + a73*rhs(ic_loc,iv_loc,3) &
                + a74*rhs(ic_loc,iv_loc,4) &
                + a75*rhs(ic_loc,iv_loc,5) &
                + a76*rhs(ic_loc,iv_loc,6))
        enddo
     enddo

     call cgyro_field_c
     call cgyro_rhs(7)

!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           h_x(ic_loc, iv_loc) = h0_x(ic_loc,iv_loc) &
                + deltah2*(b1*rhs(ic_loc, iv_loc, 1) &
                + b3*rhs(ic_loc,iv_loc,3) &
                + b4*rhs(ic_loc,iv_loc,4) &
                + b5*rhs(ic_loc,iv_loc,5) &
                + b6*rhs(ic_loc,iv_loc,6) &
                + b7*rhs(ic_loc,iv_loc,7))
        enddo
     enddo
     
     call cgyro_field_c

     !! compute error between 5th and 4th order soln
          
!$omp parallel do collapse(2)
     do iv_loc=1,nv_loc
        do ic_loc=1,nc
           rhs(ic_loc,iv_loc,1)= deltah2*((b1-b1h)*rhs(ic_loc, iv_loc, 1) &
                + (b3-b3h)*rhs(ic_loc,iv_loc,3) &
                + (b4-b4h)*rhs(ic_loc,iv_loc,4) &
                + (b5-b5h)*rhs(ic_loc,iv_loc,5) &
                + (b6-b6h)*rhs(ic_loc,iv_loc,6))
        enddo
     enddo

     ! b*p gives another soln based on a8*

     error_sum = 0.
     error_x = 0.
     error_x(1) = sum(abs(rhs(:,:,1)))
     error_x(2) = sum(abs(h_x))     
     
     call MPI_ALLREDUCE(error_x, error_sum, 2, MPI_DOUBLE_PRECISION,&
          MPI_SUM, MPI_COMM_WORLD, i_err)
     
     error_x = error_sum
     delta_x = error_x(1)
     tau = tol*max(error_x(2), 1.)

     rel_error = error_x(1)/(error_x(2)+EPS)
     var_error = sqrt(total_local_error + rel_error*rel_error)

     !! method 1 local error
     
     !!if ( error_x(1) .lt. tau ) then
        
     !! method 2 "variance" error
     
     if ( var_error .lt. tol ) then
        
        !!        if ( i_proc == 0 ) &
        !! write(*,*) " paper BS5(4) **** local error mode ", rel_error, " variance error", var_err
        !! or

        ! paper
        if ( i_proc == 0 ) &
             write(*,*) " paper dt ", deltah2, " BS5(4) **** var error mode ", rel_error, " variance error", var_error

!$omp parallel do collapse(2)
        do iv_loc=1,nv_loc
           do ic_loc=1,nc
              h0_old(ic_loc,iv_loc) = h0_x(ic_loc, iv_loc)
           enddo
        enddo

        converged = converged + 1
        conv = 1
        
        total_delta_x_step = total_delta_x_step + deltah2
        total_local_error = total_local_error + rel_error*rel_error

        scale_x = 0.95*max((tol/(delta_x + EPS )*1./delta_t)**(.2), &
             (tol/(delta_x + EPS )*1./delta_t)**(.25))

        scale_x = min(5., scale_x)
        deltah2 = deltah2*max(1., scale_x)
        local_max_error = max(local_max_error, rel_error)
     else
        conv = 0
        deltah2 = .5*deltah2
        if (i_proc .eq. 0 ) then
           write(*,*) " bs5 ***  error backing up, not converged ", &
                " total delta_x step ", total_delta_x_step, " new deltah2 ", deltah2
        endif
        flush(6)
     endif
     
     deltah2 = min(deltah2, delta_x_max)
     deltah2 = max(delta_x_min, deltah2)

     iiter = iiter + 1

     if ( iiter .gt. rk_MAX ) then
        if ( i_proc .eq. 0 ) then
           write(*,*) " bs5  max count  exceeded, stopping ", iiter
           write(*,*) " bs5  step size ", deltah2
           flush(6)           
        endif
        stop
     endif
     
  enddo
  
  call timer_lib_out('str')
  call cgyro_field_c

  delta_t_gk = delta_t_last
  
  if ( delta_t_last_step .lt.  1.e-4*delta_t_last )  & 
       delta_t_gk = delta_t_last + delta_t_last_step

  total_local_error = var_error

  !!
  !! if ( i_proc == 0 ) &
  !! write(*,*) i_proc , " paper bst deltah2_min, max converged ", deltah2_min, deltah2_max
  !!  
  ! Filter special spectral components
  
  call cgyro_filter

  if(allocated(h0_old)) deallocate(h0_old)
  
end subroutine cgyro_step_gk_bs5
