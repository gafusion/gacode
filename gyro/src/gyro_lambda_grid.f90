!---------------------------------------------------------
! mgyro_lambda_grid.f90
!
! PURPOSE:
!  Lambda grid (and tau/theta spacing) setup.
!---------------------------------------------------------

subroutine gyro_lambda_grid

  use gyro_globals
  use math_constants

  !------------------------------------------
  implicit none
  !
  real :: s_tp
  !
  real, dimension(:), allocatable :: s_temp
  real, dimension(:), allocatable :: w_temp
  !------------------------------------------

  if (n_lambda == 1) then
     lambda(:,:) = 0.0
     w_lambda(:,:) = 0.5
     return
  endif

  ! Fast startup option for flat profiles

  do i=1,n_x

     if (i == 1 .or. flat_profile_flag == 0) then

        call gyro_to_geo(i)
        
        call gyro_banana_init(nint_ORB_s)
        
        ! Passing

        call gyro_banana_getlambda(lambda_tp(i),lambda_max(i))
       
        ! two regions

        call gyro_banana_integrate_tau(lambda_tp(i),s_tp)
        
        ! s_tp => s at trapped-passing boundary

        allocate(s_temp(n_pass))
        allocate(w_temp(n_pass))
        call gauss_legendre(0.0,s_tp,s_temp,w_temp,n_pass)
        s_lambda(1:n_pass) = s_temp
        w_lambda(i,1:n_pass) = w_temp
        deallocate(s_temp)
        deallocate(w_temp)

        allocate(s_temp(n_trap))
        allocate(w_temp(n_trap))
        call gauss_legendre(s_tp,1.0,s_temp,w_temp,n_trap)
        s_lambda(n_pass+1:n_lambda) = s_temp
        w_lambda(i,n_pass+1:n_lambda) = w_temp
        deallocate(s_temp)
        deallocate(w_temp)

        print *,'s2'
        do k=1,n_lambda
           call gyro_banana_s2lambda(s_lambda(k),lambda(i,k))
        enddo

        ! Two signs of velocity

        do k=1,n_lambda
           w_lambda(i,k) = 0.5*w_lambda(i,k)
        enddo

     else

        w_lambda(i,:) = w_lambda(1,:)
        lambda(i,:)   = lambda(1,:)
        lambda_max(i) = lambda_max(1)
        lambda_tp(i)  = lambda_tp(1)

     endif

  enddo ! i

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_lambda_grid done]'
  endif

end subroutine gyro_lambda_grid
