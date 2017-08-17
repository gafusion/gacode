!----------------------------------------------------------
! gyro_timestep_error.f90
!
! PURPOSE:
!  This routine calculates the instantaeous timestep error.
!----------------------------------------------------------

subroutine gyro_timestep_error

  use mpi
  use gyro_globals
  use gyro_pointers
  use ompdata

  !---------------------------------------------------
  implicit none
  !
  real :: rk_error_loc(n_kinetic,2)
  real :: rk_error(n_kinetic,2)
  real :: tol = 1e-20
  !---------------------------------------------------
  rk_error_loc(:,:) = 0.0

!$omp parallel private(is,p_nek_loc,i,m) &
!$omp reduction(+:rk_error_loc)
  do is=1,n_kinetic
     do p_nek_loc=1,n_nek_loc_1
        do i = ibeg, iend
           do m=1,n_stack
              rk_error_loc(is,1) = rk_error_loc(is,1) + abs(h_err(m,i,p_nek_loc,is))
              rk_error_loc(is,2) = rk_error_loc(is,2) + abs(h(m,i,p_nek_loc,is))
           enddo
        enddo
     enddo
  enddo
!$omp end parallel

  call MPI_ALLREDUCE(rk_error_loc,&
       rk_error,&
       n_kinetic*2,&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       GYRO_COMM_WORLD,&
       i_err)

  if (step == 0 .or. sum(abs(rk_error(:,2))) < tol) then
     rk_error(:,1) = 0.0
     rk_error(:,2) = 1.0
  endif

  time_error(:) = rk_error(:,1)/rk_error(:,2)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_timestep_error done]'
  endif

end subroutine gyro_timestep_error
