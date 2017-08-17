!-----------------------------------------------------
! gyro_field_time_derivative.f90
!
! PURPOSE:
!  Compute partial time derivative of 
! 
!  1. h_cap (non-Doppler)
!  2. field_blend (Doppler)
!  3. gyro_uv (Doppler)
!-------------------------------------------------------

subroutine gyro_field_time_derivative

  use mpi
  use gyro_globals
  use math_constants
  use ompdata

  !---------------------------------------------------
  implicit none
  !
  real :: minus_n_omega
  !---------------------------------------------------

!$omp parallel private(minus_n_omega,is,i,m)
  do is=1,n_kinetic
     do i=ibeg,iend
        do m=1,n_stack

           h_cap(m,i,:,is) = h(m,i,:,is) &
                +z(is)*alpha_s(is,i)*gyro_u(m,i,:,is)

        enddo ! m
     enddo ! i
  enddo ! is

  do i=ibeg,iend

     ! Use 3-point rule for time derivative

     ! Note sign: minus_n_omega = -n*w0

     minus_n_omega = omega_eb_s(i)

     field_blend_dot(:,i,:) = (1.5*field_blend(:,i,:)-2.0*field_blend_old(:,i,:) &
          +0.5*field_blend_old2(:,i,:))/dt &
          +i_c*minus_n_omega*field_blend(:,i,:)

     gyro_uv_dot(:,i,:,:,:)=(1.5*gyro_uv(:,i,:,:,:)-2.0*gyro_uv_old(:,i,:,:,:) &
          +0.5*gyro_uv_old2(:,i,:,:,:))/dt &
          +i_c*minus_n_omega*gyro_uv(:,i,:,:,:)

     h_cap_dot(:,i,:,:)=(1.5*h_cap(:,i,:,:)-2.0*h_cap_old(:,i,:,:) &
          +0.5*h_cap_old2(:,i,:,:))/dt &
          +i_c*minus_n_omega*h_cap(:,i,:,:)

  enddo
!$omp end parallel

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_field_time_derivative called]'
  endif

end subroutine gyro_field_time_derivative
