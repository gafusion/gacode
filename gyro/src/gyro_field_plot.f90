!---------------------------------------------------------
! gyro_field_plot.f90
!
! PURPOSE:
!  Generate interpolation of fields suitable for plotting.
!---------------------------------------------------------

subroutine gyro_field_plot

  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  integer :: j_plot
  !  
  complex, dimension(n_theta_plot,n_x) :: phi_tmp
  complex, dimension(n_blend,n_x) :: c0_tmp
  !---------------------------------------------------

  phi_plot = (0.0,0.0)

  !---------------------------------------------------------
  ! Phi, A_parallel, B_parallel
  !
  do ix=1,n_field

     c0_tmp(:,:) = field_blend(:,:,ix)

     do i=1,n_x
        do j_plot=1,n_theta_plot

           phi_tmp(j_plot,i) = sum(c0_tmp(:,i)*blend_plot(:,j_plot,i))

        enddo ! j_plot
     enddo ! i

     phi_plot(:,:,ix) = phi_tmp(:,:)

  enddo
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! E_parallel 
  !
  if (eparallel_plot_flag == 1) then

     do i=1,n_x
        do j_plot=1,n_theta_plot

           if (n_field > 1) then

              phi_tmp(j_plot,i) = -sum(field_blend_dot(:,i,2)*blend_plot(:,j_plot,i)) &
                   -sum(field_blend(:,i,1)*blend_prime_plot(:,j_plot,i)) &
                   /rmaj_s(i)/q_s(i)/g_theta_plot(i,j_plot)

           else

              phi_tmp(j_plot,i) = &
                   -sum(field_blend(:,i,1)*blend_prime_plot(:,j_plot,i)) &
                   /rmaj_s(i)/q_s(i)/g_theta_plot(i,j_plot)

           endif


        enddo ! j_plot
     enddo ! i

     phi_plot(:,:,n_field+1) = phi_tmp(:,:)

  endif
  !---------------------------------------------------------

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_field_plot done]' 
  endif

end subroutine gyro_field_plot
