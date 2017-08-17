!------------------------------------------------------------------------
! gyro_moments_plot.f90
!
! PURPOSE:
!  Make density, energy and V_parallel fluctuation moments (with 
!  same functional form as phi_plot) for diagnostic output file.
! 
!  This routine assumes phi_plot has been computed gyro_field_plot.
! 
!     density moment: moments_plot(:,:,:,1)
!      energy moment: moments_plot(:,:,:,2)
!  V_parallel moment: moments_plot(:,:,:,3)
!
! NOTES:
!  See http://fusion.gat.com/theory/Gyrousermanual for precise
!  definitions.
!--------------------------------------------------------------------

subroutine gyro_moments_plot

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants
  use ompdata

  !---------------------------------------------------
  implicit none
  !
  integer :: j_plot
  !
  complex, dimension(n_blend,n_x,3) :: vel_sum_loc
  complex, dimension(n_blend,n_x,3) :: vel_sum
  !
  complex, dimension(n_stack,i1_buffer:i2_buffer) :: cap_h
  complex, dimension(n_theta_plot,n_x,3) :: mom_tmp
  !---------------------------------------------------

  if (alltime_index == 0) then
     moments_plot(:,:,:,:) = (0.0,0.0)
  endif
  moments_zero_plot(:,:,:) = 0.0


  do is=1,n_kinetic

!$omp parallel private(p_nek_loc,ie,k,ck,m0,p_nek,i,m,i_diff,j)
     vel_sum_loc(:,ibeg:iend,:) = (0.0,0.0)

     ! First compute gyro_h = <h+z*alpha*<U>>
     !
     ! Note that gyro_h is really a temporary variable 
     ! here; only gyro_h(:,:,:,1) is used.

     p_nek_loc = 0

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)
        k  = nek_k(p_nek)

        ck = class(k)

!$omp barrier
!$omp single
        cap_h(:,i1_buffer:0) = (0.0,0.0)
        cap_h(:,n_x+1:i2_buffer) = (0.0,0.0)
!$omp end single
        do i = ibeg, iend
           do m=1,n_stack
              cap_h(m,i) = h(m,i,p_nek_loc,is)+&
                   z(is)*alpha_s(is,i)*gyro_u(m,i,p_nek_loc,is)
           enddo ! m
        enddo ! i
!$omp barrier  ! ensure all cap_h values are available

        if (is <= n_gk) then
           do i = ibeg, iend
              gyro_h(:,i,p_nek_loc,1) = (0.0,0.0)
              do m=1,n_stack
                 m0 = m_phys(ck,m)
                 do i_diff=-m_gyro,m_gyro-i_gyro
                    gyro_h(m,i,p_nek_loc,1) = gyro_h(m,i,p_nek_loc,1)+ &
                         w_gyro0(m0,i_diff,i,p_nek_loc,is)*cap_h(m,i_loop(i+i_diff))
                 enddo ! i_diff
              enddo ! m
           enddo ! i
        else
           ! Implicit electrons (gyroave=1)
           do i = ibeg, iend
              do m=1,n_stack
                 gyro_h(m,i,p_nek_loc,1) = cap_h(m,i)
              enddo ! m
           enddo ! i
        endif
        
        ! Subtract potential for all moments to give
        !  moment(i) = Int[ d3v w(i) { <cap_h> - e (delta_phi)/T } ]
        !    where w(i)=[1,E,v_par] 
        do i = ibeg, iend
           do m=1,n_stack
              gyro_h(m,i,p_nek_loc,1) = gyro_h(m,i,p_nek_loc,1)-&
                   z(is)*alpha_s(is,i)*field_tau(m,i,p_nek_loc,1)
           enddo ! m
        enddo ! i

        !----------------------------------------------------
        ! Now, compute blending projections:
        !
        ! vel_sum(j,i) -> FV[ (F*_j) gyro_h ]
        !
        do i = ibeg, iend
           do m=1,n_stack
              m0 = m_phys(ck,m)
              do j=1,n_blend
                 ! 1: density moment
                 vel_sum_loc(j,i,1) = vel_sum_loc(j,i,1)+&
                      gyro_h(m,i,p_nek_loc,1)*cs_blend(j,m0,i,p_nek_loc)

                 ! 2: energy moment
                 vel_sum_loc(j,i,2) = vel_sum_loc(j,i,2)+&
                      energy(ie)*tem_s(is,i)*&
                      gyro_h(m,i,p_nek_loc,1)*cs_blend(j,m0,i,p_nek_loc)

                 ! 3: v_parallel moment
                 vel_sum_loc(j,i,3) = vel_sum_loc(j,i,3)+&
                      v_para(m,i,p_nek_loc,is)*&
                      gyro_h(m,i,p_nek_loc,1)*cs_blend(j,m0,i,p_nek_loc)
              enddo ! j
           enddo ! m
        enddo ! i
     enddo ! p_nek
!$omp end parallel
     !--------------------------------------------------------------

     call MPI_ALLREDUCE(vel_sum_loc,&
          vel_sum,&
          size(vel_sum),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     !---------------------------------------------------
     ! Solve for blending coefficients:
     !
!$omp parallel do private(info) schedule(static)
     do i=1,n_x

        call ZGETRS('N',&
             n_blend,&
             3,&
             ff_mm(:,:,i,1),&
             n_blend,&
             ff_mm_piv(:,i,1),&
             vel_sum(:,i,:),&
             n_blend,&
             info)

     enddo ! i
     !---------------------------------------------------

     !---------------------------------------------------------------
     ! Theta-interpolate moments_plot from vel_sum projections: 
     !
     do ix=1,3
        do i=1,n_x
           do j_plot=1,n_theta_plot
              mom_tmp(j_plot,i,ix) = sum(vel_sum(:,i,ix)*blend_plot(:,j_plot,i))
           enddo ! j_plot
        enddo ! i
     enddo ! ix
     !----------------------------------------------------------------

     !---------------------------------------------------------------
     ! Integrate to obtain flux-surface averages
     !
     do ix=1,2
        do i=1,n_x
           do j=1,n_blend
              moments_zero_plot(i,is,ix) = &
                   moments_zero_plot(i,is,ix)+c_fluxave(j,i)*real(vel_sum(j,i,ix))
           enddo ! j
        enddo ! i
     enddo ! ix
     !----------------------------------------------------------------

     !----------------------------------------------------------------
     ! Compute moments_plot with averaging:
     !
     moments_plot(:,:,is,:) = moments_plot(:,:,is,:)+&
          w_time(alltime_index+1)*mom_tmp(:,:,:)
     !----------------------------------------------------------------

  enddo ! is

end subroutine gyro_moments_plot
