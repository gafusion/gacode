!-------------------------------------------------------------------------
! gyro_phi_kp_squared.f90
!
! PURPOSE:
!  Calculate
!
!  (1) phi_squared   : (nq/r) F |phi|^2 
!  (2) k_perp_squared: F < k_perp^2 rhos_norm**2 |phi|^2 > / F < |phi|^2 >
!
! Here, F is a flux-surface average and <> is a radial average.
!-------------------------------------------------------------------------

subroutine gyro_phi_kp_squared

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants
  use ompdata

  !---------------------------------------------------
  implicit none
  !
  complex, dimension(n_x,2) :: moment
  complex, dimension(n_x,2) :: moment_glob
  !
  complex :: inqr
  real :: k_perp_squared_loc(1)
  !
  complex, dimension(i1_buffer:i2_buffer,n_stack) :: k0phi
  complex k2phi
  !
  real, dimension(n_field) :: ave_phi_0_loc
  real, dimension(n_field) :: ave_phi_n_loc
  real, dimension(n_field) :: ave_phi_0
  real, dimension(n_field) :: ave_phi_n
  !--------------------------------------------------  

!$omp parallel private(p_nek_loc,ie,k,k2phi,inqr,ip,p_nek,m,i,i_diff)
  moment(ibeg:iend,:) = 0.0

  !-----------------------------------------------------------
  p_nek_loc = 0
  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek) 
     k  = nek_k(p_nek)

!$omp barrier
!$omp single
     do m=1,n_stack
        k0phi(i1_buffer:0,m)     = 0.0
        k0phi(1:n_x,m)           = field_tau(m,:,p_nek_loc,1)
        k0phi(n_x+1:i2_buffer,m) = 0.0
     end do
!$omp end single
!$omp barrier  ! ensure all k0phi values are available

        ! k_perp^2 phi

     do i = ibeg, iend
        do m=1,n_stack
           k2phi = 0.0
           inqr = i_c*n_1(in_1)*q_s(i)/r_s(i) 
           do i_diff=-m_dx,m_dx-i_dx
              ip = i+i_diff
              k2phi = k2phi-(inqr**2*qrat_t(i,k,m)**2* &
                   (1.0+captheta_t(i,k,m)**2)*w_d0(i_diff)+ &
                   2.0*inqr*qrat_t(i,k,m)*captheta_t(i,k,m)* &
                   grad_r_t(i,k,m)*dr_eodr(i)*w_d1(i_diff)+ &
                   (grad_r_t(i,k,m)*dr_eodr(i))**2*w_d2(i_diff))*&
                   k0phi(i_loop(ip),m)
           enddo

           moment(i,1) = moment(i,1)+w_p(ie,i,k)*k0phi(i,m)*conjg(k0phi(i,m))
           moment(i,2) = moment(i,2)+w_p(ie,i,k)*k2phi*conjg(k0phi(i,m))
        enddo
     enddo

  enddo ! p_nek
!$omp end parallel
  !
  !------------------------------------------------------------

  !------------------------------------------------------------
  ! ALLREDUCE to obtain full velocity-space sums.  The nonlinear
  ! fluxes will be correct on all processors.
  !
  call MPI_ALLREDUCE(moment, &
       moment_glob, &
       size(moment_glob), &
       MPI_DOUBLE_COMPLEX, &
       MPI_SUM, &
       NEW_COMM_1, &
       i_err)
  !------------------------------------------------------------

  phi_squared(:) = real(krho_i(in_1,:)*moment_glob(:,1))

  if (abs(sum(moment_glob(:,1))) > 0.0) then
     k_perp_squared_loc(1) = abs(sum(moment_glob(:,2))/sum(moment_glob(:,1))*rhos_norm**2) 
  else
     k_perp_squared_loc(1) = 0.0
  endif

  call collect_real(k_perp_squared_loc(1),k_perp_squared(:))

  ! Now, get RMS field averages

  ave_phi_0_loc(:) = 0.0
  ave_phi_n_loc(:) = 0.0

  if (n_1(in_1) /= 0) then

     ! nonzero modes

     do i=1,n_x
        do j=1,n_theta_int
           ave_phi_n_loc(:) = ave_phi_n_loc(:)+&
                2.0*real(phi(j,i,:)*conjg(phi(j,i,:)))
        enddo
     enddo

  else

     ! zero modes

     do i=1,n_x
        do j=1,n_theta_int
           ave_phi_0_loc(:) = ave_phi_0_loc(:)+&
                real(phi(j,i,:)*conjg(phi(j,i,:)))
        enddo
     enddo

  endif

  call MPI_ALLREDUCE(ave_phi_0_loc,&
       ave_phi_0,&
       n_field,&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_2,&
       i_err)

  call MPI_ALLREDUCE(ave_phi_n_loc,&
       ave_phi_n,&
       n_field,&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_2,&
       i_err)

  ave_phi(1,:) = sqrt(ave_phi_0(:)/(n_theta_int*n_x))
  ave_phi(2,:) = sqrt(ave_phi_n(:)/(n_theta_int*n_x))

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_phi_kp_squared called]'
  endif

end subroutine gyro_phi_kp_squared
