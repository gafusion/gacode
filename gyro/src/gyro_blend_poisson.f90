!----------------------------------------------------------
! gyro_blend_poisson.f90
!
! PURPOSE:
!  Generate matrix L_P of paper.
!
!  i_elec = 0: 
!         = 1:
!----------------------------------------------------------

subroutine gyro_blend_poisson(i_elec)

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  integer, intent(in) :: i_elec
  !
  real :: omega_c
  real :: rho_gyro
  real :: a_gyro
  real :: u_gyro
  real :: v_gyro
  !
  complex :: inqr
  ! 
  complex, dimension(-m_gyro:m_gyro-i_gyro,n_blend,n_blend) :: vel_sum_loc
  complex, dimension(-m_gyro:m_gyro-i_gyro,n_blend,n_blend) :: vel_sum_glob
  complex, dimension(n_gk,-m_gyro:m_gyro-i_gyro) :: f_x
  complex, dimension(-m_dx:m_dx-i_dx) :: grad_perp_phi
  complex, dimension(n_gk,-m_gyro:m_gyro-i_gyro,n_stack,n_nek_loc_1) :: f_save
  !
  real, dimension(-m_gyro:m_gyro-i_gyro) :: trace_1
  real, dimension(n_gk,-m_gyro:m_gyro-i_gyro) :: trace_2
  real, dimension(-m_gyro:m_gyro-i_gyro) :: trace_1_glob
  real, dimension(n_gk,-m_gyro:m_gyro-i_gyro) :: trace_2_glob
  !
  complex :: cprod 
  !---------------------------------------------------

  !------------------------------------------------------------------
  ! Now, compute matrix of blending projections:
  !
  do i=1,n_x

     vel_sum_loc = (0.0,0.0)

     p_nek_loc = 0

     inqr = i_c*n_1(in_1)*q_s(i)/r_s(i) 

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)
        k  = nek_k(p_nek)

        ck = class(k)

        ! Now, compute m-projections of the RHS:

        do m=1,n_stack

           m0 = m_phys(ck,m)

           grad_perp_phi(:) =-lambda_debye**2*( &
                inqr**2*qrat_t(i,k,m)**2* &
                (1.0+captheta_t(i,k,m)**2)*w_d0(:)+ &
                2.0*inqr*qrat_t(i,k,m)*captheta_t(i,k,m)* &
                grad_r_t(i,k,m)*w_d1(:)+ &
                (grad_r_t(i,k,m))**2*w_d2(:))

           ! Fast option for flat profiles

           if (i == 1 .or. flat_profile_flag == 0) then

              do is=1,n_gk

                 if (gyro_method == 1) then

                    !---------------------------------------------------------------
                    ! Prepare argument of Bessel function
                    !  
                    omega_c = abs(z(is))*b_unit_s(i)*mu(is)**2

                    if (kill_gyro_b_flag == 0) then
                       omega_c = omega_c*b0_t(i,k,m)
                    endif

                    rho_gyro = rhos_norm*mu(is)*sqrt(tem_s(is,i))/omega_c
                    !
                    a_gyro = grad_r_t(i,k,m)/x_length
                    v_gyro = qrat_t(i,k,m)*n_1(in_1)*q_s(i)/r_s(i)
                    u_gyro = v_gyro*captheta_t(i,k,m)
                    !---------------------------------------------------------------

                    f_x(is,:) = (0.0,0.0) 

                    ! I_0
                    call gyro_bessel_operator(rho_gyro,&
                         a_gyro,&
                         u_gyro,&
                         v_gyro,&
                         f_x(is,:),&
                         7)

                 else 

                    !---------------------------------------------------------------
                    ! Prepare argument of Bessel function
                    !  
                    omega_c = abs(z(is))*b_unit_s(i)*mu(is)**2

                    if (kill_gyro_b_flag == 0) then
                       omega_c = omega_c*b0_t(i,k,m)
                    endif

                    rho_gyro = rhos_norm*v_perp(m,i,p_nek_loc,is)/omega_c
                    !
                    a_gyro = grad_r_t(i,k,m)/x_length
                    v_gyro = qrat_t(i,k,m)*n_1(in_1)*q_s(i)/r_s(i)
                    u_gyro = v_gyro*captheta_t(i,k,m)
                    !---------------------------------------------------------------

                    f_x(is,:) = (0.0,0.0) 

                    ! J_0^2
                    call gyro_bessel_operator(rho_gyro,&
                         a_gyro,&
                         u_gyro,&
                         v_gyro,&
                         f_x(is,:),&
                         2)

                 endif

              enddo ! is
              f_save(:,:,m,p_nek_loc) = f_x(:,:)
           else
              f_x(:,:) = f_save(:,:,m,p_nek_loc)
           endif

           do j=1,n_blend
              do jp=1,n_blend

                 do i_diff=-m_dx,m_dx-i_dx

                    !-----------------------------------------------
                    ! This will stagnate at i=n_x for i+i_diff > n_x
                    ! if boundary_method=2:
                    ip = i_cyc(i+i_diff)
                    !-----------------------------------------------

                    vel_sum_loc(i_diff,j,jp) = vel_sum_loc(i_diff,j,jp)+&
                         grad_perp_phi(i_diff)* &
                         cs_blend(j,m0,i,p_nek_loc)* &
                         c_blend(jp,m0,ip,p_nek_loc)

                 enddo ! i_diff

                 do i_diff=-m_gyro,m_gyro-i_gyro

                    !-----------------------------------------------
                    ! This will stagnate at i=n_x for i+i_diff > n_x
                    ! if boundary_method=2:
                    ip = i_cyc(i+i_diff)
                    !-----------------------------------------------

                    cprod = cs_blend(j,m0,i,p_nek_loc)*c_blend(jp,m0,ip,p_nek_loc)

                    ! FV[ (F*_j) (alpha_i (1-R) + alpha_e) (F_jp) ]

                    select case (electron_method) 
                    case (1,2)
                       vel_sum_loc(i_diff,j,jp) = vel_sum_loc(i_diff,j,jp)+&
                            (sum(alpha_s(1:n_gk,i)*z(1:n_gk)**2* &
                            (w_g0(i_diff)-fakefield_flag*f_x(1:n_gk,i_diff))) + &
                            i_elec*alpha_s(indx_e,i)*w_g0(i_diff))*cprod
                    case (3)
                       ! is=1 are electrons, is>1 are ions.
                       vel_sum_loc(i_diff,j,jp) = vel_sum_loc(i_diff,j,jp)+&
                            (alpha_s(1,i)*z(1)**2* &
                            (w_g0(i_diff)-fakefield_flag*f_x(1,i_diff)) + &
                            sum(alpha_s(2:n_spec,i)*z(2:n_spec)**2)*w_g0(i_diff))*cprod
                    case (4)
                       vel_sum_loc(i_diff,j,jp) = vel_sum_loc(i_diff,j,jp)+&
                            sum(alpha_s(1:n_gk,i)*z(1:n_gk)**2* &
                            (w_g0(i_diff)-fakefield_flag*f_x(1:n_gk,i_diff)))*cprod
                    end select

                 enddo ! i_diff

              enddo ! jp
           enddo ! j

        enddo ! m 

     enddo ! p_nek_loc 
     !-------------------------------------------------------------------

     call MPI_ALLREDUCE(vel_sum_loc,&
          vel_sum_glob,&
          size(vel_sum_glob),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     ap_mm(i,:,:,:) = vel_sum_glob(:,:,:)

     if (electron_method == 1 .and. n_1(in_1) == 0) then

        !---------------------
        ! ADIABATIC ELECTRONS:
        !---------------------

        do j=1,n_blend
           do jp=1,n_blend

              ap_mm(i,0,j,jp) = ap_mm(i,0,j,jp)-&
                   alpha_s(indx_e,i)*c_fluxave(j,i)*c_fluxave(jp,i)

           enddo ! jp
        enddo ! j

     endif

  enddo ! i 

  !---------------------------------------------------------------------------
  ! Compute measure of quality of double-gyroaverage truncation:
  !
  i = ir_norm

  trace_1(:) = 0.0
  trace_2(:,:) = 0.0

  p_nek_loc = 0

  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)
     k  = nek_k(p_nek)

     ck = class(k)

     ! Now, compute m-projections of the RHS:

     do m=1,n_stack

        m0 = m_phys(ck,m)

        do is=1,n_gk

           !---------------------------------------------------------------
           ! Prepare argument of Bessel function
           !  
           omega_c = abs(z(is))*b_unit_s(i)*mu(is)**2

           if (kill_gyro_b_flag == 0) then
              omega_c = omega_c*b0_t(i,k,m0)
           endif

           rho_gyro = rhos_norm*mu(is)*sqrt(tem_s(is,i))/omega_c
           !
           a_gyro = grad_r_t(i,k,m)/x_length
           v_gyro = qrat_t(i,k,m)*n_1(in_1)*q_s(i)/r_s(i)
           u_gyro = v_gyro*captheta_t(i,k,m)
           !---------------------------------------------------------------

           ! I_0
           call gyro_bessel_operator(rho_gyro,&
                a_gyro,&
                u_gyro,&
                v_gyro,&
                f_x(is,:),&
                7)

        enddo ! is

        do i_diff=-m_gyro,m_gyro-i_gyro

           !-----------------------------------------------
           ! This will stagnate at i=n_x for i+i_diff > n_x
           ! if boundary_method=2:
           ip = i_cyc(i+i_diff)
           !-----------------------------------------------

           do j=1,n_blend

              trace_1(i_diff) = trace_1(i_diff)+&
                   real(w_g0(i_diff)* &
                   cs_blend(j,m0,i,p_nek_loc)*c_blend(j,m0,i,p_nek_loc))

              trace_2(:,i_diff) = trace_2(:,i_diff)+&
                   real(f_x(:,i_diff)* &
                   cs_blend(j,m0,i,p_nek_loc)*c_blend(j,m0,i,p_nek_loc))

           enddo ! j

        enddo ! i_diff

     enddo ! m 

  enddo ! p_nek_loc 
  !-------------------------------------------------------------------

  ! We only want the n=0 trace.  

  if (n_1(in_1) /= 0) then
     trace_1(:) = 0.0
     trace_2(:,:) = 0.0
  endif

  call MPI_ALLREDUCE(trace_1,&
       trace_1_glob,&
       size(trace_1_glob),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  call MPI_ALLREDUCE(trace_2,&
       trace_2_glob,&
       size(trace_2_glob),&
       MPI_DOUBLE_PRECISION,&
       MPI_SUM,&
       NEW_COMM_1,&
       i_err)

  if (n_1(in_1) == 0) then
     gyro_trace(:,:) = trace_2_glob(:,:)/trace_1_glob(0)
  else
     gyro_trace(:,:) = 0.0
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_poisson_blend done]'
  endif

end subroutine gyro_blend_poisson
