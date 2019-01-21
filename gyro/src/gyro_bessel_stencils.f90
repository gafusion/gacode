!-----------------------------------------------------------------
! gyro_bessel_stencils.f90 
!
! PURPOSE:
!  Driver for creation of (tau-space) gyroaverage 
!  stencils w_gyro*.
!-----------------------------------------------------------------

subroutine gyro_bessel_stencils

  use gyro_globals
  use gyro_pointers
  use math_constants

  !-----------------------------------------------------------------
  implicit none
  !
  integer :: p
  !
  real :: rho_gyro
  real :: a_gyro
  real :: u_gyro
  real :: v_gyro
  !
  real :: omega_c
  !
  complex, dimension(-m_gyro:m_gyro-i_gyro) :: g
  !
  complex, dimension(:,:), allocatable :: w_temp0
  complex, dimension(:,:), allocatable :: w_temp1
  complex, dimension(:,:), allocatable :: w_temp2
  complex, dimension(:,:), allocatable :: w_temp3
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  ! Set dimensions and functions
  !
  do p=-n_x/2,n_x/2-1
     do m=-n_x/2,n_x/2
        z_gyro(m,p) = exp(-i_c*p*m*pi/(n_x/2))
     enddo
  enddo
  !-----------------------------------------------------------------

  allocate(w_temp0(n_theta(2),-m_gyro:m_gyro-i_gyro))
  allocate(w_temp2(n_theta(2),-m_gyro:m_gyro-i_gyro))
  if (n_field == 3) then
     allocate(w_temp1(n_theta(2),-m_gyro:m_gyro-i_gyro))
     allocate(w_temp3(n_theta(2),-m_gyro:m_gyro-i_gyro))
  endif

  !-----------------------------------------------------------------
  ! LOCAL gyroaverage:
  !
  do is=1,n_gk
     p_nek_loc = 0
     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)  
        k  = nek_k(p_nek)   

        ck = class(k)

        do i=1,n_x

           if (i == 1 .or. flat_profile_flag == 0) then

              do m0=1,n_theta(ck)

                 !---------------------------------------------------------------
                 ! Prepare argument of Bessel function
                 !  
                 omega_c = abs(z(is))*b_unit_s(i)*mu(is)**2

                 if (kill_gyro_b_flag == 0) then
                    omega_c = omega_c*b0_t(i,k,m0)
                 endif

                 rho_gyro = rhos_norm*v_perp(m0,i,p_nek_loc,is)/omega_c
                 !
                 a_gyro = grad_r_t(i,k,m0)/x_length
                 u_gyro = qrat_t(i,k,m0)*n_1(in_1)*q_s(i)/r_s(i)*captheta_t(i,k,m0)
                 v_gyro = qrat_t(i,k,m0)*n_1(in_1)*q_s(i)/r_s(i)
                 !----------------------------------------------------------------

                 ! G0a = J0
                 call gyro_bessel_operator(rho_gyro,&
                      a_gyro,&
                      u_gyro,&
                      v_gyro,&
                      g,&
                      1)

                 w_temp0(m0,:) = g(:)

                 ! G2a = -(i/2)*k_x*rho*[ J_0(z)+J_2(z) ]
                 call gyro_bessel_operator(rho_gyro,&
                      a_gyro,&
                      u_gyro,&
                      v_gyro,&
                      g, &
                      3)

                 w_temp2(m0,:) = g(:)

                 if (n_field == 3) then

                    ! G1a = (1/2)*[ J_0(z)+J_2(z) ] 
                    call gyro_bessel_operator(rho_gyro,&
                         a_gyro,&
                         u_gyro,&
                         v_gyro,&
                         g,&
                         4)

                    w_temp1(m0,:) = g(:)

                    ! G3a = i*k_x*rho*[ J_0(z)-J_1(z)/z ] / z^2 
                    call gyro_bessel_operator(rho_gyro,&
                         a_gyro,&
                         u_gyro,&
                         v_gyro,&
                         g,&
                         8)

                    w_temp3(m0,:) = g(:)

                 endif

              enddo ! m0

              do m=1,n_stack
                 m0 = m_phys(ck,m)
                 w_gyro0(m,:,i,p_nek_loc,is) = w_temp0(m0,:)
                 w_gyro2(m,:,i,p_nek_loc,is) = w_temp2(m0,:)
                 if (n_field == 3) then
                    w_gyro1(m,:,i,p_nek_loc,is) = w_temp1(m0,:)
                    w_gyro3(m,:,i,p_nek_loc,is) = w_temp3(m0,:)
                 endif
              enddo

           else
              w_gyro0(:,:,i,p_nek_loc,is) = w_gyro0(:,:,1,p_nek_loc,is)
              w_gyro2(:,:,i,p_nek_loc,is) = w_gyro2(:,:,1,p_nek_loc,is)
              if (n_field == 3) then
                 w_gyro1(:,:,i,p_nek_loc,is) = w_gyro1(:,:,1,p_nek_loc,is)
                 w_gyro3(:,:,i,p_nek_loc,is) = w_gyro3(:,:,1,p_nek_loc,is)
              endif
           endif

        enddo ! i
     enddo ! p_nek
  enddo ! is
  !
  !----------------------------------------------------------

  deallocate(w_temp0)
  deallocate(w_temp2)
  if (n_field == 3) then
     deallocate(w_temp1)
     deallocate(w_temp3)
  endif

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_bessel_stencils done]'
  endif

end subroutine gyro_bessel_stencils
