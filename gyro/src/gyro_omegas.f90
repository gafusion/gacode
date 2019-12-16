!---------------------------------------------------------
! gyro_omegas.f90
!
! PURPOSE:
!  Compute normalized drift coefficients, velocities, 
!  etc., of the gyrokinetic equation.  
!---------------------------------------------------------

subroutine gyro_omegas

  use gyro_globals
  use gyro_pointers
  use math_constants

  !-------------------------------------
  implicit none
  !
  integer :: p
  !
  real :: temp
  real :: e_temp
  real :: e_temp_mach
  real :: e_temp_p
  !-------------------------------------

  ! v_theta: note that 1/g_theta factor inside omega(i,k)

  do i=1,n_x
     do ie=1,n_energy
        do k=1,n_lambda    
           do is=1,n_kinetic
              v_theta(i,ie,k,is) = & 
                   omega(i,k)*sqrt(2.0*energy(ie))/(q_s(i)*rmaj_s(i))
           enddo ! is
        enddo ! k
     enddo ! ie
  enddo ! i

  omega_d1 = 0.0
  omega_dr = 0.0

  p_nek_loc = 0

  do p_nek=1+i_proc_1,n_nek_1,n_proc_1

     p_nek_loc = p_nek_loc+1

     ie = nek_e(p_nek)  
     k  = nek_k(p_nek)   

     ck = class(k)

     do i=1,n_x

        do m=1,n_stack

           m0 = m_phys(ck,m)
           p  = p_phys(ck,m) 

           !-----------------------------------------------------
           ! STANDARD CURVATURE/GRAD-B DRIFTS
           !-----------------------------------------------------

           do is=1,n_kinetic

              e_temp = (1.0/z(is))*(2.0*tem_s(is,i)/rmaj_s(i))*&
                   energy(ie)*(1.0-0.5*lambda(i,k)*b0_t(i,k,m0))/b0_t(i,k,m0)

              e_temp_p = e_temp*(1.0-lambda(i,k)*b0_t(i,k,m0))/ &
                   (1.0-0.5*lambda(i,k)*b0_t(i,k,m0))

              ! krho_i = n*q_s(r)/r_s(i)*rhos_norm/b_unit_s(i)

              omega_d1(m,i,p_nek_loc,is) = krho_i(in_1,i)*qrat_t(i,k,m0)* &
                   ( e_temp*(cos_t(i,k,m0)+captheta_t(i,k,m0)*sin_t(i,k,m0)) + &
                   e_temp_p*cos_p_t(i,k,m0) )

              ! omega_dr with nonuniform grid effect:

              omega_dr(m,i,p_nek_loc,is) = rhos_norm/b_unit_s(i)*  &
                   grad_r_t(i,k,m0)*e_temp*sin_t(i,k,m0)

           enddo ! is

        enddo ! m0

        do m=1,n_stack

           m0 = m_phys(ck,m)
           p  = p_phys(ck,m)

           temp = sigma(ck,p)* &
                sqrt(abs(energy(ie)*(1.0-lambda(i,k)*b0_t(i,k,m))))

           !-------------------------------------------------
           ! Parallel and perp. velocity for each species:
           !
           do is=1,n_kinetic
              v_para(m,i,p_nek_loc,is) = &
                   mu(is)*sqrt(2.0*tem_s(is,i))*temp
              v_perp(m,i,p_nek_loc,is) = &
                   mu(is)*sqrt(2.0*tem_s(is,i)*&
                   energy(ie)*lambda(i,k)*b0_t(i,k,m))
           enddo
           !-------------------------------------------------

           !-------------------------------------------------------- 
           ! Standard part of omega_star
           !
           do is=1,n_kinetic
              omega_star(m,i,p_nek_loc,is) = krho_i(in_1,i)*&
                   (dlnndr_s(is,i)+(energy(ie)-1.5)*dlntdr_s(is,i))*&
                   den_s(is,i)
           enddo
           !--------------------------------------------------------

           !-----------------------------------------------------
           ! ROTATION SHEAR (KELVIN-HELMHOLTZ) DRIVE 
           !-----------------------------------------------------

           ! m_a  v_para_a      1     R Bt   a gamma_p
           ! --- ----------  ------- ------ -----------
           ! m_i  c_s_hat    T_a_hat  R0 B    c_s_hat

           do is=1,n_kinetic

              omega_star(m,i,p_nek_loc,is) = omega_star(m,i,p_nek_loc,is)+&
                   krho_i(in_1,i)*gamma_p_s(i)*bigr_t(i,k,m0)/rmaj_s(i)*&
                   bt_t(i,k,m0)/b0_t(i,k,m0)*v_para(m,i,p_nek_loc,is)/ &
                   (tem_s(is,i)*mu(is)**2)*&
                   den_s(is,i)

           enddo ! is


           !-----------------------------------------------------
           ! CORIOLIS DRIFT 
           !-----------------------------------------------------

           ! m_a  v_para_a                       2
           ! --- ---------- M sqrt(Te_hat) -----------------  
           ! m_i  c_s_hat                   z_a B_hat R0_hat

           do is=1,n_kinetic

              e_temp_mach = 2.0/(z(is)*mu(is)**2*rmaj_s(i))*&
                   v_para(m,i,p_nek_loc,is)*mach_s(i)*&
                   sqrt(tem_s(indx_e,i))/b0_t(i,k,m0)

              omega_d1(m,i,p_nek_loc,is) = omega_d1(m,i,p_nek_loc,is)+ &
                   krho_i(in_1,i)*qrat_t(i,k,m0)*e_temp_mach* &
                   (ucos_t(i,k,m0)+captheta_t(i,k,m0)*usin_t(i,k,m0))

              omega_dr(m,i,p_nek_loc,is) = omega_dr(m,i,p_nek_loc,is) + &
                   rhos_norm/b_unit_s(i)*grad_r_t(i,k,m0)*&
                   e_temp_mach*usin_t(i,k,m0)

           enddo ! is

        enddo ! m

     enddo ! i

  enddo ! p_nek_loc


  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_omegas done]'
  endif

end subroutine gyro_omegas
