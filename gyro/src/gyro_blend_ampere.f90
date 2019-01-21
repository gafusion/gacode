!------------------------------------------------------------
! gyro_blend_ampere.f90
!
! PURPOSE:
!  Generate matrix L_A of paper by pushing the factor 
!  [-nabla_perp A] into a velocity integral. 
!------------------------------------------------------------

subroutine gyro_blend_ampere

  use mpi
  use gyro_globals
  use gyro_pointers
  use math_constants

  !---------------------------------------------------
  implicit none
  !
  real :: betae_eff
  !
  real :: omega_c
  real :: rho_gyro
  real :: a_gyro
  real :: u_gyro
  real :: v_gyro
  !
  complex :: inqr
  !
  complex, dimension(-mg_dx:mg_dx-ig_dx,n_blend,n_blend) :: vel_sum_loc
  complex, dimension(-mg_dx:mg_dx-ig_dx,n_blend,n_blend) :: vel_sum_glob
  complex, dimension(n_gk,-mg_dx:mg_dx-ig_dx) :: f_x
  complex, dimension(-m_dx:m_dx-i_dx) :: grad_perp_ap
  complex, dimension(n_gk,-mg_dx:mg_dx-ig_dx) :: ion_current
  complex, dimension(n_gk,-mg_dx:mg_dx-ig_dx,n_stack,n_nek_loc_1) :: f_save
  !---------------------------------------------------

  betae_eff = betae_unit_norm*ampere_scale

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

           grad_perp_ap(:) =-2.0*rhos_norm**2/betae_eff*( &
                inqr**2*qrat_t(i,k,m)**2* &
                (1.0+captheta_t(i,k,m)**2)*w_d0(:)+ &
                2.0*inqr*qrat_t(i,k,m)*captheta_t(i,k,m)* &
                grad_r_t(i,k,m)*w_d1(:)+ &
                (grad_r_t(i,k,m))**2*w_d2(:))

           ! Fast option for flat profiles

           if (i == 1 .or. flat_profile_flag == 0) then

              do is=1,n_gk       

                 if (gyro_method == 1) then 

                    ion_current(is,:) = &
                         v_para(m,i,p_nek_loc,is)**2* &
                         alpha_s(is,i)*z(is)**2*w_gd0(:)

                 else

                    !--------------------------------------------
                    ! Prepare argument of Bessel function
                    !
                    omega_c = abs(z(is))*b_unit_s(i)*mu(is)**2
                    !
                    if (kill_gyro_b_flag == 0) then
                       omega_c = omega_c*b0_t(i,k,m)
                    endif
                    !
                    rho_gyro = rhos_norm*v_perp(m,i,p_nek_loc,is)/omega_c
                    !
                    a_gyro = grad_r_t(i,k,m)/x_length
                    v_gyro = qrat_t(i,k,m)*n_1(in_1)*q_s(i)/r_s(i)
                    u_gyro = v_gyro*captheta_t(i,k,m)
                    !--------------------------------------------

                    f_x(is,:) = (0.0,0.0) 

                    ! J_0^2
                    call gyro_bessel_operator(rho_gyro,&
                         a_gyro,&
                         u_gyro,&
                         v_gyro,&
                         f_x(is,:),&
                         2)

                    ion_current(is,:) = &
                         v_para(m,i,p_nek_loc,is)**2* &
                         alpha_s(is,i)*z(is)**2*f_x(is,:)

                 endif

              enddo ! is
              f_save(:,:,m,p_nek_loc) = ion_current(:,:)
           else
              ion_current(:,:) = f_save(:,:,m,p_nek_loc)
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
                         grad_perp_ap(i_diff)* &
                         cs_blend(j,m0,i,p_nek_loc)* &
                         c_blend(jp,m0,ip,p_nek_loc)

                 enddo ! i_diff

                 do i_diff=-mg_dx,mg_dx-ig_dx

                    !-----------------------------------------------
                    ! This will stagnate at i=n_x for i+i_diff > n_x
                    ! if boundary_method=2:
                    ip = i_cyc(i+i_diff)
                    !------------------------------------------------

                    vel_sum_loc(i_diff,j,jp) = vel_sum_loc(i_diff,j,jp)+&
                         fakefield_flag*sum(ion_current(:,i_diff))* &
                         cs_blend(j,m0,i,p_nek_loc)* &
                         c_blend(jp,m0,ip,p_nek_loc)

                 enddo ! i_diff

              enddo ! jp
           enddo ! j

        enddo ! m

     enddo ! p_nek_loc

     call MPI_ALLREDUCE(vel_sum_loc,&
          vel_sum_glob,&
          n_blend*n_blend*(2*mg_dx-ig_dx+1),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     aa_mm(i,:,:,:) = vel_sum_glob(:,:,:)

  enddo ! i

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_ampere_blend done]'
  endif

end subroutine gyro_blend_ampere
