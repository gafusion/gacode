!------------------------------------------------
! gyro_blend_ampereperp.f90
!
! PURPOSE:
!  Generate matrices L_BB and LBP (LPB=-2*LBP) of paper.
!------------------------------------------------

subroutine gyro_blend_ampereperp

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
  complex, dimension(-m_gyro:m_gyro-i_gyro,n_blend,n_blend) :: vel_sum_loc1
  complex, dimension(n_gk,-m_gyro:m_gyro-i_gyro) :: f_x1
  complex, dimension(-m_gyro:m_gyro-i_gyro,n_blend,n_blend) :: vel_sum_loc2
  complex, dimension(n_gk,-m_gyro:m_gyro-i_gyro) :: f_x2
  complex, dimension(-m_gyro:m_gyro-i_gyro,n_blend,n_blend) :: vel_sum_glob
  complex, dimension(n_gk,-m_gyro:m_gyro-i_gyro,n_stack,n_nek_loc_1) :: f1_save
  complex, dimension(n_gk,-m_gyro:m_gyro-i_gyro,n_stack,n_nek_loc_1) :: f2_save
  !---------------------------------------------------

  betae_eff = betae_unit_norm*ampere_scale

  do i=1,n_x

     vel_sum_loc1 = (0.0,0.0)
     vel_sum_loc2 = (0.0,0.0)

     p_nek_loc = 0

     do p_nek=1+i_proc_1,n_nek_1,n_proc_1

        p_nek_loc = p_nek_loc+1

        ie = nek_e(p_nek)
        k  = nek_k(p_nek)

        ck = class(k)

        ! Now, compute m-projections of the RHS:

        do m=1,n_stack

           m0 = m_phys(ck,m)

           ! Fast option for flat profiles

           if (i == 1 .or. flat_profile_flag == 0) then

              do is=1,n_gk

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
                 !---------------------------------------------------------------

                 f_x1(is,:) = (0.0,0.0) 
                 f_x2(is,:) = (0.0,0.0)

                 ! G_perp^2
                 call gyro_bessel_operator(rho_gyro,&
                      a_gyro,&
                      u_gyro,&
                      v_gyro,&
                      f_x1(is,:),&
                      5)

                 ! G_perp * G
                 call gyro_bessel_operator(rho_gyro,&
                      a_gyro,&
                      u_gyro,&
                      v_gyro,&
                      f_x2(is,:),&
                      6)

              enddo ! is
              f1_save(:,:,m,p_nek_loc) = f_x1(:,:)
              f2_save(:,:,m,p_nek_loc) = f_x2(:,:)
           else
              f_x1(:,:) = f1_save(:,:,m,p_nek_loc)
              f_x2(:,:) = f2_save(:,:,m,p_nek_loc)
           endif

           do j=1,n_blend
              do jp=1,n_blend

                 do i_diff=-m_gyro,m_gyro-i_gyro

                    !-----------------------------------------------
                    ! This will stagnate at i=n_x for i+i_diff > n_x
                    ! if boundary_method=2:
                    ip = i_cyc(i+i_diff)
                    !-----------------------------------------------

                    ! FV[ (F*_j) (2*z_i*n_i*T_i*V[(ene*lambda*G_perp)^2]
                    !             + 1/beta_e) (F_jp) ]
                    !
                    !  = L_B 

                    vel_sum_loc1(i_diff,j,jp) = vel_sum_loc1(i_diff,j,jp)+&
                         fakefield_flag*(sum(2.0*den_s(1:n_gk,i)*tem_s(1:n_gk,i)&
                         * f_x1(1:n_gk,i_diff)*energy(ie)**2)) &
                         * lambda(i,k)**2 &
                         * cs_blend(j,m0,i,p_nek_loc) &
                         * c_blend(jp,m0,ip,p_nek_loc) &
                         + b_unit_s(i)**2*1.0/betae_eff * w_g0(i_diff) &
                         * cs_blend(j,m0,i,p_nek_loc) &
                         * c_blend(jp,m0,ip,p_nek_loc)              

                    ! FV[ (F*_j) (z_i*n_i*V[(ene*lambda*G*G_perp)]) (F_jp) ]
                    !
                    !  = L_BP 

                    vel_sum_loc2(i_diff,j,jp) = vel_sum_loc2(i_diff,j,jp)+&
                         fakefield_flag*(sum(z(1:n_gk)*den_s(1:n_gk,i)&
                         * f_x2(1:n_gk,i_diff)*energy(ie))) &
                         * lambda(i,k) &
                         * cs_blend(j,m0,i,p_nek_loc) &
                         * c_blend(jp,m0,ip,p_nek_loc)

                 enddo ! i_diff

              enddo ! jp
           enddo ! j

        enddo ! m 
     enddo ! p_nek_loc 
     !-------------------------------------------------------------------

     call MPI_ALLREDUCE(vel_sum_loc1,&
          vel_sum_glob,&
          n_blend*n_blend*(2*m_gyro-i_gyro+1),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     ab_mm(i,:,:,:) = vel_sum_glob(:,:,:)

     call MPI_ALLREDUCE(vel_sum_loc2,&
          vel_sum_glob,&
          n_blend*n_blend*(2*m_gyro-i_gyro+1),&
          MPI_DOUBLE_COMPLEX,&
          MPI_SUM,&
          NEW_COMM_1,&
          i_err)

     abp_mm(i,:,:,:) = vel_sum_glob(:,:,:)

  enddo ! i 

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[make_ampereperp_blend done]'
  endif

end subroutine gyro_blend_ampereperp
