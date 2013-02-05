!----------------------------------------------------------
! gyro_alloc_distrib
!
! PURPOSE:
!  Create and destroy orbit arrays.
!
! NOTES:
!  flag=0: deallocate
!  flag=1: allocate
!
!  Need to know n_lambda,n_stack before this is called.
!
!  See: gyro_banana_operators, gyro_geometry_arrays.
!----------------------------------------------------------

subroutine gyro_alloc_distrib(flag)

  use gyro_globals 
  use gyro_pointers

  !-----------------------------------------
  implicit none
  !
  integer, intent(in) :: flag
  !-----------------------------------------

  if (flag == 1 .and. allocated(omega_d1)) then
     if (i_proc == 0) then
        print *,'WARNING: already allocated arrays in gyro_alloc_distrib'
     endif
     return
  endif
  if (flag == 0 .and. .not.allocated(omega_d1)) then
     if (i_proc == 0) then
        print *,'WARNING: cannot deallocate arrays in gyro_alloc_distrib'
     endif
     return
  endif

  if (flag == 1) then

     ! Drifts
     allocate(omega_d1(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(omega_dr(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(omega_star(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(v_theta(n_x,n_energy,n_lambda,n_kinetic))
     allocate(v_para(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(v_perp(n_stack,n_x,n_nek_loc_1,n_kinetic))

     ! Pointers
     allocate(nek_n(n_nek_1))
     allocate(nek_e(n_nek_1))
     allocate(nek_k(n_nek_1))
     allocate(ine_i(n_ine_1))
     allocate(ine_n(n_ine_1))
     allocate(ine_e(n_ine_1))
     allocate(ki_k(n_ki_1))
     allocate(ki_i(n_ki_1))

     ! Blending arrays
     allocate(cs_blend(n_blend,n_theta(2),n_x,n_nek_loc_1))
     allocate(c_blend(n_blend,n_stack,n_x,n_nek_loc_1))
     allocate(cs_blend_prime(n_blend,n_theta(2),n_x,n_nek_loc_1))

     if (electron_method == 2) then
        ! Implicit advection
        allocate(o_f(n_blend,n_stack,n_x,n_nek_loc_1))
        allocate(o_fv(n_blend,n_stack,n_x,n_nek_loc_1))
        allocate(o_advect(n_stack,n_stack,n_x,n_nek_loc_1))
        allocate(imp(n_x,n_blend,n_blend,8))
     endif

     if (collision_flag == 1) then

        ! Collision arrays
        if (linsolve_method == 3) then
           allocate(d1_rbf(n_rbf,n_rbf))
        else
           allocate(d_rbf(n_rbf,n_rbf,n_ine_loc_1,n_coll))
        endif

     endif

     ! Gyroaverage arrays
     allocate(z_gyro(-n_x/2:n_x/2,-n_x/2:n_x/2-1))

     ! G0a,G2a in Technical manual
     allocate(w_gyro0(n_stack,-m_gyro:m_gyro-i_gyro,n_x,n_nek_loc_1,n_gk))
     allocate(w_gyro2(n_stack,-m_gyro:m_gyro-i_gyro,n_x,n_nek_loc_1,n_gk))

     if (n_field == 3) then
        ! G1a,G3a in Technical manual
        allocate(w_gyro1(n_stack,-m_gyro:m_gyro-i_gyro,n_x,n_nek_loc_1,n_gk))
        allocate(w_gyro3(n_stack,-m_gyro:m_gyro-i_gyro,n_x,n_nek_loc_1,n_gk))
     endif

  else 

     deallocate(omega_d1)
     deallocate(omega_dr)
     deallocate(omega_star)
     deallocate(v_theta)
     deallocate(v_para)
     deallocate(v_perp)

     deallocate(nek_n)
     deallocate(nek_e)
     deallocate(nek_k)
     deallocate(ine_i)
     deallocate(ine_n)
     deallocate(ine_e)
     deallocate(ki_k)
     deallocate(ki_i)

     deallocate(cs_blend_prime)
     deallocate(cs_blend)
     deallocate(c_blend)

     if (allocated(o_f)) deallocate(o_f)
     if (allocated(o_fv)) deallocate(o_fv)
     if (allocated(o_advect)) deallocate(o_advect)
     if (allocated(imp)) deallocate(imp)

     if (allocated(d_rbf)) deallocate(d_rbf) 
     if (allocated(d1_rbf)) deallocate(d1_rbf) 

     deallocate(z_gyro)
     deallocate(w_gyro0)
     deallocate(w_gyro2)
     if (allocated(w_gyro1)) deallocate(w_gyro1)
     if (allocated(w_gyro3)) deallocate(w_gyro3)

  endif

end subroutine gyro_alloc_distrib
