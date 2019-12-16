!----------------------------------------------------
! gyro_alloc_big.f90
!
! PURPOSE:
!  Allocation for most of the large arrays (typically 
!  physical fields) needed for time stepping.  Other 
!  arrays are allocated as they are needed.
!
! NOTES:
!
! The tentative "priority" for ordering indices in
! array is:
!
!  j         m         i     p_nek_loc   is
!  n_blend   n_stack   n_x   n_nek_loc   n_spec
!------------------------------------------------------

subroutine gyro_alloc_big(flag)

  use gyro_globals
  use gyro_pointers

  implicit none

  integer, intent(in) :: flag

  if (flag == 1 .and. allocated(field_blend)) then
     if (i_proc == 0) then
        print *,'WARNING: (GYRO) already allocated arrays in gyro_alloc_big'
     endif
     return
  endif
  if (flag == 0 .and. .not.allocated(field_blend)) then
     if (i_proc == 0) then
        print *,'WARNING: (GYRO) cannot deallocate arrays in gyro_alloc_big'
     endif
     return
  endif

  if (flag == 1) then

     allocate(field_blend(n_blend,n_x,n_field))
     allocate(field_blend_old(n_blend,n_x,n_field))
     allocate(field_blend_old2(n_blend,n_x,n_field))
     allocate(field_blend_dot(n_blend,n_x,n_field))
     allocate(field_tau(n_stack,n_x,n_nek_loc_1,n_field))
     allocate(field_tau_old(n_stack,n_x,n_nek_loc_1,n_field))
     allocate(field_tau_old2(n_stack,n_x,n_nek_loc_1,n_field))
     allocate(phi_squared(n_x))
     allocate(field_fluxave(n_x,3))
     allocate(ave_phi(2,n_field))
     allocate(h_err(n_stack,n_x,n_nek_loc_1,n_kinetic))

     if (.not.allocated(h)) allocate(h(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(h_old(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(h_0(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(h_cap(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(h_cap_old(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(h_cap_old2(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(h_cap_dot(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(rhs(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(rhs_dr(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(rhs_dt(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(entropy(n_kinetic,n_entro))
    
     if (collision_flag == 1) then
        allocate(f_coll(n_stack,n_x,n_nek_loc_1))
        allocate(fb_coll(n_stack,n_x,n_nek_loc_1))
        allocate(h_c(n_stack,n_lambda,n_ine_loc_1)) 
     endif

     allocate(h_tran(nv1_SSUB,msplit_SSUB,n_n,n_kinetic))
     allocate(gyro_h(n_stack,n_x,n_nek_loc_1,n_kinetic))

     if (n_field == 3) then
        allocate(gyro_h_aperp(n_stack,n_x,n_nek_loc_1,n_kinetic))
     endif

     allocate(gyro_uv(n_stack,n_x,n_nek_loc_1,n_kinetic,n_field))
     allocate(kyro_uv(n_stack,n_x,n_nek_loc_1,n_kinetic,n_field))
     allocate(gyro_uv_old2(n_stack,n_x,n_nek_loc_1,n_kinetic,n_field))
     allocate(gyro_uv_old(n_stack,n_x,n_nek_loc_1,n_kinetic,n_field))
     allocate(gyro_uv_dot(n_stack,n_x,n_nek_loc_1,n_kinetic,n_field))
     allocate(gyro_u(n_stack,n_x,n_nek_loc_1,n_kinetic))
     allocate(gyro_u_tran(nv1_SSUB,msplit_SSUB,n_n,n_kinetic))
     allocate(phi(n_theta_int,n_x,n_field))
     allocate(vel_sum_p(n_blend,n_x))
     allocate(vel_sum_a(n_blend,n_x))
     allocate(vel_sum_aperp(n_blend,n_x))
     allocate(phi_plot(n_theta_plot,n_x,n_field+eparallel_plot_flag))

     allocate(moments_plot(n_theta_plot,n_x,n_kinetic,3))
     allocate(moments_zero_plot(n_x,n_kinetic,6))

     allocate(kxkyspec(n_x))
     allocate(k_perp_squared(n_n))

     allocate(nonlinear_flux(n_x,n_kinetic,n_field,p_moment))
     allocate(nonlinear_flux_momparts(n_kinetic,3))
     allocate(nonlinear_flux_excparts(n_kinetic,2))
     allocate(gbflux_i(n_kinetic,n_field,p_moment,n_x))
     allocate(gbflux(n_kinetic,n_field,p_moment))
     allocate(gbflux_mom(n_kinetic,3))
     allocate(gbflux_exc(n_kinetic,2))
     allocate(gbflux_n(n_kinetic,n_field,p_moment))
     allocate(gbflux_vec(n_kinetic,n_field,p_moment,n_x))

     allocate(time_error(n_kinetic))
     allocate(omega_linear(n_n,2))

     !------------------------------------------------------------
     ! Source-related arrays
     !
     allocate(h0_eq(n_kinetic,n_energy,n_x))
     if (source_method == 3) then
        allocate(h_source(n_stack,n_x,n_nek_loc_1,n_kinetic))
     endif
      !------------------------------------------------------------

  else 

     deallocate(field_blend)
     deallocate(field_blend_old)
     deallocate(field_blend_old2)
     deallocate(field_blend_dot)
     deallocate(field_tau)
     deallocate(field_tau_old)
     deallocate(field_tau_old2)
     deallocate(phi_squared)
     deallocate(field_fluxave)
     deallocate(ave_phi)
     deallocate(h_err)

     deallocate(h)
     deallocate(h_old)
     deallocate(h_0)
     deallocate(h_cap)
     deallocate(h_cap_old)
     deallocate(h_cap_old2)
     deallocate(h_cap_dot)
     deallocate(rhs)
     deallocate(rhs_dr)
     deallocate(rhs_dt)
     deallocate(entropy)

     if (allocated(f_coll)) deallocate(f_coll)
     if (allocated(fb_coll)) deallocate(fb_coll)
     if (allocated(h_c)) deallocate(h_c) 

     deallocate(h_tran)
     deallocate(gyro_h)
     if (allocated(gyro_h_aperp)) deallocate(gyro_h_aperp)

     deallocate(gyro_uv)
     if (allocated(kyro_uv)) deallocate(kyro_uv)
     deallocate(gyro_uv_old2)
     deallocate(gyro_uv_old)
     deallocate(gyro_uv_dot)
     deallocate(gyro_u)
     deallocate(gyro_u_tran)
     deallocate(phi)
     deallocate(vel_sum_p)
     deallocate(vel_sum_a)
     deallocate(vel_sum_aperp)
     deallocate(phi_plot)

     deallocate(moments_plot)
     deallocate(moments_zero_plot)

     deallocate(kxkyspec)
     deallocate(k_perp_squared)

     deallocate(nonlinear_flux)  
     deallocate(nonlinear_flux_momparts)
     deallocate(nonlinear_flux_excparts)

     deallocate(gbflux_i)
     deallocate(gbflux)
     deallocate(gbflux_mom)
     deallocate(gbflux_exc)
     deallocate(gbflux_n)
     deallocate(gbflux_vec)

     deallocate(time_error)
     deallocate(omega_linear)
     deallocate(h0_eq)
     if (allocated(h_source)) deallocate(h_source)

  endif

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[allocate_big done]' 
  endif

end subroutine gyro_alloc_big
