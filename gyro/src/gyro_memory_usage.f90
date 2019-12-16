!------------------------------------------------
! gyro_memory_usage.f90
!
! PURPOSE:
!  Manually sum amount of memory used to allocate
!  large arrays.
!------------------------------------------------

subroutine gyro_memory_usage(data_file,io)

  use gyro_globals
  use gyro_collision_private
  use gyro_pointers

  !---------------------
  implicit none
  !---------------------

  integer, intent(in) :: io
  character (len=*), intent(in) :: data_file

  select case (output_flag)

  case (1)

     if (i_proc == 0) then

        open(unit=io,file=data_file,status='replace')

        !---------------------------------------------
        ! gyro_alloc_distrib
        !---------------------------------------------

        write(io,*) 'gyro_alloc_distrib:'
        write(io,*)

        call alloc_add(io,size(omega_d1),8,'omega_d1')
        call alloc_add(io,size(omega_dr),8,'omega_dr')
        call alloc_add(io,size(omega_star),8,'omega_star')

        call alloc_add(io,size(v_theta),8,'v_theta')
        call alloc_add(io,size(v_theta),8,'v_theta')
        call alloc_add(io,size(v_para),8,'v_para')

        call alloc_add(io,size(nek_n),4,'nek_n')
        call alloc_add(io,size(nek_e),4,'nek_e')
        call alloc_add(io,size(nek_k),4,'nek_k')
        call alloc_add(io,size(ine_i),4,'ine_i')
        call alloc_add(io,size(ine_n),4,'ine_n')
        call alloc_add(io,size(ine_e),4,'ine_e')

        call alloc_add(io,size(cs_blend),16,'cs_blend')
        call alloc_add(io,size(c_blend),16,'c_blend')
        call alloc_add(io,size(cs_blend_prime),16,'cs_blend_prime')

        if (electron_method == 2) then
           call alloc_add(io,size(o_f),16,'o_f')
           call alloc_add(io,size(o_fv),16,'o_fv')
           call alloc_add(io,size(o_advect),16,'o_advect')
           call alloc_add(io,size(imp),16,'imp')
        endif

        if (collision_flag == 1) then
           if (linsolve_method == 3) then
              call alloc_add(io,size(d1_rbf),8,'d1_rbf')
           else
              call alloc_add(io,size(d_rbf),8,'d_rbf')
           endif
        endif

        call alloc_add(io,size(z_gyro),16,'z_gyro')

        call alloc_add(io,size(w_gyro0),16,'w_gyro0')
        call alloc_add(io,size(w_gyro2),16,'w_gyro2')
        if (n_field == 3) then
           call alloc_add(io,size(w_gyro1),16,'w_gyro1')
           call alloc_add(io,size(w_gyro3),16,'w_gyro3')
        endif

        !---------------------------------------------
        ! gyro_alloc_big
        !---------------------------------------------

        write(io,*)
        write(io,*) 'gyro_alloc_big:'
        write(io,*)

        call alloc_add(io,size(field_blend),16,'field_blend')
        call alloc_add(io,size(field_blend_old),16,'field_blend_old')
        call alloc_add(io,size(field_blend_dot),16,'field_blend_dot')
        call alloc_add(io,size(phi_squared),8,'phi_squared')
        call alloc_add(io,size(field_fluxave),8,'field_fluxave')
        call alloc_add(io,size(h_err),16,'h_err')
        call alloc_add(io,size(h_0),16,'h_0')
        call alloc_add(io,size(h_old),16,'h_old')
        call alloc_add(io,size(h),16,'h')
        call alloc_add(io,size(h_cap),16,'h_cap')
        call alloc_add(io,size(h_source),16,'h_source')
        call alloc_add(io,size(rhs),16,'rhs')
        call alloc_add(io,size(rhs_dr),16,'rhs_dr')
        call alloc_add(io,size(rhs_dt),16,'rhs_dt')

        if (collision_flag == 1) then
           call alloc_add(io,size(f_coll),16,'f_coll')
           call alloc_add(io,size(fb_coll),16,'fb_coll')
           call alloc_add(io,size(h_c),16,'h_c') 
        endif

        call alloc_add(io,size(h_tran),16,'h_tran') 
        call alloc_add(io,size(gyro_h),16,'gyro_h')
        if (n_field == 3) then
           call alloc_add(io,size(gyro_h_aperp),16,'gyro_h_aperp')
        endif
        call alloc_add(io,size(field_tau),16,'field_tau') 

        call alloc_add(io,size(gyro_uv),16,'gyro_uv')
        call alloc_add(io,size(gyro_uv_old2),16,'gyro_uv_old2')
        call alloc_add(io,size(gyro_uv_old),16,'gyro_uv_old')
        call alloc_add(io,size(gyro_uv_dot),16,'gyro_uv_dot')
        call alloc_add(io,size(gyro_u),16,'gyro_u')
        call alloc_add(io,size(gyro_u_tran),16,'gyro_u_tran')
        call alloc_add(io,size(phi),16,'phi')
        call alloc_add(io,size(vel_sum_p),16,'vel_sum_p')
        call alloc_add(io,size(vel_sum_a),16,'vel_sum_a')
        call alloc_add(io,size(vel_sum_aperp),16,'vel_sum_aperp')
        call alloc_add(io,size(phi_plot),16,'phi_plot')

        call alloc_add(io,size(moments_plot),16,'moments_plot')
        call alloc_add(io,size(moments_zero_plot),8,'moments_zero_plot')

        call alloc_add(io,size(kxkyspec),8,'kxkyspec')
 
        call alloc_add(io,size(nonlinear_flux),8,'nonlinear_flux')

        call alloc_add(io,size(gbflux_i),8,'gbflux_i')
        call alloc_add(io,size(gbflux),8,'gbflux')
        call alloc_add(io,size(gbflux_n),8,'gbflux_n')

        if (transport_method > 1) then
           call alloc_add(io,size(gbflux_vec),8,'gbflux_vec')
        endif

        call alloc_add(io,size(time_error),8,'time_error')
        call alloc_add(io,size(entropy),8,'entropy')

        call alloc_add(io,size(h0_eq),8,'h0_eq')

        !---------------------------------------------
        ! gyro_alloc_orbit
        !---------------------------------------------

        write(io,*)
        write(io,*) 'gyro_alloc_orbit:'
        write(io,*)

        call alloc_add(io,size(theta_t),8,'theta_t')
        call alloc_add(io,size(tau),8,'tau')
        call alloc_add(io,size(m_cyc),4,'m_cyc')
        call alloc_add(io,size(p_cyc),16,'p_cyc')
        call alloc_add(io,size(m_phys),4,'m_phys')
        call alloc_add(io,size(p_phys),4,'p_phys')
        call alloc_add(io,size(omega),8,'omega')
        call alloc_add(io,size(m_map),4,'m_map')
        call alloc_add(io,size(theta_int),8,'theta_int')
        call alloc_add(io,size(theta_plot),8,'theta_plot')

        call alloc_add(io,size(b0_t),8,'b0_t')
        call alloc_add(io,size(g_theta_t),8,'g_theta_t')
        call alloc_add(io,size(grad_r_t),8,'grad_r_t')
        call alloc_add(io,size(qrat_t),8,'qrat_t')
        call alloc_add(io,size(cos_t),8,'cos_t')
        call alloc_add(io,size(cos_p_t),8,'cos_p_t')  
        call alloc_add(io,size(captheta_t),8,'captheta_t')  
        call alloc_add(io,size(sin_t),8,'sin_t') 
        call alloc_add(io,size(bt_t),8,'bt_t')
        call alloc_add(io,size(bp_t),8,'bp_t')
        call alloc_add(io,size(bigr_t),8,'bigr_t')
        call alloc_add(io,size(b0_plot),8,'b0_plot')
        call alloc_add(io,size(g_theta_plot),8,'g_theta_plot')  

        call alloc_add(io,size(c_fluxave),8,'c_fluxav')
        call alloc_add(io,size(ff_mm),16,'ff_mm')
        call alloc_add(io,size(ff_mm_piv),4,'ff_mm_piv')
        call alloc_add(io,size(ff2_mm),16,'ff2_mm')
        call alloc_add(io,size(ff2_mm_piv),4,'ff2_mm_piv')
        call alloc_add(io,size(blend_plot),16,'blend_plo')
        call alloc_add(io,size(blend_prime_plot),16,'blend_prime_plot')

        if (collision_flag == 1) then
           call alloc_add(io,size(nu_total),8,'nu_total')
           call alloc_add(io,size(xi),8,'xi')
        endif

        if (n_field > 1) then
           call alloc_add(io,size(coll_vel),16,'coll_vel')
        endif
        if (n_field > 2) then
           call alloc_add(io,size(coll_vel_perp1),16,'coll_vel_perp1')
           call alloc_add(io,size(coll_vel_perp2),16,'coll_vel_perp2')
        endif

        !---------------------------------------------
        ! Sparse field arrays (see gyro_cleanup)
        !---------------------------------------------

        write(io,*)
        write(io,*) 'Sparse field arrays (no alloc block):'
        write(io,*)

        if (allocated(m_poisson)) call alloc_add(io,size(m_poisson),16,'m_poisson')
        if (allocated(indx_poisson)) call alloc_add(io,size(indx_poisson),4,'indx_poisson')
        if (allocated(m_ampere)) call alloc_add(io,size(m_ampere),16,'m_ampere')
        if (allocated(indx_ampere)) call alloc_add(io,size(indx_ampere),4,'indx_ampere')
        if (allocated(m_maxwell)) call alloc_add(io,size(m_maxwell),16,'m_maxwell')
        if (allocated(indx_maxwell)) call alloc_add(io,size(indx_maxwell),4,'indx_maxwell')
        if (allocated(m_poissonaperp)) call alloc_add(io,size(m_poissonaperp),16,'m_poissonaperp')
        if (allocated(indx_poissonaperp)) call alloc_add(io,size(indx_poissonaperp),4,'indx_poissonaperp')

        write(io,*) '---------------------'
        write(io,'(f9.3,a,3x,a)') total_memory/1e6,' MB'

        close(io)

     endif

  end select

  if (i_proc == 0 .and. debug_flag == 1) then
     print *,'[gyro_memory_usage done]'
  endif

end subroutine gyro_memory_usage

!------------------------------------------------
! alloc_add.f90
!
! PURPOSE:
!  Primitive allocation addition routine.
!------------------------------------------------

subroutine alloc_add(io,n_size,bytes,name)

  use gyro_globals

  implicit none
  !
  integer, intent(in) :: io
  integer, intent(in) :: n_size
  integer, intent(in) :: bytes
  character (len=*), intent(in) :: name
  !
  real :: this_memory

  select case (output_flag)

  case (1)

     if (i_proc == 0) then

        this_memory  = 1.0*n_size*bytes
        total_memory = total_memory+this_memory 

        write(io,10) this_memory/1e6,' MB',name

     endif

  end select

10 format(f9.3,a,3x,a)

end subroutine alloc_add
