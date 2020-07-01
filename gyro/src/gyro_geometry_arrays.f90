!-----------------------------------------------------
! gyro_geometry_arrays.f90 [caller BigScience]
!
! PURPOSE:
!  Definition of arrays which require calls to 
!  general geometry library.
!-----------------------------------------------------

subroutine gyro_geometry_arrays

  use gyro_globals
  use geo

  !--------------------------------------
  implicit none
  integer :: p
  !--------------------------------------


  do i=1,n_x

     call gyro_to_geo(i)
     call geo_interp(n_lambda*n_stack,theta_t(i,:,:),.true.)

     p = 0
     do m=1,n_stack
        do k=1,n_lambda

           p = p+1

           b0_t(i,k,m)      = GEO_b(p)
           g_theta_t(i,k,m) = GEO_g_theta(p)
           grad_r_t(i,k,m)  = GEO_grad_r(p)
           qrat_t(i,k,m)    = GEO_gq(p)

           if (geo_gradbcurv_flag == 1) then

              ! No explicit pressure effect.  

              ! This is also known as the grad B -> curvature 
              ! rule, where we remark that curvature has no 
              ! pressure dependence (but grad B does!).

              cos_t(i,k,m)   = GEO_gcos1(p)
              cos_p_t(i,k,m) = 0.0

           else 

              ! Retain full grad-B and curvature effects.

              cos_t(i,k,m)   = GEO_gcos1(p)+GEO_gcos2(p)
              cos_p_t(i,k,m) = -GEO_gcos2(p)

           endif

           bt_t(i,k,m)   = GEO_bt(p)
           bp_t(i,k,m)   = GEO_bp(p)
           bigr_t(i,k,m) = GEO_bigr(p)

           captheta_t(i,k,m) = GEO_captheta(p)
           sin_t(i,k,m)      = GEO_gsin(p)

           usin_t(i,k,m) = GEO_usin(p)
           ucos_t(i,k,m) = GEO_ucos(p)

        enddo ! k
     enddo ! m

     call geo_interp(n_theta_plot,theta_plot,.false.)
     do j=1,n_theta_plot
        b0_plot(i,j)      = GEO_b(j)
        g_theta_plot(i,j) = GEO_g_theta(j)
     enddo ! j

  enddo !i

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_geometry_arrays done]'
  endif

end subroutine gyro_geometry_arrays
