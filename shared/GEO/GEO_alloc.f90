!------------------------------------------------
! GEO_alloc.f90
!
! PURPOSE:
!  Allocate (or deallocate) arrays.  This must 
!  be called before GEO_do.
!
! NOTES:
!
!  - require n_theta > 8
!  - should choose n_theta to be odd
!
!  flag=1: allocate
!  flag=0: deallocate
!------------------------------------------------

subroutine GEO_alloc(flag)

  use GEO_interface

  !-------------------------------------------
  implicit none
  !
  integer, intent(in) :: flag
  !-------------------------------------------

  if (flag == 1) then
     if (GEO_ntheta_in < 9) then 
        print *,'Need more point in GEO_alloc.' 
        stop
     endif
     if (modulo(GEO_ntheta_in,2) == 0) then 
        print *,'GEO_ntheta_in must be odd in GEO_alloc.' 
        stop
     endif
     allocate(GEOV_b(GEO_ntheta_in))
     allocate(GEOV_bt(GEO_ntheta_in))
     allocate(GEOV_bp(GEO_ntheta_in))
     allocate(GEOV_dbdt(GEO_ntheta_in))
     allocate(GEOV_dbdt2(GEO_ntheta_in))
     allocate(GEOV_gsin(GEO_ntheta_in))
     allocate(GEOV_gcos1(GEO_ntheta_in))
     allocate(GEOV_gcos2(GEO_ntheta_in))
     allocate(GEOV_grad_r(GEO_ntheta_in))
     allocate(GEOV_jac_r(GEO_ntheta_in))
     allocate(GEOV_g_theta(GEO_ntheta_in))
     allocate(GEOV_gq(GEO_ntheta_in))
     allocate(GEOV_captheta(GEO_ntheta_in))
     allocate(GEOV_nu(GEO_ntheta_in))
     allocate(GEOV_theta(GEO_ntheta_in))
     allocate(GEOV_l_r(GEO_ntheta_in))
     allocate(GEOV_l_t(GEO_ntheta_in))
     allocate(GEOV_nsin(GEO_ntheta_in))
     allocate(GEOV_usin(GEO_ntheta_in))
     allocate(GEOV_ucos(GEO_ntheta_in))
     allocate(GEOV_bigr(GEO_ntheta_in))
     allocate(GEOV_bigr_r(GEO_ntheta_in))
     allocate(GEOV_bigr_t(GEO_ntheta_in))
     allocate(GEOV_theta_nc(GEO_ntheta_in))
     allocate(GEOV_theta_s(GEO_ntheta_in))
     allocate(GEOV_chi2(GEO_ntheta_in))

     allocate(GEO_fourier_in(8,0:GEO_nfourier_in))

  else
     deallocate(GEOV_b)
     deallocate(GEOV_bt)
     deallocate(GEOV_bp)
     deallocate(GEOV_dbdt)
     deallocate(GEOV_dbdt2)
     deallocate(GEOV_gsin)
     deallocate(GEOV_gcos1)
     deallocate(GEOV_gcos2)
     deallocate(GEOV_jac_r)
     deallocate(GEOV_grad_r)
     deallocate(GEOV_g_theta)
     deallocate(GEOV_gq)
     deallocate(GEOV_captheta)
     deallocate(GEOV_nu)
     deallocate(GEOV_theta)
     deallocate(GEOV_l_r)
     deallocate(GEOV_l_t)
     deallocate(GEOV_nsin)
     deallocate(GEOV_usin)
     deallocate(GEOV_ucos)
     deallocate(GEOV_bigr)
     deallocate(GEOV_bigr_r)
     deallocate(GEOV_bigr_t)
     deallocate(GEOV_theta_nc)
     deallocate(GEOV_theta_s)
     deallocate(GEOV_chi2)

     deallocate(GEO_fourier_in)

  endif

end subroutine GEO_alloc
