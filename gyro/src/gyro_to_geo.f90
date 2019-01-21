!-----------------------------------------------------
! gyro_to_geo.f90
!
! PURPOSE:
!  Pack array and call GEO geometry package.
!-----------------------------------------------------

subroutine gyro_to_geo(i0)

  use gyro_globals
  use geo
  
  !--------------------------------
  implicit none
  !
  integer, intent(in) :: i0
  !--------------------------------

  GEO_signb_in   = 1.0
  GEO_rmin_in    = r_s(i0)
  GEO_rmaj_in    = rmaj_s(i0)
  GEO_drmaj_in   = drmaj_s(i0)
  GEO_zmag_in    = zmag_s(i0)
  GEO_dzmag_in   = dzmag_s(i0)
  GEO_q_in       = q_s(i0)
  GEO_s_in       = shat_s(i0)
  GEO_kappa_in   = kappa_s(i0)
  GEO_s_kappa_in = s_kappa_s(i0)
  GEO_delta_in   = delta_s(i0)
  GEO_s_delta_in = s_delta_s(i0)
  GEO_zeta_in    = zeta_s(i0)
  GEO_s_zeta_in  = s_zeta_s(i0)
  GEO_beta_star_in = beta_star_s(i0)

  GEO_fourier_in(:,:) = a_fourier_geo_s(:,:,i0)

end subroutine gyro_to_geo
