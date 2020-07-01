!-----------------------------------------------------
! gyro_to_geo.f90
!
! PURPOSE:
!  Pack array and call geo geometry package.
!-----------------------------------------------------

subroutine gyro_to_geo(i0)

  use gyro_globals
  use geo
  
  !--------------------------------
  implicit none
  !
  integer, intent(in) :: i0
  !--------------------------------

  geo_signb_in   = 1.0
  geo_rmin_in    = r_s(i0)
  geo_rmaj_in    = rmaj_s(i0)
  geo_drmaj_in   = drmaj_s(i0)
  geo_zmag_in    = zmag_s(i0)
  geo_dzmag_in   = dzmag_s(i0)
  geo_q_in       = q_s(i0)
  geo_s_in       = shat_s(i0)
  geo_kappa_in   = kappa_s(i0)
  geo_s_kappa_in = s_kappa_s(i0)
  geo_delta_in   = delta_s(i0)
  geo_s_delta_in = s_delta_s(i0)
  geo_zeta_in    = zeta_s(i0)
  geo_s_zeta_in  = s_zeta_s(i0)
  geo_beta_star_in = beta_star_s(i0)

  ! New geometry not implemented in GYRO
  geo_shape_cos0_in = 0.0
  geo_shape_s_cos0_in = 0.0
  geo_shape_cos1_in = 0.0
  geo_shape_s_cos1_in = 0.0
  geo_shape_cos2_in = 0.0
  geo_shape_s_cos2_in = 0.0
  geo_shape_cos3_in = 0.0
  geo_shape_s_cos3_in = 0.0
  geo_shape_sin3_in = 0.0
  geo_shape_s_sin3_in = 0.0

  geo_fourier_in(:,:) = a_fourier_geo_s(:,:,i0)

end subroutine gyro_to_geo
