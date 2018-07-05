subroutine vgen_getgeo

  use vgen_globals
  use EXPRO_interface
  use geo

  implicit none

  integer :: i, j
  real :: r_min
  integer :: n_theta=4
  real, dimension(4) :: theta=(/ 0.0, 1.57079632, 3.14159265, -1.57079632/)

  GEO_nfourier_in = EXPRO_nfourier
  GEO_signb_in    = EXPRO_signb

  open(unit=1,file='out.vgen.geo',status='replace')
  r_min = EXPRO_rmin(EXPRO_n_exp)

  do i=2,EXPRO_n_exp
     
     ! Parameters to be passed to GEO library   
     !
     ! NOTE: dp/dr set to zero without loss of generality.
     ! 
     GEO_rmin_in      = EXPRO_rmin(i)/r_min
     GEO_rmaj_in      = EXPRO_rmaj(i)/r_min
     GEO_drmaj_in     = EXPRO_drmaj(i)
     GEO_zmag_in      = EXPRO_zmag(i)/r_min
     GEO_dzmag_in     = EXPRO_dzmag(i)
     GEO_q_in         = EXPRO_q(i)
     GEO_s_in         = EXPRO_s(i)
     GEO_kappa_in     = EXPRO_kappa(i)
     GEO_s_kappa_in   = EXPRO_skappa(i)
     GEO_delta_in     = EXPRO_delta(i)
     GEO_s_delta_in   = EXPRO_sdelta(i)
     GEO_zeta_in      = EXPRO_zeta(i)
     GEO_s_zeta_in    = EXPRO_szeta(i)
     GEO_beta_star_in = 0.0
     !
     if (EXPRO_ctrl_numeq_flag == 0) then
        ! Call GEO with model shape
        GEO_model_in = 0
     else
        ! Call GEO with general (numerical) shape
        GEO_model_in = 1
        GEO_fourier_in(1:4,:) = EXPRO_geo(:,:,i)/r_min
        GEO_fourier_in(5:8,:) = EXPRO_dgeo(:,:,i)
     endif

     call geo_interp(n_theta,theta,.true.)
     
     write(1,'(e16.8)',advance='no') EXPRO_rho(i)
     do j=1,n_theta
        write(1,'(e16.8)',advance='no') theta(j)
        write(1,'(e16.8)',advance='no') GEO_b(j)  * EXPRO_bunit(i)
        write(1,'(e16.8)',advance='no') GEO_bp(j) * EXPRO_bunit(i)
        write(1,'(e16.8)',advance='no') GEO_bt(j) * EXPRO_bunit(i)
        write(1,'(e16.8)',advance='no') GEO_bigr(j) * r_min
     enddo
     write(1,*)
  enddo

  ! Clean up
  close(1)

end subroutine vgen_getgeo
