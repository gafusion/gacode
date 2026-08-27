!
!--------------------------------------------------------------
!
      SUBROUTINE qlgyro_sat_mg
!
!***********************************************************************
!  questions  should be addressed to
!  Gary Staebler 858-455-3466  or email: gary.staebler@gat.com
!***********************************************************************
!
!  TGLF multiscale saturation rule: sat_rule_in = 1 option
!  April 8, 2016: G.M. Staebler, J. Candy, N. T. Howard, and C. Holland
!                  Physics of Plasmas, 23 (2016) 062518
!  June 22, 2017: Retuned after coding error to Laplacian terms in Ampere 
!                 and Poisson equations was fixed
!
!***********************************************************************
!

        USE qlgyro_globals
        USE tglf_interface
        

      implicit none
      integer :: i_kypx0, px0_max_index, i_ky_local, i_px0_local, i_species, i_field, i_moment
      real, dimension(n_ky, n_px0) :: ql_metric_ky_px0, species_sum_flux
      real, dimension(n_ky) :: ql_metric_ky
      real :: ql_metric
      real :: px0_max, max_gam
      real :: facnorm, b, a

      call calculate_ql_metric(ql_metric_ky_px0)

      ! Normalise fluxes to total heat flux
      species_sum_flux = sum(sum(qlgyro_flux_spectrum_out(2, :, :, :, :), dim=2), dim=1)
      
      do i_ky_local=1, n_ky
         do i_px0_local=1, n_px0
            if (species_sum_flux(i_ky_local, i_px0_local) .eq. 0.0) then
               ql_metric_ky_px0(i_ky_local, i_px0_local) = 0.0
               species_sum_flux(i_ky_local, i_px0_local) = 1.0
            end if
         end do
      end do

      ! Re-scale fluxes
      do i_moment=1, 3
         do i_species=1, n_species
            do i_field=1, n_field
               qlgyro_flux_spectrum_out(i_moment, i_species, i_field, :, :) = &
               qlgyro_flux_spectrum_out(i_moment, i_species, i_field, :, :) &
               * ql_metric_ky_px0 / species_sum_flux
            end do
         end do
      end do

      ! DO PX0 average here!
      if (n_px0 .eq. 1) then
         ! CGYRO flux
         tglf_flux_spectrum_out = qlgyro_flux_spectrum_out
         ql_metric_ky(:) = ql_metric_ky_px0(:, 1)
      else
         do i_ky_local=1, n_ky
            max_gam = maxval(qlgyro_eigenvalue_spectrum_out(1, i_ky_local, :))

            if (max_gam < abs(gamma_exb) / 10) then
               ! Set all to 0 
               tglf_flux_spectrum_out(:, :, :, i_ky_local, 1) = qlgyro_flux_spectrum_out(:, :, :, i_ky_local, 1) * 0.0
               ql_metric_ky(i_ky_local) = 0.0
            else
               px0_max = min(abs(gamma_exb / (s * max_gam * 2 * pi)), 0.5)
               px0_max = max(px0_max, px0_spectrum(2))

               do i_px0_local = 1, n_px0
                  if (px0_spectrum(i_px0_local) .le. px0_max) then
                     px0_max_index = i_px0_local
                  else
                     exit
                  end if
               end do
               
               call trapezoid_flux(px0_spectrum, px0_max_index, qlgyro_flux_spectrum_out(:, :, :, i_ky_local, :) &
               , tglf_flux_spectrum_out(:, :, :, i_ky_local, 1))
               call trapezoid_metric_px0(px0_spectrum, px0_max_index, ql_metric_ky_px0(i_ky_local, :), ql_metric_ky(i_ky_local))
            end if
         end do 
      end if
 
      call trapezoid_ky(tglf_ky_spectrum_out, ql_metric_ky, ql_metric)

      a = 2.1

      ! Convert from using rho_s(B0) to rho_s(Bunit) by (bgs2**-2)
      facnorm = bgs2**(-2)

      b = exp(3.0)
      tglf_flux_spectrum_out = tglf_flux_spectrum_out * facnorm * b * ql_metric ** (a-1)
      tglf_eigenvalue_spectrum_out = qlgyro_eigenvalue_spectrum_out

    END SUBROUTINE qlgyro_sat_mg


    subroutine calculate_ql_metric(ql_metric_local)

      use qlgyro_globals
      implicit none

      integer :: i_field, i_thetab
      
      real, dimension(3, n_thetab, n_ky, n_px0) :: field_sq, numer, denom
      real, dimension(3, n_ky, n_px0) :: max_field, numer_trapz, denom_trapz
      real, intent(out), dimension(n_ky, n_px0) :: ql_metric_local
      
      field_sq = abs(qlgyro_field_spectrum_out)**2

      numer = 0.0
      denom = 0.0
      do i_thetab=1, n_thetab
         do i_field=1, n_field
            numer(i_field, i_thetab, :, :) = field_sq(i_field, i_thetab, :, :) * qlgyro_jacobian(i_thetab, :, :)
            denom(i_field, i_thetab, :, :) = numer(i_field, i_thetab, :, :) * qlgyro_k_perp(i_thetab, :, :)**2
         end do
      end do
      
      max_field = maxval(field_sq, dim=2)
      
      call trapezoid_field(qlgyro_theta_ballooning, numer, numer_trapz)
      call trapezoid_field(qlgyro_theta_ballooning, denom, denom_trapz)

      ql_metric_local = 0.0
      do i_field=1, n_field
         ql_metric_local = ql_metric_local + (qlgyro_eigenvalue_spectrum_out(1, :, :) * numer_trapz(i_field, :, :)&
         / denom_trapz(i_field, :, :) * (max_field(i_field, :, :) / max_field(1, :, :)) ** 0.5)

      end do

    end subroutine calculate_ql_metric


    subroutine trapezoid_field(thetab, integrand, integral)

      use qlgyro_globals
      implicit none

      real, dimension(n_thetab), intent(in) :: thetab
      real, dimension(3, n_thetab, n_ky, n_px0), intent(in) :: integrand
      real, dimension(3, n_ky, n_px0), intent(out) :: integral
      integer :: i_thetab
      real :: dthetab
      
      integral = 0.0

      do i_thetab=1, n_thetab - 1
         dthetab = thetab(i_thetab+1) - thetab(i_thetab)
         integral = integral + dthetab * (integrand(:, i_thetab+1, :, :) + integrand(:, i_thetab, :, :)) / 2
      end do

    end subroutine trapezoid_field
      
    
    subroutine trapezoid_flux(px0, px0_max, integrand, integral)

      use qlgyro_globals
      implicit none

      real, dimension(n_px0), intent(in) :: px0
      real, dimension(5, n_species, n_field, n_px0), intent(in) :: integrand
      real, dimension(5, n_species, n_field), intent(out) :: integral
      integer, intent(in) :: px0_max
      integer :: i_px0_local
      real :: dpx0
      
      integral = 0.0

      do i_px0_local=1, px0_max - 1
         dpx0 = px0(i_px0_local+1) - px0(i_px0_local)
         integral = integral + dpx0 * (integrand(:, :, :, i_px0_local+1) + integrand(:, :, :, i_px0_local)) / 2
      end do

      integral = integral / px0(px0_max)

    end subroutine trapezoid_flux
    
    subroutine trapezoid_metric_px0(px0, px0_max, integrand, integral)

      use qlgyro_globals
      implicit none

      real, dimension(n_px0), intent(in) :: px0
      real, dimension(n_px0), intent(in) :: integrand
      real, intent(out) :: integral
      integer, intent(in) :: px0_max
      integer :: i_px0_local
      real :: dpx0
      
      integral = 0.0
      
      do i_px0_local=1, px0_max - 1
         dpx0 = px0(i_px0_local+1) - px0(i_px0_local)
         integral = integral + dpx0 * (integrand(i_px0_local+1) + integrand(i_px0_local)) / 2
      end do

      integral = integral / px0(px0_max)

    end subroutine trapezoid_metric_px0
    
    subroutine trapezoid_ky(ky_local, integrand, integral)

      use qlgyro_globals
      implicit none

      real, dimension(n_ky), intent(in) :: ky_local
      real, dimension(n_ky), intent(in) :: integrand
      real, intent(out) :: integral
      integer :: i_ky_local
      real :: dky
      
      integral = 0.0

      do i_ky_local=1, n_ky - 1
         dky = ky_local(i_ky_local+1) - ky_local(i_ky_local)
         integral = integral + dky * (integrand(i_ky_local+1) + integrand(i_ky_local)) / 2
      end do

    end subroutine trapezoid_ky
    
