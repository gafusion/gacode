!
!--------------------------------------------------------------
!
      SUBROUTINE qlgyro_sat1
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
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_kyspectrum
      USE tglf_xgrid
      USE tglf_sgrid
      USE tglf_interface
      USE qlgyro_globals

      ! CGYRO eigenvalues and fluxes - need to set these equal when
      ! n_PX0 = 1
      tglf_eigenvalue_spectrum_out = qlgyro_eigenvalue_spectrum_out
      tglf_flux_spectrum_out = qlgyro_flux_spectrum_out

      eigenvalue_spectrum_out(:, :nky, :1) = tglf_eigenvalue_spectrum_out
      QL_flux_spectrum_out(:, :ns, :, :nky, :1) = tglf_flux_spectrum_out

      call get_multiscale_spectrum

      tglf_flux_spectrum_out = flux_spectrum_out(:, :ns, :, :nky, :1)
      tglf_field_spectrum_out = field_spectrum_out

    
    END SUBROUTINE qlgyro_sat1
