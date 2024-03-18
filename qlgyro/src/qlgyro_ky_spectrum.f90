      SUBROUTINE qlgyro_ky_spectrum(kymin)
!
      USE tglf_dimensions
      USE tglf_global
      USE tglf_species
      USE tglf_kyspectrum
      USE tglf_interface
      IMPLICIT NONE

      real, intent(in) :: kymin
      real, dimension(12) :: n_tor
      integer :: i_ky

      call get_ky_spectrum

      tglf_nky_in = nky
      if (.not. allocated(tglf_ky_spectrum_out)) then
         allocate(tglf_ky_spectrum_out(nky), tglf_dky_spectrum_out(nky))
      end if

      tglf_ky_spectrum_out = ky_spectrum(:nky)
      tglf_dky_spectrum_out = dky_spectrum(:nky)

    END SUBROUTINE qlgyro_ky_spectrum
