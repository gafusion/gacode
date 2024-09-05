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
      integer :: mg_ky_grid=0

      if (tglf_kygrid_model_in .eq. -1) then
         kygrid_model_in = 0
         mg_ky_grid = 1
      end if

      call get_ky_spectrum

      tglf_nky_in = nky
      if (.not. allocated(tglf_ky_spectrum_out)) then
         allocate(tglf_ky_spectrum_out(nky), tglf_dky_spectrum_out(nky))
      end if

      if (mg_ky_grid .eq. 1) then

         n_tor = (/2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 120.0, 140./)
         ky_spectrum = 0.0
         ky_spectrum(1) = kymin * n_tor(1)
         dky_spectrum(1) = kymin * n_tor(1)
         do i_ky=2,nky
            ky_spectrum(i_ky) = kymin * n_tor(i_ky)
            dky_spectrum(i_ky) = ky_spectrum(i_ky) - ky_spectrum(i_ky -1)
         end do
      end if

      tglf_ky_spectrum_out = ky_spectrum(:nky)
      tglf_dky_spectrum_out = dky_spectrum(:nky)

    END SUBROUTINE qlgyro_ky_spectrum
