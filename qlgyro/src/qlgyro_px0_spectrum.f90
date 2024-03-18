      SUBROUTINE qlgyro_px0_spectrum
!

      USE qlgyro_globals
      IMPLICIT NONE

      integer :: i
      
      if (.not.allocated(px0_spectrum)) allocate(px0_spectrum(n_px0))

      ! Set px0 values to scan through (currently 0 -> 0.5)
      if (n_px0 .eq. 1) then
         px0_spectrum(1) = 0.0
      else if (px0grid_model .eq. 2) then
         do i=1, n_px0 - 2
            px0_spectrum(i) = ((i-1)*0.1) * 0.5
         end do
         px0_spectrum(n_px0-1) = 0.25
         px0_spectrum(n_px0) = 0.5
      else if (px0grid_model .eq. 3) then
         do i=1, n_px0 - 2
            px0_spectrum(i) = (i-1) * 0.016
         end do
         px0_spectrum(n_px0-1) = 0.2
         px0_spectrum(n_px0) = 0.5
      else
         do i=1, n_px0
            px0_spectrum(i) = (i-1.0) / (n_px0-1.0) * 0.5
         end do
      end if
      
    END SUBROUTINE qlgyro_px0_spectrum
