!---------------------------------------------------------
! cgyro_freq.f90
!
! PURPOSE:
!  Comput estimates of linear growth rates
!---------------------------------------------------------

subroutine cgyro_freq

  use cgyro_globals

  implicit none

  real :: total_weight,dfr,dfi
  real, dimension(nc) :: mode_weight
  complex, dimension(nc) :: freq_loc

  if (i_time == 0 .or. n == 0) then

     freq = 0.0
     freq_err = 0.0

  else

     if (px0 < 0) then

        !--------------------------------------------------
        ! Standard method: sum all wavenumbers at a given n
        !--------------------------------------------------

        ! Use potential to compute frequency
        mode_weight(:) = abs(field_old(1,:))

     else

        !--------------------------------------------------
        ! Alternate method: use fixed wavenumber
        !--------------------------------------------------

        mode_weight = 0.0        

        do ic=1,nc
           if (px(ir_c(ic)) == n*(px0+box_size)) then
              mode_weight(ic) = abs(field_old(1,ic))
           endif
           if (px(ir_c(ic)) == n*(px0-box_size)) then
              mode_weight(ic) = abs(field_old(1,ic))
           endif
        enddo

     endif

     ! Define local frequencies
     freq_loc(:) = (i_c/delta_t)*log(field_old(1,:)/field_old2(1,:))

     total_weight = sum(mode_weight(:))

     freq = sum(freq_loc(:)*mode_weight(:))/total_weight

     ! Fractional Frequency Error
     dfr = sum(abs(real(freq_loc(:)-freq))*mode_weight(:))
     dfi = sum(abs(aimag(freq_loc(:)-freq))*mode_weight(:))

     freq_err = (dfr+i_c*dfi)/total_weight/abs(freq)

     if (n_toroidal == 1 .and. abs(freq_err) < freq_tol) signal=1

  endif


end subroutine cgyro_freq
