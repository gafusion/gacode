!---------------------------------------------------------
! cgyro_freq.f90
!
! PURPOSE:
!  Compute estimates of linear eigenvalue w where:
!
!                   -iwt
!            phi ~ e
!---------------------------------------------------------

subroutine cgyro_freq

  use cgyro_globals
  use cgyro_io

  implicit none

  real :: mw
  real :: total_weight,dfr,dfi
  real, dimension(nc) :: mode_weight
  integer :: itor
  complex, dimension(nc) :: freq_loc
  complex :: fl,myfr,df,total_weighted_freq

  if (i_time == 0) then

    freq(:) = 0.0
    freq_err(:) = 0.0

  else

    !--------------------------------------------------
    ! Standard method: sum all wavenumbers at a given n
    !--------------------------------------------------

    do itor=nt1,nt2

     total_weight = 0.0
     total_weighted_freq = (0.0,0.0)
     do ic=1,nc
        ! Use potential to compute frequency
        mw = abs(field_old(1,ic,itor))
        mode_weight(ic) =  mw
        total_weight = total_weight + mw
        ! Define local frequencies
        if (abs(field_old(1,ic,itor)) > 1e-12 .and. abs(field_old2(1,ic,itor)) > 1e-12) then
           fl = (i_c/delta_t)*log(field_old(1,ic,itor)/field_old2(1,ic,itor))
        else
           fl = 0.0
        endif
        freq_loc(ic) = fl
        total_weighted_freq = total_weighted_freq + fl*mw
     enddo

     myfr = total_weighted_freq/total_weight
     freq(itor) = myfr
     
     ! Fractional Frequency Error
     dfr = 0.0
     dfi = 0.0
     do ic=1,nc
        mw = mode_weight(ic)
        df = freq_loc(ic)-myfr
        dfr = dfr + abs(real(df))*mw
        dfi = dfi + abs(aimag(df))*mw
     enddo

     if (abs(myfr) > 1e-12) then
        ! Trap a division-by-zero error and halt
        freq_err(itor) = (dfr+i_c*dfi)/total_weight/abs(myfr)
     else
        call cgyro_error('Underflow in calculation of frequency error')
     endif
     
    enddo

    if (n_toroidal == 1 .and. abs(freq_err(nt1)) < freq_tol) signal = 1

  endif


end subroutine cgyro_freq
