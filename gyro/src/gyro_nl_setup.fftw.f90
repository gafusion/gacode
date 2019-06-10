!-----------------------------------------------------
! gyro_nl_setup.fftw.f90
!
! PURPOSE:
! Allocate (and define) selected arrays for
! use with direct and FFT methods.
!
! NOTES:
! FFTW2 routine.
!-------------------------------------------------------

subroutine gyro_nl_setup

  use gyro_globals
  use gyro_nl_private

  !------------------------------------
  implicit none
  !------------------------------------

  include 'fftw_f77.i'

  !----------------------------------------------------
  do nn=-n_max,n_max
     i_p(nn) = abs(nn)+1
     if (nn < 0) then
        n_p(nn) = -n(i_p(nn))
     else
        n_p(nn) = n(i_p(nn))
     endif
  enddo
  !
  ! Nonlinear coupling coefficient:
  !
  c_nl_i(:) = -rhos_norm*q_s(:)/(r_s(:)*b_unit_s(:))
  !----------------------------------------------------

  !----------------------------------------------------
  ! FFTW2 arrays and plans
  !
  n_max_d = (3*n_max)/2+1
  n_fft = 2*n_max_d+2
  !
  allocate( v_fft(0:n_fft-1,6*n_x))
  allocate(vt_fft(0:n_fft-1,6*n_x))

  call rfftw_f77_create_plan(plan_b,&
       n_fft,&
       FFTW_BACKWARD,&
       FFTW_MEASURE)

  call rfftw_f77_create_plan(plan_f,&
       n_fft,&
       FFTW_FORWARD,&
       FFTW_MEASURE)
  !----------------------------------------------------

end subroutine gyro_nl_setup
