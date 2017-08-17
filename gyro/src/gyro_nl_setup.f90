!-----------------------------------------------------
! gyro_nl_setup.fftw3.f90
!
! PURPOSE:
!  Allocate (and define) selected arrays for 
!  use with direct and FFT methods.
!
! NOTES:
! FFTW3 routine.  
!-------------------------------------------------------

subroutine gyro_nl_setup

  use gyro_globals
  use gyro_nl_private
  use ompdata

  !------------------------------------
  implicit none
  !------------------------------------

  include 'fftw3.f'
  integer :: ierr

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
  !
  ! ... and add the nonuniform grid effect:
  !
  c_nl_i(:) = dr_eodr(:)*c_nl_i(:)
  !----------------------------------------------------

  !----------------------------------------------------
  ! FFTW3 arrays and plans
  !
  n_max_d = (3*n_max)/2+1
  n_fft   = 2*n_max_d+2
  !
  allocate( v_fft3(0:n_fft-1,n_x,6))
  allocate(vt_fft3(0:n_fft-1,n_x,6))

  call dfftw_init_threads(ierr)

  call dfftw_plan_with_nthreads(n_omp)

  call dfftw_plan_many_r2r(plan_b, 1, n_fft, 6*n_x,   &
                           v_fft3,  0, 1, n_fft,       &
                           vt_fft3, 0, 1, n_fft,       &
                           FFTW_HC2R, FFTW_MEASURE)

  call dfftw_plan_many_r2r(plan_f, 1, n_fft, 6*n_x,   &
                           vt_fft3, 0, 1, n_fft,       &
                           v_fft3,  0, 1, n_fft,       &
                           FFTW_R2HC, FFTW_MEASURE)
  !----------------------------------------------------

end subroutine gyro_nl_setup
