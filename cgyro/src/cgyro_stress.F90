!---------------------------------------------------------------------------
! cgyro_stress.f90
!
! PURPOSE:
!  Compute Reynolds/Maxwell stresses from the fields as a
!  function of (kx,theta,species,ky,field).
!
! NOTES:
!
!  The stress terms are just the nonlinear term evaluated for each
!  field individually
!
!  To trigger output of the stresses
!    set STRESS_PRINT_FLAG=1

!---------------------------------------------------------------------------

subroutine cgyro_stress

  use timer_lib
  use mpi
  use cgyro_globals
  use cgyro_nl_comm

  implicit none

  integer :: i_field

  ! Set up stress for each field
  ! Not sure if this can be with async given the NL calc set up
  do i_field=1, n_field
     ! Set up first half h_x
     call cgyro_nl_fftw_comm1_async

     ! Set up fields
     call cgyro_nl_fftw_comm2_async_stress(i_field)

     ! Set up second half h_x (if split)
     call cgyro_nl_fftw_comm3_async

     ! Calculate nonlinear term
     call cgyro_nl_fftw
     call cgyro_nl_fftw_comm_test

     ! Set stress array from nonlinear output
     call cgyro_nl_fftw_comm1_r_stress(i_field)
     call cgyro_nl_fftw_comm_test
  end do
#if defined(OMPGPU)
!$omp target update from(stress)
#elif defined(_OPENACC)
!$acc update host(stress)
#endif

end subroutine cgyro_stress
