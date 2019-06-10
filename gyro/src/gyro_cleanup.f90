!-----------------------------------------------------------
! gyro_cleanup.f90
!
! PURPOSE:
!  Deallocate and clean up.
!-----------------------------------------------------------

subroutine gyro_cleanup

  use mpi
  use gyro_globals
  use gyro_interface

  implicit none

  ! No need for cleanup in test mode
  if (gyrotest_flag == 1) return

  if (allocated(r)) then

     ! GYRO arrays
     call gyro_alloc_profile_sim(0)
     call gyro_alloc_big(0)
     call gyro_alloc_orbit(0)
     call gyro_alloc_velocity(0)
     call gyro_alloc_distrib(0)
     if (nonlinear_flag == 1) call gyro_alloc_nl(0)

     ! Sparse field arrays 
     if (allocated(m_poisson)) deallocate(m_poisson)
     if (allocated(indx_poisson)) deallocate(indx_poisson)
     if (allocated(m_ampere)) deallocate(m_ampere)
     if (allocated(indx_ampere)) deallocate(indx_ampere)
     if (allocated(m_maxwell)) deallocate(m_maxwell)
     if (allocated(indx_maxwell)) deallocate(indx_maxwell)
     if (allocated(m_poissonaperp)) deallocate(m_poissonaperp)
     if (allocated(indx_poissonaperp)) deallocate(indx_poissonaperp)

     ! Blending coefficients
     call BLEND_cleanup

     call MPI_COMM_FREE(NEW_COMM_1,i_err)
     call MPI_COMM_FREE(NEW_COMM_2,i_err)

  endif

end subroutine gyro_cleanup
