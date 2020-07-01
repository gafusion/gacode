subroutine prgen_get_chi(nx,q,psi,rho,torfluxa)

  implicit none

  integer, intent(in) :: nx
  real, intent(in), dimension(nx) :: q
  real, intent(in), dimension(nx) :: psi
  real, intent(inout), dimension(nx) :: rho
  real, intent(inout) :: torfluxa

  real, dimension(:), allocatable :: chi_t
  integer :: i 

  !---------------------------------------------------------
  ! Compute rho by integrating d(chi_t) = q d(psi_p)
  ! using the trapezoidal rule
  !
  allocate(chi_t(nx))  
  !
  chi_t(1) = 0.0
  do i=2,nx
     chi_t(i) = chi_t(i-1) + 0.5*(q(i)+q(i-1))*(psi(i)-psi(i-1))
  enddo
  !
  torfluxa = chi_t(nx)
  !
  ! Normalized root of chi_t:
  !
  rho(:) = sqrt(chi_t(:)/chi_t(nx))
  ! Prevent any roundoff
  rho(nx) = 1.0
  deallocate(chi_t)
  !---------------------------------------------------------

end subroutine prgen_get_chi
