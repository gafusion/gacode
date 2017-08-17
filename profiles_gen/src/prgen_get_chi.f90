subroutine prgen_get_chi(nx,q,kappa,rmin,psi,rho,bref,arho)

  implicit none

  integer, intent(in) :: nx
  real, intent(in), dimension(nx) :: q
  real, intent(in), dimension(nx) :: kappa
  real, intent(in), dimension(nx) :: rmin
  real, intent(in), dimension(nx) :: psi
  real, intent(inout), dimension(nx) :: rho
  real, intent(inout) :: bref
  real, intent(inout) :: arho

  real, dimension(:), allocatable :: chi_t
  integer :: i 

  !---------------------------------------------------------
  ! Compute rho by integrating d(chi_t) = q d(psi_p)
  ! using the trapezoidal rule
  !
  ! NOTE: chi_t/bref)
  !
  allocate(chi_t(nx))  
  !
  chi_t(1) = 0.0
  do i=2,nx
     chi_t(i) = chi_t(i-1) + 0.5*(q(i)+q(i-1))*(psi(i)-psi(i-1))
  enddo
  !
  ! Even though bref is arbitrary, choose physically sensible values 
  ! rather than, say, bref=1.  This is exact in the limit B is constant 
  ! and the flux surfaces are elliptical with fixed kappa.
  ! 
  bref = 2*chi_t(nx)/(kappa(nx)*rmin(nx)**2)
  arho = sqrt(kappa(nx))*rmin(nx)
  !
  ! Normalized root of chi_t:
  !
  rho(:) = sqrt(chi_t(:)/chi_t(nx))
  deallocate(chi_t)
  !---------------------------------------------------------

end subroutine prgen_get_chi
