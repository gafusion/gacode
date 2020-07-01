subroutine fluxfit_moments(r,z,n,s)

  implicit none

  ! Input variables
  integer, intent(in) :: n
  real, intent(in), dimension(n) :: r,z
  real, intent(inout), dimension(3) :: s

  ! Work variables
  integer :: i,ip
  real :: dz,dr
  real :: r0,z0

  s(:) = 0.0
  
  do i=1,n-1

     ip = i+1

     dz = z(ip)-z(i)
     dr = r(ip)-r(i)
     r0 = 0.5*(r(ip)+r(i))
     z0 = 0.5*(z(ip)+z(i))
     
     ! A = Int R dZ = -Int Z dR

     s(1) = s(1)+dz*r0    ! Area
     !s(1) = s(1)-dr*z0    ! Area
     s(2) = s(2)-dr*z0*r0 ! Centroid (R)
     s(3) = s(3)+dz*r0*z0 ! Centroid (Z)
     !s(3) = s(3)-dr*z0*z0 ! Centroid (Z)
 
  enddo

end subroutine fluxfit_moments
