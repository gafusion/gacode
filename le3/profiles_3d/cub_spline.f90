!---------------------------------------------------------------
! cub_spline.f90
!
! PURPOSE:
!  Take known points x(1:n),y(1:n) and perform cubic spline
!  interpolation at the node points xi(1:ni) to give yi(1:ni).
!
!  Use 'natural' cubic spline interpolation, where natural 
!  means the second derivative is zero at the boundaries.
!
!  INPUT  : x(1:n),y(1:n),n,xi(1:ni),ni
!
!  OUTPUT : yi(1:ni)
!
!  The only requirements are 
!  
!  1. Ordered data: x(i+1) > x(i) and xi(i+1) > xi(i).
!  2. Bounded data: xi(1) < x(1) and xi(ni) > x(n).
!----------------------------------------------------------------

subroutine cub_spline(x,y,n,xi,yi,ni)

  !-------------------------------------------------------------
  implicit none
  !
  integer :: i,ii
  real :: x0 
  !
  integer, intent(in) :: n
  real, intent(in), dimension(n) :: x,y
  !
  integer, intent(in) :: ni
  real, dimension(ni) :: xi,yi
  !
  ! LAPACK working variables
  !
  integer :: info
  integer, dimension(n) :: ipiv 
  real, dimension(n) :: c,z
  real, dimension(n-1) :: zl,zu,h
  real, dimension(n-2) :: zu2
  real, dimension(n-1) :: b,d
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Check to see that interpolated point is inside data interval
  !
  if (xi(ni) > x(n)) then
     print *,'ERROR: (cub_spline) Data above upper bound'
     print *,'xi(ni) > x(n)',xi(ni),x(n) 
  endif 
  if (xi(1) < x(1)) then
     print *,'ERROR: (cub_spline) Data below lower bound'
     print *,'xi(1) < x(1)',xi(1),x(1) 
  endif 
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Define coefficients of spline matrix 
  !
  do i=1,n-1 
     h(i)  = x(i+1)-x(i)
     zl(i) = h(i)
     zu(i) = h(i)
  enddo
  zl(n-1) = 0.0
  zu(1)   = 0.0

  z(1) = 1.0
  c(1) = 0.0
  do i=2,n-1
     z(i) = 2.0*(h(i-1)+h(i))
     c(i) = 3.0*((y(i+1)-y(i))/h(i)-(y(i)-y(i-1))/h(i-1))
  enddo
  z(n) = 1.0
  c(n) = 0.0
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Solve the system using LAPACK
  !
  call DGTTRF(n,zl,z,zu,zu2,ipiv,info)
  call DGTTRS('N',n,1,zl,z,zu,zu2,ipiv,c,n,info)
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Find remaining polynomial coefficients:
  !
  c(n) = 0.0
  do i=1,n-1
     b(i) = (y(i+1)-y(i))/h(i)-h(i)*(2.0*c(i)+c(i+1))/3.0
     d(i) = (c(i+1)-c(i))/(3.0*h(i))
  enddo
  !-------------------------------------------------------------

  !print *,c(:)
  
  !-------------------------------------------------------------
  ! Using known polynomial coefficients, perform interpolation.
  !
  !  S(x) = y(i) + b(i) [x-x(i)] + c(i) [x-x(i)]^2 
  !                                        + d(i) [x-x(i)]^3
  !
  i  = 1
  ii = 1
  do while (ii <= ni)
     x0 = xi(ii)
     if (x0 <= x(i+1)) then
        yi(ii) = y(i)+(x0-x(i))*(b(i)+(x0-x(i))*(c(i)+(x0-x(i))*d(i)))
        ii = ii+1
     else
        i = i+1
     endif
  enddo
  !-------------------------------------------------------------

end subroutine cub_spline

