!-------------------------------------------------------
! bound_deriv.f90
!
! PURPOSE:
!  Compute the finite-difference derivative of a 
!  function (on a grid which may be unequally-spaced)
!  using the Lagrange interpolating polynomial.  
!-------------------------------------------------------

subroutine bound_deriv(df,f,r,n)

  implicit none

  integer, intent(in) :: n

  real, intent(inout), dimension(n) :: df
  real, intent(in), dimension(n) :: f
  real, intent(in),dimension(n) :: r

  real :: r1,r2,r3,ra
  real :: f1,f2,f3

  integer :: i

  do i=1,n

     if (i == 1) then

        ! Left boundary

        ra = r(1)

        r1 = r(1)
        r2 = r(2)
        r3 = r(3)
        f1 = f(1)
        f2 = f(2)
        f3 = f(3)

     else if (i == n) then

        ! Right boundary

        ra = r(n)

        r1 = r(n-2)
        r2 = r(n-1)
        r3 = r(n)
        f1 = f(n-2)
        f2 = f(n-1)
        f3 = f(n)

     else

        ! Interior

        ra = r(i)

        r1 = r(i-1)
        r2 = r(i)
        r3 = r(i+1)
        f1 = f(i-1)
        f2 = f(i)
        f3 = f(i+1)

     endif

     ! Derivative of Lagrange interpolating polynomial:

     df(i) = ((ra-r1)+(ra-r2))/(r3-r1)/(r3-r2)*f3 &
          + ((ra-r1)+(ra-r3))/(r2-r1)/(r2-r3)*f2 &
          + ((ra-r2)+(ra-r3))/(r1-r2)/(r1-r3)*f1

  enddo

end subroutine bound_deriv
