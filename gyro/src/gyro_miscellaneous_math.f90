!------------------------------------------------
! invert_p32.f90
!
! PURPOSE:
!  Calculate x given s, where:
!
!                         x
!                  2      /          -t
! s = p32(x) = --------   | sqrt(t) e   dt
!              sqrt(pi)   /
!                         0
!
! Speed is not esential, so we use bisection 
! for robsutness.
!
! REVISIONS:
! 26 Sept 01: jc 
!  Documented
! 25 May 05: jc
!  Moved into gyro/math subdir and renamed.
!------------------------------------------------

subroutine invert_p32(s,x,x_max)

  implicit none

  real, parameter :: eps_x = 1e-10

  real, intent(in) :: s
  real, intent(in) :: x_max
  real, intent(inout) :: x

  real :: dx
  real :: s1
  real :: s2
  real :: residual

  real, external :: p32


  residual = 1.0

  dx = x_max/2.0-eps_x
  x  = x_max/2.0

  do while (residual > eps_x)

     s1 = p32(x)
     s2 = p32(x+dx)

     dx = 0.5*dx

     if (s1 <= s .and. s2 >= s) then
        x = x+dx
     else
        x = x-dx
     endif

     residual = abs(s2-s1)

  enddo

end subroutine invert_p32

!-----------------------------------------------------
! p32.f90
!
! PURPOSE:
!  Compute the incomplete Gamma (P) function.
!
! NOTES:
!                     x
!              2      /          -t
! p32(x) = --------   | sqrt(t) e   dt == P(3/2,x)
!          sqrt(pi)   /
!                     0
! 
! REVISIONS
! 10 Jan 01: jeff.candy@gat.com
!  Created.
!---------------------------------------------

real function p32(x)

  use math_constants

  implicit none

  real, intent(in) :: x
 
  p32 = erf(sqrt(x))-sqrt(x)*exp(-x)/(0.5*sqrt(pi))

end function p32

!--------------------------------------------------------------
! pascal.f90 
!
! PURPOSE:
!  Compute coefficients in Pascal's triangle, with added 
!  alternating sign.
!
! NOTES:
!  Equivalently, these are the coefficients of the 
!  centered n-th derivative operator.  Below, we have 
!  supressed the alternating sign.
!
!  n    s(1)  s(2)  s(3) ...
! ---   ---------------------------------------------------------
!  3    1.0   2.0   1.0
!  4    1.0   3.0   3.0   1.0
!  5    1.0   4.0   6.0   4.0   1.0
!  6    1.0   5.0  10.0  10.0   5.0   1.0
!  7    1.0   6.0  15.0  20.0  15.0   6.0   1.0
!  8    1.0   7.0  21.0  35.0  35.0  21.0   7.0   1.0
!  9    1.0   8.0  28.0  56.0  70.0  56.0  28.0   8.0   1.0
! 
! REVISIONS
! 06 Feb 01: jeff.candy@gat.com
!  Created. 
!--------------------------------------------------------------

subroutine pascal(n,s)

  implicit none

  integer, intent(in) :: n
  integer :: i
  integer :: m
  real, dimension(n) :: s
  real, dimension(n) :: s0

  s = 0 

  if (n < 3) then
     print *,'PASCAL: error in n'
     stop
  endif

  if (n == 3) then

     s(1) = 1.0  
     s(2) = -2.0
     s(3) = 1.0

  else

     s0(1) = 1.0  
     s0(2) = 2.0
     s0(3) = 1.0

     do m=4,n

        s(1) = 1.0
        do i=2,m-1
           s(i) = s0(i)+s0(i-1)
        enddo
        s(m) = 1.0

        s0 = s

     enddo

     ! Ensure middle element is negative for n odd.
 
     do m=1,n
        s(m) = (-1)**((n+1)/2-m-1)*s(m)
     enddo

  endif

end subroutine pascal
