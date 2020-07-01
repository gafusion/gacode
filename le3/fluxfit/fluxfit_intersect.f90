subroutine fluxfit_intersect(x,y,n,y0,xp,xm)

  implicit none

  ! Input variables
  integer, intent(in) :: n
  real, dimension(n), intent(in) :: x,y
  real, intent(inout) :: y0
  real, intent(inout) :: xp,xm

  ! Internal variables
  integer :: i,i0
  real :: x_c
  real :: x0
  real :: x1,x2,x3
  real :: y1,y2,y3
  real :: d1,d2,d3
  real :: d12,d13,d23
  real :: xa,xb
  real :: ya,yb


  x_c = 0.5*(minval(x)+maxval(x))

  do i=1,n-1

     ya = y(i)
     yb = y(i+1)
     xa = x(i)
     xb = x(i+1)

     if ((ya-y0)*(yb-y0) < 0.0) then

        ! Linear approximation (unused)
        ! x0 = (y0-ya)*xb/(yb-ya)+(y0-yb)*xa/(ya-yb)

        if (abs(ya-y0) < abs(yb-y0)) then
           if (i == 1) then
              i0 = n-1
           else
              i0 = i-1
           endif
           y1 = y(i0)
           y2 = ya
           y3 = yb
           x1 = x(i0)
           x2 = xa
           x3 = xb
        else
           if (i == n-1) then
              i0 = 2
           else
              i0 = i+2
           endif
           y1 = ya
           y2 = yb
           y3 = y(i0)
           x1 = xa
           x2 = xb
           x3 = x(i0)
        endif

        d12 = y1-y2
        d13 = y1-y3
        d23 = y2-y3
        d3 = d13*d23
        d2 = -d12*d23
        d1 = d12*d13

        ! Quadratic approximation:
 
        x0 = (y0-y1)*(y0-y2)*x3/d3+(y0-y1)*(y0-y3)*x2/d2+(y0-y2)*(y0-y3)*x1/d1

        if (x0 > x_c) then
           xp = x0
        else
           xm = x0
        endif

     endif

  enddo

end subroutine fluxfit_intersect
