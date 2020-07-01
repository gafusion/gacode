! High accuracy estimate of minimum x given initial guess as 
! minimum of array x.

subroutine fluxfit_minmax(x,y,n,x0,y_x0,mode)

  implicit none

  ! Input variables
  integer, intent(in) :: n
  real, dimension(n), intent(in) :: x,y
  real, intent(inout) :: x0,y_x0
  character (len=3), intent(in) :: mode

  ! Internal variables
  integer :: i,im,ip
  real :: x1,x2,x3
  real :: y1,y2,y3
  real :: y0
  real :: d1,d2,d3
  real :: d12,d13,d23

  if (trim(mode) == 'min') then
     ! Minimum x
     i = minloc(x,1)
  else
     ! Maximum x
     i = maxloc(x,1)
  endif

  ip = modulo(i,n-1)+1
  im = modulo(i-2,n-1)+1
  x1 = x(im) 
  x2 = x(i)
  x3 = x(ip)
  y1 = y(im)
  y2 = y(i)
  y3 = y(ip)

  if (abs(x1-x2) <= epsilon(0.0) .and. abs(x2-x3) <= epsilon(0.0)) then
     print '(10(1pe13.6,1x))',x(:)
     print '(a)','ERROR: (fluxfit_minmax) x1=x2=x3 error.'
     stop
  endif

  d12 = y1-y2
  d13 = y1-y3
  d23 = y2-y3
  d3 = d13*d23
  d2 = -d12*d23
  d1 = d12*d13
  y0 = (y1+y2)*x3/d3+(y1+y3)*x2/d2+(y2+y3)*x1/d1
  y0 = y0/(2*x3/d3+2*x2/d2+2*x1/d1)
  x0 = (y0-y1)*(y0-y2)*x3/d3+(y0-y1)*(y0-y3)*x2/d2+(y0-y2)*(y0-y3)*x1/d1

  y_x0 = y0

end subroutine fluxfit_minmax
