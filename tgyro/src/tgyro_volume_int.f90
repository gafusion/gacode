!---------------------------------------------------
! tgyro_volume_int.f90
!
! Integrate source to obtain integrated power:
!
! s    : source (input)
! volp : derivative of volume, V' (input)
! p    : power (output; erg/s if source in erg/cm^3/s)
!
!        r
!        /
! p(r) = | dx V'(x) s(x)
!        /
!        0
!
!---------------------------------------------------

subroutine tgyro_volume_int(s,p)

  use tgyro_globals, only : n_r,r,volp,use_trap

  implicit none

  integer :: i

  real :: x1,x2,x3
  real :: f1,f2,f3

  real, dimension(n_r), intent(in) :: s 
  real, dimension(n_r), intent(inout) :: p 

  ! Integrated power in erg/s

  if (use_trap == 1) then

     ! Integration using trapezoidal rule

     p(1) = 0.0
     do i=2,n_r
        p(i) = p(i-1)+0.5*(s(i-1)*volp(i-1)+s(i)*volp(i))*(r(i)-r(i-1))
     enddo

  else

     ! Integration using quadratic polynomial interpolation

     ! NOTE: this was introduced to significantly improve 
     !       accuracy for sparse (n_r < 10) grids.

     p(1) = 0.0
     do i=2,n_r

        if (i == 2) then

           ! Use points 1-3 for integral from 1-2

           x1 = r(i-1)
           x2 = r(i)
           x3 = r(i+1)

           f1 = s(i-1)*volp(i-1)
           f2 = s(i)*volp(i)
           f3 = s(i+1)*volp(i+1)

           p(i) = p(i-1)+(x2-x1)*(3*x3-x2-2*x1)*f1/6/(x3-x1) &
                + (x2-x1)*(3*x3-2*x2-x1)*f2/6/(x3-x2) &
                - (x2-x1)**3*f3/6/(x3-x1)/(x3-x2)
        else

           ! Use points 1-3 for integral from 2-3

           x1 = r(i-2)
           x2 = r(i-1)
           x3 = r(i)

           f1 = s(i-2)*volp(i-2) 
           f2 = s(i-1)*volp(i-1)
           f3 = s(i)*volp(i)

           p(i) = p(i-1)+(x3-x2)*(2*x3+x2-3*x1)*f3/6/(x3-x1) &
                + (x3-x2)*(x3+2*x2-3*x1)*f2/6/(x2-x1) &
                - (x3-x2)**3*f1/6/(x2-x1)/(x3-x1)
        endif

     enddo

  endif

end subroutine tgyro_volume_int

subroutine tgyro_volume_ave(f,r,volp,fave,n)

  implicit none

  integer :: i
  integer, intent(in) :: n
  real, intent(in), dimension(n) :: f,r,volp
  real, intent(inout) :: fave
  real :: vol

  fave = 0.0
  vol  = 0.0
  do i=2,n
     fave = fave+0.5*(f(i-1)*volp(i-1)+f(i)*volp(i))*(r(i)-r(i-1))
     vol  = vol+0.5*(volp(i-1)+volp(i))*(r(i)-r(i-1))
  enddo
  fave = fave/vol

end subroutine tgyro_volume_ave

!--------------------------------------------------------------------
! math_scaleint.f90
!
! PURPOSE:
!  Given the linear (z=f') or negative logarithmic (z=-f'/f) gradient 
!  of a function f, find the backward integral:
!
! type=log
!                    x*
!                    /
!  f(x) = f(x*) exp[ | z(s) ds ]
!                    /
!                    x
!
! type=lin
!                 x*
!                 /
!  f(x) = f(x*) - | z(s) ds
!                 /
!                 x
!
!
!  In other words, We obtain f by analytic integration of the piecewise
!  linear gradient function z(x).
!
! NOTES:
!        z(n) : (input) gradients
!        x(n) : (input) radii
!           n : (input) dimension of input arrays
!          fb : (input) f(x*) = boundary condition 
!          x0 : (input) desired location of profile
!           f : (output) f(x0)
!---------------------------------------------------------------------

subroutine math_scaleint(z,x,n,fb,x0,f,type)

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: z(n),x(n)
  real, intent(in) :: x0,fb
  real, intent(inout) :: f
  character(len=3), intent(in) :: type

  integer :: i,i0
  real :: xs,dx

  xs = x0

  if (xs > x(n)) then
     print *,'WARNING: (tgyro_scaleint) x0 > x(n)'
     xs = x(n)
     f = fb
     return
  endif

  if (xs < 0.0) then
     print *,'WARNING: (tgyro_scaleint) x0 < 0'
     xs = 0.0
  endif

  ! Find i0 such that x(i0-1) < xs < x(i0)
  do i0=n,2,-1
     if (xs >= x(i0-1)) exit
  enddo

  ! Integrate up to nearest gridpoint with trapezoidal integration
  ! then perform final subinterval with analytic integration.
  ! NOTE: the analytic method reduces to trapezoidal rule at meshpoints

  dx = x(i0)-x(i0-1)
  f  = fb

  if (type == 'log') then
     do i=n,i0+1,-1
        f = f*exp(0.5*(z(i)+z(i-1))*(x(i)-x(i-1)))
     enddo
     ! Analytic integral of linear function
     f = f*exp(0.5/dx*(z(i0)*(dx**2-(xs-x(i0-1))**2)+z(i0-1)*(xs-x(i0))**2))
  else
     do i=n,i0+1,-1
        f = f-0.5*(z(i)+z(i-1))*(x(i)-x(i-1))
     enddo
     ! Analytic integral of linear function
     f = f-0.5/dx*(z(i0)*(dx**2-(xs-x(i0-1))**2)+z(i0-1)*(xs-x(i0))**2)
  endif

end subroutine math_scaleint

!--------------------------------------------------------------------
! math_scaleintv.f90
!
! PURPOSE:
!  Given the linear (z=f') or negative logarithmic (z=-f'/f) gradient 
!  of a function f, find the backward integral:
!
! type=log
!                    x*
!                    /
!  f(x) = f(x*) exp[ | z(s) ds ]
!                    /
!                    x
!
! type=lin
!                 x*
!                 /
!  f(x) = f(x*) - | z(s) ds
!                 /
!                 x
!
!
!  In other words, We obtain f by analytic integration of the piecewise
!  linear gradient function z(x).
!
! NOTES:
!        z(n) : (input) gradients
!        x(n) : (input) radii
!           n : (input) dimension of input arrays
!           f : (inout) f defined at n 
!---------------------------------------------------------------------

subroutine math_scaleintv(z,x,n,f,type)

  implicit none

  integer, intent(in) :: n
  real, intent(in) :: z(n),x(n)
  real, intent(inout), dimension(n) :: f
  character(len=3), intent(in) :: type

  integer :: i

  if (type == 'log') then
     do i=n,2,-1
        f(i-1) = f(i)*exp(0.5*(z(i)+z(i-1))*(x(i)-x(i-1)))
     enddo
  else
     do i=n,2,-1
        f(i-1) = f(i)-0.5*(z(i)+z(i-1))*(x(i)-x(i-1))
     enddo
  endif

end subroutine math_scaleintv
