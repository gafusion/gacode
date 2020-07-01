subroutine fluxfit_f_model(t,r,z)

  use fluxfit_globals

  implicit none

  integer :: n

  real, intent(in) :: t
  real, intent(inout) :: r,z

  select case (model)

  case (1)

     ! Elevation, shift, kappa, delta, squareness

     r = c(3)+c(1)*cos(t+asin(c(5))*sin(t))
     z = c(2)+c(4)*c(1)*sin(t+c(6)*sin(2*t))

  case (2)

     r = ar(0)/2.0
     z = az(0)/2.0
     do n=1,ns
        r = r+ar(n)*cos(n*t)+br(n)*sin(n*t)
        z = z+az(n)*cos(n*t)+bz(n)*sin(n*t)
     enddo

  end select

end subroutine fluxfit_f_model

