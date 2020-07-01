subroutine fluxfit_fourier

  use fluxfit_globals

  implicit none

  integer :: i,n
  real :: y
  real :: dtheta

  ar(:) = 0.0
  br(:) = 0.0
  az(:) = 0.0
  bz(:) = 0.0

  do i=1,nd-1
     dtheta = theta(i+1)-theta(i)
     if (dtheta > pi) dtheta = dtheta-2*pi
     if (dtheta < -pi) dtheta = dtheta+2*pi
      do n=0,ns
        y = 0.5*(cos(n*theta(i+1))*rd(i+1)+cos(n*theta(i))*rd(i))
        ar(n) = ar(n)+dtheta*y/pi
        y = 0.5*(sin(n*theta(i+1))*rd(i+1)+sin(n*theta(i))*rd(i))
        br(n) = br(n)+dtheta*y/pi
        y = 0.5*(cos(n*theta(i+1))*zd(i+1)+cos(n*theta(i))*zd(i))
        az(n) = az(n)+dtheta*y/pi
        y = 0.5*(sin(n*theta(i+1))*zd(i+1)+sin(n*theta(i))*zd(i))
        bz(n) = bz(n)+dtheta*y/pi
     enddo
  enddo

end subroutine fluxfit_fourier
