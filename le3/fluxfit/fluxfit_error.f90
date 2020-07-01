subroutine fluxfit_error(err)

  use fluxfit_globals

  implicit none

  ! Output
  real, intent(inout) :: err

  ! Number of points in search arc
  integer, parameter :: nt=9

  ! Iteration tolerance
  real, parameter :: tol=1e-9

  ! Internal variables
  integer :: i,p,k
  real, dimension(:), allocatable :: d
  real, dimension(:), allocatable :: tt
  real :: dmin,dmin0,q
  real :: r,z,t

  allocate(d(nt))
  allocate(tt(nt))

  err = 0.0
  do i=1,nd-1
  
     ! For given flux-surface datapoint, find minimum distance from 
     ! continuous model curve.

     ! 1. First, get an initial guess by looping around the contour

     if (i < nd/2) then
        do k=1,nt 
           tt(k) = -0.5*pi+(k-1)*2.0*pi/(nt-1)
           call fluxfit_f_model(tt(k),r,z)
           d(k) = sqrt((r-rd(i))**2+(z-zd(i))**2)
        enddo
     else
        do k=1,nt 
           tt(k) = 0.5*pi+(k-1)*2.0*pi/(nt-1)
           call fluxfit_f_model(tt(k),r,z)
           d(k) = sqrt((r-rd(i))**2+(z-zd(i))**2)
        enddo
     endif
   
     call fluxfit_minmax(d,tt,nt,dmin,t,"min")

     ! 2. Now refine the guess

     dmin0 = 0.0
     p = 0

     do while (abs(dmin-dmin0) > tol)
        dmin0 = dmin
        p = p+1
        do k=1,nt
           q = (nt*1.0)**p
           tt(k) = t-pi/q+(k-1.0)*pi/(nt-1)/q
           call fluxfit_f_model(tt(k),r,z)
           d(k) = sqrt((r-rd(i))**2+(z-zd(i))**2)
        enddo
        k = minloc(d,1)
        t    = tt(k)
        dmin = d(k)
     enddo
     err = err+dmin/(nd-1)/rmin
  enddo

  deallocate(d)
  deallocate(tt)

end subroutine fluxfit_error
