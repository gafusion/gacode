module locpargen_globals

  integer :: is,ise
  real :: r0
  real :: rho0
  real :: psi0
  real :: a
  real, dimension(1) :: x,y
  integer :: hasgeo
  integer :: qnflag
  integer :: appendflag
  
  character(len=1) :: tag(5)

  real :: btccw,ipccw
  real :: cc,loglam,nu_ee,pi,betae_unit

end module locpargen_globals
