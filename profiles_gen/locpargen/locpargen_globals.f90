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
  
  character(len=2) :: tag(9)

  real :: btccw,ipccw
  real :: pi
  real :: cc,loglam,nu_ee
  real :: betae_unit,lambda_star

end module locpargen_globals
