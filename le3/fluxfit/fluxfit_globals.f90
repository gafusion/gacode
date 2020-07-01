module fluxfit_globals

  integer :: model
  integer :: nd
  integer :: nc
  integer :: ns

  real, parameter :: pi=3.14159265358979323846

  real :: rmin,rmaj
  real :: r_c,z_c

  ! Decimated flux-surface data
  real, dimension(:), allocatable :: rd
  real, dimension(:), allocatable :: zd

  real, dimension(:), allocatable :: theta

  character (len=11), dimension(:), allocatable :: tag

  ! Model 1 coefficients
  real, dimension(:), allocatable :: c

  ! Model 2 (Fourier) coefficients
  real, dimension(:), allocatable :: ar,br
  real, dimension(:), allocatable :: az,bz
  
end module fluxfit_globals
