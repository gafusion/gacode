module vgen_globals
  
  character(len=80) :: path
  
  integer :: i_proc
  integer :: n_proc
  integer :: i_err

  integer :: i_loc,n_loc
  integer, dimension(:), allocatable :: i_glob
  
  ! vgen inputs
  integer :: er_method
  integer :: erspecies_indx
  integer :: vel_method
  integer :: nth_min
  integer :: nth_max
  integer :: nn_flag
  
  real, parameter :: pi=3.1415926535897932
  real, parameter :: mass_deuterium = 3.3452   ! (x 1e-27 kg)
  real, parameter :: temp_norm_fac   = 1602.2
  real, parameter :: charge_norm_fac = 1.6022
  
  real :: dens_norm, temp_norm, mass_norm, vth_norm, jbs_norm

  real, dimension(11) :: zfac

  real, dimension(:), allocatable :: vtor_measured

  character(len=8) :: fmt='(I2.2)' 
  character(len=2), dimension(100) :: tag

  real, dimension(:), allocatable :: jbs_neo
  real, dimension(:), allocatable :: jbs_sauter
  real, dimension(:), allocatable :: jbs_koh
  real, dimension(:), allocatable :: jbs_nclass
  real, dimension(:), allocatable :: pflux_sum
  real, dimension(:), allocatable :: jtor_neo
  real, dimension(:), allocatable :: jtor_sauter
  
  integer :: n_ions

  integer, parameter :: timing_flag = 0

end module vgen_globals
