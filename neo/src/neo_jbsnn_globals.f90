module jbsnn_globals
  
  character(len=80) :: path
  
  integer :: i_proc
  integer :: n_proc
  integer :: i_err
  integer :: i_loc,n_loc
  integer, dimension(:), allocatable :: i_glob
  
  ! jbsnn inputs 
  integer :: plasma_spec
  integer :: impurities_spec
  integer :: nth_min
  integer :: nth_max
  
  ! Needed constants
  real, parameter :: pi=3.1415926535897932
  real, parameter :: mass_deuterium = 3.3452   ! (x 1e-27 kg)
  real, parameter :: temp_norm_fac   = 1602.2
  real, parameter :: charge_norm_fac = 1.6022
  real :: dens_norm, temp_norm, mass_norm, vth_norm, jbs_norm

  character(len=8) :: fmt='(I2.2)' 
  character(len=2), dimension(100) :: tag
  real, dimension(:), allocatable :: jbsnn_neo
  real, dimension(:), allocatable :: jbsnn_sauter
  integer :: n_ions
  integer, parameter :: timing_flag = 0

end module jbsnn_globals
