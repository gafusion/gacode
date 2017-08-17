module create_globals

  integer,parameter :: nions_max=5

  ! Input variables

  integer :: nx      ! used only if input.profiles created from scratch
  real :: exm_b_ref
  real :: exm_arho
  real :: exm_a_meters
  real :: exm_rmaj
  real :: exm_q
  real :: exm_kappa
  real :: exm_delta
  real :: exm_te_axis
  real :: exm_ne_axis
  real :: exm_mu_e
  real :: exm_z_eff
  real :: exm_w0
  real :: exm_zeta
  real :: exm_zmag
  real, dimension(nions_max) :: exm_ni_axis
  real, dimension(nions_max) :: exm_ti_axis
  integer, dimension(nions_max) :: exm_z
  real, dimension(nions_max) :: exm_mu
  integer :: exm_pivot
  
  ! Models
  integer :: exm_te_model
  real :: exm_alte
  integer :: exm_ne_model
  real :: exm_alne
  integer, dimension(nions_max) :: exm_ti_model
  real, dimension(nions_max) :: exm_alti
  integer, dimension(nions_max) :: exm_ni_model
  real, dimension(nions_max) :: exm_alni
  integer :: exm_z_eff_model

  ! Set flags
  integer :: set_exm_b_ref
  integer :: set_exm_arho
  integer :: set_exm_rho
  integer :: set_exm_rmin
  integer :: set_exm_rmaj
  integer :: set_exm_q
  integer :: set_exm_kappa
  integer :: set_exm_delta
  integer :: set_exm_te
  integer :: set_exm_ne
  integer :: set_exm_z_eff
  integer :: set_exm_w0
  integer :: set_exm_zeta
  integer :: set_exm_zmag
  integer, dimension(nions_max) :: set_exm_ni
  integer, dimension(nions_max) :: set_exm_ti
  

  ! parameters not presently controllable from input
  real :: exm_flow_mom  = 0.0
  real :: exm_pow_e     = 0.0
  real :: exm_pow_i     = 0.0
  real :: exm_pow_ei    = 0.0
  real :: exm_flow_beam = 0.0
  real :: exm_flow_wall = 0.0
  real :: exm_ptot      = 0.0
  real, dimension(nions_max) :: exm_vpol = 0.0
  real, dimension(nions_max) :: exm_vtor = 0.0

  integer :: set_exm_flow_mom  
  integer :: set_exm_pow_e     
  integer :: set_exm_pow_i    
  integer :: set_exm_pow_ei    
  integer :: set_exm_flow_beam 
  integer :: set_exm_flow_wall 
  integer :: set_exm_ptot     
  integer :: set_exm_polflux   
  integer, dimension(nions_max) :: set_exm_vpol    
  integer, dimension(nions_max) :: set_exm_vtor    

  character(len=80),dimension(nions_max) :: exm_ni_data
  character(len=80) :: exm_ne_data

end module create_globals
