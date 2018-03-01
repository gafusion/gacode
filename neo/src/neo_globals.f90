module neo_globals

  !---------------------------------------------------------------
  ! local MPI and OpenMP variables
  ! 
  integer :: i_proc
  integer :: n_proc
  integer :: NEO_COMM_WORLD
  integer :: n_omp
  !---------------------------------------------------------------

  real, parameter :: pi=3.1415926535897932

  !---------------------------------------------------------------
  ! Input parameters:
  !
  real :: rmin_1_in
  real :: rmin_2_in
  real :: rmaj_in
  real :: dphi0dr_in
  real :: epar0_in
  real :: omega_rot_in
  real :: omega_rot_deriv_in
  real :: q_in
  real :: rho_in
  real :: shear_in
  real :: shift_in
  real :: kappa_in
  real :: s_kappa_in
  real :: delta_in
  real :: s_delta_in
  real :: zeta_in
  real :: s_zeta_in
  real :: zmag_in
  real :: s_zmag_in
  real :: beta_star_in
  integer :: geo_ny_in
  real, dimension(8,0:32) :: geo_yin_in
  !
  real :: profile_delta_scale
  real :: profile_zeta_scale
  real :: profile_zmag_scale

  real, dimension(11) :: profile_dlnndr_scale
  real, dimension(11) :: profile_dlntdr_scale
  ! 
  integer :: ipccw_in
  integer :: btccw_in
  !
  real :: te_ade_in
  real :: ne_ade_in
  real :: dlntdre_ade_in
  real :: dlnndre_ade_in
  !
  real, dimension(11) :: z_in 
  real, dimension(11) :: mass_in
  real, dimension(11) :: dens_in
  real, dimension(11) :: temp_in
  real, dimension(11) :: dlnndr_in
  real, dimension(11) :: dlntdr_in
  !
  integer, dimension(11) :: aniso_model_in
  real, dimension(11)    :: temp_para_in
  real, dimension(11)    :: dlntdr_para_in
  real, dimension(11)    :: temp_perp_in
  real, dimension(11)    :: dlntdr_perp_in
  !
  real :: nu_1_in
  character(len=16) :: rbf_dir
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Grid dimensions:
  !
  integer :: n_species
  integer :: n_energy
  integer :: n_xi
  integer :: n_theta
  integer :: n_radial
  integer :: is_ele
  integer :: matsz_scalefac
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Models:
  !
  integer :: sim_model
  integer :: equilibrium_model
  integer :: collision_model 
  integer :: profile_model
  integer :: profile_erad0_model
  integer :: profile_equilibrium_model
  integer :: rotation_model
  integer :: adiabatic_ele_model
  integer :: spitzer_model
  real    :: epar0_spitzer
  integer :: coll_uncoupledei_model
  integer :: coll_uncoupledaniso_model
  real    :: sign_q
  real    :: sign_bunit
  integer :: subroutine_flag  ! only used for neo_read_input
  integer :: laguerre_method
  integer :: write_cmoments_flag
  integer :: use_slu ! use the SuperLU solver
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! 3D parameters
  !
  integer :: threed_model
  integer :: threed_exb_model
  real    :: threed_exb_dphi0dr
  integer :: threed_drift_model
  real    :: threed_hyperxi
  !
  ! 3D grid dimensions (out.le3.geoscalar)
  integer :: tpmatsize
  integer :: n_tptheta
  integer :: n_tpvarphi
  integer :: indx_c00
  !
  ! LAPACK VARIABLES
  integer :: bw
  integer :: npb
  integer :: ldab
  integer :: nb
  integer :: n_row
  real, dimension(:,:), allocatable :: ab
  integer, dimension(:), allocatable :: ipiv
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Output mode:
  !
  integer :: silent_flag
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Charge and mass
  !
  real, dimension(:), allocatable :: z        ! charge (ns)
  real, dimension(:), allocatable :: mass     ! m/m_0  (ns)
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Path to INPUT, read in the get_inputpath subroutine
  character(len=80) :: path
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Profile functions
  !
  real, dimension(:), allocatable :: r
  real, dimension(:), allocatable :: rmaj
  real, dimension(:), allocatable :: q
  real, dimension(:), allocatable :: rho
  real, dimension(:), allocatable :: shear
  real, dimension(:), allocatable :: shift
  real, dimension(:), allocatable :: kappa
  real, dimension(:), allocatable :: s_kappa
  real, dimension(:), allocatable :: delta
  real, dimension(:), allocatable :: s_delta
  real, dimension(:), allocatable :: zeta
  real, dimension(:), allocatable :: s_zeta
  real, dimension(:), allocatable :: zmag
  real, dimension(:), allocatable :: s_zmag
  real, dimension(:), allocatable :: beta_star
  integer :: geo_numeq_flag
  integer :: geo_ny
  real, dimension(:,:,:), allocatable :: geo_yin
  !
  real, dimension(:), allocatable :: te_ade, ne_ade, dlnndre_ade, dlntdre_ade
  ! ele temp and dens -- used only if ade (for Poisson calculation)
  !
  ! (species-dependent)
  !
  real, dimension(:,:), allocatable :: dens   ! n/n_0  (ns,nr)
  real, dimension(:,:), allocatable :: temp   ! T/T_0  (ns,nr)
  real, dimension(:,:), allocatable :: vth    ! vth/vth_0=sqrt(temp/mass) 
  real, dimension(:,:), allocatable :: dlnndr ! a/Ln   (ns,nr)
  real, dimension(:,:), allocatable :: dlntdr ! a/LT   (ns,nr)
  real, dimension(:,:), allocatable :: nu     ! nu / (vth_0/a) (ns,nr)
  !
  real, dimension(:), allocatable :: dphi0dr  ! dphi0/dr (e/T_0)*a (nr) 
  ! erad0 = -dphi0/dr * |grad r(theta=0)|
  real, dimension(:), allocatable :: epar0    ! inductive Efield (e/T_0)*a (nr)
  real, dimension(:), allocatable :: omega_rot  ! omega_rot / (vth_0/a) (nr)
  real, dimension(:), allocatable :: omega_rot_deriv ! domega_rot/d(r/a) (nr)
  
  integer, dimension(:), allocatable :: aniso_model    ! 1=isotropic,2=aniso
  real, dimension(:,:), allocatable :: temp_para       ! Tpara/T_0  (ns,nr)
  real, dimension(:,:), allocatable :: dlntdr_para     ! a/LTpara   (ns,nr)
  real, dimension(:,:), allocatable :: temp_perp       ! Tperp/T_0  (ns,nr)
  real, dimension(:,:), allocatable :: dlntdr_perp     ! a/LTperp   (ns,nr)
  real, dimension(:,:), allocatable :: vth_para        ! vthpara/vth_0=sqrt(temp_para/mass)

  !---------------------------------------------------------------

  real, dimension(:), allocatable :: theta

  ! Numerical/work arrays and dimensions

  integer :: n_max

  real, dimension(:), allocatable :: g
  real, dimension(:), allocatable :: g_loc
  real, dimension(:), allocatable :: amat
  integer, dimension(:), allocatable :: amat_indx

  integer, dimension(:), allocatable :: is_indx, ie_indx, ix_indx, it_indx

  real, dimension(:,:), allocatable :: driftx, &
       driftxrot1, driftxrot2, driftxrot3
  ! driftx   = vdrift dot grad r

  ! normalizations for experimental profiles
  real  :: a_meters
  real, dimension(:), allocatable :: dens_norm
  real, dimension(:), allocatable :: temp_norm
  real, dimension(:), allocatable :: vth_norm
  real, dimension(:), allocatable :: b_unit

  real, dimension(:), allocatable :: rhoN_torflux
  real, dimension(:), allocatable :: psiN_polflux
  real :: psiN_polflux_a

  ! theta derivative variables
  integer, dimension(-2:2) :: cderiv
  integer, dimension(:), allocatable :: thcyc

  ! error checking
  integer :: error_status = 0
  character(len=80) :: error_message = '(NEO) completed successfully.'

  ! output file
  character(len=80)  :: runfile_neoout = 'out.neo.run'
  integer :: io_neoout = 12

  ! output vectors

  ! (n_species_max, transport coeff)
  ! transport coeff: 1-> gamma, 2-> Q, 3->Pi, 4-> Q-omega*Pi
  !                  5-> vpol,  6-> vtor
  real, dimension(11,6) :: neo_dke_out=0.0

  ! species-independent transport coeff (jpar, jtor)
  real, dimension(2)    :: neo_dke_1d_out= 0.0

  ! (n_species_max, transport coeff)  
  ! gyro-viscous fluxes: 1-> gamma, 2-> Q, 3->Pi, 4-> Q-omega*Pi
  real, dimension(11,4) :: neo_gv_out=0.0

  ! pure plasma theory
  ! Gamma_HH, Qi_HH, Qe_HH, Qi_CH, jpar_S, jpar_K, jpar_N, jtor_S
  real, dimension(8)   :: neo_th_out=0.0

  ! nclass (n_species_max, transport coeff)
  ! transport coeff: 1-> gamma, 2-> Q, 3-> vpol,  4-> vtor, 5-> visc coeff
  real, dimension(11,5) :: neo_nclass_out=0.0

  ! species-independent transport coeff (currently just jpar)
  real                 :: neo_nclass_1d_out= 0.0

  ! (n_species_max, transport coeff)
  ! multi-species theory: 1-> gamma_HS, 2-> Q_HS
  real, dimension(11,2)   :: neo_thHS_out=0.0

  ! geometry parameters
  real, dimension(5) :: neo_geo_out=0.0
  
end module neo_globals
