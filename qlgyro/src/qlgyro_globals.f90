!-----------------------------------------------------
! qlgyro_globals.f90
!
! PURPOSE:
!  Fundamental module containing most of the shared 
!  variables and arrays for qlgyro
!
! NOTES:
!  Variables are split into three major categories:
!  1. Shared 
!  2. Specific to LOCAL method
!  3. Specific to GLOBAL method 
!-----------------------------------------------------

module qlgyro_globals

  !===============================================================
  ! SHARED parameters for both transport methods:
  !
  ! MPI related integers
  !
  integer :: i_proc_global
  integer :: n_proc_global
  integer :: block
  integer :: color
  integer :: worker
  integer :: adjoint
  integer :: workeradj
  integer :: worker_index
  !
  integer :: qlgyro_comm
  integer :: qlgyro_comm_rank
  integer :: qlgyro_adj
  integer :: qlgyro_adj_rank
  integer :: qlgyro_rad
  integer :: qlgyro_rad_rank
  integer :: GYRO_COMM_WORLD
  integer :: CGYRO_COMM_WORLD
  integer :: QLGYRO_COMM_WORLD
  !
  integer :: ierr
  !
  ! Containers for simulation path information
  !
  integer :: n_inst
  integer :: n_worker
  !
  character(len=80) :: path
  integer :: procs
  real, dimension(:), allocatable :: inputrads

  ! Code to run, 0 for GYRO, 1 for CGYRO
  integer :: code
  !
  character(80) :: lpath, linput
  integer :: lproc
  character(len=5) :: lcode
  !
  ! Control variables
  !
  integer :: transport_method
  integer :: tag_removal = 0
  integer :: gyrotest_flag
  integer :: gyro_restart_method
  !
  integer :: error_flag
  character(len=80) :: error_msg
  !
  ! Ky run status string
  integer :: statstr
  !
  ! Automatically set box size
  integer :: auto_box_size = 0
  real :: kx_max_box = 1e10
  
  !===============================================================

  !===============================================================
  ! Variables for LOCAL transport mode
  !
  integer :: shot
  integer, parameter :: n_ion_max = 9
  integer, parameter :: n_evolve_max = 5
  character(len=14) :: runfile='out.qlgyro.run'
  character(len=1), dimension(n_ion_max) :: &
       ion_tag=(/'1','2','3','4','5','6','7','8','9'/)
  !
  ! Component fluxes
  !
  ! neo    -> neoclassical
  ! tur    -> turbulent
  ! tot    -> total
  ! target -> target (source) 

  ! - particle fluxes
  real  :: pflux_e_tur
  real  :: pflux_i_tur
  
  ! - momentum fluxes
  real  :: mflux_e_tur
  real  :: mflux_i_tur

  ! - energy fluxes
  real  :: eflux_e_tur
  real  :: eflux_i_tur

  ! - exchange power densities
  real :: expwd_e_tur
  real :: expwd_i_tur

  
  ! Moving gyroBohm diffusivity
  real :: chi_gb

  ! Moving gyroBohm fluxes
  real :: gamma_gb
  real :: pi_gb
  real :: q_gb
  real :: s_gb

  ! Collision frequencies
  real :: nue
  real, dimension(:), allocatable :: nui
  real :: nui_HH
  real :: nue_HH
  real :: nue_star

  real :: z_eff

  ! Formulary exchange rate
  real :: nu_exch

  ! Alpha heating coefficients
  real :: frac_ae
  real :: frac_ai
  real :: e_cross

  ! Electron and ion temperatures
  real :: te
  real :: dlntedr
  real, dimension(:), allocatable :: ti
  real, dimension(:), allocatable :: dlntidr

  ! Electron and ion densities
  real :: ne
  real :: dlnnedr
  real, dimension(:), allocatable :: ni
  real, dimension(:), allocatable :: dlnnidr

  real, dimension(n_ion_max) :: mi
  real, dimension(n_ion_max) :: zi_vec
  real, dimension(n_ion_max) :: mi_vec
  integer, dimension(n_ion_max) :: therm_flag

  
  ! Rotation parameters
  real :: w0
  Real :: w0p
  real :: gamma_e
  real :: gamma_p
  real :: u00
  real :: w0_norm
  integer :: qlgyro_rotation_flag=1

  real :: pr
  real :: ptot
  real :: pext
  real :: dpext
  real :: dlnpdr
  real :: dlnptotdr
  real :: beta_unit
  real :: betae_unit
  real :: c_s
  real :: v_i
  real :: rho_s
  real :: rho_i

  real :: a_min
  real :: r_min
  real :: r_maj
  real :: q
  real :: s
  real :: kappa
  real :: delta
  real :: s_kappa
  real :: s_delta
  real :: shift
  real :: zmag
  real :: dzmag
  real :: zeta
  real :: s_zeta
  real :: shape_sin3
  real :: shape_ssin3
  real :: shape_cos0
  real :: shape_scos0
  real :: shape_cos1
  real :: shape_scos1
  real :: shape_cos2
  real :: shape_scos2
  real :: shape_cos3
  real :: shape_scos3

  real :: bunit
  real :: bgs2
  real :: volp
  real :: vol
  real :: ave_grad_r
  real :: er
  real :: f_rot


  ! Physical constants
  real :: pi
  real :: e_alpha
  real :: e
  real :: bk
  real :: me
  real :: mp
  real :: malpha
  real :: c
  real :: aspect_rat
  real :: mu_0
  
  real :: gamma_exb
  real :: gamma_p0

  !---------------------------------------------------------
  ! Equilibrium normalisations
  !
  real :: a_meters
  !

  real :: rho_star_norm
  real :: csda_norm
  real :: temp_norm
  real :: dens_norm
  real :: lambda_star
  ! Proton mass in kg:
  !
  real, parameter :: kg_proton = 1.6726e-27

  !
  real :: b_ref
  !
  ! Geometry
  !
  integer :: n_fourier_geo
  real, dimension(:,:), allocatable :: a_fourier_geo
  !
  ! Orientation
  !
  integer :: signb
  integer :: signq
  !
  integer :: qlgyro_num_equil_flag=1

  integer :: nspec_max=6

  ! Grad r correction term
  real, dimension(:), allocatable :: sat_geo_spectrum
  real, dimension(:), allocatable :: kxrms_spectrum

  real, dimension(:, :, :), allocatable :: qlgyro_eigenvalue_spectrum_out
  real, dimension(:, :, :, :, :), allocatable :: qlgyro_flux_spectrum_out
  complex, dimension(:, :, :, :), allocatable :: qlgyro_field_spectrum_out
  real, dimension(:), allocatable :: qlgyro_theta_ballooning

  real, dimension(:, :, :), allocatable :: qlgyro_k_perp, qlgyro_jacobian
  
  integer :: sat_rule, restart_mode
  !
  !
  ! Iteration variables (global)
  !
  integer :: n_evolve
  integer :: p_max
  integer :: i_ion
  integer :: n_ion=1
  integer :: n_ky
  integer :: n_px0
  integer :: n_kypx0
  integer :: n_modes=1
  integer :: n_field=1
  integer :: n_theta
  integer :: n_species
  integer :: n_thetab
  integer :: i_tor

  ! No. of parallel runs
  integer :: n_parallel=1

  ! No. of GYRO runs to attempt to converge
  integer :: n_runs=1

  ! Color. of cores running each ky value
  integer, dimension(:), allocatable :: ky_color
  real, dimension(:), allocatable :: px0_spectrum
  integer, dimension(:), allocatable :: i_ky, i_px0
  integer :: kygrid_model=1
  integer :: px0grid_model=1
  real :: ky_in
  integer :: xnu_model
  
  integer :: flux_method
  integer, dimension(:), allocatable :: flux_method_vec
  integer :: i_tran=1
  character(30) :: iter_path="./"
  character(30) :: worker_str
  integer :: i_bc
  integer :: flux_counter
  integer :: i_ash
  integer :: evolve_indx(5)

  integer :: elec_ind, ion_start, ion_end, loc_n_ion
  !
  integer :: use_trap
  !---------------------------------------------------------
  integer :: debug_flag=0
  
  ! runpath for each KY
  character(len=120) :: runpath
  
  ! TGLF specific parameters
  real :: qlgyro_tglf_nn_max_error=1.0e-3
  integer :: qlgyro_tglf_revision=2
  integer ::qlgyro_tglf_dump_flag=0
end module qlgyro_globals
