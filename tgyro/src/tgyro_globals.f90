!-----------------------------------------------------
! tgyro_globals.f90
!
! PURPOSE:
!  Fundamental module containing most of the shared 
!  variables and arrays for tgyro
!
! NOTES:
!  Variables are split into three major categories:
!  1. Shared 
!  2. Specific to LOCAL method
!  3. Specific to GLOBAL method 
!-----------------------------------------------------

module tgyro_globals

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
  integer :: gyro_comm
  integer :: gyro_comm_rank
  integer :: gyro_adj
  integer :: gyro_adj_rank
  integer :: gyro_rad
  integer :: gyro_rad_rank
  !
  integer :: ierr
  !
  ! Containers for simulation path information
  !
  integer :: n_inst
  integer :: n_worker
  !
  character(len=80), dimension(:), allocatable :: paths
  integer, dimension(:), allocatable :: procs
  real, dimension(:), allocatable :: inputrads
  character(len=5), dimension(:), allocatable :: code
  !
  character(80) :: lpath, linput
  integer :: lproc
  character(len=5) :: lcode
  !
  ! Control variables
  !
  integer :: transport_method
  integer :: gyrotest_flag
  integer :: gyro_restart_method
  !
  integer :: error_flag
  character(len=80) :: error_msg
  !===============================================================

  !===============================================================
  ! Variables for LOCAL transport mode
  !
  integer :: shot
  integer, parameter :: n_ion_max = 9
  integer, parameter :: n_evolve_max = 5
  character(len=14) :: runfile='out.tgyro.run'
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
  real, dimension(:), allocatable   :: pflux_e_neo
  real, dimension(:,:), allocatable :: pflux_i_neo
  real, dimension(:), allocatable   :: pflux_e_tur
  real, dimension(:,:), allocatable :: pflux_i_tur
  real, dimension(:), allocatable   :: pflux_e_tot
  real, dimension(:,:), allocatable :: pflux_i_tot
  real, dimension(:), allocatable   :: pflux_e_target
  real, dimension(:), allocatable   :: pflux_he_target

  ! - momentum fluxes
  real, dimension(:), allocatable   :: mflux_e_neo
  real, dimension(:,:), allocatable :: mflux_i_neo
  real, dimension(:), allocatable   :: mflux_e_tur
  real, dimension(:,:), allocatable :: mflux_i_tur
  real, dimension(:), allocatable   :: mflux_tot
  real, dimension(:), allocatable   :: mflux_target

  ! - energy fluxes
  real, dimension(:), allocatable   :: eflux_e_neo
  real, dimension(:,:), allocatable :: eflux_i_neo
  real, dimension(:), allocatable   :: eflux_e_tur
  real, dimension(:,:), allocatable :: eflux_i_tur
  real, dimension(:), allocatable   :: eflux_e_tot
  real, dimension(:), allocatable   :: eflux_i_tot
  real, dimension(:), allocatable   :: eflux_e_target
  real, dimension(:), allocatable   :: eflux_i_target

  ! - exchange power densities
  real, dimension(:), allocatable :: expwd_e_tur
  real, dimension(:,:), allocatable :: expwd_i_tur

  ! Moving gyroBohm diffusivity
  real, dimension(:), allocatable :: chi_gb

  ! Moving gyroBohm fluxes
  real, dimension(:), allocatable :: gamma_gb
  real, dimension(:), allocatable :: pi_gb
  real, dimension(:), allocatable :: q_gb
  real, dimension(:), allocatable :: s_gb

  ! Collision frequencies
  real, dimension(:), allocatable :: nue
  real, dimension(:,:), allocatable :: nui
  real, dimension(:), allocatable :: nui_HH
  real, dimension(:), allocatable :: nue_HH
  real, dimension(:), allocatable :: nue_star

  real, dimension(:), allocatable :: z_eff

  ! Formulary exchange rate
  real, dimension(:), allocatable :: nu_exch

  ! Alpha heating coefficients
  real, dimension(:), allocatable :: frac_ae
  real, dimension(:), allocatable :: frac_ai
  real, dimension(:), allocatable :: e_cross

  ! Electron and ion temperatures
  real, dimension(:), allocatable :: te
  real, dimension(:), allocatable :: dlntedr
  real, dimension(:,:), allocatable :: ti
  real, dimension(:,:), allocatable :: dlntidr

  ! Electron and ion densities
  real, dimension(:), allocatable :: ne
  real, dimension(:), allocatable :: dlnnedr
  real, dimension(:,:), allocatable :: ni
  real, dimension(:,:), allocatable :: dlnnidr

  ! Rotation parameters
  real, dimension(:), allocatable :: w0
  real, dimension(:), allocatable :: w0p
  real, dimension(:), allocatable :: gamma_eb
  real, dimension(:), allocatable :: gamma_p
  real, dimension(:), allocatable :: u00
  real :: w0_norm

  real, dimension(:), allocatable :: pr
  real, dimension(:), allocatable :: ptot
  real, dimension(:), allocatable :: pext
  real, dimension(:), allocatable :: dpext
  real, dimension(:), allocatable :: dlnpdr
  real, dimension(:), allocatable :: dlnptotdr
  real, dimension(:), allocatable :: beta_unit
  real, dimension(:), allocatable :: betae_unit
  real, dimension(:), allocatable :: c_s
  real, dimension(:), allocatable :: v_i
  real, dimension(:), allocatable :: rho_s
  real, dimension(:), allocatable :: rho_i

  real, dimension(:), allocatable :: rho
  real, dimension(:), allocatable :: r
  real, dimension(:), allocatable :: r_maj
  real, dimension(:), allocatable :: q
  real, dimension(:), allocatable :: s
  real, dimension(:), allocatable :: kappa
  real, dimension(:), allocatable :: delta
  real, dimension(:), allocatable :: s_kappa
  real, dimension(:), allocatable :: s_delta
  real, dimension(:), allocatable :: shift
  real, dimension(:), allocatable :: zmag
  real, dimension(:), allocatable :: dzmag
  real, dimension(:), allocatable :: zeta
  real, dimension(:), allocatable :: s_zeta

  real, dimension(:), allocatable :: b_unit
  real, dimension(:), allocatable :: volp
  real, dimension(:), allocatable :: vol
  real, dimension(:), allocatable :: ave_grad_r
  real, dimension(:), allocatable :: er
  real, dimension(:), allocatable :: f_rot

  real, dimension(:), allocatable :: p_i_in
  real, dimension(:), allocatable :: p_e_in
  real, dimension(:), allocatable :: p_exch_in
  real, dimension(:), allocatable :: p_i
  real, dimension(:), allocatable :: p_e
  real, dimension(:), allocatable :: p_e_fus
  real, dimension(:), allocatable :: p_i_fus
  real, dimension(:), allocatable :: p_e_fus_in
  real, dimension(:), allocatable :: p_i_fus_in
  real, dimension(:), allocatable :: p_e_aux_in
  real, dimension(:), allocatable :: p_i_aux_in
  real, dimension(:), allocatable :: p_brem
  real, dimension(:), allocatable :: p_sync
  real, dimension(:), allocatable :: p_line
  real, dimension(:), allocatable :: p_exch
  real, dimension(:), allocatable :: p_brem_in
  real, dimension(:), allocatable :: p_sync_in
  real, dimension(:), allocatable :: p_line_in
  real, dimension(:), allocatable :: p_expwd
  real, dimension(:), allocatable :: s_alpha_i
  real, dimension(:), allocatable :: s_alpha_e
  real, dimension(:), allocatable :: sn_alpha
  real, dimension(:), allocatable :: s_brem
  real, dimension(:), allocatable :: s_sync
  real, dimension(:), allocatable :: s_exch
  real, dimension(:), allocatable :: s_expwd
  real, dimension(:), allocatable :: f_b_in
  real, dimension(:), allocatable :: f_w_in
  real, dimension(:), allocatable :: f_he_fus
  real, dimension(:), allocatable :: mf_in

  real, dimension(:), allocatable :: res
  real, dimension(:), allocatable :: res0
  real, dimension(:), allocatable :: relax

  integer, dimension(:), allocatable :: therm_vec

  real, dimension(n_ion_max) :: mi
  real, dimension(n_ion_max) :: n_ratio,t_ratio

  ! Physical constants
  real :: pi
  real :: e_alpha
  real :: e
  real :: k
  real :: me
  real :: mp
  real :: malpha
  real :: c
  real :: aspect_rat
  real :: mu_0
  !
  real :: b_ref
  real :: r_min
  !
  integer, dimension(:,:), allocatable :: pmap
  character(len=1), dimension(:), allocatable :: b_flag
  !
  ! Geometry
  !
  integer :: n_fourier_geo
  real, dimension(:,:,:), allocatable :: a_fourier_geo
  !
  ! Orientation
  !
  integer :: signb
  integer :: signq
  !
  ! Control variables
  !
  integer :: tgyro_mode
  integer :: tgyro_relax_iterations
  real :: loc_nu_scale
  real :: loc_dx
  real :: loc_dx_gyro
  real :: loc_dx_max
  real :: loc_relax
  integer :: loc_lock_profile_flag
  integer :: loc_evolve_grad_only_flag
  integer :: loc_restart_flag
  integer :: loc_scenario
  integer :: loc_neo_method
  integer :: loc_n_ion
  real, dimension(n_ion_max) :: zi_vec
  real, dimension(n_ion_max) :: mi_vec
  integer, dimension(n_ion_max) :: therm_flag
  integer, dimension(0:n_ion_max) :: evo_e
  real, dimension(0:n_ion_max) :: evo_c
  real :: loc_betae_scale
  integer :: loc_chang_hinton
  real :: loc_me_multiplier
  integer :: loc_bc_offset
  integer :: tgyro_tglf_revision
  integer :: tgyro_tglf_dump_flag
  integer :: tgyro_glf23_revision
  integer :: tgyro_glf23_dump_flag
  integer :: loc_ti_feedback_flag
  integer :: loc_te_feedback_flag
  integer :: loc_er_feedback_flag
  integer :: loc_he_feedback_flag
  integer :: loc_zeff_flag
  integer :: loc_pflux_method
  integer :: loc_residual_method
  integer :: loc_num_equil_flag
  integer :: tgyro_neo_gv_flag
  integer :: tglf_q_low_flag
  integer :: tgyro_iteration_method
  integer :: tgyro_rotation_flag
  real :: tgyro_rmin
  real :: tgyro_rmax
  integer :: tgyro_expwd_flag
  real :: tgyro_input_den_scale
  real :: tgyro_input_te_scale
  real :: tgyro_input_ti_scale
  real :: tgyro_input_w0_scale
  real :: tgyro_input_paux_scale
  real :: tgyro_input_dlntdr_scale
  integer :: tgyro_er_bc
  integer :: tgyro_noturb_flag
  integer :: tgyro_use_rho
  integer :: tgyro_dt_method
  integer :: tgyro_gyro_restart_flag
  integer :: tgyro_write_profiles_flag
  integer :: tgyro_neo_n_theta
  integer :: tgyro_ptot_flag
  integer :: tgyro_ped_model
  real :: tgyro_rped
  real :: tgyro_neped
  real :: tgyro_zeffped
  real :: tgyro_ped_ratio
  real :: tgyro_tglf_nn_max_error
  integer :: tgyro_quickfast_flag
  !
  ! Iteration variables (global)
  !
  integer :: n_evolve
  integer :: p_max
  integer :: i_r
  integer :: n_r
  integer :: flux_method
  integer, dimension(:), allocatable :: flux_method_vec
  integer :: i_tran
  integer :: i_bc
  integer :: flux_counter
  integer :: i_ash
  integer :: evolve_indx(5)
  !
  ! Global TGYRO variables
  integer :: igmin
  integer :: igmax
  real :: length
  real :: dlength
  !---------------------------------------------------------

end module tgyro_globals
