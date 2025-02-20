!-----------------------------------------------------------------
! cgyro_globals.f90
!
! PURPOSE:
!  CGYRO global variables.  The idea is to have a primary, large
!  module containing all essential CGYRO arrays and scalars.
!-----------------------------------------------------------------

module cgyro_globals

  use, intrinsic :: iso_c_binding
#if defined(_OPENACC) || defined(OMPGPU)
#define CGYRO_GPU_FFT
#endif

  use, intrinsic :: iso_fortran_env

  integer :: sbeta_const_flag
  real :: sbeta_h
  
  !---------------------------------------------------------------
  ! Input parameters:
  !
  integer :: n_energy
  integer :: n_xi
  integer :: n_theta
  integer :: n_radial
  integer :: n_toroidal
  integer :: n_field
  real    :: e_max
  real    :: alpha_poly
  integer :: e_fix
  integer :: delta_t_method
  real    :: delta_t
  real    :: error_tol
  real    :: max_time
  integer :: print_step
  integer :: restart_step
  integer :: restart_preservation_mode
  real    :: freq_tol
  real    :: up_radial
  real    :: up_theta
  real    :: up_alpha
  integer :: nup_radial
  integer :: nup_theta
  integer :: nup_alpha
  integer :: n_wave
  integer :: constant_stream_flag
  integer :: explicit_trap_flag
  real    :: ky
  integer :: box_size
  real    :: ipccw
  real    :: btccw
  integer :: silent_flag
  integer :: profile_model
  integer :: equilibrium_model
  integer :: collision_model
  integer :: collision_mom_restore
  integer :: collision_ene_restore
  integer :: collision_ene_diffusion
  integer :: collision_kperp
  integer :: collision_field_model
  integer :: collision_ion_model
  integer :: collision_precision_mode
  integer :: collision_test_mode
  integer :: collision_field_max_l
  integer :: collision_test_max_l
  real    :: z_eff
  integer :: z_eff_method
  integer :: zf_test_mode 
  integer :: nonlinear_flag 
  real :: temp_ae
  real :: dens_ae
  real :: mass_ae
  real :: dlntdr_ae   
  real :: dlnndr_ae   
  real :: lambda_star
  integer :: test_flag
  integer :: h_print_flag
  integer :: moment_print_flag
  integer :: gflux_print_flag
  integer :: field_print_flag
  integer :: triad_print_flag
  real :: amp0
  real :: amp
  real :: gamma_e
  real :: gamma_p
  real :: mach
  integer :: rotation_model
  integer :: mpi_rank_order
  integer :: velocity_order
  integer :: hiprec_flag
  integer :: udsymmetry_flag
  integer :: shear_method
  integer :: global_flag
  integer :: n_global
  real    :: nu_global
  integer :: psym_flag
  integer :: theta_plot
  integer :: gpu_bigmem_flag
  integer :: upwind_single_flag
  real :: px0
  integer :: stream_term
  real :: stream_factor
  integer :: exch_flag
  real :: res_weight_power
  !
  ! Geometry input
  !
  real :: rmin
  real :: rmaj
  real :: q
  real :: s
  real :: shift
  real :: kappa
  real :: s_kappa
  real :: delta
  real :: s_delta
  real :: zeta
  real :: s_zeta
  real :: zmag
  real :: dzmag
  integer, parameter :: n_shape=6
  real, dimension(0:n_shape) :: shape_sin
  real, dimension(0:n_shape) :: shape_s_sin
  real, dimension(0:n_shape) :: shape_cos
  real, dimension(0:n_shape) :: shape_s_cos
  real :: betae_unit
  !
  ! Species parameters
  !
  integer :: n_species
  real :: nu_ee
  real, dimension(11) :: z
  real, dimension(11) :: mass
  real, dimension(11) :: dens
  real, dimension(11) :: temp
  real, dimension(11) :: dlnndr
  real, dimension(11) :: dlntdr
  real, dimension(11) :: sdlnndr
  real, dimension(11) :: sdlntdr
  real :: sbeta

  integer :: subroutine_flag  ! only used for cgyro_read_input

  ! Re-scaling parameters for experimental profiles
  integer :: quasineutral_flag
  real :: lambda_star_scale
  real :: gamma_e_scale
  real :: gamma_p_scale
  real :: mach_scale
  real :: beta_star_scale, betae_unit_scale
  real :: nu_ee_scale
  real :: zf_scale
  real, dimension(11) :: dlnndr_scale
  real, dimension(11) :: dlntdr_scale

  real :: lambda_debye
  real :: rhos

  ! Global conversion variables
  real :: b_unit, b_gs2
  real :: a_meters
  real :: dens_norm, temp_norm, vth_norm, mass_norm, rho_star_norm
  real :: gamma_gb_norm, q_gb_norm, pi_gb_norm
  
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! MPI/OpenMP variables and pointers
  ! 
  integer :: n_omp
  !
  integer :: i_err
  integer :: i_proc
  integer :: i_proc_1
  integer :: i_proc_2
  integer :: i_proc_3
  integer :: i_proc_restart_io
  integer :: n_proc
  integer :: n_proc_1
  integer :: n_proc_2
  integer :: n_proc_restart_io
  integer :: i_group_1
  integer :: i_group_2
  integer :: i_group_3
  integer :: i_group_restart_io
  integer :: CGYRO_COMM_WORLD
  integer :: NEW_COMM_1
  integer :: NEW_COMM_2
  integer :: NEW_COMM_3
  integer :: nv1,nv2,nc1,nc2
  integer :: nsplit,nsplitA,nsplitB
  integer :: ns1,ns2
  integer, dimension(:), allocatable :: recv_status
  integer :: fA_req, fB_req, eA_req, eB_req, g_req
  logical :: fA_req_valid, fB_req_valid, g_req_valid
  ! Thetas present in the process after NL AllToAll
  integer :: n_jtheta
  !
  ! Pointers
  integer :: nv,iv
  integer :: nv_loc,iv_loc
  integer :: nc,ic
  integer :: nc_loc,ic_loc
  integer :: ns_loc
  integer, dimension(:), allocatable :: ie_v
  integer, dimension(:), allocatable :: ix_v
  integer, dimension(:), allocatable :: is_v
  integer, dimension(:), allocatable :: ir_c
  integer, dimension(:), allocatable :: it_c
  integer, dimension(:,:), allocatable :: ic_c
  integer, dimension(:,:,:), allocatable :: iv_v
  integer, dimension(:), allocatable :: ica_c,icb_c
  !
  integer :: nt1,nt2,nt_loc
  integer :: n_toroidal_procs
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Constants
  !
  real, parameter    :: pi   = 3.1415926535897932
  complex, parameter :: i_c  = (0.0,1.0)
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! I/O and error management variables
  !
  character(len=80) :: path
  character(len=14) :: runfile_info    = 'out.cgyro.info'
  character(len=13) :: runfile_mpi     = 'out.cgyro.mpi'
  character(len=16) :: runfile_memory  = 'out.cgyro.memory'
  character(len=15) :: runfile_hosts   = 'out.cgyro.hosts'
  character(len=17) :: runfile_restart = 'bin.cgyro.restart'
  character(len=13) :: runfile_restart_tag = 'out.cgyro.tag'
  character(len=15) :: runfile_grids   = 'out.cgyro.grids'
  character(len=14) :: runfile_prec    = 'out.cgyro.prec'
  character(len=14) :: runfile_time    = 'out.cgyro.time'
  character(len=16) :: runfile_timers  = 'out.cgyro.timing'
  character(len=18) :: runfile_startups= 'out.cgyro.startups'
  character(len=14) :: runfile_freq    = 'out.cgyro.freq'
  character(len=14) :: binfile_freq    = 'bin.cgyro.freq'
  character(len=12) :: binfile_hb      = 'bin.cgyro.hb'
  character(len=17) :: binfile_ky_flux = 'bin.cgyro.ky_flux'
  character(len=18) :: binfile_ky_cflux = 'bin.cgyro.ky_cflux'
  character(len=15) :: binfile_triad   = 'bin.cgyro.triad'
  character(len=15), dimension(3) :: binfile_fieldb = &
       (/'bin.cgyro.phib ','bin.cgyro.aparb','bin.cgyro.bparb'/)
  character(len=16), dimension(3) :: binfile_kxky = &
       (/'bin.cgyro.kxky_n','bin.cgyro.kxky_e','bin.cgyro.kxky_v'/)
  character(len=19), dimension(3) :: binfile_kxky_field = &
       (/'bin.cgyro.kxky_phi ','bin.cgyro.kxky_apar','bin.cgyro.kxky_bpar'/)
  character(len=20), dimension(3) :: binfile_lky_flux = &
       (/'bin.cgyro.lky_flux_n','bin.cgyro.lky_flux_e','bin.cgyro.lky_flux_v'/)
  integer, parameter :: io=1
  ! Restart tags
  character(len=8) :: fmt='(I2.2)'
  character(len=6), dimension(100) :: rtag
  !
  ! error checking
  integer :: error_status = 0
  character(len=80) :: error_message
  !
  integer :: io_control
  integer :: signal
  integer :: restart_flag

  logical :: printout=.true.

  integer :: mpiio_small_stripe_factor   ! optional striping for data files
  integer :: mpiio_stripe_factor         ! optional striping for restart file
  character(len=2) :: mpiio_small_stripe_str
  character(len=3) :: mpiio_stripe_str
  !
  ! Precision for IO (these are set in cgyro_init_manager)
  !
  ! BYTE=4 (standard single precision)
  ! BYTE=8 (double precision for DMD analysis)
  !
  integer :: BYTE
  character(len=8)  :: fmtstr   
  integer           :: fmtstr_len 
  character(len=15) :: fmtstrn  
  character(len=9)  :: fmtstr_prec ='(es18.12)'
  !----------------------------------------------------

  !---------------------------------------------------------------
  ! Time stepping (see detailed parameters in cgyro_step)
  integer :: i_time
  integer :: n_time
  integer :: i_current
  real    :: t_current
  real, dimension(:), allocatable    :: gtime
  complex, dimension(:), allocatable :: freq
  complex :: freq_err
  integer(KIND=8) :: kernel_start_time, kernel_exit_time, kernel_count_rate, kernel_count_max
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Physics variables
  !
  integer :: ae_flag
  integer :: is_ele
  real :: dens_ele
  real :: temp_ele
  real :: mass_ele
  real :: dlnndr_ele
  real :: dlntdr_ele
  !
  real, dimension(11) :: vth  
  real, dimension(11) :: nu
  real :: rho
  real :: length
  integer :: sign_qs
  ! k_theta = k_theta_base * my_toroidal
  real :: k_theta_base
  ! omega_eb = omega_eb_base * my_toroidal
  real :: omega_eb_base
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Numerical/work arrays and dimensions
  !
  ! Velocity space 
  integer, dimension(:), allocatable :: px
  real, dimension(:), allocatable :: energy, vel, vel2, w_e
  real, dimension(:), allocatable :: xi, w_xi
  real, dimension(:,:), allocatable :: w_exi
  real, dimension(:,:), allocatable :: xi_deriv_mat, xi_lor_mat
  real, dimension(:,:), allocatable :: e_deriv1_mat, e_deriv1_rot_mat
  !
  ! Parallel streaming
  real :: d_theta
  real, dimension(:), allocatable :: theta
  real, dimension(:), allocatable :: uderiv
  real, dimension(:), allocatable :: cderiv
  complex, dimension(:,:), allocatable :: thfac_itor
  !
  ! Wavenumber advection
  integer :: source_flag
  real, parameter :: tau_ave=50.0
  real, dimension(:), allocatable :: c_wave
  complex, dimension(:,:,:), allocatable :: source
  real :: sa
  !
  ! Distributions
  complex, dimension(:,:,:,:), allocatable :: rhs
  complex, dimension(:,:,:), allocatable :: h_x
  complex, dimension(:,:,:), allocatable :: g_x
  complex, dimension(:,:,:), allocatable :: h0_x
  complex, dimension(:,:,:), allocatable :: h0_old
  complex, dimension(:,:,:,:), allocatable :: fA_nl,fB_nl
  complex, dimension(:,:,:,:), allocatable :: eA_nl,eB_nl
  complex, dimension(:,:,:,:), allocatable :: g_nl
  complex, dimension(:,:,:), allocatable :: fpackA,fpackB
  complex, dimension(:,:,:), allocatable :: epackA,epackB
  complex, dimension(:,:,:,:), allocatable :: gpack
  complex, dimension(:,:,:), allocatable :: omega_cap_h
  complex, dimension(:,:,:), allocatable :: omega_h
  complex, dimension(:,:,:), allocatable :: diss_r
  complex, dimension(:,:,:,:), allocatable :: omega_s,omega_ss
  complex, dimension(:,:,:), allocatable :: omega_sbeta
  complex, dimension(:,:,:), allocatable :: cap_h_c
  complex, dimension(:,:,:), allocatable :: cap_h_c_dot
  complex, dimension(:,:,:), allocatable :: cap_h_c_old
  complex, dimension(:,:,:), allocatable :: cap_h_c_old2
  complex, dimension(:,:,:), allocatable :: cap_h_ct
  complex, dimension(:,:,:), allocatable :: cap_h_v
  real, dimension(:,:,:,:), allocatable :: jvec_c
  real, dimension(:,:,:,:,:), allocatable :: jvec_c_nl ! used by NL only
  real, dimension(:,:,:,:), allocatable :: jvec_v
  real, dimension(:,:,:,:), allocatable :: dvjvec_c
  real, dimension(:,:,:,:), allocatable :: dvjvec_v
  real, dimension(:,:,:,:), allocatable :: jxvec_c
  real, dimension(:,:,:), allocatable :: upfac1,upfac2
  !
  ! Fields
  real, dimension(:,:,:), allocatable :: fcoef
  real, dimension(:,:,:), allocatable :: gcoef
  complex, dimension(:,:,:), allocatable :: field, field_v
  complex, dimension(:,:,:), allocatable :: field_dot
  complex, dimension(:,:,:), allocatable :: field_loc, field_loc_v
  complex, dimension(:,:,:), allocatable :: field_old
  complex, dimension(:,:,:), allocatable :: field_old2
  complex, dimension(:,:,:), allocatable :: field_old3
  complex, dimension(:,:,:,:,:), allocatable :: moment_loc
  complex, dimension(:,:,:,:,:), allocatable :: moment
  !
  ! Nonlinear fluxes (f=standard,c=central,g=global)
  real, dimension(:,:,:,:), allocatable :: cflux_loc
  real, dimension(:,:,:,:), allocatable :: cflux
  complex, dimension(:,:,:,:,:), allocatable :: gflux_loc
  complex, dimension(:,:,:,:,:), allocatable :: gflux
  complex, dimension(:,:,:,:), allocatable :: triad,triad_loc,triad_loc_old
  real, dimension(:,:), allocatable :: cflux_tave, gflux_tave
  real :: tave_min, tave_max
  integer :: tave_step
  integer :: nflux
  !
  ! Nonlinear plans
#ifndef CGYRO_GPU_FFT
  ! CPU-FFTW plans
  type(C_PTR) :: plan_r2c
  type(C_PTR) :: plan_c2r
  !
#else
  ! GPU-FFTW plans

#if defined(HIPGPU)
  type(C_PTR) :: hip_plan_r2c_manyA,hip_plan_r2c_manyB
  type(C_PTR) :: hip_plan_c2r_manyA,hip_plan_c2r_manyB,hip_plan_c2r_manyG
#elif defined(MKLGPU)
  INTEGER*8 :: dfftw_plan_r2c_manyA,dfftw_plan_r2c_manyB
  INTEGER*8 :: dfftw_plan_c2r_manyA,dfftw_plan_c2r_manyB,dfftw_plan_c2r_manyG
#else
  integer(c_int) :: cu_plan_r2c_manyA,cu_plan_r2c_manyB
  integer(c_int) :: cu_plan_c2r_manyA,cu_plan_c2r_manyB,cu_plan_c2r_manyG
#endif

  complex, dimension(:,:,:),allocatable, target :: fxmany,fymany,gxmany,gymany
  real, dimension(:,:,:), allocatable, target :: uxmany,uymany
  real, dimension(:,:,:), allocatable, target :: vxmany,vymany,uvmany
#endif
  ! 
  ! 2D FFT dimensions 
  integer :: nx,ny
  integer :: nx0,ny0
  integer :: nx2,ny2
  !
  ! 2D FFT work arrays
#ifndef CGYRO_GPU_FFT
  real, dimension(:,:,:), allocatable :: vxmany
  real, dimension(:,:,:), allocatable :: vymany
  real, dimension(:,:,:), allocatable :: uxmany
  real, dimension(:,:,:), allocatable :: uymany
  real, dimension(:,:,:), allocatable :: uv
  complex, dimension(:,:,:),allocatable :: fx
  complex, dimension(:,:,:),allocatable :: fy
  complex, dimension(:,:,:),allocatable :: gx
  complex, dimension(:,:,:),allocatable :: gy
#endif
  !
  ! Work arrays
  real, dimension(2) :: integration_error
  ! Upwind work arrays
  complex, dimension(:,:,:),allocatable :: upwind_res_loc
  complex, dimension(:,:,:),allocatable :: upwind_res
  complex(KIND=REAL32), dimension(:,:,:),allocatable :: upwind32_res_loc
  complex(KIND=REAL32), dimension(:,:,:),allocatable :: upwind32_res
  !
  ! LAPACK work arrays 
  real, dimension(:), allocatable :: work  
  integer, dimension(:), allocatable :: i_piv
  integer :: info
  !
  ! Field solve variables
  real, dimension(:), allocatable :: sum_den_h
  real, dimension(:,:), allocatable :: sum_den_x,sum_cur_x
  real, dimension(:), allocatable :: vfac
  !
  ! n=0 test variables
  real, dimension(:,:,:), allocatable :: hzf,xzf 
  !
  ! Collision operator
  integer :: n_low_energy
  real, dimension(:,:,:,:), allocatable :: cmat ! only used if collision_precision_mode=0 & 64
  real(KIND=REAL32), dimension(:,:,:,:), allocatable :: cmat_fp32 ! only used if collision_precision_mode=1 & 32
  real(KIND=REAL32), dimension(:,:,:,:,:,:), allocatable :: cmat_stripes ! only used if collision_precision_mode=1
  real(KIND=REAL32), dimension(:,:,:,:,:,:), allocatable :: cmat_e1 ! only used if collision_precision_mode=1
  real, dimension(:,:,:,:,:,:), allocatable :: cmat_simple ! only used in collision_model=5
  ! 
  ! Equilibrium/geometry arrays
  integer :: it0
  integer, dimension(:), allocatable :: itp
  real :: bigr_th0
  real :: bigr_r_th0
  real, dimension(:,:), allocatable :: thetab
  real, dimension(:), allocatable   :: w_theta
  real, dimension(:), allocatable   :: g_theta
  real, dimension(:), allocatable   :: g_theta_geo
  real, dimension(:,:), allocatable :: k_perp
  real, dimension(:,:), allocatable :: k_x
  real, dimension(:), allocatable   :: bmag
  real, dimension(:), allocatable   :: btor
  real, dimension(:), allocatable   :: bpol
  real, dimension(:), allocatable   :: bigr
  real, dimension(:), allocatable   :: bigr_r
  real, dimension(:), allocatable   :: captheta
  real, dimension(:,:,:), allocatable :: omega_stream
  real, dimension(:,:,:), allocatable :: omega_trap
  real, dimension(:,:), allocatable :: omega_rdrift
  real, dimension(:,:), allocatable :: omega_adrift
  real, dimension(:,:), allocatable :: omega_aprdrift
  real, dimension(:,:), allocatable :: omega_cdrift
  real, dimension(:,:), allocatable :: omega_cdrift_r
  real, dimension(:),   allocatable :: omega_gammap
  integer, parameter                :: n_beta_star=2
  real, dimension(0:n_beta_star)    :: beta_star
  real                              :: beta_star_fac
  real                              :: mach_one_fac
  ! for centrifugal rotation
  real, dimension(:,:), allocatable :: lambda_rot
  real, dimension(:,:), allocatable :: dlambda_rot
  real, dimension(:,:), allocatable :: dens_rot,dens2_rot
  real, dimension(:),   allocatable :: dens_ele_rot
  real, dimension(:),   allocatable :: dens_avg_rot
  real, dimension(:),   allocatable :: dlnndr_avg_rot
  real, dimension(:,:), allocatable :: omega_rot_trap
  real, dimension(:,:), allocatable :: omega_rot_u
  real, dimension(:,:), allocatable :: omega_rot_drift
  real, dimension(:,:), allocatable :: omega_rot_drift_r
  real, dimension(:),   allocatable :: omega_rot_edrift
  real, dimension(:),   allocatable :: omega_rot_edrift_r
  real, dimension(:,:), allocatable :: omega_rot_star
  !
  ! Number of gridpoints for Miller geometry integration grid
  integer, parameter :: geo_ntheta=1001 
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  integer :: geo_ny_in
  real, dimension(8,0:32) :: geo_yin_in
  !
  integer :: geo_numeq_flag
  integer :: geo_ny
  real, dimension(:,:), allocatable :: geo_yin
  !---------------------------------------------------------------

  real :: total_memory
  real :: small

end module cgyro_globals
