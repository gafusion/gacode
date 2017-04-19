!-----------------------------------------------------------------
! cgyro_globals.f90
!
! PURPOSE:
!  CGYRO global variables.  The idea is to have a primary, large
!  module containing all essential CGYRO arrays and scalars.
!-----------------------------------------------------------------
#ifdef _OPENACC
#include "precision_m.f90"
#include "cufft_m.f90"
#endif

module cgyro_globals

  use, intrinsic :: iso_c_binding

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
  integer :: e_method
  real    :: delta_t
  real    :: max_time
  integer :: print_step
  integer :: restart_step
  real    :: freq_tol
  integer :: restart_mode
  real    :: up_radial
  real    :: up_theta
  real    :: up_alpha
  real    :: up_wave
  integer :: nup_radial
  integer :: nup_theta
  integer :: nup_alpha
  integer :: nup_wave
  integer :: constant_stream_flag
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
  real    :: collision_ele_scale
  integer :: zf_test_flag 
  integer :: nonlinear_flag 
  integer :: nonlinear_method
  real :: te_ade
  real :: ne_ade
  real :: dlntdre_ade   
  real :: dlnndre_ade   
  real :: masse_ade
  real :: lambda_star
  integer :: test_flag
  integer :: h_print_flag
  integer :: moment_print_flag
  integer :: kxkyflux_print_flag
  integer :: globalflux_print_flag
  real :: amp0
  real :: amp
  real :: gamma_e
  real :: gamma_p
  real :: mach
  integer :: rotation_model
  real :: error_tol
  integer :: mpi_rank_order
  integer :: hiprec_flag
  integer :: udsymmetry_flag
  integer :: shear_method
  integer :: n_global
  integer :: psym_flag
  integer :: profile_shear_flag
  real :: adapt_tol
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
  real :: betae_unit
  !
  ! Species parameters
  !
  integer :: n_species
  real :: nu_ee
  real, dimension(6) :: z
  real, dimension(6) :: mass
  real, dimension(6) :: dens
  real, dimension(6) :: temp
  real, dimension(6) :: dlnndr
  real, dimension(6) :: dlntdr
  real, dimension(6) :: sdlnndr
  real, dimension(6) :: sdlntdr

  integer :: subroutine_flag  ! only used for cgyro_read_input

  ! Re-scaling parameters for experimental profiles
  real :: lambda_star_scale
  real :: gamma_e_scale
  real :: gamma_p_scale
  real :: mach_scale
  real :: q_scale
  real :: s_scale
  real :: shift_scale
  real :: kappa_scale, s_kappa_scale
  real :: delta_scale, s_delta_scale
  real :: zeta_scale, s_zeta_scale
  real :: beta_star_scale, betae_unit_scale
  real :: nu_ee_scale
  real, dimension(6) :: dlnndr_scale
  real, dimension(6) :: dlntdr_scale

  real :: lambda_debye
  real :: rhos

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
  integer :: n_proc
  integer :: n_proc_1
  integer :: n_proc_2
  integer :: i_group_1
  integer :: i_group_2
  integer :: CGYRO_COMM_WORLD
  integer :: NEW_COMM_1
  integer :: NEW_COMM_2
  integer :: nv1,nv2,nc1,nc2
  integer :: nsplit
  integer, dimension(:), allocatable :: recv_status
  !
  ! Pointers
  integer :: nv,iv
  integer :: nv_loc,iv_loc
  integer :: nc,ic
  integer :: nc_loc,ic_loc
  integer, dimension(:), allocatable :: ie_v
  integer, dimension(:), allocatable :: ix_v
  integer, dimension(:), allocatable :: is_v
  integer, dimension(:), allocatable :: ir_c
  integer, dimension(:), allocatable :: it_c
  integer, dimension(:,:), allocatable :: ic_c
  integer, dimension(:,:,:), allocatable :: iv_v
  integer, dimension(:), allocatable :: ica_c,icb_c
  !
  integer :: n
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Constants
  !
  real, parameter    :: pi   = 3.1415926535897932
  complex, parameter :: i_c  = (0.0,1.0)
  real, parameter    :: num1 = 1.0
  real, parameter    :: num0 = 0.0
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! I/O and error management variables
  !
  character(len=80) :: path
  character(len=18) :: runfile_info    = 'out.cgyro.info'
  character(len=18) :: runfile_mpi     = 'out.cgyro.mpi'
  character(len=18) :: runfile_extend  = 'out.cgyro.extend'
  character(len=18) :: runfile_memory  = 'out.cgyro.memory'
  character(len=18) :: runfile_restart = 'out.cgyro.restart'
  character(len=18) :: runfile_restart_tag = 'out.cgyro.tag'
  character(len=18) :: runfile_hb      = 'out.cgyro.hb'
  character(len=18) :: runfile_grids   = 'out.cgyro.grids'
  character(len=18) :: runfile_prec    = 'out.cgyro.prec'
  character(len=18) :: runfile_time    = 'out.cgyro.time'
  character(len=18) :: runfile_timers  = 'out.cgyro.timing'
  character(len=18) :: runfile_freq    = 'out.cgyro.freq'
  character(len=18) :: runfile_kxky_phi = 'out.cgyro.kxky_phi'
  character(len=21) :: runfile_kxky_flux = 'out.cgyro.kxky_flux_e'
  character(len=19) :: runfile_ky_flux = 'out.cgyro.ky_flux'
  character(len=15), dimension(3)  :: runfile_fieldb = &
       (/'out.cgyro.phib ','out.cgyro.aparb','out.cgyro.bparb'/)
  character(len=16), dimension(2)  :: runfile_kxky = &
       (/'out.cgyro.kxky_n','out.cgyro.kxky_e'/)
  character(len=20), dimension(2)  :: runfile_lky_flux = &
       (/'out.cgyro.lky_flux_n','out.cgyro.lky_flux_e'/)
  integer, parameter :: io=1
  !
  ! error checking
  integer :: error_status = 0
  character(len=80) :: error_message
  !
  integer :: io_control
  integer :: signal
  integer :: restart_flag
  !
  ! Standard precision for IO (there are optionally reset to higher precision later)
  character(len=8)  :: fmtstr    ='(es11.4)'
  character(len=14) :: fmtstr2   ='(2(es11.4,1x))'
  character(len=15) :: fmtstrn   ='(10(es11.4,1x))'
  character(len=9)  :: fmtstr_hi ='(es18.12)'
  !----------------------------------------------------

  !---------------------------------------------------------------
  ! Time stepping
  integer :: i_time
  integer :: n_time
  integer :: i_current
  real :: t_current
  complex :: freq
  complex :: freq_err
  real :: gtime
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
  real, dimension(6) :: vth  
  real, dimension(6) :: nu
  real :: rho
  real :: k_theta
  real :: length
  real :: omega_eb
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Numerical/work arrays and dimensions
  !
  ! Velocity space 
  integer, dimension(:), allocatable :: indx_xi, px
  real, dimension(:), allocatable :: energy, vel, w_e
  real, dimension(:), allocatable :: xi, w_xi
  real, dimension(:,:), allocatable :: xi_deriv_mat, xi_lor_mat
  real, dimension(:,:), allocatable :: e_deriv1_mat, e_deriv2_mat
  !
  ! Parallel streaming
  real :: d_theta
  real, dimension(:), allocatable :: theta
  real, dimension(:), allocatable :: uderiv
  real, dimension(:), allocatable :: cderiv
  integer, dimension(:,:), allocatable :: icd_c
  complex, dimension(:,:), allocatable :: dtheta
  complex, dimension(:,:), allocatable :: dtheta_up
  !
  ! Wavenumber advection
  real, dimension(:), allocatable :: der_wave
  real, dimension(:), allocatable :: dis_wave
  !
  ! Distributions
  complex, dimension(:,:,:), allocatable :: rhs
  complex, dimension(:,:), allocatable :: h_x
  complex, dimension(:,:), allocatable :: g_x
  complex, dimension(:,:), allocatable :: h0_x
  complex, dimension(:,:), allocatable :: psi
  complex, dimension(:,:), allocatable :: chi
  complex, dimension(:,:,:), allocatable :: f_nl
  complex, dimension(:,:,:), allocatable :: g_nl
  complex, dimension(:,:), allocatable :: omega_cap_h
  complex, dimension(:,:), allocatable :: omega_h
  complex, dimension(:,:,:), allocatable :: omega_s,omega_ss
  complex, dimension(:,:), allocatable :: cap_h_c
  complex, dimension(:,:), allocatable :: cap_h_ct
  complex, dimension(:,:), allocatable :: cap_h_v
  complex, dimension(:,:), allocatable :: cap_h_v_prime
  real, dimension(:,:,:), allocatable :: jvec_c
  real, dimension(:,:,:), allocatable :: jvec_v
  real, dimension(:,:,:), allocatable :: dvjvec_c
  real, dimension(:,:,:), allocatable :: dvjvec_v
  real, dimension(:,:,:), allocatable :: jxvec_c
  real, dimension(:,:), allocatable :: upfac1,upfac2
  !
  ! Fields
  real, dimension(:,:), allocatable :: fcoef
  real, dimension(:,:), allocatable :: gcoef
  real, dimension(:,:), allocatable :: res_norm
  complex, dimension(:,:), allocatable :: field
  complex, dimension(:,:), allocatable :: field_loc
  complex, dimension(:,:), allocatable :: field_old
  complex, dimension(:,:), allocatable :: field_old2
  complex, dimension(:,:), allocatable :: field_old3
  complex, dimension(:,:,:), allocatable :: moment_loc
  complex, dimension(:,:,:), allocatable :: moment
  !
  ! Nonlinear fluxes 
  real, dimension(:,:), allocatable :: flux_loc
  real, dimension(:,:), allocatable :: flux
  real, dimension(:,:,:), allocatable :: fflux_loc
  real, dimension(:,:,:), allocatable :: fflux
  complex, dimension(:,:,:), allocatable :: gflux_loc
  complex, dimension(:,:,:), allocatable :: gflux
  !
  ! Nonlinear plans
  type(C_PTR) :: plan_r2c
  type(C_PTR) :: plan_c2r
  !
  ! GPU-FFTW plans
#ifdef _OPENACC
  integer(c_int) :: cu_plan_r2c_many
  integer(c_int) :: cu_plan_c2r_many
  complex, dimension(:,:,:),allocatable :: fxmany,fymany,gxmany,gymany
  real, dimension(:,:,:), allocatable :: uxmany,uymany
  real, dimension(:,:,:), allocatable :: vxmany,vymany,uvmany
#endif
  ! 
  ! 2D FFT dimensions 
  integer :: nx,ny
  integer :: nx0,ny0
  !
  ! 2D FFT work arrays
  real, dimension(:,:), allocatable :: ux
  real, dimension(:,:), allocatable :: uy
  real, dimension(:,:), allocatable :: vx
  real, dimension(:,:), allocatable :: vy
  real, dimension(:,:), allocatable :: uv
  complex, dimension(:,:),allocatable :: fx
  complex, dimension(:,:),allocatable :: fy
  complex, dimension(:,:),allocatable :: gx
  complex, dimension(:,:),allocatable :: gy
  !
  ! Work arrays
  complex, dimension(:,:), allocatable :: f_balloon
  real, dimension(2) :: integration_error
  !
  ! LAPACK work arrays 
  real, dimension(:), allocatable :: work  
  integer, dimension(:), allocatable :: i_piv
  integer :: info
  !
  ! Field solve variables
  real, dimension(:), allocatable :: sum_den_h
  real, dimension(:), allocatable :: sum_den_x, sum_cur_x
  real, dimension(:), allocatable :: vfac
  !
  ! n=0 test variables
  real, dimension(:,:,:), allocatable :: hzf, xzf 
  real, dimension(:), allocatable :: pvec_outr, pvec_outi
  !
  ! Collision operator
  real, dimension(:,:,:), allocatable :: cmat
  real, dimension(:,:,:,:,:), allocatable :: cmat_simple
  ! 
  ! Equilibrium/geometry arrays
  integer :: it0
  real :: bigR_th0
  real :: bigR_r_th0
  real, dimension(:,:), allocatable :: thetab
  real, dimension(:), allocatable   :: w_theta
  real, dimension(:), allocatable   :: g_theta
  real, dimension(:), allocatable   :: g_theta_geo
  real, dimension(:), allocatable   :: k_perp
  real, dimension(:), allocatable   :: k_x
  real, dimension(:), allocatable   :: bmag
  real, dimension(:), allocatable   :: btor
  real, dimension(:), allocatable   :: bpol
  real, dimension(:), allocatable   :: bigR
  real, dimension(:), allocatable   :: bigR_r
  real, dimension(:,:), allocatable :: omega_stream
  real, dimension(:,:), allocatable :: omega_trap
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
  real, dimension(:,:), allocatable :: dens_rot
  real, dimension(:),   allocatable :: dens_ele_rot
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
  
end module cgyro_globals
