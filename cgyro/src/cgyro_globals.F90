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
  integer :: hyper_flag
  integer :: hyper_order
  real    :: hyper_coeff
  integer :: n_wave
  integer :: constant_stream_flag
  integer :: explicit_trap_flag
  real    :: ky
  integer :: box_size
  real    :: ipccw
  real    :: btccw
  integer :: silent_flag
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
  integer :: z_eff_method
  real    :: z_eff
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
  integer :: momentum_print_flag
  real :: amp0
  real :: amp
  real :: gamma_e
  real :: gamma_p
  real :: mach
  integer :: rotation_model
  integer :: nt_loc
  integer :: mpi_rank_order
  integer :: velocity_order
  integer :: hiprec_flag
  integer :: udsymmetry_flag
  integer :: shear_method
  integer :: global_flag
  integer :: n_global
  real    :: nu_global
  integer :: theta_plot
  integer :: gpu_bigmem_flag
  integer :: upwind_single_flag
  integer :: nl_single_flag
  real :: px0
  integer :: exch_flag
  integer :: nl_min
  real :: dealias
  integer :: dealias_order
  integer :: scale_index
  real :: scale_stream
  real :: scale_trap
  real :: scale_drift
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
  real :: beta_star_scale
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

  real :: lambda_debye
  real :: b_gs2
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! MPI/OpenMP variables and pointers
  ! 
  integer :: n_omp
  logical :: have_COMM_4 = .FALSE.
  !
  integer :: i_err
  integer :: i_proc
  integer :: i_proc_1
  integer :: i_proc_2
  integer :: i_proc_3
  integer :: i_proc_4
  integer :: n_proc
  integer :: n_proc_1
  integer :: n_proc_2
  integer :: n_proc_4
  integer :: CGYRO_COMM_WORLD
  integer :: NEW_COMM_1  ! simple Linear
  integer :: NEW_COMM_2  ! non-linear
  integer :: NEW_COMM_3  ! simple linear by species
  integer :: CGYRO_COMM_WORLD_4
  integer :: NEW_COMM_4  ! aggregate linear, coll
  integer :: nv1,nv2,nc1,nc2,nc_cl1,nc_cl2
  integer :: nsplit,nsplitA,nsplitB
  integer :: ns1,ns2
  integer, dimension(:), allocatable :: recv_status
  integer :: fA_req, fB_req, g_req
  logical :: fA_req_valid, fB_req_valid, g_req_valid
  ! Thetas present in the process after NL AllToAll
  integer :: n_jtheta
  !
  ! Pointers
  integer :: nv,iv
  integer :: nv_loc,iv_loc
  integer :: nc,ic
  integer :: nc_loc,ic_loc
  integer :: nc_loc_coll   ! nc local slice for coll purposes, equals nc_loc for single simulation
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
  integer :: nt1,nt2
  integer :: n_toroidal_procs
  ! For multi-simulation handling
  integer :: n_sim          ! how many simulations is coll processing at once
  integer :: i_sim          ! 0-based order in n_sim
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
  character(len=21) :: binfile_ky_flux_mom = 'bin.cgyro.ky_flux_mom'
  character(len=22) :: binfile_ky_cflux_mom = 'bin.cgyro.ky_cflux_mom'
  character(len=15), dimension(4) :: binfile_fieldb = &
       (/'bin.cgyro.phib ','bin.cgyro.aparb','bin.cgyro.bparb','bin.cgyro.eparb'/)
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
  integer :: i_at_restart
  real    :: t_at_restart
  !---------------------------------------------------------------

  !---------------------------------------------------------------
  ! Numerical/work arrays and dimensions
  !
  ! Velocity space 
  integer :: px_zero ! Hostorical: used to be px(ir) := ir - px_zero
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
  real :: dsrc
  !
  ! Distributions
  complex, dimension(:,:,:,:), allocatable :: rhs
  complex, dimension(:,:,:), allocatable :: h_x
  complex, dimension(:,:,:), allocatable :: g_x
  complex, dimension(:,:,:), allocatable :: h0_x
  complex, dimension(:,:,:), allocatable :: h0_old
  complex, dimension(:,:,:,:), allocatable :: fA_nl,fB_nl
  complex(KIND=REAL32), dimension(:,:,:,:), allocatable :: fA_nl32,fB_nl32
  complex, dimension(:,:,:,:), allocatable :: g_nl
  complex(KIND=REAL32), dimension(:,:,:,:), allocatable :: g_nl32
  complex, dimension(:,:,:), allocatable :: fpackA,fpackB
  complex(KIND=REAL32), dimension(:,:,:), allocatable :: fpackA32,fpackB32
  complex, dimension(:,:,:,:), allocatable :: gpack
  complex(KIND=REAL32), dimension(:,:,:,:), allocatable :: gpack32
  complex, dimension(:,:,:), allocatable :: omega_cap_h
  complex, dimension(:,:,:), allocatable :: omega_h
  complex, dimension(:,:,:,:), allocatable :: omega_s,omega_ss
  complex, dimension(:,:,:), allocatable :: omega_sbeta
  complex, dimension(:,:,:), allocatable :: cap_h_c
  complex, dimension(:,:,:), allocatable :: cap_h_c_dot
  complex, dimension(:,:,:), allocatable :: cap_h_c_old
  complex, dimension(:,:,:), allocatable :: cap_h_c_old2
  complex, dimension(:,:,:), allocatable :: cap_h_ct
  complex, dimension(:,:,:,:), allocatable :: cap_h_v
  real, dimension(:,:,:,:), allocatable :: jvec_c
  real, dimension(:,:,:,:,:), allocatable :: jvec_c_nl ! used by NL only
  real(KIND=REAL32), dimension(:,:,:,:,:), allocatable :: jvec_c_nl32 ! used by NL only
  real, dimension(:,:,:,:,:), allocatable :: jvec_v
  real, dimension(:,:,:,:), allocatable :: dvjvec_c
  real, dimension(:,:,:,:), allocatable :: dvjvec_v
  real, dimension(:,:,:,:), allocatable :: jxvec_c
  real, dimension(:,:,:), allocatable :: upfac1,upfac2
  real, dimension(:,:,:,:), allocatable :: jmvec_c
  !
  ! Fields
  real, dimension(:,:,:), allocatable :: fcoef
  real, dimension(:,:,:), allocatable :: gcoef
  complex, dimension(:,:,:), allocatable :: field
  complex, dimension(:,:,:), allocatable :: field_dot
  complex, dimension(:,:,:), allocatable :: field_loc
  complex, dimension(:,:,:), allocatable :: field_old
  complex, dimension(:,:,:), allocatable :: field_old2
  complex, dimension(:,:,:), allocatable :: field_old3
  complex, dimension(:,:), allocatable :: epar
  complex, dimension(:,:,:,:,:), allocatable :: moment_loc
  complex, dimension(:,:,:,:,:), allocatable :: moment
  complex, dimension(:,:,:,:), allocatable :: field_v
  complex, dimension(:,:,:,:), allocatable :: field_loc_v
  !
  ! Nonlinear fluxes (f=standard,c=central,g=global)
  real, dimension(:,:,:,:), allocatable :: cflux_loc
  real, dimension(:,:,:,:), allocatable :: cflux
  complex, dimension(:,:,:,:,:), allocatable :: gflux_loc
  complex, dimension(:,:,:,:,:), allocatable :: gflux
  real, dimension(:,:), allocatable :: cflux_tave, gflux_tave
  real, dimension(:,:,:,:), allocatable :: cflux_mom_loc
  real, dimension(:,:,:,:), allocatable :: cflux_mom
  complex, dimension(:,:,:,:,:), allocatable :: gflux_mom_loc
  complex, dimension(:,:,:,:,:), allocatable :: gflux_mom
  real :: tave_min, tave_max
  integer :: tave_step
  integer :: nflux
  integer :: nflux_mom
  ! NL dealias variables
  !integer, dimension(:,:), allocatable :: dealias_pvec_count
  !integer, dimension(:,:,:), allocatable :: dealias_pvec
  integer, dimension(:,:,:,:), allocatable :: dealias_raw_ir
  integer, dimension(:,:,:,:), allocatable :: dealias_raw_it
  complex, dimension(:,:,:,:), allocatable :: dealias_raw_ph
  complex, dimension(:,:,:,:), allocatable :: inraw_dealias
  complex, dimension(:,:,:,:), allocatable :: outraw_dealias
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
  type(C_PTR)    :: plan_r2c_manyA,plan_r2c_manyB
  type(C_PTR)    :: plan_c2r_manyA,plan_c2r_manyB,plan_c2r_manyG
#elif defined(MKLGPU)
  INTEGER*8      :: plan_r2c_manyA,plan_r2c_manyB
  INTEGER*8      :: plan_c2r_manyA,plan_c2r_manyB,plan_c2r_manyG
#else
  integer(c_int) :: plan_r2c_manyA,plan_r2c_manyB
  integer(c_int) :: plan_c2r_manyA,plan_c2r_manyB,plan_c2r_manyG
#endif

  complex, dimension(:,:,:),allocatable, target :: fxmany,fymany,gxmany,gymany
  real, dimension(:,:,:), allocatable, target :: uxmany,uymany
  real, dimension(:,:,:), allocatable, target :: vxmany,vymany,uvmany

  complex(KIND=REAL32), dimension(:,:,:),allocatable, target :: fxmany32,fymany32,gxmany32,gymany32
  real(KIND=REAL32), dimension(:,:,:), allocatable, target :: uxmany32,uymany32
  real(KIND=REAL32), dimension(:,:,:), allocatable, target :: vxmany32,vymany32,uvmany32
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

  real(KIND=REAL32), dimension(:,:,:), allocatable :: vxmany32
  real(KIND=REAL32), dimension(:,:,:), allocatable :: vymany32
  real(KIND=REAL32), dimension(:,:,:), allocatable :: uxmany32
  real(KIND=REAL32), dimension(:,:,:), allocatable :: uymany32
  real(KIND=REAL32), dimension(:,:,:), allocatable :: uv32
  complex(KIND=REAL32), dimension(:,:,:),allocatable :: fx32
  complex(KIND=REAL32), dimension(:,:,:),allocatable :: fy32
  complex(KIND=REAL32), dimension(:,:,:),allocatable :: gx32
  complex(KIND=REAL32), dimension(:,:,:),allocatable :: gy32
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

  ! Unlike most other arrays, last dimension is simply nt_loc
  ! due to limis of pointer handling
  real, dimension(:,:,:,:), pointer :: cmat ! only used if collision_precision_mode=0 & 64
  real(KIND=REAL32), dimension(:,:,:,:), pointer :: cmat_fp32 ! only used if collision_precision_mode=1 & 32
  !
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

  real :: total_memory
  real :: small

contains

!
! Avoid having two copies of cmat, one in CPU and one in GPU memory
! On GPU systems under OMPGPU, keep it only in GPU memory (omp_target_alloc)
! ON CPU-only systems and OpenACC, keep the old logic
!

#ifdef OMPGPU
subroutine allocate_cmat_gpu_fp32(nv,nc_loc_coll,nt1,nt2)
    use iso_c_binding
    use omp_lib
    implicit none

    ! ----------------------
    integer, intent(in)          :: nv,nc_loc_coll,nt1,nt2

    ! Pointer to array, only need raw pointer locally, can recover
    type(c_ptr) :: c_ptr_buffer
    integer :: dev_id
    ! Use c_sizeof to avoid integer wraparound
    integer(c_size_t) :: total_bytes
    integer(c_size_t) :: nv2
    integer :: nt_loc

    nv2 = nv*nv
    nt_loc = nt2-nt1+1

    total_bytes = (nv2*nc_loc_coll)*nt_loc*4 ! fp32 =  4 bytes

    ! OMP supports multi-GPU setups, we only support 1-GPU ones
    dev_id = omp_get_default_device()
    c_ptr_buffer = omp_target_alloc(total_bytes, dev_id)

    if (.not. c_associated(c_ptr_buffer)) then
         ! Catastrophic error, do not even try to process it cleanly
         write(*,*) "Error: allocate_cmat_fp32 failed to allocate memory."
        stop
    end if

    call c_f_pointer(c_ptr_buffer, cmat_fp32, [nv, nv, nc_loc_coll, nt_loc])

end subroutine allocate_cmat_gpu_fp32

subroutine allocate_cmat_gpu(nv,nc_loc_coll,nt1,nt2)
    use iso_c_binding
    use omp_lib
    implicit none

    ! ----------------------
    integer, intent(in)          :: nv,nc_loc_coll,nt1,nt2

    ! Pointer to array, only need raw pointer locally, can recover
    type(c_ptr) :: c_ptr_buffer
    integer :: dev_id
    ! Use c_sizeof to avoid integer wraparound
    integer(c_size_t) :: total_bytes
    integer(c_size_t) :: nv2
    integer :: nt_loc

    nv2 = nv*nv
    nt_loc = nt2-nt1+1

    total_bytes = (nv2*nc_loc_coll)*nt_loc*8 ! fp64 =  8 bytes

    ! OMP supports multi-GPU setups, we only support 1-GPU ones
    dev_id = omp_get_default_device()
    c_ptr_buffer = omp_target_alloc(total_bytes, dev_id)

    if (.not. c_associated(c_ptr_buffer)) then
         ! Catastrophic error, do not even try to process it cleanly
         write(*,*) "Error: allocate_cmat failed to allocate memory."
        stop
    end if

    call c_f_pointer(c_ptr_buffer, cmat, [nv, nv, nc_loc_coll, nt_loc])

end subroutine allocate_cmat_gpu

subroutine deallocate_cmat_gpu_fp32
    use iso_c_binding
    use omp_lib
    implicit none
    type(c_ptr) :: c_ptr_buffer
    integer :: dev_id

    if (associated(cmat_fp32)) then
      c_ptr_buffer = c_loc(cmat_fp32(1,1,1,1))

      ! OMP supports multi-GPU setups, we only support 1-GPU ones
      ! Note: Must be the same as the one used in alloc
      dev_id = omp_get_default_device()
      call omp_target_free(c_ptr_buffer, dev_id)
    endif

end subroutine
        
subroutine deallocate_cmat_gpu
    use iso_c_binding
    use omp_lib
    implicit none
    type(c_ptr) :: c_ptr_buffer
    integer :: dev_id

    if (associated(cmat)) then
      c_ptr_buffer = c_loc(cmat(1,1,1,1))

      ! OMP supports multi-GPU setups, we only support 1-GPU ones
      ! Note: Must be the same as the one used in alloc
      dev_id = omp_get_default_device()
      call omp_target_free(c_ptr_buffer, dev_id)
    endif

end subroutine

subroutine copy_into_cmat_gpu_fp32(amat, ic_loc, itor)
    use iso_c_binding
    use omp_lib
    implicit none
    
    real(KIND=REAL32), target, intent(in)    :: amat(nv, nv)
    integer, intent(in) :: ic_loc, itor
    
    integer(c_size_t)   :: bytes_to_copy, dst_offset
    integer             :: host_id, device_id
    integer             :: ierr
    
    ! Size of the (nv, nv) block in bytes
    bytes_to_copy = int(nv, c_size_t) * int(nv, c_size_t) * c_sizeof(amat(1,1))
    
    host_id   = omp_get_initial_device()
    device_id = omp_get_default_device()
    
    ! 3. Calculate destination offset (0-based bytes)
    ! We are targeting: cmat_fp32(1, 1, ic_loc, itor-nt1+1)
    ! Offset = [(dim4_idx - 1) * (size_dim3 * size_dim2 * size_dim1) + (dim3_idx - 1) * (size_dim2 * size_dim1)] * element_size
    dst_offset = ( int(itor - nt1, c_size_t) * int(nc_loc_coll, c_size_t) * int(nv * nv, c_size_t) + &
                   int(ic_loc - 1, c_size_t) * int(nv * nv, c_size_t) ) * c_sizeof(amat(1,1))

    ierr = omp_target_memcpy( &
        c_loc(cmat_fp32),    & ! Destination (device pointer)
        c_loc(amat),        & ! Source (host pointer)
        bytes_to_copy,      & ! Length in bytes
        dst_offset,         & ! Offset in destination
        0_c_size_t,         & ! Offset in source
        device_id,          & ! Destination device
        host_id             & ! Source device
    )

    if (ierr /= 0) then
        ! Catastrophic error, do not even try to process it cleanly
        write(*,*) "Error in omp_target_memcpy: ", ierr
        stop
    end if

end subroutine

subroutine copy_into_cmat_gpu(amat, ic_loc, itor)
    use iso_c_binding
    use omp_lib
    implicit none
    
    real, target, intent(in)    :: amat(nv, nv)
    integer, intent(in) :: ic_loc, itor
    
    integer(c_size_t)   :: bytes_to_copy, dst_offset
    integer             :: host_id, device_id
    integer             :: ierr
    
    ! Size of the (nv, nv) block in bytes
    bytes_to_copy = int(nv, c_size_t) * int(nv, c_size_t) * c_sizeof(amat(1,1))
    
    host_id   = omp_get_initial_device()
    device_id = omp_get_default_device()
    
    ! 3. Calculate destination offset (0-based bytes)
    ! We are targeting: cmat(1, 1, ic_loc, itor-nt1+1)
    ! Offset = [(dim4_idx - 1) * (size_dim3 * size_dim2 * size_dim1) + (dim3_idx - 1) * (size_dim2 * size_dim1)] * element_size
    dst_offset = ( int(itor - nt1, c_size_t) * int(nc_loc_coll, c_size_t) * int(nv * nv, c_size_t) + &
                   int(ic_loc - 1, c_size_t) * int(nv * nv, c_size_t) ) * c_sizeof(amat(1,1))

    ierr = omp_target_memcpy( &
        c_loc(cmat),        & ! Destination (device pointer)
        c_loc(amat),        & ! Source (host pointer)
        bytes_to_copy,      & ! Length in bytes
        dst_offset,         & ! Offset in destination
        0_c_size_t,         & ! Offset in source
        device_id,          & ! Destination device
        host_id             & ! Source device
    )

    if (ierr /= 0) then
        ! Catastrophic error, do not even try to process it cleanly
        write(*,*) "Error in omp_target_memcpy: ", ierr
        stop
    end if

end subroutine

subroutine copy_from_cmat_gpu(amat, ic_loc, itor)
    use iso_c_binding
    use omp_lib
    implicit none
    
    real, target, intent(inout) :: amat(nv, nv)
    integer, intent(in) :: ic_loc, itor
    
    integer(c_size_t)   :: bytes_to_copy, src_offset
    integer             :: host_id, device_id
    integer             :: ierr
    
    ! Size of the (nv, nv) block in bytes
    bytes_to_copy = int(nv, c_size_t) * int(nv, c_size_t) * c_sizeof(amat(1,1))
    
    host_id   = omp_get_initial_device()
    device_id = omp_get_default_device()
    
    ! 3. Calculate cmat offset (0-based bytes)
    ! We are targeting: cmat(1, 1, ic_loc, itor-nt1+1)
    ! Offset = [(dim4_idx - 1) * (size_dim3 * size_dim2 * size_dim1) + (dim3_idx - 1) * (size_dim2 * size_dim1)] * element_size
    src_offset = ( int(itor - nt1, c_size_t) * int(nc_loc_coll, c_size_t) * int(nv * nv, c_size_t) + &
                   int(ic_loc - 1, c_size_t) * int(nv * nv, c_size_t) ) * c_sizeof(amat(1,1))

    ierr = omp_target_memcpy( &
        c_loc(amat),        & ! Destination (host pointer)
        c_loc(cmat),        & ! Source (device pointer)
        bytes_to_copy,      & ! Length in bytes
        0_c_size_t,         & ! Offset in destination
        src_offset,         & ! Offset in source
        host_id,            & ! Destination device
        device_id           & ! Source device
    )

    if (ierr /= 0) then
        ! Catastrophic error, do not even try to process it cleanly
        write(*,*) "Error in omp_target_memcpy: ", ierr
        stop
    end if

end subroutine

subroutine copy_from_cmat_gpu_all(cmat_cpu)
    use iso_c_binding
    use omp_lib
    implicit none
    
    real, target, intent(inout) :: cmat_cpu(nv, nv, nc_loc_coll, nt_loc)
    
    integer(c_size_t)   :: bytes_to_copy
    integer             :: host_id, device_id
    integer             :: ierr
    
    ! Size of the cmat block in bytes
    bytes_to_copy = int(nv, c_size_t) * int(nv, c_size_t) * &
                    int(nc_loc_coll, c_size_t) * int(nt_loc, c_size_t) * &
                    c_sizeof(cmat_cpu(1,1,1,1))
    
    host_id   = omp_get_initial_device()
    device_id = omp_get_default_device()
    
    ierr = omp_target_memcpy( &
        c_loc(cmat_cpu),    & ! Destination (host pointer)
        c_loc(cmat),        & ! Source (device pointer)
        bytes_to_copy,      & ! Length in bytes
        0_c_size_t,         & ! Offset in destination
        0_c_size_t,         & ! Offset in source
        host_id,            & ! Destination device
        device_id           & ! Source device
    )

    if (ierr /= 0) then
        ! Catastrophic error, do not even try to process it cleanly
        write(*,*) "Error in omp_target_memcpy: ", ierr
        stop
    end if

end subroutine

#endif
! =============  end OMPGPU =============

subroutine allocate_cmat_fp32(nv,nc_loc_coll,nt1,nt2)
    implicit none

    ! ----------------------
    integer, intent(in)          :: nv,nc_loc_coll,nt1,nt2

    integer :: nt_loc
    nt_loc = nt2-nt1+1

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! we keep only the GPU buffer with BIGMEM
      call allocate_cmat_gpu_fp32(nv,nc_loc_coll,nt1,nt2)
    else
#else
    if (.TRUE.) then
#endif
      allocate(cmat_fp32(nv,nv,nc_loc_coll,nt_loc))
#if defined(_OPENACC)
!$acc enter data create(cmat_fp32) if (gpu_bigmem_flag > 0)
#endif
    endif

end subroutine allocate_cmat_fp32

subroutine allocate_cmat(nv,nc_loc_coll,nt1,nt2)
    implicit none

    ! ----------------------
    integer, intent(in)          :: nv,nc_loc_coll,nt1,nt2

    integer :: nt_loc
    nt_loc = nt2-nt1+1

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! we keep only the GPU buffer with BIGMEM
      call allocate_cmat_gpu(nv,nc_loc_coll,nt1,nt2)
    else
#else
    if (.TRUE.) then
#endif
      allocate(cmat(nv,nv,nc_loc_coll,nt_loc))
#if defined(_OPENACC)
!$acc enter data create(cmat) if (gpu_bigmem_flag > 0)
#endif
    endif
end subroutine allocate_cmat

subroutine deallocate_cmat_fp32
    implicit none

    if (associated(cmat_fp32)) then
#ifdef OMPGPU
      if (gpu_bigmem_flag > 0) then
        call deallocate_cmat_gpu_fp32
      else
#else
      if (.TRUE.) then
#endif
        deallocate(cmat_fp32)
      endif
    endif

end subroutine
        
subroutine deallocate_cmat
    implicit none

    if (associated(cmat)) then
#ifdef OMPGPU
      if (gpu_bigmem_flag > 0) then
        call deallocate_cmat_gpu
      else
#else
      if (.TRUE.) then
#endif
        deallocate(cmat)
      endif
    endif

end subroutine

subroutine copy_into_cmat_fp32(amat,ic_loc,itor)
    implicit none

    ! ----------------------
    real(KIND=REAL32), target, intent(in)    :: amat(nv,nv)
    integer, intent(in) :: ic_loc,itor

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! cmat_fp32 in GPU memory only
      call copy_into_cmat_gpu_fp32(amat,ic_loc,itor)
    else
#else
    if (.TRUE.) then
#endif
      cmat_fp32(:,:,ic_loc,itor-nt1+1) = amat(:,:)
    endif

end subroutine

subroutine copy_into_cmat(amat,ic_loc,itor)
    implicit none

    ! ----------------------
    real, target, intent(in)    :: amat(nv,nv)
    integer, intent(in) :: ic_loc,itor

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! cmat in GPU memory only
      call copy_into_cmat_gpu(amat,ic_loc,itor)
    else
#else
    if (.TRUE.) then
#endif
      cmat(:,:,ic_loc,itor-nt1+1) = amat(:,:)
    endif

end subroutine

subroutine copy_from_cmat(amat,ic_loc,itor)
    implicit none

    ! ----------------------
    real, target, intent(inout)    :: amat(nv,nv)
    integer, intent(in) :: ic_loc,itor

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! cmat in GPU memory only
      call copy_from_cmat_gpu(amat,ic_loc,itor)
    else
#else
    if (.TRUE.) then
#endif
       amat(:,:) = cmat(:,:,ic_loc,itor-nt1+1)
    endif

end subroutine

subroutine copy_from_cmat_all(cmat_dest)
    implicit none

    ! ----------------------
    real, target, intent(inout) :: cmat_dest(nv,nv,nc_loc_coll,nt_loc)

#ifdef OMPGPU
    if (gpu_bigmem_flag > 0) then
      ! cmat in GPU memory only
      call copy_from_cmat_gpu_all(cmat_dest)
    else
#else
    if (.TRUE.) then
#endif
      cmat_dest(:,:,:,:) = cmat(:,:,:,:)
    endif

end subroutine

end module cgyro_globals
