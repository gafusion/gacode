!-------------------------------------------------------------------------
! tglf_interface.f90
!
! PURPOSE:
!  Provides interface description for TGLF.
!
! CALLING SEQUENCE:
!  call tglf_init(...)
!  set tglf_*_in variables
!  call tglf_run(...)
!  get tglf_*_out variables
!
! NOTES:
!  tglf_path:        working directory 
!  tglf_dump_flag:   TRUE = call dump routines, FALSE = do NOT call dump routines 
!  tglf_dump_suffix: appended to name of dump file 
!
!  tglf_dump_local and tglf_dump_global are called by tglf_run
!
! NAMING CONVENTION:
!  Interface parameter names are formed from TGLF variable names by 
!  prepending 'tglf_' and appending '_in'
!
! INITIALIZATION:
!  We initialize input arrays/vectors to zero unlike original TGLF module
!
!-------------------------------------------------------------------------

module tglf_interface

  use tglf_max_dimensions

  !  NOTE: tglf_max_dimensions defines 
  !   nsm=6
  !   maxmodes=4

  implicit none

  ! CONTROL PARAMETERS
  character (len=256)  :: tglf_path_in        = ''
  logical              :: tglf_dump_flag_in   = .false.

  ! INPUT PARAMETERS

  logical :: tglf_use_transport_model_in = .true.
  integer :: tglf_geometry_flag_in = 1
  integer :: tglf_write_wavefunction_flag_in=0

  ! Data passed to: put_signs
  real    :: tglf_sign_bt_in        = 1.0
  real    :: tglf_sign_it_in        = 1.0

  ! Data passed to: put_rare_switches
  real    :: tglf_theta_trapped_in  = 0.7
  real    :: tglf_park_in           = 1.0
  real    :: tglf_ghat_in           = 1.0
  real    :: tglf_gchat_in          = 1.0
  real    :: tglf_wd_zero_in        = 0.1
  real    :: tglf_linsker_factor_in = 0.0
  real    :: tglf_gradb_factor_in   = 0.0
  real    :: tglf_filter_in         = 0.0
  real    :: tglf_damp_psi_in       = 0.0
  real    :: tglf_damp_sig_in       = 0.0

  ! Data passed to: put_switches
  logical :: tglf_iflux_in          = .true.
  logical :: tglf_use_bper_in       = .false.
  logical :: tglf_use_bpar_in       = .false.
  logical :: tglf_use_mhd_rule_in   = .true.
  logical :: tglf_use_bisection_in  = .true.
  integer :: tglf_ibranch_in        = -1
  integer :: tglf_nmodes_in         = 2
  integer :: tglf_nbasis_max_in     = 4
  integer :: tglf_nbasis_min_in     = 1
  integer :: tglf_nxgrid_in         = 16
  integer :: tglf_nky_in            = 12

  ! Data passed to: put_model_parameters
  logical :: tglf_adiabatic_elec_in = .false.
  real    :: tglf_alpha_e_in        = 1.0
  real    :: tglf_alpha_p_in        = 1.0
  real    :: tglf_alpha_quench_in   = 0.0
  real    :: tglf_xnu_factor_in     = 1.0
  real    :: tglf_debye_factor_in   = 1.0
  real    :: tglf_etg_factor_in     = 1.25
  integer :: tglf_sat_rule_in       = 0
  integer :: tglf_kygrid_model_in   = 1
  integer :: tglf_xnu_model_in      = 2
  integer :: tglf_vpar_model_in     = 0
  integer :: tglf_vpar_shear_model_in = 1

  ! Data passed to: put_species
  integer :: tglf_ns_in             = 2
  real    :: tglf_mass_in(nsm)      = 0.0
  real    :: tglf_zs_in(nsm)        = 0.0

  ! Data passed to: put_kys
  real    :: tglf_ky_in             = 0.3

  ! Data passed to: put_gaussian_width
  real    :: tglf_width_in          = 1.65
  real    :: tglf_width_min_in      = 0.3
  integer :: tglf_nwidth_in         = 21
  logical :: tglf_find_width_in     = .true.

  ! Data passed to: put_gradients
  real    :: tglf_rlns_in(nsm)      = 0.0
  real    :: tglf_rlts_in(nsm)      = 0.0
  real    :: tglf_vpar_shear_in(nsm)= 0.0
  real    :: tglf_vexb_shear_in     = 0.0

  ! Data passed to: put_profile_shear
  real    :: tglf_vns_shear_in(nsm) = 0.0
  real    :: tglf_vts_shear_in(nsm) = 0.0

  ! Data passed to: put_averages
  real    :: tglf_taus_in(nsm)    = 0.0
  real    :: tglf_as_in(nsm)      = 0.0
  real    :: tglf_vpar_in(nsm)    = 0.0
  real    :: tglf_vexb_in         = 0.0
  real    :: tglf_betae_in        = 0.0
  real    :: tglf_xnue_in         = 0.0
  real    :: tglf_zeff_in         = 1.0
  real    :: tglf_debye_in        = 0.0

  ! Data passed to: put_eikonal
  logical :: tglf_new_eikonal_in    = .true.

  ! Data passed to: put_s_alpha_geometry
  real    :: tglf_rmin_sa_in        = 0.5
  real    :: tglf_rmaj_sa_in        = 3.0
  real    :: tglf_q_sa_in           = 2.0
  real    :: tglf_shat_sa_in        = 1.0
  real    :: tglf_alpha_sa_in       = 0.0
  real    :: tglf_xwell_sa_in       = 0.0
  real    :: tglf_theta0_sa_in      = 0.0
  integer :: tglf_b_model_sa_in     = 1
  integer :: tglf_ft_model_sa_in    = 1

  ! Data passed to: put_Miller_geometry
  real    :: tglf_rmin_loc_in       = 0.5
  real    :: tglf_rmaj_loc_in       = 3.0
  real    :: tglf_zmaj_loc_in       = 0.0
  real    :: tglf_drmindx_loc_in    = 1.0
  real    :: tglf_drmajdx_loc_in    = 0.0
  real    :: tglf_dzmajdx_loc_in    = 0.0
  real    :: tglf_kappa_loc_in      = 1.0
  real    :: tglf_s_kappa_loc_in    = 0.0
  real    :: tglf_delta_loc_in      = 0.0
  real    :: tglf_s_delta_loc_in    = 0.0
  real    :: tglf_zeta_loc_in       = 0.0
  real    :: tglf_s_zeta_loc_in     = 0.0
  real    :: tglf_q_loc_in          = 2.0
  real    :: tglf_q_prime_loc_in    = 16.0
  real    :: tglf_p_prime_loc_in    = 0.0
  real    :: tglf_kx0_loc_in        = 0.0

  ! Data passed to Fourier_geometry
  integer :: tglf_nfourier_in        = 16
  real :: tglf_q_fourier_in          = 2.0
  real :: tglf_q_prime_fourier_in    = 16.0
  real :: tglf_p_prime_fourier_in    = 0.0
  real :: tglf_fourier_in(8,0:max_fourier)=0.0

  ! Data passed to: put_ELITE_geometry
  integer  :: tglf_n_elite_in   = 700
  real     :: tglf_q_elite_in   = 2.0
  real     :: tglf_q_prime_elite_in = 16.0
  real     :: tglf_R_elite_in(max_ELITE) 
  real     :: tglf_Z_elite_in(max_ELITE) 
  real     :: tglf_Bp_elite_in(max_ELITE) 

  ! TRANSPORT OUTPUT PARAMETERS
  real :: tglf_elec_pflux_out = 0.0
  real :: tglf_elec_eflux_out = 0.0
  real :: tglf_elec_eflux_low_out = 0.0
  real :: tglf_elec_mflux_out = 0.0
  real :: tglf_elec_expwd_out = 0.0

  real, dimension(nsm-1) :: tglf_ion_pflux_out = 0.0
  real, dimension(nsm-1) :: tglf_ion_eflux_out = 0.0
  real, dimension(nsm-1) :: tglf_ion_eflux_low_out = 0.0
  real, dimension(nsm-1) :: tglf_ion_mflux_out = 0.0
  real, dimension(nsm-1) :: tglf_ion_expwd_out = 0.0
  
  ! LINEAR OUTPUT PARAMETERS
  complex :: tglf_eigenvalue_out(maxmodes)

  ! ERROR OUTPUT
  character (len=80) :: tglf_error_msg='null'
  integer :: tglf_error_flag=0

contains

  ! Dump LOCAL INTERFACE variables
  subroutine tglf_dump_local()

    implicit none

    integer :: ierr, i

    open(unit=1,file=trim(tglf_path_in)//'out.tglf.localdump',&
         status='replace',iostat=ierr)

    if (tglf_geometry_flag_in == 2 ) then
       write(1,*) 'DUMP NOT SUPPORTED FOR FOURIER GEOMETRY'
       close(1)
       return
    endif

    write(1,*) '# input.tglf generated by tglf_dump_local()'
    write(1,*) '#'
    write(1,*) '# See https://fusion.gat.com/theory/tglfinput'
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Control parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,20) 'NS',tglf_ns_in
    write(1,10) 'USE_TRANSPORT_MODEL',tglf_use_transport_model_in
    write(1,20) 'GEOMETRY_FLAG',tglf_geometry_flag_in
    write(1,10) 'USE_BPER',tglf_use_bper_in
    write(1,10) 'USE_BPAR',tglf_use_bpar_in
    write(1,10) 'USE_MHD_RULE',tglf_use_mhd_rule_in
    write(1,10) 'USE_BISECTION',tglf_use_bisection_in
    write(1,20) 'SAT_RULE',tglf_sat_rule_in
    write(1,20) 'KYGRID_MODEL',tglf_kygrid_model_in
    write(1,20) 'XNU_MODEL',tglf_xnu_model_in
    write(1,20) 'VPAR_MODEL',tglf_vpar_model_in
    write(1,20) 'VPAR_SHEAR_MODEL',tglf_vpar_shear_model_in 
    write(1,30) 'SIGN_BT',tglf_sign_bt_in
    write(1,30) 'SIGN_IT',tglf_sign_it_in
    write(1,30) 'KY',tglf_ky_in
    write(1,10) 'NEW_EIKONAL',tglf_new_eikonal_in
    write(1,30) 'VEXB',tglf_vexb_in
    write(1,30) 'VEXB_SHEAR',tglf_vexb_shear_in
    write(1,30) 'BETAE',tglf_betae_in
    write(1,30) 'XNUE',tglf_xnue_in
    write(1,30) 'ZEFF',tglf_zeff_in
    write(1,30) 'DEBYE',tglf_debye_in
    write(1,10) 'IFLUX',tglf_iflux_in
    write(1,20) 'IBRANCH',tglf_ibranch_in
    write(1,20) 'NMODES',tglf_nmodes_in
    write(1,20) 'NBASIS_MAX',tglf_nbasis_max_in
    write(1,20) 'NBASIS_MIN',tglf_nbasis_min_in
    write(1,20) 'NXGRID',tglf_nxgrid_in
    write(1,20) 'NKY',tglf_nky_in
    write(1,10) 'ADIABATIC_ELEC',tglf_adiabatic_elec_in
    write(1,30) 'ALPHA_E',tglf_alpha_e_in
    write(1,30) 'ALPHA_P',tglf_alpha_p_in
    write(1,30) 'ALPHA_QUENCH',tglf_alpha_quench_in
    write(1,30) 'XNU_FACTOR',tglf_xnu_factor_in
    write(1,30) 'DEBYE_FACTOR',tglf_debye_factor_in
    write(1,30) 'ETG_FACTOR',tglf_etg_factor_in
    write(1,20) 'WRITE_WAVEFUNCTION_FLAG',tglf_write_wavefunction_flag_in
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Species vectors:'
    write(1,*) '#---------------------------------------------------'
    do i = 1,tglf_ns_in
        write(1,40) 'ZS_',i,tglf_zs_in(i)
    enddo
    do i = 1,tglf_ns_in
        write(1,40) 'MASS_',i,tglf_mass_in(i)
    enddo
    do i = 1,tglf_ns_in
        write(1,40) 'RLNS_',i,tglf_rlns_in(i)
    enddo
    do i = 1,tglf_ns_in
        write(1,40) 'RLTS_',i,tglf_rlts_in(i)
    enddo
    do i = 1,tglf_ns_in
        write(1,40) 'TAUS_',i,tglf_taus_in(i)
    enddo
    do i = 1,tglf_ns_in
        write(1,40) 'AS_',i,tglf_as_in(i)
    enddo
    do i = 1,tglf_ns_in
        write(1,40) 'VPAR_',i,tglf_vpar_in(i)
    enddo
    do i = 1,tglf_ns_in
        write(1,40) 'VPAR_SHEAR_',i,tglf_vpar_shear_in(i)
    enddo
    do i = 1,tglf_ns_in
        write(1,40) 'VNS_SHEAR_',i,tglf_vns_shear_in(i)
    enddo
    do i = 1,tglf_ns_in
        write(1,40) 'VTS_SHEAR_',i,tglf_vts_shear_in(i)
    enddo
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Gaussian width parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'WIDTH',tglf_width_in
    write(1,30) 'WIDTH_MIN',tglf_width_min_in
    write(1,10) 'FIND_WIDTH',tglf_find_width_in
    write(1,20) 'NWIDTH',tglf_nwidth_in
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Miller geometry parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'RMIN_LOC',tglf_rmin_loc_in
    write(1,30) 'RMAJ_LOC',tglf_rmaj_loc_in
    write(1,30) 'ZMAJ_LOC',tglf_zmaj_loc_in
    write(1,30) 'DRMINDX_LOC',tglf_drmindx_loc_in
    write(1,30) 'DRMAJDX_LOC',tglf_drmajdx_loc_in
    write(1,30) 'DZMAJDX_LOC',tglf_dzmajdx_loc_in
    write(1,30) 'Q_LOC',tglf_q_loc_in
    write(1,30) 'KAPPA_LOC',tglf_kappa_loc_in
    write(1,30) 'S_KAPPA_LOC',tglf_s_kappa_loc_in
    write(1,30) 'DELTA_LOC',tglf_delta_loc_in
    write(1,30) 'S_DELTA_LOC',tglf_s_delta_loc_in
    write(1,30) 'ZETA_LOC',tglf_zeta_loc_in
    write(1,30) 'S_ZETA_LOC',tglf_s_zeta_loc_in
    write(1,30) 'P_PRIME_LOC',tglf_p_prime_loc_in
    write(1,30) 'Q_PRIME_LOC',tglf_q_prime_loc_in
    write(1,30) 'KX0_LOC',tglf_kx0_loc_in
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# s-alpha geometry parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'RMIN_SA',tglf_rmin_sa_in
    write(1,30) 'RMAJ_SA',tglf_rmaj_sa_in
    write(1,30) 'Q_SA',tglf_q_sa_in
    write(1,30) 'SHAT_SA',tglf_shat_sa_in
    write(1,30) 'ALPHA_SA',tglf_alpha_sa_in
    write(1,30) 'XWELL_SA',tglf_xwell_sa_in
    write(1,30) 'THETA0_SA',tglf_theta0_sa_in
    write(1,20) 'B_MODEL_SA',tglf_b_model_sa_in
    write(1,20) 'FT_MODEL_SA',tglf_ft_model_sa_in
    write(1,*) ' '
    write(1,*) '#---------------------------------------------------'
    write(1,*) '# Expert parameters:'
    write(1,*) '#---------------------------------------------------'
    write(1,30) 'DAMP_PSI',tglf_damp_psi_in
    write(1,30) 'DAMP_SIG',tglf_damp_sig_in
    write(1,30) 'PARK',tglf_park_in
    write(1,30) 'GHAT',tglf_ghat_in
    write(1,30) 'GCHAT',tglf_gchat_in
    write(1,30) 'WD_ZERO',tglf_wd_zero_in
    write(1,30) 'LINSKER_FACTOR',tglf_linsker_factor_in
    write(1,30) 'GRADB_FACTOR',tglf_gradB_factor_in
    write(1,30) 'FILTER',tglf_filter_in
    write(1,30) 'THETA_TRAPPED',tglf_theta_trapped_in

    close(1)

10  format(t2,a,'=',l1)
20  format(t2,a,'=',i3)
30  format(t2,a,'=',1pe12.5)
40  format(t2,a,i1,'=',1pe12.5)

  end subroutine tglf_dump_local

  ! Dump GLOBAL MODULE variables
  subroutine tglf_dump_global()

    use tglf_global

    implicit none

    integer :: ierr

    open(unit=1,file=trim(tglf_path_in)//'out.tglf.globaldump',&
         status='replace',iostat=ierr)

    write(1,*) 'new_eikonal = ', new_eikonal_in
    write(1,*) 'find_width = ', find_width_in
    write(1,*) 'nwidth = ', nwidth_in
    write(1,*) 'width = ', width_in
    write(1,*) 'width_min = ', width_min_in
    write(1,*) 'ky = ', ky_in
    write(1,*) 'nky = ', nky_in
    write(1,*) 'nmodes = ', nmodes_in
    write(1,*) 'ns = ', ns_in
    write(1,*) 'mass = ', mass_in(:)
    write(1,*) 'zs = ', zs_in(:)
    write(1,*) 'iflux = ', iflux_in
    write(1,*) 'use_bper = ', use_bper_in
    write(1,*) 'use_bpar = ', use_bpar_in
    write(1,*) 'use_mhd_rule = ', use_mhd_rule_in
    write(1,*) 'use_bisection = ', use_bisection_in
    write(1,*) 'ibranch = ', ibranch_in
    write(1,*) 'nbasis_max = ', nbasis_max_in
    write(1,*) 'nbasis_min = ', nbasis_min_in
    write(1,*) 'nxgrid = ', nxgrid_in
    write(1,*) 'adiabatic_elec = ', adiabatic_elec_in
    write(1,*) 'damp_psi = ', damp_psi_in
    write(1,*) 'damp_sig = ', damp_sig_in
    write(1,*) 'park = ', park_in
    write(1,*) 'ghat = ', ghat_in
    write(1,*) 'gchat = ', gchat_in
    write(1,*) 'wd_zero = ', wd_zero_in
    write(1,*) 'linsker_factor = ', linsker_factor_in
    write(1,*) 'gradB_factor = ', gradB_factor_in
    write(1,*) 'filter = ', filter_in
    write(1,*) 'sat_rule = ', sat_rule_in
    write(1,*) 'alpha_quench = ', alpha_quench_in
    write(1,*) 'alpha_e = ', alpha_e_in
    write(1,*) 'alpha_p = ', alpha_p_in
    write(1,*) 'theta_trapped = ', theta_trapped_in
    write(1,*) 'xnu_factor = ', xnu_factor_in
    write(1,*) 'debye_factor = ', debye_factor_in
    write(1,*) 'etg_factor = ', etg_factor_in
    write(1,*) 'kygrid_model = ', kygrid_model_in
    write(1,*) 'xnu_model = ', xnu_model_in
    write(1,*) 'vpar_model = ', vpar_model_in
    write(1,*) 'vpar_shear_model = ', vpar_shear_model_in 
    write(1,*) 'sign_Bt = ', sign_Bt_in
    write(1,*) 'sign_It = ', sign_It_in
    write(1,*) 'rlns = ', rlns_in(:)
    write(1,*) 'rlts = ', rlts_in(:)
    write(1,*) 'vpar_shear = ', vpar_shear_in(:)
    write(1,*) 'vexb_shear = ', vexb_shear_in
    write(1,*) 'vns_shear = ', vns_shear_in(:)
    write(1,*) 'vts_shear = ', vts_shear_in(:)
    write(1,*) 'as = ', as_in(:)
    write(1,*) 'taus = ', taus_in(:)
    write(1,*) 'vpar = ', vpar_in(:)
    write(1,*) 'vexb = ', vexb_in
    write(1,*) 'betae = ', betae_in
    write(1,*) 'xnuei = ', xnue_in
    write(1,*) 'zeff = ', zeff_in
    write(1,*) 'debye = ', debye_in

    if (tglf_geometry_flag_in == 1 ) then

       write(1,*) 'rmin_loc = ', rmin_loc
       write(1,*) 'rmaj_loc = ', rmaj_loc
       write(1,*) 'zmaj_loc = ', zmaj_loc
       write(1,*) 'drmindx_loc = ', drmindx_loc
       write(1,*) 'drmajdx_loc = ', drmajdx_loc
       write(1,*) 'dzmajdx_loc = ', dzmajdx_loc
       write(1,*) 'q_loc = ', q_loc
       write(1,*) 'kappa_loc = ', kappa_loc
       write(1,*) 's_kappa_loc = ', s_kappa_loc
       write(1,*) 'delta_loc = ', delta_loc
       write(1,*) 's_delta_loc = ', s_delta_loc
       write(1,*) 'zeta_loc = ', zeta_loc
       write(1,*) 's_zeta_loc = ', s_zeta_loc
       write(1,*) 'p_prime_loc = ', p_prime_loc
       write(1,*) 'q_prime_loc = ', q_prime_loc
       write(1,*) 'kx0_loc =', kx0_loc

    elseif (tglf_geometry_flag_in == 2 ) then

       write(1,*) 'q_fourier_in = ', q_fourier_in
       write(1,*) 'q_prime_fourier_in = ', q_prime_fourier_in
       write(1,*) 'p_prime_fourier_in = ', p_prime_fourier_in
       write(1,*) 'nfourier_in = ', nfourier_in
       write(1,*) 'fourier_in = ', fourier_in(:,:)

    elseif (tglf_geometry_flag_in == 0 ) then

       write(1,*) 'rmin_sa = ', rmin_sa
       write(1,*) 'rmaj_sa = ', rmaj_sa
       write(1,*) 'q_sa = ', q_sa
       write(1,*) 'shat_sa = ', shat_sa
       write(1,*) 'alpha_sa = ', alpha_sa
       write(1,*) 'xwell_sa = ', xwell_sa
       write(1,*) 'theta0_sa = ', theta0_sa
       write(1,*) 'b_model_sa = ', b_model_sa
       write(1,*) 'ft_model_sa = ', ft_model_sa

    endif

    close(1)

  end subroutine tglf_dump_global

end module tglf_interface
