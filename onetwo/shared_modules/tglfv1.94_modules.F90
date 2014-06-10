      MODULE tglf_max_dimensionsv194
!
! global dimensions for shared arrays
!
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: nb=32
      INTEGER, PARAMETER :: nxm=4*nb-1
      INTEGER, PARAMETER :: nsm=6, nt0=40
      INTEGER, PARAMETER :: neq = 15*nsm,iar=neq*nb
      INTEGER, PARAMETER :: nkym=50
      INTEGER, PARAMETER :: maxmodes=4
      INTEGER, PARAMETER :: max_ELITE=700
      INTEGER, PARAMETER :: max_fourier = 16
      INTEGER, PARAMETER :: ms = 128  ! ms needs to be divisible by 8
      INTEGER, PARAMETER :: max_plot =6*ms/8+1
!
      END MODULE tglf_max_dimensionsv194
!
      MODULE tglf_dimensionsv194
      USE tglf_max_dimensionsv194
!
! dimensions determined by inputs
! 
      IMPLICIT NONE
      SAVE
!
      INTEGER nx,nbasis,ns0,ns
!
      END MODULE tglf_dimensionsv194
!
      MODULE tglf_globalv194
      USE tglf_max_dimensionsv194
!
!  global controls for the tglf driver routine
!
      IMPLICIT NONE
      SAVE
! global constants
      REAL :: pi
      REAL :: pi_2
      REAL :: sqrt_pi
      REAL :: sqrt_two
      COMPLEX :: xi=(0.0,1.0)
! internal flow control switches
      LOGICAL :: new_eikonal_in=.TRUE.
      LOGICAL :: new_start=.TRUE.
      LOGICAL :: new_matrix=.TRUE.
      LOGICAL :: new_geometry=.TRUE.
      LOGICAL :: new_width=.TRUE.
      LOGICAL :: new_kyspectrum=.TRUE.
      LOGICAL :: eikonal_unsaved=.TRUE.
      INTEGER :: igeo=1
      LOGICAL :: use_default_species=.TRUE.
      INTEGER,DIMENSION(8) :: trace_path=0
! Input global variables
      LOGICAL :: find_width_in=.TRUE.
      INTEGER :: nwidth_in=21
      REAL :: width_in=1.65
      REAL :: width_min_in=0.3
      REAL :: ky_in=0.3
      INTEGER :: nky_in=12
      INTEGER :: nmodes_in=2
! Input species 
      INTEGER :: ns_in=2
      REAL,DIMENSION(nsm) :: mass_in
      REAL,DIMENSION(nsm) :: zs_in
! input switches
      LOGICAL :: iflux_in=.TRUE.
      LOGICAL :: use_bper_in=.FALSE.
      LOGICAL :: use_bpar_in=.FALSE.
      LOGICAL :: use_mhd_rule_in=.TRUE.
      LOGICAL :: use_bisection_in=.TRUE.
      INTEGER :: ibranch_in=-1
      INTEGER :: nbasis_max_in=4
      INTEGER :: nbasis_min_in=1
      INTEGER :: nxgrid_in=16
      REAL :: damp_psi_in = 0.0
      REAL :: damp_sig_in = 0.0
! input rare switches
      REAL :: alpha_quench_in=0.0
      REAL :: park_in=1.0
      REAL :: ghat_in=1.0
      REAL :: gchat_in=1.0
      REAL :: wd_zero_in=0.1
      REAL :: Linsker_factor_in=0.0
      REAL :: gradB_factor_in=0.0
      REAL :: filter_in=2.0
      INTEGER :: sat_rule_in=0
      REAL :: alpha_kx_e_in=0.0
      REAL :: alpha_kx_p_in=0.0
      REAL :: alpha_kx_n_in=0.0
      REAL :: alpha_kx_t_in=0.0
! Input model paramaters
      LOGICAL :: adiabatic_elec_in=.FALSE.
      REAL :: alpha_p_in=1.0
      REAL :: alpha_e_in=1.0
      REAL :: alpha_n_in =0.0
      REAL :: alpha_t_in =0.0
      REAL :: theta_trapped_in=0.7
      REAL :: xnu_factor_in=1.0
      REAL :: debye_factor_in=1.0
      REAL :: etg_factor_in=1.25
      INTEGER :: kygrid_model_in=1
      INTEGER :: xnu_model_in=2
      INTEGER :: vpar_model_in=0
      INTEGER :: vpar_shear_model_in=0
! Input signs
      REAL :: sign_Bt_in = 1.0
      REAL :: sign_It_in = 1.0
! Input field gradients
      REAL,DIMENSION(nsm) :: rlns_in=1.0
      REAL,DIMENSION(nsm) :: rlts_in=3.0
      REAL,DIMENSION(nsm) :: vpar_shear_in=0.0
      REAL :: vexb_shear_in=0.0
! Input profile shear
      REAL,DIMENSION(nsm) :: vns_shear_in=0.0
      REAL,DIMENSION(nsm) :: vts_shear_in=0.0
! Input field averages
      REAL,DIMENSION(nsm) :: as_in
      REAL,DIMENSION(nsm) :: taus_in
      REAL,DIMENSION(nsm) :: vpar_in=0.0
      REAL :: vexb_in = 0.0
      REAL :: betae_in=0.0
      REAL :: xnue_in=0.0
      REAL :: zeff_in=1.0
      REAL :: debye_in=0.0
! Shifted circle (s-alpha) inputs
      REAL :: rmin_sa=0.5
      REAL :: rmaj_sa=3.0
      REAL :: q_sa=2.0
      REAL :: shat_sa=1.0
      REAL :: alpha_sa=0.0
      REAL :: xwell_sa=0.0
      REAL :: theta0_sa=0.0
! Shifted circle flags
      INTEGER :: b_model_sa=1
      INTEGER :: ft_model_sa=1
! Miller inputs
      REAL :: rmin_loc=0.5
      REAL :: rmaj_loc=3.0
      REAL :: zmaj_loc=0.0
      REAL :: q_loc=2.0
!      REAL :: shat_loc=1.0
      REAL :: dlnpdr_loc=0.0
      REAL :: drmindx_loc=1.0
      REAL :: drmajdx_loc=0.0
      REAL :: dzmajdx_loc=0.0
      REAL :: kappa_loc=1.0
      REAL :: s_kappa_loc=0.0
      REAL :: delta_loc=0.0
      REAL :: s_delta_loc=0.0
      REAL :: zeta_loc=0.0
      REAL :: s_zeta_loc=0.0
      REAL :: p_prime_loc=0.0 
      REAL :: q_prime_loc=16.0
      REAL :: kx0_loc = 0.0
! Fourier geometry inputs
      INTEGER :: nfourier_in        = 16
      REAL :: q_fourier_in          = 2.0
      REAL :: q_prime_fourier_in    = 16.0
      REAL :: p_prime_fourier_in    = 0.0
      REAL :: fourier_in(8,0:max_fourier) = 0.0
! ELITE geometry inputs
      INTEGER :: n_ELITE
      REAL :: R_ELITE(0:max_ELITE)
      REAL :: Z_ELITE(0:max_ELITE)
      REAL :: Bp_ELITE(0:max_ELITE)
      REAL :: q_ELITE
      REAL :: q_prime_ELITE
      REAL :: p_prime_ELITE
! global variables
      REAL :: ft=0.5
      REAL :: ft_min=0.01
      REAL :: ky=0.3
      REAL :: R_unit=3.0
      REAL :: q_unit=2.0
      REAL :: B_unit=1.0
      REAL :: Rmaj_input = 3.0
      REAL :: q_input = 2.0
      REAL :: rmin_input=0.5
      REAL,DIMENSION(maxmodes) :: gamma_reference_kx0=0.0
      REAL :: pol=1.0
      REAL :: U0=0.0
      REAL :: kx0=0.0
      REAL :: kx0_e=0.0
      REAL :: kx0_p = 0.0
      REAL :: midplane_shear=1.0,kx_shear=0.0
      REAL :: kx0_factor=1.0
! output
      COMPLEX,DIMENSION(maxmodes,3,nb) :: field_weight_out=0.0
      COMPLEX,DIMENSION(maxmodes,3,max_plot) :: plot_field_out=0.0
      REAL,DIMENSION(max_plot) :: plot_angle_out=0.0
      REAL,DIMENSION(maxmodes,nsm,3) :: particle_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm,3) :: energy_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm,3) :: stress_par_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm,3) :: stress_tor_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm,3) :: exchange_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: N_QL_out=0.0,T_QL_out=0.0
      REAL,DIMENSION(maxmodes,nsm) :: n_bar_out=0.0,t_bar_out=0.0
      REAL,DIMENSION(nsm,3) :: particle_flux_out=0.0,energy_flux_out=0.0
      REAL,DIMENSION(nsm,3) :: exchange_out=0.0
      REAL,DIMENSION(nsm,3) :: stress_par_out=0.0,stress_tor_out=0.0
      REAL,DIMENSION(maxmodes) :: gamma_out=0.0,freq_out=0.0
      REAL,DIMENSION(maxmodes) :: wd_bar_out=0.0,phi_QL_out=0.0
      REAL,DIMENSION(maxmodes) :: phi_bar_out=0.0,v_bar_out=0.0
      REAL,DIMENSION(maxmodes) :: b0_bar_out=0.0,ne_te_phase_out=0.0
      REAL,DIMENSION(maxmodes) :: kx_bar_out=0.0,kpar_bar_out=0.0
      REAL,DIMENSION(nsm) :: n_bar_sum_out=0.0,t_bar_sum_out=0.0
      REAL,DIMENSION(nsm) :: q_low_out=0.0
      REAL,DIMENSION(2,nkym,maxmodes) :: field_spectrum_out=0.0
      REAL,DIMENSION(2,nsm,nkym,maxmodes) :: intensity_spectrum_out=0.0
      REAL,DIMENSION(5,nsm,3,nkym,maxmodes) :: flux_spectrum_out=0.0
      REAL,DIMENSION(2,nkym,maxmodes) :: eigenvalue_spectrum_out=0.0
      REAL :: phi_bar_sum_out=0.0
      REAL :: v_bar_sum_out=0.0
      REAL :: gamma_nb_min_out=0.0
      REAL :: B2_ave_out=1.0
      REAL :: R2_ave_out=1.0
      REAL :: a_pol_out=1.0
      REAL :: a_tor_out=1.0
      REAL :: Bp0_out = 1.0
      REAL :: RBt_ave_out=1.0
      REAL :: DM_out = 0.25
      REAL :: DR_out = 0.0
      INTEGER :: nmodes_out
      INTEGER :: nfields_out
!
      END MODULE tglf_globalv194
!------------------------------------------------- 
      MODULE tglf_closurev194
!
!  tglf closure coefficients
!
      IMPLICIT NONE
      SAVE
! ft = 1 coefficients
      REAL :: v1_r,v2_r,v3_r,v4_r,v5_r
      REAL :: v1_i,v2_i,v3_i,v4_i,v5_i
      REAL :: v6_r,v7_r,v8_r,v9_r,v10_r
      REAL :: v6_i,v7_i,v8_i,v9_i,v10_i
      REAL :: vb1_r,vb2_r,vb3_r,vb4_r,vb5_r
      REAL :: vb1_i,vb2_i,vb3_i,vb4_i,vb5_i
      REAL :: vb6_r,vb7_r,vb8_r,vb9_r,vb10_r
      REAL :: vb6_i,vb7_i,vb8_i,vb9_i,vb10_i
! ft < 1 coefficients
      REAL :: u1_r,u2_r,u3_r,u4_r,u5_r
      REAL :: u1_i,u2_i,u3_i,u4_i,u5_i
      REAL :: u6_r,u7_r,u8_r,u9_r,u10_r
      REAL :: u6_i,u7_i,u8_i,u9_i,u10_i
      REAL :: ub1_r,ub2_r,ub3_r,ub4_r,ub5_r
      REAL :: ub1_i,ub2_i,ub3_i,ub4_i,ub5_i
      REAL :: ub6_r,ub7_r,ub8_r,ub9_r,ub10_r
      REAL :: ub6_i,ub7_i,ub8_i,ub9_i,ub10_i
! parallel coefficients
      REAL :: dpar_HP,dper_HP,bpar_HP,bper_HP
      REAL :: bper_DH,dper_DH
      COMPLEX :: dhr13,dgr13
      REAL :: d1,d3,d33
      REAL :: b1,b3,b33
!
      END MODULE tglf_closurev194
!-------------------------------------------------
      MODULE tglf_hermitev194
! hermite basis functions and x-grid
      USE tglf_max_dimensionsv194
      IMPLICIT NONE
      SAVE
!
      REAL :: x(nxm),wx(nxm),h(nb,nxm)
!
      END MODULE tglf_hermitev194  
!
      MODULE tglf_speciesv194
! species parameters
      USE tglf_max_dimensionsv194
      IMPLICIT NONE
      SAVE
!
      REAL :: ei_exch(nsm,nsm),resist(nsm,nsm)
      REAL :: zs(nsm),mass(nsm),vs(nsm)
      REAL :: rlts(nsm),rlns(nsm),vpar_shear_s(nsm)
      REAL :: as(nsm),taus(nsm),vpar_s(nsm)
      REAL :: vexb_shear_s
! 
      END MODULE tglf_speciesv194  
!-------------------------------------------------
      MODULE tglf_kyspectrumv194
!
! ky spectrum for computing total fluxes and intensities
!
      USE tglf_max_dimensionsv194
      IMPLICIT NONE
      SAVE
!
      INTEGER nky
!
      REAL :: ky_spectrum(nkym)
      REAL :: dky_spectrum(nkym)
!
      END MODULE tglf_kyspectrumv194
!
!-------------------------------------------------
      MODULE tglf_eigenv194
!
! eigenvalues, eigenvectors and fluxes
!
      USE tglf_max_dimensionsv194
      IMPLICIT NONE
      SAVE
!
      INTEGER matz,iur,nroot
      REAL :: fv1(iar),fv2(iar),fv3(iar)
      REAL :: rr(iar), ri(iar)
      REAL :: ar(iar,iar), ai(iar,iar)
      REAL :: vr(iar,iar), vi(iar,iar)
      COMPLEX :: v(iar),eigenvalue
      COMPLEX :: amat(iar,iar),bmat(iar,iar)
      COMPLEX :: alpha(iar),beta(iar)
!
      END MODULE tglf_eigenv194
!
!-------------------------------------------------
!
      MODULE tglf_xgridv194
!
! functions on the x-grid
!
      USE tglf_max_dimensionsv194
      IMPLICIT NONE
      SAVE
!
      REAL :: hxn(nsm,nxm),hxp1(nsm,nxm),hxp3(nsm,nxm)
      REAL :: hxr11(nsm,nxm),hxr13(nsm,nxm),hxr33(nsm,nxm)
      REAL :: hxw113(nsm,nxm),hxw133(nsm,nxm),hxw333(nsm,nxm)
      REAL :: gxn(nsm,nxm),gxp1(nsm,nxm),gxp3(nsm,nxm)
      REAL :: gxr11(nsm,nxm),gxr13(nsm,nxm),gxr33(nsm,nxm)
      REAL :: gxw113(nsm,nxm),gxw133(nsm,nxm),gxw333(nsm,nxm)
      INTEGER :: mask_save(nkym)
      REAL :: wdx(nxm),b0x(nxm),kxx(nxm)
      REAL :: cx_tor_par(nxm),cx_tor_per(nxm)
      REAL :: cx_par_par(nxm)
      REAL :: p0x(nxm),Bx(nxm)
      REAL :: gamma_nb_min_save(nkym)
      REAL :: width_save(nkym),ft_save(nkym)
      REAL :: R_unit_save(nkym),q_unit_save(nkym)
      REAL :: wdx_save(nkym,nxm),b0x_save(nkym,nxm)
      REAL :: cx_par_par_save(nkym,nxm)
      REAL :: cx_tor_par_save(nkym,nxm)
      REAL :: cx_tor_per_save(nkym,nxm)
      REAL :: kxx_save(nkym,nxm)
!
      END MODULE tglf_xgridv194
!
!-------------------------------------------------
!
      MODULE tglf_sgridv194
!
! functions on the s-grid
!
      USE tglf_max_dimensionsv194
      IMPLICIT NONE
      SAVE
!
!---------------------------------------------------------------
! s_grid.m
! 
! PURPOSE:
!  include file which contains the data for mercier_luc.f :
! INPUT
!  R, Z and Bp on the s-grid with ms equally space intervals of length ds
!  as defined by the polidal co-ordinate in the mercier-luc system. 
!  B_unit = B0 (r/rho) dr/drho where B0 = magnetic field on axis
!  Rmaj_s = major radius at magnetic axis
!  rmin_s = minor radius of flux surface 
!  q_s = local flux surface safety factor 
!  q_prime_s = dq/dpsi
!  p_prime_s = dp/dpsi 
!
! OUTPUT
!  costheta_geo, sintheta_geo, costheta_p_geo, pk_geo, epsl_geo, qrat_geo
!  kyoky_geo, b_geo
!
! 24 June 05: gms
!  this version for TGLF separated miller and mercier-luc components
!  in GKS and GYRO versions
!---------------------------------------------------------------
! INPUT
      REAL :: R(0:ms), Z(0:ms), Bp(0:ms)
      REAL :: ds, Ls
      REAL :: Rmaj_s, Zmaj_s,rmin_s, q_s
      REAL :: p_prime_s, q_prime_s
      REAL :: p_prime_zero_s
! OUTPUT
      REAL :: costheta_geo(0:ms),sintheta_geo(0:ms)
      REAL :: costheta_p_geo(0:ms),s_p(0:ms)
      REAL :: pk_geo(0:ms),epsl_geo(0:ms),qrat_geo(0:ms)
      REAL :: kxoky_geo(0:ms), b_geo(0:ms),t_s(0:ms)
      REAL :: S_prime(0:ms),kx_factor(0:ms),y(0:ms)
      REAL :: f,ff_prime
!
      END MODULE tglf_sgridv194 
!
!-------------------------------------------------
!
      MODULE tglf_coeffv194
!
! store the hermite basis matrix coefficients
!
      USE tglf_max_dimensionsv194
      IMPLICIT NONE
      SAVE
! ave_h
      REAL :: hn,hp1,hp3,hr11,hr13,hr33
      REAl :: hw113,hw133,hw333
      REAL :: hu1,hu3,hu33,ht1,ht3
      REAL :: hu3ht1,hu3ht3,hu33ht1,hu33ht3
      REAL :: hb1,hb3,hb33,hb1ht1,hb3ht3,hb33ht1
      REAL :: hd1,hd3,hd33,hd1hu1,hd3hu3,hd33hu1
      REAL :: hv1r,hv2r,hv3r,hv4r,hv5r
      REAL :: hv6r,hv7r,hv8r,hv9r,hv10r
      REAL :: hv1rht1,hv2rht3,hv3rht1,hv4rht3
      REAL :: hv6rhu1,hv7rhu3,hv9rhu1,hv10rhu3
      REAL :: hv1i,hv2i,hv3i,hv4i,hv5i
      REAL :: hv6i,hv7i,hv8i,hv9i,hv10i
      REAL :: hv1iht1,hv2iht3,hv3iht1,hv4iht3
      REAL :: hv6ihu1,hv7ihu3,hv9ihu1,hv10ihu3
      REAL :: gradhp1,gradhr11,gradhr13
      REAL :: gradhp1p1,gradhr11p1,gradhr13p1
      REAl :: grad_hu1,grad_hu3
      REAL :: wdhu3ht1,wdhu3ht3,wdhu33ht1,wdhu33ht3
      REAL :: modwdhu3,modwdhu33
      REAL :: modwdhu3ht1,modwdhu3ht3,modwdhu33ht1,modwdhu33ht3
      REAL :: ave_hn(nsm,nb,nb)
      REAL :: ave_hp1(nsm,nb,nb),ave_hp3(nsm,nb,nb)
      REAL :: ave_hr11(nsm,nb,nb),ave_hr13(nsm,nb,nb),ave_hr33(nsm,nb,nb)
      REAL :: ave_hw113(nsm,nb,nb),ave_hw133(nsm,nb,nb),ave_hw333(nsm,nb,nb)
      REAL :: ave_ht1(nsm,nb,nb),ave_ht3(nsm,nb,nb)
      REAL :: ave_hu1(nsm,nb,nb),ave_hu3(nsm,nb,nb),ave_hu33(nsm,nb,nb)
      REAL :: ave_hu3ht1(nsm,nb,nb),ave_hu3ht3(nsm,nb,nb)
      REAL :: ave_hu33ht1(nsm,nb,nb),ave_hu33ht3(nsm,nb,nb)
      REAL :: ave_hninv(nsm,nb,nb),ave_hp1inv(nsm,nb,nb)
      REAL :: ave_hp3inv(nsm,nb,nb)
      REAL :: ave_gradhp1(nsm,nb,nb),ave_gradhr11(nsm,nb,nb)
      REAL :: ave_gradhr13(nsm,nb,nb)
      REAL :: ave_gradhp1p1(nsm,nb,nb),ave_gradhr11p1(nsm,nb,nb)
      REAL :: ave_gradhr13p1(nsm,nb,nb)
      REAL :: ave_gradhu1(nsm,nb,nb),ave_gradhu3(nsm,nb,nb)
      REAL :: ave_hnp0(nsm,nb,nb)
      REAL :: ave_hp1p0(nsm,nb,nb),ave_hp3p0(nsm,nb,nb)
      REAL :: ave_hr11p0(nsm,nb,nb),ave_hr13p0(nsm,nb,nb)
      REAL :: ave_hr33p0(nsm,nb,nb)
      REAL :: ave_hnb0(nsm,nb,nb)
      REAL :: ave_hp1b0(nsm,nb,nb),ave_hp3b0(nsm,nb,nb)
      REAL :: ave_hr11b0(nsm,nb,nb),ave_hr13b0(nsm,nb,nb)
      REAL :: ave_hr33b0(nsm,nb,nb)
      REAL :: ave_hw113b0(nsm,nb,nb),ave_hw133b0(nsm,nb,nb)
      REAl :: ave_hw333b0(nsm,nb,nb)
      REAL :: ave_hnbp(nsm,nb,nb)
      REAL :: ave_hp1bp(nsm,nb,nb),ave_hp3bp(nsm,nb,nb)
      REAL :: ave_hr11bp(nsm,nb,nb),ave_hr13bp(nsm,nb,nb)
      REAL :: ave_hr33bp(nsm,nb,nb)
      REAL :: ave_hw113bp(nsm,nb,nb),ave_hw133bp(nsm,nb,nb)
      REAl :: ave_hw333bp(nsm,nb,nb)
      REAL :: ave_hnp0b0(nsm,nb,nb)
      REAL :: ave_gradhp1p0(nsm,nb,nb),ave_gradhr11p0(nsm,nb,nb)
      REAL :: ave_gradhr13p0(nsm,nb,nb)
      REAL :: ave_wdhu3ht1(nsm,nb,nb),ave_wdhu3ht3(nsm,nb,nb)
      REAL :: ave_wdhu33ht1(nsm,nb,nb),ave_wdhu33ht3(nsm,nb,nb)
      REAL :: ave_modwdhu3(nsm,nb,nb),ave_modwdhu33(nsm,nb,nb)
      REAL :: ave_modwdhu3ht1(nsm,nb,nb),ave_modwdhu3ht3(nsm,nb,nb)
      REAL :: ave_modwdhu33ht1(nsm,nb,nb),ave_modwdhu33ht3(nsm,nb,nb)
      REAL :: ave_c_tor_par_hp1p0(nsm,nb,nb)
      REAL :: ave_c_tor_par_hr11p0(nsm,nb,nb)
      REAL :: ave_c_tor_par_hr13p0(nsm,nb,nb)
! ave_g
      REAL :: gn,gp1,gp3,gr11,gr13,gr33
      REAL :: gw113,gw133,gw333
      REAL :: gu1,gu3,gu33,gt1,gt3
      REAL :: gu3gt1,gu3gt3,gu33gt1,gu33gt3
      REAL :: gu1r,gu2r,gu3r,gu4r,gu5r
      REAL :: gu6r,gu7r,gu8r,gu9r,gu10r
      REAL :: gu1i,gu2i,gu3i,gu4i,gu5i
      REAL :: gu6i,gu7i,gu8i,gu9i,gu10i
      REAL :: gb1,gd1,gb3,gd3,gb33,gd33
      REAL :: gb1gt1,gd1gu1,gb3gt3,gd3gu3
      REAL :: gb33gt1,gd33gu1
      REAL :: gu1rgt1,gu1igt1,gu2rgt3,gu2igt3
      REAL :: gu3rgt1,gu3igt1,gu4rgt3,gu4igt3
      REAL :: gu6rgu1,gu7rgu3,gu9rgu1,gu10rgu3
      REAL :: gu6igu1,gu7igu3,gu9igu1,gu10igu3
      REAL :: gradgp1,gradgr11,gradgr13
      REAL :: grad_gu1,grad_gu3
      REAL :: gradgp1p1,gradgr11p1,gradgr13p1
      REAL :: wdgu3gt1,wdgu3gt3,wdgu33gt1,wdgu33gt3
      REAL :: modwdgu3,modwdgu33
      REAL :: modwdgu3gt1,modwdgu3gt3,modwdgu33gt1,modwdgu33gt3
      REAL :: ave_gn(nsm,nb,nb)
      REAL :: ave_gp1(nsm,nb,nb),ave_gp3(nsm,nb,nb)
      REAL :: ave_gr11(nsm,nb,nb),ave_gr13(nsm,nb,nb),ave_gr33(nsm,nb,nb)
      REAL :: ave_gw113(nsm,nb,nb),ave_gw133(nsm,nb,nb),ave_gw333(nsm,nb,nb)
      REAL :: ave_gt1(nsm,nb,nb),ave_gt3(nsm,nb,nb)
      REAL :: ave_gu1(nsm,nb,nb),ave_gu3(nsm,nb,nb),ave_gu33(nsm,nb,nb)
      REAL :: ave_gu3gt1(nsm,nb,nb),ave_gu3gt3(nsm,nb,nb)
      REAL :: ave_gu33gt1(nsm,nb,nb),ave_gu33gt3(nsm,nb,nb)
      REAL :: ave_gninv(nsm,nb,nb),ave_gp1inv(nsm,nb,nb)
      REAL :: ave_gp3inv(nsm,nb,nb)
      REAL :: ave_gradgp1(nsm,nb,nb),ave_gradgr11(nsm,nb,nb)
      REAL :: ave_gradgr13(nsm,nb,nb)
      REAL :: ave_gradgp1p1(nsm,nb,nb),ave_gradgr11p1(nsm,nb,nb)
      REAL :: ave_gradgr13p1(nsm,nb,nb)
      REAL :: ave_gradgu1(nsm,nb,nb),ave_gradgu3(nsm,nb,nb)
      REAL :: ave_gnp0(nsm,nb,nb)
      REAL :: ave_gp1p0(nsm,nb,nb),ave_gp3p0(nsm,nb,nb)
      REAL :: ave_gr11p0(nsm,nb,nb),ave_gr13p0(nsm,nb,nb)
      REAL :: ave_gr33p0(nsm,nb,nb)
      REAL :: ave_gnb0(nsm,nb,nb)
      REAL :: ave_gp1b0(nsm,nb,nb),ave_gp3b0(nsm,nb,nb)
      REAL :: ave_gr11b0(nsm,nb,nb),ave_gr13b0(nsm,nb,nb)
      REAL :: ave_gr33b0(nsm,nb,nb)
      REAL :: ave_gw113b0(nsm,nb,nb),ave_gw133b0(nsm,nb,nb)
      REAL :: ave_gw333b0(nsm,nb,nb)
      REAL :: ave_gnbp(nsm,nb,nb)
      REAL :: ave_gp1bp(nsm,nb,nb),ave_gp3bp(nsm,nb,nb)
      REAL :: ave_gr11bp(nsm,nb,nb),ave_gr13bp(nsm,nb,nb)
      REAL :: ave_gr33bp(nsm,nb,nb)
      REAL :: ave_gw113bp(nsm,nb,nb),ave_gw133bp(nsm,nb,nb)
      REAL :: ave_gw333bp(nsm,nb,nb)
      REAL :: ave_gradgp1p0(nsm,nb,nb),ave_gradgr11p0(nsm,nb,nb)
      REAL :: ave_gradgr13p0(nsm,nb,nb)
      REAL :: ave_wdgu3gt1(nsm,nb,nb),ave_wdgu3gt3(nsm,nb,nb)
      REAL :: ave_wdgu33gt1(nsm,nb,nb),ave_wdgu33gt3(nsm,nb,nb)
      REAL :: ave_modwdgu3(nsm,nb,nb),ave_modwdgu33(nsm,nb,nb)
      REAL :: ave_modwdgu3gt1(nsm,nb,nb),ave_modwdgu3gt3(nsm,nb,nb)
      REAL :: ave_modwdgu33gt1(nsm,nb,nb),ave_modwdgu33gt3(nsm,nb,nb)
      REAL :: ave_c_tor_par_gp1p0(nsm,nb,nb)
      REAL :: ave_c_tor_par_gr11p0(nsm,nb,nb)
      REAL :: ave_c_tor_par_gr13p0(nsm,nb,nb)
! ave_wd_h
      REAL :: wdhp1p0,wdhr11p0,wdhr13p0
      REAL :: wdhp1b0,wdhr11b0,wdhr13b0
      REAL :: wdhp1bp,wdhr11bp,wdhr13bp
      REAL :: wdhu1,wdhu3,wdhu33,modwdhu1
      REAL :: wdht1,wdht3,modwdht1,modwdht3
      REAL :: ave_wdhp1p0(nsm,nb,nb)
      REAL :: ave_wdhr11p0(nsm,nb,nb),ave_wdhr13p0(nsm,nb,nb)
      REAL :: ave_wdhp1b0(nsm,nb,nb)
      REAL :: ave_wdhr11b0(nsm,nb,nb),ave_wdhr13b0(nsm,nb,nb)
      REAL :: ave_wdhp1bp(nsm,nb,nb)
      REAL :: ave_wdhr11bp(nsm,nb,nb),ave_wdhr13bp(nsm,nb,nb)
      REAL :: ave_wdhu1(nsm,nb,nb),ave_wdhu3(nsm,nb,nb)
      REAL :: ave_wdhu33(nsm,nb,nb)
      REAL :: ave_modwdhu1(nsm,nb,nb)
      REAL :: ave_wdht1(nsm,nb,nb),ave_wdht3(nsm,nb,nb)
      REAL :: ave_modwdht1(nsm,nb,nb),ave_modwdht3(nsm,nb,nb)
! ave_wd_g
      REAL :: wdgp1p0,wdgr11p0,wdgr13p0
      REAL :: wdgp1b0,wdgr11b0,wdgr13b0
      REAL :: wdgp1bp,wdgr11bp,wdgr13bp
      REAL :: wdgu1,wdgu3,wdgu33,modwdgu1
      REAL :: wdgt1,wdgt3,modwdgt1,modwdgt3
      REAL :: ave_wdgp1p0(nsm,nb,nb)
      REAL :: ave_wdgr11p0(nsm,nb,nb),ave_wdgr13p0(nsm,nb,nb)
      REAL :: ave_wdgp1b0(nsm,nb,nb)
      REAL :: ave_wdgr11b0(nsm,nb,nb),ave_wdgr13b0(nsm,nb,nb)
      REAL :: ave_wdgp1bp(nsm,nb,nb)
      REAL :: ave_wdgr11bp(nsm,nb,nb),ave_wdgr13bp(nsm,nb,nb)
      REAL :: ave_wdgu1(nsm,nb,nb),ave_wdgu3(nsm,nb,nb)
      REAL :: ave_wdgu33(nsm,nb,nb)
      REAL :: ave_modwdgu1(nsm,nb,nb)
      REAL :: ave_wdgt1(nsm,nb,nb),ave_wdgt3(nsm,nb,nb)
      REAL :: ave_modwdgt1(nsm,nb,nb),ave_modwdgt3(nsm,nb,nb)
! kpar_h
      COMPLEX :: kpar_hnp0,kpar_hp1p0,kpar_hp3p0
      COMPLEX :: kpar_hp1b0,kpar_hr11b0,kpar_hr13b0
      COMPLEX :: kpar_hnbp,kpar_hp1bp,kpar_hp3bp
      COMPLEX :: kpar_hr11bp,kpar_hr13bp
      COMPLEX :: kpar_hu1,kpar_hu3
      COMPLEX :: kpar_hb1,kpar_hb3,kpar_hb33
      COMPLEX :: kpar_hb1ht1,kpar_hb3ht3,kpar_hb33ht1
      COMPLEX :: modkpar_hd1,modkpar_hd3,modkpar_hd33
      COMPLEX :: modkpar_hd1hu1,modkpar_hd3hu3,modkpar_hd33hu1
      COMPLEX :: modkpar_hu1,modkpar_hu3
      COMPLEX :: ave_kparhnp0(nsm,nb,nb),ave_kparhp1p0(nsm,nb,nb)
      COMPLEX :: ave_kparhp3p0(nsm,nb,nb)
      COMPLEX :: ave_kparhp1b0(nsm,nb,nb),ave_kparhr11b0(nsm,nb,nb)
      COMPLEX :: ave_kparhr13b0(nsm,nb,nb)
      COMPLEX :: ave_kparhnbp(nsm,nb,nb),ave_kparhp3bp(nsm,nb,nb)
      COMPLEX :: ave_kparhp1bp(nsm,nb,nb),ave_kparhr11bp(nsm,nb,nb)
      COMPLEX :: ave_kparhr13bp(nsm,nb,nb)
      COMPLEX :: ave_kparhu1(nsm,nb,nb),ave_kparhu3(nsm,nb,nb)
      COMPLEX :: ave_kparht1(nsm,nb,nb),ave_kparht3(nsm,nb,nb)
      COMPLEX :: ave_modkparhu1(nsm,nb,nb),ave_modkparhu3(nsm,nb,nb)
! kpar_g
      COMPLEX :: kpar_gnp0,kpar_gp1p0,kpar_gp3p0
      COMPLEX :: kpar_gp1b0,kpar_gr11b0,kpar_gr13b0
      COMPLEX :: kpar_gnbp,kpar_gp3bp
      COMPLEX :: kpar_gp1bp,kpar_gr11bp,kpar_gr13bp
      COMPLEX :: kpar_gu1,kpar_gu3
      COMPLEX :: kpar_gb1,kpar_gb3,kpar_gb33
      COMPLEX :: kpar_gb1gt1,kpar_gb3gt3,kpar_gb33gt1
      COMPLEX :: modkpar_gd1,modkpar_gd3,modkpar_gd33
      COMPLEX :: modkpar_gd1gu1,modkpar_gd3gu3,modkpar_gd33gu1
      COMPLEX :: modkpar_gu1,modkpar_gu3
      COMPLEX :: ave_kpargnp0(nsm,nb,nb),ave_kpargp1p0(nsm,nb,nb)
      COMPLEX :: ave_kpargp3p0(nsm,nb,nb)
      COMPLEX :: ave_kpargp1b0(nsm,nb,nb),ave_kpargr11b0(nsm,nb,nb)
      COMPLEX :: ave_kpargr13b0(nsm,nb,nb)
      COMPLEX :: ave_kpargnbp(nsm,nb,nb),ave_kpargp3bp(nsm,nb,nb)
      COMPLEX :: ave_kpargp1bp(nsm,nb,nb),ave_kpargr11bp(nsm,nb,nb)
      COMPLEX :: ave_kpargr13bp(nsm,nb,nb)
      COMPLEX :: ave_kpargu1(nsm,nb,nb),ave_kpargu3(nsm,nb,nb)
      COMPLEX :: ave_kpargt1(nsm,nb,nb),ave_kpargt3(nsm,nb,nb)
      COMPLEX :: ave_modkpargu1(nsm,nb,nb),ave_modkpargu3(nsm,nb,nb)
! gradB_h
      REAL :: gradBhp1,gradBhp3
      REAL :: gradBhr11,gradBhr13,gradBhr33
      REAL :: gradBhu1,gradBhu3,gradBhu33
      REAL :: ave_gradBhp1(nsm,nb,nb),ave_gradBhp3(nsm,nb,nb)
      REAL :: ave_gradBhr11(nsm,nb,nb),ave_gradBhr13(nsm,nb,nb)
      REAL :: ave_gradBhr33(nsm,nb,nb)
      REAL :: ave_gradBhu1(nsm,nb,nb),ave_gradBhu3(nsm,nb,nb)
      REAL :: ave_gradBhu33(nsm,nb,nb)
! gradB_g
      REAL :: gradBgp1,gradBgp3
      REAL :: gradBgr11,gradBgr13,gradBgr33
      REAL :: gradBgu1,gradBgu3,gradBgu33
      REAL :: ave_gradBgp1(nsm,nb,nb),ave_gradBgp3(nsm,nb,nb)
      REAL :: ave_gradBgr11(nsm,nb,nb),ave_gradBgr13(nsm,nb,nb)
      REAL :: ave_gradBgr33(nsm,nb,nb)
      REAL :: ave_gradBgu1(nsm,nb,nb),ave_gradBgu3(nsm,nb,nb)
      REAL :: ave_gradBgu33(nsm,nb,nb)
!  ave_theta
      REAL :: gradB
      REAL :: ave_kx(nb,nb)
      REAL :: ave_c_tor_par(nb,nb),ave_c_tor_per(nb,nb)
      REAL :: ave_c_par_par(nb,nb)
      REAL :: ave_wd(nb,nb),ave_modwd(nb,nb)
      REAL :: ave_gradB(nb,nb),ave_lnB(nb,nb)
      REAL :: ave_b0(nb,nb),ave_b0inv(nb,nb)
      REAL :: ave_kpar(nb,nb),ave_modkpar(nb,nb)
      REAL :: ave_p0(nb,nb),ave_p0inv(nb,nb)
      REAL :: ave_bp(nb,nb),ave_bpinv(nb,nb)
      COMPLEX :: ave_kpar_eff(nsm,nb,nb),ave_modkpar_eff(nsm,nb,nb)
!
      END MODULE tglf_coeffv194
!
!-----------------------------------------------------
!
      MODULE tglf_nuei_coeffv194
!
!     electron ion collision closure coefficients
!
      IMPLICIT NONE
!
      INTEGER,PARAMETER :: nc1=16,ncb=16
      REAL :: nuei_c1(nc1)
      REAL :: nuei_cb(ncb)
!
      DATA nuei_c1 /                                      &
       0.4919, 0.7535, 0.6727, 0.8055, 1.627, 2.013, &
       0.4972, 0.7805, 1.694, 3.103, 0.0, 0.0,         &
       0.0, 0.0, 0.0, 0.0/
!moment calc
!       0.60125, 1.1284, 0.79788, 1.1284, 1.79524, 2.2568, &
!       0.60125, 1.1284, 1.7952, 2.2568, 0.0, 0.0,         &
!       0.0, 0.0, 0.0, 0.0/
!
      END MODULE tglf_nuei_coeffv194
!
!______________________________________________________
      
      MODULE tglf_tgv194
!
! input data for the tglf model
!
      USE tglf_max_dimensionsv194
      IMPLICIT NONE
      SAVE
!
      INTEGER,PARAMETER :: nsmax=nsm,elitemax=max_ELITE,fouriermax=max_fourier
      LOGICAL :: use_TM_tg=.TRUE.
      LOGICAL :: iflux_tg=.TRUE.
      LOGICAL :: use_bper_tg=.FALSE.
      LOGICAL :: use_bpar_tg=.FALSE.
      LOGICAL :: use_mhd_rule_tg=.TRUE.
      LOGICAL :: use_bisection_tg=.TRUE.
      LOGICAL :: find_width_tg=.TRUE.
      LOGICAL :: adiabatic_elec_tg=.FALSE.
      LOGICAL :: new_eikonal_tg=.TRUE.
      INTEGER :: ibranch_tg=-1
      INTEGER :: nmodes_tg=2
      INTEGER :: nky_tg=12
      INTEGER :: ns_tg=2
      INTEGER :: igeo_tg=1
      INTEGER :: nbasis_max_tg=4
      INTEGER :: nbasis_min_tg=1
      INTEGER :: nxgrid_tg=16
      INTEGER :: nwidth_tg=21 
      INTEGER :: b_model_tg=0
      INTEGER :: ft_model_tg=1     
      INTEGER :: kygrid_model_tg=1
      INTEGER :: xnu_model_tg=1
      INTEGER :: sat_rule_tg = 0
      INTEGER :: vpar_model_tg=0
      INTEGER :: vpar_shear_model_tg=0
      REAL :: ky_tg=0.3
      REAL :: width_max_tg=1.65
      REAL :: width_min_tg=0.3
      REAL :: park_tg=1.0
      REAL :: ghat_tg=1.0
      REAL :: gchat_tg=1.0
      REAL :: betae_tg=0.0
      REAL :: damp_psi_tg=0.0
      REAL :: damp_sig_tg = 0.0
      REAL :: taui_tg=1.0
      REAL :: xnue_tg=0.0
      REAL :: rmin_tg=0.5
      REAL :: rmaj_tg=3.0
      REAL :: zeff_tg=1.0
      REAL :: debye_tg=0.0
      REAL :: vexb_shear_tg=0.0
      REAL :: vexb_tg
      REAL :: alpha_quench_tg=0.0
      REAL :: alpha_p_tg=1.0
      REAL :: alpha_e_tg=1.0
      REAL :: alpha_n_tg=0.0
      REAL :: alpha_t_tg=0.0
      REAL :: wd_zero_tg=0.1
      REAL :: Linsker_factor_tg=0.0
      REAL :: gradB_factor_tg=0.0
      REAL :: theta_trapped_tg=0.7
      REAL :: xnu_factor_tg=1.0
      REAL :: debye_factor_tg=1.0
      REAL :: etg_factor_tg = 1.25
      REAL :: filter_tg = 2.0
      REAL :: alpha_kx_e_tg=0.0
      REAL :: alpha_kx_p_tg=0.0
      REAL :: alpha_kx_n_tg=0.0
      REAL :: alpha_kx_t_tg=0.0
      REAL :: sign_Bt_tg=1.0
      REAL :: sign_It_tg=1.0
      REAL,DIMENSION(nsmax) :: taus_tg=1.0
      REAL,DIMENSION(nsmax) :: as_tg=1.0
      REAL,DIMENSION(nsmax) :: vpar_tg=0.0
      REAL,DIMENSION(nsmax) :: rlns_tg=1.0
      REAL,DIMENSION(nsmax) :: rlts_tg=3.0
      REAL,DIMENSION(nsmax) :: vpar_shear_tg=0.0
      REAL,DIMENSION(nsmax) :: zs_tg
      REAL,DIMENSION(nsmax) :: mass_tg
      REAL,DIMENSION(nsmax) :: vns_shear_tg=0.0
      REAL,DIMENSION(nsmax) :: vts_shear_tg=0.0
! Shifted cicle inputs
      REAL :: theta0_tg=0.0
      REAL :: shat_tg=1.0
      REAL :: alpha_tg=0.0
      REAL :: xwell_tg=0.0
      REAL :: q_tg=2.0
! Miller inputs
      REAL :: zmaj_tg=0.0
      REAL :: delta_tg=0.0
      REAL :: kappa_tg=1.0
      REAL :: zeta_tg=0.0
      REAL :: drmindx_tg=1.0
      REAL :: drmajdx_tg=0.0
      REAL :: dzmajdx_tg=0.0
      REAL :: s_delta_tg=0.0
      REAL :: s_kappa_tg=0.0
      REAL :: s_zeta_tg=0.0
      REAL :: q_prime_tg=0.0
      REAL :: p_prime_tg=0.0
      REAL :: kx0_tg=0.0
! Fourier inputs
      INTEGER :: nfourier_tg = 4
      REAL :: fourier_tg(8,0:fouriermax)=0.0
! ELITE inputs
      INTEGER :: n_elite_tg,n_surface_tg,j_surface_tg
      REAL :: R_elite_tg(0:elitemax)
      REAL :: Z_elite_tg(0:elitemax)
      REAL :: Bp_elite_tg(0:elitemax)
!
! namelist ....................................................

      NAMELIST /tglfin_nml/ adiabatic_elec_tg, find_width_tg, new_eikonal_tg, &
        nbasis_max_tg, nbasis_min_tg, nxgrid_tg, ibranch_tg, ns_tg, &
        nmodes_tg, iflux_tg, ky_tg, width_max_tg, width_min_tg, &
        nwidth_tg, park_tg, ghat_tg, gchat_tg, &
        alpha_e_tg, alpha_n_tg, alpha_t_tg, alpha_p_tg, &
        vexb_shear_tg, vpar_shear_tg, alpha_quench_tg, igeo_tg, theta_trapped_tg, &
        theta0_tg,taus_tg,as_tg,rlns_tg,rlts_tg,mass_tg,zs_tg, &
        rmin_tg, rmaj_tg, zmaj_tg,use_bisection_tg,vpar_tg, &
        q_tg, xnue_tg, wd_zero_tg, betae_tg, shat_tg, alpha_tg, &
        xwell_tg,kappa_tg,s_kappa_tg,delta_tg,s_delta_tg,zeta_tg,s_zeta_tg, &
        drmindx_tg,drmajdx_tg,dzmajdx_tg,zeff_tg, debye_tg, use_bper_tg, &
        use_bpar_tg,use_mhd_rule_tg,q_prime_tg,damp_psi_tg,damp_sig_tg, &
        p_prime_tg, filter_tg, Linsker_factor_tg, gradB_factor_tg,  &
        b_model_tg, ft_model_tg, xnu_factor_tg, debye_factor_tg, &
        nky_tg,etg_factor_tg,use_TM_tg,kygrid_model_tg,xnu_model_tg, &
        sat_rule_tg,alpha_kx_e_tg,alpha_kx_p_tg,alpha_kx_n_tg, alpha_kx_t_tg, &
        vpar_shear_model_tg, j_surface_tg,vpar_model_tg,sign_Bt_tg,sign_It_tg, &
        vns_shear_tg,vts_shear_tg, nfourier_tg,fourier_tg,vexb_tg,kx0_tg
!
      END MODULE tglf_tgv194
!______________________________________________________
!\
! module for MPI version of TGLF
!/

      MODULE tglf_mpiv194
#ifndef ONETWO
        use mpi

        IMPLICIT NONE

        ! local communicator
        integer :: iCommTglf    = -1
        integer :: iProcTglf    = -1
        integer :: nProcTglf    = -1
        integer :: iProc0Tglf   = 0
        integer :: iGroupIDTglf = 1
#endif
      END MODULE tglf_mpiv194
!
!______________________________________________________
!
