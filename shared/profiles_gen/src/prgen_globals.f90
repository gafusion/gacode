module prgen_globals

  ! Parameters
  real, parameter :: kevdsecpmw=1.6022e-19*1e3*1e-6
  real, parameter :: pi=3.14159265358979323846

  ! List of possible data files
  character (len=70) :: date
  character (len=70) :: raw_data_file
  character (len=70) :: raw_data_type
  character (len=70) :: cer_file

  !----------------------------------------------------------
  ! Internal control variables
  !
  integer, parameter :: n_ion_max = 10
  integer :: nx
  integer :: format_type
  integer :: efit_method
  integer :: noq_flag
  integer :: nop_flag
  integer :: verbose_flag
  integer :: gmerge_flag
  character (len=70) :: efit_header
  integer :: ipccw
  integer :: btccw
  integer :: nfourier
  integer :: lump_fast_flag
  real :: dpsi_efit
  real :: dpsi_data
  integer :: n_reorder
  integer, dimension(n_ion_max) :: reorder_vec
  integer :: n_null
  integer :: n_lump
  integer, dimension(:), allocatable :: lump_vec
  integer :: n_ion_sanitize
  character(len=2), dimension(:), allocatable :: ion_sanitize
  !----------------------------------------------------------

  ! Internal physics variables
  real, dimension(:), allocatable :: rho
  real, dimension(:), allocatable :: dpsi
  real, dimension(:), allocatable :: rmin
  real, dimension(:), allocatable :: rmaj
  real, dimension(:), allocatable :: q
  real, dimension(:), allocatable :: p_tot
  real, dimension(:), allocatable :: vphi_exp
  real, dimension(:), allocatable :: vphi_carbon
  real, dimension(:), allocatable :: vpolc_exp
  real, dimension(:), allocatable :: vtorc_exp
  real, dimension(:), allocatable :: omega0
  real, dimension(:), allocatable :: kappa
  real, dimension(:), allocatable :: delta
  real, dimension(:), allocatable :: zmag
  real, dimension(:), allocatable :: zeta
  real, dimension(:), allocatable :: pow_e
  real, dimension(:), allocatable :: pow_i
  real, dimension(:), allocatable :: pow_ei
  real, dimension(:), allocatable :: pow_e_fus
  real, dimension(:), allocatable :: pow_i_fus
  real, dimension(:), allocatable :: pow_e_sync
  real, dimension(:), allocatable :: pow_e_brem
  real, dimension(:), allocatable :: pow_e_line
  real, dimension(:), allocatable :: pow_e_rad
  real, dimension(:), allocatable :: pow_e_nb
  real, dimension(:), allocatable :: pow_i_nb
  real, dimension(:), allocatable :: pow_e_rf
  real, dimension(:), allocatable :: pow_i_rf
  real, dimension(:), allocatable :: pow_e_ohm
  real, dimension(:), allocatable :: pow_e_aux
  real, dimension(:), allocatable :: pow_i_aux
  real, dimension(:), allocatable :: powe_ion_exp
  real, dimension(:), allocatable :: powi_ion_exp
  real, dimension(:), allocatable :: powe_wdot_exp
  real, dimension(:), allocatable :: powi_wdot_exp
  real, dimension(:), allocatable :: powi_cx_exp
  real, dimension(:), allocatable :: flow_wall_exp
  real, dimension(:), allocatable :: flow_beam
  real, dimension(:), allocatable :: flow_mom
  real, dimension(:), allocatable :: sion_d
  real, dimension(:), allocatable :: sbcx_d

  ! Ion name,type,mass,charge
  integer :: n_ion
  character (len=6), dimension(n_ion_max) :: ion_name
  character (len=7), dimension(n_ion_max) :: ion_type
  real, dimension(n_ion_max) :: ion_mass
  integer, dimension(n_ion_max) :: ion_z
  character (len=7) :: type_therm = '[therm]'
  character (len=7) :: type_fast  = '[fast] '

  real :: quasi_err
  real :: pow_e_err
  real :: pow_i_err
  real :: ip_tot
  !---------------------------------------------------------
  ! ONETWO variables
  !
  integer :: onetwo_ishot
  integer :: onetwo_nb
  integer :: onetwo_npsi
  integer :: onetwo_nion
  integer :: onetwo_nion_tot
  integer :: onetwo_nprim
  integer :: onetwo_nimp
  integer :: onetwo_nneu
  integer :: onetwo_nbion
  integer :: onetwo_nalp
  character (len=2), dimension(n_ion_max) :: onetwo_namep
  character (len=2), dimension(n_ion_max) :: onetwo_namei
  character (len=2), dimension(n_ion_max) :: onetwo_nameb
  character (len=2), dimension(15) :: onetwo_ion_name
  real, dimension(n_ion_max) :: onetwo_m
  integer, dimension(n_ion_max) :: onetwo_z
  real :: onetwo_time
  real :: onetwo_Rgeom
  real :: onetwo_Rmag
  real :: onetwo_R0
  real :: onetwo_kappa
  real :: onetwo_delta
  real :: onetwo_volo
  real :: onetwo_cxareao
  real :: onetwo_Btor
  !
  real, dimension(:), allocatable :: onetwo_rho_grid
  real, dimension(:), allocatable :: onetwo_volume
  real, dimension(:), allocatable :: onetwo_hcap
  real, dimension(:), allocatable :: onetwo_qbeame
  real, dimension(:), allocatable :: onetwo_qrfe
  real, dimension(:), allocatable :: onetwo_qohm
  real, dimension(:), allocatable :: onetwo_qrad
  real, dimension(:), allocatable :: onetwo_qione
  real, dimension(:), allocatable :: onetwo_qioni
  real, dimension(:), allocatable :: onetwo_dpedt
  real, dimension(:), allocatable :: onetwo_qfuse
  real, dimension(:), allocatable :: onetwo_qbeami
  real, dimension(:), allocatable :: onetwo_qrfi
  real, dimension(:), allocatable :: onetwo_qcx
  real, dimension(:), allocatable :: onetwo_dpidt
  real, dimension(:), allocatable :: onetwo_qfusi
  real, dimension(:), allocatable :: onetwo_qdelt
  real, dimension(:), allocatable :: onetwo_sbeam
  real, dimension(:), allocatable :: onetwo_sbeame
  real, dimension(:), allocatable :: onetwo_sscxl
  real, dimension(:), allocatable :: onetwo_te
  real, dimension(:), allocatable :: onetwo_ti
  real, dimension(:), allocatable :: onetwo_talp
  real, dimension(:), allocatable :: onetwo_q
  real, dimension(:), allocatable :: onetwo_angrot
  real, dimension(:), allocatable :: onetwo_zeff
  real, dimension(:), allocatable :: onetwo_ene
  real, dimension(:), allocatable :: onetwo_enalp
  real, dimension(:), allocatable :: onetwo_psi
  real, dimension(:), allocatable :: onetwo_storqueb
  real, dimension(:,:), allocatable :: onetwo_enion
  real, dimension(:,:), allocatable :: onetwo_enbeam
  real, dimension(:,:), allocatable :: onetwo_pressb
  real, dimension(:,:), allocatable :: onetwo_sion
  real, dimension(:,:), allocatable :: onetwo_srecom
  real, dimension(:,:), allocatable :: onetwo_scx
  real, dimension(:,:), allocatable :: onetwo_sbcx
  real, dimension(:,:), allocatable :: onetwo_s
  real, dimension(:,:), allocatable :: onetwo_dudt
  !
  ! psi-grid vectors
  !
  real, dimension(:), allocatable :: onetwo_rho_mhd_gridnpsi
  real, dimension(:), allocatable :: onetwo_rmajavnpsi
  real, dimension(:), allocatable :: onetwo_rminavnpsi
  real, dimension(:), allocatable :: onetwo_psivolpnpsi
  real, dimension(:), allocatable :: onetwo_elongxnpsi
  real, dimension(:), allocatable :: onetwo_triangnpsi_u
  real, dimension(:), allocatable :: onetwo_triangnpsi_l
  !
  ! conglomerate ion vectors
  !
  real, dimension(:,:), allocatable :: onetwo_enion_vec
  real, dimension(:,:), allocatable :: onetwo_Tion_vec
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! PLASMA STATE variables
  !
  integer :: plst_shot_number
  character (len=32) :: plst_tokamak_id
  integer :: plst_dim_nrho
  integer :: plst_dim_nrho_eq_geo
  integer :: plst_dp1_nspec_th
  integer :: plst_dp1_nspec_tha
  integer :: plst_dp1_nspec_all
  integer :: plst_dp1_nspec_alla
  character (len=32) :: plst_tag
  character (len=32), dimension(:), allocatable :: plst_all_name
  character (len=32), dimension(:), allocatable :: plst_alla_name
  real :: plst_b_axis_vac
  real, dimension(:), allocatable :: plst_q_all
  real, dimension(:), allocatable :: plst_m_all
  real, dimension(:,:), allocatable :: plst_ts
  real, dimension(:,:), allocatable :: plst_ns
  real, dimension(:), allocatable :: plst_nb
  real, dimension(:), allocatable :: plst_nmini
  real, dimension(:), allocatable :: plst_nfusi
  real, dimension(:), allocatable :: plst_eperpb
  real, dimension(:), allocatable :: plst_eparb
  real, dimension(:), allocatable :: plst_vol
  real, dimension(:), allocatable :: plst_rho
  real, dimension(:), allocatable :: plst_grho1
  real, dimension(:), allocatable :: plst_phit
  real, dimension(:), allocatable :: plst_psipol
  real, dimension(:), allocatable :: plst_elong
  real, dimension(:), allocatable :: plst_triang
  real, dimension(:), allocatable :: plst_iota
  real, dimension(:), allocatable :: plst_r_midp_in
  real, dimension(:), allocatable :: plst_r_midp_out
  real, dimension(:), allocatable :: plst_z_midp
  real, dimension(:), allocatable :: plst_zeff
  real, dimension(:), allocatable :: plst_epot
  real, dimension(:), allocatable :: plst_omegat
  real, dimension(:), allocatable :: plst_pe_trans
  real, dimension(:), allocatable :: plst_pi_trans
  real, dimension(:), allocatable :: plst_pbe
  real, dimension(:), allocatable :: plst_pbi
  real, dimension(:), allocatable :: plst_pbth
  real, dimension(:), allocatable :: plst_qie
  real, dimension(:), allocatable :: plst_peech
  real, dimension(:), allocatable :: plst_pohme
  real, dimension(:), allocatable :: plst_pmine
  real, dimension(:), allocatable :: plst_pmini
  real, dimension(:), allocatable :: plst_pminth
  real, dimension(:), allocatable :: plst_picth
  real, dimension(:), allocatable :: plst_pfusi
  real, dimension(:), allocatable :: plst_pfusth
  real, dimension(:), allocatable :: plst_pfuse
  real, dimension(:), allocatable :: plst_prad_br
  real, dimension(:), allocatable :: plst_prad_cy
  real, dimension(:), allocatable :: plst_prad_li
  real, dimension(:), allocatable :: plst_tq_trans
  real, dimension(:), allocatable :: plst_sn_trans
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! PEQDSK variables
  !
  real :: peqdsk_bref
  real :: peqdsk_arho
  real, dimension(:), allocatable :: peqdsk_psi
  real, dimension(:), allocatable :: peqdsk_ne
  real, dimension(:), allocatable :: peqdsk_te
  real, dimension(:), allocatable :: peqdsk_ni
  real, dimension(:), allocatable :: peqdsk_ti
  real, dimension(:,:), allocatable :: peqdsk_nz
  real, dimension(:), allocatable :: peqdsk_omegat
  real, dimension(:), allocatable :: peqdsk_omgeb
  real, dimension(:), allocatable :: peqdsk_nb
  real, dimension(:), allocatable :: peqdsk_pb
  integer :: peqdsk_fmt, peqdsk_nion, peqdsk_nimp, peqdsk_nbeams
  real, dimension(n_ion_max) :: peqdsk_z
  real, dimension(n_ion_max) :: peqdsk_m
  character (len=7), dimension(5)  :: peqdsk_type
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! null input variables
  !
  real :: null_bref
  real :: null_arho
  !---------------------------------------------------------

  !--------------------------------------------------------
  ! Some GATO variables
  !
  integer :: nsurf
  integer :: narc
  !--------------------------------------------------------

  !---------------------------------------------------------
  ! CORSICA variables
  !
  real :: corsica_time
  real :: corsica_current
  real :: corsica_loop_voltage
  real :: corsica_fusion_power
  real :: corsica_bootstrap_fraction
  real :: corsica_betat
  real :: corsica_betan
  real :: corsica_taue
  real :: corsica_q95
  real :: corsica_r
  real :: corsica_a
  real :: corsica_kappa
  real :: corsica_delta
  real :: corsica_fusion_gain
  real :: corsica_betap
  real :: corsica_li3

  integer, parameter :: corsica_nvals = 81
  real :: corsica_bref
  real :: corsica_arho
  real, allocatable :: corsica_rho(:)
  real, allocatable :: corsica_r_a(:)
  real, allocatable :: corsica_psin(:)
  real, allocatable :: corsica_vl(:)
  real, allocatable :: corsica_te(:)
  real, allocatable :: corsica_ti(:)
  real, allocatable :: corsica_ne(:)
  real, allocatable :: corsica_ndt(:)
  real, allocatable :: corsica_nz(:)
  real, allocatable :: corsica_nalpha(:)
  real, allocatable :: corsica_zeff(:)
  real, allocatable :: corsica_q(:)
  real, allocatable :: corsica_j(:)
  real, allocatable :: corsica_jbs(:)
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! UFILE variables
  !
  integer :: ufile_nion
  real :: ufile_bref
  real :: ufile_arho
  real :: ufile_z(8)
  real :: ufile_m(8)
  character (len=10) :: ufile_tok
  character (len=40) :: ufile_shot
  character (len=10) :: ufile_time
  character (len=5), dimension(8)  :: ufile_type
  real, dimension(:), allocatable :: ufile_ne
  real, dimension(:,:), allocatable :: ufile_ni
  real, dimension(:), allocatable :: ufile_te
  real, dimension(:,:), allocatable :: ufile_ti
  real, dimension(:), allocatable :: ufile_zeff
  real, dimension(:), allocatable :: ufile_pres
  real, dimension(:), allocatable :: ufile_vrot
  real, dimension(:), allocatable :: ufile_vrotm
  real, dimension(:), allocatable :: ufile_volume
  real, dimension(:), allocatable :: ufile_qnbii
  real, dimension(:), allocatable :: ufile_qnbie
  real, dimension(:), allocatable :: ufile_qicrhi
  real, dimension(:), allocatable :: ufile_qicrhe
  real, dimension(:), allocatable :: ufile_qei
  real, dimension(:), allocatable :: ufile_qrad
  real, dimension(:), allocatable :: ufile_qeche
  real, dimension(:), allocatable :: ufile_qechi
  real, dimension(:), allocatable :: ufile_qohm
  real, dimension(:), allocatable :: ufile_qwalli
  real, dimension(:), allocatable :: ufile_qwalle

  real, dimension(:), allocatable :: ufile_qfuse
  real, dimension(:), allocatable :: ufile_qfusi
  real, dimension(:), allocatable :: ufile_qlhe
  real, dimension(:), allocatable :: ufile_qlhi
  real, dimension(:), allocatable :: ufile_snbie
  real, dimension(:), allocatable :: ufile_snbii
  real, dimension(:), allocatable :: ufile_swall
  real, dimension(:), allocatable :: ufile_torq
  !---------------------------------------------------------

end module prgen_globals
