!------------------------------------------------
! EXPRO_interface.f90
!
! PURPOSE:
!  Interface to experimental profiles.
!
! NOTES:
!  For more detailed variable definitions, see
!
!  http://fusion.gat.com/theory/input.profiles
!  
!  * Control parameters (user can change these)
! 
!  EXPRO_ctrl_quasineutral_flag (0=do nothing, 1=quasineutral)
!  EXPRO_ctrl_z(1:10) (ion charges)
!  EXPRO_ctrl_numeq_flag (0=model,1=numerical)
!  EXPRO_ctrl_silent_flag (0=normal,1=silent)
!
!  * Derived quantities:
!
!  EXPRO_bunit(:)       B_unit (T)
!  EXPRO_s(:)           (r/q)(dq/dr)
!  EXPRO_drmaj(:)       dR_0/dr
!  EXPRO_dzmag(:)       dZ_0/dr
!  EXPRO_sdelta(:)      r d(delta)/dr
!  EXPRO_skappa(:)      (r/kappa) d(kappa)/dr
!  EXPRO_szeta(:)       r d(zeta)/dr
!  EXPRO_dlnnedr(:)     -dln(ne)/dr (1/m)
!  EXPRO_dlntedr(:)     -dln(Te)/dr (1/m)
!  EXPRO_dlnnidr(1:5,:) -dln(ni)/dr (1/m)
!  EXPRO_dlntidr(1:5,:) -dln(ti)/dr (1/m)
!  EXPRO_sdlnnedr(:)     -ne''/ne (1/m^2)
!  EXPRO_sdlntedr(:)     -Te''/Te (1/m^2)
!  EXPRO_sdlnnidr(1:5,:) -ni''/ni /dr (1/m^2)
!  EXPRO_sdlntidr(1:5,:) -Ti''/Ti (1/m^2)
!  EXPRO_dlnptotdr(:)   -dln(ptot)/dr (1/m)
!  EXPRO_w0p(:)         d(w0)/dr (1/s/m)
!  EXPRO_vol(:)         V (m^3)
!  EXPRO_volp(:)        dV/dr (m^2)
!  EXPRO_cs(:)          cs (m/s)
!  EXPRO_rhos(:)        rhos (m)
!  EXPRO_ni_new(:)      ni [Corrected for quasin.]
!  EXPRO_dlnnidr_new(:) -dln(ni)/dr (1/m) [Corrected for quasin.]
!  EXPRO_sdlnnidr_new(:) -ni''/ni /dr (1/m^2) [Corrected for quasin.]
!  EXPRO_grad_r0(:)     |grad r| at theta=0 
!  EXPRO_ave_grad_r(:)  Flux-surface average <|grad r|> 
!  EXPRO_drdrho(:)      dr/d(rho) [r and rho have units of length]
!  EXPRO_bp0(:)         B_pol at theta=0 (T)
!  EXPRO_bt0(:)         B_tor at theta=0 (T)
!  EXPRO_ip(:)          Plasma current (A)
!  EXPRO_gamma_e(:)     r/q d(w0)/dr (1/s)
!  EXPRO_gamma_p(:)     R_0 d(w0)/dr (1/s)
!  EXPRO_mach(:)        R_0 w0/cs
!  EXPRO_thetascale(:)  max(gsin') [measure of poloidal scale length]
!
!  * General geometry arrays:
!
!  EXPRO_geo(:,:,:)
!  EXPRO_dgeo(:,:,:)
!------------------------------------------------ 

module EXPRO_interface

  integer :: EXPRO_error=0
  integer, parameter :: EXPRO_n_ion_max=10

  ! Fundamental input.profiles scalars

  integer :: EXPRO_shot=0
  integer :: EXPRO_n_ion=0
  integer :: EXPRO_n_exp=0
  real    :: EXPRO_b_ref
  real    :: EXPRO_ip_exp
  real    :: EXPRO_rvbv
  real    :: EXPRO_arho

  ! Fundamental input.profiles arrays

  character (len=16) :: EXPRO_null_tag = '[null]'

  real, dimension(:),allocatable :: EXPRO_rho
  character (len=16) :: EXPRO_rho_tag = 'rho(-)'

  real, dimension(:),allocatable :: EXPRO_rmin
  character (len=16) :: EXPRO_rmin_tag = 'rmin(m)'

  real, dimension(:),allocatable :: EXPRO_rmaj
  character (len=16) :: EXPRO_rmaj_tag = 'rmaj(m)'

  real, dimension(:),allocatable :: EXPRO_q
  character (len=16) :: EXPRO_q_tag = 'q(-)'

  real, dimension(:),allocatable :: EXPRO_kappa
  character (len=16) :: EXPRO_kappa_tag = 'kappa(-)'

  real, dimension(:),allocatable :: EXPRO_delta
  character (len=16) :: EXPRO_delta_tag = 'delta(-)'

  real, dimension(:),allocatable :: EXPRO_z_eff
  character (len=16) :: EXPRO_z_eff_tag = 'z_eff(-)'

  real, dimension(:),allocatable :: EXPRO_w0
  character (len=16) :: EXPRO_w0_tag = 'omega0(rad/s)'

  real, dimension(:),allocatable :: EXPRO_flow_mom
  character (len=16) :: EXPRO_flow_mom_tag = 'flow_mom(N-m)'

  real, dimension(:),allocatable :: EXPRO_sbcx
  character (len=16) :: EXPRO_sbcx_tag = 'sbcx(/m^3/s)'

  real, dimension(:),allocatable :: EXPRO_sbeame
  character (len=16) :: EXPRO_sbeame_tag = 'sbeame(/m^3/s)'

  real, dimension(:),allocatable :: EXPRO_sscxl
  character (len=16) :: EXPRO_sscxl_tag = 'sscxl(/m^3/s)'

  real, dimension(:),allocatable :: EXPRO_pow_e
  character (len=16) :: EXPRO_pow_e_tag = 'pow_e(MW)'

  real, dimension(:),allocatable :: EXPRO_pow_i
  character (len=16) :: EXPRO_pow_i_tag = 'pow_i(MW)'

  real, dimension(:),allocatable :: EXPRO_pow_ei
  character (len=16) :: EXPRO_pow_ei_tag = 'pow_ei(MW)'

  real, dimension(:),allocatable :: EXPRO_zeta
  character (len=16) :: EXPRO_zeta_tag = 'zeta(-)'

  real, dimension(:),allocatable :: EXPRO_flow_beam
  character (len=16) :: EXPRO_flow_beam_tag = 'flow_beam(kW/eV)'

  real, dimension(:),allocatable :: EXPRO_flow_wall
  character (len=16) :: EXPRO_flow_wall_tag = 'flow_wall(kW/eV)'

  real, dimension(:),allocatable :: EXPRO_zmag
  character (len=16) :: EXPRO_zmag_tag = 'zmag(m)'

  real, dimension(:),allocatable :: EXPRO_ptot
  character (len=16) :: EXPRO_ptot_tag = 'ptot(Pa)'

  real, dimension(:),allocatable :: EXPRO_polflux
  character (len=16) :: EXPRO_polflux_tag = 'polflux(Wb/rad)'

  real, dimension(:),allocatable :: EXPRO_pow_e_fus
  character (len=16) :: EXPRO_pow_e_fus_tag = 'pow_e_fus(MW)'

  real, dimension(:),allocatable :: EXPRO_pow_i_fus
  character (len=16) :: EXPRO_pow_i_fus_tag = 'pow_i_fus(MW)'

  real, dimension(:),allocatable :: EXPRO_pow_e_sync
  character (len=16) :: EXPRO_pow_e_sync_tag = 'pow_e_sync(MW)'

  real, dimension(:),allocatable :: EXPRO_pow_e_brem
  character (len=16) :: EXPRO_pow_e_brem_tag = 'pow_e_brem(MW)'

  real, dimension(:),allocatable :: EXPRO_pow_e_line
  character (len=16) :: EXPRO_pow_e_line_tag = 'pow_e_line(MW)'

  real, dimension(:),allocatable :: EXPRO_pow_e_aux
  character (len=16) :: EXPRO_pow_e_aux_tag = 'pow_e_aux(MW)'

  real, dimension(:),allocatable :: EXPRO_pow_i_aux
  character (len=16) :: EXPRO_pow_i_aux_tag = 'pow_i_aux(MW)'


  real, dimension(:),allocatable :: EXPRO_ne
  character (len=16) :: EXPRO_ne_tag = 'ne(10^19/m^3)'

  real, dimension(:,:),allocatable :: EXPRO_ni
  character (len=16), dimension(EXPRO_n_ion_max) :: EXPRO_ni_tag = (/ &
       'ni_1(10^19/m^3) ',&
       'ni_2(10^19/m^3) ',&
       'ni_3(10^19/m^3) ',&
       'ni_4(10^19/m^3) ',&
       'ni_5(10^19/m^3) ',&
       'ni_6(10^19/m^3) ',&
       'ni_7(10^19/m^3) ',&
       'ni_8(10^19/m^3) ',&
       'ni_9(10^19/m^3) ',&
       'ni_10(10^19/m^3)'/)

  real, dimension(:),allocatable :: EXPRO_te
  character (len=16) :: EXPRO_te_tag = 'Te(keV)'

  real, dimension(:,:),allocatable :: EXPRO_ti
  character (len=16), dimension(EXPRO_n_ion_max) :: EXPRO_ti_tag = (/ &
       'Ti_1(keV)       ',&
       'Ti_2(keV)       ',&
       'Ti_3(keV)       ',&
       'Ti_4(keV)       ',&
       'Ti_5(keV)       ',&
       'Ti_6(keV)       ',&
       'Ti_7(keV)       ',&
       'Ti_8(keV)       ',&
       'Ti_9(keV)       ',&
       'Ti_10(keV)      '/)

  real, dimension(:,:),allocatable :: EXPRO_vtor
  character (len=16), dimension(EXPRO_n_ion_max) :: EXPRO_vtor_tag = (/ &
       'vtor_1(m/s)     ',&
       'vtor_2(m/s)     ',&
       'vtor_3(m/s)     ',&
       'vtor_4(m/s)     ',&
       'vtor_5(m/s)     ',&
       'vtor_6(m/s)     ',&
       'vtor_7(m/s)     ',&
       'vtor_8(m/s)     ',&
       'vtor_9(m/s)     ',&
       'vtor_10(m/s)    '/)

  real, dimension(:,:),allocatable :: EXPRO_vpol
  character (len=16), dimension(EXPRO_n_ion_max) :: EXPRO_vpol_tag = (/ &
       'vpol_1(m/s)     ',&
       'vpol_2(m/s)     ',&
       'vpol_3(m/s)     ',&
       'vpol_4(m/s)     ',&
       'vpol_5(m/s)     ',&
       'vpol_6(m/s)     ',&
       'vpol_7(m/s)     ',&
       'vpol_8(m/s)     ',&
       'vpol_9(m/s)     ',&
       'vpol_10(m/s)    '/)

  ! Derived quantities

  real, dimension(:),allocatable :: EXPRO_bunit
  real, dimension(:),allocatable :: EXPRO_s
  real, dimension(:),allocatable :: EXPRO_drmaj
  real, dimension(:),allocatable :: EXPRO_dzmag
  real, dimension(:),allocatable :: EXPRO_sdelta
  real, dimension(:),allocatable :: EXPRO_skappa
  real, dimension(:),allocatable :: EXPRO_szeta
  real, dimension(:),allocatable :: EXPRO_dlnnedr
  real, dimension(:),allocatable :: EXPRO_dlntedr
  real, dimension(:,:),allocatable :: EXPRO_dlnnidr
  real, dimension(:,:),allocatable :: EXPRO_dlntidr
  real, dimension(:),allocatable :: EXPRO_sdlnnedr
  real, dimension(:),allocatable :: EXPRO_sdlntedr
  real, dimension(:,:),allocatable :: EXPRO_sdlnnidr
  real, dimension(:,:),allocatable :: EXPRO_sdlntidr
  real, dimension(:),allocatable :: EXPRO_dlnptotdr

  real, dimension(:),allocatable :: EXPRO_w0p

  real, dimension(:),allocatable :: EXPRO_vol
  real, dimension(:),allocatable :: EXPRO_volp

  real, dimension(:),allocatable :: EXPRO_cs
  real, dimension(:),allocatable :: EXPRO_rhos

  real, dimension(:),allocatable :: EXPRO_ni_new
  real, dimension(:),allocatable :: EXPRO_dlnnidr_new
  real, dimension(:),allocatable :: EXPRO_sdlnnidr_new
  real, dimension(:),allocatable :: EXPRO_grad_r0
  real, dimension(:),allocatable :: EXPRO_ave_grad_r
  real, dimension(:),allocatable :: EXPRO_drdrho

  real, dimension(:),allocatable :: EXPRO_bp0
  real, dimension(:),allocatable :: EXPRO_bt0
  real, dimension(:),allocatable :: EXPRO_ip

  ! input.profiles.geo dimension and arrays

  integer :: EXPRO_nfourier
  real, dimension(:,:,:),allocatable :: EXPRO_geo
  real, dimension(:,:,:),allocatable :: EXPRO_dgeo

  ! Gyrokinetic rotation parameters

  real, dimension(:),allocatable :: EXPRO_gamma_e
  real, dimension(:),allocatable :: EXPRO_gamma_p
  real, dimension(:),allocatable :: EXPRO_mach

  ! Scale length

  real, dimension(:),allocatable :: EXPRO_thetascale

  ! Field orientation parameters

  integer :: EXPRO_signb
  integer :: EXPRO_signq

  !------------------------------------------------------------------------

  ! Control parameters (force nonsensical default -1 for usage check)

  integer :: EXPRO_ctrl_quasineutral_flag = -1
  integer :: EXPRO_ctrl_n_ion = -1
  real, dimension(EXPRO_n_ion_max) :: EXPRO_ctrl_z = 0.0
  integer :: EXPRO_ctrl_numeq_flag = -1
  integer :: EXPRO_ctrl_silent_flag = 0
 
end module EXPRO_interface
