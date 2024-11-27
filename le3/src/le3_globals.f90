module le3_globals

  real, parameter :: pi=3.141592653589793

  logical :: initialized = .false.

  ! Resolution parameters
  integer :: nt,np
  integer :: nts,nps
  integer :: m,n

  ! Geometry parameters
  integer :: equilibrium_model
  real :: rmin,rmaj,hmin,dhmindr,q,s
  real :: kappa, s_kappa, delta, s_delta, zeta, s_zeta
  real :: shift, zmag, dzmag
  real :: beta_star
  real :: iota,iota_p
  real :: dt,dp

  integer :: nts_geo, nps_geo
  real, dimension(:,:,:), allocatable :: r_geo, z_geo, rd_geo, zd_geo

  real, dimension(:), allocatable    :: t,p
  real, dimension(:,:), allocatable  :: tb
  real, dimension(:,:), allocatable  :: dtbdp,dtbdt
  real, dimension(:,:), allocatable  :: r,z
  real, dimension(:,:), allocatable  :: drdtb,dzdtb
  real, dimension(:,:), allocatable  :: drdpb,dzdpb
  real, dimension(:,:), allocatable  :: drdt,dzdt
  real, dimension(:,:), allocatable  :: drdp,dzdp
  real, dimension(:,:), allocatable  :: jac
  real, dimension(:,:), allocatable  :: rp,rt
  real, dimension(:,:), allocatable  :: zp,zt
  real, dimension(:,:), allocatable  :: fp,ft
  real, dimension(:,:), allocatable  :: as,bs,cs,ds
  real, dimension(:,:), allocatable  :: sinm,cosm
  real, dimension(:,:), allocatable  :: sinn,cosn

  real, dimension(:,:), allocatable :: gpp,gtt,gpt
  real, dimension(:,:), allocatable :: gct, gcp, gcc
  real, dimension(:,:), allocatable :: cosu,rc
  real, dimension(:,:), allocatable :: bpol,btor
  real, dimension(:,:), allocatable :: chi1
  real, dimension(:,:), allocatable :: chi1p,chi1t
  real, dimension(:,:), allocatable :: drdpt,dzdpt
  real, dimension(:,:), allocatable :: dbdt,dbdp

  real, dimension(:), allocatable :: deriv_t
  real, dimension(:), allocatable :: deriv_p

  ! "Second order" quantities
  real, dimension(:,:), allocatable :: dtheta
  real, dimension(:,:), allocatable :: dchi
  real, dimension(:,:), allocatable :: dthetap
  real, dimension(:,:), allocatable :: b1
  real, dimension(:,:), allocatable :: gc

  !--------------------------------------------------------
  ! MINPACK variables
  integer :: info
  integer :: nwork
  integer :: msize
  real, dimension(:), allocatable :: work
  !
  real :: tol
  real, dimension(:), allocatable :: xfunc
  real, dimension(:), allocatable :: yfunc
  !--------------------------------------------------------

  real, dimension(:,:), allocatable :: bmag
  real, dimension(:,:), allocatable :: bdotgrad
  real, dimension(:,:), allocatable :: bdotgradB_overB
  real, dimension(:,:), allocatable :: dgdp,vexb_dt,vexb_dp
  real, dimension(:,:), allocatable :: vdrift_x, vdrift_dt, vdrift_dp
  real, dimension(:,:,:), allocatable :: grad_perpsq_gk, vdrift_gk
  real, dimension(:,:), allocatable :: grad_cc, grad_tt, grad_pp, grad_pt, grad_ct, grad_cp
  real, dimension(:,:), allocatable :: g
  real, dimension(:,:), allocatable :: mat_stream_dt
  real, dimension(:,:), allocatable :: mat_stream_dp
  real, dimension(:,:), allocatable :: mat_trap
  real, dimension(:,:), allocatable :: mat_coll
  real, dimension(:,:), allocatable :: mat_vexb_dt, mat_vexb_dp
  real, dimension(:,:), allocatable :: mat_vdrift_dt, mat_vdrift_dp
  real :: vprime
  integer :: matsize
  integer :: indx_c00
  real, dimension(:,:), allocatable :: bk,bkp
  real, dimension(:,:), allocatable :: bk_t,bkp_t
  real, dimension(:,:), allocatable :: bk_p,bkp_p
  integer, dimension(:), allocatable :: m_indx, n_indx, itype

end module le3_globals
