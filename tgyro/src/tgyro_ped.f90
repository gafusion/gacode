!-------------------------------------------------------------------
! tgyro_ped.f90
!
! PURPOSE:
!  Manage dynamic pedestal interface including no-man's-land (NML).
!
! CRITICAL RADII:
!  Let w be the pedestal width in psi_norm:
!
!  top -> 1-2.5*w
!  ped -> 1-2.0*w
!  sym -> 1-1.0*w
!-------------------------------------------------------------------

module tgyro_ped

  implicit none

  ! EPED_NN inputs
  real :: a_in
  real :: betan_in
  real :: bt_in
  real :: delta_in
  real :: ip_in
  real :: kappa_in
  real :: m_in
  real :: nped_in
  real :: r_in
  real :: zeffped_in
  
  real :: t_axis

  ! EPED_NN outputs
  integer, parameter :: nx_nn=1001
  real :: nn_vec(nx_nn,3)
  real :: nn_w_ped

  ! Derivatives wrt psi
  real :: t_vec(nx_nn)
  real :: n_p(nx_nn)
  real :: t_p(nx_nn)

  ! Copies of EXPRO variables
  integer :: n_exp
  !
  real, dimension(:), allocatable :: rmin_exp
  real, dimension(:), allocatable :: psi_exp
  real, dimension(:), allocatable :: dpsidr_exp
  real, dimension(:), allocatable :: volp_exp
  real, dimension(:), allocatable :: ptot_exp
  real, dimension(:), allocatable :: exp_ne
  real, dimension(:), allocatable :: exp_te
  real, dimension(:,:), allocatable :: exp_ni
  real, dimension(:,:), allocatable :: exp_ti
  real, dimension(:), allocatable :: exp_w0

  ! Pedestal top scale lengths
  real :: zn_top,zt_top
  real :: dr_nml
  real, dimension(1) :: r_top,psi_top
  real, dimension(1) :: n_top,t_top,p_top
  real :: n_frac,t_frac
  
contains

  subroutine tgyro_pedestal

    use mpi
    use tgyro_globals

    implicit none

    ! Parameters interpolated at top of pedestal
    real, dimension(1) :: n_p_top,t_p_top
    real, dimension(1) :: dpsidr_top
    real :: ntsum
    integer :: i_ion

    if (tgyro_ped_model == 1) return

    !-------------------------------------------------------------------------
    ! 1. Initializations 
    !
    ! nped_in = <n>
    nped_in    = tgyro_neped*(1.0+sum(n_ratio(1:loc_n_ion)))/2.0
    zeffped_in = tgyro_zeffped
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 2. Call EPED1NN 
    !
    if (i_proc_global == 0) then

       ! Inputs stored in interface  
       call tgyro_eped_nn
       ! Outputs (nx_nn is large, fixed number of radial points)
       !  nn_vec(1:nx_nn,1) -> psi_norm
       !  nn_vec(1:nx_nn,2) -> ne [1/cm^3]
       !  nn_vec(1:nx_nn,3) -> P [Pa]

       if (tgyro_rped < 0.0) then
          psi_top(1) = 1.0 - 1.5*nn_w_ped*abs(tgyro_rped)
       else
          psi_top(1) = tgyro_rped
       endif

  endif

    ! Broadcast needed data from single-task call above 
    call MPI_BCAST(nn_vec,size(nn_vec),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(psi_top,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

    ! NOTE: Below, n_top means ne_top
    !              t_top means Te_top 

    ! n_top: in 1/cm^3
    call cub_spline(nn_vec(:,1),nn_vec(:,2),nx_nn,psi_top,n_top,1)
    ! p_top: in Pa
    call cub_spline(nn_vec(:,1),nn_vec(:,3),nx_nn,psi_top,p_top,1) 

    ! Pressure formula based on pivot assumption (includes tiny fast ion bits):
    ! 
    ! P/k = ne Te + Sum_i ni Ti
    !     = n_top t_top [ 1 + Sum_i n_ratio(i) t_ratio(i) ]
    !     = n_top t_top ntsum
    !     = p_top/k
 
    !ntsum = 1.0+sum(n_ratio(1:loc_n_ion)*t_ratio(1:loc_n_ion))

    ! t_top [eV] 
    t_top = (10.0*p_top)/(2*n_top*k) 

    ! Calculate n' (n_p) and T' (t_p):
    call bound_deriv(n_p,nn_vec(:,2),nn_vec(:,1),nx_nn)
    t_vec = (10.0*nn_vec(:,3))/(2*nn_vec(:,2)*k)
    call bound_deriv(t_p,t_vec,nn_vec(:,1),nx_nn)

    ! n_top', t_top'
    call cub_spline(nn_vec(:,1),n_p,nx_nn,psi_top,n_p_top,1)
    call cub_spline(nn_vec(:,1),t_p,nx_nn,psi_top,t_p_top,1) 
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 3. Compute z parameters at pedestal top
    ! 
    ! zn_top = -(1/n) dn/dr 
    ! zt_top = -(1/T) dT/dr 
    !
    call cub_spline(psi_exp,rmin_exp  ,n_exp,psi_top,     r_top,1)
    call cub_spline(psi_exp,dpsidr_exp,n_exp,psi_top,dpsidr_top,1)

    zn_top = -n_p_top(1)/n_top(1)*dpsidr_top(1)
    zt_top = -t_p_top(1)/t_top(1)*dpsidr_top(1)
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! 4. Integrate to obtain TGYRO pivot values
    
    ! Integration backward from r_top to r_star
    dr_nml = r_top(1)-r(n_r)
    ne(n_r)   = n_top(1)*exp(0.5*(dlnnedr(n_r)  +zn_top)*dr_nml)*n_frac
    te(n_r)   = t_top(1)*exp(0.5*(dlntedr(n_r)  +zt_top)*dr_nml)*t_frac

    ! Self-similar ion values (fixed pivot assumption; see tgyro_profile_reintegrate)
    do i_ion=1,loc_n_ion
       if (therm_flag(i_ion) == 1) then
          ni(i_ion,n_r) = n_top(1)*exp(0.5*(dlnnidr(i_ion,n_r)+zn_top)*dr_nml)*n_ratio(i_ion)*n_frac
          ti(i_ion,n_r) = t_top(1)*exp(0.5*(dlntidr(i_ion,n_r)+zt_top)*dr_nml)*t_ratio(i_ion)*t_frac
       endif
    enddo
    !-------------------------------------------------------------------------

  end subroutine tgyro_pedestal

  !---------------------------------------------------------------------------
  ! tgyro_pededtal_map.f90
  !
  ! PURPOSE:
  !  Perform critical profile interpolation from core to edge.
  !
  ! Inputs
  !
  !  z_star : scale length at r_star (pivot)
  !  z_top  : scale length at r_top  (pedestal top)
  !  f_top  : function at r_top
  !  p_vec  : profile vector from pedestal top to edge (r_star < r < a)
  !
  ! Outputs
  !
  ! i_star : first i_exp such that r > r_star
  ! f_exp  : Global profile on experimental (input.profiles) grid 
  !---------------------------------------------------------------------------

  subroutine tgyro_pedestal_map(z_star,z_top,f_top,p_vec,i_star,f_exp)

    use tgyro_globals

    implicit none

    real, intent(in) :: z_star,z_top,f_top
    real, intent(in) :: p_vec(nx_nn)
    integer, intent(inout) :: i_star
    real, intent(inout) :: f_exp(n_exp)

    integer :: i_exp,i0
    real :: x0

    ! 1. Scale-length interpolation over NML (r_star < r < r_top)
    !                    zb                         za
    ! f(r) = f(rb)*exp[ ---- ( dr^2 - (r-ra)^2 ) + ---- ( rb-r )^2 ]
    !                   2 dr                       2 dr
    !
    i_star = 0
    do i_exp=2,n_exp
       x0 = rmin_exp(i_exp)
       if (x0 > r(n_r) .and. x0 <= r_top(1)) then 
          ! Calculate i_star [first i_exp such that r > r(n_r)]
          if (i_star == 0) i_star = i_exp
          f_exp(i_exp) = f_top*exp(&
               0.5*z_top/dr_nml*(dr_nml**2-(x0-r(n_r))**2)+&
               0.5*z_star/dr_nml*(r_top(1)-x0)**2)
       endif
       if (x0 > r_top(1)) then
          i0 = i_exp
          exit
       endif
    enddo

    ! 2. Direct spline interpolation over pedestal (r_top < r < a)
    call cub_spline(nn_vec(:,1),p_vec,nx_nn,psi_exp(i0:n_exp),f_exp(i0:n_exp),n_exp-i0+1) 

  end subroutine tgyro_pedestal_map

end module tgyro_ped
