module expro_locsim_interface

  integer :: n_species_exp

  double precision :: mass_deuterium
  double precision, parameter :: temp_norm_fac   = 1602.2
  double precision, parameter :: charge_norm_fac = 1.6022

  double precision, dimension(:), allocatable :: rmin_exp

  double precision, dimension(:,:), allocatable :: temp_exp
  double precision, dimension(:,:), allocatable :: dens_exp
  double precision, dimension(:,:), allocatable :: dlntdr_exp
  double precision, dimension(:,:), allocatable :: dlnndr_exp
  double precision, dimension(:,:), allocatable :: sdlntdr_exp
  double precision, dimension(:,:), allocatable :: sdlnndr_exp
  double precision, dimension(:,:), allocatable :: sbeta_exp

  double precision, dimension(:), allocatable :: gamma_e_exp
  double precision, dimension(:), allocatable :: gamma_p_exp
  double precision, dimension(:), allocatable :: mach_exp
  
  double precision, dimension(:,:,:), allocatable :: geo_yin_exp

  ! Local values

  double precision :: shift_loc
  double precision :: q_loc
  double precision :: s_loc
  double precision :: kappa_loc
  double precision :: delta_loc
  double precision :: zeta_loc
  double precision :: shape_cos0_loc
  double precision :: shape_cos1_loc
  double precision :: shape_cos2_loc
  double precision :: shape_cos3_loc
  double precision :: shape_cos4_loc
  double precision :: shape_cos5_loc
  double precision :: shape_cos6_loc
  double precision :: shape_sin3_loc
  double precision :: shape_sin4_loc
  double precision :: shape_sin5_loc
  double precision :: shape_sin6_loc
  double precision :: s_kappa_loc
  double precision :: s_delta_loc
  double precision :: s_zeta_loc
  double precision :: shape_s_cos0_loc
  double precision :: shape_s_cos1_loc
  double precision :: shape_s_cos2_loc
  double precision :: shape_s_cos3_loc
  double precision :: shape_s_cos4_loc
  double precision :: shape_s_cos5_loc
  double precision :: shape_s_cos6_loc
  double precision :: shape_s_sin3_loc
  double precision :: shape_s_sin4_loc
  double precision :: shape_s_sin5_loc
  double precision :: shape_s_sin6_loc
  double precision :: zmag_loc
  double precision :: dzmag_loc
  double precision :: gamma_e_loc
  double precision :: gamma_p_loc
  double precision :: mach_loc
  double precision :: rmin_loc
  double precision :: rmaj_loc
  double precision :: rhos_loc
  double precision :: z_eff_loc
  double precision :: b_unit_loc
  double precision :: rho_norm_loc
  double precision :: psi_norm_loc
  double precision :: psi_a_loc
  double precision :: cs_loc
  double precision :: betae_loc
  double precision :: beta_star_loc
  double precision :: sbeta_loc

  double precision, dimension(9) :: mass_loc

  double precision, dimension(9) :: z_loc
  double precision, dimension(9) :: dens_loc
  double precision, dimension(9) :: temp_loc
  double precision, dimension(9) :: dlnndr_loc
  double precision, dimension(9) :: dlntdr_loc
  double precision, dimension(9) :: sdlnndr_loc
  double precision, dimension(9) :: sdlntdr_loc

  integer :: geo_ny_loc
  double precision, dimension(:,:), allocatable :: geo_yin_loc

end module expro_locsim_interface

subroutine expro_locsim_alloc(flag)

  use expro
  use expro_locsim_interface

  implicit none
  integer, intent(in) :: flag  

  ! flag=0: deallocate
  ! flag=1: allocate

  if (flag == 1) then

     allocate(rmin_exp(expro_n_exp))

     allocate(temp_exp(n_species_exp,expro_n_exp))
     allocate(dens_exp(n_species_exp,expro_n_exp))
     allocate(dlntdr_exp(n_species_exp,expro_n_exp))
     allocate(dlnndr_exp(n_species_exp,expro_n_exp))
     allocate(sdlntdr_exp(n_species_exp,expro_n_exp))
     allocate(sdlnndr_exp(n_species_exp,expro_n_exp))
     allocate(sbeta_exp(n_species_exp,expro_n_exp))

     allocate(gamma_e_exp(expro_n_exp))
     allocate(gamma_p_exp(expro_n_exp))
     allocate(mach_exp(expro_n_exp))

  else

     if (allocated(rmin_exp)) deallocate(rmin_exp)

     if (allocated(temp_exp)) deallocate(temp_exp)
     if (allocated(dens_exp)) deallocate(dens_exp)
     if (allocated(dlntdr_exp)) deallocate(dlntdr_exp)
     if (allocated(dlnndr_exp)) deallocate(dlnndr_exp)
     if (allocated(sdlntdr_exp)) deallocate(sdlntdr_exp)
     if (allocated(sdlnndr_exp)) deallocate(sdlnndr_exp)
     if (allocated(sbeta_exp)) deallocate(sbeta_exp)
     
     if (allocated(gamma_e_exp)) deallocate(gamma_e_exp)
     if (allocated(gamma_p_exp)) deallocate(gamma_p_exp)
     if (allocated(mach_exp)) deallocate(mach_exp)
     
     if (allocated(geo_yin_exp)) deallocate(geo_yin_exp)

  endif

end subroutine expro_locsim_alloc

!----------------------------------------------------------------
! expro_locsim_profiles.f90
!
! PURPOSE:
!  Read experimental profiles and generate local profile 
!  parameters at rmin = r/a.
!
! INPUTS:
!  path              : path to data
!  comm              : MPI communicator
!  numeq_flag        : Fourier series equilibrium (0=no,1=yes)
!  udsymmetry_flag   : enforce up-down symmetry (0=no,1=yes)
!  quasineutral_flag : enforce quasineutrality (0=no,1=yes)
!  n_species_in      : total species (e+i)
!  z                 : vector of charges (length n_species_in-1)
!  rmin              : r/a
!
! OUTPUTS:
!  btccw             : (+1 or -1) 
!  ipccw             : (+1 or -1)
!  a_meters          : minor radius i metres
!
! OUTPUTS (interface):
!  (see expro_locsim_interface)
!----------------------------------------------------------------

subroutine expro_locsim_profiles(&
     numeq_flag,&
     udsymmetry_flag,&
     quasineutral_flag,&
     n_species_in,&
     rmin,&
     btccw,&
     ipccw,&
     a_meters,&
     path,&
     comm)

  use expro
  use expro_locsim_interface

  implicit none

  integer, intent(in) :: numeq_flag
  integer, intent(in) :: udsymmetry_flag
  integer, intent(in) :: quasineutral_flag
  integer, intent(in) :: n_species_in
  integer, intent(in) :: comm
  character(len=80), intent(in) :: path
  double precision, intent(in) :: rmin
  double precision, intent(inout) :: btccw,ipccw,a_meters
  double precision :: rhostar
  
  integer :: i,j,i_ion

  rmin_loc = rmin

  n_species_exp = n_species_in

  mass_deuterium = expro_mass_deuterium*1e24

  !--------------------------------------------------------------
  ! use expro routines to read data:
  !
  expro_ctrl_quasineutral_flag = 1  ! quasi-neutrality density flag
  expro_ctrl_numeq_flag = numeq_flag
  expro_ctrl_n_ion = n_species_exp-1

  call expro_read(trim(path)//'input.gacode',comm)
  if (allocated(rmin_exp)) call expro_locsim_alloc(0)
  call expro_locsim_alloc(1)
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Transfer data from read vector to individual arrays:
  !
  btccw = -expro_signb
  ipccw = -expro_signq*expro_signb

  rmin_exp(:) = expro_rmin(:)

  if (udsymmetry_flag == 1) then
     expro_zmag(:) = 0d0   
     expro_dzmag(:) = 0d0
  endif

  ! Minor radius, a, in meters:
  a_meters = rmin_exp(expro_n_exp)

  rmin_exp(:) = rmin_exp(:)/a_meters

  ! Pack electrons into top of species vector
  temp_exp(n_species_exp,:)    = expro_te(:)
  dlntdr_exp(n_species_exp,:)  = expro_dlntedr(:)*a_meters 
  sdlntdr_exp(n_species_exp,:) = expro_sdlntedr(:)*a_meters
  dens_exp(n_species_exp,:)    = expro_ne(:)
  dlnndr_exp(n_species_exp,:)  = expro_dlnnedr(:)*a_meters 
  sdlnndr_exp(n_species_exp,:) = expro_sdlnnedr(:)*a_meters

  mass_loc(n_species_exp) = expro_masse
  z_loc(n_species_exp) = -1d0

  ! Pack ions from the bottom
  do i_ion=1,n_species_exp-1
     ! ion temps should be equal, but not enforced 
     temp_exp(i_ion,:)    = expro_ti(i_ion,:)
     dlntdr_exp(i_ion,:)  = expro_dlntidr(i_ion,:)*a_meters 
     sdlntdr_exp(i_ion,:) = expro_sdlntidr(i_ion,:)*a_meters 

     mass_loc(i_ion) = expro_mass(i_ion)
     z_loc(i_ion) = expro_z(i_ion)

     ! First species density is reset by quasi-neutrality
     if (quasineutral_flag == 1 .and. i_ion == 1) then
        dens_exp(i_ion,:)    = expro_ni_new(:)
        dlnndr_exp(i_ion,:)  = expro_dlnnidr_new(:)*a_meters
        sdlnndr_exp(i_ion,:) = expro_sdlnnidr_new(:)*a_meters
     else
        dens_exp(i_ion,:)    = expro_ni(i_ion,:)
        dlnndr_exp(i_ion,:)  = expro_dlnnidr(i_ion,:)*a_meters
        sdlnndr_exp(i_ion,:) = expro_sdlnnidr(i_ion,:)*a_meters
     endif

  enddo

  ! Rotation
  gamma_e_exp(:) = -expro_w0p(:)*(a_meters*rmin_exp(:))/expro_q(:)
  gamma_p_exp(:) = -expro_w0p(:)*expro_rmaj(:)
  mach_exp(:)    = expro_w0(:)*expro_rmaj(:)

  !------------------------------------------------------------------
  ! Use local cubic spline interpolation to get simulation 
  ! profiles from experimental (_exp) ones.
  !
  call cub_spline1(rmin_exp,expro_rmaj/a_meters,expro_n_exp,rmin,rmaj_loc)
  call cub_spline1(rmin_exp,expro_q,expro_n_exp,rmin,q_loc)
  call cub_spline1(rmin_exp,expro_s,expro_n_exp,rmin,s_loc)
  call cub_spline1(rmin_exp,expro_drmaj,expro_n_exp,rmin,shift_loc)
  call cub_spline1(rmin_exp,expro_kappa,expro_n_exp,rmin,kappa_loc)
  call cub_spline1(rmin_exp,expro_skappa,expro_n_exp,rmin,s_kappa_loc)
  call cub_spline1(rmin_exp,expro_delta,expro_n_exp,rmin,delta_loc)
  call cub_spline1(rmin_exp,expro_sdelta,expro_n_exp,rmin,s_delta_loc)
  call cub_spline1(rmin_exp,expro_zeta,expro_n_exp,rmin,zeta_loc)
  call cub_spline1(rmin_exp,expro_szeta,expro_n_exp,rmin,s_zeta_loc)
  call cub_spline1(rmin_exp,expro_shape_cos0,expro_n_exp,rmin,shape_cos0_loc)
  call cub_spline1(rmin_exp,expro_shape_scos0,expro_n_exp,rmin,shape_s_cos0_loc)
  call cub_spline1(rmin_exp,expro_shape_cos1,expro_n_exp,rmin,shape_cos1_loc)
  call cub_spline1(rmin_exp,expro_shape_scos1,expro_n_exp,rmin,shape_s_cos1_loc)
  call cub_spline1(rmin_exp,expro_shape_cos2,expro_n_exp,rmin,shape_cos2_loc)
  call cub_spline1(rmin_exp,expro_shape_scos2,expro_n_exp,rmin,shape_s_cos2_loc)
  call cub_spline1(rmin_exp,expro_shape_cos3,expro_n_exp,rmin,shape_cos3_loc)
  call cub_spline1(rmin_exp,expro_shape_scos3,expro_n_exp,rmin,shape_s_cos3_loc)
  call cub_spline1(rmin_exp,expro_shape_cos4,expro_n_exp,rmin,shape_cos4_loc)
  call cub_spline1(rmin_exp,expro_shape_scos4,expro_n_exp,rmin,shape_s_cos4_loc)
  call cub_spline1(rmin_exp,expro_shape_cos5,expro_n_exp,rmin,shape_cos5_loc)
  call cub_spline1(rmin_exp,expro_shape_scos5,expro_n_exp,rmin,shape_s_cos5_loc)
  call cub_spline1(rmin_exp,expro_shape_cos6,expro_n_exp,rmin,shape_cos6_loc)
  call cub_spline1(rmin_exp,expro_shape_scos6,expro_n_exp,rmin,shape_s_cos6_loc)
  call cub_spline1(rmin_exp,expro_shape_sin3,expro_n_exp,rmin,shape_sin3_loc)
  call cub_spline1(rmin_exp,expro_shape_ssin3,expro_n_exp,rmin,shape_s_sin3_loc)
  call cub_spline1(rmin_exp,expro_shape_sin4,expro_n_exp,rmin,shape_sin4_loc)
  call cub_spline1(rmin_exp,expro_shape_ssin4,expro_n_exp,rmin,shape_s_sin4_loc)
  call cub_spline1(rmin_exp,expro_shape_sin5,expro_n_exp,rmin,shape_sin5_loc)
  call cub_spline1(rmin_exp,expro_shape_ssin5,expro_n_exp,rmin,shape_s_sin5_loc)
  call cub_spline1(rmin_exp,expro_shape_sin6,expro_n_exp,rmin,shape_sin6_loc)
  call cub_spline1(rmin_exp,expro_shape_ssin6,expro_n_exp,rmin,shape_s_sin6_loc)
  call cub_spline1(rmin_exp,expro_zmag/a_meters,expro_n_exp,rmin,zmag_loc)
  call cub_spline1(rmin_exp,expro_dzmag,expro_n_exp,rmin,dzmag_loc)
  call cub_spline1(rmin_exp,gamma_e_exp,expro_n_exp,rmin,gamma_e_loc)
  call cub_spline1(rmin_exp,gamma_p_exp,expro_n_exp,rmin,gamma_p_loc)
  call cub_spline1(rmin_exp,mach_exp,expro_n_exp,rmin,mach_loc)
  call cub_spline1(rmin_exp,expro_rhos,expro_n_exp,rmin,rhos_loc)
  call cub_spline1(rmin_exp,expro_cs,expro_n_exp,rmin,cs_loc)
  call cub_spline1(rmin_exp,expro_z_eff,expro_n_exp,rmin,z_eff_loc)
  call cub_spline1(rmin_exp,expro_bunit,expro_n_exp,rmin,b_unit_loc)
  call cub_spline1(rmin_exp,expro_rho,expro_n_exp,rmin,rho_norm_loc)
  call cub_spline1(rmin_exp,expro_polflux,expro_n_exp,rmin,psi_norm_loc)
  psi_norm_loc = psi_norm_loc/expro_polflux(expro_n_exp)
  psi_a_loc = expro_polflux(expro_n_exp)

  ! rhos/a for curvature corrections
  rhostar = rhos_loc/a_meters

  beta_star_loc = 0d0
  sbeta_loc = 0d0
  
  do i=1,n_species_exp
     ! Note: mapping is only done for n_species (not n_species_exp)
     call cub_spline1(rmin_exp,dens_exp(i,:),expro_n_exp,rmin,dens_loc(i))
     call cub_spline1(rmin_exp,temp_exp(i,:),expro_n_exp,rmin,temp_loc(i))
     call cub_spline1(rmin_exp,dlntdr_exp(i,:),expro_n_exp,rmin,dlntdr_loc(i))
     call cub_spline1(rmin_exp,dlnndr_exp(i,:),expro_n_exp,rmin,dlnndr_loc(i))
     call cub_spline1(rmin_exp,sdlntdr_exp(i,:),expro_n_exp,rmin,sdlntdr_loc(i))
     call cub_spline1(rmin_exp,sdlnndr_exp(i,:),expro_n_exp,rmin,sdlnndr_loc(i))
     beta_star_loc = beta_star_loc+dens_loc(i)*temp_loc(i)*(dlnndr_loc(i)+dlntdr_loc(i))
     sbeta_loc = sbeta_loc+dens_loc(i)*temp_loc(i)*(sdlnndr_loc(i)+sdlntdr_loc(i)-2*dlnndr_loc(i)*dlntdr_loc(i)*rhostar)
  enddo

  ! beta calculation in CGS:
  !
  !         8*pi ( n[1e19/m^3]*1e-6*1e19 )( T[keV]*1.6022*1e-9 )
  ! beta = ------------------------------------------------------
  !                           ( 1e4*B[T] )^2
  !
  !      = 4.027e-3 n[1e19/m^3]*T[keV]/B[T]^2

  ! beta
  betae_loc = 4.027e-3*dens_loc(n_species_exp)*temp_loc(n_species_exp)/b_unit_loc**2
  ! beta_* (local)
  beta_star_loc = beta_star_loc*betae_loc/(dens_loc(n_species_exp)*temp_loc(n_species_exp))
  ! beta_S (global)
  sbeta_loc = sbeta_loc*betae_loc/(dens_loc(n_species_exp)*temp_loc(n_species_exp))

  ! Define effective curvatures used in global formulation
  sdlnndr_loc(:) = sdlnndr_loc(:)+(1.5*dlntdr_loc(:)**2+dlnndr_loc(:)**2-dlnndr_loc(:)*dlntdr_loc(:))*rhostar
  sdlntdr_loc(:) = sdlntdr_loc(:)+dlntdr_loc(:)**2*rhostar
  
  if (numeq_flag == 1 .and. expro_nfourier > 0) then

     geo_ny_loc = expro_nfourier
     allocate(geo_yin_exp(8,0:geo_ny_loc,expro_n_exp))
     if(allocated(geo_yin_loc)) deallocate(geo_yin_loc)
     allocate(geo_yin_loc(8,0:geo_ny_loc))
     geo_yin_exp(1:4,:,:) = expro_geo(:,:,:)/a_meters
     geo_yin_exp(5:8,:,:) = expro_dgeo(:,:,:)

     do i=1,8
        do j=0,geo_ny_loc
           call cub_spline1(rmin_exp,geo_yin_exp(i,j,:),expro_n_exp,rmin, &
                geo_yin_loc(i,j))
        enddo
     enddo

  endif

end subroutine expro_locsim_profiles

!---------------------------------------------------------------
! cub_spline.f90
!
! PURPOSE:
!  Take known points x(1:n),y(1:n) and perform cubic spline
!  interpolation at the node points xi(1:ni) to give yi(1:ni).
!
!  Use 'natural' cubic spline interpolation, where natural 
!  means the second derivative is zero at the boundaries.
!
!  INPUT  : x(1:n),y(1:n),n,xi(1:ni),ni
!
!  OUTPUT : yi(1:ni)
!
!  The only requirements are 
!  
!  1. Ordered data: x(i+1) > x(i) and xi(i+1) > xi(i).
!  2. Bounded data: xi(1) < x(1) and xi(ni) > x(n).
!----------------------------------------------------------------

subroutine cub_spline(x,y,n,xi,yi,ni)

  !-------------------------------------------------------------
  implicit none
  !
  integer :: i,ii
  double precision :: x0 
  !
  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: x,y
  !
  integer, intent(in) :: ni
  double precision, dimension(ni) :: xi,yi
  !
  ! LAPACK working variables
  !
  integer :: info
  integer, dimension(n) :: ipiv 
  double precision, dimension(n) :: c,z
  double precision, dimension(n-1) :: zl,zu,h
  double precision, dimension(n-2) :: zu2
  double precision, dimension(n-1) :: b,d
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Check to see that interpolated point is inside data interval
  !
  if (xi(ni) > x(n)) then
     print *,'ERROR: (cub_spline) Data above upper bound'
     print *,'xi(ni) > x(n)',xi(ni),x(n)
     print *,'y(:)',y(:)
     stop
  endif 
  if (xi(1) < x(1)) then
     print *,'ERROR: (cub_spline) Data below lower bound'
     print *,'xi(1) < x(1)',xi(1),x(1) 
     stop
  endif
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Define coefficients of spline matrix 
  !
  do i=1,n-1 
     h(i)  = x(i+1)-x(i)
     zl(i) = h(i)
     zu(i) = h(i)
  enddo
  zl(n-1) = 0d0
  zu(1)   = 0d0

  z(1) = 1d0
  c(1) = 0d0
  do i=2,n-1
     z(i) = 2d0*(h(i-1)+h(i))
     c(i) = 3d0*((y(i+1)-y(i))/h(i)-(y(i)-y(i-1))/h(i-1))
  enddo
  z(n) = 1d0
  c(n) = 0d0
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Solve the system using LAPACK
  !
  call DGTTRF(n,zl,z,zu,zu2,ipiv,info)
  call DGTTRS('N',n,1,zl,z,zu,zu2,ipiv,c,n,info)
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Find remaining polynomial coefficients:
  !
  c(n) = 0d0
  do i=1,n-1
     b(i) = (y(i+1)-y(i))/h(i)-h(i)*(2d0*c(i)+c(i+1))/3d0
     d(i) = (c(i+1)-c(i))/(3d0*h(i))
  enddo
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Using known polynomial coefficients, perform interpolation.
  !
  !  S(x) = y(i) + b(i) [x-x(i)] + c(i) [x-x(i)]^2 
  !                                        + d(i) [x-x(i)]^3
  !
  i  = 1
  ii = 1
  do while (ii <= ni)
     x0 = xi(ii)
     if (x0 <= x(i+1)) then
        yi(ii) = y(i)+(x0-x(i))*(b(i)+(x0-x(i))*(c(i)+(x0-x(i))*d(i)))
        ii = ii+1
     else
        i = i+1
     endif
  enddo
  !-------------------------------------------------------------

end subroutine cub_spline

subroutine cub_spline1(x,y,n,xi,yi)

  !-------------------------------------------------------------
  implicit none
  !
  integer :: i,ii
  !
  integer, intent(in) :: n
  double precision, intent(in), dimension(n) :: x,y
  !
  double precision, intent(in) :: xi
  double precision, intent(inout) :: yi
  !
  ! LAPACK working variables
  !
  integer :: info
  integer, dimension(n) :: ipiv 
  double precision, dimension(n) :: c,z
  double precision, dimension(n-1) :: zl,zu,h
  double precision, dimension(n-2) :: zu2
  double precision, dimension(n-1) :: b,d
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Check to see that interpolated point is inside data interval
  !
  if (xi > x(n)) then
     print '(a)','ERROR: (cub_spline1) Data above upper bound'
     print '(a)','xi(ni) > x(n)',xi,x(n) 
  endif 
  if (xi < x(1)) then
     print '(a)','ERROR: (cub_spline1) Data below lower bound'
     print '(a)','xi(1) < x(1)',xi,x(1) 
  endif 
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Define coefficients of spline matrix 
  !
  do i=1,n-1 
     h(i)  = x(i+1)-x(i)
     zl(i) = h(i)
     zu(i) = h(i)
  enddo
  zl(n-1) = 0d0
  zu(1)   = 0d0

  z(1) = 1d0
  c(1) = 0d0
  do i=2,n-1
     z(i) = 2d0*(h(i-1)+h(i))
     c(i) = 3d0*((y(i+1)-y(i))/h(i)-(y(i)-y(i-1))/h(i-1))
  enddo
  z(n) = 1d0
  c(n) = 0d0
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Solve the system using LAPACK
  !
  call DGTTRF(n,zl,z,zu,zu2,ipiv,info)
  call DGTTRS('N',n,1,zl,z,zu,zu2,ipiv,c,n,info)
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Find remaining polynomial coefficients:
  !
  c(n) = 0d0
  do i=1,n-1
     b(i) = (y(i+1)-y(i))/h(i)-h(i)*(2d0*c(i)+c(i+1))/3d0
     d(i) = (c(i+1)-c(i))/(3d0*h(i))
  enddo
  !-------------------------------------------------------------

  !-------------------------------------------------------------
  ! Using known polynomial coefficients, perform interpolation.
  !
  !  S(x) = y(i) + b(i) [x-x(i)] + c(i) [x-x(i)]^2 
  !                                        + d(i) [x-x(i)]^3
  !
  i  = 1
  ii = 1
  do while (ii <= 1)
     if (xi <= x(i+1)) then
        yi = y(i)+(xi-x(i))*(b(i)+(xi-x(i))*(c(i)+(xi-x(i))*d(i)))
        ii = ii+1
     else
        i = i+1
     endif
  enddo
  !-------------------------------------------------------------

end subroutine cub_spline1
