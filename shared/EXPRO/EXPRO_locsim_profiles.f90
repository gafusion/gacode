!----------------------------------------------------------------
! EXPRO_locsim_profiles.f90
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
!  real :: shift_loc
!  real :: q_loc
!  real :: s_loc
!  real :: kappa_loc
!  real :: delta_loc
!  real :: zeta_loc
!  real :: s_kappa_loc
!  real :: s_delta_loc
!  real :: s_zeta_loc
!  real :: zmag_loc
!  real :: dzmag_loc
!  real :: gamma_e_loc
!  real :: gamma_p_loc
!  real :: mach_loc
!  real :: rmaj_loc
!  real :: rhos_loc [m]
!  real :: z_eff_loc
!  real :: b_unit_loc
!  real :: rho_norm_loc
!  real :: psi_norm_loc
!  real :: psi_a_loc
!
!  real, dimension(9) :: dens_loc
!  real, dimension(9) :: temp_loc
!  real, dimension(9) :: dlnndr_loc
!  real, dimension(9) :: dlntdr_loc
!  real, dimension(9) :: sdlnndr_loc
!  real, dimension(9) :: sdlntdr_loc
!----------------------------------------------------------------

subroutine EXPRO_locsim_profiles(&
     path,&
     comm,&
     numeq_flag,&
     udsymmetry_flag,&
     quasineutral_flag,&
     n_species_in,&
     z,&
     rmin,&
     btccw,&
     ipccw,&
     a_meters)

  use EXPRO_locsim_interface
  use EXPRO_interface

  implicit none

  character(len=*), intent(in) :: path 
  integer, intent(in) :: comm
  integer, intent(in) :: numeq_flag
  integer, intent(in) :: udsymmetry_flag
  integer, intent(in) :: quasineutral_flag
  integer, intent(in) :: n_species_in
  real, intent(in), dimension(n_species_in-1) :: z
  real, intent(in) :: rmin
  real, intent(inout) :: btccw,ipccw,a_meters

  integer :: i,j,is,i_ion

  n_species_exp = n_species_in

  !--------------------------------------------------------------
  ! use EXPRO routines to read data:
  !
  call EXPRO_palloc(comm,path,1)
  EXPRO_ctrl_quasineutral_flag = 1  ! quasi-neutrality density flag
  EXPRO_ctrl_numeq_flag = numeq_flag

  ! Number and charge of ion species
  EXPRO_ctrl_z(:) = 0.0
  EXPRO_ctrl_n_ion = n_species_exp-1
  do is=1,EXPRO_ctrl_n_ion
     EXPRO_ctrl_z(is) = z(is)
  enddo

  call EXPRO_pread
  call EXPRO_locsim_alloc(1)
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Transfer data from read vector to individual arrays:
  !
  btccw = -EXPRO_signb
  ipccw = -EXPRO_signq*EXPRO_signb

  rmin_exp(:) = EXPRO_rmin(:)

  if (udsymmetry_flag == 1) then
     EXPRO_zmag(:) = 0.0   
     EXPRO_dzmag(:) = 0.0
  endif

  ! Minor radius, a, in meters:
  a_meters = rmin_exp(EXPRO_n_exp)

  rmin_exp(:) = rmin_exp(:)/a_meters

  ! Pack electrons into top of species vector.
  temp_exp(n_species_exp,:)    = EXPRO_te(:)
  dlntdr_exp(n_species_exp,:)  = EXPRO_dlntedr(:)*a_meters 
  sdlntdr_exp(n_species_exp,:) = EXPRO_sdlntedr(:)*a_meters**2 
  dens_exp(n_species_exp,:)    = EXPRO_ne(:)
  dlnndr_exp(n_species_exp,:)  = EXPRO_dlnnedr(:)*a_meters 
  sdlnndr_exp(n_species_exp,:) = EXPRO_sdlnnedr(:)*a_meters**2 

  ! Pack ions from the bottom
  do i_ion=1,n_species_exp-1
     ! ion temps should be equal, but not enforced 
     temp_exp(i_ion,:)    = EXPRO_ti(i_ion,:)
     dlntdr_exp(i_ion,:)  = EXPRO_dlntidr(i_ion,:)*a_meters 
     sdlntdr_exp(i_ion,:) = EXPRO_sdlntidr(i_ion,:)*a_meters**2 

     ! First species density is reset by quasi-neutrality
     if (quasineutral_flag == 1 .and. i_ion == 1) then
        dens_exp(i_ion,:)    = EXPRO_ni_new(:)
        dlnndr_exp(i_ion,:)  = EXPRO_dlnnidr_new(:)*a_meters
        sdlnndr_exp(i_ion,:) = EXPRO_sdlnnidr_new(:)*a_meters**2
     else
        dens_exp(i_ion,:)    = EXPRO_ni(i_ion,:)
        dlnndr_exp(i_ion,:)  = EXPRO_dlnnidr(i_ion,:)*a_meters
        sdlnndr_exp(i_ion,:) = EXPRO_sdlnnidr(i_ion,:)*a_meters**2
     endif
  enddo

  ! Sanity check for densities
  do i=1,n_species_exp
     if (minval(dens_exp(i,:)) <= 0.0) then
        !call cgyro_error('Nonpositive in exp. density profile')
        return
     endif
  enddo

  ! Rotation
  gamma_e_exp(:) = -EXPRO_w0p(:)*(a_meters*rmin_exp(:))/EXPRO_q(:)
  gamma_p_exp(:) = -EXPRO_w0p(:)*EXPRO_rmaj(:)
  mach_exp(:)    = EXPRO_w0(:)*EXPRO_rmaj(:)

  !------------------------------------------------------------------
  ! Use local cubic spline interpolation to get simulation 
  ! profiles from experimental (_exp) ones.
  ! 
  call cub_spline(rmin_exp,EXPRO_rmaj/a_meters,EXPRO_n_exp,rmin,rmaj_loc,1)
  call cub_spline(rmin_exp,EXPRO_q,EXPRO_n_exp,rmin,q_loc,1)
  call cub_spline(rmin_exp,EXPRO_s,EXPRO_n_exp,rmin,s_loc,1)
  call cub_spline(rmin_exp,EXPRO_drmaj,EXPRO_n_exp,rmin,shift_loc,1)
  call cub_spline(rmin_exp,EXPRO_kappa,EXPRO_n_exp,rmin,kappa_loc,1)
  call cub_spline(rmin_exp,EXPRO_skappa,EXPRO_n_exp,rmin,s_kappa_loc,1)
  call cub_spline(rmin_exp,EXPRO_delta,EXPRO_n_exp,rmin,delta_loc,1)
  call cub_spline(rmin_exp,EXPRO_sdelta,EXPRO_n_exp,rmin,s_delta_loc,1)
  call cub_spline(rmin_exp,EXPRO_zeta,EXPRO_n_exp,rmin,zeta_loc,1)
  call cub_spline(rmin_exp,EXPRO_szeta,EXPRO_n_exp,rmin,s_zeta_loc,1)
  call cub_spline(rmin_exp,EXPRO_zmag/a_meters,EXPRO_n_exp,rmin,zmag_loc,1)
  call cub_spline(rmin_exp,EXPRO_dzmag,EXPRO_n_exp,rmin,dzmag_loc,1)
  call cub_spline(rmin_exp,gamma_e_exp,EXPRO_n_exp,rmin,gamma_e_loc,1)
  call cub_spline(rmin_exp,gamma_p_exp,EXPRO_n_exp,rmin,gamma_p_loc,1)
  call cub_spline(rmin_exp,mach_exp,EXPRO_n_exp,rmin,mach_loc,1)
  call cub_spline(rmin_exp,EXPRO_rhos,EXPRO_n_exp,rmin,rhos_loc,1)
  call cub_spline(rmin_exp,EXPRO_z_eff,EXPRO_n_exp,rmin,z_eff_loc,1)
  call cub_spline(rmin_exp,EXPRO_bunit,EXPRO_n_exp,rmin,b_unit_loc,1)
  call cub_spline(rmin_exp,EXPRO_rho,EXPRO_n_exp,rmin,rho_norm_loc,1)
  call cub_spline(rmin_exp,EXPRO_polflux,EXPRO_n_exp,rmin,psi_norm_loc,1)
  psi_norm_loc  = psi_norm_loc/EXPRO_polflux(EXPRO_n_exp)
  psi_a_loc = EXPRO_polflux(EXPRO_n_exp)
  
  do i=1,n_species_exp
     ! Note: mapping is only done for n_species (not n_species_exp)
     call cub_spline(rmin_exp,dens_exp(i,:),EXPRO_n_exp,rmin,dens_loc(i),1)
     call cub_spline(rmin_exp,temp_exp(i,:),EXPRO_n_exp,rmin,temp_loc(i),1)
     call cub_spline(rmin_exp,dlntdr_exp(i,:),EXPRO_n_exp,rmin,dlntdr_loc(i),1)
     call cub_spline(rmin_exp,dlnndr_exp(i,:),EXPRO_n_exp,rmin,dlnndr_loc(i),1)
     call cub_spline(rmin_exp,sdlntdr_exp(i,:),EXPRO_n_exp,rmin,sdlntdr_loc(i),1)
     call cub_spline(rmin_exp,sdlnndr_exp(i,:),EXPRO_n_exp,rmin,sdlnndr_loc(i),1)
  enddo

  if (numeq_flag == 1) then

     if (EXPRO_nfourier <= 0) then
        !call cgyro_error('Geometry coefficients missing')
        return
     endif

     geo_ny_loc = EXPRO_nfourier
     allocate(geo_yin_exp(8,0:geo_ny_loc,EXPRO_n_exp))
     if(allocated(geo_yin_loc)) deallocate(geo_yin_loc)
     allocate(geo_yin_loc(8,0:geo_ny_loc))
     geo_yin_exp(1:4,:,:) = EXPRO_geo(:,:,:)/a_meters
     geo_yin_exp(5:8,:,:) = EXPRO_dgeo(:,:,:)

     do i=1,8
        do j=0,geo_ny_loc
           call cub_spline(rmin_exp,geo_yin_exp(i,j,:),EXPRO_n_exp,rmin, &
                geo_yin_loc(i,j),1)
        enddo
     enddo

  endif

  call EXPRO_palloc(comm,path,0)
  call EXPRO_locsim_alloc(0)

end subroutine EXPRO_locsim_profiles
