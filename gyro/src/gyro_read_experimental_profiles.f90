!----------------------------------------------------------------
! gyro_read_experimental_profiles.f90 [caller gyro_profile_init]
!
! PURPOSE:
!  Read experimental profiles using EXPRO library.
!
! NOTES: 
!  See http://fusion.gat.com/theory/input.profiles for a 
!  complete description of the input.profiles file structure.
!----------------------------------------------------------------

subroutine gyro_read_experimental_profiles

  use gyro_globals
  use gyro_profile_exp
  use EXPRO_interface

  !---------------------------------------------------
  implicit none
  !
  real :: p_total
  !---------------------------------------------------

  !---------------------------------------------------------------------
  ! Read experimental profiles using EXPRO library.
  !
  EXPRO_ctrl_quasineutral_flag = density_method-1
  EXPRO_ctrl_numeq_flag = num_equil_flag
  EXPRO_ctrl_n_ion = n_spec-1

  call EXPRO_palloc(GYRO_COMM_WORLD,path,1) 
  call EXPRO_pread
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  n_grid_exp = EXPRO_n_exp
  if (EXPRO_nfourier < 0) then
     n_fourier_geo = 0
     if (num_equil_flag == 1) then
        call catch_error('ERROR: (GYRO) Geometry coefficients missing.')
     endif
  else
     n_fourier_geo = EXPRO_nfourier
  endif
  !
  call gyro_alloc_profile_exp
  !
  ! Transfer data from read arrays to individual arrays:
  !
  bt_exp   = EXPRO_b_ref
  arho_exp = EXPRO_arho
  !
  btccw = -EXPRO_signb
  ipccw = -EXPRO_signq*EXPRO_signb
  !
  call send_message('INFO: (GYRO) Resetting IPCCW and BTCCW based on input.profiles.')
  !
  rhogrid_exp(:)    = EXPRO_rho(:)
  rmin_exp(:)       = EXPRO_rmin(:)
  rmaj_exp(:)       = EXPRO_rmaj(:)
  zmag_exp(:)       = EXPRO_zmag(:)
  ptot_exp(:)       = EXPRO_ptot(:)
  q_exp(:)          = EXPRO_q(:)
  kappa_exp(:)      = EXPRO_kappa(:)
  delta_exp(:)      = EXPRO_delta(:)
  zeta_exp(:)       = EXPRO_zeta(:)
  tem_exp(n_spec,:) = EXPRO_te(:)
  den_exp(n_spec,:) = EXPRO_ne(:)
  z_eff_exp(:)      = EXPRO_z_eff(:)
  w0_exp(:)         = EXPRO_w0(:)
  w0p_exp(:)        = EXPRO_w0p(:)
  gamma_e_exp(:)    = EXPRO_gamma_e(:)
  gamma_p_exp(:)    = EXPRO_gamma_p(:)
  mach_exp(:)       = EXPRO_mach(:)
  !
  do is=1,n_spec-1
     den_exp(is,:)  = EXPRO_ni(is,:)*n_vec(is)
     tem_exp(is,:)  = EXPRO_ti(is,:)
     if (n_vec(is) /= 1.0) then
        call send_message_real(&
             'INFO: ni'//achar(is-1+iachar("1"))//' rescaled by: ',&
             n_vec(is))
     endif
  enddo

  !--------------------------------------------------------
  ! Minor radius, a, in meters:
  !
  a_meters = rmin_exp(n_grid_exp)
  !--------------------------------------------------------

  ! Fill in computed profile quantities:

  b_unit_p(:)    = EXPRO_bunit(:)
  shat_p(:)      = EXPRO_s(:)
  s_kappa_p(:)   = EXPRO_skappa(:)
  drmaj_p(:)     = EXPRO_drmaj(:)
  s_delta_p(:)   = EXPRO_sdelta(:)
  s_zeta_p(:)    = EXPRO_szeta(:)
  dzmag_p(:)     = EXPRO_dzmag(:)

  ! The case is=1 may be reset below.

  do is=1,n_spec-1
     dlntdr_p(is,:) = a_meters*EXPRO_dlntidr(is,:)
     dlnndr_p(is,:) = a_meters*EXPRO_dlnnidr(is,:)
  enddo
  dlntdr_p(n_spec,:) = a_meters*EXPRO_dlntedr(:)
  dlnndr_p(n_spec,:) = a_meters*EXPRO_dlnnedr(:)
  dlnptotdr_p(:)     = a_meters*EXPRO_dlnptotdr(:)
  !
  dlnndr_p(:,1)  = 0.0
  dlntdr_p(:,1)  = 0.0
  dlnptotdr_p(1) = 0.0
  !----------------------------------------------------------------------------

  ! Retain up-down asymmetry from elevation in the case of Miller shape
  if (num_equil_flag == 0 .and. udsymmetry_flag == 1) then
     zmag_exp(:) = 0.0
     dzmag_p(:)  = 0.0
  endif

  ! Fill in general geometry parameters if they exist
  if (num_equil_flag == 1) then

     ! Fill in geometry arrays with original data from 
     ! EXPRO routine.  Note that EXPRO_geo is in m, 
     ! and EXPRO_dgeo is dimensionless.  All geo_p are 
     ! dimensionless:

     geo_p(1:4,:,:) = EXPRO_geo(:,:,:)/a_meters
     geo_p(5:8,:,:) = EXPRO_dgeo(:,:,:)

  endif

  ! rhos_deuterium/a (dimensionless)
  rhosda_p(:) = EXPRO_rhos(:)/a_meters

  ! cs_deuterium/a (1/s)
  csda_p(:)   = EXPRO_cs(:)/a_meters

  !--------------------------------------------------------------
  ! Management of multiple ion densities:
  !
  select case (density_method)

  case (1)

     call send_message('INFO: (GYRO) Taking densities directly from input.profiles')

  case (2) 

     call send_message('INFO: (GYRO) Offsetting main ions to force sum(ni) = ne')
     den_exp(1,:) = EXPRO_ni_new(:)
     dlnndr_p(1,:) = a_meters*EXPRO_dlnnidr_new(:)

  case default

     call catch_error('ERROR: (GYRO) density_method')

  end select
  !--------------------------------------------------------------

  !-------------------------------------------------------------------
  ! Sanity check for densities
  !
  do is=1,n_spec
     if (minval(den_exp(is,:)) <= 0.0) then
        print *,is,den_exp(is,:)
        call catch_error(&
             'ERROR: (GYRO) Nonpositive in density profile '//achar(is-1+iachar("1")))  
     endif
  enddo

  do i_exp=1,n_grid_exp

     r_p(i_exp) = rmin_exp(i_exp)/rmin_exp(n_grid_exp)
     q_p(i_exp) = q_exp(i_exp)

     if (q_scale /= 1.0) then
        q_p(i_exp) = q_exp(i_exp)*q_scale
     endif

     !-------------------------------------------------------------
     ! Specification of beta_unit (beta in terms of b_unit).
     !
     ! beta_unit is defined with respect to b_unit, as required 
     ! by GEO.  For comparison with DIII-D, beta should be 
     ! defined with respect to bt_exp.
     !
     ! den_s  -> 1/m^3
     ! tem_s  -> keV
     ! b_unit -> T
     !
     ! beta calculation in CGS:
     !
     !         8*pi ( n[1e19/m^3]*1e-6*1e19 )( T[keV]*1.6022*1e-9 )
     ! beta = ------------------------------------------------------
     !                           ( 1e4*B[T] )^2
     !
     !      = 4.027e-3 n[1e19/m^3]*T[keV]/B[T]^2
     !
     p_total = sum(den_exp(:,i_exp)*tem_exp(:,i_exp))
     !
     beta_unit_p(i_exp) = 4.027e-3*p_total/b_unit_p(i_exp)**2
     !
     ! Here, ptot_exp is in Pa:
     !
     !  beta = 8*pi ( 10 ptot[Pa] )/( 1e4*B[T] )^2 = 2.513e-6 ptot[Pa]/B[T]^2
     !
     beta_unit_ptot_p(i_exp) = 2.513e-6*ptot_exp(i_exp)/b_unit_p(i_exp)**2
     !-------------------------------------------------------------

  enddo

  if (ampere_scale /= 1.0) then
     call send_message_real('INFO: (GYRO) betae in Ampere Eq. scaled by ',ampere_scale)
  endif

  ! Deallocate EXPRO arrays
  call EXPRO_palloc(GYRO_COMM_WORLD,path,0)

  if (debug_flag == 1 .and. i_proc == 0) then
     print *,'[gyro_read_experimental_profiles done]'
  endif
  !---------------------------------------------------------------------

end subroutine gyro_read_experimental_profiles
