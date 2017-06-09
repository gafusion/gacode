!----------------------------------------------------------------
! neo_experimental_profiles.f90
!
! PURPOSE:
!  Read experimental profiles. and generate additional profiles.
!
! NOTES: 
!  See http://fusion.gat.com/theory/INPUT_profiles for a 
!  complete description of the INPUT_profiles file structure.
! 
!  The files read by this routine are *generated* from 
!  INPUT_profiles by scripts. 
!----------------------------------------------------------------

subroutine neo_experimental_profiles

  use neo_globals
  use neo_profile_exp
  use neo_allocate_profile
  use EXPRO_interface

  !---------------------------------------------------
  implicit none

  integer :: i, is, i_ion

  !--------------------------------------------------------------
  ! use EXPRO routines to read data:
  !
  call EXPRO_alloc(path,1)
  EXPRO_ctrl_quasineutral_flag = 1  ! quasi-neutrality density flag

  if (profile_equilibrium_model == 2) then
     EXPRO_ctrl_numeq_flag = 1
  else
     EXPRO_ctrl_numeq_flag = 0
  endif

  ! Determine if electrons are to be included in the neo 
  ! simulation.  Electron profiles are read even if not 
  ! to be included in the simulation (needed to re-scale 
  ! ion density/temp is not quasi-neutral).

  if (adiabatic_ele_model == 1) then
     n_species_exp = n_species + 1
  else
     if (Z(n_species) > 0.0) then
        call neo_error('ERROR: (NEO) For exp. profiles, electron species must be n_species')
        return
     endif
     n_species_exp = n_species
  endif

  ! charge of ion species
  EXPRO_ctrl_z(:) = 0.0
  EXPRO_ctrl_n_ion = 0
  if(adiabatic_ele_model == 1) then
     do is=1,n_species
        EXPRO_ctrl_z(is) = 1.0 * Z(is)
        EXPRO_ctrl_n_ion = EXPRO_ctrl_n_ion + 1
     enddo
  else
     do is=1,n_species-1
        EXPRO_ctrl_z(is) = 1.0 * Z(is)
        EXPRO_ctrl_n_ion = EXPRO_ctrl_n_ion + 1
     enddo
  endif

  call EXPRO_read

  n_grid_exp = EXPRO_n_exp
  call PROFILE_EXP_alloc(1)

  sign_bunit = EXPRO_signb
  sign_q     = EXPRO_signq

  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Transfer data from read vector to individual arrays:
  !
  rhoN_torflux_a      = EXPRO_arho
  rhoN_torflux_exp(:) = EXPRO_rho(:)
  if(abs(EXPRO_polflux(n_grid_exp)) < epsilon(0.)) then
     psiN_polflux_exp(:) = 0.0
     psiN_polflux_a      = 0.0
  else
     psiN_polflux_exp(:) = EXPRO_polflux(:) &
          / EXPRO_polflux(n_grid_exp)
     psiN_polflux_a     = EXPRO_polflux(n_grid_exp)
  endif
  rmin_exp(:)         = EXPRO_rmin(:)
  rmaj_exp(:)         = EXPRO_rmaj(:)
  q_exp(:)            = EXPRO_q(:)
  kappa_exp(:)        = EXPRO_kappa(:)
  delta_exp(:)        = EXPRO_delta(:)  * profile_delta_scale
  zeta_exp(:)         = EXPRO_zeta(:)   * profile_zeta_scale
  zmag_exp(:)         = EXPRO_zmag(:)   * profile_zmag_scale
  shift_p(:)          = EXPRO_drmaj(:)
  s_kappa_p(:)        = EXPRO_skappa(:)
  s_delta_p(:)        = EXPRO_sdelta(:) * profile_delta_scale
  s_zeta_p(:)         = EXPRO_szeta(:)  * profile_zeta_scale
  s_zmag_p(:)         = EXPRO_dzmag(:)  * profile_zmag_scale
  shear_p(:)           = EXPRO_s(:)

  r_p(:)    = rmin_exp(:)/rmin_exp(n_grid_exp)
  rmaj_p(:) = rmaj_exp(:)/rmin_exp(n_grid_exp)
  zmag_p(:) = zmag_exp(:)/rmin_exp(n_grid_exp)

  ! Minor radius, a, in meters:
  a_meters = rmin_exp(n_grid_exp)

  if(silent_flag == 0 .and. i_proc == 0 .and. &
       abs(profile_delta_scale-1.0) > epsilon(0.) ) then
     open(unit=io_neoout,file=trim(path)//runfile_neoout,&
          status='old',position='append')
     write(io_neoout,*) 'Delta and S_Delta are re-scaled'
     close(io_neoout)
  end if
  if(silent_flag == 0 .and. i_proc == 0 .and. &
       abs(profile_zeta_scale-1.0) > epsilon(0.)) then
     open(unit=io_neoout,file=trim(path)//runfile_neoout,&
          status='old',position='append')
     write(io_neoout,*) 'Zeta and S_Zeta are re-scaled'
     close(io_neoout)
  end if
  if(silent_flag == 0 .and. i_proc == 0 .and. &
       abs(profile_zmag_scale-1.0) > epsilon(0.)) then
     open(unit=io_neoout,file=trim(path)//runfile_neoout,&
          status='old',position='append')
     write(io_neoout,*) 'Zmag and S_Zmag are re-scaled'
     close(io_neoout)
  end if

  if (profile_equilibrium_model == 2) then
     if(EXPRO_nfourier <= 0) then
        call neo_error('ERROR: (NEO) Geometry coefficients missing')
        return
     endif
     geo_numeq_flag = 1
     deallocate(geo_yin_exp)
     deallocate(geo_yin)
     geo_ny = EXPRO_nfourier
     allocate(geo_yin_exp(8,0:geo_ny,n_grid_exp))
     allocate(geo_yin(8,0:geo_ny,n_radial))
     geo_yin_exp(1:4,:,:) = EXPRO_geo(:,:,:)/rmin_exp(n_grid_exp)
     geo_yin_exp(5:8,:,:) = EXPRO_dgeo(:,:,:)
  endif

  te_ade_exp(:)      = EXPRO_te(:)
  ne_ade_exp(:)      = EXPRO_ne(:)
  dlntdre_ade_exp(:) = EXPRO_dlntedr(:) * a_meters &
       * profile_dlntdr_scale(n_species_exp)
  dlnndre_ade_exp(:) = EXPRO_dlnnedr(:) * a_meters &
       * profile_dlnndr_scale(n_species_exp)

  tem_exp(n_species_exp,:)  = EXPRO_te(:)
  dlntdr_p(n_species_exp,:) = EXPRO_dlntedr(:) * a_meters &
       * profile_dlntdr_scale(n_species_exp)
  den_exp(n_species_exp,:)  = EXPRO_ne(:)
  dlnndr_p(n_species_exp,:) = EXPRO_dlnnedr(:) * a_meters &
       * profile_dlnndr_scale(n_species_exp)

  do i_ion=1,n_species_exp-1
     ! ion temps should be equal, but not enforced 
     tem_exp(i_ion,:)  = EXPRO_ti(i_ion,:)
     dlntdr_p(i_ion,:) = EXPRO_dlntidr(i_ion,:) * a_meters &
          * profile_dlntdr_scale(i_ion)

     ! first species density is re-set by quasi-neutrality
     if(i_ion == 1) then
        den_exp(i_ion,:)  = EXPRO_ni_new(:)
        dlnndr_p(i_ion,:) = EXPRO_dlnnidr_new(:) * a_meters &
             * profile_dlnndr_scale(i_ion)
     else
        den_exp(i_ion,:)  = EXPRO_ni(i_ion,:)
        dlnndr_p(i_ion,:) = EXPRO_dlnnidr(i_ion,:) * a_meters &
             * profile_dlnndr_scale(i_ion)

     endif
  enddo

  ! If desired, reset electron density gradient to ensure quasi-neutrality
  if (profile_dlnndr_scale(n_species_exp) < 0.0) then
     dlnndr_p(n_species_exp,:) = 0.0
     do i_ion=1,n_species_exp-1
        dlnndr_p(n_species_exp,:) = dlnndr_p(n_species_exp,:) + &
             EXPRO_ctrl_z(i_ion)*den_exp(i_ion,:)/den_exp(n_species_exp,:) &
             * dlnndr_p(i_ion,:)
     enddo
  endif

  ! Sanity check for temperatures
  do i=1,n_species_exp
     if (minval(den_exp(i,:)) <= 0.0) then
        call neo_error('ERROR: (NEO) Nonpositive in exp. density profile')
        return
     endif
  enddo

  ! Bunit 

  b_unit_p(:)  = EXPRO_bunit(:)

  ! Rotation and Er

  omega_rot_p(:)       = EXPRO_w0(:) 
  omega_rot_deriv_p(:) = EXPRO_w0p(:) 
  dphi0dr_p(:) = -omega_rot_p(:) * rmin_exp(:) * b_unit_p(:) / q_exp(:)

  call EXPRO_alloc(path,0)

end subroutine neo_experimental_profiles
