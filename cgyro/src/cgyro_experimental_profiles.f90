!----------------------------------------------------------------
! cgyro_experimental_profiles.f90
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

subroutine cgyro_experimental_profiles

  use cgyro_globals
  use cgyro_io
  use cgyro_experimental_globals
  use EXPRO_interface

  implicit none

  integer :: i, is, i_ion

  !--------------------------------------------------------------
  ! use EXPRO routines to read data:
  !
  call EXPRO_palloc(CGYRO_COMM_WORLD,path,1)
  EXPRO_ctrl_quasineutral_flag = 1  ! quasi-neutrality density flag

  if (equilibrium_model == 3) then
     EXPRO_ctrl_numeq_flag = 1
  else
     EXPRO_ctrl_numeq_flag = 0
  endif

  ! Determine if electrons are to be included in the 
  ! simulation.  Electron profiles are read even if not 
  ! to be included in the simulation (needed to re-scale 
  ! ion density/temp is not quasi-neutral).

  if (ae_flag == 0) then
     n_species_exp = n_species
     if (z(n_species) > 0.0) then
        call cgyro_error('For exp. profiles, electron species must be n_species')
        return
     endif
  else
     n_species_exp = n_species + 1
  endif

  ! charge of ion species
  EXPRO_ctrl_z(:) = 0.0
  EXPRO_ctrl_n_ion = 0
  if(ae_flag == 1) then
     do is=1,n_species
        EXPRO_ctrl_z(is) = z(is)
        EXPRO_ctrl_n_ion = EXPRO_ctrl_n_ion + 1
     enddo
  else
     do is=1,n_species-1
        EXPRO_ctrl_z(is) = z(is)
        EXPRO_ctrl_n_ion = EXPRO_ctrl_n_ion + 1
     enddo
  endif

  call EXPRO_pread

  n_grid_exp = EXPRO_n_exp
  call cgyro_experimental_alloc(1)

  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! Transfer data from read vector to individual arrays:
  !

  btccw = -EXPRO_signb
  ipccw = -EXPRO_signq*EXPRO_signb

  rmin_exp(:)    = EXPRO_rmin(:)
  rmaj_exp(:)    = EXPRO_rmaj(:)
  q_exp(:)       = EXPRO_q(:)
  s_exp(:)       = EXPRO_s(:)
  shift_exp(:)   = EXPRO_drmaj(:)
  kappa_exp(:)   = EXPRO_kappa(:)
  s_kappa_exp(:) = EXPRO_skappa(:)
  delta_exp(:)   = EXPRO_delta(:) 
  s_delta_exp(:) = EXPRO_sdelta(:) 
  zeta_exp(:)    = EXPRO_zeta(:)
  s_zeta_exp(:)  = EXPRO_szeta(:)
  if (udsymmetry_flag == 1 .and. equilibrium_model /= 3) then
     zmag_exp(:)    = 0.0   
     dzmag_exp(:)   = 0.0
  else
     zmag_exp(:)    = EXPRO_zmag(:)   
     dzmag_exp(:)   = EXPRO_dzmag(:)
  endif
  
  ! Minor radius, a, in meters:
  a_meters = rmin_exp(n_grid_exp)

  rmin_exp(:) = rmin_exp(:)/a_meters
  rmaj_exp(:) = rmaj_exp(:)/a_meters
  zmag_exp(:) = zmag_exp(:)/a_meters

  if (equilibrium_model == 3) then
     if (EXPRO_nfourier <= 0) then
        call cgyro_error('Geometry coefficients missing')
        return
     endif
     geo_numeq_flag = 1
     deallocate(geo_yin_exp)
     deallocate(geo_yin)
     geo_ny = EXPRO_nfourier
     allocate(geo_yin_exp(8,0:geo_ny,n_grid_exp))
     allocate(geo_yin(8,0:geo_ny))
     geo_yin_exp(1:4,:,:) = EXPRO_geo(:,:,:)/a_meters
     geo_yin_exp(5:8,:,:) = EXPRO_dgeo(:,:,:)
  endif

  te_ade_exp(:)      = EXPRO_te(:)
  ne_ade_exp(:)      = EXPRO_ne(:)
  dlntdre_ade_exp(:) = EXPRO_dlntedr(:) * a_meters
  dlnndre_ade_exp(:) = EXPRO_dlnnedr(:) * a_meters

  temp_exp(n_species_exp,:)   = EXPRO_te(:)
  dlntdr_exp(n_species_exp,:) = EXPRO_dlntedr(:) * a_meters 
  sdlntdr_exp(n_species_exp,:) = EXPRO_sdlntedr(:) * a_meters**2 

  dens_exp(n_species_exp,:)   = EXPRO_ne(:)
  dlnndr_exp(n_species_exp,:) = EXPRO_dlnnedr(:) * a_meters 
  sdlnndr_exp(n_species_exp,:) = EXPRO_sdlnnedr(:) * a_meters**2 

  do i_ion=1,n_species_exp-1
     ! ion temps should be equal, but not enforced 
     temp_exp(i_ion,:)   = EXPRO_ti(i_ion,:)
     dlntdr_exp(i_ion,:) = EXPRO_dlntidr(i_ion,:) * a_meters 
     sdlntdr_exp(i_ion,:) = EXPRO_sdlntidr(i_ion,:) * a_meters**2 

     ! first species density is re-set by quasi-neutrality [JC: do we want this?]
     if(quasineutral_flag == 1 .and. i_ion == 1) then
        dens_exp(i_ion,:)   = EXPRO_ni_new(:)
        dlnndr_exp(i_ion,:) = EXPRO_dlnnidr_new(:) * a_meters
        sdlnndr_exp(i_ion,:) = EXPRO_sdlnnidr_new(:) * a_meters**2
     else
        dens_exp(i_ion,:)   = EXPRO_ni(i_ion,:)
        dlnndr_exp(i_ion,:) = EXPRO_dlnnidr(i_ion,:) * a_meters
        sdlnndr_exp(i_ion,:) = EXPRO_sdlnnidr(i_ion,:) * a_meters**2
     endif
  enddo

  ! Sanity check for densities
  do i=1,n_species_exp
     if (minval(dens_exp(i,:)) <= 0.0) then
        call cgyro_error('Nonpositive in exp. density profile')
        return
     endif
  enddo

  ! Bunit 

  b_unit_exp(:) = EXPRO_bunit(:)

  ! Rotation

  gamma_e_exp(:) = -EXPRO_w0p(:) * (a_meters * rmin_exp(:)) / q_exp(:)

  gamma_p_exp(:) = -EXPRO_w0p(:) * (a_meters * rmaj_exp(:))

  mach_exp(:)    = EXPRO_w0(:) * (a_meters * rmaj_exp(:))

  rhos_exp(:)    = EXPRO_rhos(:)

  call EXPRO_palloc(CGYRO_COMM_WORLD,path,0)

end subroutine cgyro_experimental_profiles
