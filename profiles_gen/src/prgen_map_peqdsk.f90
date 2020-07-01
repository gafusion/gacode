!------------------------------------------------------------
! prgen_map_peqdsk.f90
!
! PURPOSE:
!  Map native peqdsk data onto input.profiles standard.  
!------------------------------------------------------------

subroutine prgen_map_peqdsk

  use prgen_globals
  use expro

  implicit none

  integer :: i
  real, dimension(nx) :: ni_d
  real, dimension(3,nx) :: ni_imp
  real, dimension(nx) :: ni_b
  real, dimension(nx) :: z_eff

  ni_d(:) = 10*peqdsk_ni(:)
  if (peqdsk_fmt == 0) then
     ! old p-file, assume missing density is carbon
     ni_imp(1,:) = 10*(peqdsk_ne(:)-peqdsk_ni(:)-peqdsk_nb(:))/6.0
     peqdsk_nimp  = 1
     peqdsk_z(1)  = 1.0
     peqdsk_m(1)  = 2.0
     peqdsk_z(2)  = 6.0
     peqdsk_m(2)  = 12.0
     if (peqdsk_nbeams == 1) then
        peqdsk_z(3)  = 1.0
        peqdsk_m(3)  = 2.0
     endif
  else
     ! new p-file
     do i=1,peqdsk_nimp
        ni_imp(i,:) = 10*peqdsk_nz(i,:)
     enddo
  endif
  ni_b(:) = 10*peqdsk_nb(:)

  ! Compute Z_eff
  z_eff(:) = peqdsk_z(1)**2 * ni_d(:)
  do i=1,peqdsk_nimp
     z_eff(:) = z_eff(:)+peqdsk_z(1+i)**2 * ni_imp(i,:)
  enddo
  z_eff(:) = z_eff(:)+peqdsk_z(1+peqdsk_nimp+1)**2 * ni_b(:)
  z_eff(:) = z_eff(:)/(10*peqdsk_ne(:))

  !---------------------------------------------------------
  ! Map profile data onto single array:
  !
  peqdsk_nion = 1+peqdsk_nimp+peqdsk_nbeams

  expro_n_exp = nx
  expro_n_ion = peqdsk_nion
  call expro_init(1)
  !
  expro_rho(:)   = rho(:)
  expro_rmin(:)  = rmin(:)
  expro_rmaj(:)  = rmaj(:)
  expro_te(:)    = peqdsk_te(:)
  expro_ne(:)    = peqdsk_ne(:)*10
  expro_z_eff(:) = z_eff(:)
  expro_ptot(:)  = p_tot(:)      
  ! COORDINATES: -ipccw accounts for DIII-D toroidal angle convention
  ! (wrt Ip direction)
  expro_w0(:) = -ipccw*1e3*peqdsk_omgeb(:) 

  expro_mass = peqdsk_m(1:expro_n_ion)
  expro_z    = peqdsk_z(1:expro_n_ion)
  do i=1,expro_n_ion
     call prgen_ion_name(nint(expro_mass(i)),nint(expro_z(i)),expro_name(i))         
  enddo

  ! ni, nc, nb -and- ti, tc, tb
  expro_ni(1,:) = ni_d(:)
  expro_ti(1,:) = peqdsk_ti(:)
  expro_type(1) = type_therm
  do i=1,peqdsk_nimp
     expro_ni(1+i,:) = ni_imp(i,:)
     expro_ti(1+i,:) = peqdsk_ti(:)
     expro_type(1+1) = type_therm
  enddo
  if (peqdsk_nbeams == 1) then
     i = 1+peqdsk_nimp+1
     expro_ni(i,:) = ni_b(:)
     ! JC: need to check for zero density here?
     expro_ti(i,:) = peqdsk_pb(:)/(peqdsk_nb(:)*10)/1.602
     expro_type(i) = type_fast
  endif

  ! vphi
  ! COORDINATES: negative sign accounts for DIII-D toroidal angle convention
  expro_vtor(:,:) = 0.0
  expro_vtor(2,:) = -1e3*peqdsk_omegat(:)*(rmaj(:)+rmin(:))

  ! vpol
  expro_vpol(:,:) = 0.0

  !---------------------------------------------------------
  ! Read the cer file and overlay
  !
  if (file_cer /= "null") then
     call prgen_read_cer
     expro_w0 = omega0(:)
     do i=1,peqdsk_nimp
        if (peqdsk_m(i+1) .le. 12.0+epsilon(0.) .and. &
             peqdsk_m(i+1) .ge. 12.0-epsilon(0.)) then
           expro_vtor(i+1,:) = vtorc_exp(:)
           expro_vpol(i+1,:) = vpolc_exp(:)
           exit
        endif
     enddo
  endif
  !---------------------------------------------------------
  
end subroutine prgen_map_peqdsk
