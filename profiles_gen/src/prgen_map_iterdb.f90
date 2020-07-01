!------------------------------------------------------------
! prgen_map_iterdb.f90
!
! PURPOSE:
!  Map native iterdb data onto input.profiles standard.
!
! NOTES:
!  - This routine is common to both text and NetCDF formats.
!  - See map_plasmastate.f90 for analogous routine for
!    plasmastate data.
!------------------------------------------------------------

subroutine prgen_map_iterdb

  use prgen_globals
  use expro

  implicit none

  integer :: i
  integer :: n0

  !----------------------------------------------------------------------
  ! shot number and time slice (ms)
  !
  expro_shot = onetwo_ishot
  expro_time = int(onetwo_time*1e3)

  !----------------------------------------------------------------------
  ! Construct ion densities and temperatures, manage naming and numbering
  !
  onetwo_enion_vec(:,:) = 0.0
  onetwo_tion_vec(:,:) = 0.0
  do i=1,onetwo_nion
     ! ni
     onetwo_enion_vec(i,:) = onetwo_enion(:,i)*1e-19
     ! Ti
     onetwo_tion_vec(i,:) = onetwo_ti(:)
  enddo

  ! Primary ions
  onetwo_ion_name(1:onetwo_nprim) = onetwo_namep(1:onetwo_nprim)
  n0 = onetwo_nprim

  ! Impurities
  onetwo_ion_name(1+n0:onetwo_nimp+n0) = onetwo_namei(1:onetwo_nimp)
  n0 = n0+onetwo_nimp

  ! Beams
  do i=1,onetwo_nbion
     ! Only accept beams if density is everywhere positive
     if (minval(onetwo_enbeam(:,i)) > epsilon(0.0)) then

        print '(a)',"INFO: (prgen_map_iterdb) Found fast ion beam species"
        print '(a)',"INFO: (prgen_map_iterdb) Modifying fast ion beam temperature to satisfy beam pressure"

        n0 = n0+1

        onetwo_enion_vec(n0,:) = onetwo_enbeam(:,i)*1e-19

        ! Ti: T[keV] = (p/n)[J]/1.6022e-16[J/keV]
        onetwo_tion_vec(n0,:) = onetwo_pressb(:,i)/onetwo_enbeam(:,i)/1.6022e-16

        onetwo_ion_name(n0) = onetwo_nameb(i)

     endif
  enddo

  ! Fast alphas
  if (minval(onetwo_enalp(:)) > epsilon(0.0)) then

     print '(a)',"INFO: (prgen_map_iterdb) Found fast alpha species"
     print '(a)',"INFO: (prgen_map_iterdb) Modifying fast alpha temperature to satisfy total pressure"

     n0 = n0+1

     onetwo_enion_vec(n0,:) = onetwo_enalp(:)*1e-19

     onetwo_tion_vec(n0,:) = (p_tot(:)-&
          (sum(onetwo_enion_vec(1:n0-1,:)*&
          onetwo_tion_vec(1:n0-1,:),dim=1)*1e19+&
          onetwo_ene(:)*onetwo_te(:))*1.6022e-16)/(onetwo_enalp(:))/1.6022e-16

     onetwo_ion_name(n0) = 'he'

  endif

  if (n0 > 10) then
     print '(a)',"ERROR: (prgen_map_iterdb) Too many ions; report to GACODE developers."
     stop
  endif

  !---------------------------------------------------------
  ! Map profile data into expro interface variables
  !
  expro_n_exp = nx
  expro_n_ion = n0
  call expro_init(1)
  !
  expro_rho  = rho
  expro_rmin = rmin(:)
  expro_rmaj = rmaj(:)
  expro_te = onetwo_te(:)
  expro_ne = onetwo_ene(:)*1e-19
  expro_z_eff = onetwo_zeff(:)
  expro_ptot = p_tot(:) ! Total pressure

  do i=1,expro_n_ion
     expro_ni(i,:) = onetwo_enion_vec(i,:)
     expro_ti(i,:) = onetwo_tion_vec(i,:)
  enddo

  ! velocities
  expro_vpol = 0.0
  expro_vtor = 0.0

  ! Look for carbon as first impurity, and insert toroidal velocity at theta=0

  if (trim(onetwo_namei(1)) == 'c') then
     ! COORDINATES: -ipccw accounts for DIII-D toroidal angle convention
     vtorc_exp(:) = -ipccw*onetwo_angrot(:)*(rmaj(:)+rmin(:))
  else
     vtorc_exp(:) = 0.0
  endif

  ! Insert carbon toroidal velocity
  do i=1,expro_n_ion
     if (i == onetwo_nprim+1) then
        expro_vtor(i,:) = vtorc_exp(:)
     endif
  enddo

  ! Use angrot as an initial approximation for omega0
  expro_w0 = -ipccw*onetwo_angrot(:)
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Read the cer file and overlay
  !
  if (file_cer /= "null") then
     call prgen_read_cer
     expro_w0 = omega0(:)
     do i=1,expro_n_ion
        if (i == onetwo_nprim+1) then
           expro_vtor(i,:) = vtorc_exp(:)
           expro_vpol(i,:) = vpolc_exp(:)
        endif
     enddo
  endif
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Ion name association
  !
  print '(a)','INFO: (prgen_map_iterdb) Found these ion species:'
  do i=1,expro_n_ion
     call onetwo_ion_zmass(onetwo_ion_name(i),onetwo_z(i),onetwo_m(i))
     print '(t6,i2,1x,a)',i,trim(onetwo_ion_name(i))
  enddo
  !
  do i=1,expro_n_ion
     expro_mass(i) = onetwo_m(i)             
     expro_z(i)    = onetwo_z(i)             
     call prgen_ion_name(nint(expro_mass(i)),nint(expro_z(i)),expro_name(i))         
     if (i > onetwo_nion) then
        expro_type(i) = type_fast
     else
        expro_type(i) = type_therm
     endif
  enddo
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Power densities (q*) : source > 0, sink < 0
  !
  ! Ohmic
  expro_qohme = 1e-6*onetwo_qohm
  ! NBI
  expro_qbeame = 1e-6*onetwo_qbeame
  expro_qbeami = 1e-6*onetwo_qbeami
  ! RF
  expro_qrfe = 1e-6*onetwo_qrfe
  expro_qrfi = 1e-6*onetwo_qrfi
  ! Fusion
  expro_qfuse = 1e-6*onetwo_qfuse
  expro_qfusi = 1e-6*onetwo_qfusi
  ! Radiated
  expro_qbrem = 1e-6*(onetwo_qrad-onetwo_qione-onetwo_qsync)
  expro_qsync = 1e-6*onetwo_qsync
  expro_qline = 0.0
  ! electron-to-ion exchange (positive for Te > Ti)
  expro_qei = -1e-6*onetwo_qdelt
  ! neutrals (ion=ionization and recombination, cx=charge exchange) 
  expro_qione = 1e-6*onetwo_qione
  expro_qioni = 1e-6*onetwo_qioni
  expro_qcxi  = 1e-6*onetwo_qcx
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Particle/momentum sources
  !
  ! particles
  expro_qpar_beam = onetwo_sbeame
  expro_qpar_wall = sion_d
  !
  ! Torque (TAM flow) (nt-m)
  expro_qmom = -ipccw*onetwo_storqueb
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Current densities
  !
  expro_johm = -ipccw*johm*1e-6
  expro_jbs = -ipccw*jbs*1e-6
  expro_jnb = -ipccw*jnb*1e-6
  expro_jrf = -ipccw*jrf*1e-6
  !---------------------------------------------------------

end subroutine prgen_map_iterdb

!------------------------------------------------------------------
! Routine to perform Volume integration
!------------------------------------------------------------------

subroutine onetwo_ion_zmass(iname,z,m)

  implicit none
  
  character(len=2) :: iname
  integer :: z
  real :: m
  
  select case (iname)
  case ('h')
     z = 1
     m = 1.0
  case ('d')
     z = 1
     m = 2.0
  case ('t')
     z = 1
     m = 3.0
  case ('he')
     z = 2
     m = 4.0
  case ('c')
     z = 6
     m = 12.0
  case ('ar')
     z = 18
     m = 40.0
  case default
     z = 0
     m = 0.0
  end select

end subroutine onetwo_ion_zmass
