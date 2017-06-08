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
  use EXPRO_interface

  implicit none

  integer :: i
  integer :: ip
  integer :: n0

  do i=1,nx
     rho(i) = (i-1)/(nx-1.0)
  enddo

  call volint(1e-6*onetwo_qbeame,pow_e_nb)
  call volint(1e-6*onetwo_qrfe,pow_e_rf)
  call volint(1e-6*onetwo_qohm,pow_e_ohm)
  call volint(1e-6*onetwo_qrad,pow_e_rad)
  call volint(1e-6*onetwo_qione,powe_ion_exp)
  call volint(1e-6*onetwo_dpedt,powe_wdot_exp)
  call volint(1e-6*onetwo_qfuse,pow_e_fus)
  call volint(1e-6*onetwo_qfusi,pow_i_fus)
  call volint(1e-6*onetwo_qbeami,pow_i_nb)
  call volint(1e-6*onetwo_qrfi,pow_i_rf)
  call volint(1e-6*onetwo_qioni,powi_ion_exp)
  call volint(1e-6*onetwo_qcx,powi_cx_exp)
  call volint(1e-6*onetwo_dpidt,powi_wdot_exp)
  ! Note sign change to keep pow_ei_exp positive for te>ti
  call volint(-1e-6*onetwo_qdelt,pow_ei)

  ! Wall ion flow (scale is undefined) (MW/KeV = Kamp)
  call volint(kevdsecpmw*sion_d,flow_wall_exp)

  ! Beam ion flow
  call volint(kevdsecpmw*(onetwo_sbeam+sbcx_d),flow_beam)

  ! Torque (TAM flow) (nt-m)
  call volint(onetwo_storqueb,flow_mom)
  ! COORDINATES: -ipccw accounts for DIII-D toroidal angle convention
  flow_mom = -ipccw*flow_mom

  ! Total transport power (MW) to electrons: pow_e

  ! pow_e_nb  (MW)     : NBI power to electrons
  ! pow_e_rf  (MW)     : RF power to electrons
  ! pow_e_ohm (MW)     : Ohmic power to electrons
  ! pow_e_rad (MW)     : [negative] Radiated power to electrons
  ! powe_ion_exp (MW)  : [negative] Neutral ionization power to electrons
  ! powe_wdot_exp (MW) : 
  ! pow_ei    (MW)     : Exchange power from electrons to ions 
  ! pow_e_fus (MW)     : Fusion power to electrons

  pow_e_aux(:) = &
       +pow_e_nb(:) &
       +pow_e_rf(:) &
       +pow_e_ohm(:) &
       +powe_ion_exp(:) 

  pow_e(:) = &
       -0.0*powe_wdot_exp(:) &
       +pow_e_aux(:) &
       +pow_e_rad(:) &
       -pow_ei(:) &
       +pow_e_fus(:)

  ! Total transport power (MW) to ions: pow_i

  ! pow_i_nb (MW)      : NBI power to ions
  ! pow_i_rf (MW)      : RF power to ions
  ! powi_ion_exp (MW)  : [positive] Neutral ionization power to ions
  ! powi_cx_exp (MW)   : [negative] Neutral charge exch power to ions
  ! powi_wdot_exp (MW) 
  ! pow_ei (MW)        : Exchange power to ions
  ! pow_i_fus (MW)     : Fusion power to ions

  pow_i_aux(:) = &
       +pow_i_nb(:) &
       +pow_i_rf(:) &
       +powi_ion_exp(:) &
       +powi_cx_exp(:) 

  pow_i(:) = &
       -0.0*powi_wdot_exp(:) &
       +pow_i_aux(:) &
       +pow_ei(:) &
       +pow_i_fus(:)

  !---------------------------------------------------------
  ! Map profile data into EXPRO interface variables
  !
  EXPRO_n_exp = nx
  call EXPRO_alloc('./',1)
  !
  EXPRO_rho  = rho
  EXPRO_rmin = rmin(:)
  EXPRO_rmaj = rmaj(:)
  ! COORDINATES: set sign of q
  EXPRO_q = abs(q(:))*ipccw*btccw
  EXPRO_kappa = kappa(:)
  EXPRO_delta = delta(:)
  EXPRO_te = onetwo_te(:)
  EXPRO_ne = onetwo_ene(:)*1e-19
  EXPRO_z_eff = onetwo_zeff(:)
  EXPRO_flow_mom = flow_mom(:)
  EXPRO_pow_e = pow_e(:)
  EXPRO_pow_i = pow_i(:)
  EXPRO_pow_ei = pow_ei(:)
  EXPRO_zeta = zeta(:)
  EXPRO_flow_beam = flow_beam(:)
  EXPRO_flow_wall = flow_wall_exp(:)
  EXPRO_sbeame = onetwo_sbeame(:)
  EXPRO_sbcx = sbcx_d(:)
  EXPRO_sscxl = onetwo_sscxl(:)
  EXPRO_zmag = zmag(:)
  EXPRO_ptot = p_tot(:) ! Total pressure
  ! COORDINATES: set sign of poloidal flux
  EXPRO_polflux = abs(dpsi(:))*(-ipccw)

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

        print '(a)',"INFO: (prgen) Found fast ion beam species"
        print '(a)',"INFO: (prgen) Modifying fast ion beam temperature to satisfy beam pressure"

        n0 = n0+1

        onetwo_enion_vec(n0,:) = onetwo_enbeam(:,i)*1e-19

        ! Ti: T[keV] = (p/n)[J]/1.6022e-16[J/keV]
        onetwo_tion_vec(n0,:) = onetwo_pressb(:,i)/onetwo_enbeam(:,i)/1.6022e-16

        onetwo_ion_name(n0) = onetwo_nameb(i)

     endif
  enddo

  ! Fast alphas
  if (minval(onetwo_enalp(:)) > epsilon(0.0)) then

     print '(a)',"INFO: (prgen) Found fast alpha species"
     print '(a)',"INFO: (prgen) Modifying fast alpha temperature to satisfy total pressure"

     n0 = n0+1

     onetwo_enion_vec(n0,:) = onetwo_enalp(:)*1e-19

     onetwo_tion_vec(n0,:) = (p_tot(:)-&
          (sum(onetwo_enion_vec(1:n0-1,:)*&
          onetwo_tion_vec(1:n0-1,:),dim=1)*1e19+&
          onetwo_ene(:)*onetwo_te(:))*1.6022e-16)/(onetwo_enalp(:))/1.6022e-16

     onetwo_ion_name(n0) = 'he'

  endif

  onetwo_nion_tot = n0
  if (n0 > 10) then
     print '(a)',"ERROR: (prgen) Too many ions; report to GACODE developers."
     stop
  endif

  ! reorder
  do i=1,10
     ip = reorder_vec(i)
     EXPRO_ni(i,:) = onetwo_enion_vec(ip,:)
     EXPRO_ti(i,:) = onetwo_tion_vec(ip,:)
  enddo

  ! velocities
  EXPRO_vpol = 0.0
  EXPRO_vtor = 0.0

  ! Look for carbon as first impurity, and insert toroidal velocity at theta=0

  allocate(vphi_carbon(nx))

  if (trim(onetwo_namei(1)) == 'c') then
     ! COORDINATES: -ipccw accounts for DIII-D toroidal angle convention
     vphi_carbon(:) = -ipccw*onetwo_angrot(:)*(rmaj(:)+rmin(:))
  else
     vphi_carbon(:) = 0.0
  endif

  ! Insert carbon toroidal velocity
  do i=1,10
     if (reorder_vec(i) == onetwo_nprim+1) then
        EXPRO_vtor(i,:) = vphi_carbon(:)
     endif
  enddo

  ! Additional powers (fusion and radiation)
  ! * for iterdb, put all radiated power in pow_e_line
  EXPRO_pow_e_fus  = pow_e_fus(:)
  EXPRO_pow_i_fus  = pow_i_fus(:)
  EXPRO_pow_e_line = pow_e_rad(:)

  ! Additional powers (external heating)
  EXPRO_pow_e_aux = pow_e_aux(:)
  EXPRO_pow_i_aux = pow_i_aux(:)
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Read the cer file and overlay
  !
  if (cer_file /= "null") then
     allocate(vpolc_exp(nx))
     allocate(vtorc_exp(nx))
     call prgen_read_cer
     EXPRO_w0 = omega0(:)
     do i=1,10
        if (reorder_vec(i) == onetwo_nprim+1) then
           EXPRO_vtor(i,:) = vtorc_exp(:)
           EXPRO_vpol(i,:) = vpolc_exp(:)
        endif
     enddo
  endif
  !---------------------------------------------------------

  !---------------------------------------------------------
  ! Ion name association
  !
  print '(a)','INFO: (prgen) Found these ion species:'
  do i=1,onetwo_nion_tot
     call onetwo_ion_zmass(onetwo_ion_name(i),onetwo_z(i),onetwo_m(i))
     print '(t6,i2,1x,a)',i,trim(onetwo_ion_name(i))
  enddo

  !
  ! Reordering
  !

  print '(a)','INFO: (prgen) Created these species'
  do i=1,n_ion_max

     ip = reorder_vec(i)

     if (ip > onetwo_nion_tot) then
        print '(t6,i2,1x,a)',i,'[null]'
     else
        if (ip > onetwo_nion) then
           ion_type(i) = type_fast
        else
           ion_type(i) = type_therm
        endif
        print '(t6,i2,1x,a,1x,a)',i,trim(onetwo_ion_name(ip)),ion_type(i)
        ion_mass(i) = onetwo_m(ip)             
        ion_z(i)    = onetwo_z(ip)             

        call prgen_ion_name(nint(ion_mass(i)),ion_z(i),ion_name(i))     
    
     endif
  enddo
  !---------------------------------------------------------

end subroutine prgen_map_iterdb

!------------------------------------------------------------------
! Routine to perform Volume integration
!------------------------------------------------------------------

subroutine volint(f,fdv)

  use prgen_globals, &
       only : onetwo_rho_grid,pi,onetwo_R0,onetwo_hcap,nx

  implicit none

  integer :: i
  real, intent(in) :: f(nx)
  real, intent(out) :: fdv(nx)
  real :: drho
  real :: dvoldr_p,dvoldr_m

  fdv(1) = 0.0

  do i=2,nx

     drho = onetwo_rho_grid(i)-onetwo_rho_grid(i-1)

     dvoldr_p = 2.0*pi*onetwo_rho_grid(i)*2.0*pi*onetwo_R0*onetwo_hcap(i)
     dvoldr_m = 2.0*pi*onetwo_rho_grid(i-1)*2.0*pi*onetwo_R0*onetwo_hcap(i-1)

     fdv(i) = fdv(i-1)+0.5*(dvoldr_p*f(i)+dvoldr_m*f(i-1))*drho

  enddo

end subroutine volint

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
