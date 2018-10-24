!------------------------------------------------------------
! locpargen.f90
!
! PURPOSE:
!  Generate list of local parameters given input radius or 
!  input rho.
!------------------------------------------------------------

program locpargen

  use EXPRO_interface
  use EXPRO_locsim_interface

  implicit none

  integer :: j1,j2
  integer :: is,ise
  real :: r0
  real :: rho0
  real :: psi0
  real :: a
  real, dimension(1) :: x,y
  integer :: hasgeo
  real :: btccw,ipccw
  real :: cc,loglam,nu_ee,pi,betae_unit
  
  character(len=1) :: tag(5)

  open(unit=1,file='input.locpargen',status='old')
  read(1,*) r0
  read(1,*) rho0
  read(1,*) psi0
  read(1,*) hasgeo
  close(1)

  EXPRO_ctrl_quasineutral_flag = 0
  ! We don't need the numerical eq. flag set for this routine.
  EXPRO_ctrl_numeq_flag = hasgeo

  call EXPRO_alloc('./',1) 
  call EXPRO_read

  print '(a)','INFO: (locpargen) Local input parameters:'
  print *

  ! Minor radius
  a = EXPRO_rmin(EXPRO_n_exp)

  ! Electron index
  ise = EXPRO_n_ion+1

  if (rho0 > 0.0) then

     ! Use local rho

     x(1) = rho0
     call cub_spline(EXPRO_rho,EXPRO_rmin/a,EXPRO_n_exp,x,y,1)
     r0 = y(1)

  else if (psi0 > 0.0) then

     ! Use local psi_N

     x(1) = psi0*EXPRO_polflux(EXPRO_n_exp)
     call cub_spline(EXPRO_polflux,EXPRO_rmin/a,EXPRO_n_exp,x,y,1)
     r0 = y(1)

  endif

  !------------------------------------------------------------
  ! Create input.geo with local parameters for general geometry
  !
  if (hasgeo == 1) call locpargen_geo
  !------------------------------------------------------------

  call EXPRO_alloc('./',0) 

  call EXPRO_locsim_profiles('./',&
       -1,&
       0,&
       0,&
       0,&
       EXPRO_n_ion+1,&
       r0,&
       btccw,&
       ipccw,&
       a)

  print 10,'# rhos/a=',rhos_loc/a

  print 10,'RMIN=',r0
  print 10,'RMAJ=',rmaj_loc

  print *
  print 10,'SHIFT=',shift_loc
  print 10,'ZMAG=',zmag_loc
  print 10,'DZMAG=',dzmag_loc
  print 10,'Q=',q_loc
  print 10,'S=',s_loc
  print 10,'KAPPA=',kappa_loc
  print 10,'S_KAPPA=',s_kappa_loc
  print 10,'DELTA=',delta_loc
  print 10,'S_DELTA=',s_delta_loc
  print 10,'ZETA=',zeta_loc
  print 10,'S_ZETA=',s_zeta_loc

  print *
  print 10,'GAMMA_E=',gamma_e_loc*a/cs_loc
  print 10,'GAMMA_P=',gamma_p_loc*a/cs_loc
  print 10,'MACH=',mach_loc/cs_loc

  ! Compute collision frequency
  !
  !          4 pi ne e^4 ne ln(Lambda)
  !  nu_ei = ------------------------- 
  !            me^(1/2) (2 Te)^(3/2)

  pi = 4.0*atan(1.0)
  cc = sqrt(2.0)*pi*charge_norm_fac**4/(4.0*pi*8.8542)**2 &
       *1e9/(sqrt(mass_deuterium)*temp_norm_fac**1.5)

  loglam = 24.0-log(sqrt(dens_loc(ise)*1e13)/(temp_loc(ise)*1e3))
  nu_ee  = cc*loglam*dens_loc(ise)/(sqrt(mass_loc(ise)/2.0)*temp_loc(ise)**1.5)

  print *
  print 10,'NU_EE=',nu_ee*a/cs_loc

  betae_unit = 4.027e-3*dens_loc(ise)*temp_loc(ise)/b_unit_loc**2
  print 10,'BETAE_UNIT=',betae_unit

  !---------------------------------------------------------
  ! Species data
  
  print *
  print 11,'N_SPECIES=',ise

  tag(:) = (/'1','2','3','4','5'/)
  do is=1,ise
     print *
     print 11,'Z_'//tag(is)//'=',int(z_loc(is))
     ! Deuteron mass normalization
     print 10,'MASS_'//tag(is)//'=',mass_loc(is)/2.0
     print 10,'DENS_'//tag(is)//'=',dens_loc(is)/dens_loc(ise)
     print 10,'TEMP_'//tag(is)//'=',temp_loc(is)/temp_loc(ise)
     print 10,'DLNNDR_'//tag(is)//'=',dlnndr_loc(is)
     print 10,'DLNTDR_'//tag(is)//'=',dlntdr_loc(is)
     print 10,'SDLNNDR_'//tag(is)//'=',sdlnndr_loc(is)
     print 10,'SDLNTDR_'//tag(is)//'=',sdlntdr_loc(is)
  enddo

10 format(a,sp,1pe12.5)
11 format(a,i0)

end program locpargen
