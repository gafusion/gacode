!------------------------------------------------------------
! locpargen.f90
!
! PURPOSE:
!  Generate list of local parameters given input radius or 
!  input rho.
!------------------------------------------------------------

program locpargen

  use locpargen_globals
  use EXPRO_interface
  use EXPRO_locsim_interface

  implicit none

  open(unit=1,file='input.locpargen',status='old')
  read(1,*) r0
  read(1,*) rho0
  read(1,*) psi0
  read(1,*) hasgeo
  read(1,*) qnflag
  read(1,*) appendflag
  close(1)

  EXPRO_ctrl_quasineutral_flag = qnflag
  ! We don't need the numerical eq. flag set for this routine.
  EXPRO_ctrl_numeq_flag = hasgeo

  call EXPRO_alloc('./',1) 
  call EXPRO_read

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

     x(1) = psi0*abs(EXPRO_polflux(EXPRO_n_exp))
     call cub_spline(abs(EXPRO_polflux),EXPRO_rmin/a,EXPRO_n_exp,x,y,1)
     r0 = y(1)

  endif

  call EXPRO_alloc('./',0) 

  call EXPRO_locsim_profiles('./',&
       -1,&
       hasgeo,&
       0,&
       qnflag,&
       EXPRO_n_ion+1,&
       r0,&
       btccw,&
       ipccw,&
       a)

  !------------------------------------------------------------
  ! Create input.geo with local parameters for general geometry
  !
  if (hasgeo == 1) call locpargen_geo
  !------------------------------------------------------------

  if (qnflag == 0) then 
     print 10,'INFO: (locpargen) Quasineutrality NOT enforced.'
  else
     print 10,'INFO: (locpargen) Quasineutrality enforced.'
  endif

  print 10,'INFO: (locpargen) rhos/a   =',rhos_loc/a
  !print 10,'rhoi/a   =',rhos_loc/a*sqrt(temp_loc(ise)/temp_loc(1))
  print 10,'INFO: (locpargen) Te [keV] =',temp_loc(ise)
  print 10,'INFO: (locpargen) Ti [keV] =',temp_loc(1)
  print 10,'INFO: (locpargen) Bunit    =',b_unit_loc
  print 10,'INFO: (locpargen) beta_*   =',beta_star_loc

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

  betae_unit = 4.027e-3*dens_loc(ise)*temp_loc(ise)/b_unit_loc**2

  tag(:) = (/'1','2','3','4','5'/)

  open(unit=1,file='out.locpargen',status='replace')
  write(1,*) q_loc
  write(1,*) r0
  write(1,*) rmaj_loc
  write(1,*) shift_loc
  write(1,*) kappa_loc
  write(1,*) q_loc*rhos_loc*sqrt(temp_loc(ise)/temp_loc(1))/a
  close(1)

  call fileopen('input.cgyro.locpargen') ; call locpargen_cgyro
  call fileopen('input.tglf.locpargen')  ; call locpargen_tglf
  call fileopen('input.neo.locpargen')   ; call locpargen_neo
  print 10,'INFO: (locpargen) Wrote input.*.locpargen'

10 format(a,sp,1pe12.5)

end program locpargen

subroutine fileopen(fname)

  use locpargen_globals

  implicit none

  character(len=*) :: fname

  if (appendflag == 0) then
     open(unit=1,file=trim(fname),status='replace')
  else
     open(unit=1,file=trim(fname),position='append')
     write(1,*) '# **************** CUT HERE **********************'
  endif

end subroutine fileopen
