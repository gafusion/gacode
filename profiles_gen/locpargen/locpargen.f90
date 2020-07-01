!------------------------------------------------------------
! locpargen.f90
!
! PURPOSE:
!  Generate list of local parameters given input radius or 
!  input rho.
!------------------------------------------------------------

program locpargen

  use locpargen_globals
  use expro
  use expro_locsim_interface

  implicit none

  open(unit=1,file='input.locpargen',status='old')
  read(1,*) r0
  read(1,*) rho0
  read(1,*) psi0
  read(1,*) hasgeo
  read(1,*) qnflag
  read(1,*) appendflag
  close(1)

  expro_ctrl_quasineutral_flag = qnflag
  ! We don't need the numerical eq. flag set for this routine.
  expro_ctrl_numeq_flag = hasgeo

  call expro_read('input.gacode') 

  ! Minor radius
  a = expro_rmin(expro_n_exp)

  ! Electron index
  ise = expro_n_ion+1

  if (rho0 > 0.0) then

     ! Use local rho

     x(1) = rho0
     call cub_spline(expro_rho,expro_rmin/a,expro_n_exp,x,y,1)
     r0 = y(1)

  else if (psi0 > 0.0) then

     ! Use local psi_N

     x(1) = psi0*abs(expro_polflux(expro_n_exp))
     call cub_spline(abs(expro_polflux),expro_rmin/a,expro_n_exp,x,y,1)
     r0 = y(1)

  endif

  call expro_locsim_profiles(&
       hasgeo,&
       0,&
       qnflag,&
       expro_n_ion+1,&
       r0,&
       btccw,&
       ipccw,&
       a)

  !------------------------------------------------------------
  ! Create input.geo with local parameters for general geometry
  !
  !if (hasgeo == 1) call locpargen_geo
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
  print 10,'INFO: ----->  n=1: ky*rhos =',q_loc/rmin_loc*rhos_loc/a

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
  
  lambda_star = 7.43 * sqrt((1e3*temp_loc(ise))/(1e13*dens_loc(ise)))/rhos_loc

  tag(:) = (/'1','2','3','4','5','6','7','8','9'/)

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
  call fileopen('input.tglf.locpargen_stack') ; call locpargen_tglf_stack
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
     write(1,*) '#---------------------------------------------------------'
  endif

end subroutine fileopen
